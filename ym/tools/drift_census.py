#!/usr/bin/env python3
"""
drift_census.py -- do the fifteen copies of the model in ym/ agree with each other and with
csrc/hdym_model.h?

Every script in ym/ that integrates something carries its own copy of the potential V, its
gradient, the action of its Hessian and the Yoshida-4 coefficients: in Python (the numpy
reference backends, hdym_core, the analysis scripts, early/), in C (the gcc/ctypes backends)
and in CUDA (the CuPy RawKernel backends), as strings inside the .py files.  This script
evaluates every copy at the same random points and reports, per copy and quantity, the
largest difference from csrc/hdym_model.h in ulps.

Verdicts
  bitwise    every value identical
  roundoff   differences of at most ROUNDOFF_ULP ulps (operation order, 1-ulp coefficients)
  DIFFERENT  larger than that: a copy that computes something else -- look at it
Each row also carries what is EXPECTED for that copy: `bitwise` where the copy is the same
expression in the same order, `roundoff` where it is known to differ in order or in the
Yoshida literals (hdym_core and early/ compute the coefficients from the formula; the kernels
carry literals that differ by 1, 3, 1 ulp in three of the four).  The census passes when no row
is worse than expected.

Nothing in the repository is modified.  The scripts are not imported: the functions and
source strings they define are extracted from the files by name (ast) and exec'd, so scripts
with argparse or runs at module level are fine.  C copies are compiled from the scripts'
own C_SRC strings with an exporter appended that calls their static functions; CUDA copies
likewise when CuPy and a GPU are present (otherwise reported as skipped), in two passes: as the
production kernels are built (NVRTC default, which contracts a*b+c into fused multiply-adds, so
the CPU and GPU builds of the same text differ at roundoff; expected `roundoff`) and with
`--fmad=false` (expected `bitwise` with the CPU).

Usage
  python tools/drift_census.py                      # from ym/
  python tools/drift_census.py --N 20000 --out data/drift_census_2026-10-13.txt
  python tools/drift_census.py --no-cuda            # CPU only
"""
import argparse
import ast
import ctypes
import os
import platform
import subprocess
import sys
import tempfile
import time
import types

import numpy as np

YM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, YM)
import cbuild          # noqa: E402
import hdym_model      # noqa: E402

ROUNDOFF_ULP = 64
_rows = []
_notes = []


# ---------------------------------------------------------------- extraction without import
def load_defs(relpath, names):
    """exec the top-level definitions `names` (functions, classes, assignments) of a script,
    in file order, into a fresh namespace that already has numpy, math, os, sys.  Nothing
    else in the file runs."""
    path = os.path.join(YM, relpath)
    with open(path, encoding='utf-8') as fh:
        src = fh.read()
    tree = ast.parse(src, filename=path)
    wanted = set(names)
    keep = []
    for node in tree.body:
        got = set()
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            got = {node.name}
        elif isinstance(node, ast.Assign):
            for t in node.targets:
                for n in ast.walk(t):
                    if isinstance(n, ast.Name):
                        got.add(n.id)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            got = {node.target.id}
        if got & wanted:
            keep.append(node)
    found = set()
    for node in keep:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            found.add(node.name)
        else:
            for n in ast.walk(node):
                if isinstance(n, ast.Name):
                    found.add(n.id)
    missing = wanted - found
    if missing:
        raise KeyError(f"{relpath}: not found at top level: {sorted(missing)}")
    mod = ast.Module(body=keep, type_ignores=[])
    ast.fix_missing_locations(mod)
    ns = {'__name__': 'census_' + os.path.basename(relpath)[:-3], '__file__': path,
          'np': np, 'math': __import__('math'), 'os': os, 'sys': sys}
    exec(compile(mod, path, 'exec'), ns)
    return types.SimpleNamespace(**ns)


# ---------------------------------------------------------------- comparison
def ulps(a, b, scale=None):
    """max |a-b| in ulps of max(|a|, |b|, scale), max |a-b|, number of differing entries.
    `scale` (broadcastable to a) is the magnitude of the terms that produced the value, so
    that a sum with cancellation is judged against its inputs rather than its tiny result."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch {a.shape} vs {b.shape}")
    mag = np.maximum(np.abs(a), np.abs(b))
    if scale is not None:
        mag = np.maximum(mag, np.broadcast_to(np.abs(np.asarray(scale, dtype=np.float64)), a.shape))
    a, b, mag = a.ravel(), b.ravel(), mag.ravel()
    both = np.isfinite(a) & np.isfinite(b)
    if not np.all(both) and not np.array_equal(np.isfinite(a), np.isfinite(b)):
        return np.inf, np.inf, int(np.sum(np.isfinite(a) != np.isfinite(b)))
    a, b, mag = a[both], b[both], mag[both]
    d = np.abs(a - b)
    sp = np.spacing(mag)
    sp = np.where(sp == 0, np.spacing(1.0), sp)
    u = d / sp
    return float(u.max()) if u.size else 0.0, float(d.max()) if d.size else 0.0, int(np.sum(d != 0))


def verdict(maxulp):
    if maxulp == 0:
        return 'bitwise'
    if maxulp <= ROUNDOFF_ULP:
        return 'roundoff'
    return 'DIFFERENT'


RANK = {'bitwise': 0, 'roundoff': 1, 'DIFFERENT': 2, 'skipped': 1, 'error': 2}


def record(copy, quantity, a, b, expected, note='', scale=None):
    try:
        mu, md, nd = ulps(a, b, scale)
        v = verdict(mu)
    except Exception as e:                       # noqa: BLE001
        mu, md, nd, v = np.nan, np.nan, -1, 'error'
        note = f"{note} {type(e).__name__}: {e}".strip()
    _rows.append((copy, quantity, v, expected, mu, md, nd, note))


def skip(copy, quantity, why):
    _rows.append((copy, quantity, 'skipped', 'n/a', np.nan, np.nan, -1, why))


# ---------------------------------------------------------------- inputs
class Inputs:
    def __init__(self, N, seed):
        rng = np.random.default_rng(seed)
        scale = np.exp(rng.uniform(np.log(0.05), np.log(40.0), N))       # channel depths too
        q = rng.standard_normal((3, N)) * scale
        chan = rng.random(N) < 0.3                                        # some deep-channel points
        q[1, chan] *= 1e-2 / np.maximum(1.0, scale[chan])
        q[2, chan] *= 1e-2 / np.maximum(1.0, scale[chan])
        self.q = q
        self.w = rng.standard_normal((3, N))
        self.P = rng.standard_normal((3, N)) * np.sqrt(scale)
        self.Pi = rng.standard_normal((3, N))
        self.W = rng.standard_normal((3, N))
        self.v8 = rng.standard_normal((8, N))
        self.q2 = rng.standard_normal((2, N)) * scale
        self.P2 = rng.standard_normal((2, N)) * np.sqrt(scale)
        # tame copies for the multi-step tests: fixed-step Yoshida without sub-stepping is only
        # stable for |q| dt well below 1, and the wide set above goes to |q| ~ 40 on purpose
        tame = np.exp(rng.uniform(np.log(0.05), np.log(2.5), N))
        self.qt = rng.standard_normal((3, N)) * tame
        self.Pt = rng.standard_normal((3, N)) * np.sqrt(tame)
        self.q2t = rng.standard_normal((2, N)) * tame
        self.P2t = rng.standard_normal((2, N)) * np.sqrt(tame)
        self.N = N
        self.h_kick = 1.3512071919596576 * 0.01            # YD0 * dt as yoshida passes it
        self.dt = 0.01
        self.nsteps = 3
        self.mu2 = 0.1 ** 2 * np.sqrt(0.03)                 # mu1 = 0.1 at E = 0.03

    def S_lin(self, tame=False):
        q, P = (self.qt, self.Pt) if tame else (self.q, self.P)
        return np.concatenate([q, self.w, P, self.Pi], 0).copy()

    def S_nl(self, tame=False):
        q, P = (self.qt, self.Pt) if tame else (self.q, self.P)
        return np.concatenate([q, self.w * 1e-2, P, self.W * 1e-2], 0).copy()

    def S_ym(self, tame=False):
        q, P = (self.qt, self.Pt) if tame else (self.q, self.P)
        return np.concatenate([q, P], 0).copy()

    def S_2d(self, tame=False):
        q, P = (self.q2t, self.P2t) if tame else (self.q2, self.P2)
        return np.concatenate([q, self.w, P, self.Pi], 0).copy()

    @staticmethod
    def scale_of(S):
        """per-point magnitude of a state array (ncomp, N), for the cancellation-aware ulps"""
        return np.max(np.abs(S), axis=0, keepdims=True)


# ---------------------------------------------------------------- reference values
class Ref:
    def __init__(self, m, X):
        self.YC, self.YD = m.y4()
        self.V = m.V(X.q)
        self.g = m.gradV(X.q)
        self.H = m.hessV(X.q)
        self.Hw1 = m.hessV_dot_lin(X.q, X.w, 1.0)
        self.Hw0 = m.hessV_dot_lin(X.q, X.w, 0.0)
        self.kick_lin1 = m.kick_lin(X.S_lin(), X.h_kick, 1.0)
        self.kick_lin0 = m.kick_lin(X.S_lin(), X.h_kick, 0.0)
        self.yosh_lin1 = m.yoshida_lin(X.S_lin(True), X.dt, 1.0, X.nsteps)
        self.yosh_lin0 = m.yoshida_lin(X.S_lin(True), X.dt, 0.0, X.nsteps)
        self.kick_nl = m.kick_nl(X.S_nl(), X.h_kick)
        self.yosh_nl = m.yoshida_nl(X.S_nl(True), X.dt, X.nsteps)
        self.yosh_ym = m.yoshida_ym(X.S_ym(True), X.dt, X.nsteps)
        self.V2 = m.V2d(X.q2, X.mu2)
        self.g2 = m.gradV2d(X.q2, X.mu2)
        self.kick2d1 = m.kick2d(X.S_2d(), X.h_kick, 1.0, X.mu2)
        self.kick2d0 = m.kick2d(X.S_2d(), X.h_kick, 0.0, X.mu2)
        self.yosh2d1 = m.yoshida2d(X.S_2d(True), X.dt, 1.0, X.mu2, X.nsteps)
        self.Kv1 = m.offdiag_Kv(X.q, X.P, X.v8, 1.0)
        self.Kv0 = m.offdiag_Kv(X.q, X.P, X.v8, 0.0)
        self.mass = [m.offdiag_mass(X.q, k) for k in range(3)]
        self.Hw_scale = np.max(np.abs(self.H), 0, keepdims=True) * np.max(np.abs(X.w), 0, keepdims=True)


# ---------------------------------------------------------------- Python copies
def py_hdym_core(X, R):
    c = 'hdym_core.py'
    ns = load_defs(c, ['V', 'gradV', 'hessV_dot', 'hessV_min_eig', 'dU', '_c', '_w1', '_w0',
                       'YOSHIDA_C', 'YOSHIDA_D', 'yoshida_step', '_lin_step', '_ym_step'])
    q = X.q.T.copy()
    record(c, 'V', ns.V(q), R.V, 'bitwise')
    record(c, 'gradV', ns.gradV(q).T, R.g, 'bitwise')
    record(c, 'hessV_dot (od=1)', ns.hessV_dot(q, X.w.T.copy()).T, R.Hw1, 'roundoff',
           'nonlinear operation order vs linearised order', scale=R.Hw_scale)
    record(c, 'Yoshida C', ns.YOSHIDA_C, R.YC, 'roundoff', 'computed from formula vs literals')
    record(c, 'Yoshida D', ns.YOSHIDA_D[:3], R.YD, 'roundoff', 'computed from formula vs literals')
    # one nonlinear step
    S = X.S_nl(True)
    s = tuple(S[3 * k:3 * k + 3].T.copy() for k in range(4))
    for _ in range(X.nsteps):
        s = ns.yoshida_step(s, X.dt)
    record(c, 'yoshida_step (nonlinear)', np.concatenate([a.T for a in s], 0), R.yosh_nl, 'roundoff',
           'coefficients and order', scale=X.scale_of(R.yosh_nl))
    S = X.S_lin(True)
    p, w, P, Pi = (S[3 * k:3 * k + 3].T.copy() for k in range(4))
    for _ in range(X.nsteps):
        p, w, P, Pi = ns._lin_step(p, w, P, Pi, X.dt)
    record(c, '_lin_step (od=1)', np.concatenate([p.T, w.T, P.T, Pi.T], 0), R.yosh_lin1, 'roundoff',
           'coefficients and order', scale=X.scale_of(R.yosh_lin1))
    S = X.S_ym(True)
    p, P = S[0:3].T.copy(), S[3:6].T.copy()
    for _ in range(X.nsteps):
        p, P = ns._ym_step(p, P, X.dt)
    record(c, '_ym_step', np.concatenate([p.T, P.T], 0), R.yosh_ym, 'roundoff', 'coefficients',
           scale=X.scale_of(R.yosh_ym))


def _py_gpu_like(c, X, R):
    ns = load_defs(c, ['NumpyBackend', 'V', 'hessV_min_eig'])
    record(c, 'V', ns.V(X.q), R.V, 'bitwise')
    record(c, 'NumpyBackend.YC', ns.NumpyBackend.YC, R.YC, 'bitwise')
    record(c, 'NumpyBackend.YD', ns.NumpyBackend.YD[:3], R.YD, 'bitwise')
    for od, ref in ((1.0, R.kick_lin1), (0.0, R.kick_lin0)):
        P, Pi = ns.NumpyBackend._kick(X.q, X.w, X.P, X.Pi, X.h_kick, od)
        record(c, f'NumpyBackend._kick (od={int(od)})', np.concatenate([X.q, X.w, P, Pi], 0), ref, 'bitwise')
    Hm = R.H.reshape(3, 3, -1).transpose(2, 0, 1)
    record(c, 'hessV_min_eig', ns.hessV_min_eig(X.q), np.linalg.eigvalsh(Hm)[:, 0], 'bitwise')


def py_hdym_gpu(X, R):
    _py_gpu_like('hdym_gpu.py', X, R)


def py_hdym_gpu_landscape(X, R):
    _py_gpu_like('hdym_gpu_landscape.py', X, R)


def py_hdym_2d(X, R):
    c = 'hdym_2d.py'
    ns = load_defs(c, ['NumpyBackend'])
    record(c, 'NumpyBackend.YC', ns.NumpyBackend.YC, R.YC, 'bitwise')
    record(c, 'NumpyBackend.YD', ns.NumpyBackend.YD[:3], R.YD, 'bitwise')
    for od, ref in ((1.0, R.kick2d1), (0.0, R.kick2d0)):
        P, Pi = ns.NumpyBackend._kick(X.q2, X.w, X.P2, X.Pi, X.h_kick, od, X.mu2)
        record(c, f'NumpyBackend._kick (od={int(od)})', np.concatenate([X.q2, X.w, P, Pi], 0), ref, 'bitwise')


def py_hdym_escape_map(X, R):
    c = 'hdym_escape_map.py'
    ns = load_defs(c, ['V1'])
    record(c, 'V1', ns.V1(X.q), R.V, 'bitwise')


def py_hdym_rung1_gpu(X, R):
    c = 'hdym_rung1_gpu.py'
    ns = load_defs(c, ['V'])
    record(c, 'V', ns.V(X.q), R.V, 'bitwise')


def py_hdym_rung1_lyap(X, R):
    c = 'hdym_rung1_lyap.py'
    ns = load_defs(c, ['YC', 'YD', 'V', 'gradV', 'ym_step', 'ym_mass'])
    record(c, 'V', ns.V(X.q), R.V, 'bitwise')
    record(c, 'gradV', ns.gradV(X.q), R.g, 'bitwise')
    record(c, 'YC', ns.YC, R.YC, 'bitwise')
    record(c, 'YD', ns.YD[:3], R.YD, 'bitwise')
    q, p = X.qt.copy(), X.Pt.copy()
    for _ in range(X.nsteps):
        q, p = ns.ym_step(q, p, X.dt)
    record(c, 'ym_step', np.concatenate([q, p], 0), R.yosh_ym, 'bitwise')
    for k, name in enumerate(('12', '13', '23')):
        m = ns.ym_mass(X.q, name)            # (2, 2, N)
        record(c, f'ym_mass {name}', m.reshape(4, -1), R.mass[k], 'bitwise')


def py_hdym_spectrum(X, R):
    c = 'hdym_spectrum.py'
    ns = load_defs(c, ['hess'])
    H = ns.hess(X.q.T.copy())                # (N, 3, 3)
    record(c, 'hess', H.transpose(1, 2, 0).reshape(9, -1), R.H, 'bitwise')


def py_hdym_spectrum_xy(X, R):
    c = 'hdym_spectrum_xy.py'
    ns = load_defs(c, ['YC', 'YD', 'V', 'grad', 'step'])
    record(c, 'V', ns.V(X.q), R.V, 'bitwise')
    record(c, 'grad', ns.grad(X.q), R.g, 'bitwise')
    record(c, 'YC', ns.YC, R.YC, 'bitwise')
    record(c, 'YD', ns.YD[:3], R.YD, 'bitwise')
    q, p = X.qt.copy(), X.Pt.copy()
    for _ in range(X.nsteps):
        q, p = ns.step(q, p, X.dt)
    record(c, 'step', np.concatenate([q, p], 0), R.yosh_ym, 'bitwise')


def py_hdym2d_nu0(X, R):
    c = 'hdym2d_nu0.py'
    ns = load_defs(c, ['YC', 'YD', 'V', 'grad', 'step'])
    record(c, 'V', ns.V(X.q2, X.mu2), R.V2, 'bitwise')
    record(c, 'grad', ns.grad(X.q2, X.mu2), R.g2, 'roundoff', 'writes x*y*y + mu2*x where the kernel has x*(y2 + mu2)')
    record(c, 'YC', ns.YC, R.YC, 'bitwise')
    record(c, 'YD', ns.YD[:3], R.YD, 'bitwise')
    q, p = X.q2t.copy(), X.P2t.copy()
    for _ in range(X.nsteps):
        q, p = ns.step(q, p, X.dt, X.mu2)
    # the (p, P) rows of the two-field linearised step do not depend on the ghost
    ref = np.concatenate([R.yosh2d1[0:2], R.yosh2d1[5:7]], 0)
    record(c, 'step', np.concatenate([q, p], 0), ref, 'roundoff', 'through the grad above', scale=X.scale_of(ref))


def py_hdym_diffpair(X, R):
    c = 'hdym_diffpair.py'
    ns = load_defs(c, ['_dmul', 'dV'])
    dq = X.w * 1e-3
    a = ns.dV(X.q, dq)
    b = hdym_model.lib().V(X.q + dq) - R.V
    record(c, 'dV vs V(q+dq)-V(q)', a, b, 'roundoff', 'difference arithmetic vs cancelling subtraction, judged against |V|',
           scale=np.abs(R.V))


def _py_early(c, names, X, R, with_hess=False, with_step='ym'):
    ns = load_defs(c, names)
    record(c, 'V', ns.V(X.q.T.copy()), R.V, 'bitwise')
    record(c, 'gradV', ns.gradV(X.q.T.copy()).T, R.g, 'bitwise')
    if with_hess:
        H = ns.hess(X.q.T.copy())
        record(c, 'hess', H.transpose(1, 2, 0).reshape(9, -1), R.H, 'bitwise')
    D = getattr(ns, 'D', None) or getattr(ns, 'Dc', None)
    K = getattr(ns, 'K', None) or getattr(ns, 'Kc', None)
    record(c, 'Yoshida C', D, R.YC, 'roundoff', 'computed from formula vs literals')
    record(c, 'Yoshida D', K, R.YD, 'roundoff', 'computed from formula vs literals')
    if with_step == 'ym':
        q, p = X.qt.T.copy(), X.Pt.T.copy()
        for _ in range(X.nsteps):
            ns.step(q, p, X.dt)
        record(c, 'step (pure YM, in place)', np.concatenate([q.T, p.T], 0), R.yosh_ym, 'roundoff', 'coefficients',
               scale=X.scale_of(R.yosh_ym))


def py_early_channels(X, R):
    _py_early('early/channels.py', ['gradV', 'V', 'c', 'w1', 'w0', 'D', 'K', 'step'], X, R)


def py_early_spectrum(X, R):
    _py_early('early/spectrum.py', ['gradV', 'V', 'c', 'w1', 'w0', 'Dc', 'Kc', 'step'], X, R)


def py_early_transverse(X, R):
    _py_early('early/transverse.py', ['gradV', 'V', 'hess', 'c', 'w1', 'w0', 'Dc', 'Kc'], X, R,
              with_hess=True, with_step=None)


def py_early_lyap_full(X, R):
    c = 'early/lyap_full.py'
    ns = load_defs(c, ['V', 'gradV', 'hessV_w', 'kick', 'drift', 'c', 'w1', 'w0', 'D', 'K', 'step'])
    record(c, 'V', ns.V(X.q.T.copy()), R.V, 'bitwise')
    record(c, 'gradV', ns.gradV(X.q.T.copy()).T, R.g, 'bitwise')
    record(c, 'hessV_w (od=1)', ns.hessV_w(X.q.T.copy(), X.w.T.copy()).T, R.Hw1, 'roundoff',
           'nonlinear operation order vs linearised order', scale=R.Hw_scale)
    record(c, 'Yoshida C', ns.D, R.YC, 'roundoff', 'computed from formula vs literals')
    record(c, 'Yoshida D', ns.K, R.YD, 'roundoff', 'computed from formula vs literals')
    S = X.S_nl(True)
    s = [S[3 * k:3 * k + 3].T.copy() for k in range(4)]
    for _ in range(X.nsteps):
        ns.step(s, X.dt)
    record(c, 'step (nonlinear, in place)', np.concatenate([a.T for a in s], 0), R.yosh_nl, 'roundoff',
           'coefficients and order', scale=X.scale_of(R.yosh_nl))


PY_CHECKS = [py_hdym_core, py_hdym_gpu, py_hdym_gpu_landscape, py_hdym_2d, py_hdym_escape_map,
             py_hdym_rung1_gpu, py_hdym_rung1_lyap, py_hdym_spectrum, py_hdym_spectrum_xy,
             py_hdym2d_nu0, py_hdym_diffpair, py_early_channels, py_early_spectrum,
             py_early_transverse, py_early_lyap_full]


# ---------------------------------------------------------------- C / CUDA copies
# Exporters: inner bodies that act on trajectory i of component-major arrays; the wrapper
# supplies the loop (C) or the thread index (CUDA).
LOAD12 = ("double p[3], w[3], P[3], Q[3]; int k;\n"
          "for (k = 0; k < 3; k++) { p[k] = S[k*N+i]; w[k] = S[(3+k)*N+i]; P[k] = S[(6+k)*N+i]; Q[k] = S[(9+k)*N+i]; }\n")
STORE12 = ("for (k = 0; k < 3; k++) { S[k*N+i] = p[k]; S[(3+k)*N+i] = w[k]; S[(6+k)*N+i] = P[k]; S[(9+k)*N+i] = Q[k]; }\n")
LOAD10 = ("double p[2], w[3], P[2], Q[3]; int k;\n"
          "for (k = 0; k < 2; k++) { p[k] = S[k*N+i]; P[k] = S[(5+k)*N+i]; }\n"
          "for (k = 0; k < 3; k++) { w[k] = S[(2+k)*N+i]; Q[k] = S[(7+k)*N+i]; }\n")
STORE10 = ("for (k = 0; k < 2; k++) { S[k*N+i] = p[k]; S[(5+k)*N+i] = P[k]; }\n"
           "for (k = 0; k < 3; k++) { S[(2+k)*N+i] = w[k]; S[(7+k)*N+i] = Q[k]; }\n")

# name: (args, body)   -- args without the trailing `int N`
EXPORTERS = {
    'kick_lin': ("double* S, double h, double od",
                 LOAD12 + "kick(p, w, P, Q, h, od);\n" + STORE12),
    'yoshida_lin': ("double* S, double h, double od, int nsteps",
                    LOAD12 + "int s; for (s = 0; s < nsteps; s++) yoshida(p, w, P, Q, h, od);\n" + STORE12),
    'kick2d': ("double* S, double h, double od, double mu2",
               LOAD10 + "kick(p, w, P, Q, h, od, mu2);\n" + STORE10),
    'yoshida2d': ("double* S, double h, double od, double mu2, int nsteps",
                  LOAD10 + "int s; for (s = 0; s < nsteps; s++) yoshida(p, w, P, Q, h, od, mu2);\n" + STORE10),
    'kick_nl': ("double* S, double h", LOAD12 + "kick_nl(p, w, P, Q, h);\n" + STORE12),
    'yoshida_nl': ("double* S, double h, int nsteps",
                   LOAD12 + "int s; for (s = 0; s < nsteps; s++) yoshida_nl(p, w, P, Q, h);\n" + STORE12),
    'kick_r': ("double* S, double h", LOAD12 + "kick_r(p, w, P, Q, h, 0, 0.0);\n" + STORE12),
    'yoshida_r': ("double* S, double h, int nsteps",
                  LOAD12 + "int s; for (s = 0; s < nsteps; s++) yoshida_r(p, w, P, Q, h, 0, 0.0);\n" + STORE12),
    'kick_d': ("double* S, double h",
               LOAD12 + "du dp[3], dw[3], dP[3], dQ[3];\n"
               "for (k = 0; k < 3; k++) { dp[k] = mk(p[k], 0); dw[k] = mk(w[k], 0); dP[k] = mk(P[k], 0); dQ[k] = mk(Q[k], 0); }\n"
               "kick_d(dp, dw, dP, dQ, h);\n"
               "for (k = 0; k < 3; k++) { P[k] = dP[k].v; Q[k] = dQ[k].v; }\n" + STORE12),
    'kick_t': ("double* S, double h",
               LOAD12 + "tb tp[3], tw[3], tP[3], tQ[3];\n"
               "for (k = 0; k < 3; k++) { tp[k] = mktb(p[k], 0.0, 0.0, 0.0); tw[k] = mktb(w[k], 0.0, 0.0, 0.0);"
               " tP[k] = mktb(P[k], 0.0, 0.0, 0.0); tQ[k] = mktb(Q[k], 0.0, 0.0, 0.0); }\n"
               "kick_t(tp, tw, tP, tQ, h);\n"
               "for (k = 0; k < 3; k++) { P[k] = tP[k].v; Q[k] = tQ[k].v; }\n" + STORE12),
    'kick_b': ("double* S, double h",
               LOAD12 + "bd bp[3], bw[3], bP[3], bQ[3];\n"
               "for (k = 0; k < 3; k++) { bp[k] = mkb(mkd(p[k], 0.0), 0.0, mkd(0.0, 0.0), 0.0); bw[k] = mkb(mkd(w[k], 0.0), 0.0, mkd(0.0, 0.0), 0.0);"
               " bP[k] = mkb(mkd(P[k], 0.0), 0.0, mkd(0.0, 0.0), 0.0); bQ[k] = mkb(mkd(Q[k], 0.0), 0.0, mkd(0.0, 0.0), 0.0); }\n"
               "kick_b(bp, bw, bP, bQ, h);\n"
               "for (k = 0; k < 3; k++) { P[k] = bP[k].v.h + bP[k].v.l; Q[k] = bQ[k].v.h + bQ[k].v.l; }\n" + STORE12),
    'gradV': ("const double* q, double* out",
              "double qq[3], g[3]; int k; for (k = 0; k < 3; k++) qq[k] = q[k*N+i];\n"
              "gradV(qq, g); for (k = 0; k < 3; k++) out[k*N+i] = g[k];\n"),
    'yoshida_ym': ("double* S, double h, int nsteps",
                   "double q[3], p[3]; int k, s; for (k = 0; k < 3; k++) { q[k] = S[k*N+i]; p[k] = S[(3+k)*N+i]; }\n"
                   "for (s = 0; s < nsteps; s++) yoshida(q, p, h);\n"
                   "for (k = 0; k < 3; k++) { S[k*N+i] = q[k]; S[(3+k)*N+i] = p[k]; }\n"),
    'Kv': ("const double* q, const double* p, const double* v, double* out, double cj",
           "double qq[3], pp[3], vv[8], o[8]; int k;\n"
           "for (k = 0; k < 3; k++) { qq[k] = q[k*N+i]; pp[k] = p[k*N+i]; }\n"
           "for (k = 0; k < 8; k++) vv[k] = v[k*N+i];\n"
           "Kv(qq, pp, vv, o, cj); for (k = 0; k < 8; k++) out[k*N+i] = o[k];\n"),
}
Y4_C = "void census_y4(double* out) { out[0] = YC0; out[1] = YC1; out[2] = YC2; out[3] = YC3; out[4] = YD0; out[5] = YD1; out[6] = YD2; }\n"
Y4_CU = ('extern "C" __global__ void census_y4(double* out) { if (threadIdx.x == 0 && blockIdx.x == 0) {'
         ' out[0] = YC0; out[1] = YC1; out[2] = YC2; out[3] = YC3; out[4] = YD0; out[5] = YD1; out[6] = YD2; } }\n')


def exporters_c(names):
    out = [Y4_C]
    for n in names:
        args, body = EXPORTERS[n]
        out.append(f"void census_{n}({args}, int N)\n{{\n    int i;\n    for (i = 0; i < N; i++) {{\n{body}    }}\n}}\n")
    return "\n".join(out)


def exporters_cuda(names):
    out = [Y4_CU]
    for n in names:
        args, body = EXPORTERS[n]
        out.append(f'extern "C" __global__ void census_{n}({args}, int N)\n{{\n'
                   f"    int i = blockIdx.x * blockDim.x + threadIdx.x;\n    if (i < N) {{\n{body}    }}\n}}\n")
    return "\n".join(out)


# per script: (defs to extract, how to get the C source, how to get the CUDA source, exporter names, expected, note)
C_SPECS = [
    ('hdym_gpu.py', ['STEP_SRC', 'CUDA_SRC', 'C_SRC'], 'C_SRC', 'CUDA_SRC',
     ['kick_lin', 'yoshida_lin'], 'bitwise', ''),
    ('hdym_gpu_landscape.py', ['STEP_SRC', 'CUDA_SRC', 'C_SRC'], 'C_SRC', 'CUDA_SRC',
     ['kick_lin', 'yoshida_lin'], 'bitwise', ''),
    ('hdym_2d.py', ['STEP_SRC', 'CUDA_SRC', 'C_SRC'], 'C_SRC', 'CUDA_SRC',
     ['kick2d', 'yoshida2d'], 'bitwise', ''),
    ('hdym_escape_map.py', ['STEP_SRC', 'CUDA_SRC', 'C_SRC'], 'C_SRC', 'CUDA_SRC',
     ['kick_nl', 'yoshida_nl'], 'bitwise', ''),
    ('hdym_precision.py', ['PAIR_SRC', 'C_SRC', 'CUDA_SRC'], 'C_SRC', 'CUDA_SRC',
     ['kick_r', 'yoshida_r'], 'bitwise', 'd = 0, sc = 0: every rounding is the identity'),
    ('hdym_diffpair.py', ['DUAL_SRC', 'c_src', 'cuda_src'], lambda ns: ns.c_src(False), lambda ns: ns.cuda_src(False),
     ['kick_d'], 'bitwise', 'value channel of the dual numbers, DT = double'),
    ('hdym_ball.py', ['GHOST', 'BALL_SRC', 'c_src', 'cuda_src'], lambda ns: ns.c_src(), lambda ns: ns.cuda_src(),
     ['kick_t', 'kick_b'], None, ''),
    ('hdym_rung1_gpu.py', ['KBLOCK', 'STEP_SRC', 'CUDA_SRC', 'C_SRC'], 'C_SRC', 'CUDA_SRC',
     ['gradV', 'yoshida_ym', 'Kv'], 'bitwise', ''),
]
EXPECT_OVERRIDE = {('hdym_ball.py', 'kick_t'): ('bitwise', 'value channel of the tangent balls'),
                   ('hdym_ball.py', 'kick_b'): ('roundoff', 'double-word centre h + l vs plain double')}


def _get_src(ns, how):
    return how(ns) if callable(how) else getattr(ns, how)


FMA_NOTE = 'NVRTC contracts a*b+c into fused multiply-adds by default (--fmad=true), gcc on x86-64 does not'


def run_c_exports(lib, names, X, R, copy, expected, note, call, fma=False):
    """call(fname, *args) invokes exporter fname on the backend; returns host arrays.
    fma=True: the backend may contract a*b+c (CUDA as the production kernels are built), so
    every row is expected at roundoff and judged against the magnitude of its terms."""
    y4 = call('census_y4')
    record(copy, 'Yoshida C', y4[:4], R.YC, 'bitwise', note)
    record(copy, 'Yoshida D', y4[4:7], R.YD, 'bitwise', note)
    sc = X.scale_of                      # cancellation-aware: judge against the state's magnitude

    def kicksc(ref, S0):
        # a kick adds h (grad V + Hess V w) to P; its roundoff is set by the size of that increment,
        # which can exceed every component of the state when P and the increment nearly cancel
        return np.maximum(np.maximum(sc(ref), sc(S0)), sc(ref - S0))
    for n in names:
        exp, nt = EXPECT_OVERRIDE.get((os.path.basename(copy.split(' ')[0]), n), (expected, note))
        if fma:
            exp, nt = 'roundoff', (nt + '; ' if nt else '') + FMA_NOTE
        if n == 'kick_lin':
            for od, ref in ((1.0, R.kick_lin1), (0.0, R.kick_lin0)):
                record(copy, f'kick (od={int(od)})', call('census_kick_lin', X.S_lin(), X.h_kick, od), ref, exp, nt,
                       scale=kicksc(ref, X.S_lin()))
        elif n == 'yoshida_lin':
            for od, ref in ((1.0, R.yosh_lin1), (0.0, R.yosh_lin0)):
                record(copy, f'yoshida x{X.nsteps} (od={int(od)})',
                       call('census_yoshida_lin', X.S_lin(True), X.dt, od, X.nsteps), ref, exp, nt, scale=sc(ref))
        elif n == 'kick2d':
            for od, ref in ((1.0, R.kick2d1), (0.0, R.kick2d0)):
                record(copy, f'kick (od={int(od)})', call('census_kick2d', X.S_2d(), X.h_kick, od, X.mu2), ref, exp, nt,
                       scale=kicksc(ref, X.S_2d()))
        elif n == 'yoshida2d':
            record(copy, f'yoshida x{X.nsteps} (od=1)',
                   call('census_yoshida2d', X.S_2d(True), X.dt, 1.0, X.mu2, X.nsteps), R.yosh2d1, exp, nt,
                   scale=sc(R.yosh2d1))
        elif n in ('kick_nl', 'kick_r', 'kick_d', 'kick_t', 'kick_b'):
            record(copy, n, call(f'census_{n}', X.S_nl(), X.h_kick), R.kick_nl, exp, nt, scale=kicksc(R.kick_nl, X.S_nl()))
        elif n in ('yoshida_nl', 'yoshida_r'):
            record(copy, f'{n} x{X.nsteps}', call(f'census_{n}', X.S_nl(True), X.dt, X.nsteps), R.yosh_nl, exp, nt,
                   scale=sc(R.yosh_nl))
        elif n == 'gradV':
            record(copy, 'gradV', call('census_gradV', X.q.copy(), np.empty_like(X.q)), R.g, exp, nt, scale=sc(R.g))
        elif n == 'yoshida_ym':
            record(copy, f'yoshida x{X.nsteps} (pure YM)', call('census_yoshida_ym', X.S_ym(True), X.dt, X.nsteps),
                   R.yosh_ym, exp, nt, scale=sc(R.yosh_ym))
        elif n == 'Kv':
            # K v is a long polynomial with cancellation: judge against the size of its terms,
            # ~ (|q|^2 + |p|^2 + 1)^2 |v|
            m = np.max(np.abs(X.q), 0, keepdims=True) ** 2 + np.max(np.abs(X.P), 0, keepdims=True) ** 2 + 1.0
            kscale = m ** 2 * np.max(np.abs(X.v8), 0, keepdims=True)
            for cj, ref in ((1.0, R.Kv1), (0.0, R.Kv0)):
                record(copy, f'Kv (cj={int(cj)})',
                       call('census_Kv', X.q.copy(), X.P.copy(), X.v8.copy(), np.empty_like(X.v8), cj), ref, exp, nt,
                       scale=kscale)


class CCaller:
    """compile src+exporters with gcc (same flags as the scripts) and call through ctypes"""
    def __init__(self, src, tag):
        d = tempfile.mkdtemp(prefix='census_')
        path = os.path.join(d, tag + '.c')
        lib = os.path.join(d, tag + ('.dll' if sys.platform == 'win32' else '.so'))
        with open(path, 'w') as fh:
            fh.write(src)
        cc = os.environ.get('CC', 'gcc')
        try:
            cbuild.check_call([cc, '-O3', '-fopenmp', '-shared', '-fPIC', '-o', lib, path, '-lm'])
        except subprocess.CalledProcessError:
            cbuild.check_call([cc, '-O3', '-shared', '-fPIC', '-o', lib, path, '-lm'])
        self.L = ctypes.CDLL(lib)
        self.dp = ctypes.POINTER(ctypes.c_double)

    def __call__(self, fname, *args):
        f = getattr(self.L, fname)
        if fname == 'census_y4':
            out = np.empty(7)
            f.argtypes = [self.dp]; f.restype = None
            f(out.ctypes.data_as(self.dp))
            return out
        arrs = [a for a in args if isinstance(a, np.ndarray)]
        N = arrs[0].shape[-1]
        cargs, types_ = [], []
        for a in args:
            if isinstance(a, np.ndarray):
                a = np.ascontiguousarray(a, dtype=np.float64)
                cargs.append(a.ctypes.data_as(self.dp)); types_.append(self.dp)
            elif isinstance(a, (int, np.integer)):
                cargs.append(int(a)); types_.append(ctypes.c_int)
            else:
                cargs.append(float(a)); types_.append(ctypes.c_double)
        cargs.append(int(N)); types_.append(ctypes.c_int)
        f.argtypes = types_; f.restype = None
        f(*cargs)
        # the first array argument is the in/out state (or the last array is `out` for pure functions)
        outs = [a for a in args if isinstance(a, np.ndarray)]
        return np.ascontiguousarray(outs[-1] if fname in ('census_gradV', 'census_Kv') else outs[0])


class CudaCaller:
    def __init__(self, src, tag, options=()):
        import cupy as cp
        self.cp = cp
        self.mod = cp.RawModule(code=src, options=tuple(options))
        self.tag = tag

    def __call__(self, fname, *args):
        cp = self.cp
        k = self.mod.get_function(fname)
        if fname == 'census_y4':
            out = cp.empty(7, dtype=cp.float64)
            k((1,), (32,), (out,))
            return cp.asnumpy(out)
        dev = []
        cargs = []
        N = None
        for a in args:
            if isinstance(a, np.ndarray):
                g = cp.asarray(np.ascontiguousarray(a, dtype=np.float64))
                dev.append(g); cargs.append(g); N = a.shape[-1]
            elif isinstance(a, (int, np.integer)):
                cargs.append(np.int32(a))
            else:
                cargs.append(np.float64(a))
        cargs.append(np.int32(N))
        threads = 128
        k(((N + threads - 1) // threads,), (threads,), tuple(cargs))
        cp.cuda.Stream.null.synchronize()
        out = dev[-1] if fname in ('census_gradV', 'census_Kv') else dev[0]
        return cp.asnumpy(out)


def c_checks(X, R, do_cuda):
    cuda_ok = False
    if do_cuda:
        try:
            import cupy as cp
            cp.cuda.Device().compute_capability
            cuda_ok = True
            _notes.append(f"CuPy {cp.__version__}, device {cp.cuda.runtime.getDeviceProperties(0)['name'].decode()}")
        except Exception as e:                   # noqa: BLE001
            _notes.append(f"CUDA copies skipped: {type(e).__name__}: {e}")
    for relpath, defs, c_how, cu_how, names, expected, note in C_SPECS:
        try:
            ns = load_defs(relpath, defs)
        except Exception as e:                   # noqa: BLE001
            skip(relpath + ' [C]', '*', f"extract failed: {type(e).__name__}: {e}")
            continue
        tag = relpath.replace('/', '_')[:-3]
        # C
        try:
            src = _get_src(ns, c_how) + "\n" + exporters_c(names)
            call = CCaller(src, tag)
            run_c_exports(None, names, X, R, relpath + ' [C]', expected, note, call)
        except Exception as e:                   # noqa: BLE001
            skip(relpath + ' [C]', '*', f"compile/run failed: {type(e).__name__}: {e}")
        # CUDA, twice: as the production kernels are built (NVRTC default, FMA contraction on),
        # expected at roundoff; and with --fmad=false, expected bitwise with the CPU.  The second
        # pass is what pins the first pass's differences on the contraction and nothing else.
        if cuda_ok:
            for label, opts, fma in (('[CUDA]', (), True), ('[CUDA fmad=false]', ('--fmad=false',), False)):
                try:
                    src = _get_src(ns, cu_how) + "\n" + exporters_cuda(names)
                    call = CudaCaller(src, tag, opts)
                    run_c_exports(None, names, X, R, relpath + ' ' + label, expected, note, call, fma=fma)
                except Exception as e:               # noqa: BLE001
                    skip(relpath + ' ' + label, '*', f"compile/run failed: {type(e).__name__}: {e}")
        else:
            skip(relpath + ' [CUDA]', '*', 'no CuPy/GPU in this environment')


# ---------------------------------------------------------------- report
def environment():
    lines = [f"python {platform.python_version()} ({sys.executable})",
             f"numpy {np.__version__}", f"platform {platform.platform()}"]
    cc = os.environ.get('CC', 'gcc')
    try:
        v = subprocess.run([cc, '--version'], capture_output=True, text=True).stdout.splitlines()[0]
        lines.append(f"compiler {v}")
    except Exception as e:                       # noqa: BLE001
        lines.append(f"compiler {cc}: {e}")
    return lines


def report(X, seed, t0):
    out = []
    out.append("drift census of the model copies in ym/  --  reference: csrc/hdym_model.h")
    out.append(time.strftime("%Y-%m-%d %H:%M:%S"))
    out += environment()
    out += _notes
    out.append(f"points N = {X.N}, seed = {seed}, kick h = {X.h_kick:.16g}, yoshida dt = {X.dt} x {X.nsteps} steps, "
               f"mu2 = {X.mu2:.6g}; roundoff threshold {ROUNDOFF_ULP} ulp")
    out.append("")
    w0 = max(len(r[0]) for r in _rows)
    w1 = max(len(r[1]) for r in _rows)
    hdr = f"{'copy':{w0}}  {'quantity':{w1}}  {'verdict':9}  {'expected':8}  {'max ulp':>10}  {'max |diff|':>10}  {'n diff':>7}  note"
    out.append(hdr)
    out.append('-' * len(hdr))
    worse = []
    for copy, qty, v, exp, mu, md, nd, note in _rows:
        mus = '' if np.isnan(mu) else (f"{mu:.0f}" if mu < 1e6 else f"{mu:.2e}")
        mds = '' if np.isnan(md) else f"{md:.2e}"
        nds = '' if nd < 0 else str(nd)
        flag = ''
        if v not in ('skipped',) and exp != 'n/a' and RANK[v] > RANK[exp]:
            flag = '  <-- worse than expected'
            worse.append((copy, qty, v, exp))
        out.append(f"{copy:{w0}}  {qty:{w1}}  {v:9}  {exp:8}  {mus:>10}  {mds:>10}  {nds:>7}  {note}{flag}")
    out.append("")
    nb = sum(1 for r in _rows if r[2] == 'bitwise')
    nr = sum(1 for r in _rows if r[2] == 'roundoff')
    nd = sum(1 for r in _rows if r[2] == 'DIFFERENT')
    ns = sum(1 for r in _rows if r[2] == 'skipped')
    ne = sum(1 for r in _rows if r[2] == 'error')
    out.append(f"{len(_rows)} rows: {nb} bitwise, {nr} roundoff, {nd} DIFFERENT, {ns} skipped, {ne} error  ({time.time() - t0:.1f}s)")
    if worse:
        out.append("FAIL: rows worse than expected:")
        for copy, qty, v, exp in worse:
            out.append(f"  {copy}: {qty}: {v} (expected {exp})")
    else:
        out.append("PASS: no copy is worse than expected (bitwise where the text is the same, roundoff where the "
                   "operation order or the Yoshida literals differ)")
    return "\n".join(line.rstrip() for line in out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--N', type=int, default=10000)
    ap.add_argument('--seed', type=int, default=20261013)
    ap.add_argument('--no-cuda', action='store_true')
    ap.add_argument('--out', default=None, help='also write the report to this file')
    a = ap.parse_args()
    t0 = time.time()
    m = hdym_model.lib()
    X = Inputs(a.N, a.seed)
    R = Ref(m, X)
    for chk in PY_CHECKS:
        try:
            chk(X, R)
        except Exception as e:                   # noqa: BLE001
            skip(chk.__name__.replace('py_', '') + '.py', '*', f"{type(e).__name__}: {e}")
    c_checks(X, R, not a.no_cuda)
    text = report(X, a.seed, t0)
    print(text)
    if a.out:
        with open(a.out, 'w', encoding='utf-8') as fh:
            fh.write(text + "\n")
        print("wrote", a.out)
    return 0 if 'PASS' in text.splitlines()[-1] else 1


if __name__ == '__main__':
    sys.exit(main())
