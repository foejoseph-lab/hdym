"""offdiag_escape.py -- nonlinear escape times of the FULL nine-component Lee-Wick system seeded with
off-diagonal perturbations of the Dahlqvist-Russberg orbit: the nonlinear counterpart of the typed
Floquet scan (dr_offdiag_floquet_typed.py, B001).

Model: ym/hdym_full_derive.py (all nine A^a_i, temporal gauge, g = M = 1; its CHECK 3 shows the
linearisation about the diagonal background is exactly rung1_system.pkl's K, so this is the same system
the Floquet scans analysed, now without the linearisation).  Integrator: ym/csrc/hdym_full_escape.c
(Dormand-Prince 5(4), adaptive, OpenMP).  Nothing here is symplectic or certified; the Ostrogradsky
energy and the Gauss constraints are recorded per trajectory as the drift diagnostics.

Initial data.  The diagonal sector is placed ON the DR orbit at energy E (same point and scaling as
dr_offdiag_floquet.orbit) with q2 = -grad V + w, q3 = -Hess V q1 + w1 (qk = k-th time derivative), where
the diagonal ghost (w, w1) is a random vector of amplitude delta_diag (default = delta; --delta-diag 0
switches it off, but then the integrator's truncation error seeds it at ~rtol anyway: the ghost-free
manifold is Floquet-unstable on the orbit for E >~ 0.315 and in tongues below, and an unseeded run at
E = 1 escapes through the diagonal ghost at t ~ 120 from a 1e-10 truncation seed).  So the experiment is
a race between the two sectors seeded at equal amplitude; the escaping sector and the growth curves say
who won.  The diagonal ghost's own Floquet rate on the orbit is computed per E (6-dim Hill monodromy) and
stored as rate_diag.  The off-diagonal state (u and its first three derivatives for the six components,
24 numbers) is a random vector of norm delta, projected onto the linearised Gauss constraint (G of
rung1_system.pkl at the DR point); the full nonlinear Gauss residual is then O(delta^2) and is printed.
--block restricts the seed to one block.

Sections
  --section E       nE energies (log-spaced in [Emin, Emax]) x nseed random seeds each.  Observable:
                    t_esc(E, seed), the escaping sector, and the sector norms at nchk checkpoints (the
                    nonlinear growth curves).  Compare with the linear prediction
                    t_lin = ln(R / delta) / rate from the typed log (drawn by --plot if the log is present).
  --section orient  fixed E; a 2D grid over the sphere angles (theta, phi) of the seed direction in the
                    span of three fixed random Gauss-projected directions e1, e2, e3 (seeded by --seed):
                    seed = delta (cos theta e1 + sin theta cos phi e2 + sin theta sin phi e3).  The
                    escape-time map in the style of hdym_escape_map's --section ghost, for the off-diagonal
                    sector.

Linear transients can be large: near a Krein band edge the two block-12 modes are nearly parallel and their beat
amplifies the seed by up to ~10^3 without instability (E = 0.05: seed 1e-4 -> 0.23 at t ~ 2150 and back, energy
conserved).  The growth curves (chk) show it; the "driver" of a runaway is therefore read at norm 1, not at the seed.
A blow-up too sharp to follow to R (step underflow or overflow with the state already beyond |A| = 10) counts as
an escape at the last resolved time and is flagged status = 1 ("singular"; the Lee-Wick runaway can reach a
finite-time singularity, cf. hdym_core.confirm_runaway); status = 2 is a real integrator failure and fails the run.

Output: <results>/<tag>.npz, written after every row (E value, or theta row) so a killed run keeps what
it has (a partial file carries rows_done < rows).  Log: CHECK lines first (build, integrator
cross-check against scipy DOP853 on one trajectory, Gauss residual of the seeds), then one line per row,
then DONE <seconds>.  --plot <tag> draws the figures from an existing npz (no run).

Usage (Spyder, from the repo root, kernel's own interpreter; the derivation must have been run once):
    import sys; %cd C:\\GitHub\\hdym
    !"{sys.executable}" -B -u ym\\hdym_full_derive.py > islands\\results\\hdym_full_derive.txt 2>&1
    !"{sys.executable}" -B -u islands\\code\\offdiag_escape.py --section E --nE 100 --nseed 32 --tag esc_E > islands\\results\\esc_E.txt 2>&1
    !"{sys.executable}" -B -u islands\\code\\offdiag_escape.py --section orient --E 0.412 --n 96 --tag esc_orient_0412 > islands\\results\\esc_orient_0412.txt 2>&1
    !"{sys.executable}" -B -u islands\\code\\offdiag_escape.py --plot esc_E

Run parameters and why: delta = delta_diag = 1e-4.  Convergence test (E = 0.549, both seeds, same initial
state, escape time): delta 1e-4 -> 151.64 / 151.58 / 151.59 at rtol 1e-10 / 1e-11 / 1e-12 against 151.61 from
scipy DOP853 at 1e-12 (converged, 4e-4); delta 1e-5 -> 198.0 / 200.2 / 200.2 vs 200.1 (needs 1e-11);
delta 1e-6 -> 231.3 / 238.8 / 240.9 vs 238.3 (not converged at any practical tolerance: the ghost-free manifold is
unstable and the integrator's truncation error competes with the seed).  Hence delta = 1e-4 with rtol 1e-11,
atol 1e-13 (the escape-time CHECK must pass to 1e-3); ln(1/delta) = 9.2, ln(R/delta) = 18.4 e-foldings.  The
paper's regime (ghost fraction eps = 0.01, |w| ~ 0.1) is larger still.  R = 1e4 (hdym_escape_map); tmax = 3000
(hdym_island at E = 0.3); hmax 0.25; nchk 60.  Cost: ~0.1-1 s per trajectory single-threaded.
Added 2026-10-19.
"""
import argparse
import ctypes
import os
import pickle
import subprocess
import sys
import tempfile
import time

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'ym'))
import dr_offdiag_floquet as base            # noqa: E402  (orbit, PKL, BLOCKS)
import cbuild                                # noqa: E402

OD_M = [1, 2, 3, 5, 6, 7]                    # m = 3 i + a of the off-diagonal components, rung1 ulist order
BLOCK_M = {'12': [1, 3], '13': [2, 6], '23': [5, 7]}
SECTOR_NAMES = ['diagonal', 'block 12', 'block 13', 'block 23']


# ---------------------------------------------------------------- the kernel
class Kernel:
    def __init__(self):
        src = os.path.join(ROOT, 'ym', 'csrc', 'hdym_full_escape.c')
        hdr = os.path.join(ROOT, 'ym', 'csrc', 'hdym_full_rhs.h')
        if not os.path.exists(hdr):
            raise SystemExit("ym/csrc/hdym_full_rhs.h missing: run ym/hdym_full_derive.py first")
        d = tempfile.mkdtemp()
        lib = os.path.join(d, 'hdym_full_escape' + ('.dll' if sys.platform == 'win32' else '.so'))
        cc = os.environ.get('CC', 'gcc')
        inc = '-I' + os.path.join(ROOT, 'ym', 'csrc')
        try:
            cbuild.check_call([cc, '-O3', '-fopenmp', '-shared', '-fPIC', inc, '-o', lib, src, '-lm'])
            omp = 'with OpenMP'
        except subprocess.CalledProcessError:
            cbuild.check_call([cc, '-O3', '-shared', '-fPIC', inc, '-o', lib, src, '-lm'])
            omp = 'single-threaded (no OpenMP)'
        self.lib = ctypes.CDLL(lib)
        dp, ip = ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int)
        self.lib.escape.argtypes = [dp, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, dp, ip, ip, dp, dp]
        self.lib.integrate_to.argtypes = [dp, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        self.lib.omp_threads.restype = ctypes.c_int
        print(f"CHECK kernel built with {cc} {omp}, {self.lib.omp_threads()} threads: PASS", flush=True)

    def escape(self, S0, cj, tmax, R, rtol, atol, hmax, nchk):
        S = np.array(S0, dtype=np.float64, order='C', copy=True)   # (N, 36); the kernel overwrites it with the final state
        N = S.shape[0]
        tesc = np.zeros(N)
        sector = np.zeros(N, dtype=np.int32)
        status = np.zeros(N, dtype=np.int32)
        chk = np.zeros((N, nchk, 6))
        diag = np.zeros((N, 4))
        dp, ip = ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int)
        self.lib.escape(S.ctypes.data_as(dp), N, cj, tmax, R, rtol, atol, hmax, nchk,
                        tesc.ctypes.data_as(dp), sector.ctypes.data_as(ip), status.ctypes.data_as(ip), chk.ctypes.data_as(dp), diag.ctypes.data_as(dp))
        return tesc, sector, chk, diag, S, status

    def integrate_to(self, s0, cj, T, rtol, atol, hmax):
        s = np.ascontiguousarray(s0, dtype=np.float64).copy()
        self.lib.integrate_to(s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), cj, T, rtol, atol, hmax)
        return s


# ---------------------------------------------------------------- the model's numpy side (checks and seeds)
class Model:
    def __init__(self):
        d = pickle.load(open(os.path.join(ROOT, 'ym', 'full_system.pkl'), 'rb'))
        args = list(d['state']) + [d['cj']]
        self.R = sp.lambdify(args, d['R'], 'numpy')
        self.G = sp.lambdify(args, d['G'], 'numpy')
        self.E = sp.lambdify(args, d['E'], 'numpy')
        r1 = pickle.load(open(base.PKL, 'rb'))
        self.Glin = sp.lambdify(list(r1['X']) + list(r1['P']) + [r1['cj']], r1['G'], 'numpy')   # 3 x 24, rung1 state order

    def rhs(self, t, s, cj=1.0):
        return np.concatenate([s[9:36], np.array(self.R(*s, cj), dtype=float)])

    def gauss(self, s, cj=1.0):
        return np.array(self.G(*s, cj), dtype=float)

    def energy(self, s, cj=1.0):
        return float(self.E(*s, cj))


def hessV(q):
    x, y, z = q
    return np.array([[y * y + z * z, 2 * x * y, 2 * x * z], [2 * x * y, z * z + x * x, 2 * y * z], [2 * x * z, 2 * y * z, x * x + y * y]])


def diag_ghost_rate(E):
    """Floquet rate of the colour-diagonal ghost on the DR orbit at E: the Hill system w2 = -(1 + Hess V(q(t))) w,
    6-dim monodromy over one period (as dr_orbit_ghost.py's coarse scan; float)."""
    from scipy.integrate import solve_ivp
    q0, p0, T = base.orbit(E)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        M = s[6:].reshape(6, 6)
        A = np.zeros((6, 6))
        A[:3, 3:] = np.eye(3)
        A[3:, :3] = -(np.eye(3) + hessV(q))
        return np.concatenate([p, -base.gradV(q), (A @ M).ravel()])
    sT = solve_ivp(rhs, (0.0, T), np.concatenate([q0, p0, np.eye(6).ravel()]), **base.RT).y[:, -1]
    mu = np.linalg.eigvals(sT[6:].reshape(6, 6))
    return float(np.log(np.abs(mu).max()) / T)


def seed_diag_ghost(s, delta_diag, rng):
    """add a diagonal ghost of amplitude delta_diag: w = q2 + grad V(q), so q2 += w and q3 += w1, random
    orientation, |w|^2 + |w1|^2 = delta_diag^2 (the (w, W) seeding of hdym_core.full_initial, scaled)."""
    if delta_diag <= 0:
        return s
    v = rng.normal(size=6)
    v *= delta_diag / np.linalg.norm(v)
    s = s.copy()
    for i in range(3):
        s[18 + 4 * i] += v[i]
        s[27 + 4 * i] += v[3 + i]
    return s


def diagonal_state(E):
    """s (36,) on the DR orbit at energy E, ghost-free, no off-diagonal part."""
    q0, p0, T = base.orbit(E)
    g = base.gradV(q0)
    H = hessV(q0)
    s = np.zeros(36)
    for i in range(3):
        s[4 * i] = q0[i]
        s[9 + 4 * i] = p0[i]
        s[18 + 4 * i] = -g[i]
        s[27 + 4 * i] = -(H @ p0)[i]
    return s, T


def od_index():
    """indices into s (36) of the 24 off-diagonal state entries, in rung1 state order (k major, ulist minor)."""
    return np.array([9 * k + m for k in range(4) for m in OD_M])


def gauss_projector(model, s_diag, block):
    """P (24 x 24) projecting an off-diagonal state vector onto ker G_lin, restricted to one block if asked."""
    q = s_diag[[0, 4, 8]]
    p = s_diag[[9, 13, 17]]
    G = np.array(model.Glin(*q, *p, 1.0), dtype=float)      # 3 x 24
    mask = np.ones(24, dtype=bool)
    if block != 'all':
        keep = [j for j, m in enumerate(OD_M) if m in BLOCK_M[block]]
        mask[:] = False
        for k in range(4):
            for j in keep:
                mask[6 * k + j] = True
    Gm = G[:, mask]
    Gm = Gm[np.abs(Gm).sum(axis=1) > 0]                       # rows that touch this block
    P = np.eye(24)
    if len(Gm):
        Pm = np.eye(mask.sum()) - Gm.T @ np.linalg.solve(Gm @ Gm.T, Gm)
        P = np.zeros((24, 24))
        P[np.ix_(mask, mask)] = Pm
    else:
        P = np.diag(mask.astype(float))
    return P


def seeded_state(s_diag, vec24):
    s = s_diag.copy()
    s[od_index()] += vec24
    return s


def pre_escape_drift(ch, dg, lin=1.0):
    """|E(t_k) - E(0)| / |E(0)| at the last checkpoint still in the linear regime (every sector norm and the ghost
    amplitude below `lin`): the integrator's accuracy where it matters for the escape time.  0 if no such checkpoint."""
    E0 = dg[:, 0]
    out = np.zeros(len(E0))
    for n in range(len(E0)):
        okk = ~np.isnan(ch[n, :, 5]) & (np.nanmax(ch[n, :, :5], axis=1) < lin)
        if okk.any():
            k = np.where(okk)[0][-1]
            out[n] = abs(ch[n, k, 5] - E0[n]) / abs(E0[n])
    return out


def driver_sector(ch, tesc=None):
    """who STARTED the runaway: among {block 12, block 13, block 23, diagonal ghost} (chk columns 1..4, all of which
    start at the seed size), the largest at the last checkpoint before escape; returned as 1, 2, 3, or 4
    (= diagonal ghost); -1 for a trajectory that did not escape.  Column 0 (the diagonal positions) is the orbit
    itself, not a deviation, and is not used; `sector` from the kernel says who was largest at |A| = R, a different
    question.  Caveat: linear transients can be large (near a Krein band edge the two block-12 modes are nearly
    parallel and their beat amplifies a 1e-4 seed to 0.2 while Floquet-stable, E = 0.05 in the 19 Oct smoke run),
    so "largest before escape" is the honest reading, not "first past a threshold"."""
    N = ch.shape[0]
    out = np.full(N, -1, dtype=int)
    for n in range(N):
        if tesc is not None and not tesc[n] > 0:
            continue
        valid = ~np.isnan(ch[n, :, 1])
        if valid.any():
            k = np.where(valid)[0][-1]
            out[n] = 1 + int(np.argmax(ch[n, k, 1:5]))
    return out


# ---------------------------------------------------------------- checks
def cross_check(kern, model, E, delta, rng, rtol, atol, hmax):
    from scipy.integrate import solve_ivp
    s_diag, T = diagonal_state(E)
    P = gauss_projector(model, s_diag, 'all')
    v = P @ rng.normal(size=24)
    v *= delta / np.linalg.norm(v)
    s0 = seeded_state(s_diag, v)
    Tc = 2 * T
    sc = kern.integrate_to(s0, 1.0, Tc, rtol, atol, hmax)
    sol = solve_ivp(model.rhs, (0.0, Tc), s0, method='DOP853', rtol=1e-12, atol=1e-14)
    sp_ = sol.y[:, -1]
    rel = np.abs(sc - sp_).max() / np.abs(sp_).max()
    od = od_index()
    rel_od = np.abs(sc[od] - sp_[od]).max() / max(np.abs(sp_[od]).max(), 1e-300)
    print(f"CHECK integrator vs scipy DOP853 at E={E}, t = 2T = {Tc:.2f}: relative difference {rel:.1e} (whole state), "
          f"{rel_od:.1e} (off-diagonal part, |u| = {np.abs(sp_[od]).max():.2e}) (tol 1e-6) {'PASS' if max(rel, rel_od) < 1e-6 else 'FAIL'}", flush=True)
    E0, E1 = model.energy(s0), model.energy(sc)
    G0, G1 = np.linalg.norm(model.gauss(s0)), np.linalg.norm(model.gauss(sc))
    print(f"CHECK energy drift over 2T: {abs(E1 - E0) / abs(E0):.1e} (tol 1e-8) {'PASS' if abs(E1 - E0) / abs(E0) < 1e-8 else 'FAIL'}", flush=True)
    print(f"CHECK Gauss residual of the seed: {G0:.1e} (O(delta^2) = {delta ** 2:.0e}; tol 100 delta^2) "
          f"{'PASS' if G0 < 100 * delta ** 2 else 'FAIL'}; after 2T: {G1:.1e}", flush=True)
    # escape time against scipy, one trajectory with both seeds at an energy where the diagonal ghost is fast
    s_d, _ = diagonal_state(0.549)
    Pd = gauss_projector(model, s_d, 'all')
    v = Pd @ rng.normal(size=24)
    v *= delta / np.linalg.norm(v)
    s1 = seed_diag_ghost(seeded_state(s_d, v), delta, rng)
    te, se, ch, dg, _, st = kern.escape(s1[None, :], 1.0, 600.0, 1e4, rtol, atol, hmax, 24)

    def ev(t, s):
        return 1e4 - np.abs(s[:18]).max()
    ev.terminal = True
    sol = solve_ivp(model.rhs, (0.0, 600.0), s1, method='DOP853', rtol=1e-12, atol=1e-14, events=ev)
    t_sp = sol.t_events[0][0] if len(sol.t_events[0]) else np.nan
    dt_rel = abs(te[0] - t_sp) / t_sp
    print(f"CHECK escape time vs scipy DOP853 (rtol 1e-12) at E=0.549, both seeds {delta:.0e}: kernel {te[0]:.3f}, scipy {t_sp:.3f}, "
          f"relative {dt_rel:.1e} (tol 1e-3) {'PASS' if dt_rel < 1e-3 else 'FAIL'}", flush=True)
    # the diagonal orbit alone must close: ghost-free, no seed
    s_end = kern.integrate_to(s_diag, 1.0, T, rtol, atol, hmax)
    clo = np.abs(s_end[[0, 4, 8, 9, 13, 17]] - s_diag[[0, 4, 8, 9, 13, 17]]).max()
    print(f"CHECK unseeded DR orbit closes after one period in the full system: {clo:.1e} (tol 1e-7) {'PASS' if clo < 1e-7 else 'FAIL'}", flush=True)
    return max(rel, rel_od) < 1e-6 and G0 < 100 * delta ** 2 and clo < 1e-7 and dt_rel < 1e-3


# ---------------------------------------------------------------- sections
def run_E(a, kern, model, rng, out):
    Es = np.logspace(np.log10(a.Emin), np.log10(a.Emax), a.nE)
    N = a.nseed
    tesc = np.full((a.nE, N), np.nan)
    sector = np.full((a.nE, N), -1, dtype=np.int32)
    status = np.full((a.nE, N), -1, dtype=np.int32)
    chk = np.full((a.nE, N, a.nchk, 6), np.nan)
    diag = np.full((a.nE, N, 4), np.nan)
    Tper = np.zeros(a.nE)
    rate_diag = np.full(a.nE, np.nan)
    t0 = time.time()
    for i, E in enumerate(Es):
        s_diag, T = diagonal_state(E)
        Tper[i] = T
        P = gauss_projector(model, s_diag, a.block)
        S0 = np.empty((N, 36))
        for n in range(N):
            v = P @ rng.normal(size=24)
            v *= a.delta / np.linalg.norm(v)
            S0[n] = seed_diag_ghost(seeded_state(s_diag, v), a.delta_diag, rng)
        te, se, ch, dg, _, st = kern.escape(S0, a.cj, a.tmax, a.R, a.rtol, a.atol, a.hmax, a.nchk)
        tesc[i], sector[i], chk[i], diag[i], status[i] = te, se, ch, dg, st
        rate_diag[i] = diag_ghost_rate(E)
        esc = te > 0
        drift = pre_escape_drift(ch, dg)
        secs = ",".join(f"{SECTOR_NAMES[g]}:{(se[esc] == g).sum()}" for g in range(4) if (se[esc] == g).any())
        fail = (te == -2).sum()
        sing = (st == 1).sum()
        print(f"E={E:.5f} T={T:7.3f}  escaped {esc.sum():3d}/{N}{' (%d singular)' % sing if sing else ''}  t_esc median {np.median(te[esc]) if esc.any() else float('nan'):8.1f} "
              f"min {te[esc].min() if esc.any() else float('nan'):8.1f}  sectors [{secs}]  diag-ghost rate {rate_diag[i]:.4f}  max |dE/E| {drift.max():.1e}  "
              f"max |G|_end {dg[:, 3].max():.1e}{'  INTEGRATOR-FAIL x%d' % fail if fail else ''}  ({time.time() - t0:.0f}s)", flush=True)
        np.savez(out, section='E', E=Es, T=Tper, tesc=tesc, sector=sector, status=status, chk=chk, diag=diag, rate_diag=rate_diag, rows_done=i + 1, rows=a.nE,
                 delta=a.delta, delta_diag=a.delta_diag, R=a.R, tmax=a.tmax, rtol=a.rtol, atol=a.atol, block=a.block, cj=a.cj, nchk=a.nchk, seed=a.seed)
    return tesc, sector


def run_orient(a, kern, model, rng, out):
    s_diag, T = diagonal_state(a.E)
    P = gauss_projector(model, s_diag, a.block)
    basis = []
    for _ in range(3):
        v = P @ rng.normal(size=24)
        for b in basis:
            v -= (v @ b) * b
        basis.append(v / np.linalg.norm(v))
    e1, e2, e3 = basis
    rng_diag = np.random.default_rng(a.seed + 1)      # one diagonal-ghost seed orientation per cell, reproducible
    th = np.linspace(0, np.pi, a.n)
    ph = np.linspace(0, 2 * np.pi, a.n, endpoint=False)
    tesc = np.full((a.n, a.n), np.nan)
    sector = np.full((a.n, a.n), -1, dtype=np.int32)
    status = np.full((a.n, a.n), -1, dtype=np.int32)
    diag = np.full((a.n, a.n, 4), np.nan)
    t0 = time.time()
    for i, t_ in enumerate(th):
        S0 = np.empty((a.n, 36))
        for j, p_ in enumerate(ph):
            v = a.delta * (np.cos(t_) * e1 + np.sin(t_) * np.cos(p_) * e2 + np.sin(t_) * np.sin(p_) * e3)
            S0[j] = seed_diag_ghost(seeded_state(s_diag, v), a.delta_diag, rng_diag)
        te, se, ch, dg, _, st = kern.escape(S0, a.cj, a.tmax, a.R, a.rtol, a.atol, a.hmax, a.nchk)
        tesc[i], sector[i], diag[i], status[i] = te, se, dg, st
        esc = te > 0
        print(f"theta={t_:.4f} ({i + 1}/{a.n}): escaped {esc.sum()}/{a.n} ({(st == 1).sum()} singular, {(st == 2).sum()} failed), t_esc in [{te[esc].min() if esc.any() else float('nan'):.1f}, "
              f"{te[esc].max() if esc.any() else float('nan'):.1f}], max |dE/E| {pre_escape_drift(ch, dg).max():.1e}  ({time.time() - t0:.0f}s)", flush=True)
        np.savez(out, section='orient', E=a.E, T=T, theta=th, phi=ph, basis=np.array(basis), tesc=tesc, sector=sector, status=status, diag=diag,
                 rows_done=i + 1, rows=a.n, delta=a.delta, delta_diag=a.delta_diag, R=a.R, tmax=a.tmax, rtol=a.rtol, atol=a.atol, block=a.block, cj=a.cj, seed=a.seed)
    return tesc, sector


# ---------------------------------------------------------------- plots (from an npz; no run)
def linear_rates(results_dir):
    """rate per (E, block) from dr_offdiag_floquet_typed.txt, if present: {block: (E array, rate array, type array)}."""
    import re
    f = os.path.join(results_dir, 'dr_offdiag_floquet_typed.txt')
    if not os.path.exists(f):
        return None
    out = {b: [] for b in base.BLOCKS}
    E = T = None
    fine = False
    for l in open(f):
        if l.startswith('== Fine'):
            fine = True
        if not fine:
            continue
        m = re.match(r"E=([0-9.]+) T=\s*([0-9.]+)", l)
        if m:
            E, T = float(m.group(1)), float(m.group(2))
            continue
        m = re.match(r"\s+blk(\d\d)\s+([0-9.]+)\s+(\S+)", l)
        if m and E is not None:
            out[m.group(1)].append((E, np.log(max(float(m.group(2)), 1.0)) / T, m.group(3)))
    return {b: np.array(v, dtype=object) for b, v in out.items()}


def plot(tag, results_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = np.load(os.path.join(results_dir, tag + '.npz'), allow_pickle=True)
    sec = str(d['section'])
    cols = {0: '#52514e', 1: '#2a78d6', 2: '#eb6834', 3: '#1baf7a'}
    if sec == 'E':
        Es, tesc, sector = d['E'], d['tesc'], d['sector']
        nd = int(d['rows_done'])
        fig, ax = plt.subplots(figsize=(10, 5.5), facecolor='#fcfcfb')
        ax.set_facecolor('#fcfcfb')
        chk = d['chk']
        cols5 = dict(cols)
        cols5[4] = '#9085e9'
        names5 = SECTOR_NAMES + ['diagonal ghost']
        for i in range(nd):
            te = tesc[i]
            drv = driver_sector(chk[i], te)
            for g in range(5):
                m = (drv == g) & (te > 0)
                if m.any():
                    ax.scatter(np.full(m.sum(), Es[i]), te[m], s=9, color=cols5[g], alpha=0.7, lw=0)
            m = te == -1
            if m.any():
                ax.scatter([Es[i]], [float(d['tmax']) * 1.15], s=14, marker='^', color='#c3c2b7', lw=0)
        lr = linear_rates(results_dir)
        if lr:
            ln = np.log(float(d['R']) / float(d['delta']))
            for b, c in (('12', 1), ('13', 2), ('23', 3)):
                arr = lr[b]
                E_ = np.array([r[0] for r in arr], dtype=float)
                rate = np.array([r[1] for r in arr], dtype=float)
                ok = rate > 0
                ax.plot(E_[ok], ln / rate[ok], '.', ms=2.5, color=cols[c], alpha=0.5)
                ax.plot(E_[ok], np.log(1.0 / float(d['delta'])) / rate[ok], '.', ms=1.2, color=cols[c], alpha=0.35)
            ax.plot([], [], '.', color='#0b0b0b', label='linear: ln(R/δ)/rate and ln(1/δ)/rate per block (typed log)')
        if 'rate_diag' in d.files:
            rd = d['rate_diag'][:nd]
            okd = rd > 1e-6
            dd_ = max(float(d['delta_diag']), 1e-10)
            ax.plot(Es[:nd][okd], np.log(float(d['R']) / dd_) / rd[okd], '-', color=cols5[4], lw=1,
                    label='linear, diagonal ghost: ln(R/δ_diag)/rate_diag (and ln(1/δ_diag), dotted)')
            ax.plot(Es[:nd][okd], np.log(1.0 / dd_) / rd[okd], ':', color=cols5[4], lw=1)
        for g in (1, 2, 3, 4):
            ax.scatter([], [], color=cols5[g], label=f'runaway led by {names5[g]} (largest before escape)')
        ax.scatter([], [], marker='^', color='#c3c2b7', label=f'not escaped by tmax = {float(d["tmax"]):.0f}')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('E')
        ax.set_ylabel('escape time')
        ax.set_title(f'Nonlinear escape times, off-diagonal seed δ = {float(d["delta"]):.0e} on the DR orbit, block {str(d["block"])} '
                     f'({nd}/{int(d["rows"])} rows)', fontsize=10, loc='left')
        ax.legend(frameon=False, fontsize=8)
        ax.grid(True, color='#eeede9')
        fig.savefig(os.path.join(results_dir, tag + '.png'), dpi=150, bbox_inches='tight')
        print(f"wrote {tag}.png")
    else:
        te, se = d['tesc'], d['sector']
        fig, axs = plt.subplots(1, 2, figsize=(12, 5), facecolor='#fcfcfb')
        im = axs[0].imshow(np.where(te > 0, np.log10(np.where(te > 0, te, 1)), np.nan), origin='lower', aspect='auto',
                           extent=[0, 2, 0, 1], cmap='viridis')
        plt.colorbar(im, ax=axs[0], label='log10 t_esc')
        axs[0].set_xlabel('φ/π')
        axs[0].set_ylabel('θ/π')
        axs[0].set_title(f'escape time, E = {float(d["E"])}, δ = {float(d["delta"]):.0e}', fontsize=10)
        cm = matplotlib.colors.ListedColormap(['#c3c2b7', cols[0], cols[1], cols[2], cols[3]])
        axs[1].imshow(se + 1, origin='lower', aspect='auto', extent=[0, 2, 0, 1], cmap=cm, vmin=0, vmax=4)
        axs[1].set_xlabel('φ/π')
        axs[1].set_title('escaping sector (grey: not escaped)', fontsize=10)
        fig.savefig(os.path.join(results_dir, tag + '.png'), dpi=150, bbox_inches='tight')
        print(f"wrote {tag}.png")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--section', choices=['E', 'orient'], default='E')
    ap.add_argument('--plot', default=None, help='tag of an existing npz to draw; no run')
    ap.add_argument('--tag', default='offdiag_escape')
    ap.add_argument('--results', default=os.path.join(HERE, '..', 'results'))
    ap.add_argument('--Emin', type=float, default=0.05)
    ap.add_argument('--Emax', type=float, default=1.0)
    ap.add_argument('--nE', type=int, default=100)
    ap.add_argument('--nseed', type=int, default=32)
    ap.add_argument('--E', type=float, default=0.412)
    ap.add_argument('--n', type=int, default=96)
    ap.add_argument('--block', choices=['all', '12', '13', '23'], default='all')
    ap.add_argument('--delta', type=float, default=1e-4)
    ap.add_argument('--delta-diag', type=float, default=1e-4, help='diagonal ghost seed amplitude; 0 = none (then truncation error seeds it at ~rtol)')
    ap.add_argument('--R', type=float, default=1e4)
    ap.add_argument('--tmax', type=float, default=3000.0)
    ap.add_argument('--rtol', type=float, default=1e-11)
    ap.add_argument('--atol', type=float, default=1e-13)
    ap.add_argument('--hmax', type=float, default=0.25)
    ap.add_argument('--nchk', type=int, default=60)
    ap.add_argument('--cj', type=float, default=1.0)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    results = os.path.normpath(a.results)
    os.makedirs(results, exist_ok=True)
    if a.plot:
        plot(a.plot, results)
        return
    t0 = time.time()
    print(f"offdiag_escape.py  python {sys.version.split()[0]}  numpy {np.__version__}  sympy {sp.__version__}  ({sys.executable})")
    print("run parameters: " + " ".join(f"{k}={v}" for k, v in vars(a).items() if k not in ('plot', 'results')))
    print("orbit: DR, X0=%r PX0=%r T_HALF=%r (E=1/2, scaled); diagonal ghost not seeded; seed Gauss-projected (linear G)" % (base.X0, base.PX0, base.T_HALF))
    kern = Kernel()
    model = Model()
    rng = np.random.default_rng(a.seed)
    ok = cross_check(kern, model, a.E if a.section == 'orient' else 0.412, a.delta, rng, a.rtol, a.atol, a.hmax)
    out = os.path.join(results, a.tag + '.npz')
    if a.section == 'E':
        tesc, sector = run_E(a, kern, model, rng, out)
    else:
        tesc, sector = run_orient(a, kern, model, rng, out)
    fails = int((tesc == -2).sum())
    print(f"\nCHECK no integrator failures: {fails} {'PASS' if fails == 0 else 'FAIL'}")
    print(f"{'ALL CHECKS PASS' if ok and fails == 0 else 'SOME CHECK FAILED'}")
    print(f"DONE {time.time() - t0:.0f}")


if __name__ == "__main__":
    main()
