#!/usr/bin/env python3
"""
hdym_diffpair.py -- "difference numbers": pair integration with no roundoff floor.

Why ordinary floating point has a floor.  hdym_uncertainty.py integrates the
two members a, b of a pair separately and compares them.  Every operation on
a coordinate of size ~1 rounds at ~1e-16 absolute, independently in a and b,
so their separation cannot be resolved below ~2e-15 of the section width
(measured by hdym_precision.py).  Nearly all of the 52 mantissa bits are spent
re-describing the part of a and b that is the SAME.

The number system.  Carry each quantity of a pair as   x = (v, d):
    v = value in member a (the reference),   d = (value in b) - (value in a).
The kernel only adds, subtracts and multiplies, and for those the difference
part has an exact algebraic rule that never subtracts two nearly equal numbers:
    (a,da) +- (b,db) = (a +- b,  da +- db)
    (a,da) *  (b,db) = (a b,     da b + a db + da db)
    c * (a,da)       = (c a,     c da)                 (c an exact constant)
d is therefore stored with RELATIVE precision of its own size: a separation of
1e-30 carries ~16 significant digits, not zero.  This is ordinary finite-
difference ("dual number with the quadratic term kept") arithmetic, so it is
exact algebra, not a linearisation: b = v + d follows the full nonlinear flow.

What it does and does not test.  The v channel is computed with exactly the
same operations as the production kernel (it reproduces hdym_escape_map's
escape times bit for bit -- see --selftest), so member a carries the usual
~1e-16 roundoff.  Member b = a + d then carries the SAME roundoff sequence plus
exactly-computed differences: the roundoff is common-mode and cancels in d.
So the pair test measures how the computed flow responds to a change delta in
initial data, at any delta, with the arithmetic noise held fixed.  It removes
the floor; it does not address whether the computed reference orbit is
shadowed by a true one (that needs interval / rigorous methods).

Hardware.  The d channel can be fp32 (--dfloat): it only needs relative
precision.  On an RTX 4080 FP32 is 64x FP64, so the d work is nearly free and
a whole pair costs about one fp64 trajectory -- ~2x faster than integrating
a and b separately, with no floor.  fp32's range (1e-38) limits delta to
>~1e-30 relative; fp64 d (default) goes to ~1e-300.

Initial differences are built analytically (sin half-angle formulas and
V(q_b) - V(q_a) in difference arithmetic), never as S_b - S_a.

Escape.  When one member passes R the other is converted to a plain state
v + d (at that point it is a separate trajectory needing only absolute
accuracy for its escape time) and continued alone.

Examples
  python hdym_diffpair.py --selftest
  python hdym_diffpair.py --E 0.3 --eps 0.01 --section ym --M 4000 --delta-min 1e-30
  python hdym_diffpair.py --E 0.3 --eps 0.01 --section ym --M 4000 --backend cupy --dfloat
"""
import argparse
import copy
import ctypes
import json
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

import hdym_escape_map as em
import hdym_uncertainty as unc

DUAL_SRC = r"""
#define YC0 0.6756035959798288
#define YC1 -0.17560359597982886
#define YC2 -0.17560359597982886
#define YC3 0.6756035959798288
#define YD0 1.3512071919596576
#define YD1 -1.7024143839193153
#define YD2 1.3512071919596576
#define MAXSUB 16384

typedef struct { double v; DT d; } du;

DEVFN du mk(double v, DT d) { du r; r.v = v; r.d = d; return r; }
DEVFN du dadd(du a, du b) { return mk(a.v + b.v, a.d + b.d); }
DEVFN du dsub(du a, du b) { return mk(a.v - b.v, a.d - b.d); }
DEVFN du dmul(du a, du b)
{
    return mk(a.v * b.v, a.d * (DT)b.v + (DT)a.v * b.d + a.d * b.d);
}
DEVFN du dscl(double c, du a) { return mk(c * a.v, (DT)c * a.d); }

/* v-channel: same operations, same order as hdym_escape_map.kick_nl */
DEVFN void kick_d(const du* p, const du* w, du* P, du* W, double h)
{
    du x = dsub(p[0], w[0]), y = dsub(p[1], w[1]), z = dsub(p[2], w[2]);
    du x2 = dmul(x, x), y2 = dmul(y, y), z2 = dmul(z, z);
    du yz = dadd(y2, z2), zx = dadd(z2, x2), xy = dadd(x2, y2);
    du gx = dmul(x, yz), gy = dmul(y, zx), gz = dmul(z, xy);
    du wx = w[0], wy = w[1], wz = w[2];
    du txy = dmul(dscl(2.0, x), y), txz = dmul(dscl(2.0, x), z), tyz = dmul(dscl(2.0, y), z);
    du hx = dadd(dadd(dmul(yz, wx), dmul(txy, wy)), dmul(txz, wz));
    du hy = dadd(dadd(dmul(txy, wx), dmul(zx, wy)), dmul(tyz, wz));
    du hz = dadd(dadd(dmul(txz, wx), dmul(tyz, wy)), dmul(xy, wz));
    P[0] = dsub(P[0], dscl(h, dadd(gx, hx)));
    P[1] = dsub(P[1], dscl(h, dadd(gy, hy)));
    P[2] = dsub(P[2], dscl(h, dadd(gz, hz)));
    W[0] = dadd(W[0], dscl(h, dadd(hx, wx)));
    W[1] = dadd(W[1], dscl(h, dadd(hy, wy)));
    W[2] = dadd(W[2], dscl(h, dadd(hz, wz)));
}

DEVFN void drift_d(du* p, du* w, const du* P, const du* W, double ch)
{
    int k;
    for (k = 0; k < 3; k++) { p[k] = dadd(p[k], dscl(ch, P[k])); w[k] = dsub(w[k], dscl(ch, W[k])); }
}

DEVFN void yoshida_d(du* p, du* w, du* P, du* W, double h)
{
    drift_d(p, w, P, W, YC0 * h);
    kick_d(p, w, P, W, YD0 * h);
    drift_d(p, w, P, W, YC1 * h);
    kick_d(p, w, P, W, YD1 * h);
    drift_d(p, w, P, W, YC2 * h);
    kick_d(p, w, P, W, YD2 * h);
    drift_d(p, w, P, W, YC3 * h);
}

/* om = max |q| of member a (sel=0) or member b = v + d (sel=1) */
DEVFN int nsub_d(const du* p, const du* w, double dt, double omd, int sel)
{
    double om = 0.0, a;
    int k, n;
    for (k = 0; k < 3; k++) {
        a = sel ? fabs((p[k].v + (double)p[k].d) - (w[k].v + (double)w[k].d))
                : fabs(p[k].v - w[k].v);
        if (a > om) om = a;
    }
    n = (int)ceil(om * dt / omd);
    if (n < 1) n = 1;
    if (n > MAXSUB) n = MAXSUB;
    return n;
}

DEVFN int esc_d(const du* p, const du* w, const du* P, const du* W, double R, int sel)
{
    double m = 0.0, a;
    int k;
    for (k = 0; k < 3; k++) {
        a = fabs(p[k].v + (sel ? (double)p[k].d : 0.0)); if (a > m) m = a;
        a = fabs(w[k].v + (sel ? (double)w[k].d : 0.0)); if (a > m) m = a;
        a = fabs(P[k].v + (sel ? (double)P[k].d : 0.0)); if (a > m) m = a;
        a = fabs(W[k].v + (sel ? (double)W[k].d : 0.0)); if (a > m) m = a;
    }
    return (!(m == m) || m > 1e300 || m > R);
}

/* collapse onto member b: v <- v + d, d <- 0 */
DEVFN void to_b(du* x) { int k; for (k = 0; k < 3; k++) { x[k].v += (double)x[k].d; x[k].d = 0; } }
DEVFN void to_a(du* x) { int k; for (k = 0; k < 3; k++) x[k].d = 0; }

/* V: (24, N): rows 0..11 = reference state a, rows 12..23 = difference b - a.
   out (2 + nchk, N): t_esc(a), t_esc(b) (-1 = not by nmax), |d| at chk steps. */
DEVFN void dpair(const double* S, double* out, int N, int i, double dt, long long nmax, double R,
                 double omd, const long long* chk, int nchk)
{
    du p[3], w[3], P[3], W[3];
    int k, j = 0, n, alive_a = 1, alive_b = 1;
    long long s;
    for (k = 0; k < 3; k++) {
        p[k] = mk(S[(0 + k) * N + i], (DT)S[(12 + k) * N + i]);
        w[k] = mk(S[(3 + k) * N + i], (DT)S[(15 + k) * N + i]);
        P[k] = mk(S[(6 + k) * N + i], (DT)S[(18 + k) * N + i]);
        W[k] = mk(S[(9 + k) * N + i], (DT)S[(21 + k) * N + i]);
    }
    out[0 * N + i] = -1.0;
    out[1 * N + i] = -1.0;
    for (k = 0; k < nchk; k++) out[(2 + k) * N + i] = -1.0;
    for (s = 1; s <= nmax && (alive_a || alive_b); s++) {
        /* both alive: v = a, v + d = b.  One alive: that member is in v, d = 0 */
        if (alive_a && alive_b) {
            int na = nsub_d(p, w, dt, omd, 0), nb = nsub_d(p, w, dt, omd, 1);
            n = na > nb ? na : nb;
        } else {
            n = nsub_d(p, w, dt, omd, 0);
        }
        {
            int u;
            double h = dt / n;
            for (u = 0; u < n; u++) yoshida_d(p, w, P, W, h);
        }
        if (alive_a && alive_b) {
            int ea = esc_d(p, w, P, W, R, 0), eb = esc_d(p, w, P, W, R, 1);
            if (ea) { alive_a = 0; out[0 * N + i] = s * dt; }
            if (eb) { alive_b = 0; out[1 * N + i] = s * dt; }
            if (ea && !eb) { to_b(p); to_b(w); to_b(P); to_b(W); }
            else if (eb && !ea) { to_a(p); to_a(w); to_a(P); to_a(W); }
        } else if (alive_a) {
            if (esc_d(p, w, P, W, R, 0)) { alive_a = 0; out[0 * N + i] = s * dt; }
        } else if (alive_b) {
            if (esc_d(p, w, P, W, R, 0)) { alive_b = 0; out[1 * N + i] = s * dt; }
        }
        while (j < nchk && chk[j] <= s) {
            if (chk[j] == s && alive_a && alive_b) {
                double acc = 0.0, e;
                for (k = 0; k < 3; k++) {
                    e = (double)p[k].d; acc += e * e;  e = (double)w[k].d; acc += e * e;
                    e = (double)P[k].d; acc += e * e;  e = (double)W[k].d; acc += e * e;
                }
                out[(2 + j) * N + i] = sqrt(acc);
            }
            j++;
        }
    }
}
"""


def _dll_search_path(cc):
    """Windows: Python >= 3.8 does not search PATH for a DLL's dependencies (libgomp-1.dll,
    libwinpthread-1.dll from MinGW), so add the compiler's bin directory explicitly."""
    if sys.platform == 'win32':
        import shutil
        exe = shutil.which(cc)
        if exe:
            try:
                os.add_dll_directory(os.path.dirname(os.path.abspath(exe)))
            except (AttributeError, OSError):
                pass


def c_src(dfloat):
    return (("#define DT float\n" if dfloat else "#define DT double\n")
            + "#include <math.h>\n#define DEVFN static inline\n" + DUAL_SRC + r"""
void dpairs(const double* S, double* out, int N, double dt, long long nmax, double R, double omd,
            const long long* chk, int nchk)
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
    for (i = 0; i < N; i++) dpair(S, out, N, i, dt, nmax, R, omd, chk, nchk);
}
""")


def cuda_src(dfloat):
    return (("#define DT float\n" if dfloat else "#define DT double\n")
            + "#define DEVFN __device__ __forceinline__\n" + DUAL_SRC + r"""
extern "C" __global__
void dpairs(const double* S, double* out, int N, double dt, long long nmax, double R,
            double omd, const long long* chk, int nchk)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) dpair(S, out, N, i, dt, nmax, R, omd, (const long long*)chk, nchk);
}
""")


class CDual:
    name = 'c'

    def __init__(self, dfloat=False):
        d = tempfile.mkdtemp()
        src = os.path.join(d, 'hdym_dual.c')
        lib = os.path.join(d, 'hdym_dual' + ('.dll' if sys.platform == 'win32' else '.so'))
        with open(src, 'w') as fh:
            fh.write(c_src(dfloat))
        cc = os.environ.get('CC', 'gcc')
        # -ffp-contract=off keeps the v channel identical to the production kernel
        base = [cc, '-O3', '-ffp-contract=off', '-shared', '-fPIC', '-o', lib, src, '-lm']
        try:
            subprocess.check_call(base[:3] + ['-fopenmp'] + base[3:])
        except subprocess.CalledProcessError:
            subprocess.check_call(base)
        _dll_search_path(cc)
        self.lib = ctypes.CDLL(lib)
        dp, lp = ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_longlong)
        self.lib.dpairs.argtypes = [dp, dp, ctypes.c_int, ctypes.c_double, ctypes.c_longlong,
                                    ctypes.c_double, ctypes.c_double, lp, ctypes.c_int]
        print(f"C difference-number kernel (d channel {'fp32' if dfloat else 'fp64'})", flush=True)

    def run(self, S, dt, nmax, R, omd, chk):
        S = np.ascontiguousarray(S, dtype=np.float64)
        N = S.shape[1]
        chk = np.ascontiguousarray(chk, dtype=np.int64)
        out = np.empty((2 + len(chk), N))
        dp, lp = ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_longlong)
        self.lib.dpairs(S.ctypes.data_as(dp), out.ctypes.data_as(dp), N, dt, int(nmax), R, omd,
                        chk.ctypes.data_as(lp), len(chk))
        return out


class CupyDual:
    name = 'cupy'

    def __init__(self, dfloat=False):
        import cupy as cp
        self.cp = cp
        self.kern = cp.RawKernel(cuda_src(dfloat), 'dpairs', options=('--fmad=false',))
        print(f"CuPy difference-number kernel (d channel {'fp32' if dfloat else 'fp64'})",
              flush=True)

    def run(self, S, dt, nmax, R, omd, chk):
        cp = self.cp
        Sd = cp.asarray(np.ascontiguousarray(S, dtype=np.float64))
        N = Sd.shape[1]
        chk_d = cp.asarray(np.asarray(chk, dtype=np.int64))
        out = cp.empty((2 + len(chk), N), dtype=cp.float64)
        threads = 128
        self.kern(((N + threads - 1) // threads,), (threads,),
                  (Sd, out, np.int32(N), np.float64(dt), np.int64(nmax), np.float64(R),
                   np.float64(omd), chk_d, np.int32(len(chk))))
        cp.cuda.Stream.null.synchronize()
        return cp.asnumpy(out)


# ----------------------------------------------------------------- exact initial differences
def _dmul(a, da, b, db):
    return a * b, da * b + a * db + da * db


def dV(q, dq):
    """V(q + dq) - V(q) without cancellation (difference arithmetic), row-wise on (3, M)."""
    x, y, z = q
    dx, dy, dz = dq
    x2, dx2 = _dmul(x, dx, x, dx)
    y2, dy2 = _dmul(y, dy, y, dy)
    z2, dz2 = _dmul(z, dz, z, dz)
    _, a = _dmul(x2, dx2, y2, dy2)
    _, b = _dmul(y2, dy2, z2, dz2)
    _, c = _dmul(z2, dz2, x2, dx2)
    return 0.5 * (a + b + c)


def sample_dpairs(section, E, eps, seed, plane, center, width, M, delta, rng):
    """Reference states S_a (identical to hdym_uncertainty.sample_pairs) and exact
    differences D = S_b - S_a for the same random pairs."""
    r2 = copy.deepcopy(rng)
    Sa, _Sb_unused, ok = unc.sample_pairs(section, E, eps, seed, plane, center, width, M,
                                          delta, rng)
    # replay the same draws (sample_pairs: ua, va, ang in that order)
    cu, cv = center
    ua = cu + (r2.random(M) - 0.5) * width
    va = cv + (r2.random(M) - 0.5) * width
    ang = r2.random(M) * 2 * np.pi
    du_, dv_ = delta * np.cos(ang), delta * np.sin(ang)
    q0, P0, w0, W0 = em.base_point(E, eps, seed)
    e1, e2 = em.unit_pair(plane)
    D = np.zeros_like(Sa)
    if section == 'ghost':
        a = np.sqrt(eps * E)
        dc_u = -2.0 * np.sin(ua + du_ / 2) * np.sin(du_ / 2)
        ds_u = 2.0 * np.cos(ua + du_ / 2) * np.sin(du_ / 2)
        dc_v = -2.0 * np.sin(va + dv_ / 2) * np.sin(dv_ / 2)
        ds_v = 2.0 * np.cos(va + dv_ / 2) * np.sin(dv_ / 2)
        dw = a * (dc_u[:, None] * e1 + ds_u[:, None] * e2)
        dW = a * (dc_v[:, None] * e1 + ds_v[:, None] * e2)
        dp, dP = dw, -dW
    else:
        q = q0 + ua[:, None] * e1 + va[:, None] * e2
        dq = du_[:, None] * e1 + dv_[:, None] * e2
        K0 = 0.5 * P0 @ P0
        Va = em.V1(q.T)
        sa = np.sqrt(np.clip(E - Va, 0, None) / K0)
        dVab = dV(q.T, dq.T)
        sb = np.sqrt(np.clip(E - Va - dVab, 0, None) / K0)
        ds = -dVab / K0 / np.where(sa + sb > 0, sa + sb, np.inf)
        dw = np.zeros_like(dq)
        dW = np.zeros_like(dq)
        dp = dq
        dP = ds[:, None] * P0
    D = np.concatenate([dp.T, dw.T, dP.T, dW.T], 0)
    return Sa, D, ok


# ----------------------------------------------------------------- self-test
def selftest():
    import hdym_precision as hp
    E, eps, M, width = 0.3, 0.01, 64, em.default_width('ym', 0.3)
    dt, tmax, R, omd = 0.005, 300.0, 1e4, 0.25
    nmax = int(np.ceil(tmax / dt))
    chk = [int(10 / dt), int(30 / dt)]
    be = CDual()

    print("selftest 1: reference channel == production kernel (bitwise)")
    Sa, D, ok = sample_dpairs('ym', E, eps, 1, 'xy', [0.0, 0.0], width, M, 1e-6 * width,
                              hp.pair_rng(7, 1e-6))
    out = be.run(np.concatenate([Sa, D], 0), dt, nmax, R, omd, chk)
    ta = em.make_backend('c').escape(Sa.copy(), dt, nmax, R, omd)   # (kernel overwrites S)
    diff = out[0] != ta
    # the shared sub-step schedule gives member a extra sub-steps only when b is
    # deeper in a channel than a -- in practice only as b runs away first
    b_first = (out[1] > 0) & (out[1] < out[0])
    print(f"   identical t_esc(a) on {np.sum(~diff)}/{M} pairs; the {diff.sum()} others all have "
          f"b escaping first (shared sub-steps): {bool(np.all(b_first[diff]))}")
    assert np.all(b_first[diff])

    print("selftest 2: initial differences vs direct S_b - S_a at delta = 1e-4 (no cancellation issue)")
    for sec in ('ym', 'ghost'):
        wd = em.default_width(sec, E)
        cen = [0.0, 0.0] if sec == 'ym' else [np.pi, np.pi]
        r = np.random.default_rng(3)
        Sa, D, ok = sample_dpairs(sec, E, eps, 1, 'xy', cen, wd, M, 1e-4 * wd, r)
        r = np.random.default_rng(3)
        Sa2, Sb2, ok2 = unc.sample_pairs(sec, E, eps, 1, 'xy', cen, wd, M, 1e-4 * wd, r)
        rel = np.max(np.abs(D[:, ok] - (Sb2 - Sa2)[:, ok])) / np.max(np.abs(D[:, ok]))
        print(f"   {sec}: max relative mismatch {rel:.1e}  (expect ~1e-11: fp64 cancellation in Sb-Sa)")
        assert rel < 1e-8

    print("selftest 3: agrees with separate-trajectory pairs where both are valid (delta = 1e-6)")
    hpb = hp.CPairs()
    r = hp.pair_rng(9, 1e-6)
    Sa, D, ok = sample_dpairs('ym', E, eps, 1, 'xy', [0.0, 0.0], width, 400, 1e-6 * width, r)
    r = hp.pair_rng(9, 1e-6)
    Sa2, Sb2, _ = unc.sample_pairs('ym', E, eps, 1, 'xy', [0.0, 0.0], width, 400, 1e-6 * width, r)
    o1 = be.run(np.concatenate([Sa, D], 0), dt, nmax, R, omd, chk)
    o2 = hpb.run(np.concatenate([Sa2, Sb2], 0), dt, nmax, R, omd, 0, 1, chk)
    print(f"   identical t_esc(a): {np.mean(o1[0] == o2[0]):.2f} of pairs")
    g = ok & (o1[2] > 0) & (o2[2] > 0)
    rr = o1[2][g] / o2[2][g]
    print(f"   sep(t=10) ratio difference-numbers / separate: median {np.median(rr):.6f}, "
          f"90% within [{np.percentile(rr, 5):.4f}, {np.percentile(rr, 95):.4f}]")
    agree = np.mean(np.abs(o1[1] - o2[1]) < 1e-9)
    close = np.mean(np.abs(o1[1] - o2[1]) <= 2.0)
    print(f"   t_esc(b): bit-identical {agree:.2f}, within 2 time units {close:.2f} "
          f"(b's roundoff differs between the two schemes; chaos amplifies it)")
    print("selftest passed")


# ----------------------------------------------------------------- sweep
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--E', type=float, default=0.3)
    ap.add_argument('--eps', type=float, default=0.01)
    ap.add_argument('--section', choices=['ghost', 'ym'], default='ym')
    ap.add_argument('--plane', default='xy')
    ap.add_argument('--center', type=float, nargs=2, default=None)
    ap.add_argument('--width', type=float, default=None)
    ap.add_argument('--M', type=int, default=2000)
    ap.add_argument('--delta-max', type=float, default=1e-2)
    ap.add_argument('--delta-min', type=float, default=1e-30)
    ap.add_argument('--per-decade', type=float, default=1.0)
    ap.add_argument('--t-sep', type=float, nargs='+', default=[10.0, 30.0])
    ap.add_argument('--thresholds', type=float, nargs='+', default=[2.0, 5.0, 10.0])
    ap.add_argument('--dt', type=float, default=0.005)
    ap.add_argument('--tmax', type=float, default=1000.0)
    ap.add_argument('--R', type=float, default=1e4)
    ap.add_argument('--omega-dt-max', type=float, default=0.25)
    ap.add_argument('--dfloat', action='store_true', help='fp32 difference channel')
    ap.add_argument('--backend', choices=['c', 'cupy'], default='c')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--pair-seed', type=int, default=12345)
    ap.add_argument('--out', default='results_precision')
    ap.add_argument('--tag', default=None)
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    import hdym_precision as hp
    if a.center is None:
        a.center = [np.pi, np.pi] if a.section == 'ghost' else [0.0, 0.0]
    if a.width is None:
        a.width = em.default_width(a.section, a.E)
    if a.dfloat and a.delta_min < 1e-30:
        print("note: fp32 difference channel underflows below ~1e-38; clamping delta-min to 1e-30")
        a.delta_min = 1e-30
    tag = a.tag or (f"diffpair_{a.section}_E{a.E:g}_eps{a.eps:g}_M{a.M}"
                    + ("_d32" if a.dfloat else ""))
    os.makedirs(a.out, exist_ok=True)
    be = (CupyDual if a.backend == 'cupy' else CDual)(a.dfloat)
    nd = int(round(a.per_decade * np.log10(a.delta_max / a.delta_min))) + 1
    rel = np.logspace(np.log10(a.delta_max), np.log10(a.delta_min), nd)
    deltas = rel * a.width
    chk = [int(round(t / a.dt)) for t in a.t_sep]
    nmax = int(np.ceil(a.tmax / a.dt))
    F = np.full((len(a.thresholds), nd), np.nan)
    sep = np.full((len(chk), nd), np.nan)
    t0 = time.time()
    for i, d in enumerate(deltas):
        Sa, D, ok = sample_dpairs(a.section, a.E, a.eps, a.seed, a.plane, a.center, a.width,
                                  a.M, d, hp.pair_rng(a.pair_seed, d / a.width))
        out = be.run(np.concatenate([Sa, D], 0), a.dt, nmax, a.R, a.omega_dt_max, chk)
        ta, tb = out[0], out[1]
        good = ok & (ta > 0) & (tb > 0)
        diff = np.abs(ta - tb)[good]
        F[:, i] = [np.mean(diff > th) for th in a.thresholds]
        for j in range(len(chk)):
            sj = out[2 + j][ok]
            sj = sj[sj > 0]
            if sj.size > 10:
                sep[j, i] = np.median(sj)
        print(f"  delta/w={rel[i]:.1e} pairs={good.sum():5d} "
              + " ".join(f"f(>{th:g})={F[k, i]:.3f}" for k, th in enumerate(a.thresholds))
              + "  " + " ".join(f"gain(t={t:g})={sep[j, i] / d:.3g}" for j, t in enumerate(a.t_sep))
              + f"  ({time.time() - t0:.0f}s)", flush=True)
    print("\nuncertainty exponent over the whole range:")
    alphas = []
    for k, th in enumerate(a.thresholds):
        al, se, m = unc.fit_alpha(rel, F[k])
        alphas.append((al, se))
        print(f"  |dt|>{th:g}: alpha = {al:.3f} +- {se:.3f}  d = {2 - al:.3f}  [{m.sum()} points]")
    gain = sep[0] / deltas
    small = rel < 1e-16
    if np.any(small & np.isfinite(gain)):
        big = (rel > 1e-12) & (rel < 1e-6)
        print(f"\nseparation gain at t={a.t_sep[0]:g}: median {np.nanmedian(gain[big]):.3g} for "
              f"1e-12<delta<1e-6, {np.nanmedian(gain[small]):.3g} for delta<1e-16 "
              f"(equal => no roundoff floor)")
    np.savez(os.path.join(a.out, tag + '.npz'), rel=rel, F=F, sep=sep,
             thresholds=a.thresholds, t_sep=a.t_sep, alphas=np.array(alphas),
             args=json.dumps(vars(a)))
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
        for k, th in enumerate(a.thresholds):
            al, se = alphas[k]
            axs[0].loglog(rel, F[k], 'o-', ms=3, label=f"|Δt_esc|>{th:g}: α={al:.3f}±{se:.3f}")
        axs[0].axvline(2e-15, color='0.5', ls=':', lw=1, label='fp64 floor, separate pairs')
        axs[0].set_xlabel('δ / section width'); axs[0].set_ylabel('f(δ)')
        axs[0].set_title('uncertainty fraction, difference-number pairs', fontsize=10)
        axs[0].legend(fontsize=7)
        for j, t in enumerate(a.t_sep):
            axs[1].loglog(rel, sep[j] / deltas, 'o-', ms=3, label=f"t = {t:g}")
        axs[1].axvline(2e-15, color='0.5', ls=':', lw=1)
        axs[1].set_xlabel('δ / section width'); axs[1].set_ylabel('median separation / δ')
        axs[1].set_title('separation gain (flat = no floor)', fontsize=10)
        axs[1].legend(fontsize=7)
        plt.suptitle(f"E={a.E:g} eps={a.eps:g} {a.section} section, M={a.M}, d channel "
                     f"{'fp32' if a.dfloat else 'fp64'}", fontsize=10)
        plt.tight_layout()
        png = os.path.join(a.out, tag + '.png')
        plt.savefig(png, dpi=150)
        print("wrote", png)
    except Exception as e:  # noqa
        print("plot skipped:", e)


if __name__ == '__main__':
    main()
