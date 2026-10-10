#!/usr/bin/env python3
"""
hdym_precision.py -- mantissa-precision sweep of the escape-time landscape.

Question: is the fractal structure of t_esc at delta ~ 1e-14 of the section
width (Fig. 3c) a property of the flow, or roundoff in fractal clothing?

Method: rerun the uncertainty-exponent experiment of hdym_uncertainty.py with
the arithmetic emulated at p stored mantissa bits (round-to-nearest-even after
EVERY floating-point operation, fp64 exponent range kept, so p is the only
knob).  p = 52 is ordinary double; 23 has fp32's mantissa, 10 fp16's, 7
bf16's.  Exponent range is deliberately NOT emulated: fp8/fp16 would overflow
at the escape radius R = 1e4 and you would be measuring exponent limits.

Two observables per (p, delta), on the same pairs for every p (common random
numbers: the pair set at a given delta does not depend on p):

  f(delta)   fraction of pairs whose escape times differ by > threshold --
             the paper's observable.  With alpha ~ 0.05 a noise floor shows up
             only as a weak flattening, so f alone cannot locate it well.

  G(delta)   median phase-space separation of the pair at fixed early times
             t1 (default 10 and 30), before escape and before saturation.
             Above the noise floor  sep(t1) = g * delta;  below it sep(t1)
             stops shrinking.  Fitting  sep = g * sqrt(delta^2 + delta_f^2)
             gives the floor delta_f directly, in the same section units as
             delta, with no dependence on alpha.  This is the sharp test.

Expected (random-walk roundoff injection, amplified at the YM rate):
      delta_f(p)  ~  C * 2^-(p+1) / sqrt(2 lambda_YM dt)
i.e. log10 delta_f falls by log10(2) = 0.301 per bit.  The summary fits that
line over the bits where the floor is resolved, extrapolates to p = 52, and
compares with the paper's lowest point (1e-14 of the section width).

Schedule.  The production kernel sub-steps each trajectory by its own channel
depth (n = ceil(max|q| dt / omega_dt_max)).  Two members of a pair can then
take different numbers of sub-steps on the same step: a jump of order the
truncation error, far bigger than delta, that would fake a floor which does
NOT move with p.  --schedule shared (default) gives both members the larger
n; --schedule independent reproduces the production kernel (at p = 52 it is
bit-for-bit hdym_escape_map's C kernel; --selftest checks this).  Running
both is itself a diagnostic.

Examples
  python hdym_precision.py --selftest
  python hdym_precision.py --E 0.3 --eps 0.01 --section ym --M 2000 \
         --bits 52 40 32 28 23 20 16 12 10 8 6 --backend cupy
  python hdym_precision.py --E 0.3 --eps 0.01 --section ym --M 2000 \
         --bits 52 23 16 --schedule independent --backend cupy --tag indep
  python hdym_precision.py --replot results_precision/prec_ym_E0.3_eps0.01_M2000_shared.npz

Notes
  * FMA contraction is disabled (gcc -ffp-contract=off, NVRTC --fmad=false):
    a fused multiply-add rounds once where the emulation assumes twice.
  * On the RTX 4080 the kernel is FP64 throughout (the rounding is done on
    fp64 bit patterns), so cost is ~the production kernel times ~3-5 for the
    rounding.  Frozen low-p trajectories run to tmax: keep tmax modest.
  * Constants (Yoshida weights, dt, initial state) are also rounded to p bits.
"""
import argparse
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

# ----------------------------------------------------------------- kernel source
PAIR_SRC = r"""
#define YC0 0.6756035959798288
#define YC1 -0.17560359597982886
#define YC2 -0.17560359597982886
#define YC3 0.6756035959798288
#define YD0 1.3512071919596576
#define YD1 -1.7024143839193153
#define YD2 1.3512071919596576
#define MAXSUB 16384

/* round x to (52 - d) stored mantissa bits, round-to-nearest-even.
   d = 0 is plain fp64.  Carries into the exponent are handled by the integer
   add; +-inf is preserved; NaN may become inf (caught by the escape test). */
DEVFN double chop(double x, int d)
{
    if (d <= 0) return x;
    U64 u = AS_U64(x);
    U64 half = ((U64)1) << (d - 1);
    U64 lsb = (u >> d) & (U64)1;
    u += half - (U64)1 + lsb;
    u &= ~((((U64)1) << d) - (U64)1);
    return AS_DBL(u);
}
/* fixed point: round to the nearest multiple of 1/sc (sc = 2^F exactly, so the
   scaling is exact), ties to even.  sc <= 0 selects floating point (chop).   */
DEVFN double rnd(double x, int d, double sc)
{
    if (sc > 0.0) return rint(x * sc) / sc;
    return chop(x, d);
}
#define RA(a, b) rnd((a) + (b), d, sc)
#define RS(a, b) rnd((a) - (b), d, sc)
#define RM(a, b) rnd((a) * (b), d, sc)

/* same algebra and operation order as hdym_escape_map.STEP_SRC, every op rounded */
DEVFN void kick_r(const double* p, const double* w, double* P, double* W, double h, int d, double sc)
{
    double x = RS(p[0], w[0]), y = RS(p[1], w[1]), z = RS(p[2], w[2]);
    double x2 = RM(x, x), y2 = RM(y, y), z2 = RM(z, z);
    double yz = RA(y2, z2), zx = RA(z2, x2), xy = RA(x2, y2);
    double gx = RM(x, yz), gy = RM(y, zx), gz = RM(z, xy);
    double wx = w[0], wy = w[1], wz = w[2];
    double txy = RM(RM(2.0, x), y), txz = RM(RM(2.0, x), z), tyz = RM(RM(2.0, y), z);
    double hx = RA(RA(RM(yz, wx), RM(txy, wy)), RM(txz, wz));
    double hy = RA(RA(RM(txy, wx), RM(zx, wy)), RM(tyz, wz));
    double hz = RA(RA(RM(txz, wx), RM(tyz, wy)), RM(xy, wz));
    P[0] = RS(P[0], RM(h, RA(gx, hx)));
    P[1] = RS(P[1], RM(h, RA(gy, hy)));
    P[2] = RS(P[2], RM(h, RA(gz, hz)));
    W[0] = RA(W[0], RM(h, RA(hx, wx)));
    W[1] = RA(W[1], RM(h, RA(hy, wy)));
    W[2] = RA(W[2], RM(h, RA(hz, wz)));
}

DEVFN void drift_r(double* p, double* w, const double* P, const double* W, double ch, int d, double sc)
{
    int k;
    for (k = 0; k < 3; k++) { p[k] = RA(p[k], RM(ch, P[k])); w[k] = RS(w[k], RM(ch, W[k])); }
}

DEVFN void yoshida_r(double* p, double* w, double* P, double* W, double h, int d, double sc)
{
    drift_r(p, w, P, W, RM(rnd(YC0, d, sc), h), d, sc);
    kick_r(p, w, P, W, RM(rnd(YD0, d, sc), h), d, sc);
    drift_r(p, w, P, W, RM(rnd(YC1, d, sc), h), d, sc);
    kick_r(p, w, P, W, RM(rnd(YD1, d, sc), h), d, sc);
    drift_r(p, w, P, W, RM(rnd(YC2, d, sc), h), d, sc);
    kick_r(p, w, P, W, RM(rnd(YD2, d, sc), h), d, sc);
    drift_r(p, w, P, W, RM(rnd(YC3, d, sc), h), d, sc);
}

DEVFN int nsub_of(const double* p, const double* w, double dt, double omd)
{
    double om = 0.0, a;
    int k, n;
    for (k = 0; k < 3; k++) { a = fabs(p[k] - w[k]); if (a > om) om = a; }
    n = (int)ceil(om * dt / omd);
    if (n < 1) n = 1;
    if (n > MAXSUB) n = MAXSUB;
    return n;
}

DEVFN void step_n(double* p, double* w, double* P, double* W, double dt, int n, int d, double sc)
{
    int s;
    double h = rnd(dt / n, d, sc);
    for (s = 0; s < n; s++) yoshida_r(p, w, P, W, h, d, sc);
}

DEVFN int escaped(const double* p, const double* w, const double* P, const double* W, double R)
{
    double m = 0.0, a;
    int k;
    for (k = 0; k < 3; k++) {
        a = fabs(p[k]); if (a > m) m = a;
        a = fabs(w[k]); if (a > m) m = a;
        a = fabs(P[k]); if (a > m) m = a;
        a = fabs(W[k]); if (a > m) m = a;
    }
    return (!(m == m) || m > 1e300 || m > R);
}

/* Pair i: member a = S rows 0..11, member b = S rows 12..23 (layout (24, N)).
   out (2 + nchk, N): row 0 t_esc(a), row 1 t_esc(b)  (-1 = not escaped by nmax),
   rows 2.. phase-space distance |a - b| at step chk[j]  (-1 = a member had escaped). */
DEVFN void pair_traj(const double* S, double* out, int N, int i, double dt, long long nmax,
                     double R, double omd, int d, double sc, int shared, const long long* chk, int nchk)
{
    double pa[3], wa[3], Pa[3], Wa[3], pb[3], wb[3], Pb[3], Wb[3];
    int k, j = 0, na, nb, alive_a = 1, alive_b = 1;
    long long s;
    for (k = 0; k < 3; k++) {
        pa[k] = rnd(S[(0 + k) * N + i], d, sc);  wa[k] = rnd(S[(3 + k) * N + i], d, sc);
        Pa[k] = rnd(S[(6 + k) * N + i], d, sc);  Wa[k] = rnd(S[(9 + k) * N + i], d, sc);
        pb[k] = rnd(S[(12 + k) * N + i], d, sc); wb[k] = rnd(S[(15 + k) * N + i], d, sc);
        Pb[k] = rnd(S[(18 + k) * N + i], d, sc); Wb[k] = rnd(S[(21 + k) * N + i], d, sc);
    }
    out[0 * N + i] = -1.0;
    out[1 * N + i] = -1.0;
    for (k = 0; k < nchk; k++) out[(2 + k) * N + i] = -1.0;
    for (s = 1; s <= nmax && (alive_a || alive_b); s++) {
        na = alive_a ? nsub_of(pa, wa, dt, omd) : 0;
        nb = alive_b ? nsub_of(pb, wb, dt, omd) : 0;
        if (shared && alive_a && alive_b) { if (nb > na) na = nb; nb = na; }
        if (alive_a) {
            step_n(pa, wa, Pa, Wa, dt, na, d, sc);
            if (escaped(pa, wa, Pa, Wa, R)) { alive_a = 0; out[0 * N + i] = s * dt; }
        }
        if (alive_b) {
            step_n(pb, wb, Pb, Wb, dt, nb, d, sc);
            if (escaped(pb, wb, Pb, Wb, R)) { alive_b = 0; out[1 * N + i] = s * dt; }
        }
        while (j < nchk && chk[j] <= s) {
            if (chk[j] == s && alive_a && alive_b) {
                double acc = 0.0, e;
                for (k = 0; k < 3; k++) {
                    e = pa[k] - pb[k]; acc += e * e;  e = wa[k] - wb[k]; acc += e * e;
                    e = Pa[k] - Pb[k]; acc += e * e;  e = Wa[k] - Wb[k]; acc += e * e;
                }
                out[(2 + j) * N + i] = sqrt(acc);
            }
            j++;
        }
    }
}
"""

C_SRC = (r"""#include <math.h>
#include <string.h>
#include <stdint.h>
#define DEVFN static inline
typedef uint64_t U64;
static inline U64 AS_U64(double x) { U64 u; memcpy(&u, &x, 8); return u; }
static inline double AS_DBL(U64 u) { double x; memcpy(&x, &u, 8); return x; }
""" + PAIR_SRC + r"""
void pairs(const double* S, double* out, int N, double dt, long long nmax, double R, double omd,
           int d, double sc, int shared, const long long* chk, int nchk)
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
    for (i = 0; i < N; i++)
        pair_traj(S, out, N, i, dt, nmax, R, omd, d, sc, shared, chk, nchk);
}
""")

CUDA_SRC = (r"""#define DEVFN __device__ __forceinline__
typedef unsigned long long U64;
#define AS_U64(x) ((U64)__double_as_longlong(x))
#define AS_DBL(u) __longlong_as_double((long long)(u))
""" + PAIR_SRC + r"""
extern "C" __global__
void pairs(const double* S, double* out, int N, double dt, long long nmax, double R, double omd,
           int d, double sc, int shared, const long long* chk, int nchk)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) pair_traj(S, out, N, i, dt, nmax, R, omd, d, sc, shared, (const long long*)chk, nchk);
}
""")


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


class CPairs:
    name = 'c'

    def __init__(self):
        d = tempfile.mkdtemp()
        src = os.path.join(d, 'hdym_pairs.c')
        lib = os.path.join(d, 'hdym_pairs' + ('.dll' if sys.platform == 'win32' else '.so'))
        with open(src, 'w') as fh:
            fh.write(C_SRC)
        cc = os.environ.get('CC', 'gcc')
        base = [cc, '-O3', '-ffp-contract=off', '-shared', '-fPIC', '-o', lib, src, '-lm']
        try:
            subprocess.check_call(base[:3] + ['-fopenmp'] + base[3:])
            omp = 'with OpenMP'
        except subprocess.CalledProcessError:
            subprocess.check_call(base)
            omp = 'single-threaded'
        _dll_search_path(cc)
        self.lib = ctypes.CDLL(lib)
        dp = ctypes.POINTER(ctypes.c_double)
        lp = ctypes.POINTER(ctypes.c_longlong)
        self.lib.pairs.argtypes = [dp, dp, ctypes.c_int, ctypes.c_double, ctypes.c_longlong,
                                   ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_double,
                                   ctypes.c_int, lp, ctypes.c_int]
        print(f"C pair kernel compiled with {cc}, {omp}, -ffp-contract=off", flush=True)

    def run(self, S, dt, nmax, R, omd, drop, shared, chk, sc=0.0):
        S = np.ascontiguousarray(S, dtype=np.float64)
        N = S.shape[1]
        chk = np.ascontiguousarray(chk, dtype=np.int64)
        out = np.empty((2 + len(chk), N))
        dp = ctypes.POINTER(ctypes.c_double)
        lp = ctypes.POINTER(ctypes.c_longlong)
        self.lib.pairs(S.ctypes.data_as(dp), out.ctypes.data_as(dp), N, dt, int(nmax), R, omd,
                       int(drop), float(sc), int(shared), chk.ctypes.data_as(lp), len(chk))
        return out


class CupyPairs:
    name = 'cupy'

    def __init__(self):
        import cupy as cp
        self.cp = cp
        self.kern = cp.RawKernel(CUDA_SRC, 'pairs', options=('--fmad=false',))
        dev = cp.cuda.Device()
        free, total = dev.mem_info
        print(f"CuPy {cp.__version__} on device {dev.id}, {free / 1e9:.1f}/{total / 1e9:.1f} GB "
              f"free, NVRTC --fmad=false", flush=True)

    def run(self, S, dt, nmax, R, omd, drop, shared, chk, sc=0.0):
        cp = self.cp
        Sd = cp.asarray(np.ascontiguousarray(S, dtype=np.float64))
        N = Sd.shape[1]
        chk_d = cp.asarray(np.asarray(chk, dtype=np.int64))
        out = cp.empty((2 + len(chk), N), dtype=cp.float64)
        threads = 128
        blocks = (N + threads - 1) // threads
        self.kern((blocks,), (threads,),
                  (Sd, out, np.int32(N), np.float64(dt), np.int64(nmax), np.float64(R),
                   np.float64(omd), np.int32(drop), np.float64(sc), np.int32(shared), chk_d,
                   np.int32(len(chk))))
        cp.cuda.Stream.null.synchronize()
        return cp.asnumpy(out)


def make_backend(name):
    return {'c': CPairs, 'cupy': CupyPairs}[name]()


# ----------------------------------------------------------------- analysis helpers
def pair_rng(pair_seed, delta_rel):
    """Same pairs at the same delta for every p (common random numbers)."""
    key = int(round(-np.log10(delta_rel) * 1000))
    return np.random.default_rng([pair_seed, key])


def fit_floor(deltas, sep_med, sat=1e-3):
    """Fit median separation  sep = g * sqrt(delta^2 + delta_f^2)  in log space,
    using only unsaturated points (sep < sat).  Returns (g, delta_f, resolved).
    resolved = False means every used delta is >> delta_f; delta_f is then an
    upper bound (a third of the smallest delta)."""
    m = np.isfinite(sep_med) & (sep_med > 0) & (sep_med < sat)
    if m.sum() < 3:
        return np.nan, np.nan, False
    x, y = deltas[m], sep_med[m]
    order = np.argsort(x)
    x, y = x[order], y[order]
    # profile likelihood over delta_f on a log grid; g solved in closed form
    grid = np.logspace(np.log10(x[0]) - 3, np.log10(x[-1]), 600)
    best = (np.inf, np.nan, np.nan)
    for df in np.concatenate([[0.0], grid]):
        r = np.log(np.sqrt(x * x + df * df))
        lg = np.mean(np.log(y) - r)
        res = np.sum((np.log(y) - r - lg) ** 2)
        if res < best[0] - 1e-12:
            best = (res, np.exp(lg), df)
    _, g, df = best
    resolved = df > x[0] / 3.0
    if not resolved:
        df = x[0] / 3.0
    return g, df, bool(resolved)


def predicted_floor(bits, lam, dt):
    """Random-walk estimate: unit roundoff u = 2^-(p+1), injected each step and
    amplified at rate lam:  delta_f ~ u / sqrt(2 lam dt)  (coordinate units)."""
    return 2.0 ** -(bits + 1) / np.sqrt(2.0 * lam * dt)


# ----------------------------------------------------------------- self-test
def selftest():
    """At p = 52 with the independent schedule the pair kernel must reproduce
    hdym_escape_map's C kernel bit for bit; the shared schedule and p < 52 must
    run; chop must agree with numpy's float32 rounding at p = 23 (normal range)."""
    print("selftest 1: chop(p=23) vs float32 rounding")
    be = CPairs()
    rng = np.random.default_rng(0)
    x = rng.standard_normal(20000) * 10.0 ** rng.integers(-20, 20, 20000)
    # emulate chop in numpy for the comparison
    u = x.view(np.uint64).copy()
    dbits = 52 - 23
    half = np.uint64(1 << (dbits - 1))
    lsb = (u >> np.uint64(dbits)) & np.uint64(1)
    u = (u + half - np.uint64(1) + lsb) & ~np.uint64((1 << dbits) - 1)
    ok = np.array_equal(u.view(np.float64), x.astype(np.float32).astype(np.float64))
    print("   numpy chop == float32 rounding:", ok)
    assert ok

    print("selftest 2: p=52 independent schedule == hdym_escape_map C kernel")
    E, eps, M = 0.3, 0.01, 64
    width = em.default_width('ym', E)
    r = pair_rng(7, 1e-6)
    Sa, Sb, okp = unc.sample_pairs('ym', E, eps, 1, 'xy', [0.0, 0.0], width, M, 1e-6 * width, r)
    S = np.concatenate([Sa, Sb], 0)
    dt, tmax, R, omd = 0.005, 300.0, 1e4, 0.25
    nmax = int(np.ceil(tmax / dt))
    out = be.run(S, dt, nmax, R, omd, 0, 0, [int(10 / dt)])
    ref = em.make_backend('c')
    ta = ref.escape(Sa, dt, nmax, R, omd)
    tb = ref.escape(Sb, dt, nmax, R, omd)
    same = np.array_equal(out[0], ta) and np.array_equal(out[1], tb)
    print(f"   bitwise-identical escape times on {M} pairs: {same}  "
          f"(escaped {np.mean(ta > 0):.2f})")
    assert same, "pair kernel at p=52/independent does not match the production kernel"

    print("selftest 3: shared schedule and reduced precision run")
    for bits in (52, 23, 10):
        o = be.run(S, dt, nmax, R, omd, 52 - bits, 1, [int(10 / dt)])
        print(f"   p={bits:2d}: escaped {np.mean(o[0] > 0):.2f}, median t_esc "
              f"{np.median(o[0][o[0] > 0]) if np.any(o[0] > 0) else np.nan:.1f}, "
              f"median sep(t=10) {np.median(o[2][o[2] > 0]) if np.any(o[2] > 0) else np.nan:.2e}")
    print("selftest 4: fixed point keeps the state on the 2^-F grid")
    for F in (30, 12):
        o = be.run(S, dt, int(np.ceil(20.0 / dt)), R, omd, 0, 1, [int(10 / dt)], 2.0 ** F)
        print(f"   F={F}: escaped by t=20 {np.mean(o[0] > 0):.2f}, median sep(t=10) "
              f"{np.median(o[2][o[2] > 0]) if np.any(o[2] > 0) else np.nan:.2e}")
    x = rng.standard_normal(1000) * 3
    g = _round_np(x, 0, 2.0 ** 20)
    assert np.all(g * 2.0 ** 20 == np.rint(g * 2.0 ** 20))
    print("selftest passed")


# ----------------------------------------------------------------- main sweep
def run_sweep(a):
    be = make_backend(a.backend)
    width = a.width
    lam = a.lam_ym if a.lam_ym is not None else 0.38 * a.E ** 0.25
    q_scale = a.E ** 0.25
    chk_steps = [int(round(t / a.dt)) for t in a.t_sep]
    nmax = int(np.ceil(a.tmax / a.dt))
    thresholds = list(a.thresholds)
    results = []
    t0 = time.time()
    fixed = (a.mode == 'fixed')
    for bits in a.bits:
        drop = 0 if fixed else 52 - bits
        sc = 2.0 ** bits if fixed else 0.0
        # smallest delta: a few quanta (fixed) or ulps of the coordinates (float)
        quantum = 2.0 ** -bits * (1.0 if fixed else q_scale)
        dmin_abs = max(a.delta_min * width, a.ulps * quantum)
        dmax_abs = a.delta_max * width
        nd = max(3, int(round(a.per_decade * np.log10(dmax_abs / dmin_abs))) + 1)
        deltas = np.logspace(np.log10(dmax_abs), np.log10(dmin_abs), nd)
        F = np.full((len(thresholds), nd), np.nan)
        sepmed = np.full((len(chk_steps), nd), np.nan)
        npairs = np.zeros(nd, int)
        frozen = np.zeros(nd)
        identical = np.zeros(nd)
        medt = np.full(nd, np.nan)
        for i, d in enumerate(deltas):
            r = pair_rng(a.pair_seed, d / width)
            Sa, Sb, ok = unc.sample_pairs(a.section, a.E, a.eps, a.seed, a.plane, a.center,
                                          width, a.M, d, r)
            S = np.concatenate([Sa, Sb], 0)
            out = be.run(S, a.dt, nmax, a.R, a.omega_dt_max, drop,
                         1 if a.schedule == 'shared' else 0, chk_steps, sc)
            ta, tb = out[0], out[1]
            # pairs that the initial rounding collapsed onto the same state
            same0 = np.all(_round_np(Sa, drop, sc) == _round_np(Sb, drop, sc), axis=0)
            good = ok & (ta > 0) & (tb > 0) & ~same0
            frozen[i] = np.mean(ok & ((ta < 0) | (tb < 0)))
            identical[i] = np.mean(same0 & ok)
            npairs[i] = good.sum()
            if good.sum():
                diff = np.abs(ta - tb)[good]
                medt[i] = np.median(np.concatenate([ta[good], tb[good]]))
                F[:, i] = [np.mean(diff > th) for th in thresholds]
            for j in range(len(chk_steps)):
                sj = out[2 + j][ok & ~same0]
                sj = sj[sj > 0]
                if sj.size >= max(10, a.M // 10):
                    sepmed[j, i] = np.median(sj)
            print(f"  p={bits:2d} delta/w={d / width:.2e} pairs={good.sum():5d} "
                  f"f(>{thresholds[-1]:g})={F[-1, i]:.3f} "
                  + " ".join(f"sep(t={t:g})={sepmed[j, i]:.2e}" for j, t in enumerate(a.t_sep))
                  + f" frozen={frozen[i]:.2f} med t_esc={medt[i]:.0f} ({time.time() - t0:.0f}s)",
                  flush=True)
        floors = [fit_floor(deltas, sepmed[j]) for j in range(len(chk_steps))]
        # alpha above the floor (use the larger of the per-checkpoint floors)
        fv = [fl[1] for fl in floors if np.isfinite(fl[1])]
        df_fit = max(fv) if fv else np.nan
        above = deltas > 3.0 * df_fit if np.isfinite(df_fit) else np.ones(nd, bool)
        alphas = []
        for k in range(len(thresholds)):
            al, se, _ = unc.fit_alpha(deltas[above] / width, F[k, above])
            alphas.append((al, se))
        res = dict(bits=bits, deltas=deltas, F=F, sepmed=sepmed, npairs=npairs, frozen=frozen,
                   identical=identical, medt=medt, floors=floors, alphas=alphas,
                   pred=predicted_floor(bits, lam, a.dt))
        results.append(res)
        fl_txt = ", ".join(f"t={t:g}: {fl[1] / width:.2e}{'' if fl[2] else ' (upper bound)'}"
                           for t, fl in zip(a.t_sep, floors))
        al_txt = ", ".join(f">{th:g}: {al:.3f}+-{se:.3f}" for th, (al, se) in zip(thresholds, alphas))
        print(f"p={bits:2d}: floor/width [{fl_txt}]; predicted {res['pred'] / width:.2e}; "
              f"alpha above floor [{al_txt}]", flush=True)
    return results, width, lam


def _round_np(x, drop, sc=0.0):
    if sc > 0:
        return np.rint(np.asarray(x, dtype=np.float64) * sc) / sc
    if drop <= 0:
        return x
    u = np.ascontiguousarray(x, dtype=np.float64).view(np.uint64).copy()
    half = np.uint64(1 << (drop - 1))
    lsb = (u >> np.uint64(drop)) & np.uint64(1)
    u = (u + half - np.uint64(1) + lsb) & ~np.uint64((1 << drop) - 1)
    return u.view(np.float64)


def summarize(results, width, lam, a):
    """Fit log10(floor) = c - p*log10(2) over resolved floors; extrapolate to p=52."""
    summary = dict(E=a.E, eps=a.eps, section=a.section, mode=a.mode, schedule=a.schedule, dt=a.dt,
                   width=width, lam_ym=lam, t_sep=list(a.t_sep), per_bits=[])
    for r in results:
        summary['per_bits'].append(dict(
            bits=r['bits'],
            floor_over_width=[fl[1] / width for fl in r['floors']],
            floor_resolved=[fl[2] for fl in r['floors']],
            predicted_over_width=r['pred'] / width,
            alpha=[al for al, _ in r['alphas']],
            alpha_se=[se for _, se in r['alphas']],
            frozen_max=float(np.nanmax(r['frozen'])),
            median_t_esc=float(np.nanmedian(r['medt']))))
    fits = []
    for j, t in enumerate(a.t_sep):
        P, L = [], []
        for r in results:
            g, df, resolved = r['floors'][j]
            if resolved and r['bits'] >= a.fit_min_bits:
                P.append(r['bits'])
                L.append(np.log10(df / width))
        P, L = np.array(P, float), np.array(L)
        fit = dict(t_sep=t, n=len(P))
        if len(P) >= 2:
            c_fixed = np.mean(L + P * np.log10(2.0))
            fit.update(c_fixed=c_fixed,
                       extrap52_fixed=10 ** (c_fixed - 52 * np.log10(2.0)),
                       scatter_fixed=float(np.std(L - (c_fixed - P * np.log10(2.0)))))
            C_ratio = 10 ** np.mean([np.log10(10 ** l / (predicted_floor(p, lam, a.dt) / width))
                                     for p, l in zip(P, L)])
            fit['C_vs_prediction'] = C_ratio
        if len(P) >= 3:
            s, c = np.polyfit(P, L, 1)
            fit.update(slope_free=s, slope_expected=-np.log10(2.0),
                       extrap52_free=10 ** (c + 52 * s))
        fits.append(fit)
    summary['fits'] = fits
    print("\n=== summary ===")
    print(f"{'bits':>4} " + " ".join(f"{'floor/w t=' + format(t, 'g'):>18}" for t in a.t_sep)
          + f" {'predicted/w':>12} {'alpha(>' + format(a.thresholds[-1], 'g') + ')':>12} "
            f"{'frozen':>7} {'med t_esc':>9}")
    for pb in summary['per_bits']:
        fl = " ".join(f"{v:>16.2e}{' ' if res else '<'} " for v, res in
                      zip(pb['floor_over_width'], pb['floor_resolved']))
        print(f"{pb['bits']:>4} {fl} {pb['predicted_over_width']:>12.2e} "
              f"{pb['alpha'][-1]:>12.3f} {pb['frozen_max']:>7.2f} {pb['median_t_esc']:>9.1f}")
    print("('<' = not resolved: upper bound)")
    for f in fits:
        if f['n'] >= 2:
            print(f"t_sep={f['t_sep']:g}: {f['n']} resolved floors; fixed slope -log10(2)/bit "
                  f"-> p=52 floor/width ~ {f['extrap52_fixed']:.1e} (scatter "
                  f"{f['scatter_fixed']:.2f} dex); C vs random-walk prediction = "
                  f"{f['C_vs_prediction']:.2g}"
                  + (f"; free slope {f['slope_free']:.3f}/bit (expect -0.301), p=52 -> "
                     f"{f['extrap52_free']:.1e}" if 'slope_free' in f else ''))
        else:
            print(f"t_sep={f['t_sep']:g}: fewer than 2 resolved floors; no extrapolation")
    print("paper's lowest point: 1e-14 of the section width")
    return summary


def plot(results, width, summary, a, png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap('viridis')
    bits = [r['bits'] for r in results]
    col = {b: cmap(i / max(1, len(bits) - 1)) for i, b in enumerate(sorted(bits))}
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    ax = axs[0]
    for r in results:
        ax.loglog(r['deltas'] / width, r['F'][-1], 'o-', ms=3, color=col[r['bits']],
                  label=f"p={r['bits']}")
    ax.set_xlabel('delta / section width')
    ax.set_ylabel(f"f(delta), |dt_esc| > {a.thresholds[-1]:g}")
    ax.set_title('uncertainty fraction vs mantissa bits', fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    ax = axs[1]
    j = len(a.t_sep) - 1
    for r in results:
        g = r['sepmed'][j] / r['deltas']
        ax.loglog(r['deltas'] / width, g, 'o-', ms=3, color=col[r['bits']])
        gg, df, res = r['floors'][j]
        if np.isfinite(gg):
            xx = np.logspace(np.log10(r['deltas'].min()), np.log10(r['deltas'].max()), 100)
            ax.loglog(xx / width, gg * np.sqrt(xx ** 2 + df ** 2) / xx, '-', lw=0.8,
                      color=col[r['bits']], alpha=0.6)
    ax.set_xlabel('delta / section width')
    ax.set_ylabel(f"median sep(t={a.t_sep[j]:g}) / delta")
    ax.set_title('separation gain: flat above the floor, ~1/delta below', fontsize=10)
    ax = axs[2]
    for j, t in enumerate(a.t_sep):
        pb = summary['per_bits']
        P = np.array([x['bits'] for x in pb])
        Fl = np.array([x['floor_over_width'][j] for x in pb])
        R = np.array([x['floor_resolved'][j] for x in pb])
        mk = 'o' if j == len(a.t_sep) - 1 else 's'
        cj = f"C{j}"
        ax.semilogy(P[R], Fl[R], mk, color=cj, label=f"floor, t={t:g}")
        ax.semilogy(P[~R], Fl[~R], 'v', mfc='none', color='gray',
                    label='upper bound' if j == 0 else None)
        f = summary['fits'][j]
        if f['n'] >= 2:
            pp = np.arange(min(P.min(), a.fit_min_bits), 53)
            ax.semilogy(pp, 10 ** (f['c_fixed'] - pp * np.log10(2.0)), '--', lw=0.8,
                        color=cj, label=f"fit -0.301/bit, t={t:g}")
    pp = np.arange(min(bits), 53)
    ax.semilogy(pp, [predicted_floor(p, summary['lam_ym'], a.dt) / width for p in pp], ':',
                color='k', label='u/sqrt(2 lam dt)')
    ax.axhline(1e-14, color='C3', lw=0.8, label="paper's lowest delta (Fig. 3c)")
    ax.set_xlabel('fractional bits F (quantum 2^-F)' if getattr(a, 'mode', 'float') == 'fixed'
                  else 'stored mantissa bits p')
    ax.set_ylabel('floor delta_f / section width')
    ax.set_title('noise floor vs precision (slope -0.301/bit expected)', fontsize=10)
    ax.legend(fontsize=7)
    plt.suptitle(f"{getattr(a, 'mode', 'float')} point, E={a.E:g} eps={a.eps:g} {a.section} section, "
                 f"{a.schedule} schedule, "
                 f"dt={a.dt:g}, M={a.M}", fontsize=10)
    plt.tight_layout()
    plt.savefig(png, dpi=150)
    print("wrote", png)


def save(results, width, summary, a, path):
    arr = {}
    for r in results:
        b = r['bits']
        arr[f'p{b}_deltas'] = r['deltas']
        arr[f'p{b}_F'] = r['F']
        arr[f'p{b}_sepmed'] = r['sepmed']
        arr[f'p{b}_npairs'] = r['npairs']
        arr[f'p{b}_frozen'] = r['frozen']
        arr[f'p{b}_identical'] = r['identical']
        arr[f'p{b}_medt'] = r['medt']
    np.savez(path, bits=np.array([r['bits'] for r in results]), width=width,
             summary=json.dumps(summary, default=float), args=json.dumps(vars(a)), **arr)
    print("wrote", path)


def load(path):
    z = np.load(path, allow_pickle=True)
    a = argparse.Namespace(**json.loads(str(z['args'])))
    if not hasattr(a, 'mode'):
        a.mode = 'float'
    width = float(z['width'])
    lam = a.lam_ym if a.lam_ym is not None else 0.38 * a.E ** 0.25
    results = []
    for b in z['bits']:
        b = int(b)
        r = dict(bits=b, deltas=z[f'p{b}_deltas'], F=z[f'p{b}_F'], sepmed=z[f'p{b}_sepmed'],
                 npairs=z[f'p{b}_npairs'], frozen=z[f'p{b}_frozen'],
                 identical=z[f'p{b}_identical'], medt=z[f'p{b}_medt'],
                 pred=predicted_floor(b, lam, a.dt))
        r['floors'] = [fit_floor(r['deltas'], r['sepmed'][j]) for j in range(len(a.t_sep))]
        df = np.nanmax([fl[1] for fl in r['floors']])
        above = r['deltas'] > 3.0 * df if np.isfinite(df) else np.ones(len(r['deltas']), bool)
        r['alphas'] = [unc.fit_alpha(r['deltas'][above] / width, r['F'][k, above])[:2]
                       for k in range(len(a.thresholds))]
        results.append(r)
    return results, width, lam, a


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--replot', default=None, help='npz from an earlier run: refit and replot')
    ap.add_argument('--E', type=float, default=0.3)
    ap.add_argument('--eps', type=float, default=0.01)
    ap.add_argument('--section', choices=['ghost', 'ym'], default='ym')
    ap.add_argument('--plane', default='xy')
    ap.add_argument('--center', type=float, nargs=2, default=None)
    ap.add_argument('--width', type=float, default=None)
    ap.add_argument('--M', type=int, default=2000, help='pairs per delta')
    ap.add_argument('--bits', type=int, nargs='+',
                    default=[52, 44, 36, 32, 28, 23, 20, 16, 12, 10, 8, 6],
                    help='stored mantissa bits (52 = fp64, 23 = fp32 mantissa, 10 = fp16, 7 = bf16)')
    ap.add_argument('--delta-max', type=float, default=1e-3, help='largest delta / width')
    ap.add_argument('--delta-min', type=float, default=1e-17,
                    help='smallest delta / width (also clamped to --ulps ulps of the coordinates)')
    ap.add_argument('--ulps', type=float, default=4.0)
    ap.add_argument('--per-decade', type=float, default=2.0, help='delta points per decade')
    ap.add_argument('--t-sep', type=float, nargs='+', default=[10.0, 30.0],
                    help='times at which the pair separation is recorded')
    ap.add_argument('--thresholds', type=float, nargs='+', default=[5.0, 10.0],
                    help='|dt_esc| thresholds; the last one is plotted')
    ap.add_argument('--schedule', choices=['shared', 'independent'], default='shared')
    ap.add_argument('--mode', choices=['float', 'fixed'], default='float',
                    help='float: --bits = stored mantissa bits, fp64 exponent range. '
                         'fixed: --bits = fractional bits F, every result rounded to a '
                         'multiple of 2^-F (word length F + 1 + ceil(log2 R) bits)')
    ap.add_argument('--dt', type=float, default=0.005)
    ap.add_argument('--tmax', type=float, default=1000.0,
                    help='frozen (low-p) trajectories run to tmax; median t_esc at E=0.3 is ~95')
    ap.add_argument('--R', type=float, default=1e4)
    ap.add_argument('--omega-dt-max', type=float, default=0.25)
    ap.add_argument('--lam-ym', type=float, default=None,
                    help='lambda_YM for the predicted floor (default 0.38 E^1/4)')
    ap.add_argument('--fit-min-bits', type=int, default=16,
                    help='only floors at p >= this enter the extrapolation (low p is not YM)')
    ap.add_argument('--backend', choices=['cupy', 'c'], default='c')
    ap.add_argument('--seed', type=int, default=1, help='base-point seed (as in the maps)')
    ap.add_argument('--pair-seed', type=int, default=12345)
    ap.add_argument('--out', default='results_precision')
    ap.add_argument('--tag', default=None)
    a = ap.parse_args()

    if a.selftest:
        selftest()
        return
    if a.replot:
        results, width, lam, a2 = load(a.replot)
        summary = summarize(results, width, lam, a2)
        plot(results, width, summary, a2, a.replot.replace('.npz', '.png'))
        return

    if a.center is None:
        a.center = [np.pi, np.pi] if a.section == 'ghost' else [0.0, 0.0]
    if a.width is None:
        a.width = em.default_width(a.section, a.E)
    a.bits = sorted(set(a.bits), reverse=True)
    for b in a.bits:
        if not 1 <= b <= 52:
            raise SystemExit("--bits must be in 1..52")
    if a.mode == 'fixed' and max(a.bits) > 52 - int(np.ceil(np.log2(a.R))):
        print(f"note: fixed point with F > {52 - int(np.ceil(np.log2(a.R)))} is limited by the "
              f"fp64 carrier for |x| near R; such F are exact only while |x| < 2^(52-F)")
    os.makedirs(a.out, exist_ok=True)
    tag = a.tag or (f"prec_{a.section}_E{a.E:g}_eps{a.eps:g}_M{a.M}_{a.schedule}"
                    + ("_fixed" if a.mode == 'fixed' else ""))
    print(f"{tag}: bits {a.bits}, M={a.M} pairs/delta, dt={a.dt}, tmax={a.tmax}, "
          f"width={a.width:.3g}, schedule={a.schedule}", flush=True)
    results, width, lam = run_sweep(a)
    summary = summarize(results, width, lam, a)
    path = os.path.join(a.out, tag + '.npz')
    save(results, width, summary, a, path)
    with open(os.path.join(a.out, tag + '_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=1, default=float)
    try:
        plot(results, width, summary, a, path.replace('.npz', '.png'))
    except Exception as e:                    # noqa
        print("plot skipped:", e)


if __name__ == '__main__':
    main()
