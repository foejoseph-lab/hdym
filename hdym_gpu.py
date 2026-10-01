#!/usr/bin/env python3
"""
hdym_gpu.py (v2) -- large-ensemble transverse (linearised-ghost) Lyapunov
exponent with the ENTIRE integrator inside one CUDA kernel.

Why v2: v1 built the Yoshida step from ~170 small CuPy element-wise kernels.
Under Windows' WDDM driver each launch costs ~100 us, so v1 ran at ~29 ms/step
on an RTX 4080 SUPER (GPU '97% busy' while actually launch-starved).  Here one
thread integrates one trajectory for `renorm_every` steps per launch:
launches per step drop from ~170 to ~0.02.

Backends
  --backend cupy   : CUDA RawKernel (production)
  --backend c      : the SAME C source compiled with gcc via ctypes, OpenMP
                     over all cores (how the kernel was tested without a GPU;
                     needs gcc -- on Windows e.g. MinGW-w64 / MSYS2 -- and
                     `set CC=...` if it is not called gcc)
  --backend numpy  : pure-numpy reference (slow, no compiler needed)

Usage
  python hdym_gpu.py --E 0.2154 --N 4096 --target-abs 1e-4 --t-block 100   # check: lam ~ 7.0e-3
  python hdym_gpu.py --E 0.015 0.02 0.027 0.035 --N 200000 --target-abs 3e-7

Everything else (estimator, stopping rule, outputs) is as in v1 / hdym_core:
lambda from the slope of the per-trajectory log-growth over the trailing half
of the checkpoints; stop when standard error and drift are both below
max(target_abs, target_rel*|lambda|) for n_consec consecutive blocks.
State layout is (12, N): p[3], w[3], P[3], Pi[3] rows, so neighbouring
threads read neighbouring addresses.
"""
import argparse
import csv
import ctypes
import os
import sys
import tempfile
import time

import numpy as np

# ---------------------------------------------------------------- the integrator (C)
# Shared between CUDA and plain C.  One call advances trajectory i by `nsteps`
# steps of size dt (with per-trajectory channel sub-stepping), renormalising the
# ghost vector every `renorm_every` steps and (if accumulate) adding log|w| to
# logsum[i].  Yoshida-4 coefficients are baked in as literals.
STEP_SRC = r"""
#define YC0 0.6756035959798288
#define YC1 -0.17560359597982886
#define YC2 -0.17560359597982886
#define YC3 0.6756035959798288
#define YD0 1.3512071919596576
#define YD1 -1.7024143839193153
#define YD2 1.3512071919596576

DEVFN void kick(const double* p, const double* w, double* P, double* Pi, double h, double od)
{
    double x = p[0], y = p[1], z = p[2];
    double x2 = x * x, y2 = y * y, z2 = z * z;
    /* grad V */
    double gx = x * (y2 + z2), gy = y * (z2 + x2), gz = z * (x2 + y2);
    /* Hess V . w */
    double wx = w[0], wy = w[1], wz = w[2];
    /* od = 1: full Hessian (physical).  od = 0: off-diagonal entries dropped,
       i.e. three independent scalar oscillators driven by K_xx, K_yy, K_zz
       (the norm then grows at the largest of the three rates) -- a diagnostic
       for which entries of Hess V drive the ghost at a given energy. */
    double hx = (y2 + z2) * wx + od * (2.0 * x * y * wy + 2.0 * x * z * wz);
    double hy = od * (2.0 * x * y * wx + 2.0 * y * z * wz) + (z2 + x2) * wy;
    double hz = od * (2.0 * x * z * wx + 2.0 * y * z * wy) + (x2 + y2) * wz;
    P[0] -= h * gx;  P[1] -= h * gy;  P[2] -= h * gz;
    Pi[0] -= h * (wx + hx);  Pi[1] -= h * (wy + hy);  Pi[2] -= h * (wz + hz);
}

DEVFN void yoshida(double* p, double* w, double* P, double* Pi, double h, double od)
{
    int k;
    for (k = 0; k < 3; k++) { p[k] += YC0 * h * P[k]; w[k] += YC0 * h * Pi[k]; }
    kick(p, w, P, Pi, YD0 * h, od);
    for (k = 0; k < 3; k++) { p[k] += YC1 * h * P[k]; w[k] += YC1 * h * Pi[k]; }
    kick(p, w, P, Pi, YD1 * h, od);
    for (k = 0; k < 3; k++) { p[k] += YC2 * h * P[k]; w[k] += YC2 * h * Pi[k]; }
    kick(p, w, P, Pi, YD2 * h, od);
    for (k = 0; k < 3; k++) { p[k] += YC3 * h * P[k]; w[k] += YC3 * h * Pi[k]; }
}

DEVFN void step_traj(double* S, double* logsum, int N, int i, double dt, int nsteps,
                     int renorm_every, double omega_dt_max, int accumulate, double od)
{
    double p[3], w[3], P[3], Pi[3];
    int k, s, n, sub;
    for (k = 0; k < 3; k++) {
        p[k]  = S[(0 + k) * N + i];
        w[k]  = S[(3 + k) * N + i];
        P[k]  = S[(6 + k) * N + i];
        Pi[k] = S[(9 + k) * N + i];
    }
    double acc = 0.0;
    for (s = 1; s <= nsteps; s++) {
        double om = fabs(p[0]);
        if (fabs(p[1]) > om) om = fabs(p[1]);
        if (fabs(p[2]) > om) om = fabs(p[2]);
        n = (int)ceil(om * dt / omega_dt_max);
        if (n < 1) n = 1;
        if (n == 1) {
            yoshida(p, w, P, Pi, dt, od);
        } else {
            double h = dt / n;
            for (sub = 0; sub < n; sub++) yoshida(p, w, P, Pi, h, od);
        }
        if (s % renorm_every == 0) {
            if (od == 0.0) {
                /* three independent scalar oscillators: renormalise each
                   (w_k, Pi_k) pair separately and accumulate the MEAN of the
                   three log-growths = the single-mode exponent, not the max
                   of three finite-time exponents. */
                double lsum = 0.0;
                for (k = 0; k < 3; k++) {
                    double nk = sqrt(w[k]*w[k] + Pi[k]*Pi[k]);
                    lsum += log(nk);
                    w[k] /= nk;  Pi[k] /= nk;
                }
                if (accumulate) acc += lsum / 3.0;
            } else {
                double nrm = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]
                                + Pi[0]*Pi[0] + Pi[1]*Pi[1] + Pi[2]*Pi[2]);
                if (accumulate) acc += log(nrm);
                for (k = 0; k < 3; k++) { w[k] /= nrm; Pi[k] /= nrm; }
            }
        }
    }
    if (accumulate) logsum[i] += acc;
    for (k = 0; k < 3; k++) {
        S[(0 + k) * N + i] = p[k];
        S[(3 + k) * N + i] = w[k];
        S[(6 + k) * N + i] = P[k];
        S[(9 + k) * N + i] = Pi[k];
    }
}
"""

CUDA_SRC = ("#define DEVFN __device__ __forceinline__\n" + STEP_SRC + r"""
extern "C" __global__
void advance(double* S, double* logsum, int N, double dt, int nsteps,
             int renorm_every, double omega_dt_max, int accumulate, double od)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) step_traj(S, logsum, N, i, dt, nsteps, renorm_every, omega_dt_max, accumulate, od);
}
""")

C_SRC = ("#include <math.h>\n#define DEVFN static inline\n" + STEP_SRC + r"""
void advance(double* S, double* logsum, int N, double dt, int nsteps,
             int renorm_every, double omega_dt_max, int accumulate, double od)
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
    for (i = 0; i < N; i++)
        step_traj(S, logsum, N, i, dt, nsteps, renorm_every, omega_dt_max, accumulate, od);
}
""")


# ---------------------------------------------------------------- backends
class CupyBackend:
    name = 'cupy'

    def __init__(self):
        import cupy as cp
        self.cp = cp
        self.kern = cp.RawKernel(CUDA_SRC, 'advance')
        dev = cp.cuda.Device()
        free, total = dev.mem_info
        print(f"CuPy {cp.__version__} on device {dev.id}, {free / 1e9:.1f}/{total / 1e9:.1f} GB free",
              flush=True)

    def alloc(self, S_host):
        self.S = self.cp.asarray(S_host)
        self.logsum = self.cp.zeros(S_host.shape[1])
        self.N = S_host.shape[1]

    def advance(self, dt, nsteps, renorm_every, omega_dt_max, accumulate, od=1.0):
        N = self.N
        threads = 128
        blocks = (N + threads - 1) // threads
        self.kern((blocks,), (threads,),
                  (self.S, self.logsum, np.int32(N), np.float64(dt), np.int32(nsteps),
                   np.int32(renorm_every), np.float64(omega_dt_max), np.int32(accumulate),
                   np.float64(od)))

    def logsum_host(self):
        return self.cp.asnumpy(self.logsum)

    def S_host(self):
        return self.cp.asnumpy(self.S)


class CBackend:
    name = 'c'

    def __init__(self):
        import subprocess
        d = tempfile.mkdtemp()
        src = os.path.join(d, 'hdym_step.c')
        lib = os.path.join(d, 'hdym_step' + ('.dll' if sys.platform == 'win32' else '.so'))
        with open(src, 'w') as fh:
            fh.write(C_SRC)
        cc = os.environ.get('CC', 'gcc')
        try:
            subprocess.check_call([cc, '-O3', '-fopenmp', '-shared', '-fPIC', '-o', lib, src, '-lm'])
            omp = 'with OpenMP'
        except subprocess.CalledProcessError:
            subprocess.check_call([cc, '-O3', '-shared', '-fPIC', '-o', lib, src, '-lm'])
            omp = 'single-threaded (no OpenMP)'
        self.lib = ctypes.CDLL(lib)
        dp = ctypes.POINTER(ctypes.c_double)
        self.lib.advance.argtypes = [dp, dp, ctypes.c_int, ctypes.c_double, ctypes.c_int,
                                     ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_double]
        print(f"C backend compiled with {cc}, {omp}", flush=True)

    def alloc(self, S_host):
        self.S = np.ascontiguousarray(S_host)
        self.logsum = np.zeros(S_host.shape[1])
        self.N = S_host.shape[1]

    def advance(self, dt, nsteps, renorm_every, omega_dt_max, accumulate, od=1.0):
        dp = ctypes.POINTER(ctypes.c_double)
        self.lib.advance(self.S.ctypes.data_as(dp), self.logsum.ctypes.data_as(dp),
                         self.N, dt, nsteps, renorm_every, omega_dt_max, accumulate, od)

    def logsum_host(self):
        return self.logsum.copy()

    def S_host(self):
        return self.S.copy()


class NumpyBackend:
    """Reference implementation (vectorised numpy).  Slow; for checking."""
    name = 'numpy'
    YC = (0.6756035959798288, -0.17560359597982886, -0.17560359597982886, 0.6756035959798288)
    YD = (1.3512071919596576, -1.7024143839193153, 1.3512071919596576, 0.0)

    def __init__(self):
        print("numpy reference backend", flush=True)

    def alloc(self, S_host):
        self.S = S_host.copy()
        self.logsum = np.zeros(S_host.shape[1])
        self.N = S_host.shape[1]

    @staticmethod
    def _kick(p, w, P, Pi, h, od):
        x, y, z = p
        x2, y2, z2 = x * x, y * y, z * z
        g = np.array([x * (y2 + z2), y * (z2 + x2), z * (x2 + y2)])
        wx, wy, wz = w
        H = np.array([(y2 + z2) * wx + od * (2 * x * y * wy + 2 * x * z * wz),
                      od * (2 * x * y * wx + 2 * y * z * wz) + (z2 + x2) * wy,
                      od * (2 * x * z * wx + 2 * y * z * wy) + (x2 + y2) * wz])
        return P - h * g, Pi - h * (w + H)

    def _yoshida(self, p, w, P, Pi, h, od):
        for c, d in zip(self.YC, self.YD):
            p = p + c * h * P
            w = w + c * h * Pi
            if d != 0.0:
                P, Pi = self._kick(p, w, P, Pi, d * h, od)
        return p, w, P, Pi

    def advance(self, dt, nsteps, renorm_every, omega_dt_max, accumulate, od=1.0):
        S = self.S
        p, w, P, Pi = S[0:3], S[3:6], S[6:9], S[9:12]
        for s in range(1, nsteps + 1):
            om = np.max(np.abs(p), axis=0)
            nsub = np.maximum(1, np.ceil(om * dt / omega_dt_max)).astype(int)
            out = self._yoshida(p, w, P, Pi, dt, od)
            for n in np.unique(nsub):
                if n == 1:
                    continue
                m = nsub == n
                sub = (p[:, m], w[:, m], P[:, m], Pi[:, m])
                for _ in range(n):
                    sub = self._yoshida(*sub, dt / n, od)
                for a, b in zip(out, sub):
                    a[:, m] = b
            p, w, P, Pi = out
            if s % renorm_every == 0:
                nrm = np.sqrt(np.sum(w * w, 0) + np.sum(Pi * Pi, 0))
                if accumulate:
                    self.logsum += np.log(nrm)
                w = w / nrm
                Pi = Pi / nrm
        self.S = np.concatenate([p, w, P, Pi], 0)

    def logsum_host(self):
        return self.logsum.copy()

    def S_host(self):
        return self.S.copy()


# ---------------------------------------------------------------- physics helpers (host)
def V(q):
    x, y, z = q
    return 0.5 * (x * x * y * y + y * y * z * z + z * z * x * x)


def hessV_min_eig(q):
    x, y, z = q
    Hm = np.empty((q.shape[1], 3, 3))
    Hm[:, 0, 0] = y * y + z * z
    Hm[:, 1, 1] = z * z + x * x
    Hm[:, 2, 2] = x * x + y * y
    Hm[:, 0, 1] = Hm[:, 1, 0] = 2 * x * y
    Hm[:, 0, 2] = Hm[:, 2, 0] = 2 * x * z
    Hm[:, 1, 2] = Hm[:, 2, 1] = 2 * y * z
    return np.linalg.eigvalsh(Hm)[:, 0]


def initial_state(N, E, seed):
    rng = np.random.default_rng(seed)
    q = rng.standard_normal((3, N))
    P = rng.standard_normal((3, N))
    E0 = 0.5 * np.sum(P * P, 0) + V(q)
    lam = (E / E0) ** 0.25
    q *= lam
    P *= lam ** 2
    w = rng.standard_normal((3, N))
    Pi = rng.standard_normal((3, N))
    nrm = np.sqrt(np.sum(w * w, 0) + np.sum(Pi * Pi, 0))
    w /= nrm
    Pi /= nrm
    return np.concatenate([q, w, P, Pi], 0)


# ---------------------------------------------------------------- adaptive run
def run_energy(be, E, N, dt, seed, target_abs, target_rel, t_block, t_burn, tmax, wall,
               n_consec, renorm_every, omega_dt_max, log, od=1.0, t_burn_ym=None,
               trend_check=True):
    t0 = time.time()
    if t_burn_ym is not None:
        t_burn = t_burn_ym * E ** -0.25      # YM relaxation time scales as E^-1/4
    S0 = initial_state(N, E, seed)
    E0 = 0.5 * np.sum(S0[6:9] ** 2, 0) + V(S0[0:3])
    be.alloc(S0)
    nburn = int(np.ceil(t_burn / dt / renorm_every)) * renorm_every
    nblock = int(np.ceil(t_block / dt / renorm_every)) * renorm_every
    be.advance(dt, nburn, renorm_every, omega_dt_max, 0, od)
    k = 0
    T = 0.0
    chk_T, chk_S, hist = [], [], []
    tach_num, tach_den = 0.0, 0
    lam_prev, good, converged = None, 0, False
    lam_i = lam = se = np.nan
    lam_extrap = np.nan
    while True:
        be.advance(dt, nblock, renorm_every, omega_dt_max, 1, od)
        k += nblock
        T = k * dt
        S = be.logsum_host()
        if not np.all(np.isfinite(S)):
            raise FloatingPointError(f"non-finite log-growth at E={E}, T={T}")
        chk_T.append(T)
        chk_S.append(S)
        Sh = be.S_host()
        mn = hessV_min_eig(Sh[0:3, :min(N, 2048)])
        tach_num += float(np.mean(1.0 + mn <= 0.0))
        tach_den += 1
        Ta = np.array(chk_T)
        i0 = max(0, int(len(Ta) / 2) - 1)
        if len(Ta) - i0 >= 3:
            Ts = Ta[i0:] - Ta[i0:].mean()
            Ss = np.array(chk_S[i0:])
            lam_i = (Ts[:, None] * (Ss - Ss.mean(0))).sum(0) / (Ts * Ts).sum()
            lam = float(lam_i.mean())
            se = float(lam_i.std(ddof=1) / np.sqrt(N))
            drift = np.inf if lam_prev is None else abs(lam - lam_prev)
            lam_prev = lam
            tol = max(target_abs, target_rel * abs(lam))
            hist.append((T, lam, se, drift))
            # Trend check: consecutive trailing-half slopes share most of their
            # data, so `drift` is tiny even while the estimate slides like
            # a + b/T.  Fit the trailing half of the history to a + b/T and
            # require the extrapolated shift |lam - a| < tol as well.
            trend = 0.0
            if trend_check and len(hist) >= 4:
                h = np.array(hist)[len(hist) // 2:]
                A = np.stack([np.ones(len(h)), 1.0 / h[:, 0]], 1)
                lam_extrap = float(np.linalg.lstsq(A, h[:, 1], rcond=None)[0][0])
                trend = abs(lam - lam_extrap)
            el = time.time() - t0
            if log:
                log(f"    E={E:.4g} T={T:9.0f} lam={lam:+.3e} se={se:.1e} drift={drift:.1e} "
                    f"trend={trend:.1e} tol={tol:.1e} ({el:.0f}s, {el / (k + nburn) * 1e3:.3f} ms/step)")
            good = good + 1 if (se < tol and drift < tol and trend < tol) else 0
            if good >= n_consec:
                converged = True
                break
        if T >= tmax or (wall is not None and time.time() - t0 > wall):
            break
    Sh = be.S_host()
    E1 = 0.5 * np.sum(Sh[6:9] ** 2, 0) + V(Sh[0:3])
    return dict(E=E, lam_mean=lam, lam_se=se, lam_extrap=lam_extrap, lam_i=lam_i,
                converged=converged, T=T, t_burn=t_burn,
                tach_frac=tach_num / max(1, tach_den),
                ym_relerr_max=float(np.max(np.abs(E1 - E0) / E0)),
                history=np.array(hist), seconds=time.time() - t0,
                ms_per_step=(time.time() - t0) / (k + nburn) * 1e3)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--E', type=float, nargs='+', required=True)
    ap.add_argument('--N', type=int, default=100000)
    ap.add_argument('--backend', choices=['cupy', 'c', 'numpy'], default='cupy')
    ap.add_argument('--dt0', type=float, default=0.01)
    ap.add_argument('--target-abs', type=float, default=3e-7)
    ap.add_argument('--target-rel', type=float, default=0.05)
    ap.add_argument('--t-block', type=float, default=500.0)
    ap.add_argument('--t-burn', type=float, default=50.0, help='absolute time')
    ap.add_argument('--t-burn-ym', type=float, default=None,
                    help='burn-in in YM units (t_burn = this * E^-1/4); overrides --t-burn. '
                         'ym_initial is not microcanonical and channel depths are heavy-'
                         'tailed, so ~200-500 YM units is a safer default than 50 absolute.')
    ap.add_argument('--ghost', choices=['vec', 'diag'], default='vec',
                    help="vec: physical 3-vector ghost.  diag: off-diagonal Hessian entries "
                         "dropped (three scalar oscillators driven by K_xx, K_yy, K_zz).  "
                         "The ratio vec/diag says whether mode mixing (interior, degenerate "
                         "modes) or single-entry parametric driving controls lambda at that E.")
    ap.add_argument('--no-trend-check', action='store_true')
    ap.add_argument('--n-consec', type=int, default=3)
    ap.add_argument('--tmax', type=float, default=5e4)
    ap.add_argument('--wall', type=float, default=None, help='seconds per energy')
    ap.add_argument('--renorm-every', type=int, default=50)
    ap.add_argument('--omega-dt-max', type=float, default=0.25)
    ap.add_argument('--seed', type=int, default=777)
    ap.add_argument('--out', default='results_gpu')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()
    be = {'cupy': CupyBackend, 'c': CBackend, 'numpy': NumpyBackend}[a.backend]()
    os.makedirs(a.out, exist_ok=True)
    rows = []
    log = None if a.quiet else (lambda m: print(m, flush=True))
    for i, E in enumerate(a.E):
        dt = a.dt0 / max(1.0, E ** 0.25)
        r = run_energy(be, E, a.N, dt, a.seed + i, a.target_abs, a.target_rel, a.t_block,
                       a.t_burn, a.tmax, a.wall, a.n_consec, a.renorm_every,
                       a.omega_dt_max, log, od=1.0 if a.ghost == 'vec' else 0.0,
                       t_burn_ym=a.t_burn_ym, trend_check=not a.no_trend_check)
        verdict = ('zero@res' if abs(r['lam_mean']) < a.target_abs else
                   'positive' if r['lam_mean'] > 2 * r['lam_se'] else 'unresolved')
        status = 'converged' if r['converged'] else 'CAP HIT'
        print(f"E={E:9.4g} {a.ghost} lam={r['lam_mean']:+.3e}+-{r['lam_se']:.1e} "
              f"(extrap {r['lam_extrap']:+.3e}) T={r['T']:8.0f} "
              f"{status:9s} {verdict:10s} tach={r['tach_frac']:.4f} "
              f"ymerr={r['ym_relerr_max']:.1e} ({r['seconds']:.0f}s, "
              f"{r['ms_per_step']:.3f} ms/step)", flush=True)
        np.save(os.path.join(a.out, f'lam_i_E{E:.4g}.npy'), r['lam_i'])
        np.save(os.path.join(a.out, f'history_E{E:.4g}.npy'), r['history'])
        rows.append(dict(E=E, N=a.N, backend=be.name, ghost=a.ghost, lam_mean=r['lam_mean'],
                         lam_se=r['lam_se'], lam_extrap=r['lam_extrap'], t_burn=r['t_burn'],
                         lam_p10=float(np.percentile(r['lam_i'], 10)),
                         lam_p50=float(np.median(r['lam_i'])),
                         lam_p90=float(np.percentile(r['lam_i'], 90)),
                         T=r['T'], status=status, verdict=verdict,
                         tach_frac=r['tach_frac'], ym_relerr_max=r['ym_relerr_max'],
                         seconds=r['seconds'], ms_per_step=r['ms_per_step']))
        with open(os.path.join(a.out, 'transverse_gpu.csv'), 'w', newline='') as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            wr.writeheader()
            wr.writerows(rows)
    print("wrote", os.path.join(a.out, 'transverse_gpu.csv'))


if __name__ == '__main__':
    main()
