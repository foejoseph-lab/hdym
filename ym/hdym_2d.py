#!/usr/bin/env python3
"""
hdym_2d.py -- transverse (linearised-ghost) Lyapunov exponent of the TWO-FIELD model

    V(x, y) = 1/2 x^2 y^2 + 1/2 mu^2 (x^2 + y^2),        mu^2 = mu1^2 * sqrt(E),

i.e. the z = 0 invariant plane of SU(2) Yang-Mills mechanics, A^a_i = diag(x, y, 0), with a
small mass regulator that closes the hyperbolic channels.  The regulator is needed because
the pure 2-d x^2 y^2 model has a LOG-DIVERGENT density of states (channel cross-section
~ 1/|x| instead of 1/x^2): the microcanonical entry flux vanishes as T -> infinity and the
exponent drifts like 1/ln T.  Scaling mu^2 with sqrt(E) keeps the model exactly scale-free,
so every energy is the rescaled E = 1 system, as in 3-d.

Prediction (NOTE_2d.md; same derivation as the paper with ONE transverse action):

    lambda_diag(E) = (4 pi^2 / (15 sqrt2)) E^2 / Z2(mu1) * (1 - mu1^2 / (2 sqrt E))^{5/2}

with Z2(mu1) = 2 pi * area{V_1 < 1} the E = 1 density of states (Z2 = 189.7 at mu1 = 0.1,
275.3 at mu1 = 0.03).  Exponent 2 = 3/4 (kinematic) + 1/4 (one action) + 1 (driving power).
The last factor is the regulator correction to the passage integral (12% at E = 0.01,
mu1 = 0.1).  Validity: E >> mu1^4 / 4.

Ghost components.  w_x and w_y are driven resonantly by K_xx = y^2 + mu^2 and K_yy = x^2 + mu^2
inside the x- and y-channels (--ghost diag: off-diagonal 2xy dropped, each renormalised
separately, lambda_diag = mean of the two).  --ghost vec keeps the 2x2 coupling.  A THIRD
component w_z is always integrated: it is the z-ghost of the 3-d model at z = 0, driven by
K_zz = x^2 + y^2 + mu^2, which has NO resonant piece (all its channel power is detuned by
1/D from the parametric resonance).  The paper claims such detuned power is inert;
lambda_z / lambda_diag -> 0 as E -> 0 is the test.

Backends: --backend cupy (RTX, production) | c (gcc + ctypes, OpenMP) | numpy (reference).
State layout (10, N): p[2], w[3], P[2], Pi[3].

Usage (Spyder):
  %run hdym_2d.py --E 0.1 --N 100000 --target-rel 0.02 --target-abs 1e-6 --t-block 2000
  %run hdym_2d.py --E 0.03 --N 200000 --target-rel 0.02 --target-abs 1e-7 --t-block 10000
  %run hdym_2d.py --E 0.01 --N 200000 --target-rel 0.02 --target-abs 2e-8 --t-block 20000
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
STEP_SRC = r"""
#define YC0 0.6756035959798288
#define YC1 -0.17560359597982886
#define YC2 -0.17560359597982886
#define YC3 0.6756035959798288
#define YD0 1.3512071919596576
#define YD1 -1.7024143839193153
#define YD2 1.3512071919596576

DEVFN void kick(const double* p, const double* w, double* P, double* Pi, double h, double od, double mu2)
{
    double x = p[0], y = p[1];
    double x2 = x * x, y2 = y * y;
    /* grad V,  V = x2 y2 / 2 + mu2 (x2 + y2) / 2 */
    double gx = x * (y2 + mu2), gy = y * (x2 + mu2);
    /* Hess V . w for (wx, wy); od = 0 drops the off-diagonal 2xy entries.
       wz: the z-ghost of the 3-d model at z = 0, K_zz = x2 + y2 + mu2 (decoupled exactly). */
    double wx = w[0], wy = w[1], wz = w[2];
    double hx = (y2 + mu2) * wx + od * 2.0 * x * y * wy;
    double hy = od * 2.0 * x * y * wx + (x2 + mu2) * wy;
    double hz = (x2 + y2 + mu2) * wz;
    P[0] -= h * gx;  P[1] -= h * gy;
    Pi[0] -= h * (wx + hx);  Pi[1] -= h * (wy + hy);  Pi[2] -= h * (wz + hz);
}

DEVFN void yoshida(double* p, double* w, double* P, double* Pi, double h, double od, double mu2)
{
    int k;
    for (k = 0; k < 2; k++) p[k] += YC0 * h * P[k];
    for (k = 0; k < 3; k++) w[k] += YC0 * h * Pi[k];
    kick(p, w, P, Pi, YD0 * h, od, mu2);
    for (k = 0; k < 2; k++) p[k] += YC1 * h * P[k];
    for (k = 0; k < 3; k++) w[k] += YC1 * h * Pi[k];
    kick(p, w, P, Pi, YD1 * h, od, mu2);
    for (k = 0; k < 2; k++) p[k] += YC2 * h * P[k];
    for (k = 0; k < 3; k++) w[k] += YC2 * h * Pi[k];
    kick(p, w, P, Pi, YD2 * h, od, mu2);
    for (k = 0; k < 2; k++) p[k] += YC3 * h * P[k];
    for (k = 0; k < 3; k++) w[k] += YC3 * h * Pi[k];
}

/* logsum has 3 rows (x, y, z) of N: per-component accumulated log-growth. */
DEVFN void step_traj(double* S, double* logsum, int N, int i, double dt, int nsteps,
                     int renorm_every, double omega_dt_max, int accumulate, double od, double mu2)
{
    double p[2], w[3], P[2], Pi[3];
    int k, s, n, sub;
    for (k = 0; k < 2; k++) { p[k] = S[(0 + k) * N + i]; P[k] = S[(5 + k) * N + i]; }
    for (k = 0; k < 3; k++) { w[k] = S[(2 + k) * N + i]; Pi[k] = S[(7 + k) * N + i]; }
    double acc[3] = {0.0, 0.0, 0.0};
    for (s = 1; s <= nsteps; s++) {
        double om = fabs(p[0]);
        if (fabs(p[1]) > om) om = fabs(p[1]);
        n = (int)ceil(om * dt / omega_dt_max);
        if (n < 1) n = 1;
        if (n == 1) {
            yoshida(p, w, P, Pi, dt, od, mu2);
        } else {
            double h = dt / n;
            for (sub = 0; sub < n; sub++) yoshida(p, w, P, Pi, h, od, mu2);
        }
        if (s % renorm_every == 0) {
            if (od == 0.0) {
                for (k = 0; k < 2; k++) {
                    double nk = sqrt(w[k]*w[k] + Pi[k]*Pi[k]);
                    acc[k] += log(nk);
                    w[k] /= nk;  Pi[k] /= nk;
                }
            } else {
                double nrm = sqrt(w[0]*w[0] + w[1]*w[1] + Pi[0]*Pi[0] + Pi[1]*Pi[1]);
                acc[0] += log(nrm);  acc[1] += log(nrm);
                for (k = 0; k < 2; k++) { w[k] /= nrm; Pi[k] /= nrm; }
            }
            {
                double nz = sqrt(w[2]*w[2] + Pi[2]*Pi[2]);
                acc[2] += log(nz);
                w[2] /= nz;  Pi[2] /= nz;
            }
        }
    }
    if (accumulate) for (k = 0; k < 3; k++) logsum[k * N + i] += acc[k];
    for (k = 0; k < 2; k++) { S[(0 + k) * N + i] = p[k]; S[(5 + k) * N + i] = P[k]; }
    for (k = 0; k < 3; k++) { S[(2 + k) * N + i] = w[k]; S[(7 + k) * N + i] = Pi[k]; }
}
"""

CUDA_SRC = ("#define DEVFN __device__ __forceinline__\n" + STEP_SRC + r"""
extern "C" __global__
void advance(double* S, double* logsum, int N, double dt, int nsteps,
             int renorm_every, double omega_dt_max, int accumulate, double od, double mu2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) step_traj(S, logsum, N, i, dt, nsteps, renorm_every, omega_dt_max, accumulate, od, mu2);
}
""")

C_SRC = ("#include <math.h>\n#define DEVFN static inline\n" + STEP_SRC + r"""
void advance(double* S, double* logsum, int N, double dt, int nsteps,
             int renorm_every, double omega_dt_max, int accumulate, double od, double mu2)
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
    for (i = 0; i < N; i++)
        step_traj(S, logsum, N, i, dt, nsteps, renorm_every, omega_dt_max, accumulate, od, mu2);
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
        self.N = S_host.shape[1]
        self.logsum = self.cp.zeros((3, self.N))

    def advance(self, dt, nsteps, renorm_every, omega_dt_max, accumulate, od, mu2):
        N = self.N
        threads = 128
        blocks = (N + threads - 1) // threads
        self.kern((blocks,), (threads,),
                  (self.S, self.logsum, np.int32(N), np.float64(dt), np.int32(nsteps),
                   np.int32(renorm_every), np.float64(omega_dt_max), np.int32(accumulate),
                   np.float64(od), np.float64(mu2)))

    def logsum_host(self):
        return self.cp.asnumpy(self.logsum)

    def S_host(self):
        return self.cp.asnumpy(self.S)


class CBackend:
    name = 'c'

    def __init__(self):
        import subprocess
        d = tempfile.mkdtemp()
        src = os.path.join(d, 'hdym2d_step.c')
        lib = os.path.join(d, 'hdym2d_step' + ('.dll' if sys.platform == 'win32' else '.so'))
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
                                     ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_double,
                                     ctypes.c_double]
        print(f"C backend compiled with {cc}, {omp}", flush=True)

    def alloc(self, S_host):
        self.S = np.ascontiguousarray(S_host)
        self.N = S_host.shape[1]
        self.logsum = np.zeros((3, self.N))

    def advance(self, dt, nsteps, renorm_every, omega_dt_max, accumulate, od, mu2):
        dp = ctypes.POINTER(ctypes.c_double)
        self.lib.advance(self.S.ctypes.data_as(dp), self.logsum.ctypes.data_as(dp),
                         self.N, dt, nsteps, renorm_every, omega_dt_max, accumulate, od, mu2)

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
        self.N = S_host.shape[1]
        self.logsum = np.zeros((3, self.N))

    @staticmethod
    def _kick(p, w, P, Pi, h, od, mu2):
        x, y = p
        x2, y2 = x * x, y * y
        g = np.array([x * (y2 + mu2), y * (x2 + mu2)])
        wx, wy, wz = w
        H = np.array([(y2 + mu2) * wx + od * 2 * x * y * wy,
                      od * 2 * x * y * wx + (x2 + mu2) * wy,
                      (x2 + y2 + mu2) * wz])
        return P - h * g, Pi - h * (w + H)

    def _yoshida(self, p, w, P, Pi, h, od, mu2):
        for c, d in zip(self.YC, self.YD):
            p = p + c * h * P
            w = w + c * h * Pi
            if d != 0.0:
                P, Pi = self._kick(p, w, P, Pi, d * h, od, mu2)
        return p, w, P, Pi

    def advance(self, dt, nsteps, renorm_every, omega_dt_max, accumulate, od, mu2):
        S = self.S
        p, w, P, Pi = S[0:2], S[2:5], S[5:7], S[7:10]
        for s in range(1, nsteps + 1):
            om = np.max(np.abs(p), axis=0)
            nsub = np.maximum(1, np.ceil(om * dt / omega_dt_max)).astype(int)
            out = self._yoshida(p, w, P, Pi, dt, od, mu2)
            for n in np.unique(nsub):
                if n == 1:
                    continue
                m = nsub == n
                sub = (p[:, m], w[:, m], P[:, m], Pi[:, m])
                for _ in range(n):
                    sub = self._yoshida(*sub, dt / n, od, mu2)
                for a, b in zip(out, sub):
                    a[:, m] = b
            p, w, P, Pi = out
            if s % renorm_every == 0:
                if od == 0.0:
                    nrm = np.sqrt(w[:2] ** 2 + Pi[:2] ** 2)          # (2, N)
                    if accumulate:
                        self.logsum[:2] += np.log(nrm)
                    w[:2] /= nrm
                    Pi[:2] /= nrm
                else:
                    nrm = np.sqrt(np.sum(w[:2] ** 2, 0) + np.sum(Pi[:2] ** 2, 0))
                    if accumulate:
                        self.logsum[0] += np.log(nrm)
                        self.logsum[1] += np.log(nrm)
                    w[:2] /= nrm
                    Pi[:2] /= nrm
                nz = np.sqrt(w[2] ** 2 + Pi[2] ** 2)
                if accumulate:
                    self.logsum[2] += np.log(nz)
                w[2] /= nz
                Pi[2] /= nz
        self.S = np.concatenate([p, w, P, Pi], 0)

    def logsum_host(self):
        return self.logsum.copy()

    def S_host(self):
        return self.S.copy()


# ---------------------------------------------------------------- physics helpers (host)
def V(q, mu2):
    x, y = q
    return 0.5 * x * x * y * y + 0.5 * mu2 * (x * x + y * y)


def hessV_min_eig(q, mu2):
    x, y = q
    x2, y2 = x * x, y * y
    tr = x2 + y2
    disc = np.sqrt((y2 - x2) ** 2 + 16 * x2 * y2)
    return 0.5 * (tr - disc) + mu2


def Z2(mu1):
    """E = 1 density of states, Z2 = int delta(1 - H) d^2q d^2p = 2 pi * area{V_1 < 1}."""
    if mu1 <= 0.0:
        return np.inf
    from scipy.integrate import quad
    f = lambda x: np.sqrt((2 - mu1 * mu1 * x * x) / (x * x + mu1 * mu1))
    v, _ = quad(f, 0, np.sqrt(2) / mu1, limit=500, points=[1e-3, 1e-2, 0.1, 1, 10])
    return 2 * np.pi * 4 * v


def predicted_lambda(E, mu1):
    """lambda_diag = (4 pi^2 / 15 sqrt2) E^2 / Z2 * (1 - mu1^2/(2 sqrt E))^{5/2}."""
    Z = Z2(mu1)
    c = max(0.0, 1.0 - 0.5 * mu1 * mu1 / np.sqrt(E))
    # ghost mass is 1 + mu^2: resonance at u = 2 sqrt(1 + mu^2) E^{-1/4} and lambda = S_K(2 w0)/(8 w0^2)
    # -> factor (1 + mu^2)^{-7/2} (1% at E = 0.1, mu1 = 0.1).
    mu2 = mu1 * mu1 * np.sqrt(E)
    return 4 * np.pi ** 2 / (15 * np.sqrt(2)) * E ** 2 / Z * c ** 2.5 * (1 + mu2) ** -3.5


def initial_state(N, E, mu1, seed):
    """Random phase points on the E = 1 shell of V_1 (mu = mu1), then the exact scaling
    q -> E^{1/4} q, P -> E^{1/2} P, mu^2 = mu1^2 sqrt(E) takes them to the E shell."""
    rng = np.random.default_rng(seed)
    q = rng.standard_normal((2, N))
    P = rng.standard_normal((2, N))
    mu2_1 = mu1 * mu1
    s = np.ones(N)
    for _ in range(80):                       # fixed-point iteration for the amplitude
        E1 = 0.5 * s ** 4 * np.sum(P * P, 0) + V(s * q, mu2_1)
        s *= (1.0 / E1) ** 0.25 * 0.5 + 0.5
    q = q * s * E ** 0.25
    P = P * s ** 2 * E ** 0.5
    w = rng.standard_normal((3, N))
    Pi = rng.standard_normal((3, N))
    nrm = np.sqrt(w * w + Pi * Pi)
    w /= nrm
    Pi /= nrm
    return np.concatenate([q, w, P, Pi], 0)


# ---------------------------------------------------------------- adaptive run
def run_energy(be, E, mu1, N, dt, seed, target_abs, target_rel, t_block, t_burn, tmax, wall,
               n_consec, renorm_every, omega_dt_max, log, od, t_burn_ym=None, trend_check=True):
    t0 = time.time()
    mu2 = mu1 * mu1 * np.sqrt(E)
    if t_burn_ym is not None:
        t_burn = t_burn_ym * E ** -0.25
    S0 = initial_state(N, E, mu1, seed)
    E0 = 0.5 * np.sum(S0[5:7] ** 2, 0) + V(S0[0:2], mu2)
    be.alloc(S0)
    nburn = int(np.ceil(t_burn / dt / renorm_every)) * renorm_every
    nblock = int(np.ceil(t_block / dt / renorm_every)) * renorm_every
    be.advance(dt, nburn, renorm_every, omega_dt_max, 0, od, mu2)
    k = 0
    T = 0.0
    chk_T, chk_S, hist = [], [], []
    tach_num, tach_den = 0.0, 0
    lam_prev, good, converged = None, 0, False
    lam_i = np.full((3, N), np.nan)
    lam = se = lam_extrap = np.nan
    lam3 = np.full(3, np.nan)
    while True:
        be.advance(dt, nblock, renorm_every, omega_dt_max, 1, od, mu2)
        k += nblock
        T = k * dt
        S = be.logsum_host()
        if not np.all(np.isfinite(S)):
            raise FloatingPointError(f"non-finite log-growth at E={E}, T={T}")
        chk_T.append(T)
        chk_S.append(S)
        Sh = be.S_host()
        mn = hessV_min_eig(Sh[0:2, :min(N, 2048)], mu2)
        tach_num += float(np.mean(1.0 + mn <= 0.0))
        tach_den += 1
        Ta = np.array(chk_T)
        i0 = max(0, int(len(Ta) / 2) - 1)
        if len(Ta) - i0 >= 3:
            Ts = Ta[i0:] - Ta[i0:].mean()
            Ss = np.array(chk_S[i0:])                       # (nchk, 3, N)
            lam_i = (Ts[:, None, None] * (Ss - Ss.mean(0))).sum(0) / (Ts * Ts).sum()
            lam3 = lam_i.mean(1)
            lam_diag_i = 0.5 * (lam_i[0] + lam_i[1])      # per-trajectory single-mode rate
            lam = float(lam_diag_i.mean())
            se = float(lam_diag_i.std(ddof=1) / np.sqrt(N))
            drift = np.inf if lam_prev is None else abs(lam - lam_prev)
            lam_prev = lam
            tol = max(target_abs, target_rel * abs(lam))
            hist.append((T, lam, se, drift, lam3[0], lam3[1], lam3[2]))
            trend = 0.0
            if trend_check and len(hist) >= 4:
                h = np.array(hist)[len(hist) // 2:]
                A = np.stack([np.ones(len(h)), 1.0 / h[:, 0]], 1)
                lam_extrap = float(np.linalg.lstsq(A, h[:, 1], rcond=None)[0][0])
                trend = abs(lam - lam_extrap)
            el = time.time() - t0
            if log:
                log(f"    E={E:.4g} T={T:9.0f} lam_xy={lam:+.3e} se={se:.1e} drift={drift:.1e} "
                    f"trend={trend:.1e} tol={tol:.1e} | lam_x={lam3[0]:+.2e} lam_y={lam3[1]:+.2e} "
                    f"lam_z={lam3[2]:+.2e} ({el:.0f}s, {el / (k + nburn) * 1e3:.3f} ms/step)")
            good = good + 1 if (se < tol and drift < tol and trend < tol) else 0
            if good >= n_consec:
                converged = True
                break
        if T >= tmax or (wall is not None and time.time() - t0 > wall):
            break
    Sh = be.S_host()
    E1 = 0.5 * np.sum(Sh[5:7] ** 2, 0) + V(Sh[0:2], mu2)
    return dict(E=E, mu1=mu1, lam_mean=lam, lam_se=se, lam_extrap=lam_extrap, lam_i=lam_i,
                lam_x=float(lam3[0]), lam_y=float(lam3[1]), lam_z=float(lam3[2]),
                lam_z_se=float(lam_i[2].std(ddof=1) / np.sqrt(N)) if np.all(np.isfinite(lam_i[2])) else np.nan,
                converged=converged, T=T, t_burn=t_burn,
                tach_frac=tach_num / max(1, tach_den),
                ym_relerr_max=float(np.max(np.abs(E1 - E0) / E0)),
                history=np.array(hist), seconds=time.time() - t0,
                ms_per_step=(time.time() - t0) / (k + nburn) * 1e3)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--E', type=float, nargs='+', required=True)
    ap.add_argument('--mu1', type=float, default=0.1,
                    help='regulator at E = 1; mu^2 = mu1^2 sqrt(E).  0 = pure x^2y^2 (log-divergent Z2; '
                         'lambda then drifts as 1/ln T -- for demonstrating that only).')
    ap.add_argument('--N', type=int, default=100000)
    ap.add_argument('--backend', choices=['cupy', 'c', 'numpy'], default='cupy')
    ap.add_argument('--dt0', type=float, default=0.01)
    ap.add_argument('--target-abs', type=float, default=3e-7)
    ap.add_argument('--target-rel', type=float, default=0.02)
    ap.add_argument('--t-block', type=float, default=500.0)
    ap.add_argument('--t-burn', type=float, default=50.0, help='absolute time')
    ap.add_argument('--t-burn-ym', type=float, default=300.0,
                    help='burn-in in YM units (t_burn = this * E^-1/4); default 300.  '
                         'Channel visits last up to 2 sqrt2 * sqrt2/mu1 YM units (40 at mu1 = 0.1, '
                         '130 at 0.03), so use >= 500 for mu1 <= 0.03.')
    ap.add_argument('--ghost', choices=['vec', 'diag'], default='diag',
                    help="diag: w_x, w_y independent scalar oscillators (single-mode law, the prediction). "
                         "vec: physical 2x2 coupled ghost.  w_z is integrated separately in both cases.")
    ap.add_argument('--no-trend-check', action='store_true')
    ap.add_argument('--n-consec', type=int, default=3)
    ap.add_argument('--tmax', type=float, default=5e5)
    ap.add_argument('--wall', type=float, default=None, help='seconds per energy')
    ap.add_argument('--renorm-every', type=int, default=50)
    ap.add_argument('--omega-dt-max', type=float, default=0.25)
    ap.add_argument('--seed', type=int, default=777)
    ap.add_argument('--out', default='results_2d')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()
    be = {'cupy': CupyBackend, 'c': CBackend, 'numpy': NumpyBackend}[a.backend]()
    os.makedirs(a.out, exist_ok=True)
    Z = Z2(a.mu1)
    print(f"mu1 = {a.mu1}: Z2 = {Z:.2f}, nu0 = 4pi/Z2 = {4 * np.pi / Z:.4f}, "
          f"C2 = 4pi^2/(15 sqrt2 Z2) = {4 * np.pi ** 2 / (15 * np.sqrt(2)) / Z:.4e}", flush=True)
    rows = []
    log = None if a.quiet else (lambda m: print(m, flush=True))
    for i, E in enumerate(a.E):
        dt = a.dt0 / max(1.0, E ** 0.25)
        pred = predicted_lambda(E, a.mu1)
        print(f"--- E = {E:.4g}: predicted lambda_diag = {pred:.4e} "
              f"(regulator factor {(max(0, 1 - 0.5 * a.mu1 ** 2 / np.sqrt(E))) ** 2.5:.3f})", flush=True)
        r = run_energy(be, E, a.mu1, a.N, dt, a.seed + i, a.target_abs, a.target_rel, a.t_block,
                       a.t_burn, a.tmax, a.wall, a.n_consec, a.renorm_every, a.omega_dt_max, log,
                       od=1.0 if a.ghost == 'vec' else 0.0, t_burn_ym=a.t_burn_ym,
                       trend_check=not a.no_trend_check)
        verdict = ('zero@res' if abs(r['lam_mean']) < a.target_abs else
                   'positive' if r['lam_mean'] > 2 * r['lam_se'] else 'unresolved')
        status = 'converged' if r['converged'] else 'CAP HIT'
        print(f"E={E:9.4g} mu1={a.mu1} {a.ghost} lam_xy={r['lam_mean']:+.3e}+-{r['lam_se']:.1e} "
              f"(extrap {r['lam_extrap']:+.3e}) pred={pred:.3e} RATIO={r['lam_mean'] / pred:.3f} | "
              f"lam_z={r['lam_z']:+.2e}+-{r['lam_z_se']:.1e} (z/xy={r['lam_z'] / r['lam_mean']:.3f}) "
              f"T={r['T']:8.0f} {status:9s} {verdict:10s} tach={r['tach_frac']:.4f} "
              f"ymerr={r['ym_relerr_max']:.1e} ({r['seconds']:.0f}s, {r['ms_per_step']:.3f} ms/step)",
              flush=True)
        tag = f'E{E:.4g}_mu{a.mu1:g}'
        np.save(os.path.join(a.out, f'lam_i_{tag}.npy'), r['lam_i'])
        np.save(os.path.join(a.out, f'history_{tag}.npy'), r['history'])
        rows.append(dict(E=E, mu1=a.mu1, N=a.N, backend=be.name, ghost=a.ghost,
                         lam_mean=r['lam_mean'], lam_se=r['lam_se'], lam_extrap=r['lam_extrap'],
                         pred=pred, ratio=r['lam_mean'] / pred,
                         lam_x=r['lam_x'], lam_y=r['lam_y'], lam_z=r['lam_z'], lam_z_se=r['lam_z_se'],
                         t_burn=r['t_burn'], T=r['T'], status=status, verdict=verdict,
                         tach_frac=r['tach_frac'], ym_relerr_max=r['ym_relerr_max'],
                         seconds=r['seconds'], ms_per_step=r['ms_per_step']))
        # one CSV for all runs: append a row per energy (header written once), so separate
        # invocations at the same mu1 do not overwrite each other
        path = os.path.join(a.out, 'transverse_2d.csv')
        new = not os.path.exists(path)
        with open(path, 'a', newline='') as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[-1].keys()))
            if new:
                wr.writeheader()
            wr.writerow(rows[-1])
    print("appended to", os.path.join(a.out, 'transverse_2d.csv'))


if __name__ == '__main__':
    main()
