#!/usr/bin/env python3
"""
hdym_rung1_gpu.py -- Rung 1 (a): Lyapunov exponents of the OFF-DIAGONAL ghost sector along
diagonal Yang-Mills backgrounds, large ensemble, GPU.

The linearised full homogeneous SU(2) Lee-Wick system about A = diag(x, y, z) (hdym_rung1.py)
splits into three equivalent 8-dim blocks; this integrates block (u12, u21):

    state  v = (u12, u21, u12', u21', u12'', u21'', u12''', u21''')
    v' = K(q, p) v ,   q'' = -grad V(q)   (ghost-free background, pure YM mechanics)

with the two leading Lyapunov exponents from two vectors (Gram-Schmidt).  The block has one
gauge zero mode, one conserved Gauss charge, one YM-type dof (neutral: the off-diagonal
directions are the 3 gauge + 3 spatial-rotation directions of the SVD A = O_c diag O_s) and
two ghost dof.  The two positive exponents are therefore the off-diagonal ghost rates
(g1 >= g2); compare with the diagonal-sector lambda_diag = 0.0120 E^{9/4} and the physical
three-vector lambda_perp ~ 2.2 lambda_diag.  Question: do they scale as E^{9/4}, or as
E^{5/4} (u^-2 tail of the xy entry driving a degenerate ghost pair)?

Background: Yoshida-4 in two half-steps; tangent: RK4 with K at the three background points.
Sub-stepping as in hdym_gpu.py.  CPU results (hdym_rung1_lyap.py, N = 18, T = 3000 YM units):
E = 1: g1 = 0.0126(9), g2 = 0.0070(5); E = 0.1: g1 <= 0.002 (noise floor).

Usage (Spyder):
  %run hdym_rung1_gpu.py --E 1 0.3 --N 20000 --t-block 2000 --target-rel 0.02 --target-abs 1e-4
  %run hdym_rung1_gpu.py --E 0.1 0.03 --N 100000 --t-block 10000 --target-rel 0.02 --target-abs 1e-6
"""
import argparse
import csv
import ctypes
import os
import sys
import tempfile
import time

import numpy as np

HIST_SHAPE = (64, 50)          # (depth D in [0, 8], action J in [0, 1]) bins; must match HIST_* in the C source
HIST_EDGES = (np.linspace(0, 8, 65), np.linspace(0, 1, 51))

KBLOCK = r"""
    double x0 = pow(z, 2);
    double x1 = v0*x0;
    double x2 = pow(z, 4);
    double x3 = 4*pz*z;
    double x4 = 2*pow(pz, 2);
    double x5 = 2*x0;
    double x6 = pow(x, 2);
    double x7 = cj*x6;
    double x8 = pow(y, 3);
    double x9 = pow(x, 3);
    double x10 = x9*y;
    double x11 = 2*cj;
    double x12 = px*py;
    double x13 = pow(y, 2);
    double x14 = x13*x6;
    double x15 = x13*x7;
    double x16 = v1*x0;
    a6 = 2*cj*pow(px, 2)*v0 + 2*cj*px*v3*y + cj*v1*x*x0*y + cj*v1*x9*y + cj*v5*x*y + 2*px*py*v1 - px*v2*x*x11 + 2*px*v3*y + 2*py*v3*x + 2*v0*x0*x13 + 2*v0*x0*x6 - v0*x14 - v0*x15 - v0*x2 - v0*x4 - v1*x*x8 + v1*x*y - v1*x10 - v1*x11*x12 - v2*x3 - v4*x5 - v4*x7 - v4 + 2*v5*x*y - x1*x7 - x1;
    a7 = 2*cj*pow(py, 2)*v1 + 2*cj*py*v2*x + cj*v0*x*x0*y + cj*v0*x*x8 + cj*v4*x*y - cj*v5*x13 - cj*x13*x16 + 2*px*py*v0 + 2*px*v2*y + 2*py*v2*x - py*v3*x11*y - v0*x*x8 + v0*x*y - v0*x10 - v0*x11*x12 + 2*v1*x0*x13 + 2*v1*x0*x6 - v1*x14 - v1*x15 - v1*x2 - v1*x4 - v3*x3 + 2*v4*x*y - v5*x5 - v5 - x16;
"""

STEP_SRC = r"""
#define HIST_ND 64
#define HIST_NJ 50
#define HIST_DMAX 8.0
#define HIST_JMAX 1.0
#define YC0 0.6756035959798288
#define YC1 -0.17560359597982886
#define YC2 -0.17560359597982886
#define YC3 0.6756035959798288
#define YD0 1.3512071919596576
#define YD1 -1.7024143839193153
#define YD2 1.3512071919596576

DEVFN void gradV(const double* q, double* g)
{
    double x = q[0], y = q[1], z = q[2];
    g[0] = x * (y*y + z*z);  g[1] = y * (z*z + x*x);  g[2] = z * (x*x + y*y);
}

DEVFN void yoshida(double* q, double* p, double h)
{
    double g[3]; int k;
    for (k = 0; k < 3; k++) q[k] += YC0 * h * p[k];
    gradV(q, g); for (k = 0; k < 3; k++) p[k] -= YD0 * h * g[k];
    for (k = 0; k < 3; k++) q[k] += YC1 * h * p[k];
    gradV(q, g); for (k = 0; k < 3; k++) p[k] -= YD1 * h * g[k];
    for (k = 0; k < 3; k++) q[k] += YC2 * h * p[k];
    gradV(q, g); for (k = 0; k < 3; k++) p[k] -= YD2 * h * g[k];
    for (k = 0; k < 3; k++) q[k] += YC3 * h * p[k];
}

/* K(q,p) v for the (u12, u21) block */
DEVFN void Kv(const double* q, const double* p, const double* v, double* out, double cj)
{
    double x = q[0], y = q[1], z = q[2], px = p[0], py = p[1], pz = p[2];
    double v0 = v[0], v1 = v[1], v2 = v[2], v3 = v[3], v4 = v[4], v5 = v[5], v6 = v[6], v7 = v[7];
    double a6, a7;
""" + KBLOCK + r"""
    out[0] = v2; out[1] = v3; out[2] = v4; out[3] = v5; out[4] = v6; out[5] = v7;
    out[6] = a6; out[7] = a7;
}

DEVFN void rk4_pair(double* q, double* p, double* v, double* w, double h, double cj)
{
    /* background to half step and full step; tangent RK4 */
    double qh[3], ph[3], q1[3], p1[3];
    double k1[8], k2[8], k3[8], k4[8], tmp[8];
    double l1[8], l2[8], l3[8], l4[8];
    int k;
    for (k = 0; k < 3; k++) { qh[k] = q[k]; ph[k] = p[k]; }
    yoshida(qh, ph, 0.5 * h);
    for (k = 0; k < 3; k++) { q1[k] = qh[k]; p1[k] = ph[k]; }
    yoshida(q1, p1, 0.5 * h);
    Kv(q, p, v, k1, cj);  Kv(q, p, w, l1, cj);
    for (k = 0; k < 8; k++) tmp[k] = v[k] + 0.5 * h * k1[k];
    Kv(qh, ph, tmp, k2, cj);
    for (k = 0; k < 8; k++) tmp[k] = w[k] + 0.5 * h * l1[k];
    Kv(qh, ph, tmp, l2, cj);
    for (k = 0; k < 8; k++) tmp[k] = v[k] + 0.5 * h * k2[k];
    Kv(qh, ph, tmp, k3, cj);
    for (k = 0; k < 8; k++) tmp[k] = w[k] + 0.5 * h * l2[k];
    Kv(qh, ph, tmp, l3, cj);
    for (k = 0; k < 8; k++) tmp[k] = v[k] + h * k3[k];
    Kv(q1, p1, tmp, k4, cj);
    for (k = 0; k < 8; k++) tmp[k] = w[k] + h * l3[k];
    Kv(q1, p1, tmp, l4, cj);
    for (k = 0; k < 8; k++) {
        v[k] += h / 6.0 * (k1[k] + 2.0 * k2[k] + 2.0 * k3[k] + k4[k]);
        w[k] += h / 6.0 * (l1[k] + 2.0 * l2[k] + 2.0 * l3[k] + l4[k]);
    }
    for (k = 0; k < 3; k++) { q[k] = q1[k]; p[k] = p1[k]; }
}

/* S layout (22, N): q[3], p[3], v[8], w[8].  logsum (2, N). */
DEVFN void step_traj(double* S, double* logsum, int N, int i, double dt, int nsteps,
                     int renorm_every, double omega_dt_max, int accumulate, double cj, double a_split,
                     double* hist_g, double* hist_t, double qscale, double jscale)
{
    double q[3], p[3], v[8], w[8];
    int k, s, n, sub;
    for (k = 0; k < 3; k++) { q[k] = S[k * N + i]; p[k] = S[(3 + k) * N + i]; }
    for (k = 0; k < 8; k++) { v[k] = S[(6 + k) * N + i]; w[k] = S[(14 + k) * N + i]; }
    double acc1 = 0.0, acc2 = 0.0, acc_in = 0.0, acc_out = 0.0, n_in = 0.0, n_out = 0.0;
    for (s = 1; s <= nsteps; s++) {
        double om = fabs(q[0]);
        if (fabs(q[1]) > om) om = fabs(q[1]);
        if (fabs(q[2]) > om) om = fabs(q[2]);
        n = (int)ceil(om * dt / omega_dt_max);
        if (n < 1) n = 1;
        if (n == 1) {
            rk4_pair(q, p, v, w, dt, cj);
        } else {
            double h = dt / n;
            for (sub = 0; sub < n; sub++) rk4_pair(q, p, v, w, h, cj);
        }
        if (s % renorm_every == 0) {
            double nv = 0.0, dot = 0.0, nw = 0.0;
            for (k = 0; k < 8; k++) nv += v[k] * v[k];
            nv = sqrt(nv);
            for (k = 0; k < 8; k++) v[k] /= nv;
            for (k = 0; k < 8; k++) dot += w[k] * v[k];
            for (k = 0; k < 8; k++) w[k] -= dot * v[k];
            for (k = 0; k < 8; k++) nw += w[k] * w[k];
            nw = sqrt(nw);
            for (k = 0; k < 8; k++) w[k] /= nw;
            acc1 += log(nv);  acc2 += log(nw);
            {
                double mx = fabs(q[0]);
                if (fabs(q[1]) > mx) mx = fabs(q[1]);
                if (fabs(q[2]) > mx) mx = fabs(q[2]);
                if (mx < a_split) { acc_in += log(nv); n_in += 1.0; } else { acc_out += log(nv); n_out += 1.0; }
                if (accumulate && hist_g != 0) {
                    /* growth landscape: depth D = max|q| E^-1/4, transverse action J = E_perp/|q_i| E^-3/4 */
                    int i0 = 0; if (fabs(q[1]) > fabs(q[0])) i0 = 1; if (fabs(q[2]) > fabs(q[i0])) i0 = 2;
                    int j0 = (i0 + 1) % 3, k0 = (i0 + 2) % 3;
                    double qi = fabs(q[i0]);
                    double Eperp = 0.5 * (p[j0]*p[j0] + p[k0]*p[k0]) + 0.5 * qi*qi * (q[j0]*q[j0] + q[k0]*q[k0]) + 0.5 * q[j0]*q[j0]*q[k0]*q[k0];
                    double D = qi * qscale;
                    double Jv = (qi > 0.0) ? Eperp / qi * jscale : 0.0;
                    int bd = (int)(D / HIST_DMAX * HIST_ND); if (bd >= HIST_ND) bd = HIST_ND - 1;
                    int bj = (int)(Jv / HIST_JMAX * HIST_NJ); if (bj >= HIST_NJ) bj = HIST_NJ - 1;
                    ATOMIC_ADD(&hist_g[bd * HIST_NJ + bj], log(nv));
                    ATOMIC_ADD(&hist_t[bd * HIST_NJ + bj], 1.0);
                }
            }
        }
    }
    if (accumulate) { logsum[i] += acc1; logsum[N + i] += acc2; logsum[2 * N + i] += acc_in; logsum[3 * N + i] += acc_out;
                      logsum[4 * N + i] += n_in; logsum[5 * N + i] += n_out; }
    for (k = 0; k < 3; k++) { S[k * N + i] = q[k]; S[(3 + k) * N + i] = p[k]; }
    for (k = 0; k < 8; k++) { S[(6 + k) * N + i] = v[k]; S[(14 + k) * N + i] = w[k]; }
}
"""

CUDA_SRC = ("#define DEVFN __device__ __forceinline__\n#define ATOMIC_ADD(p, v) atomicAdd((p), (v))\n" + STEP_SRC + r"""
extern "C" __global__
void advance(double* S, double* logsum, int N, double dt, int nsteps,
             int renorm_every, double omega_dt_max, int accumulate, double cj, double a_split,
             double* hist_g, double* hist_t, double qscale, double jscale)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) step_traj(S, logsum, N, i, dt, nsteps, renorm_every, omega_dt_max, accumulate, cj, a_split,
                         hist_g, hist_t, qscale, jscale);
}
""")

C_SRC = ("#include <math.h>\n#define DEVFN static inline\n#define ATOMIC_ADD(p, v) do { _Pragma(\"omp atomic\") *(p) += (v); } while (0)\n" + STEP_SRC + r"""
void advance(double* S, double* logsum, int N, double dt, int nsteps,
             int renorm_every, double omega_dt_max, int accumulate, double cj, double a_split,
             double* hist_g, double* hist_t, double qscale, double jscale)
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
    for (i = 0; i < N; i++)
        step_traj(S, logsum, N, i, dt, nsteps, renorm_every, omega_dt_max, accumulate, cj, a_split,
                  hist_g, hist_t, qscale, jscale);
}
""")


class CupyBackend:
    name = 'cupy'

    def __init__(self):
        import cupy as cp
        self.cp = cp
        self.kern = cp.RawKernel(CUDA_SRC, 'advance')
        dev = cp.cuda.Device()
        free, total = dev.mem_info
        print(f"CuPy {cp.__version__} on device {dev.id}, {free / 1e9:.1f}/{total / 1e9:.1f} GB free", flush=True)

    def alloc(self, S_host):
        self.S = self.cp.asarray(S_host)
        self.N = S_host.shape[1]
        self.logsum = self.cp.zeros((6, self.N))
        self.hist_g = self.cp.zeros(HIST_SHAPE)
        self.hist_t = self.cp.zeros(HIST_SHAPE)

    def advance(self, dt, nsteps, renorm_every, omega_dt_max, accumulate, cj=1.0, a_split=1e9, qscale=1.0, jscale=1.0):
        N = self.N
        threads = 128
        blocks = (N + threads - 1) // threads
        self.kern((blocks,), (threads,),
                  (self.S, self.logsum, np.int32(N), np.float64(dt), np.int32(nsteps),
                   np.int32(renorm_every), np.float64(omega_dt_max), np.int32(accumulate), np.float64(cj), np.float64(a_split),
                   self.hist_g, self.hist_t, np.float64(qscale), np.float64(jscale)))

    def hist_host(self):
        return self.cp.asnumpy(self.hist_g), self.cp.asnumpy(self.hist_t)

    def logsum_host(self):
        return self.cp.asnumpy(self.logsum)

    def S_host(self):
        return self.cp.asnumpy(self.S)


class CBackend:
    name = 'c'

    def __init__(self):
        import subprocess
        d = tempfile.mkdtemp()
        src = os.path.join(d, 'rung1_step.c')
        lib = os.path.join(d, 'rung1_step' + ('.dll' if sys.platform == 'win32' else '.so'))
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
                                     ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_double, ctypes.c_double,
                                     dp, dp, ctypes.c_double, ctypes.c_double]
        print(f"C backend compiled with {cc}, {omp}", flush=True)

    def alloc(self, S_host):
        self.S = np.ascontiguousarray(S_host)
        self.N = S_host.shape[1]
        self.logsum = np.zeros((6, self.N))
        self.hist_g = np.zeros(HIST_SHAPE)
        self.hist_t = np.zeros(HIST_SHAPE)

    def advance(self, dt, nsteps, renorm_every, omega_dt_max, accumulate, cj=1.0, a_split=1e9, qscale=1.0, jscale=1.0):
        dp = ctypes.POINTER(ctypes.c_double)
        self.lib.advance(self.S.ctypes.data_as(dp), self.logsum.ctypes.data_as(dp),
                         self.N, dt, nsteps, renorm_every, omega_dt_max, accumulate, cj, a_split,
                         self.hist_g.ctypes.data_as(dp), self.hist_t.ctypes.data_as(dp), qscale, jscale)

    def hist_host(self):
        return self.hist_g.copy(), self.hist_t.copy()

    def logsum_host(self):
        return self.logsum.copy()

    def S_host(self):
        return self.S.copy()


def V(q):
    x, y, z = q
    return 0.5 * (x * x * y * y + y * y * z * z + z * z * x * x)


def initial_state(N, E, seed):
    rng = np.random.default_rng(seed)
    q = rng.standard_normal((3, N))
    p = rng.standard_normal((3, N))
    lam = (E / (0.5 * np.sum(p * p, 0) + V(q))) ** 0.25
    q *= lam
    p *= lam ** 2
    v = rng.standard_normal((8, N))
    w = rng.standard_normal((8, N))
    v /= np.sqrt(np.sum(v * v, 0))
    w -= np.sum(w * v, 0) * v
    w /= np.sqrt(np.sum(w * w, 0))
    return np.concatenate([q, p, v, w], 0)


def run_energy(be, E, N, dt, seed, target_abs, target_rel, t_block, t_burn_ym, tmax, wall,
               n_consec, renorm_every, omega_dt_max, log, cj=1.0, a_split_ym=3.0):
    t0 = time.time()
    t_burn = t_burn_ym * E ** -0.25
    a_split = a_split_ym * E ** 0.25        # channel mouth in physical units
    S0 = initial_state(N, E, seed)
    E0 = 0.5 * np.sum(S0[3:6] ** 2, 0) + V(S0[0:3])
    be.alloc(S0)
    nburn = int(np.ceil(t_burn / dt / renorm_every)) * renorm_every
    nblock = int(np.ceil(t_block / dt / renorm_every)) * renorm_every
    be.advance(dt, nburn, renorm_every, omega_dt_max, 0, cj, a_split, E ** -0.25, E ** -0.75)
    k, T = 0, 0.0
    chk_T, chk_S, hist = [], [], []
    lam_prev, good, converged = None, 0, False
    lam_i = np.full((2, N), np.nan)
    lam = se = lam_extrap = np.nan
    lam2 = np.nan
    while True:
        be.advance(dt, nblock, renorm_every, omega_dt_max, 1, cj, a_split, E ** -0.25, E ** -0.75)
        k += nblock
        T = k * dt
        S = be.logsum_host()
        if not np.all(np.isfinite(S)):
            raise FloatingPointError(f"non-finite log-growth at E={E}, T={T}")
        chk_T.append(T)
        chk_S.append(S)
        Ta = np.array(chk_T)
        i0 = max(0, int(len(Ta) / 2) - 1)
        if len(Ta) - i0 >= 3:
            Ts = Ta[i0:] - Ta[i0:].mean()
            Ss = np.array(chk_S[i0:])
            lam_i = (Ts[:, None, None] * (Ss - Ss.mean(0))).sum(0) / (Ts * Ts).sum()
            lam = float(lam_i[0].mean())
            lam2 = float(lam_i[1].mean())
            se = float(lam_i[0].std(ddof=1) / np.sqrt(N))
            drift = np.inf if lam_prev is None else abs(lam - lam_prev)
            lam_prev = lam
            tol = max(target_abs, target_rel * abs(lam))
            hist.append((T, lam, se, drift, lam2))
            trend = 0.0
            if len(hist) >= 4:
                h = np.array(hist)[len(hist) // 2:]
                A = np.stack([np.ones(len(h)), 1.0 / h[:, 0]], 1)
                lam_extrap = float(np.linalg.lstsq(A, h[:, 1], rcond=None)[0][0])
                trend = abs(lam - lam_extrap)
            el = time.time() - t0
            if log:
                log(f"    E={E:.4g} T={T:9.0f} g1={lam:+.3e} se={se:.1e} drift={drift:.1e} trend={trend:.1e} "
                    f"tol={tol:.1e} | g2={lam2:+.3e} ({el:.0f}s, {el / (k + nburn) * 1e3:.3f} ms/step)")
            good = good + 1 if (se < tol and drift < tol and trend < tol) else 0
            if good >= n_consec:
                converged = True
                break
        if T >= tmax or (wall is not None and time.time() - t0 > wall):
            break
    Sh = be.S_host()
    E1 = 0.5 * np.sum(Sh[3:6] ** 2, 0) + V(Sh[0:3])
    Sl = be.logsum_host()
    tot_in, tot_out = float(Sl[2].sum()), float(Sl[3].sum())
    frac_growth_channel = tot_out / (tot_in + tot_out) if (tot_in + tot_out) != 0 else np.nan
    frac_time_channel = float(Sl[5].sum() / (Sl[4].sum() + Sl[5].sum()))
    hg, ht = be.hist_host()
    return dict(E=E, hist_g=hg, hist_t=ht, frac_growth_channel=frac_growth_channel, frac_time_channel=frac_time_channel, g1=lam, g1_se=se, g1_extrap=lam_extrap, g2=lam2,
                g2_se=float(lam_i[1].std(ddof=1) / np.sqrt(N)), lam_i=lam_i, converged=converged, T=T,
                ym_relerr_max=float(np.max(np.abs(E1 - E0) / E0)), history=np.array(hist),
                seconds=time.time() - t0, ms_per_step=(time.time() - t0) / (k + nburn) * 1e3)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--E', type=float, nargs='+', required=True)
    ap.add_argument('--N', type=int, default=20000)
    ap.add_argument('--backend', choices=['cupy', 'c'], default='cupy')
    ap.add_argument('--dt0', type=float, default=0.01)
    ap.add_argument('--target-abs', type=float, default=1e-5)
    ap.add_argument('--target-rel', type=float, default=0.02)
    ap.add_argument('--t-block', type=float, default=1000.0)
    ap.add_argument('--t-burn-ym', type=float, default=300.0)
    ap.add_argument('--n-consec', type=int, default=3)
    ap.add_argument('--tmax', type=float, default=5e5)
    ap.add_argument('--wall', type=float, default=None)
    ap.add_argument('--renorm-every', type=int, default=20)
    ap.add_argument('--omega-dt-max', type=float, default=0.25)
    ap.add_argument('--cj', type=float, default=1.0,
                    help='coefficient of the (D^mu F_mu0)^2 Gauss-charge term of the Lee-Wick action: 1 = the theory; 0 = diagnostic (drops the j^2 term)')
    ap.add_argument('--seed', type=int, default=4242)
    ap.add_argument('--out', default='results_rung1')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()
    be = {'cupy': CupyBackend, 'c': CBackend}[a.backend]()
    os.makedirs(a.out, exist_ok=True)
    log = None if a.quiet else (lambda m: print(m, flush=True))
    for i, E in enumerate(a.E):
        dt = a.dt0 / max(1.0, E ** 0.25)
        ref_diag = 0.0120 * E ** 2.25
        print(f"--- E = {E:.4g}: diagonal-sector references lambda_diag = {ref_diag:.3e}, "
              f"lambda_perp(vec) ~ {2.2 * ref_diag:.3e} (low E) ; E^(5/4) from E=1 CPU value: {0.0126 * E ** 1.25:.3e}", flush=True)
        r = run_energy(be, E, a.N, dt, a.seed + i, a.target_abs, a.target_rel, a.t_block, a.t_burn_ym,
                       a.tmax, a.wall, a.n_consec, a.renorm_every, a.omega_dt_max, log, cj=a.cj)
        status = 'converged' if r['converged'] else 'CAP HIT'
        print(f"E={E:9.4g} cj={a.cj:g} offdiag-ghost g1={r['g1']:+.3e}+-{r['g1_se']:.1e} (extrap {r['g1_extrap']:+.3e}) "
              f"g2={r['g2']:+.3e}+-{r['g2_se']:.1e} | g1/lam_diag={r['g1'] / ref_diag:.2f} "
              f"T={r['T']:8.0f} {status} ymerr={r['ym_relerr_max']:.1e} | growth in channels (max|q| > 3 E^1/4): "
              f"{100 * r['frac_growth_channel']:.1f}% of log-growth in {100 * r['frac_time_channel']:.1f}% of time "
              f"({r['seconds']:.0f}s, {r['ms_per_step']:.3f} ms/step)", flush=True)
        tag = f'E{E:.4g}_cj{a.cj:g}'
        np.save(os.path.join(a.out, f'lam_i_{tag}.npy'), r['lam_i'])
        np.save(os.path.join(a.out, f'history_{tag}.npy'), r['history'])
        np.savez(os.path.join(a.out, f'landscape_{tag}.npz'), hist_g=r['hist_g'], hist_t=r['hist_t'],
                 D_edges=HIST_EDGES[0], J_edges=HIST_EDGES[1], E=E, cj=a.cj, g1=r['g1'])
        row = dict(E=E, cj=a.cj, N=a.N, backend=be.name, g1=r['g1'], g1_se=r['g1_se'], g1_extrap=r['g1_extrap'],
                   g2=r['g2'], g2_se=r['g2_se'], ref_lam_diag=ref_diag, frac_growth_channel=r['frac_growth_channel'],
                   frac_time_channel=r['frac_time_channel'], T=r['T'], status=status,
                   ym_relerr_max=r['ym_relerr_max'], seconds=r['seconds'], ms_per_step=r['ms_per_step'])
        path = os.path.join(a.out, 'rung1_offdiag.csv')
        new = not os.path.exists(path)
        with open(path, 'a', newline='') as fh:
            wr = csv.DictWriter(fh, fieldnames=list(row.keys()))
            if new:
                wr.writeheader()
            wr.writerow(row)
    print("appended to", os.path.join(a.out, 'rung1_offdiag.csv'))


if __name__ == '__main__':
    main()
