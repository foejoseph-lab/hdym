#!/usr/bin/env python3
"""
hdym_escape_map.py -- escape-time landscape t_esc(initial condition) on a 2D
section of the 12-dimensional Lee-Wick Yang-Mills phase space (units g = M = 1).

Sections (fixed energy E and ghost fraction eps):

  --section ghost   Yang-Mills point (q0, P0) fixed; ghost orientation varied:
                    w = a (cos th e1 + sin th e2),  W = a (cos ph e1 + sin ph e2),
                    a = sqrt(eps E)  so that |w|^2 + |W|^2 = 2 eps E as in
                    hdym_core.full_initial.  Axes: th, ph in [0, 2 pi).
  --section ym      Ghost (direction and size) fixed; the YM point varied over a
                    patch  q = q0 + (u e1 + v e2),  P = s P0  with s chosen so
                    that 1/2|P|^2 + V(q) = E.  Cells with V(q) >= E are NaN.
                    Axes: u, v in [-width/2, width/2] around --center.

Both sections use the same base YM point for a given --seed, so ghost/ym maps
at the same (E, eps, seed) are two slices through the same neighbourhood.

Zoom: --center cu cv --width w re-grids a sub-region with the same base point,
same dt, same R (the ordering of escape times is the observable, so keep those
fixed between a map and its zoom).

Backends
  --backend cupy   CUDA RawKernel (RTX 4080: ~1 h for 1024^2 at E = 0.3)
  --backend c      same C source via gcc + OpenMP (256^2 in minutes on 12 cores)
  --backend numpy  hdym_core.evolve_ensemble, no sub-stepping, slow; checks only

Outputs <out>/<tag>.npz (t_esc grid, axes, all parameters, base point) and
<out>/<tag>.png (log10 t_esc).

Examples
  python hdym_escape_map.py --E 0.3 --eps 0.01 --section ym    --n 256  --backend c
  python hdym_escape_map.py --E 0.3 --eps 0.01 --section ghost --n 256  --backend c
  python hdym_escape_map.py --E 0.3 --eps 0.01 --section ym    --n 1024 --backend cupy
  python hdym_escape_map.py --E 0.3 --eps 0.01 --section ym    --n 1024 --backend cupy \
         --center 0.12 -0.05 --width 0.04 --tag zoom10x
"""
import argparse
import ctypes
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

# ---------------------------------------------------------------- kernel (C / CUDA)
STEP_SRC = r"""
#define YC0 0.6756035959798288
#define YC1 -0.17560359597982886
#define YC2 -0.17560359597982886
#define YC3 0.6756035959798288
#define YD0 1.3512071919596576
#define YD1 -1.7024143839193153
#define YD2 1.3512071919596576
#define MAXSUB 16384

/* full nonlinear kick:  q = p - w,
   Pdot = -(grad V(q) + Hess V(q) w),   Wdot = +(Hess V(q) w + w)            */
DEVFN void kick_nl(const double* p, const double* w, double* P, double* W, double h)
{
    double x = p[0] - w[0], y = p[1] - w[1], z = p[2] - w[2];
    double x2 = x * x, y2 = y * y, z2 = z * z;
    double gx = x * (y2 + z2), gy = y * (z2 + x2), gz = z * (x2 + y2);
    double wx = w[0], wy = w[1], wz = w[2];
    double hx = (y2 + z2) * wx + 2.0 * x * y * wy + 2.0 * x * z * wz;
    double hy = 2.0 * x * y * wx + (z2 + x2) * wy + 2.0 * y * z * wz;
    double hz = 2.0 * x * z * wx + 2.0 * y * z * wy + (x2 + y2) * wz;
    P[0] -= h * (gx + hx);  P[1] -= h * (gy + hy);  P[2] -= h * (gz + hz);
    W[0] += h * (hx + wx);  W[1] += h * (hy + wy);  W[2] += h * (hz + wz);
}

/* Yoshida-4 for H = 1/2|P|^2 - 1/2|W|^2 + U(p,w):  pdot = P, wdot = -W     */
DEVFN void yoshida_nl(double* p, double* w, double* P, double* W, double h)
{
    int k;
    for (k = 0; k < 3; k++) { p[k] += YC0 * h * P[k]; w[k] -= YC0 * h * W[k]; }
    kick_nl(p, w, P, W, YD0 * h);
    for (k = 0; k < 3; k++) { p[k] += YC1 * h * P[k]; w[k] -= YC1 * h * W[k]; }
    kick_nl(p, w, P, W, YD1 * h);
    for (k = 0; k < 3; k++) { p[k] += YC2 * h * P[k]; w[k] -= YC2 * h * W[k]; }
    kick_nl(p, w, P, W, YD2 * h);
    for (k = 0; k < 3; k++) { p[k] += YC3 * h * P[k]; w[k] -= YC3 * h * W[k]; }
}

/* integrate trajectory i until max|component| > R (or non-finite) or nmax
   steps; tesc[i] = escape time, or -1 if it did not escape by nmax*dt.
   Sub-stepping uses the local frequency scale max|q| (channel depth), which
   bounds both the normal-sector and the ghost-sector frequencies.            */
DEVFN void escape_traj(double* S, double* tesc, int N, int i, double dt, long nmax,
                       double R, double omega_dt_max)
{
    double p[3], w[3], P[3], W[3];
    int k, n, sub;
    long s;
    for (k = 0; k < 3; k++) {
        p[k] = S[(0 + k) * N + i];
        w[k] = S[(3 + k) * N + i];
        P[k] = S[(6 + k) * N + i];
        W[k] = S[(9 + k) * N + i];
    }
    double t = -1.0;
    for (s = 1; s <= nmax; s++) {
        double om = 0.0, a;
        for (k = 0; k < 3; k++) { a = fabs(p[k] - w[k]); if (a > om) om = a; }
        n = (int)ceil(om * dt / omega_dt_max);
        if (n < 1) n = 1;
        if (n > MAXSUB) n = MAXSUB;
        if (n == 1) {
            yoshida_nl(p, w, P, W, dt);
        } else {
            double h = dt / n;
            for (sub = 0; sub < n; sub++) yoshida_nl(p, w, P, W, h);
        }
        double m = 0.0;
        int bad = 0;
        for (k = 0; k < 3; k++) {
            a = fabs(p[k]); if (a > m) m = a;
            a = fabs(w[k]); if (a > m) m = a;
            a = fabs(P[k]); if (a > m) m = a;
            a = fabs(W[k]); if (a > m) m = a;
        }
        if (!(m == m) || m > 1e300) bad = 1;      /* NaN / overflow */
        if (bad || m > R) { t = s * dt; break; }
    }
    tesc[i] = t;
    for (k = 0; k < 3; k++) {
        S[(0 + k) * N + i] = p[k];
        S[(3 + k) * N + i] = w[k];
        S[(6 + k) * N + i] = P[k];
        S[(9 + k) * N + i] = W[k];
    }
}
"""

CUDA_SRC = ("#define DEVFN __device__ __forceinline__\n" + STEP_SRC + r"""
extern "C" __global__
void escape(double* S, double* tesc, int N, double dt, long nmax, double R, double omega_dt_max)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) escape_traj(S, tesc, N, i, dt, nmax, R, omega_dt_max);
}
""")

C_SRC = ("#include <math.h>\n#define DEVFN static inline\n" + STEP_SRC + r"""
void escape(double* S, double* tesc, int N, double dt, long nmax, double R, double omega_dt_max)
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (i = 0; i < N; i++)
        escape_traj(S, tesc, N, i, dt, nmax, R, omega_dt_max);
}
""")


class CupyBackend:
    name = 'cupy'

    def __init__(self):
        import cupy as cp
        self.cp = cp
        self.kern = cp.RawKernel(CUDA_SRC, 'escape')
        dev = cp.cuda.Device()
        free, total = dev.mem_info
        print(f"CuPy {cp.__version__} on device {dev.id}, {free / 1e9:.1f}/{total / 1e9:.1f} GB free",
              flush=True)

    def escape(self, S_host, dt, nmax, R, omega_dt_max):
        cp = self.cp
        S = cp.asarray(np.ascontiguousarray(S_host))
        N = S.shape[1]
        tesc = cp.full(N, -1.0)
        threads = 128
        blocks = (N + threads - 1) // threads
        self.kern((blocks,), (threads,),
                  (S, tesc, np.int32(N), np.float64(dt), np.int64(nmax), np.float64(R),
                   np.float64(omega_dt_max)))
        cp.cuda.Stream.null.synchronize()
        return cp.asnumpy(tesc)


class CBackend:
    name = 'c'

    def __init__(self):
        d = tempfile.mkdtemp()
        src = os.path.join(d, 'hdym_escape.c')
        lib = os.path.join(d, 'hdym_escape' + ('.dll' if sys.platform == 'win32' else '.so'))
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
        self.lib.escape.argtypes = [dp, dp, ctypes.c_int, ctypes.c_double, ctypes.c_long,
                                    ctypes.c_double, ctypes.c_double]
        print(f"C backend compiled with {cc}, {omp}", flush=True)

    def escape(self, S_host, dt, nmax, R, omega_dt_max):
        S = np.ascontiguousarray(S_host, dtype=np.float64)
        N = S.shape[1]
        tesc = np.full(N, -1.0)
        dp = ctypes.POINTER(ctypes.c_double)
        self.lib.escape(S.ctypes.data_as(dp), tesc.ctypes.data_as(dp), N, dt, int(nmax), R,
                        omega_dt_max)
        return tesc


class NumpyBackend:
    """hdym_core.evolve_ensemble (no sub-stepping) -- for checking only."""
    name = 'numpy'

    def __init__(self):
        import hdym_core
        self.core = hdym_core
        print("numpy backend (hdym_core.evolve_ensemble; no channel sub-stepping)", flush=True)

    def escape(self, S_host, dt, nmax, R, omega_dt_max):
        s0 = tuple(np.ascontiguousarray(S_host[3 * k:3 * k + 3].T) for k in range(4))
        r = self.core.evolve_ensemble(s0, dt, nmax * dt, R, n_samples=10)
        t = r['t_escape']
        t = np.where(np.isfinite(t), t, -1.0)
        return t


def integrate_escape(be, S0, dt, tmax, R, omega_dt_max=0.25, chunk=None, log=None):
    """S0: (12, N) initial states.  Returns t_esc (N,), NaN where not escaped by tmax."""
    N = S0.shape[1]
    nmax = int(np.ceil(tmax / dt))
    chunk = N if chunk is None else chunk
    out = np.empty(N)
    t0 = time.time()
    for a in range(0, N, chunk):
        b = min(N, a + chunk)
        out[a:b] = be.escape(S0[:, a:b], dt, nmax, R, omega_dt_max)
        if log:
            el = time.time() - t0
            log(f"    {b}/{N} trajectories, {el:.0f}s elapsed, ~{el / b * (N - b):.0f}s left")
    out[out < 0] = np.nan
    return out


# ---------------------------------------------------------------- sections
def V1(q):
    x, y, z = q
    return 0.5 * (x * x * y * y + y * y * z * z + z * z * x * x)


def base_point(E, eps, seed):
    """One YM phase point at energy E and one ghost direction, from the seed."""
    rng = np.random.default_rng(seed)
    q = rng.standard_normal(3)
    P = rng.standard_normal(3)
    E0 = 0.5 * P @ P + V1(q)
    lam = (E / E0) ** 0.25
    q, P = q * lam, P * lam ** 2
    nw = rng.standard_normal(3)
    nW = rng.standard_normal(3)
    nrm = np.sqrt(nw @ nw + nW @ nW)
    amp = np.sqrt(2.0 * eps * E) / nrm
    return q, P, nw * amp, nW * amp


def unit_pair(plane):
    e = {'x': np.array([1., 0, 0]), 'y': np.array([0, 1., 0]), 'z': np.array([0, 0, 1.])}
    return e[plane[0]], e[plane[1]]


def section_ghost(E, eps, seed, n, plane, center, width):
    q0, P0, _, _ = base_point(E, eps, seed)
    e1, e2 = unit_pair(plane)
    a = np.sqrt(eps * E)
    cu, cv = center
    th = cu + np.linspace(-width / 2, width / 2, n, endpoint=False)
    ph = cv + np.linspace(-width / 2, width / 2, n, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing='xy')          # PH rows, TH columns
    w = a * (np.cos(TH)[..., None] * e1 + np.sin(TH)[..., None] * e2)
    W = a * (np.cos(PH)[..., None] * e1 + np.sin(PH)[..., None] * e2)
    q = np.broadcast_to(q0, w.shape)
    Pq = np.broadcast_to(P0, w.shape)
    p = q + w
    P = Pq - W
    S = np.concatenate([p.reshape(-1, 3).T, w.reshape(-1, 3).T,
                        P.reshape(-1, 3).T, W.reshape(-1, 3).T], 0)
    return S, th, ph, np.ones((n, n), bool), dict(q0=q0, P0=P0)


def section_ym(E, eps, seed, n, plane, center, width):
    q0, P0, w0, W0 = base_point(E, eps, seed)
    e1, e2 = unit_pair(plane)
    cu, cv = center
    u = cu + np.linspace(-width / 2, width / 2, n, endpoint=False)
    v = cv + np.linspace(-width / 2, width / 2, n, endpoint=False)
    U, Vv = np.meshgrid(u, v, indexing='xy')
    q = q0 + U[..., None] * e1 + Vv[..., None] * e2
    Vq = V1(q.reshape(-1, 3).T).reshape(n, n)
    K0 = 0.5 * P0 @ P0
    ok = (E - Vq) > 0
    s = np.sqrt(np.clip(E - Vq, 0, None) / K0)
    Pq = s[..., None] * P0
    w = np.broadcast_to(w0, q.shape)
    W = np.broadcast_to(W0, q.shape)
    p = q + w
    P = Pq - W
    S = np.concatenate([p.reshape(-1, 3).T, w.reshape(-1, 3).T,
                        P.reshape(-1, 3).T, W.reshape(-1, 3).T], 0)
    return S, u, v, ok, dict(q0=q0, P0=P0, w0=w0, W0=W0)


def make_backend(name):
    return {'cupy': CupyBackend, 'c': CBackend, 'numpy': NumpyBackend}[name]()


def default_width(section, E):
    # ghost: full torus.  ym: a patch of size ~ q_typ ~ E^{1/4} (whole neighbourhood);
    # zoom with --width.
    return 2 * np.pi if section == 'ghost' else 1.0 * E ** 0.25


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--E', type=float, required=True)
    ap.add_argument('--eps', type=float, required=True)
    ap.add_argument('--section', choices=['ghost', 'ym'], default='ym')
    ap.add_argument('--n', type=int, default=256, help='grid is n x n')
    ap.add_argument('--plane', default='xy', help='two of x,y,z: the section directions')
    ap.add_argument('--center', type=float, nargs=2, default=None,
                    help='section centre (ghost: th ph; ym: u v).  Default: pi pi / 0 0')
    ap.add_argument('--width', type=float, default=None,
                    help='side of the square section (ghost default 2pi, ym default E^1/4)')
    ap.add_argument('--dt', type=float, default=0.005)
    ap.add_argument('--tmax', type=float, default=None,
                    help='stop integrating at tmax (default 40/lam_perp with lam ~ 0.027 E^2.24, '
                         'capped at 5000); unescaped cells are NaN')
    ap.add_argument('--R', type=float, default=1e4, help='escape radius (max |component|)')
    ap.add_argument('--omega-dt-max', type=float, default=0.25)
    ap.add_argument('--backend', choices=['cupy', 'c', 'numpy'], default='c')
    ap.add_argument('--chunk', type=int, default=None, help='trajectories per kernel call')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--tag', default=None)
    ap.add_argument('--out', default='results_escape')
    a = ap.parse_args()

    if a.center is None:
        a.center = [np.pi, np.pi] if a.section == 'ghost' else [0.0, 0.0]
    if a.width is None:
        a.width = default_width(a.section, a.E)
    if a.tmax is None:
        lam_guess = 0.027 * a.E ** 2.24
        a.tmax = min(5000.0, 40.0 / lam_guess)
    tag = a.tag or f"{a.section}_E{a.E:g}_eps{a.eps:g}_n{a.n}_s{a.seed}"
    os.makedirs(a.out, exist_ok=True)

    sec = section_ghost if a.section == 'ghost' else section_ym
    S, ax1, ax2, ok, base = sec(a.E, a.eps, a.seed, a.n, a.plane, a.center, a.width)
    N = S.shape[1]
    print(f"{tag}: {N} initial conditions, {ok.sum()} on-shell, dt={a.dt}, tmax={a.tmax:.0f}, "
          f"R={a.R:g}", flush=True)

    be = make_backend(a.backend)
    chunk = a.chunk or (N if be.name == 'cupy' else max(a.n * 8, 4096))
    t0 = time.time()
    tesc = integrate_escape(be, S, a.dt, a.tmax, a.R, a.omega_dt_max, chunk=chunk,
                            log=lambda m: print(m, flush=True))
    tesc = tesc.reshape(a.n, a.n)
    tesc[~ok] = np.nan
    el = time.time() - t0
    fin = np.isfinite(tesc)
    print(f"done in {el:.0f}s ({el / max(1, N) * 1e3:.2f} ms/trajectory).  escaped: "
          f"{fin.mean():.3f}; t_esc median {np.nanmedian(tesc):.1f}, "
          f"p5 {np.nanpercentile(tesc, 5):.1f}, p95 {np.nanpercentile(tesc, 95):.1f}", flush=True)

    npz = os.path.join(a.out, tag + '.npz')
    np.savez(npz, t_esc=tesc, ax1=ax1, ax2=ax2, ok=ok, E=a.E, eps=a.eps, section=a.section,
             plane=a.plane, center=np.array(a.center), width=a.width, dt=a.dt, tmax=a.tmax,
             R=a.R, seed=a.seed, n=a.n, backend=be.name, seconds=el, **base)
    print("wrote", npz)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 6))
        im = ax.imshow(np.log10(tesc), origin='lower', cmap='magma',
                       extent=[ax1[0], ax1[-1], ax2[0], ax2[-1]], aspect='auto',
                       interpolation='nearest')
        lab = ('theta', 'phi') if a.section == 'ghost' else (f'u ({a.plane[0]})', f'v ({a.plane[1]})')
        ax.set_xlabel(lab[0]); ax.set_ylabel(lab[1])
        ax.set_title(f"log10 t_esc   {a.section} section   E={a.E:g} eps={a.eps:g}  {a.n}^2")
        plt.colorbar(im, ax=ax, label='log10 t_esc')
        plt.tight_layout()
        png = os.path.join(a.out, tag + '.png')
        plt.savefig(png, dpi=150)
        print("wrote", png)
    except Exception as e:                       # noqa
        print("plot skipped:", e)


if __name__ == '__main__':
    main()
