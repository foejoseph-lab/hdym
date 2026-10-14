#!/usr/bin/env python3
"""
arb_smoke.py -- smoke test: can Arb (FLINT) carry the ball kernel's centre arithmetic, from
Python and from C, on this machine, and what would it cost?  Scratch tool; changes nothing.

    python ym/tools/arb_smoke.py       (any folder; also %runfile in Spyder)

Needs:  pip install python-flint   (Python side)
        MSYS2:  pacman -S mingw-w64-x86_64-flint   (C side; Linux: libflint-dev)
Each part is skipped, with the reason printed, if its library is missing.

What it does: one Yoshida-4 step (dt = 0.005, one sub-step) of the 12-d Lee-Wick system, the
nonlinear kick of hdym_ball.py, from an exact double initial state, at 200 bits:
  A. in python-flint `arb`;  B. in FLINT's C library through ctypes (same algorithm text);
  R. as reference, in mpmath at 120 digits with the same double stage constants fl(YC*h).
Checks: R lies inside every A ball and every B ball; A and B midpoints agree; prints the
radii (the one-step rounding defect) and the time per step in each, and projects the cost of
rerunning one long-lived pair (300 000 steps).
"""
import ctypes, os, shutil, subprocess, sys, tempfile, time

YC = (0.6756035959798288, -0.17560359597982886, -0.17560359597982886, 0.6756035959798288)
YD = (1.3512071919596576, -1.7024143839193153, 1.3512071919596576)
DT, PREC = 0.005, 200
X0 = [3.0e-1, -2.0e-1, 1.0e-1, 1.0e-2, -2.0e-2, 5.0e-3,      # p, w
      4.0e-1, 3.0e-1, -5.0e-1, -1.0e-2, 2.0e-2, 1.0e-2]      # P, W   (exact doubles)


def step_generic(z, num, h):
    """One Yoshida step; `num` converts a Python float to the number type."""
    p, w, P, W = z[0:3], z[3:6], z[6:9], z[9:12]

    def kick(d):
        x, y, zz = p[0] - w[0], p[1] - w[1], p[2] - w[2]
        x2, y2, z2 = x * x, y * y, zz * zz
        yz, zx, xy = y2 + z2, z2 + x2, x2 + y2
        g = (x * yz, y * zx, zz * xy)
        txy, txz, tyz = (2 * x) * y, (2 * x) * zz, (2 * y) * zz
        hv = (yz * w[0] + txy * w[1] + txz * w[2],
              txy * w[0] + zx * w[1] + tyz * w[2],
              txz * w[0] + tyz * w[1] + xy * w[2])
        hh = num(d)
        for k in range(3):
            P[k] = P[k] - hh * (g[k] + hv[k])
            W[k] = W[k] + hh * (hv[k] + w[k])

    def drift(c):
        cc = num(c)
        for k in range(3):
            p[k] = p[k] + cc * P[k]
            w[k] = w[k] - cc * W[k]
    drift(YC[0] * h); kick(YD[0] * h); drift(YC[1] * h); kick(YD[1] * h)
    drift(YC[2] * h); kick(YD[2] * h); drift(YC[3] * h)
    return p + w + P + W


# ---------------------------------------------------------------- R: mpmath reference
import mpmath as mp
mp.mp.dps = 120
ref = step_generic([mp.mpf(v) for v in X0], mp.mpf, DT)

# ---------------------------------------------------------------- A: python-flint
A = None
try:
    import flint
    flint.ctx.prec = PREC
    t0 = time.perf_counter()
    reps = 200
    for _ in range(reps):
        A = step_generic([flint.arb(v) for v in X0], flint.arb, DT)
    tA = (time.perf_counter() - t0) / reps
    insideA = all(a.contains(flint.arb(mp.nstr(r, 110))) for a, r in zip(A, ref))
    radA = max(float(a.rad()) for a in A)
    print(f"A python-flint {flint.__version__}: reference inside every ball: {insideA};  "
          f"max radius {radA:.2e};  {tA * 1e6:.0f} us per step")
except ImportError:
    print("A skipped: python-flint not installed (pip install python-flint)")

# ---------------------------------------------------------------- B: FLINT C library
C_SRC = r"""
#include <flint/arb.h>
#define N 12
static const double YC[4] = {0.6756035959798288, -0.17560359597982886, -0.17560359597982886, 0.6756035959798288};
static const double YD[3] = {1.3512071919596576, -1.7024143839193153, 1.3512071919596576};
static void kick(arb_ptr z, double d, slong prec, arb_ptr t)
{   /* t: 20 scratch balls */
    arb_ptr p = z, w = z + 3, P = z + 6, W = z + 9;
    arb_ptr x = t, y = t+1, q = t+2, x2 = t+3, y2 = t+4, z2 = t+5, yz = t+6, zx = t+7, xy = t+8,
            txy = t+9, txz = t+10, tyz = t+11, g = t+12, hv = t+15, hh = t+18, tmp = t+19;
    int k;
    arb_sub(x, p, w, prec); arb_sub(y, p+1, w+1, prec); arb_sub(q, p+2, w+2, prec);
    arb_mul(x2, x, x, prec); arb_mul(y2, y, y, prec); arb_mul(z2, q, q, prec);
    arb_add(yz, y2, z2, prec); arb_add(zx, z2, x2, prec); arb_add(xy, x2, y2, prec);
    arb_mul(g, x, yz, prec); arb_mul(g+1, y, zx, prec); arb_mul(g+2, q, xy, prec);
    arb_mul_2exp_si(tmp, x, 1); arb_mul(txy, tmp, y, prec); arb_mul(txz, tmp, q, prec);
    arb_mul_2exp_si(tmp, y, 1); arb_mul(tyz, tmp, q, prec);
    /* hv = H w, same operation order as hdym_ball's kick */
    arb_mul(hv, yz, w, prec);   arb_addmul(hv, txy, w+1, prec); arb_addmul(hv, txz, w+2, prec);
    arb_mul(hv+1, txy, w, prec); arb_addmul(hv+1, zx, w+1, prec); arb_addmul(hv+1, tyz, w+2, prec);
    arb_mul(hv+2, txz, w, prec); arb_addmul(hv+2, tyz, w+1, prec); arb_addmul(hv+2, xy, w+2, prec);
    arb_set_d(hh, d);
    for (k = 0; k < 3; k++) {
        arb_add(tmp, g+k, hv+k, prec); arb_submul(P+k, hh, tmp, prec);
        arb_add(tmp, hv+k, w+k, prec); arb_addmul(W+k, hh, tmp, prec);
    }
}
static void drift(arb_ptr z, double c, slong prec, arb_ptr t)
{
    int k;
    arb_set_d(t, c);
    for (k = 0; k < 3; k++) { arb_addmul(z+k, t, z+6+k, prec); arb_submul(z+3+k, t, z+9+k, prec); }
}
void step(const double* x0, double h, long prec, int reps, double* mid, double* rad)
{
    arb_ptr z = _arb_vec_init(N), t = _arb_vec_init(20);
    int k, r;
    for (r = 0; r < reps; r++) {
        for (k = 0; k < N; k++) arb_set_d(z + k, x0[k]);
        drift(z, YC[0]*h, prec, t); kick(z, YD[0]*h, prec, t); drift(z, YC[1]*h, prec, t);
        kick(z, YD[1]*h, prec, t); drift(z, YC[2]*h, prec, t); kick(z, YD[2]*h, prec, t);
        drift(z, YC[3]*h, prec, t);
    }
    for (k = 0; k < N; k++) {
        mid[k] = arf_get_d(arb_midref(z + k), ARF_RND_NEAR);
        rad[k] = mag_get_d(arb_radref(z + k));
    }
    _arb_vec_clear(z, N); _arb_vec_clear(t, 20);
}
"""
cc = os.environ.get('CC', 'gcc')
exe = shutil.which(cc)
if not exe:
    print("B skipped: no C compiler on the PATH")
else:
    d = tempfile.mkdtemp()
    src, lib = os.path.join(d, 'arbstep.c'), os.path.join(d, 'arbstep' + ('.dll' if sys.platform == 'win32' else '.so'))
    open(src, 'w').write(C_SRC)
    bindir = os.path.dirname(os.path.abspath(exe))
    env = os.environ.copy()
    env['PATH'] = bindir + os.pathsep + env.get('PATH', '')       # as ym/cbuild.py (Spyder DLLs)
    cmd = [cc, '-O2', '-ffp-contract=off', '-shared', '-fPIC', '-o', lib, src, '-lflint', '-lmpfr', '-lgmp']
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"B skipped: compiling against FLINT failed (exit {r.returncode}); install it with "
              f"'pacman -S mingw-w64-x86_64-flint' in MSYS2.\n{r.stdout}{r.stderr}".rstrip())
    else:
        if sys.platform == 'win32':
            os.add_dll_directory(bindir)
        L = ctypes.CDLL(lib)
        dp = ctypes.POINTER(ctypes.c_double)
        L.step.argtypes = [dp, ctypes.c_double, ctypes.c_long, ctypes.c_int, dp, dp]
        x0 = (ctypes.c_double * 12)(*X0); mid = (ctypes.c_double * 12)(); rad = (ctypes.c_double * 12)()
        L.step(x0, DT, PREC, 1, mid, rad)
        reps = 20000
        t0 = time.perf_counter(); L.step(x0, DT, PREC, reps, mid, rad); tB = (time.perf_counter() - t0) / reps
        insideB = all(abs(mp.mpf(m) - r) <= mp.mpf(rd) + mp.mpf(2) ** -52 * abs(mp.mpf(m)) for m, rd, r in zip(mid, rad, ref))
        print(f"B FLINT C library: reference inside every ball (midpoint read as double, +1 ulp): {insideB};  "
              f"max radius {max(rad):.2e};  {tB * 1e6:.1f} us per step")
        if A is not None:
            agree = all(float(a.mid()) == m for a, m in zip(A, mid))
            print(f"A vs B: midpoints (as doubles) identical: {agree}")
            print(f"cost of rerunning one long-lived pair (300 000 steps, 2 members, ~3 sub-steps/step): "
                  f"python-flint ~{tA * 300000 * 6 / 60:.1f} min,  C ~{tB * 300000 * 6 / 60:.1f} min "
                  f"(centre arithmetic only; the frame stays in C either way)")
