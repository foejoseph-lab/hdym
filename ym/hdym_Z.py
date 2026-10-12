"""hdym_Z.py -- the classical density of states Z of pure Yang-Mills mechanics at E = 1.

    Z = int delta(1 - H) d^3q d^3p,   H = |P|^2/2 + V(q),   V = (x^2 y^2 + y^2 z^2 + z^2 x^2)/2,

the number behind nu0 = 8 pi^2 / Z and the "parameter-free" law
lambda_diag = (32 pi^3 / (105 sqrt 2)) E^{9/4} / Z  (paper/main.tex, sec. "nu0 from the density
of states").  Nothing here assumes a value for Z.  The paper's 557.1 is printed only at the
end, as a comparison, after Z has been computed.

The p-integral is done exactly:  int delta(1 - |P|^2/2 - V) d^3P = 4 pi sqrt(2) sqrt(1 - V),
so  Z = 4 pi sqrt(2) int_{V < 1} sqrt(1 - V) d^3q.  The region V < 1 is unbounded (the channels
along the axes) but the integral converges because the channels narrow as 1/|x|.

Two independent evaluations, which share nothing but V:

  (A) Radial / angular.  V is homogeneous of degree 4: V(r n) = r^4 v(n).  The radial integral is
      int_0^{v^{-1/4}} r^2 sqrt(1 - r^4 v) dr = v^{-3/4} K,   K = int_0^1 s^2 sqrt(1 - s^4) ds
      = B(3/4, 3/2) / 4, so Z = 4 pi sqrt(2) K int_{S^2} v^{-3/4} dOmega.  The angular integrand
      is singular (integrably, like theta^{-1/2}) on the six axis directions; we integrate over
      the 1/48 of the sphere with n_x >= n_y >= n_z >= 0, in polar angles about the x axis, with
      theta = s^2 so the integrand is smooth.

  (B) Cartesian.  Same 1/48 region 0 <= z <= y <= x.  The z-integral is done in closed form,
      int_0^{zm} sqrt(a - b z^2) dz with a = 1 - x^2 y^2 / 2, b = (x^2 + y^2) / 2,
      zm = min(y, sqrt(a/b)); then y on [0, min(x, sqrt(2)/x)] and x on [0, inf) numerically.

Both are adaptive Gauss-Kronrod (scipy.integrate.quad) with tight tolerances; the printed
error estimates are quad's.  Runs in well under a minute.

Usage (from the repository root):
    python ym/hdym_Z.py > ym/data/Z_density_of_states.txt
Requirements: numpy, scipy.
Added 2026-10-12.
"""
import math
import sys
import numpy as np
from scipy import integrate, special

EPS = dict(epsabs=0.0, epsrel=1e-12, limit=400)
EPS_IN = dict(epsabs=1e-15, epsrel=1e-11, limit=400)    # inner y-integral of (B)
EPS_OUT = dict(epsabs=0.0, epsrel=1e-10, limit=400)   # outer x-integral of (B): its integrand is itself a quad


def v_ang(n):
    nx, ny, nz = n
    return 0.5 * (nx * nx * ny * ny + ny * ny * nz * nz + nz * nz * nx * nx)


# ---------------------------------------------------------------- (A) radial / angular
def method_A():
    K_beta = special.beta(0.75, 1.5) / 4.0
    K_quad, K_err = integrate.quad(lambda s: s * s * math.sqrt(max(0.0, 1.0 - s ** 4)), 0.0, 1.0, **EPS)

    # n = (cos th, sin th cos ph, sin th sin ph); n_y >= n_z  <=>  ph in [0, pi/4];
    # n_x >= n_y  <=>  tan th <= 1/cos ph.  dOmega = sin th dth dph;  th = s^2, dth = 2 s ds.
    def inner(ph):
        th_max = math.atan(1.0 / math.cos(ph))
        def f(s):
            th = s * s
            n = (math.cos(th), math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph))
            if s == 0.0:
                return 2.0 ** 1.75                 # limit at the axis, checked below
            return v_ang(n) ** -0.75 * math.sin(th) * 2.0 * s
        val, _ = integrate.quad(f, 0.0, math.sqrt(th_max), **EPS)
        return val

    # check of the s -> 0 limit used above: v -> th^2/2 (n_y^2 + n_z^2 = sin^2 th, n_x -> 1), so
    # v^{-3/4} sin th * 2 s -> (s^4/2)^{-3/4} s^2 * 2 s = 2^{7/4}.  Verified numerically here:
    s_small = 1e-5
    th = s_small ** 2
    n = (math.cos(th), math.sin(th) * math.cos(0.3), math.sin(th) * math.sin(0.3))
    lim_num = v_ang(n) ** -0.75 * math.sin(th) * 2.0 * s_small
    ang48, ang48_err = integrate.quad(inner, 0.0, math.pi / 4, **EPS)
    ang = 48.0 * ang48
    Z = 4.0 * math.pi * math.sqrt(2.0) * K_beta * ang
    return dict(K_beta=K_beta, K_quad=K_quad, K_err=K_err, lim=2.0 ** 1.75, lim_num=lim_num,
                ang=ang, ang_err=48.0 * ang48_err, Z=Z)


# ---------------------------------------------------------------- (B) Cartesian
def zint(x, y, a):
    """int_0^{zm} sqrt(a - b z^2) dz, b = (x^2 + y^2)/2, zm = min(y, sqrt(a/b)); a = 1 - x^2 y^2 / 2
    is passed in, computed without cancellation by the caller."""
    if a <= 0.0:
        return 0.0
    b = 0.5 * (x * x + y * y)
    if y * y * b >= a:                          # the disc edge cuts z before z = y: half-ellipse quarter
        return 0.25 * math.pi * a / math.sqrt(b)
    r = a - b * y * y                           # > 0 here
    return 0.5 * (y * math.sqrt(r) + a / math.sqrt(b) * math.asin(y * math.sqrt(b / a)))


def method_B():
    def over_y(x):
        if x == 0.0:
            return 0.0
        ymax = min(x, math.sqrt(2.0) / x)
        # kink where zm switches from y to sqrt(a/b): y^4/2 + x^2 y^2 - 1 = 0,
        # y^2 = sqrt(x^4 + 2) - x^2 = 2 / (x^2 + sqrt(x^4 + 2))
        yk = math.sqrt(2.0 / (x * x + math.sqrt(x ** 4 + 2.0)))
        # y = ymax t, t in [0, 1];  a = 1 - (x ymax t)^2 / 2, exactly 1 - t^2 = (1-t)(1+t) when ymax = sqrt2/x
        c = x * ymax / math.sqrt(2.0)               # = 1 in the channel range x >= 2^(1/4)
        def g(t):
            ct = c * t
            return ymax * zint(x, ymax * t, (1.0 - ct) * (1.0 + ct))
        pts = [yk / ymax] if 0.0 < yk < ymax else None
        val, _ = integrate.quad(g, 0.0, 1.0, points=pts, **EPS_IN)
        return val

    x_split = 2.0 ** 0.25                       # where x = sqrt(2)/x
    a1, e1 = integrate.quad(over_y, 0.0, x_split, **EPS_OUT)
    a2, e2 = integrate.quad(over_y, x_split, 50.0, **EPS_OUT)
    # x > 50: deep in the channel, z <= y <= sqrt(2)/x << x, so V = x^2 (y^2 + z^2)/2 + O(y^2 z^2)
    # and the yz cross-section of 1/48 is a 1/8 sector of the disc rho < sqrt(2)/x:
    # int sqrt(1 - x^2 rho^2 / 2) rho drho dphi over phi in [0, pi/4] = (pi/4) (2/(3 x^2)) = pi/(6 x^2).
    # Its integral from X to inf is pi/(6 X).  Done numerically too, as a check on the asymptote:
    a3_num, e3 = integrate.quad(over_y, 50.0, np.inf, **EPS_OUT)
    a3_asym = math.pi / (6.0 * 50.0)
    vol48 = a1 + a2 + a3_num
    Z = 4.0 * math.pi * math.sqrt(2.0) * 48.0 * vol48
    return dict(parts=(a1, a2, a3_num), errs=(e1, e2, e3), tail_asym=a3_asym, Z=Z)


def main():
    print(f"hdym_Z.py  python {sys.version.split()[0]}  numpy {np.__version__}  "
          f"scipy {__import__('scipy').__version__}")
    print("Z = int delta(1 - H) d^3q d^3p for V = (x^2y^2 + y^2z^2 + z^2x^2)/2, E = 1\n")

    A = method_A()
    print("(A) radial / angular")
    print(f"    K = B(3/4,3/2)/4 = {A['K_beta']:.15f}   (quad: {A['K_quad']:.15f}, est. err {A['K_err']:.1e})")
    print(f"    integrand limit at the axis: 2^(7/4) = {A['lim']:.10f}   (numerical {A['lim_num']:.10f})")
    print(f"    int_S2 v^(-3/4) dOmega = {A['ang']:.12f}   (est. err {A['ang_err']:.1e})")
    print(f"    Z_A = {A['Z']:.10f}\n")

    B = method_B()
    print("(B) Cartesian, closed-form z integral")
    print("    1/48 volume integral, x in [0, 2^(1/4)], [2^(1/4), 50], [50, inf):")
    for p, e in zip(B['parts'], B['errs']):
        print(f"        {p:.15f}   (est. err {e:.1e})")
    print(f"    tail x > 50 by the channel asymptote pi/(6*50) = {B['tail_asym']:.15f}  (vs numerical above)")
    print(f"    Z_B = {B['Z']:.10f}\n")

    Z = A['Z']
    print(f"relative difference (A - B)/A = {(A['Z'] - B['Z']) / A['Z']:.2e}")
    print(f"\nRESULT: Z = {Z:.6f}")
    nu0 = 8 * math.pi ** 2 / Z
    C = 32 * math.pi ** 3 / (105 * math.sqrt(2))
    print(f"        nu0 = 8 pi^2 / Z = {nu0:.6f}")
    print(f"        lambda_diag coefficient (32 pi^3 / 105 sqrt 2) / Z = {C / Z:.6f}   (times E^(9/4))")
    print(f"        c_j = 0 off-diagonal coefficient (4 sqrt2 pi^3 / 15) / Z = {4 * math.sqrt(2) * math.pi ** 3 / 15 / Z:.6f}   (times E^(5/4))")

    print("\ncomparison with the values quoted in the repository (not used above):")
    print(f"    paper/main.tex Z = 557.1         -> computed rounds to {Z:.1f}")
    print(f"    paper/main.tex nu0 = 0.1417      -> computed {nu0:.4f}")
    print(f"    paper/main.tex 0.0120 E^(9/4)    -> computed {C / Z:.4f}")
    print(f"    ym/notes/NOTE_rung1.md 0.0210 E^(5/4) -> computed {4 * math.sqrt(2) * math.pi ** 3 / 15 / Z:.4f}")


if __name__ == "__main__":
    main()
