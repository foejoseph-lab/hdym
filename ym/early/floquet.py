"""Floquet stability of the linearised ghost on the diagonal periodic orbit x=y, z=0 of YM mechanics.
On it V = x^4/2 (with x=y), so xdd = -2x^3 ... check: V = 1/2 (x^2y^2) at z=0 -> dV/dx = x y^2 = x^3 along x=y.
Ghost modes: constant eigenvectors of Hess V = [[y^2,2xy,0],[2xy,x^2,0],[0,0,x^2+y^2]] at x=y:
  (1,1,0)/sqrt2 -> 3x^2 ; (1,-1,0)/sqrt2 -> -x^2 ; (0,0,1) -> 2x^2.
Monodromy of wdd = -(1+K(t)) w over one period; growth rate = ln(max|mu|)/period.
"""
import numpy as np
from scipy.integrate import solve_ivp
def orbit_period(E):
    # 1/2 xdot^2 + 1/2 x^4 = E  (energy of the (x,y) motion: 1/2(xd^2+yd^2) + 1/2 x^2 y^2 = xd^2 + x^4/2 with y=x)
    # along x=y: E = xdot^2 + x^4/2  -> amplitude A=(2E)^{1/4}; period = 4 * int_0^A dx / sqrt(E - x^4/2)
    A = (2*E)**0.25
    from scipy.integrate import quad
    return 4*quad(lambda x: 1/np.sqrt(E - x**4/2), 0, A)[0], A
def growth(E, coef):
    P, A = orbit_period(E)
    def rhs(t, s):
        x, xd, w1, w1d, w2, w2d = s
        K = coef*x*x
        return [xd, -x**3, w1d, -(1+K)*w1, w2d, -(1+K)*w2]
    sol = solve_ivp(rhs, [0, P], [A, 0, 1, 0, 0, 1], rtol=1e-10, atol=1e-12)
    x, xd, w1, w1d, w2, w2d = sol.y[:, -1]
    M = np.array([[w1, w2], [w1d, w2d]])
    tr = np.trace(M); mu = np.abs(np.linalg.eigvals(M)).max()
    return tr, np.log(mu)/P, abs(x-A)+abs(xd)
Es = np.geomspace(0.05, 10, 600)
out = {}
for name, coef in [('diag', 3), ('ridge', -1), ('z', 2)]:
    rows = [growth(E, coef) for E in Es]
    tr = np.array([r[0] for r in rows]); g = np.array([r[1] for r in rows]); err = max(r[2] for r in rows)
    out[name] = (tr, g)
    unstable = np.abs(tr) > 2 + 1e-9
    # report tongues
    edges = np.flatnonzero(np.diff(unstable.astype(int)))
    tongues = []
    i = 0
    while i < len(Es):
        if unstable[i]:
            j = i
            while j+1 < len(Es) and unstable[j+1]: j += 1
            tongues.append((Es[i], Es[j] if j+1 < len(Es) else np.inf, g[i:j+1].max()))
            i = j+1
        else: i += 1
    print(f"{name} mode (K={coef}x^2): orbit closure err {err:.1e}; unstable intervals (E_lo, E_hi, max rate):")
    for t in tongues: print("   %.3f  %.3f   %.3f" % t)
    k = np.argmin(abs(Es-0.3)); print(f"   rate at E=0.3: {g[k]:.3f}")
np.savez('floquet.npz', E=Es, **{f'{n}_tr': v[0] for n, v in out.items()}, **{f'{n}_g': v[1] for n, v in out.items()})
