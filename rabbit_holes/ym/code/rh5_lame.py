"""rh5_lame.py -- the diagonal ghost on the x = y = z orbit is a finite-gap (Lame n = 2) Hill equation.

Background (RH5_nonplanar.md): on x = y = z = s(t) the motion is s'' = -2 s^3, E = 3 s'^2/2 + 3 s^4/2, solved by
    s(t) = A cn(sqrt2 A t | k^2 = 1/2),   E = 3 A^4 / 2.
Hess V there is 2 s^2 J (J = all-ones 3x3): eigenvalues 6 s^2 (isotropic, (1,1,1)) and 0, 0 (anisotropic). The
diagonal ghost on a ghost-free background obeys w'' = -(1 + Hess V) w (units g = M = 1; B003), so it splits into
    isotropic:   w'' + (1 + 6 s^2) w = 0        anisotropic (x2):   w'' + w = 0  (always stable).
With u = sqrt2 A t and cn^2 = 1 - sn^2 the isotropic mode is the Lame equation
    w_uu + (h - n(n+1) k^2 sn^2 u) w = 0,   n(n+1) k^2 = 3,  k^2 = 1/2  ->  n = 2,   h = 3 + 1/(2 A^2).
n = 2 Lame is a finite-gap potential: its band edges are the five Lame-polynomial eigenvalues
    1 + k^2, 1 + 4k^2, 4 + k^2, 2(1 + k^2) -+ 2 sqrt(1 - k^2 + k^4)        (k^2 = 1/2: 1.268, 1.5, 3, 4.5, 4.732)
with gaps (1.5, 3) and (4.5, 4.732) and no instability above 4.732. Since h > 3 for every E, only the second gap
is reachable: exactly ONE instability tongue in all of E, at h in (9/2, 3 + sqrt3), i.e. E in (1/8, 1/6) exactly.
(The unreachable gap ends at h = 3, the ghost-free limit E -> infinity: there the isotropic mode is the orbit's own
tangent s'(t), a band-edge solution.)
Prediction checked here three ways: (1) the closed-form edges against a numerical Lame discriminant; (2) a direct
Floquet scan of the full 3x3 diagonal ghost along the integrated orbit (no eigen-split assumed) over E in [0.02, 3]
(the RH5 range); (3) the tongue's edges by bisection, against the closed form. Also supplies the committed log for
RH5's "|mu| = 1.22 at E = 0.13-0.16", which was stated in the note without one.

Standalone (numpy, scipy). Usage: python rabbit_holes/ym/code/rh5_lame.py > rabbit_holes/ym/results/rh5_lame.txt
Cost: under a minute. Added 2026-10-26.
"""
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipj, ellipk

T0 = time.time()
K2 = 0.5
OK = True


def rep(name, err, tol):
    global OK
    good = bool(err <= tol)
    OK &= good
    print(f"[{'OK ' if good else 'FAIL'}] {name}: {err:.2e}", flush=True)


# ---------------------------------------------------------------- closed form
edges = np.sort([1 + K2, 1 + 4 * K2, 4 + K2, 2 * (1 + K2) - 2 * np.sqrt(1 - K2 + K2 ** 2),
                 2 * (1 + K2) + 2 * np.sqrt(1 - K2 + K2 ** 2)])
gap2 = (4 + K2, 2 * (1 + K2) + 2 * np.sqrt(1 - K2 + K2 ** 2))       # (4.5, 3 + sqrt3)


def E_of_h(h):
    """h = 3 + 1/(2 A^2), E = 3 A^4 / 2."""
    A2 = 1.0 / (2.0 * (h - 3.0))
    return 1.5 * A2 ** 2


E_lo, E_hi = sorted((E_of_h(gap2[1]), E_of_h(gap2[0])))
print("rh5_lame.py: the diagonal ghost on x = y = z as Lame n = 2 (k^2 = 1/2)")
print(f"closed-form band edges (h): {', '.join(f'{e:.6f}' for e in edges)}")
print(f"gaps (h): ({edges[1]:.6f}, {edges[2]:.6f}) unreachable (h > 3);  ({gap2[0]:.6f}, {gap2[1]:.6f}) reachable")
print(f"PREDICTION: exactly one tongue, E in ({E_lo:.6f}, {E_hi:.6f}) = (1/8, 1/6) exactly "
      f"(h = 3 + sqrt3 -> A^4 = 1/12 -> E = 1/8;  h = 9/2 -> A^4 = 1/9 -> E = 1/6)")
rep('tongue edges == (1/8, 1/6) exactly', max(abs(E_lo - 1 / 8), abs(E_hi - 1 / 6)), 1e-15)


# ---------------------------------------------------------------- (1) Lame discriminant vs closed form
def lame_trace(h, period=None):
    """Trace of the monodromy of w_uu + (h - 3 sn^2 u) w = 0 over the potential's period 2K."""
    Kk = ellipk(K2)
    P = 2 * Kk if period is None else period

    def rhs(u, y):
        sn = ellipj(u, K2)[0]
        return [y[1], -(h - 3 * sn * sn) * y[0], y[3], -(h - 3 * sn * sn) * y[2]]
    sol = solve_ivp(rhs, (0, P), [1, 0, 0, 1], rtol=1e-12, atol=1e-13, method='DOP853')
    y = sol.y[:, -1]
    return y[0] + y[3]


errs = []
for e in edges:
    # an edge is where |trace| = 2: check |tr| - 2 changes sign across it (or touches 0 at a closed double point)
    d = [abs(lame_trace(e + s)) - 2 for s in (-1e-4, 1e-4)]
    errs.append(min(abs(abs(lame_trace(e)) - 2), 1.0))
    print(f"  h = {e:.6f}: |tr| - 2 at h -+ 1e-4: {d[0]:+.2e}, {d[1]:+.2e}   at h: {abs(lame_trace(e)) - 2:+.2e}")
rep('closed-form band edges: max ||tr(h_edge)| - 2|', max(errs), 1e-7)
inside = [abs(lame_trace(h)) - 2 for h in (4.6, 4.7, 2.0, 2.5)]
rep('inside both gaps |tr| > 2 (min of |tr| - 2 at h = 4.6, 4.7, 2.0, 2.5, negated)', -min(inside), 0)
above = max(abs(lame_trace(h)) - 2 for h in np.linspace(4.75, 40, 200))
rep('no instability anywhere above the last edge (h in [4.75, 40]): max |tr| - 2', above, 0)


# ---------------------------------------------------------------- (2) direct Floquet of the full 3x3 diagonal ghost
def orbit_and_monodromy(E):
    """Integrate x = y = z = s from the turning point together with the 3x3 ghost w'' = -(1 + Hess V(q)) w,
    over one orbit period (found from the quarter period). Returns (T, multipliers, max |s - A cn|)."""
    A = (2 * E / 3) ** 0.25
    om = np.sqrt(2) * A
    T = 4 * ellipk(K2) / om                            # period of s
    def hess(q):
        x, y, z = q
        return np.array([[y * y + z * z, 2 * x * y, 2 * x * z], [2 * x * y, x * x + z * z, 2 * y * z],
                         [2 * x * z, 2 * y * z, x * x + y * y]])

    def rhs(t, Y):
        q, p = Y[:3], Y[3:6]
        W = Y[6:].reshape(6, 6)                        # columns: fundamental solutions (w, w')
        x, y, z = q
        g = np.array([x * (y * y + z * z), y * (x * x + z * z), z * (x * x + y * y)])
        M = np.eye(3) + hess(q)
        dW = np.vstack([W[3:], -M @ W[:3]])
        return np.concatenate([p, -g, dW.ravel()])
    Y0 = np.concatenate([[A, A, A], [0, 0, 0], np.eye(6).ravel()])
    sol = solve_ivp(rhs, (0, T), Y0, rtol=1e-12, atol=1e-13, method='DOP853', dense_output=True)
    Mono = sol.y[6:, -1].reshape(6, 6)
    mu = np.linalg.eigvals(Mono)
    tt = np.linspace(0, T, 50)
    cn = ellipj(om * tt, K2)[1]
    err_orbit = np.abs(sol.sol(tt)[0] - A * cn).max() / A
    closure = np.abs(sol.y[:6, -1] - Y0[:6]).max() / A
    return T, mu, err_orbit, closure


Es = np.geomspace(0.02, 3.0, 300)
worst_orbit = worst_clo = 0
unstable = []
print("\nscan of the full 3x3 diagonal ghost, E in [0.02, 3], 300 log-spaced (the RH5 grid range):")
for E in Es:
    T, mu, eo, clo = orbit_and_monodromy(E)
    worst_orbit, worst_clo = max(worst_orbit, eo), max(worst_clo, clo)
    m = np.abs(mu).max()
    if m > 1 + 1e-6:
        unstable.append((E, m))
        print(f"  E={E:.5f} T={T:8.4f}  max|mu| = {m:.4f}   (predicted tongue {E_lo:.5f}-{E_hi:.5f})")
rep('orbit == A cn(sqrt2 A t | 1/2) along the integration (rel)', worst_orbit, 1e-9)
rep('orbit closure after one period (rel)', worst_clo, 1e-9)
inside_pred = all(E_lo <= E <= E_hi for E, _ in unstable)
rep('every unstable grid point lies inside the predicted tongue (count outside)',
    sum(not (E_lo <= E <= E_hi) for E, _ in unstable), 0)
rep('the predicted tongue is seen on the grid (1 if no unstable point)', 0 if unstable else 1, 0)
if unstable:
    print(f"  peak max|mu| on the grid = {max(m for _, m in unstable):.4f} at E = "
          f"{max(unstable, key=lambda r: r[1])[0]:.5f}  (RH5_nonplanar.md quoted 1.22 at E = 0.13-0.16)")


# ---------------------------------------------------------------- (3) tongue edges by bisection
def unstable_at(E):
    return np.abs(orbit_and_monodromy(E)[1]).max() > 1 + 1e-9


def bisect(a, b, ua):
    for _ in range(40):
        c = 0.5 * (a + b)
        if unstable_at(c) == ua:
            a = c
        else:
            b = c
    return 0.5 * (a + b)


mid = np.sqrt(E_lo * E_hi)
if unstable_at(mid):
    lo = bisect(mid, E_lo * 0.9, True)
    hi = bisect(mid, E_hi * 1.1, True)
    print(f"\ntongue edges by bisection: E in ({lo:.8f}, {hi:.8f}); closed form ({E_lo:.8f}, {E_hi:.8f})")
    rep('lower tongue edge vs closed form (rel)', abs(lo / E_lo - 1), 1e-6)
    rep('upper tongue edge vs closed form (rel)', abs(hi / E_hi - 1), 1e-6)
    peak = max(np.abs(orbit_and_monodromy(E)[1]).max() for E in np.geomspace(lo, hi, 41)[1:-1])
    print(f"  peak max|mu| inside the tongue (41 points) = {peak:.4f}")
else:
    rep('predicted tongue unstable at its geometric centre', 1, 0)

print('\nall checks passed' if OK else '\nSOME CHECKS FAILED')
print(f"# END rh5_lame status={'OK' if OK else 'FAIL'} wall={time.time() - T0:.0f}s", flush=True)
