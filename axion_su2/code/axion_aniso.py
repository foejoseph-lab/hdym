"""axion_aniso.py -- stage 3, first run: the anisotropic sector on the isotropic axion-SU(2) background.

Equation (axion_derive.py, CHECK 6; the two permutation-odd directions, degenerate, w0 drops out):
    eps'''' / M^2 + eps'' + kappa(t) Q(t) eps = 0,       kappa = (g/f) phidot   (flat space),
a Pais-Uhlenbeck (PU) oscillator with one time-dependent coefficient, integrated as the 4-d
first-order system y = (eps, eps', eps'', eps'''), never decomposed into two oscillators
(Kuwabara, Yumibayashi, Harada, arXiv:1503.03657).  Units g = M = 1 (time in 1/M, fields in M/g),
so with M fixed at 1 the scan parameter of the physical case is omega = m_phi / M.

Scales of the physical case (Fujita, Murai, Nakayama, Yin, arXiv:2312.06889, end of inflation,
their Fig. 5 benchmark): m_phi >> sqrt(kappa Q) >> g Q >> H, about 5e10 : 1e8 : 5e5 : 1e5 GeV.
The axion oscillates ~1e5 times per oscillation of Q, so over one axion period Q is frozen and
kappa(t) Q = c1 cos(omega t) exactly: the periodic drive is the physical fast problem, and the slow
Q cycle (and later the Hubble damping) sweeps c1.  The general periodic drive  c0 + c1 cos(omega t)
is kept as the control (c0 != 0 is the toy generalisation, never quoted as physics).

The result.  Two second-order effects, derived in the docstrings of Parts B and C and measured in D.
(i) With Q frozen, the periodic drive c1 cos(omega t): splitting eps into a slow part s and the fast
O(c1) response f, the response goes through the ghost pole, f = -c1 s cos(omega t)/(omega^2 (omega^2 - 1)),
and the period average gives  s'' = + c1^2 s / (2 omega^2 (omega^2 - 1))  (M = 1): for omega < 1 the
usual Kapitza stabilisation (omega -> 0 is the Hill limit), for omega > 1 tachyonic.  (ii) Q is not
frozen: its own response to kappa, Q'' = kappa Q^2, is the jitter dQ = phi_amp X^2 sin(omega t)/omega
about the slow amplitude X, and kappa dQ has the period average -phi_amp^2 X^2 / 2, a direct tachyonic
term in c(t), independent of M.  With c1 = omega phi_amp X the two add to
    s'' = + (phi_amp^2 X^2 / 2) * omega^2 / (omega^2 - 1) * s ,
so as M -> infinity (omega -> 0) they cancel exactly and the direction is marginal, as without the
axion (stage 1); for omega < 1 (M > m_phi) it is stable; for omega > 1 (M < m_phi) it is tachyonic
with rate  mu = phi_amp X omega / sqrt(2 (omega^2 - 1)), i.e. in physical units
    mu = g Q phi_amp / sqrt(2 (1 - r^2)),   r = M / m_phi,
not suppressed at small M.  On top of it: the ghost branch's parametric tongue at omega = 2 (ghost-pair
rate c1/4, below the band rate there) and the ghost / free-branch combination resonance at omega near
1, where the law diverges.

Parts:
  A  checks on the periodic drive: c1 = 0 gives the frozen PU multipliers exp(+-i omega_pm T);
     det = 1 and reciprocal pairs (PU is Hamiltonian); M -> infinity (here omega -> 0 with
     c0, c1 = O(omega^2)) reproduces the Hill/Mathieu equation eps'' + (c0 + c1 cos omega t) eps = 0,
     both against a direct 2-d integration and against scipy's Mathieu characteristic values.
  B  the map of the 4-d PU equation at c0 = 0.  Multipliers come as reciprocal pairs {lambda, 1/lambda};
     each pair is read off x = lambda + 1/lambda from the characteristic polynomial (the free branch
     sits at a Jordan block, x = 2, where eigenvalues are sqrt-sensitive; x itself is not).  B1 the
     unstable intervals in omega; B2 the frozen-coefficient tachyon (kappa Q < 0 half of every period)
     against the Kapitza threshold c1/omega^2 = 0.454 (Mathieu b_1(q) = 0 at q = 0.908); B3 the
     second-order law above against the measured rates and, for omega < 1, the measured frequency
     shift; B4 the combination resonance near omega = 1, its peak rate and width vs c1.
  C  the physical point: Fujita's numbers converted to g = M = 1; the full law in r = M/m_phi per
     Hubble time at the end of inflation; the frozen-Q fast problem (ghost-pole term only) measured
     along r, linear in c1 down to the physical c1 ~ 1e-5, with the tongue at M = m_phi/2 and the
     resonance at M near m_phi resolved; the Q-cycle average and an adiabatic estimate of the total
     e-folds through matter-dominated reheating (a quadrature of flat-space rates with Q ~ 1/a and
     phidot ~ a^-3/2; the FRW integration proper is the next run, not this one).
  D  direct two-frequency check at reduced scale separation omega / (g Q_amp) = 50, 200 (physical
     1e5): Q(t) from Q'' + 2 Q^3 = kappa Q^2 (Fujita (13) at H = 0) with the sinusoidal kappa, the
     4-d transfer matrix over 30 Q periods, growth slope of ln sigma_max per Q cycle (fitted over the
     last 20 cycles, where the free branch's linear drift ~ ln t is negligible) against the full law
     integrated along the slow amplitude (times the turning-point factor of the reduced slow equation
     s'' = mu(t)^2 s, which it also integrates) and against the frozen-Q (ghost-pole only) integral;
     omega = 0.8 (M > m_phi) and the M -> infinity Hill equation on the same background as the two
     controls that must show no growth; the background itself against the Kapitza-averaged quartic
     Q'' = -(2 + phi_amp^2) Q^3 with the slow amplitude X = Q_amp (1 + 2 phi_amp^2/(2 + phi_amp^2))^(1/4)
     that the jitter velocity at t = 0 gives.

Assumptions stated in the output: V(phi) quadratic for the oscillating axion (2312.06889 does not
specify V), phi_amp = phi_end/f = 1 (not stated there; everything in C scales linearly with it),
m = 0 (their dark-matter mass is irrelevant at this epoch), the M -> infinity background Q(t)
(the ghost-sourced part w0 = kappa Q^2/M^2 changes Q at relative order kappa Q/M^2).

Float only (DOP853, rtol 1e-11); nothing here is certified.  Usage (from the repository root, ~6 min on the cloud CPU):
    python axion_su2/code/axion_aniso.py > axion_su2/results/axion_aniso.txt
Added 2026-10-14.
"""
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.special import ellipk, mathieu_a, mathieu_b

RT = dict(method='DOP853', rtol=1e-11, atol=1e-13)
T0 = time.time()


def stamp(msg=''):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ------------------------------------------------------------------ the 4-d PU system, M = 1
def pu_monodromy(cfun, T, extra_rhs=None, y_extra0=None):
    """Transfer matrix of eps'''' = -(eps'' + c(t) eps) over [0, T]; c may depend on an extra
    background state (extra_rhs, y_extra0) integrated alongside.  Returns (Phi, extra_final)."""
    nx = 0 if y_extra0 is None else len(y_extra0)

    def rhs(t, y):
        xb = y[:nx]
        Y = y[nx:].reshape(4, 4)                       # columns = the four fundamental solutions
        c = cfun(t, xb)
        dY = np.empty_like(Y)
        dY[0] = Y[1]
        dY[1] = Y[2]
        dY[2] = Y[3]
        dY[3] = -(Y[2] + c * Y[0])
        out = dY.ravel() if nx == 0 else np.concatenate([extra_rhs(t, xb), dY.ravel()])
        return out
    y0 = np.eye(4).ravel() if nx == 0 else np.concatenate([y_extra0, np.eye(4).ravel()])
    sol = solve_ivp(rhs, (0.0, T), y0, **RT)
    y = sol.y[:, -1]
    return y[nx:].reshape(4, 4), y[:nx]


def periodic(c0, c1, omega):
    """Monodromy over one period of c(t) = c0 + c1 cos(omega t); multipliers; max rate per unit t."""
    T = 2.0 * np.pi / omega
    Phi, _ = pu_monodromy(lambda t, _x: c0 + c1 * np.cos(omega * t), T)
    lam = np.linalg.eigvals(Phi)
    return Phi, lam, np.max(np.log(np.abs(lam))) / T


def pairs(Phi, T):
    """Reciprocal pairs of a 4x4 symplectic monodromy: x = lambda + 1/lambda from the characteristic
    polynomial lambda^4 - A lambda^3 + B lambda^2 - A lambda + 1 (A = tr, B = sum of principal 2x2
    minors): x^2 - A x + (B - 2) = 0.  Returns [(x, rate)] for the two pairs, rate = |ln|lambda|| / T,
    the pair nearer x = 2 (the free Yang-Mills branch, a Jordan block at c1 = 0) first.  Complex x
    means the two pairs have merged into a quadruplet (Krein collision); both then carry the rate."""
    A = np.trace(Phi)
    B = 0.5 * (A * A - np.trace(Phi @ Phi))
    disc = A * A - 4.0 * (B - 2.0)
    out = []
    for sgn in (1, -1):
        x = 0.5 * (A + sgn * np.sqrt(disc + 0j))
        lam = 0.5 * (x + np.sqrt(x * x - 4.0 + 0j))
        out.append((x, abs(np.log(abs(lam))) / T))
    out.sort(key=lambda p: abs(p[0] - 2.0))
    return out


def rate_of(c0, c1, omega):
    """Max rate over the two pairs, per unit t (pair analysis, not eigenvalues)."""
    T = 2.0 * np.pi / omega
    Phi, _ = pu_monodromy(lambda t, _x: c0 + c1 * np.cos(omega * t), T)
    return max(r for _, r in pairs(Phi, T))


def second_order(c1, omega):
    """Frozen-Q (ghost-pole) law: s'' = + c1^2 s / (2 omega^2 (omega^2 - 1)); rate for omega > 1, else 0."""
    return c1 / (omega * np.sqrt(2.0 * (omega ** 2 - 1.0))) if omega > 1 else 0.0


def full_law(phiX, omega):
    """Full law with Q's own jitter: s'' = + (phi X)^2/2 * omega^2/(omega^2 - 1) s; rate for omega > 1."""
    return phiX * omega / np.sqrt(2.0 * (omega ** 2 - 1.0)) if omega > 1 else 0.0


def hill_monodromy(c0, c1, omega):
    T = 2.0 * np.pi / omega

    def rhs(t, y):
        c = c0 + c1 * np.cos(omega * t)
        return [y[1], -c * y[0], y[3], -c * y[2]]
    y = solve_ivp(rhs, (0.0, T), [1, 0, 0, 1], **RT).y[:, -1]
    return np.array([[y[0], y[2]], [y[1], y[3]]])


def mathieu_unstable(a, q, nmax=12):
    """Mathieu y'' + (a - 2q cos 2z) y = 0 unstable?  Tongues: a < a_0(q), b_n(q) < a < a_n(q)."""
    q = abs(q)
    if a < mathieu_a(0, q):
        return True
    for n in range(1, nmax + 1):
        if mathieu_b(n, q) < a < mathieu_a(n, q):
            return True
    return False


def intervals(grid, flags):
    """Maximal runs of True in flags over grid: list of (lo, hi, count)."""
    out, start = [], None
    for i, f in enumerate(flags):
        if f and start is None:
            start = i
        if (not f or i == len(flags) - 1) and start is not None:
            end = i if f else i - 1
            out.append((grid[start], grid[end], end - start + 1))
            start = None
    return out


# ------------------------------------------------------------------ Part A: checks
def part_A():
    print("\nPart A: checks on the periodic drive c(t) = c0 + c1 cos(omega t), M = 1\n")
    # A1 frozen limit c1 = 0: multipliers exp(+-i omega_pm T), omega_pm^2 = (1 +- sqrt(1 - 4 c0))/2
    print("A1  c1 = 0: multipliers against the frozen PU frequencies (max |difference|)")
    worst = 0.0
    for c0, omega in ((0.05, 0.7), (0.2, 1.3), (0.24, 2.0), (-0.3, 0.9)):
        T = 2 * np.pi / omega
        _, lam, _ = periodic(c0, 0.0, omega)
        wp2 = (1 + np.sqrt(1 - 4 * c0 + 0j)) / 2
        wm2 = (1 - np.sqrt(1 - 4 * c0 + 0j)) / 2
        ref = np.array([np.exp(s * 1j * np.sqrt(w2) * T) for w2 in (wp2, wm2) for s in (1, -1)])
        d = max(np.min(np.abs(lam - r)) for r in ref)
        worst = max(worst, d)
        print(f"    c0 = {c0:5.2f} omega = {omega:.1f}: max |lambda - exp(+-i omega_pm T)| = {d:.1e}")
    print(f"    PASS" if worst < 1e-8 else f"    FAIL ({worst:.1e})")

    # A2 Liouville / reciprocal pairs
    print("\nA2  det Phi = 1 and reciprocal pairs {lambda, 1/lambda} (PU is Hamiltonian)")
    worst_det, worst_pair = 0.0, 0.0
    for c0, c1, omega in ((0.0, 0.3, 2.0), (0.1, 0.2, 1.0), (0.0, 0.8, 0.5), (0.3, 0.1, 1.9)):
        Phi, lam, _ = periodic(c0, c1, omega)
        worst_det = max(worst_det, abs(np.linalg.det(Phi) - 1))
        worst_pair = max(worst_pair, max(np.min(np.abs(lam - 1 / l)) for l in lam))
    print(f"    max |det - 1| = {worst_det:.1e}, max |lambda - 1/lambda'| = {worst_pair:.1e}: "
          f"{'PASS' if worst_det < 1e-8 and worst_pair < 1e-7 else 'FAIL'}")

    # A3 M -> infinity: omega -> 0 with Mathieu (a, q) fixed, c0 = a omega^2/4, c1 = -q omega^2/2.
    #    Mathieu y'' + (a - 2q cos 2z) y = 0 with z = omega t/2: a = 4 c0/omega^2, q = -2 c1/omega^2.
    print("\nA3  M -> infinity (omega -> 0 at fixed Mathieu a = 4c0/omega^2, q = -2c1/omega^2):")
    print("    rate of the 4-d PU system vs the 2-d Hill equation (per unit t), and scipy's verdict")
    print(f"    {'a':>5} {'q':>5} {'omega':>6} | {'PU rate':>10} {'Hill rate':>10} {'diff':>8} | {'PU unst':>7} {'Hill':>5} {'mathieu_a/b':>11}")
    ok = True
    for a, q in ((1.0, 0.3), (0.5, 0.3), (0.0, 0.5), (0.0, 1.2), (4.0, 1.0), (-0.1, 0.2)):
        verdict = mathieu_unstable(a, q)
        for omega in (0.1, 0.03, 0.01):
            c0, c1 = a * omega ** 2 / 4, -q * omega ** 2 / 2
            T = 2 * np.pi / omega
            _, lam, rate = periodic(c0, c1, omega)
            Mh = hill_monodromy(c0, c1, omega)
            lh = np.linalg.eigvals(Mh)
            rate_h = np.max(np.log(np.abs(lh))) / T
            unst_pu, unst_h = rate > 1e-6 * omega, rate_h > 1e-6 * omega
            if omega == 0.01:
                ok = ok and (unst_pu == verdict) and (unst_h == verdict) and abs(rate - rate_h) < 1e-3 * max(rate_h, omega)
            print(f"    {a:5.2f} {q:5.2f} {omega:6.2f} | {rate:10.3e} {rate_h:10.3e} {rate - rate_h:8.1e} | "
                  f"{str(unst_pu):>7} {str(unst_h):>5} {str(verdict):>11}")
    print(f"    (the ghost branch at frequency 1 >> omega is far from its own resonances 2/n at these omega; "
          f"the difference falls with omega^2)\n    {'PASS' if ok else 'FAIL'} at omega = 0.01")


# ------------------------------------------------------------------ Part B: the map at c0 = 0
def part_B():
    print("\nPart B: the map of  eps'''' + eps'' + c1 cos(omega t) eps = 0  (c0 = 0, M = 1)\n")
    print("B1  unstable intervals in omega (grid 0.25-4.5, step 1.25e-3; rate > 1e-6 per unit t, pair analysis)")
    grid = np.linspace(0.25, 4.5, 3401)
    for c1 in (0.01, 0.03, 0.1):
        rates = np.array([rate_of(0.0, c1, om) for om in grid])
        iv = intervals(grid, rates > 1e-6)
        print(f"  c1 = {c1:5.2f}: {len(iv)} interval(s)")
        for lo, hi, n in iv:
            sel = (grid >= lo) & (grid <= hi)
            k = np.argmax(rates[sel])
            print(f"     omega in [{lo:.4f}, {hi:.4f}]  width {hi - lo:.4f}  max rate {rates[sel][k]:.4g} "
                  f"at omega = {grid[sel][k]:.4f}  ({n} pts)")
        for om in (0.5, 0.75, 1.5, 2.0, 3.0, 4.0):
            k = np.argmin(abs(grid - om))
            print(f"       rate at omega = {grid[k]:.3f}: {rates[k]:.3e}   second-order law: {second_order(c1, grid[k]):.3e}")
    print("  Reading: below omega = 1 stable except the ghost's 2/n tongues (O(c1^n)) and, at larger c1, the")
    print("  frozen-tachyon regime (B2); near omega = 1 the combination resonance (B4); above omega = 1 the whole")
    print("  range is unstable, by the second-order law (B3), with the ghost's first tongue at omega = 2 on top.")

    print("\nB2  the frozen-coefficient tachyon (kappa Q < 0 for half of every period) at small omega:")
    print("    threshold in c1/omega^2 of the 4-d system vs Mathieu b_1(q) = 0 (the a = 0 line enters the first")
    print("    tongue b_1 < a < a_1 at q = 2c1/omega^2 = 0.9080, c1/omega^2 = 0.4540)")
    for omega in (0.1, 0.03):
        lo, hi = 0.3, 0.7
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            if rate_of(0.0, mid * omega ** 2, omega) > 1e-7:
                hi = mid
            else:
                lo = mid
        print(f"    omega = {omega}: threshold c1/omega^2 = {0.5 * (lo + hi):.4f}")
    lo, hi = 0.5, 1.5
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if mathieu_b(1, mid) > 0:
            lo = mid
        else:
            hi = mid
    print(f"    Mathieu: b_1(q) = 0 at q = {0.5 * (lo + hi):.4f}, c1/omega^2 = {0.25 * (lo + hi):.4f}")
    print("    Below it the fast sign-changing coefficient is Kapitza-stabilised (omega << 1): no growth although")
    print("    the frozen Yang-Mills branch omega^2 = kappa Q is tachyonic half the time.")

    print("\nB3  the frozen-Q (ghost-pole) law.  omega > 1: measured rate of the free-branch pair vs c1/(omega sqrt(2(omega^2-1)));")
    print("    omega < 1: measured x - 2 of that pair vs -(nu T)^2, nu^2 = c1^2/(2 omega^2 (1 - omega^2)) (a stable shift).")
    print(f"    {'omega':>6} {'c1':>7} | {'rate meas':>10} {'rate law':>10} {'ratio':>7} | {'x-2 meas':>11} {'-(nu T)^2':>11}")
    worst = 0.0
    for omega in (1.2, 1.5, 2.0, 3.0, 5.0, 8.0, 0.9, 0.7, 0.55, 0.3):
        for c1 in (0.001, 0.01, 0.1):
            T = 2 * np.pi / omega
            Phi, _ = pu_monodromy(lambda t, _x: c1 * np.cos(omega * t), T)
            (x, r), _ = pairs(Phi, T)
            law = second_order(c1, omega)
            nu2 = c1 ** 2 / (2 * omega ** 2 * (1 - omega ** 2)) if omega < 1 else np.nan
            ratio = r / law if law > 0 else np.nan
            if omega >= 1.5 and c1 <= 0.01:
                worst = max(worst, abs(ratio - 1))
            print(f"    {omega:6.2f} {c1:7.3f} | {r:10.3e} {law:10.3e} {ratio:7.3f} | {x.real - 2:11.3e} "
                  f"{-nu2 * T * T if omega < 1 else np.nan:11.3e}")
    print(f"    worst |ratio - 1| for omega >= 1.5, c1 <= 0.01: {worst:.1e}  ({'PASS' if worst < 0.02 else 'FAIL'}; "
          f"corrections O(c1^2) and O(1/(omega^2-1)) near omega = 1)")
    print("    (omega = 1/n exactly is avoided: there the ghost pair also sits at x = 2, the monodromy has a")
    print("    quadruple multiplier 1 and the pair analysis is eps^(1/4)-sensitive; the tongues there are O(c1^(2/omega)))")

    print("\nB4  the combination resonance near omega = 1 (ghost at 1, free branch at 0): peak rate and width vs c1")
    print(f"    {'c1':>7} {'bubble below 1':>22} {'width':>8} {'peak rate':>10} {'at omega':>9} | {'law there':>10}")
    fine = np.linspace(0.85, 1.25, 1601)
    rows = []
    for c1 in (0.001, 0.003, 0.01, 0.03):
        rates = np.array([rate_of(0.0, c1, om) for om in fine])
        above = rates > 1e-6
        iv = [v for v in intervals(fine, above) if v[1] < 1.0]
        bub = iv[-1] if iv else (np.nan, np.nan, 0)
        kpk = np.argmax(rates)
        rows.append((c1, bub[0], bub[1], rates[kpk], fine[kpk]))
        print(f"    {c1:7.3f} [{bub[0]:.4f}, {bub[1]:.4f}] {bub[1] - bub[0]:8.4f} {rates[kpk]:10.4e} {fine[kpk]:9.4f} | "
              f"{second_order(c1, fine[kpk]) if fine[kpk] > 1 else np.nan:10.3e}")
    c = np.array([r[0] for r in rows])
    w = np.array([r[2] - r[1] for r in rows])
    pk = np.array([r[3] for r in rows])
    pos = np.array([r[4] for r in rows]) - 1
    ok = np.isfinite(w) & (w > 0)
    print(f"    bubble width ~ c1^{np.polyfit(np.log(c[ok]), np.log(w[ok]), 1)[0]:.2f}, "
          f"peak rate ~ c1^{np.polyfit(np.log(c), np.log(pk), 1)[0]:.2f}, peak position - 1 ~ "
          f"c1^{np.polyfit(np.log(c), np.log(pos), 1)[0]:.2f}")
    print("    (the bubble below 1 is a Krein collision, complex x; above 1 the resonance merges into the")
    print("    second-order band, which the law overestimates as 1/sqrt(omega^2 - 1) -> infinity)")


# ------------------------------------------------------------------ Part C: the physical point
GEV = dict(MPl=2.4e18, g=1e-12, mQ=5.0, mphi=5e10, f_over_MPl=0.004)


def part_C(phi_amp=1.0):
    print("\nPart C: Fujita et al. 2312.06889 end-of-inflation data in g = M = 1 units\n")
    MPl, g, mQ, mphi = GEV['MPl'], GEV['g'], GEV['mQ'], GEV['mphi']
    H = g * MPl / mQ ** 2                                    # their (33)
    Qend = MPl / mQ                                          # their (32)
    gQ = g * Qend
    xi = mQ + 1 / mQ                                         # their (20)-(21) at m_Q = 5
    q_F = gQ / mphi                                          # their q = gQ/m_phi at the end
    h = H / mphi
    f = GEV['f_over_MPl'] * MPl
    print(f"  their Fig. 5 benchmark: g = {g:.0e}, m_Q,end = {mQ:g}, m_phi = {mphi:.1e} GeV, f = {GEV['f_over_MPl']} M_Pl "
          f"(Lambda_Q = Q/f = {Qend / f:.0f})")
    print(f"  H_end = g M_Pl/m_Q^2 = {H:.2e} GeV,  Q_end = M_Pl/m_Q = {Qend:.2e} GeV,  g Q_end = m_Q H_end = {gQ:.2e} GeV")
    print(f"  their variables at the end: q = gQ/m_phi = {q_F:.2e},  h = H/m_phi = {h:.2e},  xi = m_Q + 1/m_Q = {xi:.2f}, "
          f"dphi/dx = 2 xi h = {2 * xi * h:.1e}")
    print(f"  assumed: V(phi) quadratic after inflation, phi_amp = phi/f oscillation amplitude = {phi_amp:g} "
          f"(not stated in 2312.06889; everything below is linear in it), m = 0")
    print(f"  oscillating stage: kappa = (g/f) phidot = -g m_phi phi_amp sin(m_phi t);  "
          f"kappa_amp Q_end = m_phi^2 q phi_amp = ({np.sqrt(mphi ** 2 * q_F * phi_amp):.1e} GeV)^2")
    print(f"  scales: m_phi : sqrt(kappa Q) : gQ : H = 1 : {np.sqrt(q_F * phi_amp):.1e} : {q_F:.1e} : {h:.1e}")
    print("  conversion to g = M = 1 (time 1/M, fields M/g): omega = m_phi/M, Q = q m_phi/M, H/M = h m_phi/M,")
    print("    kappa/(gM) = 2 xi H/M, the Hill coefficient kappa Q/M^2 = q (dphi/dx) (m_phi/M)^2 = c1 cos(omega t)")
    print(f"    with c1 = omega^2 q phi_amp = {q_F * phi_amp:.2e} omega^2  (frozen Q = Q_end; the Q cycle scales it by cn).")
    print("  The jitter of Q (its response to kappa): dQ/Q = phi_amp g Q / m_phi = phi_amp q = "
          f"{q_F * phi_amp:.1e}, and kappa dQ averages to -(phi_amp g Q)^2 / 2.")
    c1_of = lambda om: q_F * phi_amp * om ** 2
    lawH_full = lambda r: mQ * phi_amp / np.sqrt(2 * (1 - r * r)) if r < 1 else 0.0
    lawH_pole = lambda r: mQ * phi_amp * r / np.sqrt(2 * (1 - r * r)) if r < 1 else 0.0

    print("\nC1  the law in these variables, r = M/m_phi < 1, per Hubble time at the end of inflation:")
    print("    full:        mu = g Q phi_amp / sqrt(2 (1 - r^2)),     mu/H_end = m_Q phi_amp / sqrt(2 (1 - r^2))")
    print("    ghost pole:  mu = g Q phi_amp r / sqrt(2 (1 - r^2)),   the part the frozen-Q fast problem carries")
    print("    jitter:      mu = g Q phi_amp / sqrt(2),               the M-independent part (cancelled by the")
    print("                 ghost-pole term as M -> infinity: without the ghost the direction is marginal)")
    print(f"    {'r = M/m_phi':>11} {'full/H_end':>10} {'pole/H_end':>10} | {'frozen-Q rate (M=1)':>19} {'pole law (M=1)':>14}")
    for r in (0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.97, 1.03, 1.1, 1.5, 2.5):
        om = 1.0 / r
        rate = rate_of(0.0, c1_of(om), om)
        print(f"    {r:11.2f} {lawH_full(r):10.3f} {lawH_pole(r):10.3f} | {rate:19.3e} {second_order(c1_of(om), om):14.3e}")
    print("    (frozen-Q rates measured by the pair analysis at the physical c1; r = 0.5 carries the ghost tongue,")
    print("     r = 0.97 the combination resonance; the full law is checked in Part D, not here: a frozen Q")
    print("     cannot carry the jitter term.  The scan avoids r = integer, where omega = 1/n is degenerate.)")

    print("\nC2  linearity of the frozen-Q rate in c1 at r = 0.3 and 0.7: rate / c1 at 1000x, 100x, 10x and 1x the physical c1")
    for r in (0.3, 0.7):
        om = 1.0 / r
        vals = [rate_of(0.0, k * c1_of(om), om) / (k * c1_of(om)) for k in (1000, 100, 10, 1)]
        print(f"    r = {r}: rate/c1 = " + "  ".join(f"{v:.5f}" for v in vals) + f"   pole law {second_order(1.0, om):.5f}")

    print("\nC3  the ghost tongue at r = 1/2 (omega = 2), frozen Q: relative half-width q phi_amp, ghost-pair rate c1/4")
    delta = np.linspace(-2, 2, 17) * q_F * phi_amp
    M_over_H = mphi / 2 / H
    print(f"    {'(omega-2)/2 / (q phi)':>22} {'rate (M=1)':>12} {'rate/H_end':>11}")
    for d in delta:
        om = 2 * (1 + d)
        rr = rate_of(0.0, c1_of(om), om)
        print(f"    {d / (q_F * phi_amp):22.2f} {rr:12.4e} {rr * M_over_H:11.3f}")
    print(f"    the pole band at r = 1/2 gives {lawH_pole(0.5):.3f} per Hubble time, the tongue's ghost pair c1/4 = "
          f"{q_F * phi_amp * M_over_H:.3f} (a separate pair: the larger counts, they do not add);")
    print(f"    the full band, {lawH_full(0.5):.3f} per Hubble time, is above both: the tongue is not the leading feature.")

    print("\nC4  the combination resonance near r = 1 (omega near 1) at the physical c1, frozen Q, fine scan in omega")
    fine = 1.0 + np.linspace(-2e-3, 3e-3, 51)
    rates = np.array([rate_of(0.0, c1_of(om), om) for om in fine])
    kpk = np.argmax(rates)
    unst = intervals(fine, rates > 1e-6)
    print(f"    peak rate {rates[kpk]:.3e} (M = 1) at omega = {fine[kpk]:.5f}, i.e. r = {1 / fine[kpk]:.5f}; "
          f"{rates[kpk] * mphi / fine[kpk] / H:.1f} per Hubble time")
    for lo, hi, n in unst:
        print(f"    unstable (> 1e-6) for omega in [{lo:.5f}, {hi:.5f}] ({n} pts)")
    print(f"    rate at omega = 1 - 1e-3: {rates[np.argmin(abs(fine - 0.999))]:.3e}; at 1 + 2e-3: {rates[np.argmin(abs(fine - 1.002))]:.3e}")
    print("    (indicative only: the full law diverges at r = 1 and the jitter term is not in a frozen-Q scan)")

    from scipy.special import ellipj
    K = ellipk(0.5)
    tau = np.linspace(0, 4 * K, 4001)
    _, cn, _, _ = ellipj(tau, 0.5)
    mean_abs_cn = np.trapezoid(np.abs(cn), tau) / (4 * K)
    print("\nC5  average over the Q cycle and the adiabatic estimate through matter-dominated reheating")
    print(f"    rates are linear in |X|: <|cn|> over a period = {mean_abs_cn:.4f}; the reduced slow equation")
    print("    s'' = mu(t)^2 s grows 1.12x faster per cycle than exp(integral mu dt) (Part D).  With Q ~ 1/a and")
    print("    phidot ~ a^-3/2 the rate falls as a^-5/2 = (t_end/t)^(5/3) while H = 2/(3t): total e-folds =")
    print("    (rate/H)_end, i.e.  N = 1.12 <|cn|> m_Q phi_amp / sqrt(2 (1 - r^2))  (flat-space rates, the FRW run")
    print(f"    is next, not this).  phi_amp = {phi_amp:g}:")
    print("    " + "  ".join(f"r={r}: N={1.12 * mean_abs_cn * lawH_full(r):.2f}" for r in (0.0, 0.3, 0.5, 0.7, 0.9, 0.97)))
    print("    For r > 1 (M > m_phi) the sector is stable at this order (the jitter tachyon is over-compensated);")
    print(f"    the frozen tachyon of the fast problem is at c1/omega^2 = q phi_amp = {q_F * phi_amp:.1e}, five orders")
    print("    below the Kapitza threshold 0.454.")


# ------------------------------------------------------------------ Part D: direct two-frequency check
def part_D():
    print("\nPart D: direct integration with the oscillating background Q(t), reduced scale separation\n")
    print("  background Q'' + 2 Q^3 = kappa Q^2 with kappa = -omega phi_amp sin(omega t) (Fujita (13) at H = 0,")
    print("  M -> infinity), started at Q = Q_amp, Q' = 0; the anisotropic 4-d transfer matrix alongside for")
    print("  30 Q periods; growth slope of ln sigma_max per Q cycle fitted over cycles 10-30.  Predictions per")
    print("  cycle: 'full' = integral of the full law along the slow amplitude |X(t)| (X = the background with its")
    print("  jitter removed by a running average over the drive period), 'slow eq' = the reduced s'' = mu(t)^2 s")
    print("  integrated with that mu (its turning-point factor), 'pole' = the frozen-Q integral (ghost-pole term")
    print("  only, mu from the pair analysis along |Q(t)|).  'Hill' = the M -> infinity equation eps'' + kappa Q eps")
    print("  = 0 on the same background (the no-ghost control).")
    K = ellipk(0.5)
    print(f"\n  {'omega':>6} {'ratio':>6} {'phi':>4} {'T_Q meas':>9} {'T_Q Kap':>8} | {'slope':>7} {'full':>7} {'slow eq':>7} "
          f"{'pole':>7} | {'Hill':>7} | {'ln s_max end':>12}")
    for omega, ratio, phi_amp in ((2.0, 50, 2.0), (2.0, 200, 2.0), (2.0, 200, 1.0), (3.0, 200, 2.0),
                                  (1.5, 200, 2.0), (0.8, 200, 2.0), (0.6, 200, 2.0)):
        Q_amp = omega / ratio
        lam_eff = 2.0 + phi_amp ** 2
        X_amp = Q_amp * (1.0 + 2.0 * phi_amp ** 2 / lam_eff) ** 0.25
        T_kap = 4 * K / (np.sqrt(lam_eff) * X_amp)
        kappa = lambda t: -omega * phi_amp * np.sin(omega * t)

        def bg(t, x):
            return np.array([x[1], -2.0 * x[0] ** 3 + kappa(t) * x[0] ** 2])
        n_cyc = 30
        Tend = (n_cyc + 1) * T_kap
        sb = solve_ivp(lambda t, x: bg(t, x), (0.0, Tend), [Q_amp, 0.0], dense_output=True, **RT)
        tt = np.linspace(0, Tend, 200001)
        Qs = sb.sol(tt)[0]
        zc = tt[:-1][(Qs[:-1] > 0) & (Qs[1:] <= 0)]
        T_meas = np.mean(np.diff(zc)) if len(zc) > 2 else np.nan
        T_run = n_cyc * T_meas
        nx = 2

        def rhs(t, y):
            xb = y[:nx]
            Y = y[nx:].reshape(4, 4)
            c = kappa(t) * xb[0]
            dY = np.empty_like(Y)
            dY[0] = Y[1]
            dY[1] = Y[2]
            dY[2] = Y[3]
            dY[3] = -(Y[2] + c * Y[0])
            return np.concatenate([bg(t, xb), dY.ravel()])
        y0 = np.concatenate([[Q_amp, 0.0], np.eye(4).ravel()])
        t_eval = np.arange(1, n_cyc + 1) * T_meas
        sol = solve_ivp(rhs, (0.0, T_run), y0, t_eval=t_eval, **RT)
        lns = np.array([np.log(np.linalg.svd(sol.y[nx:, k].reshape(4, 4), compute_uv=False)[0]) for k in range(n_cyc)])
        slope = np.polyfit(np.arange(10, n_cyc + 1), lns[9:], 1)[0]

        # slow amplitude X(t): running average of Q over one drive period
        Td = 2 * np.pi / omega
        tt2 = np.linspace(0, T_meas, 4001)
        X_t = np.array([np.mean(sb.sol(np.linspace(t - Td / 2, t + Td / 2, 41))[0]) for t in tt2])
        mu_full = np.array([full_law(phi_amp * abs(x), omega) for x in X_t])
        full = np.trapezoid(mu_full, tt2)
        # reduced slow equation with that mu(t), over the 30 cycles (mu periodic with T_meas)
        mu_interp = interp1d(tt2, mu_full, kind='linear')

        def slow_rhs(t, y):
            m2 = mu_interp(t % T_meas) ** 2
            return [y[1], m2 * y[0], y[3], m2 * y[2]]
        ss = solve_ivp(slow_rhs, (0.0, T_run), [1, 0, 0, 1], t_eval=t_eval, **RT)
        lns2 = [np.log(np.linalg.svd(np.array([[ss.y[0, k], ss.y[2, k]], [ss.y[1, k], ss.y[3, k]]]),
                                     compute_uv=False)[0]) for k in range(n_cyc)]
        slow = np.polyfit(np.arange(10, n_cyc + 1), lns2[9:], 1)[0]
        # frozen-Q (pole) integral
        c1_grid = np.linspace(0.0, omega * phi_amp * X_amp * 1.1, 31)
        mu_grid = np.array([rate_of(0.0, c, omega) for c in c1_grid])
        mu_of = interp1d(c1_grid, mu_grid, kind='linear')
        c1_t = omega * phi_amp * np.abs(sb.sol(tt2)[0])
        pole = np.trapezoid(mu_of(np.minimum(c1_t, c1_grid[-1])), tt2)
        # Hill (M -> infinity) on the same background
        def hill_rhs(t, y):
            xb = y[:2]
            c = kappa(t) * xb[0]
            return np.concatenate([bg(t, xb), [y[3], -c * y[2], y[5], -c * y[4]]])
        hs = solve_ivp(hill_rhs, (0.0, T_run), [Q_amp, 0.0, 1, 0, 0, 1], t_eval=t_eval, **RT)
        lnh = [np.log(np.linalg.svd(np.array([[hs.y[2, k], hs.y[4, k]], [hs.y[3, k], hs.y[5, k]]]),
                                    compute_uv=False)[0]) for k in range(n_cyc)]
        hill = np.polyfit(np.arange(10, n_cyc + 1), lnh[9:], 1)[0]
        print(f"  {omega:6.2f} {ratio:6d} {phi_amp:4.1f} {T_meas:9.2f} {T_kap:8.2f} | {slope:7.3f} {full:7.3f} {slow:7.3f} "
              f"{pole:7.3f} | {hill:7.3f} | {lns[-1]:12.2f}", flush=True)
    print("\n  T_Q meas: measured Q period (zero-crossing spacing); T_Q Kap: the Kapitza-averaged quartic with the")
    print("  slow amplitude X_amp = Q_amp (1 + 2 phi^2/(2 + phi^2))^(1/4) set by the jitter velocity at t = 0.")
    print("  omega = 0.8 and 0.6 (M > m_phi) and the Hill column (no ghost) are the controls: residual drift")
    print("  ~ ln(t) only.  'slow eq' is the prediction to compare with 'slope'; 'pole' alone falls short by the")
    print("  jitter term, 'full' by the turning-point factor.")


if __name__ == "__main__":
    print(f"axion_aniso.py  python {sys.version.split()[0]}  numpy {np.__version__}  ({sys.executable})")
    print("units g = M = 1; eps'''' + eps'' + c(t) eps = 0 as the 4-d first-order system; DOP853 rtol 1e-11")
    part_A()
    stamp("Part A done")
    part_B()
    stamp("Part B done")
    part_C()
    stamp("Part C done")
    part_D()
    stamp("Part D done")
    print(f"\n({time.time() - T0:.0f} s)")
