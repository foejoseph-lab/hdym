"""axion_twins.py -- the twin suite: the massless anisotropic direction coupled by the axion drive to heavy
modes of either Krein sign, on the isotropic axion-SU(2) background.

The Pais-Uhlenbeck equation of axion_aniso.py, eps''''/M^2 + eps'' + c(t) eps = 0 with c = kappa Q, is the
Lagrangian  L = 1/2 s'^2 - (1/2 h'^2 - 1/2 M^2 h^2) - 1/2 c(t) (s + h)^2  with eps = s + h: a massless
direction s (the anisotropic Yang-Mills direction, free without the axion) and the ghost h at M, both seeing
the drive through the one combination s + h.  The suite keeps that Lagrangian and changes only the heavy
sector:
    L = 1/2 s'^2 + sum_i sigma_i (1/2 h_i'^2 - 1/2 M_i^2 h_i^2) - 1/2 c(t) (s + sum_i h_i)^2 ,
sigma_i = -1 ghost, +1 healthy.  Equations:  s'' = -c (s + sum h),   h_i'' + M_i^2 h_i = -sigma_i c (s + sum h).
Cases (units g = M = 1, M = the Lee-Wick mass of axion_aniso.py; omega = m_phi/M):
    ghost            {(-1, M)}                      the run of 14 Oct, as reference
    healthy twin     {(+1, M)}                      the same Lagrangian, the heavy mode's sign flipped
    ghost + healthy  {(-1, M), (+1, M_h)}, M_h = 2M, M/2, M
                     M_h = 2M is the sixth-order Lee-Wick operator with poles 0, M, 2M and Krein signs +, -, +
                     (bianchi9's two-pole theory, which there stopped at the beta' = 0 singularity; here the
                     system is linear and regular); M_h = M is the tuned point; M_h = M/2 the reversed order
    heavy mode alone h'' + (M^2 + sigma c) h = 0,   no massless direction: Mathieu, narrow tongues at omega = 2M/n

[Withdrawn 15 Oct (README, History): the "full law" and the sign rule "tachyonic iff D(m_phi^2) > 0, healthy
twin unstable on the other side" below are not physics.  They were measured in Part D on the M -> infinity
background Q'' = -2Q^3 + kappa Q^2, on which the jitter carries only the massless pole -1/w; on the self-consistent
background the jitter goes through the full Pi(w) = -1/w + D(w) and the two second-order terms cancel identically,
for any heavy sector (axion_dressed.py A1).  Parts A-C stand as statements about the periodic problem, and the
sign rule is correct for the pole term alone (README, "The twin suite").  The law is kept below as written.]

The second-order law (derived as in axion_aniso.py, now for any heavy sector).  With Q frozen the drive is
c = c1 cos(omega t), c1 = omega phi X (phi = phi_amp/f, X the slow amplitude of Q).  The fast O(c1) responses are
f = c1 s cos(omega t)/omega^2 for s itself and h_i = -sigma_i c1 s cos(omega t)/(M_i^2 - omega^2) through each
pole; averaging -c (f + sum h_i) over the drive period gives the frozen-Q slow equation
    s'' = -(c1^2/2) [ 1/omega^2 - sum_i sigma_i/(M_i^2 - omega^2) ] s                      (pole terms).
Q is not frozen: its jitter under kappa, dQ = phi X^2 sin(omega t)/omega, makes kappa dQ average to -phi^2 X^2/2,
which cancels the 1/omega^2 (Kapitza) term exactly, leaving
    s'' = + (phi^2 X^2 omega^2 / 2) * D(omega^2) * s,     D(w) = sum_i sigma_i/(M_i^2 - w)                (full law):
the massless direction is tachyonic iff D(m_phi^2) > 0, i.e. iff the heavy part of the propagator of the
higher-derivative operator, evaluated at the drive frequency, is positive; rate
    mu = phi X omega sqrt(D(omega^2)/2),     in physical units  mu = g Q phi m_phi sqrt(D(m_phi^2)/2).
Predictions: ghost, D = 1/(omega^2 - 1): tachyonic for M < m_phi (the 14 Oct result).  Healthy twin,
D = 1/(1 - omega^2): tachyonic for M > m_phi, stable for M < m_phi -- the mirror image.  Ghost + healthy
at M_h: D = 1/(omega^2 - 1) + 1/(M_h^2 - omega^2) = (M_h^2 - 1)/((omega^2 - 1)(M_h^2 - omega^2)): for M_h > M
tachyonic iff M < m_phi < M_h (the sixth-order theory's window closes above its top pole); for M_h < M
tachyonic iff m_phi is outside [M_h, M]; for M_h = M, D = 0 and the direction is marginal at every m_phi.
No heavy mode (pure Yang-Mills): D = 0, marginal (the Hill control of axion_aniso.py).

Measurement.  Frozen-Q: the 2n-dim monodromy over one drive period; multipliers come in reciprocal pairs,
read off x = lambda + 1/lambda from the palindromic characteristic polynomial (Faddeev-LeVerrier
coefficients, then the degree-n polynomial in x), which is well conditioned at the massless Jordan block
where eigenvalues are not.  Direct: the two-frequency integration on the driven background Q'' + 2Q^3 =
kappa Q^2 (Fujita (13) at H = 0, M -> infinity), 30 Q periods, growth slope of ln sigma_max per Q cycle over
cycles 10-30, against the reduced slow equation s'' = mu(t)^2 s with mu from the full law along the slow
amplitude X(t), as in axion_aniso.py Part D.  omega = 1/n and the combination frequencies (M_i +- M_j)/n are
degenerate points of the pair analysis and are avoided in spot tables.

Parts:
  A  checks: ghost case reproduces axion_aniso.pairs at three (c1, omega); c1 = 0 gives the frozen
     frequencies for every case; palindromic coefficients symmetric to roundoff; M_h = M cancels (rate 0
     off resonance at c1 = 0.01); heavy mode alone against Mathieu (rate c1/(4M) at omega = 2M).
  B  frozen-Q map, omega in [0.25, 4.5], c1 = 0.01, for every case: unstable intervals; spot rates against
     the pole-term law.  (The frozen-Q map is the pole terms only; on 15 Oct the "full law" of Part D was
     taken as the physics -- withdrawn, see above.)
  C  the full law per Hubble time at the end of inflation, Fujita's benchmark, along r = M/m_phi for every
     case, with the sign rule (as a statement about the pole term; phi_amp = 1 here, replaced 16 Oct by (A11)).
  D  [withdrawn 15 Oct, see above] direct two-frequency integration at m_phi/gQ = 100, phi = 2, 30 Q cycles: slope per cycle vs the slow
     equation, for the healthy twin at omega = 0.8, 0.6 (tachyonic side), 1.5, 2 (stable side); ghost +
     healthy at M_h = 2 for omega = 0.8 (stable), 1.5 (window), 3 (closed); M_h = 1 at omega = 0.8, 1.5
     (marginal); the ghost at omega = 2 under the same conditions as reference.

Float only (DOP853, rtol 1e-11); nothing here is certified.  Usage (from the repository root, ~12 min on the
cloud CPU):
    python axion_su2/code/axion_twins.py > axion_su2/results/axion_twins.txt
Added 2026-10-15.
"""
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.special import ellipk

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import axion_aniso as ax                                     # noqa: E402

RT = dict(method='DOP853', rtol=1e-11, atol=1e-13)
T0 = time.time()
CASES = {
    'ghost':            [(-1.0, 1.0)],
    'healthy twin':     [(+1.0, 1.0)],
    'ghost+healthy(2)': [(-1.0, 1.0), (+1.0, 2.0)],
    'ghost+healthy(1/2)': [(-1.0, 1.0), (+1.0, 0.5)],
    'ghost+healthy(1)': [(-1.0, 1.0), (+1.0, 1.0)],
}


def stamp(msg=''):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ------------------------------------------------------------------ the linear system
def kmat(c, modes):
    """First-order matrix for y = (s, s', h_1, h_1', ...), s'' = -c (s + sum h), h_i'' = -M_i^2 h_i - sigma_i c (s + sum h)."""
    n = 1 + len(modes)
    K = np.zeros((2 * n, 2 * n))
    K[0, 1] = 1.0
    K[1, 0::2] = -c
    for i, (sig, M) in enumerate(modes):
        r = 2 * (i + 1)
        K[r, r + 1] = 1.0
        K[r + 1, 0::2] = -sig * c
        K[r + 1, r] += -M * M
    return K


def monodromy(c1, omega, modes, c0=0.0):
    T = 2.0 * np.pi / omega
    n2 = 2 * (1 + len(modes))

    def rhs(t, y):
        K = kmat(c0 + c1 * np.cos(omega * t), modes)
        return (K @ y.reshape(n2, n2)).ravel()
    y = solve_ivp(rhs, (0.0, T), np.eye(n2).ravel(), **RT).y[:, -1]
    return y.reshape(n2, n2), T


def char_coeffs(A):
    """Faddeev-LeVerrier: coefficients a[0..n] of det(lambda I - A) = sum a[k] lambda^k."""
    n = A.shape[0]
    a = np.zeros(n + 1)
    a[n] = 1.0
    Mk = np.zeros_like(A)
    for k in range(1, n + 1):
        Mk = A @ Mk + a[n - k + 1] * np.eye(n)
        a[n - k] = -np.trace(A @ Mk) / k
    return a


def pairs_general(Phi, T):
    """Reciprocal pairs of a 2n x 2n symplectic monodromy via x = lambda + 1/lambda.  Returns
    ([(x, rate)] sorted by |x - 2| (massless pair first), palindromic asymmetry)."""
    a = char_coeffs(Phi)
    n2 = len(a) - 1
    n = n2 // 2
    asym = np.max(np.abs(a - a[::-1])) / np.max(np.abs(a))
    a = 0.5 * (a + a[::-1])
    # q(x) = a_n + sum_{m=1}^n a_{n+m} p_m(x),  p_m = lambda^m + lambda^-m,  p_0 = 2, p_1 = x, p_m = x p_{m-1} - p_{m-2}
    p_prev, p_cur = np.array([2.0]), np.array([1.0, 0.0])          # polynomials in x, highest power first
    q = np.array([a[n]])
    for m in range(1, n + 1):
        pm = p_cur if m == 1 else np.polysub(np.polymul([1.0, 0.0], p_cur), p_prev)
        if m > 1:
            p_prev, p_cur = p_cur, pm
        q = np.polyadd(q, a[n + m] * pm)
    xs = np.roots(q)
    out = []
    for x in xs:
        lam = 0.5 * (x + np.sqrt(x * x - 4.0 + 0j))
        out.append((x, abs(np.log(abs(lam))) / T))
    out.sort(key=lambda p: abs(p[0] - 2.0))
    return out, asym


def rate_general(c1, omega, modes, which='max'):
    Phi, T = monodromy(c1, omega, modes)
    pr, _ = pairs_general(Phi, T)
    return max(r for _, r in pr) if which == 'max' else pr[0][1]


def D_of(w, modes):
    """D(w) = sum_i sigma_i/(M_i^2 - w); inf at a pole (the law is not defined there)."""
    return sum(sig / (M * M - w) if M * M != w else np.inf for sig, M in modes)


def law_full(phiX, omega, modes):
    """Full law: s'' = + (phi X omega)^2/2 * D(omega^2) * s; rate if tachyonic else 0.
    Withdrawn 15 Oct: holds on the M -> infinity background only (Part D); on the dressed one the terms cancel."""
    D = D_of(omega ** 2, modes)
    A = 0.5 * (phiX * omega) ** 2 * D if np.isfinite(D) else np.nan
    return np.sqrt(A) if A > 0 else 0.0


def law_pole(c1, omega, modes):
    """Frozen-Q (pole terms only): s'' = -(c1^2/2) [1/omega^2 - D(omega^2)] s; rate if tachyonic else 0."""
    D = D_of(omega ** 2, modes)
    A = -0.5 * c1 ** 2 * (1.0 / omega ** 2 - D) if np.isfinite(D) else np.nan
    return np.sqrt(A) if A > 0 else 0.0


# ------------------------------------------------------------------ Part A
def part_A():
    print("\nPart A: checks\n")
    print("A1  ghost case against axion_aniso.pairs (massless-pair x and rate)")
    worst = 0.0
    for c1, om in ((0.01, 1.5), (0.1, 2.5), (0.03, 0.7)):
        Phi, T = monodromy(c1, om, CASES['ghost'])
        pr, asym = pairs_general(Phi, T)
        Phi2, _ = ax.pu_monodromy(lambda t, _x: c1 * np.cos(om * t), T)
        pr2 = ax.pairs(Phi2, T)
        d = max(abs(pr[0][0] - pr2[0][0]), abs(pr[0][1] - pr2[0][1]))
        worst = max(worst, d)
        print(f"    c1 = {c1}, omega = {om}: x = {pr[0][0].real:.10f} vs {pr2[0][0].real:.10f}, rate {pr[0][1]:.6e} vs {pr2[0][1]:.6e}, "
              f"palindromic asymmetry {asym:.1e}")
    print(f"    {'PASS' if worst < 1e-7 else 'FAIL'} ({worst:.1e})")

    print("\nA2  c1 = 0: pairs at the frozen frequencies exp(+-i M_i T) and the massless x = 2, every case (omega = 1.3)")
    worst = 0.0
    for name, modes in CASES.items():
        Phi, T = monodromy(0.0, 1.3, modes)
        pr, asym = pairs_general(Phi, T)
        xs = sorted([p[0].real for p in pr])
        ref = sorted([2.0] + [2 * np.cos(M * T) for _, M in modes])
        d = max(abs(np.array(xs) - np.array(ref)))
        worst = max(worst, d)
        print(f"    {name:20s}: x = {np.round(xs, 8)}  ref {np.round(ref, 8)}  asym {asym:.1e}")
    print(f"    {'PASS' if worst < 1e-7 else 'FAIL'} ({worst:.1e})")

    print("\nA3  M_h = M (ghost + healthy at the same mass): pole terms cancel, only the Kapitza term is left (stable)")
    for om in (0.7, 1.5, 2.5, 3.3):
        r = rate_general(0.01, om, CASES['ghost+healthy(1)'])
        print(f"    omega = {om}: rate {r:.2e}  (law_pole {law_pole(0.01, om, CASES['ghost+healthy(1)']):.2e}, D = {D_of(om**2, CASES['ghost+healthy(1)']):.1e})")

    print("\nA4  heavy mode alone, h'' + (M^2 + sigma c1 cos omega t) h = 0 (2-d Mathieu): rate at omega = 2M vs c1/(4M); off the tongue")
    for sig in (-1.0, +1.0):
        for om, lab in ((2.0, 'omega = 2M (first tongue)'), (1.7, 'omega = 1.7M'), (0.9, 'omega = 0.9M')):
            T = 2 * np.pi / om
            Mh = ax.hill_monodromy(1.0, sig * 0.01, om)
            lam = np.linalg.eigvals(Mh)
            r = np.max(np.log(np.abs(lam))) / T
            print(f"    sigma = {sig:+.0f}, {lab:26s}: rate {r:.3e}   (c1/4 = {0.0025:.3e})")
    print("    narrow tongues only: a heavy mode by itself, of either sign, is not destabilised by the drive away from omega = 2M/n")


# ------------------------------------------------------------------ Part B
def part_B():
    print("\nPart B: frozen-Q map, c1 = 0.01, omega in [0.25, 4.5] (step 2.5e-3); unstable intervals (rate > 1e-6) and spot rates\n")
    grid = np.linspace(0.25, 4.5, 1701)
    spots = (0.5, 0.8, 0.9, 1.5, 2.5, 3.5)
    for name, modes in CASES.items():
        rates = np.array([rate_general(0.01, om, modes) for om in grid])
        iv = ax.intervals(grid, rates > 1e-6)
        print(f"  {name}: D(w) = " + " + ".join(f"({sig:+.0f})/({M*M:g} - w)" for sig, M in modes) + f";  {len(iv)} interval(s)")
        for lo, hi, n in iv:
            sel = (grid >= lo) & (grid <= hi)
            k = np.argmax(rates[sel])
            print(f"     omega in [{lo:.4f}, {hi:.4f}]  width {hi - lo:.4f}  max rate {rates[sel][k]:.3e} at {grid[sel][k]:.4f}  ({n} pts)")
        print("     spot: " + "  ".join(f"omega={om}: {rates[np.argmin(abs(grid - om))]:.2e}/{law_pole(0.01, om, modes):.2e}" for om in spots)
              + "   (measured/pole law)")
        stamp(name)
    print("  Frozen Q carries only the pole terms: for the healthy twin they are tachyonic in 1/sqrt2 < omega < 1 and stable")
    print("  elsewhere; the jitter term (Part D) turns that into tachyonic for all omega < 1.")


# ------------------------------------------------------------------ Part C
def part_C(phi_amp=1.0):
    mQ = ax.GEV['mQ']
    print("\nPart C: the full law per Hubble time at the end of inflation (Fujita's benchmark, m_Q,end = 5, phi_amp = 1):")
    print("        mu/H_end = m_Q phi sqrt(D(m_phi^2)/2) / r,  r = M/m_phi, D in units of M = 1 (omega = 1/r); 0 = stable\n")
    rs = (0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 1.1, 1.5, 1.8, 2.5, 4.0)
    print(f"    {'r = M/m_phi':>12} | " + " ".join(f"{name:>18s}" for name in CASES))
    for r in rs:
        om = 1.0 / r
        vals = []
        for modes in CASES.values():
            D = D_of(om ** 2, modes)
            vals.append(mQ * phi_amp * np.sqrt(D / 2) / r if D > 0 else 0.0)
        print(f"    {r:12.2f} | " + " ".join(f"{v:18.3f}" for v in vals))
    print("    sign rule: tachyonic iff D(m_phi^2) = sum_i sigma_i/(M_i^2 - m_phi^2) > 0.  ghost: m_phi > M.  healthy twin: m_phi < M.")
    print("    ghost + healthy at 2M (the sixth-order operator): M < m_phi < 2M.  at M/2: m_phi < M/2 or m_phi > M.  at M: never.")
    print("    (m_phi/gQ ~ 1e5 at the end of inflation, so these are flat-space rates at the frozen end-of-inflation Q,")
    print("     before the Q-cycle average <|cn|> = 0.599 and the turning-point factor 1.12 of axion_aniso.py C5.)")


# ------------------------------------------------------------------ Part D
def part_D():
    print("\nPart D: direct two-frequency integration, m_phi/gQ = 100, phi = 2, 30 Q periods; slope of ln sigma_max per Q cycle")
    print("        (cycles 10-30) against the slow equation s'' = mu(t)^2 s with mu from the full law along X(t)\n")
    K = ellipk(0.5)
    phi_amp, ratio = 2.0, 100
    runs = (('ghost', 2.0), ('healthy twin', 0.8), ('healthy twin', 0.6), ('healthy twin', 1.5), ('healthy twin', 2.0),
            ('ghost+healthy(2)', 0.8), ('ghost+healthy(2)', 1.5), ('ghost+healthy(2)', 3.0),
            ('ghost+healthy(1)', 0.8), ('ghost+healthy(1)', 1.5))
    print(f"  {'case':>18} {'omega':>6} {'T_Q meas':>9} | {'slope':>7} {'slow eq':>8} {'full int':>8} | {'D(w)':>8} {'ln s_max end':>12}")
    for name, omega in runs:
        modes = CASES[name]
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
        T_meas = np.mean(np.diff(zc))
        T_run = n_cyc * T_meas
        n2 = 2 * (1 + len(modes))

        def rhs(t, y):
            xb = y[:2]
            Kc = kmat(kappa(t) * xb[0], modes)
            return np.concatenate([bg(t, xb), (Kc @ y[2:].reshape(n2, n2)).ravel()])
        y0 = np.concatenate([[Q_amp, 0.0], np.eye(n2).ravel()])
        t_eval = np.arange(1, n_cyc + 1) * T_meas
        sol = solve_ivp(rhs, (0.0, T_run), y0, t_eval=t_eval, **RT)
        lns = np.array([np.log(np.linalg.svd(sol.y[2:, k].reshape(n2, n2), compute_uv=False)[0]) for k in range(n_cyc)])
        slope = np.polyfit(np.arange(10, n_cyc + 1), lns[9:], 1)[0]
        # slow amplitude and the slow equation
        Td = 2 * np.pi / omega
        tt2 = np.linspace(0, T_meas, 4001)
        X_t = np.array([np.mean(sb.sol(np.linspace(t - Td / 2, t + Td / 2, 41))[0]) for t in tt2])
        mu_t = np.array([law_full(phi_amp * abs(x), omega, modes) for x in X_t])
        full = np.trapezoid(mu_t, tt2)
        mu_interp = interp1d(tt2, mu_t, kind='linear')

        def slow_rhs(t, y):
            m2 = mu_interp(t % T_meas) ** 2
            return [y[1], m2 * y[0], y[3], m2 * y[2]]
        ss = solve_ivp(slow_rhs, (0.0, T_run), [1, 0, 0, 1], t_eval=t_eval, **RT)
        lns2 = [np.log(np.linalg.svd(np.array([[ss.y[0, k], ss.y[2, k]], [ss.y[1, k], ss.y[3, k]]]), compute_uv=False)[0])
                for k in range(n_cyc)]
        slow = np.polyfit(np.arange(10, n_cyc + 1), lns2[9:], 1)[0]
        print(f"  {name:>18} {omega:6.2f} {T_meas:9.2f} | {slope:7.3f} {slow:8.3f} {full:8.3f} | {D_of(omega**2, modes):8.3f} {lns[-1]:12.2f}",
              flush=True)
    print("\n  A slope at the drift level (~0.05 per cycle, the ln t of a free direction over cycles 10-30) is 'stable'.")
    print("  The ghost row repeats axion_aniso.py Part D at the new ratio (there 3.286 at ratio 200).")


if __name__ == "__main__":
    print(f"axion_twins.py  python {sys.version.split()[0]}  numpy {np.__version__}  ({sys.executable})")
    print("units g = M = 1; massless direction + heavy modes (sigma_i, M_i), coupling -c(t)(s + sum h)^2/2; DOP853 rtol 1e-11")
    part_A()
    stamp("Part A done")
    part_B()
    stamp("Part B done")
    part_C()
    stamp("Part C done")
    part_D()
    stamp("Part D done")
    print(f"\n({time.time() - T0:.0f} s)")
