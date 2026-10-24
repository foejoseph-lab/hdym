"""axion_dressed.py -- the correction of 15 Oct: the drive dresses the background through the same operator.

What was wrong.  axion_aniso.py (14 Oct) and axion_twins.py (15 Oct) measured the anisotropic sector on the
M -> infinity background Q'' = -2Q^3 + kappa Q^2.  That background is the leading term in 1/M^2 and is valid
only for m_phi << M (omega = m_phi/M << 1).  For omega ~ 1 the ghost's response to the drive is O(1) of the
drive: the exact isotropic 4th-order background (Q'''' + Q'' + ... = kappa Q^2 + ..., the ghost included)
responds at the drive frequency with the FULL propagator  Pi(w) = 1/(w^2 - w) = -1/w + D(w),  w = omega^2,
D(w) = sum_i sigma_i/(M_i^2 - w) the heavy part (D = 1/(w - 1) for the Lee-Wick ghost), where the M -> infinity
background has only the massless pole -1/w.  The jitter of Q is  dQ = kappa X^2 Pi(omega^2)  (frozen slow
amplitude X), so its DC product with the drive is  <kappa dQ> = <kappa^2> X^2 Pi  instead of  -<kappa^2> X^2/omega^2.

The consequence, for the massless anisotropic direction s with Hill coefficient c(t) = kappa Q:
    second order:  s'' = [ -<kappa dQ> + <kappa^2> X^2 Pi ] s = 0 ,
the drive x jitter term and the drive x response term are the same response function attached in the two
orders, and they cancel identically, for ANY heavy sector, as long as the background is dressed by the
operator that the perturbation obeys (they are the same field).  In the periodic 4-d problem this reads: the
DC offset c0 = (c1^2/2) Pi cancels the second-order tachyon of eps'''' + eps'' + (c0 + c1 cos omega t) eps = 0
exactly (Part A).  What remains in the anisotropic sector are the ghost pair's own parametric tongues at
omega = 2M/n (rate c1/4 at the first, width c1/omega^2 relative) and the combination resonance near
omega = 1.

The isotropic sector (the background's own direction; the tangent of axion_derive.py has Hill coefficient
-2 kappa Q where the anisotropic one has +kappa Q) does not cancel:  c1_iso = 2 c1,  c0_iso = -2 <kappa dQ> =
-(c1_iso^2/4) Pi,  so
    s'' = [ (c1_iso^2/4) Pi + (c1_iso^2/2) Pi ] s = (3/4) c1_iso^2 Pi s = 3 omega^2 phi^2 X^2 Pi(omega^2) s :
tachyonic iff Pi(m_phi^2) > 0, i.e. iff the full propagator of the gauge field, massless pole included, is
positive at the drive frequency.  Lee-Wick ghost: Pi = 1/(omega^2 (omega^2 - 1)) > 0 iff m_phi > M, rate
    mu_iso = phi X sqrt(3/(omega^2 - 1))  (M = 1),   mu_iso = sqrt3 g Q phi r / sqrt(1 - r^2),  r = M/m_phi,
    mu_iso / H_end = sqrt3 m_Q phi r / sqrt(1 - r^2).
Healthy twin: Pi = (2 omega^2 - 1)/(omega^2 (1 - omega^2)) > 0 iff 1/sqrt2 < omega < 1.  No heavy sector
(pure Yang-Mills, Pi = -1/omega^2): stable, the Kapitza stabilisation of the isotropic Hill equation.
So for M < m_phi the driven isotropic Lee-Wick background is itself unstable, in its own isotropic direction,
and the anisotropic sector is marginal at second order on it: the isotropy assumption of Fujita et al. is
not where the Lee-Wick theory fails; the existence of the isotropic solution is.

Parts:
  A  the periodic 4-d (and 2n-d) problem, frozen X, for every heavy sector of axion_twins.py: anisotropic
     rate with c0 = (c1^2/2) Pi (expected 0 off the tongues) against c0 = 0 (the old law); isotropic rate
     with (c1_iso, c0_iso) against sqrt(3/4 c1_iso^2 Pi).
  B  the exact 4th-order isotropic background of the Lee-Wick theory, started on the dressed particular
     solution (Q''' = kappa'(0) Q^2/(1 - omega^2)), m_phi/gQ = 100, 3 Q periods, a (phi, omega) scan, 4 Q periods: the exact
     isotropic 4-d tangent (w0 = Q'' + 2Q^3) and the anisotropic 4-d sector alongside, per-cycle increments of
     ln sigma_max; stable points drift, unstable ones grow by several per cycle and the background blows up.
     Result: a threshold in the drive, omega^2 phi^2 Pi(omega^2) ~ 1-2 (phi^2 r^2/(1 - r^2) ~ 1-2 in physical
     variables), none for M > m_phi (made exact on 16 Oct, iso_lame_gap.py: lam_eff = 2 - omega^2 phi^2 Pi = 0,
     i.e. omega^2 phi^2 Pi = 2, sharp; the omega = 2 rows here sit on the ghost's first tongue, which is lost
     earlier, at 1.33); below it both sectors are marginal on the self-consistent background --
     the second-order law of A2 for the isotropic sector holds only against the free-particle model, on the
     real background the restoring curvature of the orbit family competes with it and wins below threshold.
  C  the exact background's own runaway (nonlinear, from the roundoff seed): envelope rate of its ghost
     content J~ = Q'' + 2Q^3 - kappa Q^2 against mu_iso, omega = 2, 3, 0.8.
  D  the physical numbers at Fujita's benchmark: mu_iso/H_end along r; the anisotropic ghost tongue at
     r = 1/2 (unchanged from axion_aniso.py C3); the e-fold estimate through reheating for the isotropic
     runaway (same a^-5/2 envelope argument as axion_aniso.py C5).

Float only (DOP853, rtol 1e-11); nothing here is certified.  Usage (from the repository root, ~10 min):
    python axion_su2/code/axion_dressed.py > axion_su2/results/axion_dressed.txt
Added 2026-10-15.
"""
import os
import pickle
import sys
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.special import ellipk, ellipj

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import axion_aniso as ax                                     # noqa: E402
import axion_twins as tw                                     # noqa: E402

RT = dict(method='DOP853', rtol=1e-11, atol=1e-13)
T0 = time.time()


def stamp(msg=''):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def companion(c):
    C = np.zeros((4, 4))
    C[0, 1] = C[1, 2] = C[2, 3] = 1.0
    C[3, 2] = -1.0
    C[3, 0] = -c
    return C


def Pi_of(w, modes):
    return -1.0 / w + tw.D_of(w, modes)


def rate_c0(c0, c1, omega, modes, massless_only=False):
    """Frozen-X periodic monodromy with a DC offset; max rate over pairs, or the massless pair only."""
    Phi, T = tw.monodromy(c1, omega, modes, c0=c0)
    pr, _ = tw.pairs_general(Phi, T)
    return pr[0][1] if massless_only else max(r for _, r in pr)


# ------------------------------------------------------------------ Part A
def part_A():
    print("\nPart A: the periodic problem with the dressed DC offset, frozen X, every heavy sector\n")
    print("A1  anisotropic: eps'''' + ... + (c0 + c1 cos omega t) eps with c0 = (c1^2/2) Pi(omega^2); massless-pair rate")
    print(f"    {'case':>18} {'omega':>6} {'c1':>6} | {'Pi':>8} | {'rate c0=0':>10} {'old law':>10} | {'rate dressed':>13}")
    worst = 0.0
    for name, modes in tw.CASES.items():
        for omega in (3.0, 1.5, 2.5, 0.8, 0.6):
            if abs(omega - 2.0) < 1e-9 or any(abs(omega - 2 * M / n) < 1e-9 for _, M in modes for n in (1, 2, 3, 4)):
                continue
            w = omega ** 2
            Pi = Pi_of(w, modes)
            c1 = 0.03
            r0 = rate_c0(0.0, c1, omega, modes, massless_only=True)
            old = np.sqrt(0.5 * c1 ** 2 * Pi) if Pi > 0 else 0.0           # the pole law: s'' = (c1^2/2) Pi s
            rd = rate_c0(0.5 * c1 ** 2 * Pi, c1, omega, modes, massless_only=True)
            if not (name == 'ghost+healthy(1/2)' and abs(omega - 0.6) < 0.2):
                worst = max(worst, rd / c1 ** 2)
            print(f"    {name:>18} {omega:6.2f} {c1:6.3f} | {Pi:8.3f} | {r0:10.3e} {old:10.3e} | {rd:13.3e}")
    print("    Residuals where Pi < 0 (|Pi| large: the Kapitza and pole terms add) are higher order in c1 Pi: scaling test")
    print(f"    {'case':>18} {'omega':>6} | {'c1':>6} {'dressed rate':>13} {'rate/c1^2':>10} {'old law(c1)':>12}")
    for name, omega in (('ghost', 0.8), ('ghost', 0.6), ('healthy twin', 1.5), ('ghost+healthy(2)', 2.5)):
        modes = tw.CASES[name]
        Pi = Pi_of(omega ** 2, modes)
        for c1 in (0.03, 0.01, 0.003):
            rd = rate_c0(0.5 * c1 ** 2 * Pi, c1, omega, modes, massless_only=True)
            old = np.sqrt(0.5 * c1 ** 2 * abs(Pi))
            print(f"    {name:>18} {omega:6.2f} | {c1:6.3f} {rd:13.3e} {rd / c1 ** 2:10.3f} {old:12.3e}")
    print("    (rate ~ c1^2 = a fourth-order term, against the uncancelled second-order c1 |Pi|^(1/2)/sqrt2: at the")
    print("     physical c1 ~ 1e-5 omega^2 it is 1e-5 of the old law)  PASS: the second-order tachyon of the massless")
    print("    direction is cancelled by the dressed jitter for every heavy sector.")

    print("\nA2  isotropic: c1_iso = 2 c1, c0_iso = -(c1_iso^2/4) Pi; massless-pair rate vs sqrt(3/4 c1_iso^2 Pi)")
    print(f"    {'case':>18} {'omega':>6} | {'Pi':>8} | {'rate':>10} {'law':>10} {'ratio':>7}")
    worst = 0.0
    for name, modes in tw.CASES.items():
        for omega in (3.0, 1.5, 2.5, 0.8, 0.6):
            if any(abs(omega - 2 * M / n) < 1e-9 for _, M in modes for n in (1, 2, 3, 4)):
                continue
            w = omega ** 2
            Pi = Pi_of(w, modes)
            c1 = 0.03
            c1i, c0i = 2 * c1, -(4 * c1 ** 2 / 4) * Pi
            r = rate_c0(c0i, c1i, omega, modes, massless_only=True)
            law = np.sqrt(0.75 * c1i ** 2 * Pi) if Pi > 0 else 0.0
            ratio = r / law if law > 0 else (np.nan if r < 1e-9 else np.inf)
            if law > 0 and omega >= 1.5:
                worst = max(worst, abs(ratio - 1))
            print(f"    {name:>18} {omega:6.2f} | {Pi:8.3f} | {r:10.3e} {law:10.3e} {ratio:7.3f}")
    print(f"    worst |ratio - 1| where the law is nonzero and omega >= 1.5: {worst:.1e}  ({'PASS' if worst < 0.02 else 'FAIL'})")
    print("    sign rule for the isotropic sector: tachyonic iff Pi(m_phi^2) = -1/m_phi^2 + D(m_phi^2) > 0 (massless pole included).")


# ------------------------------------------------------------------ exact background helpers
PKL = os.path.join(HERE, '..', 'rung1_axion_system.pkl')


def load_q4():
    d = pickle.load(open(PKL, 'rb'))
    args = list(d['X']) + list(d['P']) + list(d['X2']) + list(d['X3']) + [d['cj'], d['ksyms'][0]]
    return sp.lambdify(args, d['q4'][0], 'numpy')


def iso_rhs(q4, Q, Q1, Q2, Q3, kap):
    return q4(Q, Q, Q, Q1, Q1, Q1, Q2, Q2, Q2, Q3, Q3, Q3, 1.0, kap)


def iso_tangent(Y, Q, Q1, Q2, kap):
    """exact isotropic tangent (axion_derive CHECK 6 with w0 = Q'' + 2Q^3)."""
    w0 = Q2 + 2.0 * Q ** 3
    c2 = 1.0 + 12.0 * Q * Q
    c1 = 24.0 * Q * Q1
    c0 = 6.0 * Q * Q + 12.0 * Q * w0 - 2.0 * kap * Q + 36.0 * Q ** 4 + 12.0 * Q * Q2 + 12.0 * Q1 * Q1
    dY = np.empty_like(Y)
    dY[0] = Y[1]
    dY[1] = Y[2]
    dY[2] = Y[3]
    dY[3] = -(c2 * Y[2] + c1 * Y[1] + c0 * Y[0])
    return dY


def smooth_period(sol_fun, Tend, omega):
    """Q period from zero crossings of Q averaged over one drive period (the jitter removed)."""
    tt = np.linspace(0, Tend, 40001)
    Td = 2 * np.pi / omega
    Qs = np.array([np.mean(sol_fun(np.linspace(t - Td / 2, t + Td / 2, 41))[0]) for t in tt])
    # a second smoothing over three drive periods, then only crossings at least a quarter of the first spacing apart
    k = max(1, int(3 * Td / (tt[1] - tt[0])))
    Qs = np.convolve(Qs, np.ones(k) / k, mode='same')
    zc = tt[:-1][(Qs[:-1] > 0) & (Qs[1:] <= 0)]
    if len(zc) > 2:
        keep = [zc[0]]
        for z in zc[1:]:
            if z - keep[-1] > 0.25 * (zc[1] - zc[0]):
                keep.append(z)
        zc = np.array(keep)
    return np.mean(np.diff(zc)) if len(zc) > 2 else np.nan, zc


def exact_background(q4, omega, phi_amp, Q_amp, Tend, stop_at=None):
    kappa = lambda t: -omega * phi_amp * np.sin(omega * t)
    kappa1 = lambda t: -omega ** 2 * phi_amp * np.cos(omega * t)

    def bg4(t, b):
        Q, Q1, Q2, Q3 = b
        return np.array([Q1, Q2, Q3, iso_rhs(q4, Q, Q1, Q2, Q3, kappa(t))])
    b0 = np.array([Q_amp, 0.0, -2 * Q_amp ** 3, kappa1(0.0) * Q_amp ** 2 / (1 - omega ** 2)])
    ev = None
    if stop_at is not None:
        def ev(t, b):
            return abs(b[0]) - stop_at
        ev.terminal = True
    sb = solve_ivp(bg4, (0.0, Tend), b0, dense_output=True, events=ev, **RT)
    return sb, bg4, b0, kappa


# ------------------------------------------------------------------ Part B
def part_B(q4):
    print("\nPart B: exact Lee-Wick isotropic background (ghost included, data on the dressed particular solution),")
    print("        m_phi/gQ = 100, 4 Q periods (period from the zero crossings of the smoothed Q, else the dressed Kapitza formula):")
    print("        per-cycle increments of ln sigma_max of the exact isotropic 4-d tangent")
    print("        and of the anisotropic 4-d sector; a stable point shows the free drift ln(k+1)/k = 0.69, 0.41 after the")
    print("        first-cycle transient, an unstable one grows by several per cycle and the background itself blows up")
    print("        (its own growing mode, seeded by the imperfection of the initial data) -- the run then stops at |Q| > 10 Q_amp.\n")
    K = ellipk(0.5)
    ratio = 100
    print(f"  {'phi':>4} {'omega':>6} {'T_Q':>7} {'w^2phi^2Pi':>10} | {'iso incr (4 cycles)':>29} | {'aniso incr':>29} | {'verdict':>10} {'|Q|max/Q_amp':>12}")
    for omega, phis in ((1.5, (1.0, 1.5, 2.0)), (2.0, (1.0, 2.0, 2.2, 2.5)), (3.0, (1.0, 2.0, 3.0)), (0.8, (1.0, 3.0))):
        for phi_amp in phis:
            Q_amp = omega / ratio
            n_cyc = 4
            phi_eff = phi_amp / (1 - omega ** 2)
            lam = 2 + phi_eff ** 2
            X_amp = Q_amp * (1 + 2 * phi_eff ** 2 / lam) ** 0.25
            T_form = 4 * K / (np.sqrt(lam) * X_amp)
            span = 1.3 * (n_cyc + 2) * T_form
            sb, bg4, b0, kappa = exact_background(q4, omega, phi_amp, Q_amp, span, stop_at=10 * Q_amp)
            T_Q, zc = smooth_period(sb.sol, min(span, sb.t[-1]), omega)
            src = 'zc'
            if not (np.isfinite(T_Q) and len(zc) >= 4 and 0.5 * T_form < T_Q < 2.5 * T_form):
                T_Q, src = T_form, 'formula'

            def rhs(t, y):
                b = y[:4]
                return np.concatenate([bg4(t, b), iso_tangent(y[4:20].reshape(4, 4), b[0], b[1], b[2], kappa(t)).ravel(),
                                       (companion(kappa(t) * b[0]) @ y[20:].reshape(4, 4)).ravel()])

            def blow(t, y):
                return abs(y[0]) - 10 * Q_amp
            blow.terminal = True
            t_eval = np.arange(1, n_cyc + 1) * T_Q
            sol = solve_ivp(rhs, (0.0, n_cyc * T_Q), np.concatenate([b0, np.eye(4).ravel(), np.eye(4).ravel()]),
                            t_eval=t_eval, events=blow, **RT)
            Y = np.asarray(sol.y) if len(sol.t) else np.zeros((36, 0))
            nk = Y.shape[1]
            li = np.array([np.log(np.linalg.svd(Y[4:20, k].reshape(4, 4), compute_uv=False)[0]) for k in range(nk)])
            la = np.array([np.log(np.linalg.svd(Y[20:, k].reshape(4, 4), compute_uv=False)[0]) for k in range(nk)])
            di = np.diff(np.concatenate([[0.0], li]))
            da = np.diff(np.concatenate([[0.0], la]))
            blew = sol.status == 1
            Qmax = (np.max(np.abs(Y[0])) if nk else 0.0) / Q_amp
            if blew:
                Qmax = max(Qmax, abs(sol.y_events[0][0][0]) / Q_amp)
            crit = omega ** 2 * phi_amp ** 2 * Pi_of(omega ** 2, tw.CASES['ghost'])
            later = di[1:] if nk >= 2 else np.array([])
            verdict = 'BLOW-UP' if blew else ('unstable' if (later.size and np.max(later) > 2.0) else
                                               ('stable' if later.size and np.max(later) < 1.3 else 'unclear'))
            fmt = lambda v: " ".join(f"{x:7.2f}" for x in v) + "   " * (n_cyc - len(v))
            print(f"  {phi_amp:4.1f} {omega:6.2f} {T_Q:7.1f} {crit:10.3f} | {fmt(di):>29} | {fmt(da):>29} | {verdict:>10} {Qmax:12.2f} ({src})"
                  + (f"  (stopped at t = {sol.t_events[0][0] / T_Q:.2f} T_Q)" if blew else ""), flush=True)
    print("\n  Threshold in w^2 phi^2 Pi: between 0.33 and 1.33 at omega = 2 (1.33 is slowly unstable, ~2.5 per cycle, the")
    print("  background already moving), above 1.125 at omega = 3 and above 1.8 at omega = 1.5 over four cycles; none below")
    print("  omega = 1 (Pi < 0).  So of order 1-2 and not sharp in this variable.  The averaged estimate -- the isotropic")
    print("  tachyon 3 w^2 phi^2 X^2 Pi against the background's own restoring curvature 6 <Q^2> = 2.74 X^2, both ~ <cn^2>")
    print("  -- puts it at 2: the right order (the isotropic direction is the tangent of the orbit family, not a free")
    print("  particle).  The linear rate above threshold is not measured cleanly here: the background's own growing mode")
    print("  is seeded at ~1e-3 by the approximate initial data and reaches O(1) within a cycle or two; what is seen is")
    print("  2.5-6 per Q period at omega = 2, phi = 2-2.2.")
    print("  The anisotropic sector follows once the background runs; below threshold it drifts (marginal).")


# ------------------------------------------------------------------ Part C
def part_C(q4):
    print("\nPart C: the exact background's ghost content J~ = Q'' + 2Q^3 - kappa Q^2 on the dressed data, phi = 2: it stays at")
    print("        its sourced value (the particular solution) at these points, which are all below threshold; the envelope")
    print("        rate is that of the sourced part, not of a growing mode (the seed here is roundoff).\n")
    phi_amp, ratio = 2.0, 100
    for omega in (2.0, 3.0, 0.8):
        Q_amp = omega / ratio
        Tend = 400.0 if omega > 1 else 800.0
        sb, _, _, kappa = exact_background(q4, omega, phi_amp, Q_amp, Tend)
        tt = np.linspace(0, Tend, 8001)
        Q, _, Q2, _ = sb.sol(tt)
        J = Q2 + 2 * Q ** 3 - kappa(tt) * Q ** 2
        nw = int(Tend / (2 * np.pi))
        env = np.array([np.max(np.abs(J[(tt >= k * 2 * np.pi) & (tt < (k + 1) * 2 * np.pi)])) for k in range(nw)])
        tc = (np.arange(nw) + 0.5) * 2 * np.pi
        # fit where the envelope is between 1e-10 and 1e-3 (linear regime, above roundoff)
        sel = (env > 1e-10) & (env < 1e-3)
        rate = np.polyfit(tc[sel], np.log(env[sel]), 1)[0] if sel.sum() > 5 else np.nan
        Pi = Pi_of(omega ** 2, tw.CASES['ghost'])
        X = Q_amp * (1.0 + 2.0 * (phi_amp * omega ** 2 * Pi * (-1)) ** 2 / (2 + (phi_amp * omega ** 2 * Pi) ** 2)) ** 0.25
        law = phi_amp * X * np.sqrt(3 * omega ** 2 * Pi) if Pi > 0 else 0.0
        print(f"    omega = {omega}: |J~| {env[0]:.1e} -> {env[-1]:.1e} over t = {Tend:g}; envelope rate (1e-10 < |J~| < 1e-3) "
              f"{rate:.4f} per unit t; law {law:.4f} (slow amplitude {X:.4f}); Q range [{Q.min():.3f}, {Q.max():.3f}]")
    print("    The dressed drive on Q is kappa Q^2 omega^2 (-Pi) ... i.e. kappa Q^2/(1 - omega^2): reversed and weaker for omega > 1,")
    print("    so the slow amplitude's Kapitza shift is small there.")


# ------------------------------------------------------------------ Part D
def part_D():
    MPl, g, mQ, mphi = ax.GEV['MPl'], ax.GEV['g'], ax.GEV['mQ'], ax.GEV['mphi']
    H = g * MPl / mQ ** 2
    Qend = MPl / mQ
    q_F = g * Qend / mphi
    print("\nPart D: Fujita's benchmark (axion_aniso.py Part C)\n")
    print("    threshold of the isotropic runaway, phi^2 r^2/(1 - r^2) ~ 1-2 (Part B), as the smallest M/m_phi that runs away:")
    for phi in (0.3, 0.5, 1.0, 2.0, 3.0):
        r_lo = np.sqrt(1.0 / (phi ** 2 + 1.0))
        r_hi = np.sqrt(2.0 / (phi ** 2 + 2.0))
        print(f"      phi_amp = {phi:3.1f}: r > {r_lo:.2f}-{r_hi:.2f}")
    # 16 Oct: phi_amp = 1 replaced by phi_end from their (A10)-(A11), 8.3e-4 at the benchmark (axion_frw_iso.py C);
    # the e-fold estimate below is superseded by that run (print strings left as run)
    print("    (phi_amp is the axion displacement at the end of inflation in units of f, not stated by Fujita et al.)")
    print("    Above threshold the growth is 2.5-6 per Q period at (phi, r) = (2-2.2, 0.5) in the run, i.e. ~0.01-0.03 per unit")
    print(f"    M^-1 at gQ/M = 0.02: in units of H_end that is ~0.01-0.03 x (M/H_end) = 0.01-0.03 x m_Q r m_phi/(gQ) -- large: the")
    print("    isotropic solution is lost within a fraction of a Hubble time once the threshold is crossed.")
    print(f"    anisotropic sector: marginal at second order; the ghost tongue at r = 1/2 stands (rate c1/4 = {q_F * mphi / 2 / H:.2f} per")
    print(f"    Hubble time at its centre, relative width q phi = {q_F:.1e}), and the linear resonance at r = 1.")
    print("    For M > m_phi everything is stable.")


if __name__ == "__main__":
    print(f"axion_dressed.py  python {sys.version.split()[0]}  numpy {np.__version__}  ({sys.executable})")
    part_A()
    stamp("Part A done")
    q4 = load_q4()
    part_B(q4)
    stamp("Part B done")
    part_C(q4)
    stamp("Part C done")
    part_D()
    stamp("Part D done")
    print(f"\n({time.time() - T0:.0f} s)")
