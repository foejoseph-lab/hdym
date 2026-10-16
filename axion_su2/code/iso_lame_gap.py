"""iso_lame_gap.py -- the Lamé-gap argument for the isotropic direction, and why it fails (README, correction section).

On the undriven orbit Q = A cn(sqrt2 A t | 1/2) of Q'' = -2 Q^3 the isotropic tangent obeys
    eps'' + 6 Q^2 eps = 0,
which in tau = sqrt2 A t is Lamé with n = 2, k^2 = 1/2, h = 3: the lower edge of the band [3, 4.5]
(the tangent of the orbit family: multiplier 1, Jordan block).  A constant negative shift of h lands
in the gap (1.5, 3) -> unstable for any shift.  That is the Lamé-gap argument.

The dressed second-order term is 3 kappa^2 Pi Q(t)^2 (axion_dressed.py A2), i.e. a Q^2-shaped term:
    eps'' + (6 - 3 kappa^2 Pi) Q^2 eps = 0  =  eps'' + 3 lam Q^2 eps,  lam = 2 - kappa^2 Pi,
which is exactly the tangent of the averaged family Q'' = -lam Q^3 (the Kapitza-averaged quartic).
Test the three versions numerically with the monodromy over one period:
  (a) constant shift:  eps'' + (6 Q^2 - delta <Q^2>) eps = 0, Q on the undriven orbit (gap argument);
  (b) Q^2 shift on the undriven orbit:  eps'' + (6 - delta) Q^2 eps = 0  (the inconsistent version);
  (c) Q^2 shift on the consistently shifted orbit: eps'' + 3 lam Q^2 eps = 0 on Q'' = -lam Q^3.
Report ln(max |multiplier|) per period for delta in a range.  Result: (a) and (b) unstable for every delta > 0,
(c) multiplier 1 for every lam > 0: the dressed term rescales the Q^2 coefficient and the averaged tangent stays the
tangent of the averaged orbit family (a Jordan block), so the threshold is lam_eff = 0, i.e. omega^2 phi^2 Pi = 2,
not "any positive coefficient".  threshold_offtongue() then checks this on the exact fourth-order background in the
(Q, w) form of axion_frw_iso.py, off the omega = 2 tongue (which axion_dressed.py Part B had sampled).

Float only.  Usage (from the repository root, ~2 min):
    python axion_su2/code/iso_lame_gap.py > axion_su2/results/iso_lame_gap.txt
Added 2026-10-16.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipk

RT = dict(method='DOP853', rtol=1e-11, atol=1e-13)
K = ellipk(0.5)


def period(lam, A):
    # Q'' = -lam Q^3, Q = A cn(sqrt(lam) A t | 1/2): period 4K/(sqrt(lam) A)
    return 4 * K / (np.sqrt(lam) * A)


def monodromy(lam_bg, coef, A=1.0, const=0.0):
    """Background Q'' = -lam_bg Q^3 from (A, 0); tangent eps'' + (coef Q^2 + const) eps = 0; one period."""
    T = period(lam_bg, A)

    def rhs(t, y):
        Q, P = y[0], y[1]
        E = y[2:].reshape(2, 2)
        dE = np.array([[0.0, 1.0], [-(coef * Q * Q + const), 0.0]]) @ E
        return np.concatenate([[P, -lam_bg * Q ** 3], dE.ravel()])
    sol = solve_ivp(rhs, (0, T), np.concatenate([[A, 0.0], np.eye(2).ravel()]), **RT)
    M = sol.y[2:, -1].reshape(2, 2)
    ev = np.linalg.eigvals(M)
    return np.max(np.log(np.abs(ev))), np.trace(M), T


def mean_Q2(lam, A):
    T = period(lam, A)
    sol = solve_ivp(lambda t, y: [y[1], -lam * y[0] ** 3], (0, T), [A, 0.0], dense_output=True, **RT)
    tt = np.linspace(0, T, 20001)
    return np.mean(sol.sol(tt)[0] ** 2)


def threshold_offtongue():
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import axion_frw_iso as fi
    print("\nexact Lee-Wick isotropic background, (Q, w) form, m_phi/gQ = 100, dressed data, 8 Kapitza periods")
    print("T = 4K/(sqrt(lam_eff) X) (the undriven period where lam_eff <= 0); stop at |Q| > 10 Q_amp.  omega = 2 is the")
    print("centre of the ghost pair's first parametric tongue (m_phi = 2M); 1.8 and 2.4 are off it.")
    print(f"  {'omega':>5} {'phi':>6} {'w2phi2Pi':>8} {'lam_eff':>7} {'T':>7} | |Q|max/Q_amp per cycle (8)" + " " * 20 + "| blow-up at")
    for omega in (1.8, 2.0, 2.4, 3.0):
        for crit in (1.0, 1.33, 1.7, 2.0, 2.3):
            w2 = omega ** 2
            P = 1.0 / (w2 * w2 - w2)
            phi = np.sqrt(crit / (w2 * P))
            Q_amp = omega / 100.0
            lam = 2.0 - crit
            phi_eff = -phi * w2 * P
            X = Q_amp * (1 + 2 * phi_eff ** 2 / (2 + phi_eff ** 2)) ** 0.25
            T = 4 * K / (np.sqrt(lam) * X) if lam > 0 else 4 * K / (np.sqrt(2) * X)
            rhs = fi.flat_rhs(omega, phi)
            y0 = np.array([Q_amp, 0.0, 0.0, w2 * phi * Q_amp ** 2 / (w2 - 1.0)])

            def blow(t, y):
                return abs(y[0]) - 10 * Q_amp
            blow.terminal = True
            sol = solve_ivp(rhs, (0.0, 8 * T), y0, dense_output=True, events=blow, **RT)
            per = []
            for k in range(8):
                a, b = k * T, min((k + 1) * T, sol.t[-1])
                per.append(np.max(np.abs(sol.sol(np.linspace(a, b, 4001))[0])) / Q_amp if b > a else np.nan)
            s = " ".join(f"{p:5.2f}" if np.isfinite(p) else "    -" for p in per)
            blew = sol.status == 1
            print(f"  {omega:5.2f} {phi:6.3f} {crit:8.2f} {lam:7.3f} {T:7.1f} | {s} | "
                  + (f"{sol.t_events[0][0] / T:.2f} T" if blew else "-"), flush=True)
    print("  Off the tongue the exact background is steady at omega^2 phi^2 Pi = 1.0, 1.33, 1.7 (its amplitude grows as lam_eff")
    print("  shrinks) and lost at 2.0 and 2.3: the averaged threshold.  On the tongue (omega = 2) it is lost already at 1.33.")


if __name__ == "__main__":
    A = 1.0
    q2 = mean_Q2(2.0, A)
    print(f"undriven orbit A = {A}: period {period(2.0, A):.6f}, <Q^2> = {q2:.6f} (= 0.457 A^2 expected)")
    print("\n(a) constant shift  eps'' + (6 Q^2 - delta <Q^2>) eps = 0 on the undriven orbit  [the gap argument]")
    print("(b) Q^2 shift on the undriven orbit  eps'' + (6 - delta) Q^2 eps = 0")
    print("(c) Q^2 shift on the consistent orbit  eps'' + 3 lam Q^2 eps = 0,  Q'' = -lam Q^3,  lam = 2 - delta/3")
    print(f"\n  {'delta':>6} {'w2phi2Pi':>8} | {'(a) rate/T':>10} {'tr':>8} | {'(b) rate/T':>10} {'tr':>8} | {'(c) rate/T':>10} {'tr':>8}")
    for delta in (0.0, 0.03, 0.1, 0.3, 1.0, 2.0, 3.0, 4.0, 5.0, 5.9):
        ra, tra, T = monodromy(2.0, 6.0, A, const=-delta * q2)
        rb, trb, _ = monodromy(2.0, 6.0 - delta, A)
        lam = 2.0 - delta / 3.0
        rc, trc, Tc = monodromy(lam, 3.0 * lam, A)
        # delta = 3 kappa^2 Pi  <->  w^2 phi^2 Pi = delta/3
        print(f"  {delta:6.2f} {delta / 3:8.3f} | {ra:10.4f} {tra:8.4f} | {rb:10.4f} {trb:8.4f} | {rc:10.4f} {trc:8.4f}")
    print("\n  tr = trace of the one-period monodromy (|tr| > 2 unstable, = 2 the band edge / Jordan block)")
    threshold_offtongue()
