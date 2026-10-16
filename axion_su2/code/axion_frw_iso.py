"""axion_frw_iso.py -- the isotropic threshold through reheating: does the Lee-Wick isotropic solution survive
the expansion once the drive dilutes?

Background (flat, axion_derive.py CHECK 6, auxiliary form; g = M = 1):
    Q'' = w - 2 Q^3,     w'' + (1 + 6 Q^2) w = kappa Q^2,     w = (Q'' + 2Q^3)/M^2 the ghost content,
which is the fourth-order isotropic equation of results/axion_derive.txt (checked here against the pickle).
FRW in conformal time with comoving Q (axion_derive.py CHECK 7, stage 2): M -> M a(eta) and the drive
(g/f)(dphi/deta) det A with comoving fields, i.e.
    Q'' = a^2 w - 2 Q^3,     w'' + (a^2 + 6 Q^2) w = a kappa_phys(t(eta)) Q^2,
kappa_phys = (g/f) dphi/dt.  Reheating: the axion oscillates in its quadratic potential and dominates,
H = H_end a^(-3/2), a(eta) = (1 + H_end eta/2)^2, t = (2/(3 H_end))(a^(3/2) - 1); its amplitude in units
of f is phi(a) = phi_end a^(-3/2) with phi_end = sqrt3 (H_end/m_phi)(M_Pl/f), Fujita et al. (A10)-(A11)
(the README of 14-15 Oct took phi_end = 1 as an assumption; (A11) fixes it).  phi = f phi(a) cos(m_phi t),
kappa_phys = -m_phi phi(a) sin(m_phi t) - (3/2) H phi(a) cos(m_phi t).  The gauge field's own energy is
not fed back into H (their A6-A7 do that; here a(eta) is the matter-era one).

Averaged prediction (the Lamé-gap resolution, README correction section): the drive renormalises the
quartic coupling to lam_eff(a) = 2 - phi(a)^2/(omega^2 - 1), omega = m_phi/M (Pi = 1/(omega^2(omega^2-1)),
omega^2 phi^2 Pi = phi^2/(omega^2 - 1)); the isotropic solution is lost while lam_eff < 0, i.e. for
a < a_x = (phi_end^2 / (2 (omega^2 - 1)))^(1/3), and only if the runaway completes before a_x.

Parameters per point: m_phi/H_end, r = M/m_phi, m_Q,end = g Q_end/H_end (physical Q at the end of
inflation, 5 as in Fujita et al.), f/M_Pl.  Units g = M = 1: omega = 1/r, H_end = omega/(m_phi/H_end),
Q_end = m_Q H_end.  Initial data at a = 1: Q = Q_end, Q' = H_end Q_end (physical Qdot = 0), w = 0,
w' = omega^2 phi_end Q_end^2/(omega^2 - 1) (the dressed particular solution, as axion_dressed.py B).
Run to a = A_MAX or until |Q| > 10 Q_end.  Report per e-fold of a: max |Q|/Q_end and ln of the comoving
Yang-Mills energy e = Q'^2/2 + Q^4/2 relative to its initial value.

Parts:
  0  the (Q, w) form against the pickle's fourth-order q4 (20 random points) and, at a = 1 frozen, against
     axion_dressed.py Part B (omega = 2, phi = 1 stable; phi = 2.5 blow-up at 0.65 T_Q).
  A  scan over m_phi/H_end and r at f/M_Pl = 0.004, m_Q = 5: verdict against the averaged prediction.
  B  the pure Yang-Mills control (M -> infinity: Q'' = -2Q^3 + a kappa Q^2) at the low-m_phi/H_end points of A.
  B' the slow-drive (Kapitza-threshold) instability in flat space, with and without the ghost, far from and near
     the m_phi = M resonance: the second threshold, which the ghost's dressing of the drive, kappa Q^2/(1 - omega^2),
     sharpens near resonance on either side.
  C  Fujita's benchmark numbers: phi_end from (A11), lam_eff at a = 1, and what m_phi/H_end the threshold needs.

Float only (DOP853, rtol 1e-10).  Usage (from the repository root, ~10 min):
    python axion_su2/code/axion_frw_iso.py > axion_su2/results/axion_frw_iso.txt
Added 2026-10-16.
"""
import os
import pickle
import sys
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.special import ellipk

HERE = os.path.dirname(os.path.abspath(__file__))
PKL = os.path.join(HERE, '..', 'rung1_axion_system.pkl')
RT = dict(method='DOP853', rtol=1e-10, atol=1e-13)
K = ellipk(0.5)
T0 = time.time()
A_MAX = 20.0
MPL_OVER_F = 1.0 / 0.004
MQ = 5.0


def stamp(msg=''):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def Pi(w):
    return 1.0 / (w * w - w)


# ------------------------------------------------------------------ Part 0
def part_0():
    print("\nPart 0: checks\n")
    d = pickle.load(open(PKL, 'rb'))
    args = list(d['X']) + list(d['P']) + list(d['X2']) + list(d['X3']) + [d['cj'], d['ksyms'][0]]
    q4 = sp.lambdify(args, d['q4'][0], 'numpy')
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(20):
        Q, Q1, Q2, Q3, kap = rng.normal(size=5)
        Q4 = q4(Q, Q, Q, Q1, Q1, Q1, Q2, Q2, Q2, Q3, Q3, Q3, 1.0, kap)
        w = Q2 + 2 * Q ** 3
        w2 = kap * Q * Q - (1 + 6 * Q * Q) * w
        Q4b = w2 - (6 * Q * Q * Q2 + 12 * Q * Q1 * Q1)
        worst = max(worst, abs(Q4 - Q4b) / (1 + abs(Q4)))
    print(f"    (Q, w) form vs pickle q4, 20 random points: worst relative residual {worst:.1e}  "
          f"({'PASS' if worst < 1e-12 else 'FAIL'})")
    # frozen a = 1 against axion_dressed.py Part B rows
    for omega, phi, expect in ((2.0, 1.0, 'stable'), (2.0, 2.5, 'BLOW-UP'), (3.0, 3.0, 'stable')):
        Q_amp = omega / 100.0
        w2 = omega ** 2
        rhs = flat_rhs(omega, phi)
        y0 = np.array([Q_amp, 0.0, 0.0, w2 * phi * Q_amp ** 2 / (w2 - 1.0)])
        P = Pi(w2)
        phi_eff = phi * w2 * P * (-1)
        lam = 2.0 - w2 * phi ** 2 * P
        X = Q_amp * (1 + 2 * phi_eff ** 2 / (2 + phi_eff ** 2)) ** 0.25
        T = 4 * K / (np.sqrt(lam) * X) if lam > 0 else 4 * K / (np.sqrt(2) * X)

        def blow(t, y):
            return abs(y[0]) - 10 * Q_amp
        blow.terminal = True
        sol = solve_ivp(rhs, (0, 4 * T), y0, events=blow, **RT)
        got = 'BLOW-UP' if sol.status == 1 else 'stable'
        where = f" at {sol.t_events[0][0] / T:.2f} T" if sol.status == 1 else ''
        print(f"    flat, omega = {omega}, phi = {phi}: {got}{where} (axion_dressed.py B: {expect})  "
              f"({'PASS' if got == expect else 'FAIL'})")


def flat_rhs(omega, phi):
    def rhs(t, y):
        Q, Q1, w, w1 = y
        kap = -omega * phi * np.sin(omega * t)
        return np.array([Q1, w - 2.0 * Q ** 3, w1, kap * Q * Q - (1.0 + 6.0 * Q * Q) * w])
    return rhs


# ------------------------------------------------------------------ FRW run
def frw_run(mphi_over_H, r, mQ=MQ, mpl_over_f=MPL_OVER_F, a_max=A_MAX, ghost=True):
    omega = 1.0 / r
    H = omega / mphi_over_H
    Q_end = mQ * H
    h_end = 1.0 / mphi_over_H
    phi_end = np.sqrt(3.0) * h_end * mpl_over_f
    w2 = omega ** 2

    def a_of(eta):
        return (1.0 + 0.5 * H * eta) ** 2

    def t_of(eta):
        return (2.0 / (3.0 * H)) * (a_of(eta) ** 1.5 - 1.0)

    def kappa_phys(eta):
        a = a_of(eta)
        ph = phi_end * a ** -1.5
        th = omega * t_of(eta)
        Hp = H * a ** -1.5
        return -omega * ph * np.sin(th) - 1.5 * Hp * ph * np.cos(th)

    if ghost:
        def rhs(eta, y):
            Q, Q1, w, w1 = y
            a = a_of(eta)
            a2 = a * a
            return np.array([Q1, a2 * w - 2.0 * Q ** 3, w1, a * kappa_phys(eta) * Q * Q - (a2 + 6.0 * Q * Q) * w])
        y0 = np.array([Q_end, H * Q_end, 0.0, w2 * phi_end * Q_end ** 2 / (w2 - 1.0)])
    else:
        def rhs(eta, y):
            Q, Q1 = y[:2]
            a = a_of(eta)
            return np.array([Q1, -2.0 * Q ** 3 + a * kappa_phys(eta) * Q * Q, 0.0, 0.0])
        y0 = np.array([Q_end, H * Q_end, 0.0, 0.0])

    eta_max = 2.0 * (np.sqrt(a_max) - 1.0) / H

    def blow(eta, y):
        return abs(y[0]) - 10.0 * Q_end
    blow.terminal = True
    sol = solve_ivp(rhs, (0.0, eta_max), y0, dense_output=True, events=blow, **RT)
    eta_end = sol.t[-1]
    # per e-fold of a: max |Q|/Q_end and ln e/e0
    e0 = 0.5 * (H * Q_end) ** 2 + 0.5 * Q_end ** 4
    rows = []
    edges = [2.0 * (np.exp(k) ** 0.5 - 1.0) / H for k in range(int(np.log(a_max)) + 2)]
    for k in range(len(edges) - 1):
        lo, hi = edges[k], min(edges[k + 1], eta_end)
        if hi <= lo:
            rows.append((np.nan, np.nan))
            continue
        tt = np.linspace(lo, hi, 20001)
        Y = sol.sol(tt)
        e = 0.5 * Y[1] ** 2 + 0.5 * Y[0] ** 4
        rows.append((np.max(np.abs(Y[0])) / Q_end, np.log(np.max(e) / e0)))
    blew = sol.status == 1
    a_blow = a_of(sol.t_events[0][0]) if blew else np.nan
    lam0 = 2.0 - phi_end ** 2 / (w2 - 1.0) if w2 > 1 else 2.0 - phi_end ** 2 * w2 * Pi(w2)
    a_x = (phi_end ** 2 / (2.0 * (w2 - 1.0))) ** (1.0 / 3.0) if (w2 > 1 and lam0 < 0) else np.nan
    return dict(omega=omega, H=H, Q_end=Q_end, phi_end=phi_end, lam0=lam0, a_x=a_x, rows=rows, blew=blew,
                a_blow=a_blow, nfev=sol.nfev)


def verdict(res):
    """LOST: |Q| passed 10 Q_end.  Otherwise the orbit persists; report its comoving amplitude in the last e-fold."""
    if res['blew']:
        return 'LOST'
    return 'steady'


def s_slow(mphi_over_H, r, mQ=MQ, mpl_over_f=MPL_OVER_F):
    """The slow-drive (Kapitza) parameter of the isotropic Hill coefficient 2 kappa Q: c1_iso/omega^2 for pure
    Yang-Mills, 2 phi m_Q h_end, and the ghost-dressed one, divided by |1 - omega^2| (the dressed drive on Q is
    kappa Q^2/(1 - omega^2), axion_dressed.py Part C)."""
    omega = 1.0 / r
    h = 1.0 / mphi_over_H
    phi_end = np.sqrt(3.0) * h * mpl_over_f
    s_ym = 2.0 * phi_end * mQ * h
    return s_ym, s_ym / abs(1.0 - omega ** 2)


# ------------------------------------------------------------------ Part A
def part_A():
    print(f"\nPart A: f/M_Pl = {1 / MPL_OVER_F}, m_Q,end = {MQ}, matter-era reheating, a from 1 to {A_MAX:g}; "
          f"phi_end = sqrt3 (H_end/m_phi)(M_Pl/f) from their (A11).")
    print("        Averaged prediction: lost iff lam_eff(1) = 2 - phi_end^2/(omega^2 - 1) < 0, for a < a_x (m_phi > M only).")
    print("        s_YM = 2 phi_end m_Q H_end/m_phi: the slow-drive parameter c1_iso/omega^2 of pure Yang-Mills; s_dr = s_YM/|1 - omega^2|")
    print("        the ghost-dressed one (Part B' locates their threshold).  Columns per e-fold of a (a in [1,e), [e,e^2),")
    print("        [e^2,e^3), [e^3,20]): max|Q|/Q_end, ln e/e0 with e = Q'^2/2 + Q^4/2 comoving (the jitter velocity is in e).\n")
    print(f"  {'m_phi/H':>7} {'r':>5} {'phi_end':>7} {'lam_eff(1)':>10} {'a_x':>5} {'s_YM':>5} {'s_dr':>5} | "
          + "  ".join(f"{'efold ' + str(k):>13}" for k in range(4)) + f" | {'verdict':>7} {'a_blow':>6} {'nfev':>8}")
    lost = []
    for mphi_over_H in (100.0, 200.0, 300.0, 450.0, 700.0):
        for r in (0.5, 0.7, 0.8, 0.9, 1.2):
            res = frw_run(mphi_over_H, r)
            v = verdict(res)
            if v == 'LOST':
                lost.append((mphi_over_H, r))
            s_ym, s_dr = s_slow(mphi_over_H, r)
            cols = "  ".join(f"{q:6.2f} {le:6.2f}" if np.isfinite(q) else f"{'-':>6} {'-':>6}" for q, le in res['rows'][:4])
            pred = 'lost' if res['lam0'] < 0 else 'stable'
            print(f"  {mphi_over_H:7.0f} {r:5.2f} {res['phi_end']:7.3f} {res['lam0']:10.3f} "
                  f"{res['a_x']:5.2f} {s_ym:5.2f} {s_dr:5.2f} | {cols} | {v:>7} {res['a_blow']:6.2f} {res['nfev']:8d}   (averaged: {pred})"
                  f"   [{time.time() - T0:5.0f}s]", flush=True)
    return lost


# ------------------------------------------------------------------ Part B
def part_B(lost):
    print("\nPart B: pure Yang-Mills control (no ghost, Q'' = -2Q^3 + a kappa Q^2) at every point with m_phi/H_end <= 200")
    print("        and at the r = 1.2 points of 300 and 450: the slow-drive instability exists without the ghost.\n")
    pts = [(m, r) for m in (100.0, 200.0) for r in (0.5, 0.7, 0.9, 1.2)] + [(300.0, 0.7), (300.0, 1.2), (450.0, 1.2)]
    for mphi_over_H, r in pts:
        res = frw_run(mphi_over_H, r, ghost=False)
        s_ym, _ = s_slow(mphi_over_H, r)
        cols = "  ".join(f"{q:6.2f} {le:6.2f}" if np.isfinite(q) else f"{'-':>6} {'-':>6}" for q, le in res['rows'][:4])
        print(f"  {mphi_over_H:7.0f} {r:5.2f} phi_end {res['phi_end']:6.3f} s_YM {s_ym:5.2f} | {cols} | {verdict(res):>7} "
              f"{res['a_blow']:6.2f}   [{time.time() - T0:5.0f}s]", flush=True)


# ------------------------------------------------------------------ Part B'
def part_Bp():
    print("\nPart B': the slow-drive threshold in flat space, 16 Q periods, m_phi/gQ = 40 (Q_amp = omega/40):")
    print("         phi chosen so that s = 2 phi Q_amp/(omega |1 - omega^2|) (with ghost) or s_YM = 2 phi Q_amp/omega (no ghost)")
    print("         takes the listed values; LOST = |Q| > 10 Q_amp.  omega = 0.5 (far from M) and 0.833 (r = 1.2, near it).\n")
    print(f"  {'omega':>5} {'s':>5} {'phi (ghost)':>11} {'ghost':>8} {'phi (YM)':>9} {'no ghost':>9}")
    for omega in (0.5, 0.833):
        Q_amp = omega / 40.0
        w2 = omega ** 2
        for s in (0.2, 0.3, 0.4, 0.5, 0.7):
            out = []
            for ghost in (True, False):
                phi = s * omega * (abs(1 - w2) if ghost else 1.0) / (2 * Q_amp)
                rhs = flat_rhs(omega, phi) if ghost else (lambda t, y, om=omega, ph=phi: np.array(
                    [y[1], -2 * y[0] ** 3 - om * ph * np.sin(om * t) * y[0] ** 2, 0.0, 0.0]))
                y0 = np.array([Q_amp, 0.0, 0.0, w2 * phi * Q_amp ** 2 / (w2 - 1.0) if ghost else 0.0])
                Tq = 4 * K / (np.sqrt(2) * Q_amp)

                def blow(t, y):
                    return abs(y[0]) - 10 * Q_amp
                blow.terminal = True
                sol = solve_ivp(rhs, (0, 16 * Tq), y0, events=blow, **RT)
                out.append((phi, f"LOST {sol.t_events[0][0] / Tq:.1f}T" if sol.status == 1 else "steady"))
            print(f"  {omega:5.3f} {s:5.2f} {out[0][0]:11.3f} {out[0][1]:>8} {out[1][0]:9.3f} {out[1][1]:>9}   [{time.time() - T0:5.0f}s]",
                  flush=True)


# ------------------------------------------------------------------ Part C
def part_C():
    MPl, g, mQ, mphi, f = 2.4e18, 1e-12, 5.0, 5e10, 0.004
    H = g * MPl / mQ ** 2
    h = H / mphi
    phi_end = np.sqrt(3.0) * h / f
    print("\nPart C: Fujita et al.'s benchmark (axion_aniso.py GEV: g = 1e-12, m_Q,end = 5, m_phi = 5e10 GeV, f = 0.004 M_Pl)\n")
    print(f"    H_end = {H:.2e} GeV (their (33)), h_end = H_end/m_phi = {h:.2e}, phi_end = sqrt3 h_end M_Pl/f = {phi_end:.2e} (their (A11))")
    for r in (0.9, 0.99, 0.999):
        lam = 2.0 - phi_end ** 2 * r ** 2 / (1 - r ** 2)
        print(f"    r = {r}: lam_eff(1) = 2 - phi_end^2 r^2/(1 - r^2) = {lam:.6f}  (stable; the threshold needs 1 - r^2 < "
              f"{phi_end ** 2 / 2:.1e})")
    print("    The threshold phi_end^2 r^2/(1 - r^2) > 2 with phi_end = sqrt3 (H_end/m_phi)(M_Pl/f) reads")
    print("        m_phi/H_end < sqrt(3/2) (M_Pl/f) r/sqrt(1 - r^2):")
    for fM in (0.002, 0.004, 0.04):
        s = ", ".join(f"r = {r}: {np.sqrt(1.5) / fM * r / np.sqrt(1 - r ** 2):6.0f}" for r in (0.5, 0.7, 0.9))
        print(f"        f/M_Pl = {fM}: {s}")
    print("    against the fast-drive condition m_phi >> g Q_end = m_Q H_end = 5 H_end.")
    print("    The slow-drive threshold of Part B', s_YM = 2 phi_end m_Q H_end/m_phi = 2 sqrt3 m_Q (M_Pl/f)(H_end/m_phi)^2 between")
    print("    0.2 and 0.3 (pure Yang-Mills, any M; with the ghost the same far from resonance, lower near it) reads, at 0.25,")
    print("    m_phi/H_end < sqrt(2 sqrt3 m_Q (M_Pl/f)/0.25):")
    for fM in (0.002, 0.004, 0.04):
        print(f"        f/M_Pl = {fM}: m_phi/H_end < {np.sqrt(2 * np.sqrt(3) * 5.0 / fM / 0.25):5.0f}")
    print("    Both thresholds are quadratic in H_end/m_phi; both are out of reach at the benchmark by six and nine orders of")
    print("    magnitude in that ratio squared.  The gauge field's back-reaction on the axion (their A6; EL of (g/f) phidot Q^3")
    print("    for phi: phiddot + m_phi^2 phi = -(g/f) d(Q^3)/dt), relative to m_phi^2 phi with g^2 Q^4 ~ rho_Q ~ M_Pl^2 H^2 (their")
    print("    epsilon_B = 1), is of order sqrt6 (M_Pl/f)(H_end/m_phi) = sqrt2 phi_end: negligible at the benchmark (1e-3) and O(1)")
    print("    wherever either threshold is reachable.  So in the window the prescribed-axion runs locate, the coupled system is")
    print("    required; outside it the ghost does nothing to the isotropic solution.")


if __name__ == "__main__":
    print(f"axion_frw_iso.py  python {sys.version.split()[0]}  numpy {np.__version__}  ({sys.executable})")
    part_0()
    stamp("Part 0 done")
    lost = part_A()
    stamp("Part A done")
    part_B(lost)
    stamp("Part B done")
    part_Bp()
    stamp("Part B' done")
    part_C()
    stamp("Part C done")
    print(f"\n({time.time() - T0:.0f} s)")
