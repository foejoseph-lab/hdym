"""dilaton_iso.py -- item 5b, the runs (17 Oct; column A primary, column B control -- Joseph, 17 Oct).  The
isotropic SU(2) background under a prescribed dilaton Theta(t) = exp(delta_a cos(omega t)), delta = phi/Lambda, with the
Lee-Wick ghost, in the (q, w) auxiliary form of dilaton_derive.py; flat space, g = M = 1, the regime of axion_su2/:
Q_amp = omega/100, so gQ_amp = M/100 at omega = 1 (omega = m_phi/M).  Everything is an initial-value problem in phase
space on the equations of dilaton_derive.py (checked there, all PASS) and, for column A, on
axion_su2/code/probe_core.py, which carries Theta.  "Periods" and "T" below and in the log are those of the undriven
background orbit, T = 4K(1/2)/(sqrt2 Q_amp) ~ 524/omega (262 time units at omega = 2), not the drive period 2 pi/omega.

Result (17 Oct, README): the ghost's effective mass M^2 Theta(t) is modulated and the isotropic solution is lost --
within 0.08-2.8 T, every realisation, no threshold in delta_a -- wherever omega/M sits in the ghost's parametric
tongue (omega = 2 for every delta_a; 2.4 at delta_a = 1 in A; B wider, lost also at 1.5 and 3.0 at delta_a = 1);
off the tongue steady for 60 periods, seeds identical, E_heavy bounded; M -> infinity steady everywhere.  Probes off
the tongue: spectators at m_p = 0 and 3 gQ_amp; at m_p = gQ_amp a bounded x10 ghost-mediated exchange, c <= 7e-4
against 1/4 (6e-4 until 28 Oct 2026; the log's 'both ghost 1' row at (3.0, 1.0) reads 0.0007).  Decision of 17 Oct: column A primary (the cited form, Gauss law untouched), column B the control.

  column A  (Theta on F^2 only, the cited coupling):  Q'' = w - 2Q^3,  w'' + (Theta + 6Q^2) w = -Theta-dot Q'
            (probe_core with Theta, kappa = 0; the lean isotropic rhs here is checked against it at random points).
  column B  (Theta on the whole operator):  the Euler-Lagrange equations of L'_B = Theta L' + Theta-dot w.qdot on
            x = y = z, derived with sympy at import, lambdified; checked to coincide with A at Theta = 1 and to conserve
            the energy at constant Theta.
  M -> infinity (both columns):  Theta Q'' + Theta-dot Q' + 2 Theta Q^3 = 0  (Yang-Mills with g -> g Theta in
            d tau = dt/Theta; <Theta^2>_tau = 1 for a cosine, so no inversion of the potential: the control).
  The (p, w) twin is not definable with Theta on (probe_core docstring); column B is the control here (17 Oct).

  A  The isotropic background without probes: omega in (0.5, 0.8, 1.5, 2.0, 2.4, 3.0), delta_a in (0.1, 0.3, 1.0), 8
     periods of the undriven orbit T = 4K/(sqrt2 Q_amp): |Q|max/Q_amp per cycle, the heavy sector's cycle-mean energy
     over E_bg = (3/2) Q_amp^4, lost at |Q| = 10 Q_amp; columns A, B, M -> infinity.  Then 60 periods at omega = 2.4 and
     3.0, delta_a = 0.3 and 1.0, three perturbed seeds each (the fraction statistics of the probe programme).  This
     answers the derivation's open question: is the friction-like drive -Theta-dot Q' a threshold (as the axion's
     renormalised coupling) or a secular pump?  (Neither: the w equation's homogeneous part is the anisotropic Mathieu
     equation, sourced by the drive and fed back through Q'' = w - 2Q^3 -- the ghost's first-order parametric
     resonance at m_phi = 2M/n, and nothing else.)
  B  The anisotropic sector (background-independent, dilaton_derive.py CHECK 6): the one-period monodromy of
        A: u'' + Theta(t) u = 0,       B: (Theta u')' + Theta u = 0,       (u = eps'; the constant C is a drift)
     over omega in a fine grid 0.3-3.2 and delta_a in (0.1, 0.3, 1.0): the growth rate per unit time, the tongues at
     omega = 2/n, against the Mathieu heuristics printed by the derivation (A: delta_a/4 at omega = 2, B: delta_a/2).
     This sector is the same for a healthy heavy mode (a Hill equation with positive coefficient), so it is not a ghost
     statement; recorded as the dilaton's anisotropic map.
  C  The probes on the dilaton background (column A, probe_core): boson / fermion / both at m_p/gQ_amp = 0, 1, 3 with
     the ghost and with M -> infinity, at (omega, delta_a) = (2.4, 0.3) and (3.0, 1.0), 8 periods; per cycle |Q|max,
     E_boson/E_boson(0), fermion n_-, max c = (g^2/2) f^2 against the collision value 1/4.

Controls that must be steady: every M -> infinity row.

Usage (from the repository root, ~13 min):
    python dilaton_su2/code/dilaton_iso.py > dilaton_su2/results/dilaton_iso.txt
Added 2026-10-17.
"""
import os
import sys
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.special import ellipk

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'axion_su2', 'code'))
import probe_core as pc  # noqa: E402

T0 = time.time()
RT = dict(method='DOP853', rtol=1e-9, atol=1e-16)
K = ellipk(0.5)
RATIO = 100.0
EPS_P = 0.1
NCYC = 8


def stamp(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


# ------------------------------------------------------------------ the drive
def theta_funs(omega, da):
    Th = lambda t: np.exp(da * np.cos(omega * t))
    Thd = lambda t: -da * omega * np.sin(omega * t) * np.exp(da * np.cos(omega * t))
    Thdd = lambda t: (da * da * omega * omega * np.sin(omega * t) ** 2 - da * omega * omega * np.cos(omega * t)) * np.exp(da * np.cos(omega * t))
    return Th, Thd, Thdd


# ------------------------------------------------------------------ column A: isotropic, from the core
def rhs_A(P):
    """Isotropic truncation of probe_core.rhs with Theta on (column A), lean 10-d (Q, Q', w, w', f, f', Re zeta, Im zeta)."""
    g, M, ms, mD = P.g, P.M, P.ms, P.mD
    c = pc.C_FERMION
    bos, fer = P.boson, P.fermion

    def rhs(t, y):
        Q, Q1, w, w1, f, f1 = y[:6]
        zr, zi = y[6:8], y[8:10]
        Th, Thd = P.Theta(t), P.Thetadot(t)
        F = 0.0
        if bos:
            F -= 0.5 * g * g * Q * f * f
        if fer:
            F += g * (-(c / 3.0) * 2.0 * (zr[0] * zr[1] + zi[0] * zi[1]))
        gV = 2.0 * g * g * Q ** 3
        Qdd = M * M * w - gV
        wdd = -6.0 * g * g * Q * Q * w - Th * M * M * w - Thd * Q1 + F
        out = np.zeros(10)
        out[0], out[1], out[2], out[3] = Q1, Qdd, w1, wdd
        if bos:
            out[4], out[5] = f1, -(ms * ms + 0.75 * g * g * Q * Q) * f
        if fer:
            h01 = c * g * Q
            hz_r = np.array([mD * zr[0] + h01 * zr[1], h01 * zr[0] - mD * zr[1]])
            hz_i = np.array([mD * zi[0] + h01 * zi[1], h01 * zi[0] - mD * zi[1]])
            out[6:8], out[8:10] = hz_i, -hz_r
        return out
    return rhs


def iso_pack(s):
    return np.array([s[0], s[3], s[6], s[9], s[12], s[13], s[14], s[15], s[16], s[17]])


def iso_unpack(y):
    s = np.zeros(18)
    s[0:3], s[3:6], s[6:9], s[9:12] = y[0], y[1], y[2], y[3]
    s[12], s[13] = y[4], y[5]
    s[14:18] = y[6:10]
    return s


def check_rhs_A():
    rng = np.random.default_rng(7)
    Th, Thd, _ = theta_funs(2.4, 0.7)
    err = 0.0
    for bos in (False, True):
        for fer in (False, True):
            P = pc.Params(boson=bos, fermion=fer, ms=0.7, mD=0.6, Theta=Th, Thetadot=Thd)
            for _ in range(20):
                y = rng.normal(size=10)
                if not bos:
                    y[4:6] = 0
                if not fer:
                    y[6:10] = 0
                tt = rng.uniform(0, 5)
                a = rhs_A(P)(tt, y)
                b = iso_pack(pc.rhs(tt, iso_unpack(y), P))
                err = max(err, np.max(np.abs(a - b)) / max(1.0, np.max(np.abs(b))))
    print(f"  CHECK column A isotropic rhs against probe_core.rhs with Theta on, all probe switches, random points: "
          f"max rel err {err:.1e} {'PASS' if err < 1e-12 else 'FAIL'}", flush=True)


# ------------------------------------------------------------------ column B: isotropic, from L'_B by sympy
def build_B():
    t = sp.Symbol('t', real=True)
    Q, w = sp.Function('Q')(t), sp.Function('w')(t)
    Th = sp.Function('Theta', positive=True)(t)
    V = sp.Rational(3, 2) * Q ** 4                      # V on x = y = z, g = 1
    gradV = 2 * Q ** 3                                  # component of grad V along each axis
    # L' per component (three equal components): 1/2 Q'^2 + Q' w' - V/3... careful: L' = sum_i [1/2 q_i'^2 + q_i' w_i'] - V - w.gradV + M^2/2 |w|^2
    Lp = 3 * (sp.Rational(1, 2) * Q.diff(t) ** 2 + Q.diff(t) * w.diff(t)) - V - 3 * w * gradV + sp.Rational(3, 2) * w ** 2
    LB = Th * Lp + Th.diff(t) * 3 * w * Q.diff(t)
    EL = [sp.expand(e.lhs) for e in sp.euler_equations(LB, [Q, w], t)]
    sol = sp.solve(EL, [Q.diff(t, 2), w.diff(t, 2)], dict=True)[0]
    syms = sp.symbols('Qs Q1s ws w1s Ths Thds Thdds')
    rep = {Q.diff(t, 2): 0}
    sub = {Th.diff(t, 2): syms[6], Th.diff(t): syms[5], Th: syms[4], Q.diff(t): syms[1], Q: syms[0], w.diff(t): syms[3], w: syms[2]}

    def lam(e):
        for k_ in (Th.diff(t, 2), Th.diff(t), Th, Q.diff(t), Q, w.diff(t), w):
            e = e.subs(k_, sub[k_])
        return sp.lambdify(syms, e, 'numpy')
    fQ, fw = lam(sol[Q.diff(t, 2)]), lam(sol[w.diff(t, 2)])
    # energy: Legendre transform of LB
    vel = [Q.diff(t), w.diff(t)]
    E = sp.expand(sum(v * LB.diff(v) for v in vel) - LB)
    fE = lam(E)
    return fQ, fw, fE, sol, Q, w, Th, t


FQ_B, FW_B, FE_B, SOL_B, _Q, _w, _Th, _t = build_B()


def rhs_B(omega, da):
    Th, Thd, Thdd = theta_funs(omega, da)

    def rhs(t, y):
        Q, Q1, w, w1 = y[:4]
        th, thd, thdd = Th(t), Thd(t), Thdd(t)
        out = np.zeros(10)
        out[0], out[1], out[2], out[3] = Q1, FQ_B(Q, Q1, w, w1, th, thd, thdd), w1, FW_B(Q, Q1, w, w1, th, thd, thdd)
        return out
    return rhs


def check_B():
    rng = np.random.default_rng(8)
    # at Theta = 1 column B equals column A (no probes)
    P = pc.Params()
    err = 0.0
    for _ in range(20):
        y = np.zeros(10); y[:4] = rng.normal(size=4)
        a = rhs_A(P)(0.0, y)[:4]
        b = np.array([y[1], FQ_B(*y[:4], 1.0, 0.0, 0.0), y[3], FW_B(*y[:4], 1.0, 0.0, 0.0)])
        err = max(err, np.max(np.abs(a - b)))
    print(f"  CHECK column B at Theta = 1 equals column A (random points): max abs err {err:.1e} {'PASS' if err < 1e-12 else 'FAIL'}", flush=True)
    # energy conserved at constant Theta = 1.7
    th = 1.7
    f = lambda t, y: np.concatenate([[y[1], FQ_B(*y[:4], th, 0.0, 0.0), y[3], FW_B(*y[:4], th, 0.0, 0.0)], np.zeros(6)])
    y0 = np.zeros(10); y0[:4] = [0.3, 0.1, 0.02, -0.01]
    sol = solve_ivp(f, (0, 30), y0, t_eval=np.linspace(0, 30, 31), method='DOP853', rtol=1e-12, atol=1e-14)
    E = np.array([FE_B(*yy[:4], th, 0.0, 0.0) for yy in sol.y.T])
    d = np.max(np.abs(E - E[0])) / max(1.0, abs(E[0]))
    print(f"  CHECK column B energy conserved at constant Theta over t = 30: drift {d:.1e} {'PASS' if d < 1e-9 else 'FAIL'}", flush=True)
    print(f"  column B isotropic equations (sympy): Q'' = {sp.simplify(SOL_B[_Q.diff(_t, 2)])}")
    print(f"                                        w'' = {sp.simplify(SOL_B[_w.diff(_t, 2)])}")


# ------------------------------------------------------------------ the isotropic run
def run_iso(omega, da, column, ncyc=NCYC, probe='none', mp_over_gQ=0, seed=0):
    Q_amp = omega / RATIO
    T = 4 * K / (np.sqrt(2) * Q_amp)
    Th, Thd, Thdd = theta_funs(omega, da)
    m_p = mp_over_gQ * Q_amp
    P = pc.Params(boson=probe in ('boson', 'both'), fermion=probe in ('fermion', 'both'), ms=m_p, mD=m_p, Theta=Th, Thetadot=Thd)
    y = np.zeros(10)
    y[0] = Q_amp
    E_bg = 1.5 * Q_amp ** 4
    if P.boson:
        y[4] = np.sqrt(EPS_P * E_bg / (m_p ** 2 + 0.75 * Q_amp ** 2))
    if P.fermion:
        h = pc.fermion_h(Q_amp, P)
        ev, U = np.linalg.eigh(h)
        zeta = np.sqrt(EPS_P * E_bg / ev[1]) * U[:, 1]
        y[6:8], y[8:10] = zeta.real, zeta.imag
    if seed:
        rng = np.random.default_rng(seed)
        y *= 1 + 1e-6 * rng.normal(size=10)
    if column == 'A':
        f = rhs_A(P)
    elif column == 'B':
        f = rhs_B(omega, da)
    else:   # M -> infinity: Theta Q'' + Theta' Q' + 2 Theta Q^3 = F
        base = rhs_A(P)

        def f(t, yy):
            z = yy.copy(); z[2] = z[3] = 0.0
            out = base(t, z)
            F = out[3] + Thd(t) * yy[1]          # out[3] at w = 0 is -Thd Q' + F
            out[1] = (-Thd(t) * yy[1] - 2.0 * Th(t) * yy[0] ** 3 + F) / Th(t)
            out[2] = out[3] = 0.0
            return out
        y[2] = y[3] = 0.0

    def blow(t, yy):
        return abs(yy[0]) - 10 * Q_amp
    blow.terminal = True
    E_b0 = (m_p ** 2 + 0.75 * Q_amp ** 2) * y[4] ** 2 if P.boson else 1.0
    per, eh, eb, nf, cm = [], [], [], [], []
    for k in range(ncyc):
        sol = solve_ivp(f, (k * T, (k + 1) * T), y, dense_output=True, events=blow, **RT)
        tt = np.linspace(k * T, sol.t[-1], 1001)
        Y = sol.sol(tt)
        per.append(np.max(np.abs(Y[0])) / Q_amp)
        if column == 'A':
            eh.append(np.mean([pc.energy_parts(t_, iso_unpack(Y[:, i]), P)[1] for i, t_ in enumerate(tt[::50])]) / E_bg)
        elif column == 'B':
            # heavy share of column B's energy: the terms in w', w of the Legendre transform (Theta-weighted, bookkeeping)
            eh.append(np.mean([FE_B(*Y[:4, i], Th(t_), Thd(t_), Thdd(t_)) - 1.5 * Th(t_) * Y[1, i] ** 2 - 1.5 * Th(t_) * Y[0, i] ** 4
                               for i, t_ in enumerate(tt[::50])]) / E_bg)
        else:
            eh.append(0.0)
        eb.append(np.mean((Y[5] ** 2 + (m_p ** 2 + 0.75 * Y[0] ** 2) * Y[4] ** 2)[::50]) / E_b0 if P.boson else np.nan)
        nf.append(np.mean([pc.fermion_adiabatic_number(Y[0, i], Y[6:8, i] + 1j * Y[8:10, i], P) for i in range(0, 1001, 50)]) if P.fermion else np.nan)
        cm.append(0.5 * np.max(Y[4] ** 2))
        y = sol.y[:, -1]
        if sol.status == 1:
            lost = sol.t_events[0][0] / T
            n = ncyc - k - 1
            per += [np.nan] * n; eh += [np.nan] * n; eb += [np.nan] * n; nf += [np.nan] * n; cm += [np.nan] * n
            return map(np.array, (per, eh, eb, nf, cm)), lost, T
    return map(np.array, (per, eh, eb, nf, cm)), None, T


def f2(v, fmt="{:5.2f}"):
    return " ".join(fmt.format(x) if np.isfinite(x) else " " * (len(fmt.format(0.0)) - 1) + "-" for x in v)


# ------------------------------------------------------------------ Part B: anisotropic monodromy
def aniso_rate(omega, da, column):
    Th, Thd, _ = theta_funs(omega, da)
    Td = 2 * np.pi / omega

    def rhs(t, y):
        E = y.reshape(2, 2)
        if column == 'A':        # u'' + Theta u = 0
            A = np.array([[0.0, 1.0], [-Th(t), 0.0]])
        else:                    # (Theta u')' + Theta u = 0  ->  u'' = -u - (Theta'/Theta) u'
            A = np.array([[0.0, 1.0], [-1.0, -Thd(t) / Th(t)]])
        return (A @ E).ravel()
    sol = solve_ivp(rhs, (0, Td), np.eye(2).ravel(), method='DOP853', rtol=1e-11, atol=1e-13)
    Mo = sol.y[:, -1].reshape(2, 2)
    ev = np.linalg.eigvals(Mo)
    return np.max(np.log(np.abs(ev))) / Td


if __name__ == "__main__":
    check_rhs_A()
    check_B()
    stamp("A: the isotropic background without probes, 8 periods; columns A / B / M -> infinity")
    print(f"  {'omega':>5} {'d_a':>4} {'col':>3} | |Q|max/Q_amp per cycle (8)" + " " * 22 + "| E_heavy/E_bg cycle mean" + " " * 26 + "| lost at")
    for omega in (0.5, 0.8, 1.5, 2.0, 2.4, 3.0):
        for da in (0.1, 0.3, 1.0):
            for column in ('inf', 'A', 'B'):
                (per, eh, eb, nf, cm), lost, T = run_iso(omega, da, column)
                tag = f"{lost:.2f} T" if lost is not None else "-"
                if column == 'inf' and lost is not None:
                    tag += "   CONTROL-FAIL"
                print(f"  {omega:5.2f} {da:4.1f} {column:>3} | {f2(per)} | {f2(eh, '{:6.3f}')} | {tag}", flush=True)
    stamp("A, long: 60 periods at omega = 2.4, 3.0 and delta_a = 0.3, 1.0; columns A, B, M -> infinity; three seeds for A")
    print(f"  {'omega':>5} {'d_a':>4} {'col':>3} {'seed':>4} | |Q|max/Q_amp min, max per 20-cycle block | E_heavy/E_bg max |.| per block | lost at")
    for omega in (2.4, 3.0):
        for da in (0.3, 1.0):
            for column, seeds in (('inf', (0,)), ('A', (0, 1, 2)), ('B', (0,))):
                for seed in seeds:
                    (per, eh, eb, nf, cm), lost, T = run_iso(omega, da, column, ncyc=60, seed=seed)
                    blocks = ((0, 20), (20, 40), (40, 60))
                    s1 = " / ".join(f"{np.nanmin(per[a:b]):5.2f} {np.nanmax(per[a:b]):5.2f}" if np.any(np.isfinite(per[a:b])) else "    -     -" for a, b in blocks)
                    s2 = " ".join(f"{np.nanmax(np.abs(eh[a:b])):7.3f}" if np.any(np.isfinite(eh[a:b])) else "      -" for a, b in blocks)
                    tag = f"{lost:.2f} T" if lost is not None else "-"
                    if column == 'inf' and lost is not None:
                        tag += "   CONTROL-FAIL"
                    print(f"  {omega:5.2f} {da:4.1f} {column:>3} {seed:4d} | {s1} | {s2} | {tag}", flush=True)

    stamp("B: the anisotropic sector, one-period monodromy, growth rate per unit time; columns A and B")
    print("  Mathieu heuristics (dilaton_derive.py): first tongue omega = 2M, rate delta_a/4 (A), delta_a/2 (B); tongues at 2/n")
    grid = np.concatenate([np.arange(0.3, 3.21, 0.05), [2.0 / n for n in (1, 2, 3, 4)], [1.0]])
    grid = np.unique(np.round(grid, 4))
    for da in (0.1, 0.3, 1.0):
        rA = np.array([aniso_rate(om, da, 'A') for om in grid])
        rB = np.array([aniso_rate(om, da, 'B') for om in grid])
        print(f"\n  delta_a = {da}: unstable omega (rate > 1e-4) for A: " + ", ".join(f"{om:.2f}({r:.3f})" for om, r in zip(grid, rA) if r > 1e-4))
        print(f"  delta_a = {da}: unstable omega (rate > 1e-4) for B: " + ", ".join(f"{om:.2f}({r:.3f})" for om, r in zip(grid, rB) if r > 1e-4))
        iA, iB = int(np.argmax(rA)), int(np.argmax(rB))
        print(f"  max rate: A {rA[iA]:.4f} at omega = {grid[iA]:.2f} (heuristic {da / 4:.4f} at 2.00); B {rB[iB]:.4f} at omega = {grid[iB]:.2f} (heuristic {da / 2:.4f})", flush=True)

    stamp("C: the probes on the dilaton background (column A), ghost and M -> infinity, 8 periods")
    print(f"  {'omega':>5} {'d_a':>4} {'probe':>7} {'heavy':>5} {'m_p/gQ':>6} | |Q|max/Q_amp per cycle" + " " * 26 + "| E_boson/E_boson(0)" + " " * 30 + "| fermion n_-" + " " * 38 + "| max c (1/4 collides) | lost at")
    for omega, da in ((2.4, 0.3), (3.0, 1.0)):
        for mp in (0, 1, 3):
            for probe in ('boson', 'fermion', 'both'):
                for column in ('A', 'inf'):
                    (per, eh, eb, nf, cm), lost, T = run_iso(omega, da, column, probe=probe, mp_over_gQ=mp)
                    tag = f"{lost:.2f} T" if lost is not None else "-"
                    if column == 'inf' and lost is not None:
                        tag += "   CONTROL-FAIL"
                    hv = 'ghost' if column == 'A' else 'inf'
                    print(f"  {omega:5.2f} {da:4.1f} {probe:>7} {hv:>5} {mp:6d} | {f2(per)} | {f2(eb)} | {f2(nf, '{:5.3f}')} | {np.nanmax(cm):8.4f} | {tag}", flush=True)
    stamp("C, long: the boson at m_p = gQ_amp over 60 periods, ghost and M -> infinity (is the 8-cycle pump bounded?)")
    print(f"  {'omega':>5} {'d_a':>4} {'heavy':>5} | |Q|max/Q_amp min, max per 20-cycle block | E_boson/E_boson(0) min, max per block | max c per block | lost at")
    for omega, da in ((2.4, 0.3), (3.0, 1.0)):
        for column in ('A', 'inf'):
            (per, eh, eb, nf, cm), lost, T = run_iso(omega, da, column, ncyc=60, probe='boson', mp_over_gQ=1)
            blocks = ((0, 20), (20, 40), (40, 60))
            s1 = " / ".join(f"{np.nanmin(per[a:b]):5.2f} {np.nanmax(per[a:b]):5.2f}" if np.any(np.isfinite(per[a:b])) else "    -     -" for a, b in blocks)
            s2 = " / ".join(f"{np.nanmin(eb[a:b]):6.2f} {np.nanmax(eb[a:b]):7.2f}" if np.any(np.isfinite(eb[a:b])) else "     -       -" for a, b in blocks)
            s3 = " ".join(f"{np.nanmax(cm[a:b]):7.4f}" if np.any(np.isfinite(cm[a:b])) else "      -" for a, b in blocks)
            tag = f"{lost:.2f} T" if lost is not None else "-"
            if column == 'inf' and lost is not None:
                tag += "   CONTROL-FAIL"
            print(f"  {omega:5.2f} {da:4.1f} {'ghost' if column == 'A' else 'inf':>5} | {s1} | {s2} | {s3} | {tag}", flush=True)
    print(f"DONE {time.time() - T0:.0f}", flush=True)
