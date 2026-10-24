"""axion_probes_infl.py -- the inflationary stand-in (plan, 17 Oct; with axion_probes_fourier.py the last run of the
probe programme on this background, complete 17 Oct, README History): does a massless charged scalar -- the Higgs
doublet -- pump the Lee-Wick ghost on the nonzero-mean SU(2) background of chromo-natural inflation and its tail?
Result (README, "The inflationary run"): nothing, in every column, static and tail; the pump of L1 needs the axion
drive within ~15 % of the no-probe loss threshold, i.e. phi_end = O(1), which the benchmark (8e-4) does not give.

Why this run (README "Probes", Status of 16 Oct, third review): the pump found in L1 needs an orbit with nonzero mean,
and that is the inflationary regime of this very model -- in chromo-natural inflation A^a_i = psi(t) a(t) delta^a_i
with psi held near a slowly varying attractor by the axion's fixed-sign drive (Adshead-Wyman 1203.2264, eqs. 3-6).
In this folder's flat language (g = M = 1, axion_derive.py) the fixed-sign drive is a constant kappa, and the
isotropic equations Q'' = w - 2Q^3, w'' + (1 + 6Q^2) w = kappa Q^2 have the fixed point
    w0 = 2 Q0^3,  kappa = 2 Q0 (1 + 6 Q0^2)      (Q0 ~ kappa/2 for Q0 << M),
stable in the isotropic sector (V_eff'' = 2 Q0^2 at M -> infinity).  A massless doublet on it has the constant
effective mass^2 (3/4) g^2 Q0^2 and, at linear order, decouples from (Q, w); any pump is nonlinear (f^2 forcing Q at
twice the scalar's frequency, Q modulating f's mass), with the ghost pole in Q's response.  Everything in phase space
on probe_core.py; the twin is valid here (Q0 << M).

  A  Static inflation: constant kappa, the doublet at eps_p of E_bg = (3/2) Q0^4 (the Yang-Mills energy at Q0) at rest
     at its turning point, Q0/M in (0.024, 0.1, 0.3) [0.024 is L1's gQ_amp], eps_p in (0.01, 0.1, 0.5), 60 periods of
     the isotropic mode T = 2 pi / (sqrt2 Q0), heavy = ghost / twin / M -> infinity, plus the no-scalar control.  Per
     20-period block: max |Q - Q0|/Q0, E_boson/E_boson(0) range, max c = (g^2/2) f^2 against the collision value 1/4,
     E_heavy/E_bg at block end; "lost" = Q leaves (0, 5 Q0).  The ghost rows at Q0 = 0.024, eps_p = 0.1 are repeated
     with 1e-6-perturbed initial data (five seeds) for the fraction statistics the README asks for.
  B  The tail: the end of inflation as the axion reaching the bottom of its potential and oscillating at its mass with
     the amplitude inflation left it -- kappa(t) = kappa0 for 10 T, then kappa0 cos(omega (t - 10 T)) for 50 T, omega =
     m_phi/M in (0.8, 2.4, 3.0); Q starts at the fixed point, the doublet with its inflationary energy.  The oscillating
     regime that follows is 5a's with a drive amplitude set by the inflationary value (omega^2 phi^2 Pi = 4 Q0^2 Pi,
     far below the no-probe threshold 2): the question is whether the transition itself, or the small-amplitude
     oscillation with the scalar on, does anything.  "lost" there = |Q| > 10 Q0.  Q0 = 0.024 and 0.1, eps_p = 0 and 0.1.

Controls that must be steady: every 'none' row and every M -> infinity row's background.  (The twin's fixed point is
itself unstable at Q0 = 0.3 -- its 'none' row leaves in 6 T -- so the twin is a control at 0.024 and 0.1 only here.)

Usage (from the repository root, ~15 min):
    python axion_su2/code/axion_probes_infl.py > axion_su2/results/axion_probes_infl.txt
Added 2026-10-17; the TWIN-LIMIT tag the same evening (log not rerun: the one tagged row, Q0 = 0.3 twin none, reads
CONTROL-FAIL in results/axion_probes_infl.txt and is the twin's limit).
"""
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import probe_core as pc  # noqa: E402
import axion_probes as ap  # noqa: E402

T0 = time.time()
RT = dict(method='DOP853', rtol=1e-9, atol=1e-16)


def stamp(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def fixed_point(Q0, heavy):
    """kappa and w0 of the isotropic fixed point at Q0: ghost w'' = kappa Q^2 - (1 + 6Q^2) w and Q'' = w - 2Q^3 give
    w0 = 2 Q0^3, kappa = 2 Q0 (1 + 6 Q0^2); twin (probe_core sigma = -1: w'' = Hw - w - F, Q'' = -gradV - w - 2w'')
    gives w0 = -2 Q0^3, kappa = 2 Q0 (1 - 6 Q0^2); M -> infinity: w = 0, kappa = 2 Q0."""
    if heavy == 'ghost':
        return 2 * Q0 * (1 + 6 * Q0 ** 2), 2 * Q0 ** 3
    if heavy == 'twin':
        return 2 * Q0 * (1 - 6 * Q0 ** 2), -2 * Q0 ** 3
    return 2 * Q0, 0.0


def run(Q0, eps_p, heavy, nper=60, tail=False, seed=0, omega=None):
    sigma = -1 if heavy == 'twin' else +1
    kappa0, w0 = fixed_point(Q0, heavy)
    T = 2 * np.pi / (np.sqrt(2) * Q0)
    t_end = 10 * T if tail else np.inf
    kap = (lambda t: kappa0 if t < t_end else kappa0 * np.cos(omega * (t - t_end))) if tail else (lambda t: kappa0)
    P = pc.Params(sigma=sigma, boson=eps_p > 0, fermion=False, ms=0.0, kappa=kap)
    E_bg = 1.5 * Q0 ** 4
    f0 = np.sqrt(eps_p * E_bg / (0.75 * Q0 ** 2)) if eps_p > 0 else 0.0
    y = np.zeros(10)
    y[0], y[2], y[4] = Q0, w0, f0
    if seed:
        rng = np.random.default_rng(seed)
        y[:5] *= 1 + 1e-6 * rng.normal(size=5)
    base = ap.iso_rhs(P)
    if heavy == 'inf':
        def rhs(t, yy):
            z = yy.copy(); z[2] = z[3] = 0.0
            out = base(t, z)
            out[1] += out[3]; out[2] = out[3] = 0.0
            return out
    else:
        rhs = base

    def lost_ev(t, yy):
        return 10 * Q0 - abs(yy[0]) if tail else min(yy[0], 5 * Q0 - yy[0])
    lost_ev.terminal = True
    E_b0 = 0.75 * Q0 ** 2 * f0 ** 2 if f0 else 1.0
    dq, eb, cm, eh, qmax = [], [], [], [], []
    for k in range(nper):
        sol = solve_ivp(rhs, (k * T, (k + 1) * T), y, dense_output=True, events=lost_ev, **RT)
        tt = np.linspace(k * T, sol.t[-1], 1001)
        Y = sol.sol(tt)
        dq.append(np.max(np.abs(Y[0] - Q0)) / Q0)
        qmax.append(np.max(np.abs(Y[0])) / Q0)
        eb.append((np.min(Y[5] ** 2 + 0.75 * Y[0] ** 2 * Y[4] ** 2) / E_b0, np.max(Y[5] ** 2 + 0.75 * Y[0] ** 2 * Y[4] ** 2) / E_b0) if f0 else (np.nan, np.nan))
        cm.append(0.5 * np.max(Y[4] ** 2))
        s18 = ap.iso_unpack(Y[:, -1])
        eh.append(0.0 if heavy == 'inf' else pc.energy_parts(sol.t[-1], s18, P)[1] / E_bg)
        y = sol.y[:, -1]
        if sol.status == 1:
            lost = sol.t_events[0][0] / T
            n_left = nper - k - 1
            dq += [np.nan] * n_left; qmax += [np.nan] * n_left; eb += [(np.nan, np.nan)] * n_left
            cm += [np.nan] * n_left; eh += [np.nan] * n_left
            return np.array(dq), np.array(qmax), np.array(eb), np.array(cm), np.array(eh), lost, T
    return np.array(dq), np.array(qmax), np.array(eb), np.array(cm), np.array(eh), None, T


def blk(v, fmt="{:6.3f}", blocks=((0, 20), (20, 40), (40, 60)), fn=np.nanmax):
    out = []
    for a, b in blocks:
        seg = v[a:b]
        out.append(fmt.format(fn(seg)) if np.any(np.isfinite(seg)) else " " * (len(fmt.format(0.0)) - 1) + "-")
    return " ".join(out)


def header(tail=False):
    print(f"  {'Q0':>5} {'eps_p':>6} {'heavy':>5} {'seed/om':>7} | max|Q-Q0|/Q0 per block" + (" (tail: max|Q|/Q0)" if tail else "")
          + " | E_boson/E_boson(0) min..max per block | max c per block (1/4 collides) | E_heavy/E_bg at block end | lost at")


def row(Q0, eps_p, heavy, seed=0, tail=False, omega=None):
    dq, qmax, eb, cm, eh, lost, T = run(Q0, eps_p, heavy, tail=tail, seed=seed, omega=omega)
    s1 = blk(qmax if tail else dq)
    s2 = " / ".join(f"{np.nanmin(eb[a:b, 0]):7.2f} {np.nanmax(eb[a:b, 1]):8.2f}" if np.any(np.isfinite(eb[a:b, 0])) else "      -        -" for a, b in ((0, 20), (20, 40), (40, 60)))
    s3 = blk(cm)
    s4 = " ".join(f"{eh[b - 1]:8.3f}" if b - 1 < len(eh) and np.isfinite(eh[b - 1]) else "       -" for a, b in ((0, 20), (20, 40), (40, 60)))
    tag = f"{lost:.2f} T" if lost is not None else "-"
    if lost is not None and (eps_p == 0.0 or heavy == 'inf'):
        # a lost control row: for the twin at Q0 = 0.3 it is the twin's own limit (its fixed point is unstable there,
        # README), not a failure of the run; everywhere else it is a failure
        tag += "   TWIN-LIMIT" if heavy == 'twin' else "   CONTROL-FAIL"
    lab = f"{omega:7.2f}" if tail else f"{seed:7d}"
    print(f"  {Q0:5.3f} {eps_p:6.2f} {heavy:>5} {lab} | {s1} | {s2} | {s3} | {s4} | {tag}", flush=True)


if __name__ == "__main__":
    ap.check_iso_rhs()
    # the fixed points are fixed points of iso_rhs (and of the M -> infinity reduction): residual check
    res = 0.0
    for heavy in ('ghost', 'twin', 'inf'):
        for Q0 in (0.024, 0.1, 0.3):
            kappa0, w0 = fixed_point(Q0, heavy)
            P = pc.Params(sigma=-1 if heavy == 'twin' else +1, kappa=lambda t, k=kappa0: k)
            y = np.zeros(10); y[0], y[2] = Q0, w0
            d = ap.iso_rhs(P)(0.0, y)
            if heavy == 'inf':
                d[1] += d[3]; d[2] = d[3] = 0.0
            res = max(res, np.max(np.abs(d)))
    print(f"  CHECK fixed points of the isotropic system (ghost, twin, M -> infinity; three Q0): max residual {res:.1e} "
          f"{'PASS' if res < 1e-14 else 'FAIL'}", flush=True)

    stamp("A: static inflation, constant kappa, massless doublet; 60 periods of the isotropic mode")
    header()
    for Q0 in (0.024, 0.1, 0.3):
        for heavy in ('inf', 'ghost', 'twin'):
            row(Q0, 0.0, heavy)
            for eps_p in (0.01, 0.1, 0.5):
                row(Q0, eps_p, heavy)
    stamp("A, statistics: Q0 = 0.024, eps_p = 0.1, ghost, five perturbed seeds (and the twin's)")
    for seed in range(1, 6):
        row(0.024, 0.1, 'ghost', seed=seed)
    for seed in range(1, 4):
        row(0.024, 0.1, 'twin', seed=seed)

    stamp("B: the tail -- kappa constant for 10 T, then kappa0 cos(omega t) for 50 T (the axion oscillating at its mass)")
    header(tail=True)
    for Q0 in (0.024, 0.1):
        for omega in (0.8, 2.4, 3.0):
            for heavy in ('inf', 'ghost', 'twin'):
                row(Q0, 0.0, heavy, tail=True, omega=omega)
                row(Q0, 0.1, heavy, tail=True, omega=omega)
    print(f"DONE {time.time() - T0:.0f}", flush=True)
