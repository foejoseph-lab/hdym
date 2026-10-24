"""axion_probes_long.py -- the two open items of axion_probes.py that need more than the cloud session's
20 minutes, run by Joseph locally (README, "Probes"; protocol of 16 Oct: self-checks first, every failure a
named line, DONE <seconds> at the end).  Everything in phase space on probe_core.py; g = M = 1, flat space.

  L1  The m_p = 0 exchange over long times, with back-reaction; boson-only and fermion-only columns (the fermion's
      energy is bounded by its norm, so a large probe energy is the boson's); the scalar's Hill coefficient c =
      (g^2/2) f^2 per cycle against the pole-collision value M^2/4 (the ghost tangent's light and heavy poles leave the
      real axis for c > M^2/4, axion_derive.py CHECK 6 and the review of 16 Oct), to see whether the loss coincides
      with the crossing.  The exact (Q, w) isotropic background with the
      axion drive (omega = 2.4 and 3.0, omega^2 phi^2 Pi = 1.7 and 2.0, m_phi/gQ = 100) and both probes massless
      at 10 % of the background's energy, 60 Kapitza periods (axion_probes.py ran 8): ghost / twin / M -> infinity
      side by side.  Per cycle: |Q|max/Q_amp, E_probe/E_probe(0), and the running extremes.  Questions: does the
      ghost case's swing (0.3-2.3 Q_amp over 8 cycles) stay bounded, grow secularly, or lose the orbit; does the
      twin or M -> infinity case ever swing.
  L2  The m_s ~ M case with back-reaction.  The undriven isotropic orbit (no axion, Q_amp = 1, T = 5.244) with the
      scalar doublet coupled and back-reacting (the full probe_core system on x = y = z), m_s in (0.8, 0.9, 1.0,
      1.1, 1.2), scalar at eps_p = 0.01 and 0.1 of the background energy, 100 periods: the scalar's energy and
      the background amplitude per period (does the preheating saturate, and where), and the anisotropic tangent
      eps'''' + eps'' + (f^2/2) eps = 0 driven by the actual f(t) (renormalised each period) in the ghost, twin
      and M -> infinity versions.  Control rows: the same with the scalar absent (must be steady).

Controls that must be steady (a BLOW-UP or a runaway there is a no-go, not a result): L1 'inf' rows with no probe;
L2 rows with the scalar absent.

Result (README, "The mechanism"; the first readings of 16 Oct -- "a bounded three-way exchange", "a massless probe
lowers the loss threshold from 2 to at most 1.7", "a three-body pump needing both probes", "the loss time as the
number to report" -- were withdrawn the same night): the loss coincides with c = (g^2/2) f^2 crossing M^2/4 (the
pole collision), the m_D = 0 fermion is a constant force giving Q a nonzero mean, and 'boson+F' reproduces the pump;
the loss is realisation-dependent (axion_probes_spread.py), a fraction, not a threshold or a rate.

Usage (Joseph, from Spyder, output outside the repository; ~30-60 min expected):
    %cd C:\\GitHub\\hdym
    !"{sys.executable}" -B -u axion_su2\\code\\axion_probes_long.py > C:\\hdym_runs\\axion_probes_long.txt 2>&1
Added 2026-10-16; c(t) columns and the one-probe rows added the same night; the Hessian term and the constant-force
column after the second review (the m_D = 0 fermion is a constant force; L1's coefficient carries 6 g^2 Q^2).
"""
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import probe_core as pc  # noqa: E402
import axion_probes as ap  # noqa: E402

T0 = time.time()
RT = dict(method='DOP853', rtol=1e-9, atol=1e-14)
K = ellipk(0.5)


def stamp(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def fmt(v):
    return " ".join(f"{x:5.2f}" if np.isfinite(x) else "    -" for x in v)


# ---------------------------------------------------------------- L1
def long_iso(omega, crit, heavy, column='both', mp_over_gQ=0.0, ncyc=60):
    """ap.run_iso with 60 cycles; also the scalar's Hill coefficient c = (g^2/2) f^2 per cycle (max), to be read
    against the pole-collision value M^2/4 (review of 16 Oct, point 1: the ghost tangent's poles leave the real axis
    for c > M^2/4), and the first cycle at which it crosses."""
    sigma = -1 if heavy == 'twin' else +1
    out = ap.run_iso(omega, crit, column, sigma, mp_over_gQ, ncyc=ncyc, no_ghost=(heavy == 'inf'), return_sol=True)
    phi, lam, T, per, eb, nf, eh, lost, sol = out
    cmax, hmax = [], []
    for k in range(ncyc):
        a, b = k * T, min((k + 1) * T, sol.t[-1])
        if b > a:
            Y = sol.sol(np.linspace(a, b, 2001))
            cmax.append(0.5 * np.max(Y[4] ** 2))
            hmax.append(6.0 * np.max(Y[0] ** 2))          # the Hessian term of the isotropic tangent, 6 g^2 Q^2
        else:
            cmax.append(np.nan); hmax.append(np.nan)
    cmax, hmax = np.array(cmax), np.array(hmax)
    cross_c = next((k + 1 for k in range(ncyc) if np.isfinite(cmax[k]) and cmax[k] > 0.25), None)
    cross_h = next((k + 1 for k in range(ncyc) if np.isfinite(hmax[k]) and hmax[k] > 0.25), None)
    cross_s = next((k + 1 for k in range(ncyc) if np.isfinite(cmax[k]) and cmax[k] + hmax[k] > 0.25), None)
    return per, eb, nf, lost, T, cmax, hmax, (cross_c, cross_h, cross_s)


# ---------------------------------------------------------------- L2
def coupled_scalar(m_s, eps_p, heavy, nper=100, Q_amp=1.0):
    """Undriven orbit with the scalar coupled (back-reacting) on x = y = z, plus the anisotropic tangent driven by the
    actual f(t).  State: iso_pack layout (10) + tangent (4 or 2).  Returns per-period |Q|max/Q_amp, E_scalar/E_scalar(0),
    and the tangent's ln growth per period (list), with the scalar at rest at its turning point."""
    sigma = -1 if heavy == 'twin' else +1
    P = pc.Params(sigma=sigma, boson=eps_p > 0, fermion=False, ms=m_s)
    T = 4 * K / (np.sqrt(2) * Q_amp)
    E_bg = 1.5 * Q_amp ** 4
    f0 = np.sqrt(eps_p * E_bg / (m_s ** 2 + 0.75 * Q_amp ** 2)) if eps_p > 0 else 0.0
    base = ap.iso_rhs(P)
    dim = 2 if heavy == 'inf' else 4

    def rhs(t, y):
        yb = y[:10].copy()
        if heavy == 'inf':
            yb[2] = yb[3] = 0.0
        db = base(t, yb)
        if heavy == 'inf':
            db[1] += db[3]
            db[2] = db[3] = 0.0
        f = y[4]
        c = 0.5 * f * f
        e = y[10:]
        de = np.empty(dim)
        if heavy == 'ghost':
            de[:] = [e[1], e[2], e[3], -e[2] - c * e[0]]
        elif heavy == 'twin':
            de[:] = [e[1], e[2] - 2 * c * e[0], e[3], -e[2] - c * e[0]]
        else:
            de[:] = [e[1], -c * e[0]]
        return np.concatenate([db, de])
    y = np.zeros(10 + dim)
    y[0], y[4] = Q_amp, f0
    y[10:] = 1.0 / np.sqrt(dim)
    E0 = (m_s ** 2 + 0.75 * 3 * Q_amp ** 2 / 3) * f0 ** 2 if f0 else 1.0   # fdot^2 + (ms^2 + |q|^2/4) f^2 at t = 0, |q|^2 = 3Q^2
    qmax, es, grow, cmax = [], [], [], []
    for k in range(nper):
        sol = solve_ivp(rhs, (k * T, (k + 1) * T), y, dense_output=True, **RT)
        tt = np.linspace(k * T, (k + 1) * T, 401)
        Y = sol.sol(tt)
        qmax.append(np.max(np.abs(Y[0])) / Q_amp)
        cmax.append(0.5 * np.max(Y[4] ** 2))
        Ef = np.mean(Y[5] ** 2 + (m_s ** 2 + 0.75 * Y[0] ** 2) * Y[4] ** 2)
        es.append(Ef / E0 if f0 else 0.0)
        y = sol.y[:, -1]
        n = np.linalg.norm(y[10:])
        grow.append(np.log(n))
        y[10:] /= n
        if not np.all(np.isfinite(y)) or abs(y[0]) > 10 * Q_amp:
            qmax += [np.nan] * (nper - k - 1); es += [np.nan] * (nper - k - 1); grow += [np.nan] * (nper - k - 1)
            cmax += [np.nan] * (nper - k - 1)
            break
    return np.array(qmax), np.array(es), np.array(grow), np.array(cmax)


if __name__ == "__main__":
    ap.check_iso_rhs()      # CHECK line: iso_rhs against probe_core.rhs
    stamp("L1: m_p = 0, 60 Kapitza periods, ghost / twin / M -> infinity; c = (g^2/2) f^2 and the Hessian term 6 g^2 Q^2 against 1/4")
    print("  L1 is isotropic (no anisotropic eps in the state), so the coefficient that can collide is the isotropic tangent's,")
    print("  which carries the Hessian 6 g^2 Q^2 as well as c; both are logged (max per 20-cycle block) and the first cycle at")
    print("  which c, 6g^2Q^2 and their sum exceed 1/4.  'boson+F' = boson plus the constant force the m_D = 0 fermion exerts")
    print("  (it sits in a sigma_x eigenstate for all t), no spinor dynamics: the test of whether the pump is offset-plus-boson.")
    print(f"  {'omega':>5} {'w2phi2Pi':>8} {'heavy':>5} {'col':>8} | |Q|max/Q_amp per block | E_boson/E_boson(0) per block | max c per block | max 6Q^2 per block | >1/4 first (c, 6Q^2, sum) | lost at")
    for omega in (2.4, 3.0):
        for crit in (1.7, 2.0):
            for heavy in ('inf', 'ghost', 'twin'):
                if heavy == 'inf':
                    cols = ('none', 'both')
                elif heavy == 'ghost':
                    cols = ('both', 'boson', 'fermion', 'boson+F') if crit == 1.7 else ('both', 'boson+F')
                else:
                    cols = ('both',)
                for column in cols:
                    per, eb, nf, lost, T, cmax, hmax, cross = long_iso(omega, crit, heavy, column)
                    per, eb = np.array(per), np.array(eb)
                    blocks = [(0, 20), (20, 40), (40, 60)]
                    s1 = " / ".join(f"{np.nanmin(per[a:b]):5.2f} {np.nanmax(per[a:b]):5.2f}" if np.any(np.isfinite(per[a:b])) else "    -     -" for a, b in blocks)
                    s2 = " / ".join(f"{np.nanmin(eb[a:b]):7.2f} {np.nanmax(eb[a:b]):7.2f}" if np.any(np.isfinite(eb[a:b])) else "      -       -" for a, b in blocks)
                    s3 = " ".join(f"{np.nanmax(cmax[a:b]):6.3f}" if np.any(np.isfinite(cmax[a:b])) else "     -" for a, b in blocks)
                    s4 = " ".join(f"{np.nanmax(hmax[a:b]):6.3f}" if np.any(np.isfinite(hmax[a:b])) else "     -" for a, b in blocks)
                    s5 = " ".join(f"{x if x else '-':>4}" for x in cross)
                    tag = "" if lost is None else f"{lost:.2f} T"
                    if heavy == 'inf' and column == 'none' and lost is not None:
                        tag += "   CONTROL-FAIL"
                    print(f"  {omega:5.2f} {crit:8.2f} {heavy:>5} {column:>8} | {s1} | {s2} | {s3} | {s4} | {s5} | {tag or '-'}", flush=True)
    stamp("L2: m_s ~ M with back-reaction, undriven orbit, 100 periods")
    print(f"  {'m_s':>4} {'eps_p':>6} {'heavy':>5} | |Q|max/Q_amp at periods 10, 30, 60, 100 | E_s/E_s(0) at 10, 30, 60, 100 | max c at 10, 30, 60, 100 (collision at 0.25) | tangent growth/period over 10-30, 30-60, 60-100")
    for m_s in (0.8, 0.9, 1.0, 1.1, 1.2):
        for eps_p in (0.0, 0.01, 0.1):
            for heavy in (('inf', 'ghost') if eps_p == 0.0 else ('inf', 'ghost')):
                qmax, es, grow, cmax = coupled_scalar(m_s, eps_p, heavy)
                pick = [9, 29, 59, 99]
                s1 = fmt(qmax[pick]); s2 = fmt(es[pick])
                s4 = " ".join(f"{x:6.3f}" if np.isfinite(x) else "     -" for x in cmax[pick])
                g = [np.nanmean(grow[a:b]) for a, b in ((10, 30), (30, 60), (60, 100))]
                s3 = " ".join(f"{x:7.3f}" if np.isfinite(x) else "      -" for x in g)
                tag = ""
                if eps_p == 0.0 and (not np.all(np.isfinite(qmax)) or np.nanmax(qmax) > 2.0):
                    tag = "   CONTROL-FAIL"
                print(f"  {m_s:4.1f} {eps_p:6.2f} {heavy:>5} | {s1} | {s2} | {s4} | {s3}{tag}", flush=True)
    print("  (twin rows omitted: on the Q_amp = 1 orbit the twin's heavy coefficient M^2 - Hess V is negative, so it is")
    print("   not a control there; the ghost-vs-M->infinity comparison is the control at this amplitude)")
    print(f"DONE {time.time() - T0:.0f}", flush=True)
