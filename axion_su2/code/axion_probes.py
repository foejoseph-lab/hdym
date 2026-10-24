"""axion_probes.py -- the axion background with the real-particle probes of probe_core.py (16 Oct):
the boson (doublet modulus f), the fermion (singlet two-level zeta), both, each coupled to the ghost
through the gauge field, with the twin as the sign control.  Three parts, all on the exact Lee-Wick
isotropic background in the (q, w) form, flat space, g = M = 1.  Revised after the review of 16 Oct
(twin tangent sign, probe-mass axis, energy split, fermion observable, twin initial data, Part C masses).

Two regimes, labelled because they are opposite ends of one axis:
  Part A / C:  gQ_amp = M/100 (m_phi/gQ = 100 as in iso_lame_gap.py).  The probe's coupling to the field is
               (3/4) g^2 Q^2 (scalar) and (3/2) g Q (fermion mixing), so the axis that matters is m_p / gQ_amp, not
               m_p / M: at m_p = 100 gQ_amp = omega M (2.4M, 3M; labelled "= M" in the table title, corrected 17 Oct)
               the probe is decoupled by construction (coupling ~1e-4 of m_p^2, force ~1e-5 of
               2Q^3) and reproduces the no-probe rows to every digit -- a hierarchy statement, not a ghost one.
               Rows are at m_p / gQ_amp = 0, 1, 3 and 100 (= omega M, the decoupled control); the probe decouples
               already at m_p ~ 3 gQ_amp (README, "Results" A).
  Part B:      Q_amp = 1, i.e. gQ_amp = M: the regime where the probe's induced mass is of order M.

  A  The isotropic threshold with the probes.  Rows: omega = m_phi/M in (2.4, 3.0) [off the ghost's omega = 2M
     tongue] and 0.8 [M > m_phi control]; omega^2 phi^2 Pi in (1.0, 1.7, 2.0, 2.3) (the no-probe threshold is 2.0).
     Columns: none / ghost + boson / ghost + fermion / ghost + both, then the twin and the M -> infinity control.
     Each run: 8 Kapitza periods, per cycle |Q|max/Q_amp, the boson's energy relative to its start (its own
     fdot^2 + (ms^2 + |q|^2/4) f^2), the fermion's adiabatic number n_- (weight on the negative-energy level of
     h(Q), 0 at the start, the first-quantised "pair" fraction -- its energy is bounded by construction and says
     nothing; and it is meaningful only for m_p >~ gQ_amp: at m_p = 0 the two levels cross whenever Q
     passes through zero and the adiabatic basis flips, so n_- ~ 1/2 there is level crossing, not production),
     and the heavy sector's energy E_heavy (ghost: -|wdot|^2/2 - M^2|w|^2/2; twin: +) over the background's
     energy, averaged over the cycle (with a drive on the pieces of the energy are not separately conserved:
     the split is bookkeeping of where the drive's work went); blow-up at |Q| = 10 Q_amp.  The probe starts with 10 % of the
     background's Yang-Mills energy E_bg = 1.5 Q_amp^4 (V at the turning point; the ghost sector's negative share,
     O(kappa^2 Q^4 Pi^2) ~ 1e-3 E_bg here, is not counted in that reference).  The twin's dressed initial wdot has
     the opposite sign from the ghost's (its w-equation is w'' + (M^2 - 6Q^2) w = -kappa Q^2).
  B  The anisotropic sector with the scalar: the formerly free direction obeys
        ghost  eps''''/M^2 + eps'' + c eps = 0,   twin  zeta'' = -M^2 zeta + c eps, eps'' = M^2 zeta - 2 c eps,
        M -> infinity  eps'' + c eps = 0,         c = (g^2/2) f(t)^2
     (probe_derive.py CHECK 5, 5b) with f a test field on the undriven orbit (no axion: the scalar's own drive);
     the hand-coded tangents are checked against probe_core.tangent_rhs (a finite-difference linearisation of the
     core) at random points before the scan; 30 Q periods renormalised each period; the test scalar's own growth
     reported (where it is not ~0 the background preheats it and the test-field setup is not valid: see
     axion_probes_long.py L2 for the back-reacting version).
  C  On the ghost's first tongue (omega = 2M, omega^2 phi^2 Pi = 1.0: the background steady, the ghost pair
     growing) with the ghost and the twin, at m_p / gQ_amp = 0, 1, 3, 50 and 100.  (Label correction, 17 Oct: at
     omega = 2, Q_amp = 0.02, so these are m_p = M and 2M, not M/2 and M as the table title says; the fermion at m_p = M
     has level spacing 2M = the drive frequency.  Likewise Part A's "m_p = M (= 100 gQ_amp)" is m_p = omega M.  The
     table titles in the committed log are left as run; axion_probes_fourier.py takes the row apart.)

Usage (from the repository root, ~10 min):
    python axion_su2/code/axion_probes.py > axion_su2/results/axion_probes.txt
Added 2026-10-16; revised the same night.
"""
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import probe_core as pc  # noqa: E402

T0 = time.time()
RT = dict(method='DOP853', rtol=1e-9, atol=1e-14)
K = ellipk(0.5)
RATIO = 100.0       # m_phi / g Q_amp
EPS_P = 0.1         # probe energy / background energy at the start
NCYC = 8


def stamp(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def kappa_fun(omega, phi):
    return lambda t: -omega * phi * np.sin(omega * t)


def iso_rhs(P):
    """The isotropic truncation x = y = z of probe_core.rhs (exact: it is a consistent truncation), as a lean
    10-d system (Q, Q', w, w', f, f', Re zeta, Im zeta) for speed; checked against probe_core.rhs below."""
    g, M, sig, ms, mD = P.g, P.M, P.sigma, P.ms, P.mD
    c = pc.C_FERMION
    bos, fer = P.boson, P.fermion

    def rhs(t, y):
        Q, Q1, w, w1, f, f1 = y[:6]
        zr, zi = y[6:8], y[8:10]
        ka = P.kappa(t)
        F = ka * Q * Q
        if bos:
            F -= 0.5 * g * g * Q * f * f
        if fer:
            F += g * (-(c / 3.0) * 2.0 * (zr[0] * zr[1] + zi[0] * zi[1]))
        gV = 2.0 * g * g * Q ** 3
        Hw = 6.0 * g * g * Q * Q * w
        if sig == +1:
            Qdd = M * M * w - gV
            wdd = -Hw - M * M * w + F
        else:
            wdd = Hw - M * M * w - F
            Qdd = -gV - M * M * w - 2.0 * wdd
        out = np.zeros(10)
        out[0], out[1], out[2], out[3] = Q1, Qdd, w1, wdd
        if bos:
            out[4], out[5] = f1, -(ms * ms + 0.75 * g * g * Q * Q) * f
        if fer:
            h01 = c * g * Q
            # i zeta' = h zeta,  h = [[mD, h01],[h01, -mD]]:  zeta' = -i h zeta
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


def check_iso_rhs():
    rng = np.random.default_rng(3)
    err = 0.0
    for sig in (+1, -1):
        for bos in (False, True):
            for fer in (False, True):
                P = pc.Params(sigma=sig, boson=bos, fermion=fer, ms=0.7, mD=0.6, kappa=kappa_fun(2.4, 0.3))
                for _ in range(20):
                    y = rng.normal(size=10)
                    if not bos:
                        y[4:6] = 0
                    if not fer:
                        y[6:10] = 0
                    tt = rng.uniform(0, 5)
                    a = iso_rhs(P)(tt, y)
                    b = iso_pack(pc.rhs(tt, iso_unpack(y), P))
                    err = max(err, np.max(np.abs(a - b)) / max(1.0, np.max(np.abs(b))))
    print(f"  iso_rhs against probe_core.rhs on x = y = z, all switches, random points: max rel err {err:.1e} "
          f"{'PASS' if err < 1e-12 else 'FAIL'}", flush=True)


def initial(omega, phi, Q_amp, P, m_p):
    """Isotropic data: Q at its turning point, w on the dressed particular solution (iso_lame_gap.py), the
    probes at EPS_P of the background energy 1.5 Q_amp^4 (V on x = y = z at the turning point)."""
    s = np.zeros(18)
    s[0:3] = Q_amp
    w2 = omega * omega
    s[9:12] = P.sigma * w2 * phi * Q_amp ** 2 / (w2 - 1.0)      # twin: opposite sign (docstring)
    E_bg = 1.5 * Q_amp ** 4
    if P.boson:
        meff2 = m_p ** 2 + 0.75 * Q_amp ** 2
        s[12] = np.sqrt(EPS_P * E_bg / meff2)
    if P.fermion:
        h = pc.fermion_h(Q_amp, P)
        ev, U = np.linalg.eigh(h)
        pos = U[:, 1]
        n = EPS_P * E_bg / ev[1]
        zeta = np.sqrt(n) * pos
        s[14:16], s[16:18] = zeta.real, zeta.imag
    return s


def probe_energy(t, s, P):
    q, v, w, u, f, fd, zeta = pc.unpack(s)
    E = 0.0
    if P.boson:
        E += fd * fd + (P.ms ** 2 + 0.25 * (q @ q)) * f * f
    if P.fermion:
        Q = (q[0] + q[1] + q[2]) / 3.0
        E += np.real(np.conj(zeta) @ (pc.fermion_h(Q, P) @ zeta))
    return E


def run_iso(omega, crit, column, sigma, mp_over_gQ, ncyc=NCYC, no_ghost=False, return_sol=False):
    """mp_over_gQ: the probe mass in units of g Q_amp (both probes get the same mass)."""
    w2 = omega * omega
    Pi = 1.0 / (w2 * w2 - w2)
    phi = np.sqrt(crit / (w2 * Pi)) if w2 * Pi > 0 else np.sqrt(abs(crit / (w2 * Pi)))
    Q_amp = omega / RATIO
    m_p = mp_over_gQ * Q_amp
    lam = 2.0 - crit if w2 * Pi > 0 else 2.0 + crit     # omega < 1: Pi < 0, the drive stiffens
    phi_eff = -phi * w2 * Pi
    X = Q_amp * (1 + 2 * phi_eff ** 2 / (2 + phi_eff ** 2)) ** 0.25
    T = 4 * K / (np.sqrt(lam) * X) if lam > 0 else 4 * K / (np.sqrt(2) * X)
    P = pc.Params(sigma=sigma, boson=column in ('boson', 'both', 'boson+F'), fermion=column in ('fermion', 'both'),
                  ms=m_p, mD=m_p, kappa=kappa_fun(omega, phi))
    s0 = initial(omega, phi, Q_amp, P, m_p)
    E_p0 = probe_energy(0.0, s0, P)
    f_rhs = iso_rhs(P)
    if column == 'boson+F':
        # the m_D = 0 fermion in the positive eigenstate of h(Q_amp) is a sigma_x eigenstate for all t (h = cgQ sigma_x),
        # so its current is constant: replace it by that constant force g j0 (review of 16 Oct, point 1)
        Pf = pc.Params(sigma=sigma, boson=False, fermion=True, ms=m_p, mD=m_p, kappa=kappa_fun(omega, phi))
        sf = initial(omega, phi, Q_amp, Pf, m_p)
        F0 = P.g * pc.fermion_current(sf[14:16] + 1j * sf[16:18])
        base_f = f_rhs

        def f_rhs(t, y):
            out = base_f(t, y)
            out[3] += F0 if sigma == +1 else -F0       # F enters w'' (ghost) / -F (twin), as in iso_rhs
            return out
    if no_ghost:
        # M -> infinity: Q'' = -2 Q^3 + F (the ghost-free Hill background); w frozen at zero
        base = f_rhs

        def f_rhs(t, y):
            yy = y.copy(); yy[2] = yy[3] = 0.0
            out = base(t, yy)
            out[1] = out[1] + out[3]          # at w = 0, w'' = F, so Q'' = -gradV + F
            out[2] = out[3] = 0.0
            return out
        s0[6:12] = 0.0

    def blow(t, y):
        return abs(y[0]) - 10 * Q_amp
    blow.terminal = True
    sol = solve_ivp(f_rhs, (0.0, ncyc * T), iso_pack(s0), dense_output=True, events=blow, **RT)
    E_bg = 1.5 * Q_amp ** 4
    E_b0 = (P.ms ** 2 + 0.75 * Q_amp ** 2) * s0[12] ** 2 if P.boson else 0.0
    per, eb, nf, eh = [], [], [], []
    for k in range(ncyc):
        a, b = k * T, min((k + 1) * T, sol.t[-1])
        if b > a:
            tt = np.linspace(a, b, 2001)
            Y = sol.sol(tt)
            per.append(np.max(np.abs(Y[0])) / Q_amp)
            if P.boson:
                eb.append(np.mean((Y[5] ** 2 + (P.ms ** 2 + 0.75 * Y[0] ** 2) * Y[4] ** 2)[::40]) / E_b0)
            else:
                eb.append(np.nan)
            if P.fermion:
                nf.append(np.mean([pc.fermion_adiabatic_number(Y[0, i], Y[6:8, i] + 1j * Y[8:10, i], P) for i in range(0, 2001, 40)]))
            else:
                nf.append(np.nan)
            eh.append(0.0 if no_ghost else np.mean([pc.energy_parts(t_, iso_unpack(Y[:, i]), P)[1] for i, t_ in enumerate(tt[::40])]) / E_bg)
        else:
            per.append(np.nan); eb.append(np.nan); nf.append(np.nan); eh.append(np.nan)
    blew = sol.status == 1
    lost = sol.t_events[0][0] / T if blew else None
    if return_sol:
        return phi, lam, T, per, eb, nf, eh, lost, sol
    return phi, lam, T, per, eb, nf, eh, lost


def table(rows, title):
    print(f"\n{title}")
    print(f"  {'omega':>5} {'w2phi2Pi':>8} {'col':>7} {'heavy':>5} {'m_p/gQ':>6} | |Q|max/Q_amp per cycle ({NCYC})" + " " * 15
          + "| E_boson/E_boson(0) per cycle" + " " * 15 + "| fermion n_- per cycle" + " " * 26
          + "| E_heavy/E_bg, cycle mean" + " " * 23 + "| lost at")
    f2 = lambda v: " ".join(f"{x:5.2f}" if np.isfinite(x) else "    -" for x in v)
    f3 = lambda v: " ".join(f"{x:5.3f}" if np.isfinite(x) else "    -" for x in v)
    for row in rows:
        omega, crit, column, sigma, mp = row[:5]
        no_ghost = len(row) > 5 and row[5]
        phi, lam, T, per, eb, nf, eh, lost = run_iso(omega, crit, column, sigma, mp, no_ghost=no_ghost)
        hv = 'inf' if no_ghost else ('ghost' if sigma == 1 else 'twin')
        print(f"  {omega:5.2f} {crit:8.2f} {column:>7} {hv:>5} {mp:6.0f} | {f2(per)} | {f2(eb)} | {f3(nf)} | {f3(eh)} | "
              + (f"{lost:.2f} T" if lost is not None else "-"), flush=True)


# ------------------------------------------------------------------ Part B: anisotropic sector with the scalar
_TANGENT_CHECKED = set()


def check_tangent(rhs_local, heavy):
    """The hand-coded tangent of aniso_scalar against probe_core.tangent_rhs (finite differences of the core's
    full rhs along (1, -1, 0)): at random points, and at 20 points along the actual orbit-plus-scalar trajectory
    (m_s = 0, eps_p = 0.5, one period), where c = f^2/2 swings through zero and a sign error is invisible wherever
    f is small; once per heavy type; prints a CHECK line."""
    if heavy in _TANGENT_CHECKED:
        return
    _TANGENT_CHECKED.add(heavy)
    rng = np.random.default_rng(5)
    P = pc.Params(sigma=+1 if heavy == 'ghost' else -1, boson=True, ms=0.0)
    T = 4 * K / np.sqrt(2)
    f0 = np.sqrt(0.5 * 1.5 / 0.75)
    orb = solve_ivp(lambda t, y: [y[1], -2 * y[0] ** 3, y[3], -(0.75 * y[0] ** 2) * y[2]], (0, T), [1.0, 0.0, f0, 0.0],
                    dense_output=True, **RT)
    pts = [tuple(rng.normal(size=4)) for _ in range(20)] + [tuple(orb.sol(t_)) for t_ in np.linspace(0, T, 21)[:-1]]
    err = 0.0
    for (Q, Qd, f, fd) in pts:
        e = rng.normal(size=4)
        d_local = rhs_local(0.0, np.array([Q, Qd, f, fd, *e]))[4:]
        s = np.zeros(18); s[0:3] = Q; s[3:6] = Qd; s[12], s[13] = f, fd
        ds = np.zeros(18); v = np.array([1.0, -1.0, 0.0])
        # local state: ghost (eps, d eps, d2 eps, d3 eps) with d2 eps = zeta; twin (eps, d eps, zeta, d zeta)
        ds[0:3], ds[3:6], ds[6:9], ds[9:12] = e[0] * v, e[1] * v, e[2] * v, e[3] * v
        d = pc.tangent_rhs(0.0, s, ds, P)
        d_core = np.array([d[0], d[3], d[6], d[9]])
        err = max(err, np.max(np.abs(d_local - d_core)) / max(1.0, np.max(np.abs(d_core))))
    print(f"  CHECK tangent ({heavy}) against probe_core.tangent_rhs at 20 random points and 20 along the orbit: "
          f"max rel err {err:.1e} {'PASS' if err < 1e-6 else 'FAIL'}", flush=True)


def aniso_scalar(m_s, eps_p, heavy='ghost', nper=30, Q_amp=1.0):
    """Undriven isotropic orbit Q'' = -2 Q^3 from (Q_amp, 0); scalar f'' = -(m_s^2 + 3/4 Q^2) f from its turning point
    with eps_p of the background energy; the anisotropic tangent with the scalar's mass c(t) = f^2/2 (probe_derive.py
    CHECK 5), no back-reaction (a linear tangent on a fixed orbit):
        ghost   eps'''' + eps'' + c eps = 0   (M = 1; 4-d)
        twin    the same direction of the sigma = -1 system of probe_core: zeta'' = -zeta - c eps,
                eps'' = zeta - 2 c eps        (4-d)
        inf     eps'' + c eps = 0             (M -> infinity, Hill)
    The tangent is renormalised every period; returns the mean growth of ln|state| per period over periods 10-30
    and the total over 30."""
    T = 4 * K / (np.sqrt(2) * Q_amp)
    E_bg = 1.5 * Q_amp ** 4
    f0 = np.sqrt(eps_p * E_bg / (m_s ** 2 + 0.75 * Q_amp ** 2))
    dim = 2 if heavy == 'inf' else 4

    def rhs(t, y):
        Q, Qd, f, fd = y[:4]
        e = y[4:]
        out = np.empty_like(y)
        out[0], out[1] = Qd, -2 * Q ** 3
        out[2], out[3] = fd, -(m_s ** 2 + 0.75 * Q * Q) * f
        c = 0.5 * f * f
        if heavy == 'ghost':      # (eps, eps', eps'', eps''')
            out[4], out[5], out[6], out[7] = e[1], e[2], e[3], -e[2] - c * e[0]
        elif heavy == 'twin':     # (eps, eps', zeta, zeta'): zeta'' = -zeta + c eps, eps'' = zeta - 2 c eps
            out[4], out[5], out[6], out[7] = e[1], e[2] - 2 * c * e[0], e[3], -e[2] + c * e[0]
        else:
            out[4], out[5] = e[1], -c * e[0]
        return out
    if heavy != 'inf':
        check_tangent(rhs, heavy)
    y = np.concatenate([[Q_amp, 0.0, f0, 0.0], np.ones(dim) / np.sqrt(dim)])
    lns = [0.0]
    lnf = []
    for k in range(nper):
        sol = solve_ivp(rhs, (k * T, (k + 1) * T), y, dense_output=True, **RT)
        lnf.append(np.log(np.max(np.abs(sol.sol(np.linspace(k * T, (k + 1) * T, 401))[2])) / f0))
        y = sol.y[:, -1]
        n = np.linalg.norm(y[4:])
        lns.append(lns[-1] + np.log(n))
        y[4:] /= n
    lns = np.array(lns)
    slope = np.polyfit(np.arange(10, nper + 1), lns[10:], 1)[0]
    # the scalar is a test field on the fixed orbit: if the orbit pumps it (its own parametric resonance with Q^2),
    # its growth is reported too, because then c(t) grows and every column inherits it
    fslope = np.polyfit(np.arange(10, nper), np.array(lnf)[10:], 1)[0]
    return slope, lns[-1], fslope


if __name__ == "__main__":
    check_iso_rhs()
    stamp("Part A: the isotropic threshold with the probes (regime gQ_amp = M/100; m_p in units of g Q_amp)")
    rows = [(om, cr, 'none', +1, 0) for om in (2.4, 3.0) for cr in (1.0, 1.7, 2.0, 2.3)]
    rows += [(om, cr, col, +1, mp) for om in (2.4, 3.0) for cr in (1.0, 1.7, 2.0, 2.3) for mp in (0, 1, 3) for col in ('boson', 'fermion', 'both')]
    rows += [(om, cr, 'both', +1, 100) for om in (2.4, 3.0) for cr in (1.7, 2.3)]
    rows += [(0.8, cr, col, +1, mp) for cr in (1.0, 2.3) for (col, mp) in (('none', 0), ('both', 0), ('both', 1))]
    # title left as run; "m_p = M (= 100 gQ_amp)" is m_p = omega M (2.4M, 3M) -- label correction of 17 Oct, docstring
    table(rows, "A1  ghost; columns none / boson / fermion / both at m_p/gQ_amp = 0, 1, 3; the decoupled control m_p = M (= 100 gQ_amp)")
    stamp("Part A2: the twin and the M -> infinity control")
    rows = [(om, cr, 'both', -1, mp) for om in (2.4, 3.0) for cr in (1.0, 1.7, 2.0, 2.3) for mp in (0, 1)]
    rows += [(om, cr, 'none', -1, 0) for om in (2.4, 3.0) for cr in (1.7, 2.3)]
    rows += [(om, cr, col, +1, mp, True) for om in (2.4, 3.0) for cr in (1.7, 2.0) for (col, mp) in (('none', 0), ('both', 0), ('both', 1))]
    table(rows, "A2  twin (both probes at m_p/gQ = 0, 1; none) and M -> infinity (none / both at 0, 1)")

    stamp("Part B: the anisotropic sector with the scalar's (g^2/2) f^2 on the undriven orbit (regime gQ_amp = M)")
    print("  undriven orbit Q_amp = 1 (T = %.3f), scalar at eps_p of the background energy, 30 periods; growth of" % (4 * K / np.sqrt(2)))
    print("  ln|eps| per period over periods 10-30 (the free sector drifts as ln t: ~0.05 per period here)")
    print("  2 m_eff is the scalar's f^2 frequency (m_eff^2 = m_s^2 + 3/4 <Q^2>, <Q^2> = 0.457 Q_amp^2); the ghost pair's")
    print("  first tongue is at 2 m_eff = 2M, i.e. m_eff = M")
    print("  'f rate' is the test scalar's own growth per period on the fixed orbit (its parametric resonance with Q^2);")
    print("  where it is not ~0 the scalar is being preheated by the background and every column inherits that")
    print(f"  {'m_s':>4} {'2m_eff':>6} {'eps_p':>6} {'f rate':>7} | {'ghost rate':>10} {'ln(30)':>8} | {'twin rate':>10} {'ln(30)':>8} | {'M->inf rate':>11} {'ln(30)':>8}")
    for m_s in (0.0, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0):
        meff2 = 2 * np.sqrt(m_s ** 2 + 0.75 * 0.457)
        for eps_p in (0.01, 0.1, 0.5):
            rg, lg, fr = aniso_scalar(m_s, eps_p, 'ghost')
            rt, lt, _ = aniso_scalar(m_s, eps_p, 'twin')
            rh, lh, _ = aniso_scalar(m_s, eps_p, 'inf')
            print(f"  {m_s:4.1f} {meff2:6.2f} {eps_p:6.2f} {fr:7.3f} | {rg:10.4f} {lg:8.2f} | {rt:10.4f} {lt:8.2f} | {rh:11.4f} {lh:8.2f}", flush=True)

    stamp("Part C: on the ghost's tongue omega = 2M, omega^2 phi^2 Pi = 1.0 (regime gQ_amp = M/100)")
    rows = [(2.0, 1.0, col, sg, mp) for sg in (+1, -1) for mp in (0, 1, 3, 50, 100) for col in ('boson', 'fermion', 'both')]
    # title left as run; at omega = 2, Q_amp = 0.02, so 50 and 100 gQ_amp are m_D = M and 2M, not M/2 and M (17 Oct)
    table(rows, "C   omega = 2 (tongue), ghost and twin, m_p/gQ_amp = 0, 1, 3, 50 (= M/2), 100 (= M)")
    stamp("done")
