"""axion_probes_fourier.py -- Part C's asymmetry taken apart (README "Probes", second review, point 4): on the ghost's
tongue (omega = 2M, omega^2 phi^2 Pi = 1.0) a fermion of mass m_D = M (level spacing 2M = the drive frequency) is driven into its negative-energy level by the
twin background (n_- -> 0.62 in 8 cycles) and hardly by the ghost one (0.009).  The fermion couples to Q(t) alone,
h = m_D sigma_z + (3/2) g Q sigma_x, so the difference is a property of the two Q(t) waveforms; a two-level system with
spacing Omega = 2 sqrt(m_D^2 + (3gQ/2)^2) ~ M responds to the Fourier component of Q at Omega.  Three steps:

  1  Record Q(t) from the coupled ghost and twin runs of Part C (fermion column, m_p = 50 gQ_amp = M at omega = 2 -- Part C's
     row label '50 (= M/2)' was wrong by the factor Q_amp = 0.02: the rows are M and 2M; also the no-probe runs), 8 Kapitza periods, and report the Fourier amplitude |Q^(nu)| (two-sided, normalised so a pure
     A cos(nu t) gives A) at nu = Omega (the level spacing), at the drive 2M, at the ghost's M, and the mean <Q>;
     also |Q|max.  If the twin has more power at Omega, the asymmetry is amplitude; if not, it is something else.
  2  The fermion as a test field on each recorded waveform (same initial state, no back-reaction): n_- at the end.
     It must reproduce the coupled rows (the fermion's back-reaction is ~1e-7 here) -- the check that the asymmetry
     is in Q(t) and nowhere else.
  3  The matched comparison: the test fermion on the ghost waveform rescaled so that |Q^(Omega)| equals the twin's,
     and on the twin waveform rescaled down to the ghost's.  If n_- follows the rescaling, amplitude is the whole
     story.

Usage (from the repository root, ~5 min):
    python axion_su2/code/axion_probes_fourier.py > axion_su2/results/axion_probes_fourier.txt
Added 2026-10-17.
"""
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import probe_core as pc  # noqa: E402
import axion_probes as ap  # noqa: E402

T0 = time.time()
OMEGA, CRIT, MP = 2.0, 1.0, 50
ap.check_iso_rhs()


def record(sigma, column):
    out = ap.run_iso(OMEGA, CRIT, column, sigma, MP, return_sol=True)
    phi, lam, T, per, eb, nf, eh, lost, sol = out
    tt = np.linspace(0, sol.t[-1], 40001)
    Y = sol.sol(tt)
    return tt, Y, T, nf, lost


def fourier_amp(tt, q, nu):
    dt = tt[1] - tt[0]
    return 2 * abs(np.sum(q * np.exp(-1j * nu * tt)) * dt) / (tt[-1] - tt[0])


def test_fermion(tt, Q, P, zeta0):
    """i zeta' = h(Q(t)) zeta on the recorded Q(t); returns n_- at the end (and the max along the way)."""
    Qi = lambda t: np.interp(t, tt, Q)

    def rhs(t, y):
        z = y[:2] + 1j * y[2:]
        dz = -1j * (pc.fermion_h(Qi(t), P) @ z)
        return np.concatenate([dz.real, dz.imag])
    z0 = zeta0 / np.linalg.norm(zeta0)          # n_- is a ratio: unit norm keeps the tolerances meaningful
    y0 = np.concatenate([z0.real, z0.imag])
    sol = solve_ivp(rhs, (tt[0], tt[-1]), y0, dense_output=True, method='DOP853', rtol=1e-9, atol=1e-12)
    ts = np.linspace(tt[0], tt[-1], 2001)
    Z = sol.sol(ts)
    n = [pc.fermion_adiabatic_number(Qi(t), Z[:2, i] + 1j * Z[2:, i], P) for i, t in enumerate(ts)]
    return n[-1], max(n)


Q_amp = OMEGA / ap.RATIO
m_D = MP * Q_amp
P = pc.Params(fermion=True, mD=m_D)
h0 = pc.fermion_h(Q_amp, P)
ev, U = np.linalg.eigh(h0)
zeta0 = np.sqrt(ap.EPS_P * 1.5 * Q_amp ** 4 / ev[1]) * U[:, 1]
Omega = 2 * np.sqrt(m_D ** 2 + (pc.C_FERMION * Q_amp) ** 2)
print(f"\n  omega = {OMEGA}, omega^2 phi^2 Pi = {CRIT}, Q_amp = {Q_amp}, m_D = {m_D} = M (Part C label said M/2: wrong), level spacing Omega(Q_amp) = {Omega:.4f}")
print(f"\n  1  Fourier amplitudes of the recorded Q(t) over 8 cycles (A cos(nu t) -> A), in units of Q_amp")
print(f"  {'heavy':>5} {'col':>8} | {'<Q>':>7} {'|Q|max':>7} | {'|Q^(Omega)|':>11} {'|Q^(2M)|':>9} {'|Q^(M)|':>8} | coupled n_- at cycle 8")
rec = {}
for sigma, heavy in ((+1, 'ghost'), (-1, 'twin')):
    for column in ('none', 'fermion'):
        tt, Y, T, nf, lost = record(sigma, column)
        Q = Y[0]
        rec[(heavy, column)] = (tt, Q)
        a_om = fourier_amp(tt, Q, Omega) / Q_amp
        a_2 = fourier_amp(tt, Q, 2.0) / Q_amp
        a_1 = fourier_amp(tt, Q, 1.0) / Q_amp
        print(f"  {heavy:>5} {column:>8} | {np.mean(Q) / Q_amp:7.4f} {np.max(np.abs(Q)) / Q_amp:7.3f} | {a_om:11.5f} {a_2:9.5f} {a_1:8.5f} | "
              + (f"{nf[-1]:.3f}" if column == 'fermion' else "-"), flush=True)

print(f"\n  2  The fermion as a test field on each recorded waveform (no back-reaction): n_- at the end / max")
for heavy in ('ghost', 'twin'):
    for column in ('none', 'fermion'):
        tt, Q = rec[(heavy, column)]
        n_end, n_max = test_fermion(tt, Q, P, zeta0)
        print(f"  {heavy:>5} waveform ({column:>7} run): n_- = {n_end:.3f} (max {n_max:.3f})", flush=True)

print(f"\n  3  Matched: the test fermion on the ghost waveform rescaled to the twin's |Q^(Omega)|, and the twin's to the ghost's")
tg, Qg = rec[('ghost', 'fermion')]
tw, Qw = rec[('twin', 'fermion')]
ag, aw = fourier_amp(tg, Qg, Omega), fourier_amp(tw, Qw, Omega)
print(f"  |Q^(Omega)| ghost {ag / Q_amp:.5f}, twin {aw / Q_amp:.5f} Q_amp; ratio twin/ghost = {aw / ag:.3f}")
for label, tt, Q in (("ghost x (twin/ghost)", tg, Qg * (aw / ag)), ("twin x (ghost/twin)", tw, Qw * (ag / aw)),
                     ("ghost, mean removed", tg, Qg - np.mean(Qg)), ("twin, mean removed", tw, Qw - np.mean(Qw))):
    n_end, n_max = test_fermion(tt, Q, P, zeta0)
    print(f"  {label:>22}: n_- = {n_end:.3f} (max {n_max:.3f})", flush=True)
# the pure-tone reference: Q = <Q> + A cos(Omega t) with the twin's and the ghost's A
for label, A, Qm in (("tone at twin's A", aw, np.mean(Qw)), ("tone at ghost's A", ag, np.mean(Qg))):
    n_end, n_max = test_fermion(tg, Qm + A * np.cos(Omega * tg), P, zeta0)
    print(f"  {label:>22}: n_- = {n_end:.3f} (max {n_max:.3f})", flush=True)
print(f"DONE {time.time() - T0:.0f}", flush=True)
