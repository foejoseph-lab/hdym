"""axion_probes_spread.py -- the spread of the L1 loss time (axion_probes_long.py): the 'both' and 'boson+F' rows
integrate the same ODE (the m_D = 0 fermion's current is constant to ten digits) and lost the background at 32.1 vs
27.1 T (omega = 2.4) and 47.0 vs 28.5 T (omega = 3.0) -- roundoff amplified by a pumped, ultimately unstable system.
So the loss time is not the number to report (withdrawn 16 Oct, README History); this script measures the spread:
the 'boson+F' row at omega^2 phi^2 Pi = 1.7,
initial data perturbed by 1e-6 relative (five seeds), at two tolerances (rtol 1e-9 and 1e-11), both omega.
Reports the loss time of each run and the range; then the twin and M -> infinity controls with the same perturbations.  ~15 min.
Result (README, "The mechanism"): the pumped loss is a fraction -- 8 of 20 perturbed ghost runs lose the background
within 60 T (omega = 2.4: 5 of 10 at 12-46 T; 3.0: 3 of 10 at 30-48 T), c crossing 1/4 in every lost run; twin 0 of 5,
M -> infinity 0 of 5, c <= 1e-4.
    python axion_su2/code/axion_probes_spread.py > axion_su2/results/axion_probes_spread.txt
Added 2026-10-16.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import axion_probes as ap  # noqa: E402

T0 = time.time()
ap.check_iso_rhs()
print(f"  {'omega':>5} {'rtol':>6} {'seed':>4} | lost at (T)   | c > 1/4 at cycle")
for omega in (2.4, 3.0):
    allt = []
    for rtol in (1e-9, 1e-11):
        ap.RT = dict(method='DOP853', rtol=rtol, atol=1e-16)
        for seed in range(5):
            rng = np.random.default_rng(seed)
            orig_initial = ap.initial

            def pert_initial(om, ph, Qa, P, mp, _rng=rng, _orig=orig_initial):
                s = _orig(om, ph, Qa, P, mp)
                return s * (1 + 1e-6 * _rng.normal(size=s.shape)) if seed > 0 else s
            ap.initial = pert_initial
            out = ap.run_iso(omega, 1.7, 'boson+F', +1, 0, ncyc=60, return_sol=True)
            ap.initial = orig_initial
            phi, lam, T, per, eb, nf, eh, lost, sol = out
            cmax = []
            for k in range(60):
                a, b = k * T, min((k + 1) * T, sol.t[-1])
                cmax.append(0.5 * np.max(sol.sol(np.linspace(a, b, 2001))[4] ** 2) if b > a else np.nan)
            cross = next((k + 1 for k in range(60) if np.isfinite(cmax[k]) and cmax[k] > 0.25), None)
            allt.append(lost if lost is not None else np.nan)
            print(f"  {omega:5.2f} {rtol:6.0e} {seed:4d} | {lost if lost is not None else '> 60':>12} | {cross if cross else '-':>5}", flush=True)
    allt = np.array(allt, dtype=float)
    print(f"  omega = {omega}: loss time range {np.nanmin(allt):.1f} - {np.nanmax(allt):.1f} T over {np.sum(np.isfinite(allt))} of {len(allt)} runs "
          f"(median {np.nanmedian(allt):.1f} T)", flush=True)
# the controls with the same perturbations: the twin with both probes and M -> infinity with both, omega = 2.4
print("  controls, same perturbations, omega = 2.4, rtol 1e-9: twin (both probes) and M -> infinity (both)")
print(f"  {'heavy':>5} {'seed':>4} | lost at (T)   | max c over 60 cycles | |Q|max/Q_amp max")
ap.RT = dict(method='DOP853', rtol=1e-9, atol=1e-16)
for heavy in ('twin', 'inf'):
    for seed in range(5):
        rng = np.random.default_rng(seed)
        orig_initial = ap.initial

        def pert_initial(om, ph, Qa, P, mp, _rng=rng, _orig=orig_initial):
            s = _orig(om, ph, Qa, P, mp)
            return s * (1 + 1e-6 * _rng.normal(size=s.shape)) if seed > 0 else s
        ap.initial = pert_initial
        out = ap.run_iso(2.4, 1.7, 'both', -1 if heavy == 'twin' else +1, 0, ncyc=60, return_sol=True, no_ghost=(heavy == 'inf'))
        ap.initial = orig_initial
        phi, lam, T, per, eb, nf, eh, lost, sol = out
        cm = max(0.5 * np.max(sol.sol(np.linspace(k * T, min((k + 1) * T, sol.t[-1]), 1001))[4] ** 2) for k in range(60) if min((k + 1) * T, sol.t[-1]) > k * T)
        print(f"  {heavy:>5} {seed:4d} | {lost if lost is not None else '> 60':>12} | {cm:8.4f} | {np.nanmax(per):6.2f}", flush=True)
print(f"DONE {time.time() - T0:.0f}")
