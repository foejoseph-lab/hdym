"""rh12_elapsed.py -- does hole 12's subspace memory depend on how slow the handover is?

From the per-energy files of rh12_direction.py (results/rh12_esc_E_snap_<band>_E<E><suffix>.npz, `rows`): per
trajectory d = rows[:, 3] - rows[:, 5] (angle to the unstable subspace at the right phase minus the angle of the same
end state at a wrong phase, deg; memory means d < 0, as in rh12_paired.py) and the handover time tau = rows[:, 4], the
periods elapsed between the |w| = 0.1 and |w| = 1 crossings.

Reported per band:
- Spearman rho(tau, d) pooled over the band, and within each energy (each energy spans tau ~ 0.3 to 10+ periods,
  so the within-energy test separates "slow handover" from "this energy");
- the energy-stratified test: rho of tau against d after ranking within each energy (rank-transform per energy, then
  pool), so differences between energies cannot drive it;
- median d with a 95 % bootstrap interval in tau bins [0, 1), [1, 3), [3, 10), [10, inf) periods, pooled over the band.
A slow-handover memory predicts rho < 0 within energies and d falling with the bin.

Usage (from the repo root): python rabbit_holes/ym/code/rh12_elapsed.py [--suffix _n128] > rabbit_holes/ym/results/rh12_elapsed<suffix>.txt
Seconds.  Added 2026-10-28.
"""
import argparse
import glob
import os
import re
import time

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', 'results'))
BINS = [0.0, 1.0, 3.0, 10.0, np.inf]


def boot_median(d, rng, nb=10000):
    b = np.median(rng.choice(d, size=(nb, len(d)), replace=True), axis=1)
    return np.median(d), *np.percentile(b, [2.5, 97.5])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--suffix', default='_n128')
    a = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(28)
    print(f"# START rh12_elapsed suffix={a.suffix!r}")
    print("d = angle(end state, unstable subspace, right phase) - angle(same end state, wrong phase), deg (memory: d < 0);")
    print("tau = periods between the |w| = 0.1 and |w| = 1 crossings.  Slow-handover memory predicts rho(tau, d) < 0.")
    for band in ('krein', 'param'):
        files = sorted(f for f in glob.glob(os.path.join(RES, f'rh12_esc_E_snap_{band}_E*{a.suffix}.npz'))
                       if re.search(rf'_E\d\.\d{{4}}{re.escape(a.suffix)}\.npz$', f))
        if not files:
            print(f"\n{band}: no files")
            continue
        print(f"\n== {band} band ==")
        D, TAU, RD, RT = [], [], [], []
        for f in files:
            E = float(re.search(r'_E(\d\.\d{4})', f).group(1))
            R = np.load(f)['rows']
            d, tau = R[:, 3] - R[:, 5], R[:, 4]
            ok = np.isfinite(d) & np.isfinite(tau)
            d, tau = d[ok], tau[ok]
            r = spearmanr(tau, d)
            print(f"  E = {E:.4f}: n = {len(d):4d}, tau median {np.median(tau):5.1f} [p10 {np.percentile(tau, 10):4.1f}, "
                  f"p90 {np.percentile(tau, 90):5.1f}] periods; Spearman rho(tau, d) = {r.statistic:+.3f} (p = {r.pvalue:.3g})")
            D.append(d); TAU.append(tau)
            RD.append((rankdata(d) - 0.5) / len(d)); RT.append((rankdata(tau) - 0.5) / len(tau))
        d, tau = np.concatenate(D), np.concatenate(TAU)
        r = spearmanr(tau, d)
        rs = spearmanr(np.concatenate(RT), np.concatenate(RD))
        print(f"  pooled:   n = {len(d)}, rho(tau, d) = {r.statistic:+.3f} (p = {r.pvalue:.3g});   "
              f"stratified by energy: rho = {rs.statistic:+.3f} (p = {rs.pvalue:.3g})")
        print(f"  {'tau bin (periods)':>18} {'n':>5} {'median d':>9} {'95 % interval':>17}")
        for lo, hi in zip(BINS[:-1], BINS[1:]):
            m = (tau >= lo) & (tau < hi)
            if m.sum() < 5:
                print(f"  {f'[{lo:g}, {hi:g})':>18} {m.sum():5d}   (too few)")
                continue
            med, l, h = boot_median(d[m], rng)
            print(f"  {f'[{lo:g}, {hi:g})':>18} {m.sum():5d} {med:+9.1f} [{l:+6.1f}, {h:+6.1f}]")
    print(f"# END rh12_elapsed status=OK wall={time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
