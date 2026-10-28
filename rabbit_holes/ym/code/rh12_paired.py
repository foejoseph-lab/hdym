"""rh12_paired.py -- hole 12's residual question as a paired test: does the nonlinear end state keep any memory of
its block-12 unstable Floquet subspace?

From the per-energy files rh12_direction.py writes (results/rh12_<file>_E<E><suffix>.npz, `rows`): per trajectory,
column 3 is the angle (deg) of the nonlinear end state to the unstable subspace at the right phase, column 5 the angle
of the SAME end state to the subspace at a wrong phase (shifted by a uniform 0.15-0.85 of the period).  Their
difference d = theta_right - theta_wrong is paired, so orbit-phase effects common to both cancel.  Memory of the
subspace means d < 0.

Per energy: n, median d with a 95 % bootstrap interval (10 000 resamples, fixed seed), the fraction with d < 0 and its
two-sided sign-test p-value, the Wilcoxon signed-rank p-value.  Then the same pooled over each band.  A single wrong
phase per trajectory makes d noisier than a multi-phase null would; it does not bias it.

Usage (from the repo root): python rabbit_holes/ym/code/rh12_paired.py [--suffix _n128] > rabbit_holes/ym/results/rh12_paired<suffix>.txt
Seconds.  Added 2026-10-28.
"""
import argparse
import glob
import os
import re
import time

import numpy as np
from scipy.stats import binomtest, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', 'results'))


def stats(d, rng):
    d = d[np.isfinite(d)]
    n = len(d)
    boot = np.median(rng.choice(d, size=(10000, n), replace=True), axis=1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    neg = int((d < 0).sum())
    p_sign = binomtest(neg, n, 0.5).pvalue if n else np.nan
    p_wil = wilcoxon(d).pvalue if n > 0 and np.any(d != 0) else np.nan
    return n, np.median(d), lo, hi, neg / n, p_sign, p_wil


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--suffix', default='')
    a = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(28)
    print(f"# START rh12_paired suffix={a.suffix!r}")
    print("d = angle(end state, unstable subspace at the right phase) - angle(same end state, subspace at a wrong phase), deg;")
    print("memory of the subspace means d < 0.  95 % interval: bootstrap of the median.")
    for band in ('krein', 'param'):
        pat = os.path.join(RES, f'rh12_esc_E_snap_{band}_E*{a.suffix}.npz')
        files = sorted(f for f in glob.glob(pat) if re.search(rf'_E\d\.\d{{4}}{re.escape(a.suffix)}\.npz$', f))
        if not files:
            print(f"\n{band}: no files matching {pat}")
            continue
        print(f"\n== {band} band ==")
        print(f"  {'E':>6} {'n':>4} {'median d':>9} {'95 % interval':>17} {'frac d<0':>9} {'p sign':>8} {'p Wilcoxon':>11}")
        pooled = []
        for f in files:
            E = float(re.search(r'_E(\d\.\d{4})', f).group(1))
            R = np.load(f)['rows']
            if R.shape[1] < 6:
                print(f"  {E:6.4f}: file has no null column (made before 28 Oct 2026)")
                continue
            d = R[:, 3] - R[:, 5]
            pooled.append(d)
            n, med, lo, hi, fneg, ps, pw = stats(d, rng)
            print(f"  {E:6.4f} {n:4d} {med:+9.1f} [{lo:+6.1f}, {hi:+6.1f}] {fneg:9.2f} {ps:8.3f} {pw:11.3g}")
        if pooled:
            n, med, lo, hi, fneg, ps, pw = stats(np.concatenate(pooled), rng)
            print(f"  {'pooled':>6} {n:4d} {med:+9.1f} [{lo:+6.1f}, {hi:+6.1f}] {fneg:9.2f} {ps:8.3f} {pw:11.3g}")
    print(f"# END rh12_paired status=OK wall={time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
