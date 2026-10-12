"""hdym_tailA.py -- the u^-6 tail constant A from a saved hdym_spectrum.py run.

Reads <dir>/spectrum.npz written by hdym_spectrum.py (om = u, S_diag = s(u) of the mean
diagonal Hessian entry) and prints, without assuming any value of A:
  - the local log-log slopes over the paper's u-bands (5-8, 8-12, 12-20, 20-30, 30-45),
  - u^6 s(u) per band and over u = 8-45: mean, min/max of the band means, and the spread,
  - A as the mean of u^6 s(u) over u = 8-45, with the band-to-band spread as its +-.
The paper's A = 21.6 +- 0.3 and slopes are printed only at the end, for comparison.

Usage (repository root):
    python ym/hdym_tailA.py ym/data/results_spectrum > ym/data/tailA_output.txt
Added 2026-10-12.
"""
import os
import sys
import numpy as np

BANDS = [(5, 8), (8, 12), (12, 20), (20, 30), (30, 45)]


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "results_spectrum"
    f = np.load(os.path.join(d, "spectrum.npz"))
    u, s = f["om"], f["S_diag"]
    print(f"hdym_tailA.py  file {os.path.join(d, 'spectrum.npz')}  taper {f['taper']}  "
          f"{len(f['depth'])} segments")
    ok = (u > 1) & (s > 0)
    lu, ls = np.log(u[ok]), np.log(s[ok])
    y = u ** 6 * s
    print(f"\n{'band':>8} {'slope':>7} {'mean u^6 s':>11} {'min':>9} {'max':>9}  points")
    means = []
    for u0, u1 in BANDS:
        m = (u >= u0) & (u <= u1) & (s > 0)
        mm = (lu >= np.log(u0)) & (lu <= np.log(u1))
        slope = np.polyfit(lu[mm], ls[mm], 1)[0]
        ym = y[m].mean()
        if u0 >= 8:
            means.append(ym)
        print(f"{u0:>3}-{u1:<4} {slope:7.2f} {ym:11.3f} {y[m].min():9.3f} {y[m].max():9.3f}  {m.sum()}")
    m = (u >= 8) & (u <= 45) & (s > 0)
    A = y[m].mean()
    spread = 0.5 * (max(means) - min(means))
    print(f"\nA = mean of u^6 s(u) over u = 8-45: {A:.2f}")
    print(f"    band means over 8-45: {', '.join(f'{v:.2f}' for v in means)};  half-range {spread:.2f} "
          f"({100 * spread / A:.1f} %)")
    print(f"RESULT: A = {A:.2f} +- {spread:.2f}")
    print("\ncomparison (not used above): paper/main.tex eq:A gives A = 21.6 +- 0.3, slopes "
          "-6.14, -6.05, -6.00, -5.98, -5.89; NOTE_2d.md predicts A = (7/2) A_res = 21.5")


if __name__ == "__main__":
    main()
