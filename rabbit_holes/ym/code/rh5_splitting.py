"""rh5_splitting.py -- delta(E): the phase splitting of block 12's two ghost pairs, read from B005's log (no integration).

With |psi_A|, |psi_B| the unsigned half-period angles of the two pairs (B005's tracked members k0 and k1), the
full-period angles are phi = 2|psi| mod 2 pi.  Two pairs coincide AS SETS when phi_B = -phi_A (mod 2 pi), so
   Delta(E) = (phi_A + phi_B) mod 2 pi, folded to (-pi, pi],
is their splitting: 0 = degenerate, +-pi = antipodal (one at +1 when the other is at -1).  The decoupled symmetric
model of rh7_selection_rule.py (sector cross-blocks zero) has Delta = 0 exactly (N has lambda and -lambda, angles
summing to pi).  Hole 5's "near-degenerate pairs" and the reviewer's "mirror relation" (Delta = pi) are the two
ends of this one curve.
Usage: python rabbit_holes/ym/code/rh5_splitting.py > rabbit_holes/ym/results/rh5_splitting.txt   (writes .png, .npz)
Added 2026-10-21 (late).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from rh2_ladder import parse_b005, intervals_b001   # noqa: E402

RES = os.path.join(HERE, '..', 'results')


def main():
    E, mod, ang, sgn = parse_b005()
    on = (np.abs(mod - 1) < 1e-4).all(1)
    a, b = np.abs(ang[:, 0]), np.abs(ang[:, 1])
    D = np.mod(2 * a + 2 * b + 1, 2.0) - 1.0                     # units of pi, in (-1, 1]
    print("rh5_splitting.py: Delta(E) = phi_A + phi_B mod 2pi (units of pi), on-circle energies only; 0 = pairs degenerate as sets, +-1 = antipodal")
    print(" E       |psiA| |psiB|  Delta/pi")
    for i in range(len(E)):
        if on[i] and i % 8 == 0:
            print(f"{E[i]:.4f}  {a[i]:.3f}  {b[i]:.3f}  {D[i]:+.3f}")
    print("\nsummary by range (median |Delta|/pi over on-circle energies):")
    for lo, hi in [(0.05, 0.1), (0.1, 0.2), (0.2, 0.35), (0.35, 0.5), (0.5, 0.75), (0.75, 1.0)]:
        m = on & (E >= lo) & (E < hi)
        if m.any():
            print(f"  E in [{lo:.2f}, {hi:.2f}): median |Delta| = {np.median(np.abs(D[m])):.3f} pi, range [{np.abs(D[m]).min():.3f}, {np.abs(D[m]).max():.3f}], {m.sum()} energies")
    np.savez(os.path.join(RES, 'rh5_splitting.npz'), E=E, on=on, Delta=D, psiA=a, psiB=b)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 4))
        for blk, lo, hi, ty, mu in intervals_b001():
            if blk == '12':
                ax.axvspan(lo, hi, color=('tab:red' if ty == 'CPLX' else 'tab:blue'), alpha=0.15, lw=0)
        ax.plot(E[on], np.abs(D[on]), 'k.', ms=3)
        ax.axhline(0, color='grey', lw=0.5)
        ax.axhline(1, color='grey', lw=0.5)
        ax.set_xscale('log')
        ax.set_xlabel('E')
        ax.set_ylabel('|Delta| / pi  (0: degenerate pairs, 1: antipodal)')
        ax.set_title('rabbit hole 5 / delta(E): splitting of the two ghost pairs of block 12 (B005); bands: red Krein, blue parametric', fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, 'rh5_splitting.png'), dpi=130)
        print("figure: results/rh5_splitting.png")
    except Exception as e:                                        # noqa: BLE001
        print(f"(no figure: {e})")


if __name__ == '__main__':
    main()
