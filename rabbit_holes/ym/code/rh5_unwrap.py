"""rh5_unwrap.py -- signature-oriented unwrapping of both ghost pairs' half-period angles through the bands, from B005's
log (the reviewer's recipe): within an on-circle stretch the pair's angle is that of its Krein-positive member; across
an off-circle stretch the same member's logged argument is continued (off-circle arguments are continuous in E); a
change of the positive member's identity at a fold re-entry is recorded and attributed no winding.  The slope
d psi / d ln E in windows gives an effective frequency estimate omega ~ slope / (-T / 8 pi) (branch-ambiguous:
a Floquet frequency is defined mod Omega; see rh5_branch.py for the Fourier test of the branch).
Usage: python rabbit_holes/ym/code/rh5_unwrap.py > rabbit_holes/ym/results/rh5_unwrap.txt    (writes rh5_unwrapped.npz)
Added 2026-10-21 (late).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from rh2_ladder import parse_b005   # noqa: E402

RES = os.path.join(HERE, '..', 'results')


def unwrap(a):
    out = np.zeros_like(a)
    out[0] = a[0]
    for i in range(1, len(a)):
        d = a[i] - a[i - 1]
        d = (d + 1) % 2 - 1
        out[i] = out[i - 1] + d
    return out


def main():
    E, mod, ang, sgn = parse_b005()
    T = 31.053330268585066 * (2 * E) ** -0.25
    U = np.array([unwrap(ang[:, k]) for k in range(4)])
    on = np.abs(mod - 1) < 1e-4
    res = {}
    for name, (i, j) in {'A': (0, 3), 'B': (1, 2)}.items():
        plus = np.full(len(E), -1)
        cur, events = None, []
        for n in range(len(E)):
            if on[n, i] and on[n, j]:
                p = i if sgn[n, i] > 0 else j
                if cur is not None and p != cur:
                    events.append(E[n])
                cur = p
            plus[n] = cur if cur is not None else i
        psi = np.array([U[plus[n], n] for n in range(len(E))])
        out = psi.copy()
        for n in range(1, len(E)):
            if plus[n] != plus[n - 1]:
                out[n:] = psi[n:] - psi[n] + out[n - 1] + (U[plus[n], n] - U[plus[n], n - 1])
        res[name] = out
        print(f"pair {name}: positive-member identity switches at E = {[round(float(e), 4) for e in events]}; unwrapped psi+/pi {out[0]:+.3f} -> {out[-1]:+.3f} (winding {out[-1] - out[0]:+.3f} pi; the omega = 1 line winds {(T[-1] - T[0]) / 2 / np.pi:+.3f} pi)")
    print("\n window          dpsiA+/dlnE  dpsiB+/dlnE   omega_A~   omega_B~   (slope / (-T/8pi) at the window's mean T; branch-ambiguous)")
    for lo, hi in [(0.05, 0.075), (0.08, 0.15), (0.165, 0.22), (0.27, 0.34), (0.5, 0.56), (0.73, 0.79), (0.8, 1.0)]:
        m = (E >= lo) & (E <= hi)
        sa = np.polyfit(np.log(E[m]), res['A'][m], 1)[0]
        sb = np.polyfit(np.log(E[m]), res['B'][m], 1)[0]
        k = -np.mean(T[m]) / (8 * np.pi)
        print(f" [{lo:.3f},{hi:.3f}]    {sa:+.3f}       {sb:+.3f}       {sa / k:.2f}       {sb / k:.2f}")
    np.savez(os.path.join(RES, 'rh5_unwrapped.npz'), E=E, T=T, psiA=res['A'], psiB=res['B'])


if __name__ == '__main__':
    main()
