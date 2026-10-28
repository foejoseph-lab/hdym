"""rh6_frozen_K.py -- the frozen block-12 generator K(t) (8 x 8, c_j = 1) along the DR orbit: its eigenvalues at the
turning points (t = 0, T/2: x extremal, y = 0) and the bounces (t = T/4, 3T/4), plus the eighth-periods between.

Why (28 Oct 2026): README "Smaller things" put the real pair +-1.32 of the frozen K at the turning point, RH6 put the
same numbers at the bounce, and neither had a log.  This prints both, so the place is settled by a committed log.
Frozen eigenvalues are the eigenvalues of K at fixed t (not Floquet multipliers); each block also carries a gauge zero
mode and a Gauss charge, hence the eigenvalue pair at 0.

Usage (from the repo root): python rabbit_holes/ym/code/rh6_frozen_K.py [--E 0.3] > rabbit_holes/ym/results/rh6_frozen_K.txt
Seconds.  Added 2026-10-28.
"""
import argparse
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402

FRACS = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]
NAMES = {0.0: 'turning point', 0.25: 'bounce', 0.5: 'turning point', 0.75: 'bounce'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--E', type=float, nargs='+', default=[0.3])
    ap.add_argument('--cj', type=float, default=1.0)
    a = ap.parse_args()
    t0 = time.time()
    import scipy
    print(f"# START rh6_frozen_K E={a.E} cj={a.cj}  python {sys.version.split()[0]}  numpy {np.__version__}  scipy {scipy.__version__}")
    Kf = base.load_K()
    idx = base.block_idx('12')
    for E in a.E:
        q0, p0, T = base.orbit(E)
        ts = np.array(FRACS) * T
        sol = solve_ivp(lambda t, s: np.concatenate([s[3:], -base.gradV(s[:3])]), (0.0, T), np.concatenate([q0, p0]),
                        t_eval=ts, **base.RT)
        print(f"\nE = {E:g}, c_j = {a.cj:g}, T = {T:.3f}; eigenvalues of the frozen block-12 K, sorted by |.|, as Re+Im i")
        for k, t in enumerate(ts):
            q, p = sol.y[:3, k], sol.y[3:, k]
            K = np.array(Kf(*q, *p, a.cj), dtype=float)[np.ix_(idx, idx)]
            ev = np.linalg.eigvals(K)
            ev = ev[np.argsort(-np.abs(ev))]
            real = [e.real for e in ev if abs(e.imag) < 1e-9 and abs(e.real) > 1e-3]
            label = NAMES.get(FRACS[k], '')
            print(f"  t/T = {FRACS[k]:.3f}  x = {q[0]:+.3f}  y = {q[1]:+.3f}  {label:13s} "
                  + ' '.join(f"{e.real:+.3f}{e.imag:+.3f}i" for e in ev)
                  + (f"   real pair: {', '.join(f'{r:+.3f}' for r in real)}" if real else "   no real pair"))
    print(f"# END rh6_frozen_K status=OK wall={time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
