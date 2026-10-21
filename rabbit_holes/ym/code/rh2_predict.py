"""rh2_predict.py -- rabbit hole 2: test the ladder's out-of-sample predictions.

The fit T_c(n) = 4.057 (n + 0.210) to the centres of the seven block-12 bands of B001 (rungs n = 7..13; odd n
parametric, even n Krein) predicts rung 6 (Krein) at T = 25.20, E = 1.154 and rung 14 (Krein) at T = 57.65,
E = 0.0421, rung 15 (parametric) at E = 0.0321 -- all outside B001's scan [0.05, 1].  This script scans the
half-period map N of block 12 over E in [1.0, 1.5] and [0.025, 0.05] and types every off-circle event.
Usage (from the repo root): python rabbit_holes/ym/code/rh2_predict.py > rabbit_holes/ym/results/rh2_predict.txt
Cost: ~2 s per energy; default 70 energies.  Added 2026-10-21 (rabbit-hole session).
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402
import krein_b005 as kb                      # noqa: E402

A, PHI = 4.057, 0.210
T_HALF = base.T_HALF
STILDE = np.diag(np.array(kb.STILDE['12'] * 4, dtype=float))


def scan(grid):
    Kf = base.load_K()
    rows = []
    for E in grid:
        Ph, M, _, _, T = kb.half_and_full(Kf, E, '12')
        N = STILDE @ Ph
        lam = np.linalg.eig(N)[0]
        g = lam[kb.ghost_four(lam)]
        off = np.abs(np.abs(g) - 1) > 1e-6
        noff = int(off.sum())
        if noff == 0:
            typ = 'stable'
        elif noff == 4:
            typ = 'KREIN quadruplet'
        else:
            m = g[off][np.argmax(np.abs(g[off]))]
            typ = f"FOLD at N = {'+1' if m.real > 0 else '-1'} (M = +1)" if abs(m.imag) < 1e-6 * abs(m) else f"off-circle pair at arg {np.angle(m) / np.pi:+.3f} pi"
        rung = T / A - PHI
        rows.append((E, T, rung, typ, np.abs(g).max()))
        print(f"E={E:.5f} T={T:8.4f} rung T/a - phi = {rung:6.3f} | {typ:28s} max |lambda| {np.abs(g).max():.4f} | angles/pi {np.round(np.sort(np.angle(g) / np.pi), 3)}")
    return rows


def main():
    t0 = time.time()
    print(f"rh2_predict.py: ladder T_c = {A} (n + {PHI}); predictions: rung 6 (Krein) E = {(T_HALF / (A * (6 + PHI))) ** 4 / 2:.4f}, "
          f"rung 14 (Krein) E = {(T_HALF / (A * (14 + PHI))) ** 4 / 2:.4f}, rung 15 (parametric) E = {(T_HALF / (A * (15 + PHI))) ** 4 / 2:.4f}")
    print("\n== E in [1.0, 1.5] (rung 6 expected Krein near E = 1.15) ==")
    r1 = scan(np.linspace(1.0, 1.5, 36))
    print("\n== E in [0.025, 0.05] (rung 14 Krein near 0.042, rung 15 parametric near 0.032) ==")
    r2 = scan(np.logspace(np.log10(0.025), np.log10(0.05), 36))
    for name, rows in (('high', r1), ('low', r2)):
        ev = [r for r in rows if r[3] != 'stable']
        if ev:
            print(f"\n{name}: unstable energies and their rung values:")
            for r in ev:
                print(f"  E={r[0]:.5f} rung={r[2]:.3f} {r[3]} |lambda|={r[4]:.4f}")
        else:
            print(f"\n{name}: no instability found")
    print(f"DONE {time.time() - t0:.0f}")


if __name__ == '__main__':
    main()
