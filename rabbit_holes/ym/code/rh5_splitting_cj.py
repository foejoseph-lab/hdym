"""rh5_splitting_cj.py -- Delta(E) of block 12's two ghost pairs at a chosen c_j (the reviewer's test of the width-by-
splitting reading of the c_j family).  Same quantity as rh5_splitting.py (Delta = phi_A + phi_B mod 2pi from the
half-period map N = S~ Phi_T/2), computed here rather than parsed, with the pairs matched by continuity in E.
Usage: python rabbit_holes/ym/code/rh5_splitting_cj.py --cj 0.25 [--nE 150] > rabbit_holes/ym/results/rh5_splitting_cj025.txt
Cost ~1.7 s per E.  Added 2026-10-21 (late).
"""
import argparse
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402
import krein_b005 as kb                      # noqa: E402

RES = os.path.join(HERE, '..', 'results')
IDX = base.block_idx('12')
STILDE = np.diag(np.array(kb.STILDE['12'] * 4, dtype=float))


def half_map(Kf, E, cj):
    q0, p0, T = base.orbit(E)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, cj), dtype=float)[np.ix_(IDX, IDX)]
        return np.concatenate([p, -base.gradV(q), (K @ s[6:].reshape(8, 8)).ravel()])
    sh = solve_ivp(rhs, (0.0, T / 2), np.concatenate([q0, p0, np.eye(8).ravel()]), **base.RT).y[:, -1]
    return T, STILDE @ sh[6:].reshape(8, 8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cj', type=float, default=0.25)
    ap.add_argument('--nE', type=int, default=150)
    ap.add_argument('--emin', type=float, default=0.05)
    ap.add_argument('--emax', type=float, default=1.0)
    a = ap.parse_args()
    Kf = base.load_K()
    grid = np.logspace(np.log10(a.emin), np.log10(a.emax), a.nE)
    print(f"rh5_splitting_cj.py: c_j = {a.cj}; Delta = phi_A + phi_B mod 2pi (units of pi, folded to (-1,1]); pairs = the two conjugate pairs among the four ghost eigenvalues of N")
    print(" E       |psiA|  |psiB|  |lamA| |lamB|  Delta/pi   (both on circle?)")
    out = []
    for E in grid:
        T, N = half_map(Kf, E, a.cj)
        lam = np.linalg.eigvals(N)
        g = lam[kb.ghost_four(lam)]
        # pair up conjugates
        used, pairs = set(), []
        for i in range(4):
            if i in used:
                continue
            j = min((abs(g[k] - np.conj(g[i])), k) for k in range(4) if k != i and k not in used)[1]
            used |= {i, j}
            pairs.append(g[i] if g[i].imag >= 0 else g[j])
        pairs.sort(key=lambda z: abs(np.angle(z)))
        psi = [abs(np.angle(z)) / np.pi for z in pairs]
        mods = [abs(z) for z in pairs]
        on = all(abs(m - 1) < 1e-6 for m in mods)
        D = np.mod(2 * psi[0] + 2 * psi[1] + 1, 2.0) - 1.0
        out.append((E, psi[0], psi[1], mods[0], mods[1], D, on))
        print(f"{E:.4f}  {psi[0]:.3f}   {psi[1]:.3f}   {mods[0]:.3f}  {mods[1]:.3f}  {D:+.3f}   {'on' if on else 'OFF'}")
    R = np.array([(r[0], r[1], r[2], r[3], r[4], r[5], float(r[6])) for r in out])
    np.savez(os.path.join(RES, f'rh5_splitting_cj{a.cj:g}.npz'), rows=R)
    on = R[:, 6] > 0
    for lo, hi in [(0.05, 0.2), (0.2, 0.41), (0.41, 0.7), (0.7, 1.0)]:
        m = on & (R[:, 0] >= lo) & (R[:, 0] < hi)
        if m.any():
            print(f"summary E in [{lo}, {hi}): median |Delta| = {np.median(np.abs(R[m, 5])):.3f} pi over {m.sum()} on-circle energies; off-circle energies in range: {int(((~on) & (R[:, 0] >= lo) & (R[:, 0] < hi)).sum())}")


if __name__ == '__main__':
    main()
