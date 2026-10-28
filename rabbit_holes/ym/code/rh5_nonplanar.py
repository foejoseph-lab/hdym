"""rh5_nonplanar.py -- rabbit hole 5: are Krein collisions a planar-orbit phenomenon?  Typed Floquet scan of all three
off-diagonal blocks on the NON-PLANAR periodic orbit x = y = z of Yang-Mills mechanics.

On x = y = z = s(t) the potential is V = 3 s^4 / 2 and the motion is the quartic oscillator s'' = -2 s^3 (the "diagonal"
orbit; hyperbolic in the base flow, which does not matter for the ghost's Floquet problem).  Energy E = 3 s'^2/2 +
3 s^4/2; at E = 1/2 the turning point is s0 = 3^{-1/4}; other energies by the scaling s -> lam s, t -> t / lam,
lam = (2E)^{1/4}.  Period T = 4 int_0^{s0} ds / sqrt(2(E - 3 s^4/2)/3).  z != 0 throughout, so this is the test the list
asks for: on the planar DR orbit blocks 13 and 23 are exactly block-diagonal (no cross-coupling) and only block 12
shows Krein collisions; on x = y = z the three blocks are related by the cyclic permutation symmetry and each has
cross-coupling, so if Krein collisions are "a block-12 thing" they should be absent or present in all three alike.

Output per E: for each block the ghost multipliers (drop-4 rule), typed as in dr_offdiag_floquet_typed.py (REAL+ /
REAL- / CPLX / stable), the Gauss violation of the unstable eigenvectors, and the orbit closure.  Also the half-period
symmetry check: on this orbit s(t + T/2) = -s(t), so K(t + T/2) = S K(t) S for S = the sign flip of the odd-derivative
entries? -- measured rather than assumed: the residual of S~ K(t+T/2) S~ - K(t) for the candidate S~ = diag(+1,+1,-1,-1,
+1,+1,-1,-1) per block (u, u'' even; u', u''' odd under t -> t + T/2 with s -> -s), printed.

Usage (from the repo root): python rabbit_holes/ym/code/rh5_nonplanar.py [--nE 120] [--emin 0.05] [--emax 1.0]
(the committed results/rh5_nonplanar.txt has 300 energies in [0.02, 3], i.e. --nE 300 --emin 0.02 --emax 3)
Cost ~1 s per E.  Added 2026-10-22.
"""
import argparse
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp, quad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402

RES = os.path.join(HERE, '..', 'results')


def orbit111(E):
    lam = (2.0 * E) ** 0.25
    s0 = 3.0 ** -0.25 * lam                       # turning point: 3 s0^4 / 2 = E
    # period: s'' = -2 s^3, energy per component e = s'^2/2 + s^4/2 = E/3
    e = E / 3.0
    Tq = 4 * quad(lambda s: 1.0 / np.sqrt(2 * (e - s ** 4 / 2)), 0, s0, limit=200)[0]
    q0 = np.array([s0, s0, s0])
    p0 = np.zeros(3)
    return q0, p0, Tq


def typed(M):
    g = base.ghost_multipliers(M)
    off = np.abs(np.abs(g) - 1) > 1e-6
    if not off.any():
        return 'stable', g
    if off.sum() == 4 and np.abs(g[off].imag).max() > 1e-6 * np.abs(g[off]).max():
        return 'CPLX', g
    m = g[off][np.argmax(np.abs(g[off]))]
    if abs(m.imag) > 1e-6 * abs(m):
        return 'CPLX', g
    return ('REAL+' if m.real > 0 else 'REAL-') + (' x2' if off.sum() == 4 else ''), g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nE', type=int, default=120)
    ap.add_argument('--emin', type=float, default=0.05)
    ap.add_argument('--emax', type=float, default=1.0)
    a = ap.parse_args()
    Kf = base.load_K()
    import pickle
    import sympy as sp
    d = pickle.load(open(base.PKL, 'rb'))
    Gf = sp.lambdify(list(d['X']) + list(d['P']) + [d['cj']], d['G'], 'numpy')
    print("rh5_nonplanar.py: typed Floquet scan of the three off-diagonal blocks on the x = y = z orbit")
    grid = np.logspace(np.log10(a.emin), np.log10(a.emax), a.nE)
    # half-period symmetry check at one E
    q0, p0, T = orbit111(0.3)
    sol = solve_ivp(lambda t, s: np.concatenate([s[3:], -base.gradV(s[:3])]), (0, T), np.concatenate([q0, p0]), t_eval=np.linspace(0, T, 64, endpoint=False), **base.RT)
    St = np.diag([1, 1, -1, -1, 1, 1, -1, -1.0])
    for b in base.BLOCKS:
        idx = base.block_idx(b)
        Ks = [np.array(Kf(*sol.y[:3, i], *sol.y[3:, i], 1.0), dtype=float)[np.ix_(idx, idx)] for i in range(64)]
        r = max(np.abs(St @ Ks[(i + 32) % 64] @ St - Ks[i]).max() for i in range(64)) / max(np.abs(K).max() for K in Ks)
        cross = max(np.abs(K[np.ix_([6, 7], [0, 1])] - np.diag(np.diag(K[np.ix_([6, 7], [0, 1])]))).max() for K in Ks) / max(np.abs(K).max() for K in Ks)
        print(f"  E=0.3 block {b}: half-period symmetry S~ K(t+T/2) S~ = K(t) with S~ = diag(+,+,-,-,+,+,-,-): residual {r:.1e}; component cross-coupling (off-diagonal of the u'''' rows' u block) relative {cross:.2f}")
    print(f"  orbit: E=0.3, T = {T:.4f}, closure after one period {np.abs(sol.y[:, 0] - solve_ivp(lambda t, s: np.concatenate([s[3:], -base.gradV(s[:3])]), (0, T), np.concatenate([q0, p0]), **base.RT).y[:, -1]).max():.1e}")
    print("\ncolumns: E | T | per block: type |mu|max (Gauss violation of the leading unstable eigenvector)")
    rows = []
    worst_gv = 0.0
    for E in grid:
        q0, p0, T = orbit111(E)
        # monodromies on THIS orbit (base.monodromies is hard-wired to the DR orbit): integrate here
        def rhs(t, s):
            q, p = s[:3], s[3:6]
            K = np.array(Kf(*q, *p, 1.0), dtype=float)
            out = [p, -base.gradV(q)]
            off = 6
            for b in base.BLOCKS:
                idx = base.block_idx(b)
                out.append((K[np.ix_(idx, idx)] @ s[off:off + 64].reshape(8, 8)).ravel())
                off += 64
            return np.concatenate(out)
        s0 = np.concatenate([q0, p0] + [np.eye(8).ravel()] * 3)
        sT = solve_ivp(rhs, (0.0, T), s0, **base.RT).y[:, -1]
        clo = np.abs(sT[:6] - s0[:6]).max()
        line = f"E={E:.5f} T={T:8.4f} clo={clo:.0e} |"
        rec = [E, T]
        for k, b in enumerate(base.BLOCKS):
            M = sT[6 + 64 * k:6 + 64 * (k + 1)].reshape(8, 8)
            ty, g = typed(M)
            gv = np.nan
            if ty != 'stable':
                mu, V = np.linalg.eig(M)
                j = int(np.argmax(np.abs(mu)))
                G = np.array(Gf(*q0, *p0, 1.0), dtype=float)[:, base.block_idx(b)]
                G = G[np.abs(G).sum(1) > 0]
                gv = np.linalg.norm(G @ V[:, j]) / np.linalg.norm(V[:, j]) if len(G) else 0.0
                worst_gv = max(worst_gv, gv)
            line += f" {b}: {ty:9s} {np.abs(g).max():7.4f} ({gv:.0e}) |"
            rec += [np.abs(g).max(), ty]
        rows.append(rec)
        print(line)
    print(f"\nCHECK unstable eigenvectors satisfy the linearised Gauss law: worst {worst_gv:.1e} (tol 1e-6) {'PASS' if worst_gv < 1e-6 else 'FAIL'}")
    # summary: bands per block
    for k, b in enumerate(base.BLOCKS):
        types = [r[3 + 2 * k] for r in rows]
        Es = [r[0] for r in rows]
        bands, i = [], 0
        while i < len(types):
            if types[i] != 'stable':
                j = i
                while j + 1 < len(types) and types[j + 1] == types[i]:
                    j += 1
                bands.append((Es[i], Es[j], types[i]))
                i = j + 1
            else:
                i += 1
        print(f"block {b}: {len(bands)} bands: " + "; ".join(f"[{lo:.4f},{hi:.4f}] {ty}" for lo, hi, ty in bands))
    np.savez(os.path.join(RES, 'rh5_nonplanar.npz'), rows=np.array([[r[0], r[1], r[2], r[4], r[6]] for r in rows]), types=np.array([[r[3], r[5], r[7]] for r in rows]))


if __name__ == '__main__':
    main()
