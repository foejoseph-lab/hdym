"""rh11_weakcell.py -- rabbit hole 11: B002's one in-model Krein cell (D = 3.8946, J = 0.0315, psi = 0, c_j = 1: |mu| = 1.002
at arg 0.065 pi, blocks 12 and 13) -- artefact of the frozen background or real?

The suspicion (response to the 19 Oct review): the frozen-depth background is not an exact solution, so the gauge
Jordan pair at 1 splits by ~sqrt(eps) = 0.15 and a "collision" at arg 0.065 pi ~ 0.2 rad is at the scale of that split.
Two tests, neither needing a periodic orbit:

  A. eps scaling in the frozen model.  At fixed D, J -> J/8 ... 2J (eps proportional to J).  If |mu| - 1 of the cell
     tracks the Jordan split (~ sqrt(eps)) or eps itself, it is the background's; a Krein collision of the physical
     pairs would instead switch on at a threshold in J (a collision needs the two pairs to meet) and grow with the
     coupling above it.  Printed per J: eps, the Jordan split (min |mu - 1| of the gauge pair), the four ghost-sector
     eigenvalues (drop-4 rule, which is meaningful only as a label here), type, |mu|max.
  B. the exact background.  Same initial data (x = D, px = 0, y = z = A, py = pz = 0) integrated with the real
     Yang-Mills equations for one transverse period T = 2 pi / D (depth drifts by eps, so this is a transfer matrix,
     not a monodromy; det = 1 and the flow is Hamiltonian, but eigenvalues are not invariants).  What the exact
     background restores is the constraint: the gauge zero mode is exact and the linearised Gauss law is conserved,
     so any off-circle quadruplet is physical growth over that period, not constraint-sector leakage.  Reported: the
     8 eigenvalues of the transfer matrix per block, their Gauss residuals, and -- the honest number -- the growth of
     the Gauss-clean subspace over the period, max singular value of the transfer matrix restricted to ker G(0)
     and projected onto ker G(T), compared with the same on the frozen background.

Usage (from the repo root): python rabbit_holes/ym/code/rh11_weakcell.py
Cost: seconds.  Added 2026-10-22.
"""
import os
import pickle
import sys

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402
import channel_offdiag_monodromy as CH       # noqa: E402

D0, J0, PSI, CJ = 3.8946, 0.0315, 0.0, 1.0


def typed(g):
    off = np.abs(np.abs(g) - 1) > 1e-6
    if not off.any():
        return 'stable'
    if np.abs(g[off].imag).max() > 1e-6 * np.abs(g[off]).max():
        return 'CPLX'
    return 'REAL+' if g[off].real.max() > 0 else 'REAL-'


def main():
    Kf = base.load_K()
    d = pickle.load(open(base.PKL, 'rb'))
    Gf = sp.lambdify(list(d['X']) + list(d['P']) + [d['cj']], d['G'], 'numpy')
    CH._Kf = Kf
    print(f"rh11_weakcell.py: B002's weak Krein cell D={D0}, J={J0}, psi={PSI}, c_j={CJ}")
    print("\n== A. eps scaling in the frozen model (block 12) ==")
    print(" J        eps      Jordan split |mu-1|min   ghost-sector eigenvalues (drop-4)                      type   |mu|max-1")
    for J in J0 * np.array([1 / 8, 1 / 4, 1 / 2, 1, 2, 4]):
        T, Ms, _ = CH.monodromies(D0, J, PSI, CJ)
        M = Ms['12']
        mu = np.linalg.eigvals(M)
        split = np.abs(mu - 1).min()
        g = base.ghost_multipliers(M)
        eps = 4 * np.pi ** 2 * J / D0 ** 3
        print(f" {J:.5f}  {eps:.2e}  {split:.3e}             {np.round(g, 4)}  {typed(g):6s} {np.abs(g).max() - 1:+.2e}")

    print("\n== B. exact background over one transverse period (block 12 and 13) ==")
    qp, T = CH.background(D0, J0, PSI)
    q0, p0 = qp(0.0)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, CJ), dtype=float)
        out = [p, -base.gradV(q)]
        off = 6
        for b in base.BLOCKS:
            idx = base.block_idx(b)
            out.append((K[np.ix_(idx, idx)] @ s[off:off + 64].reshape(8, 8)).ravel())
            off += 64
        return np.concatenate(out)
    sT = solve_ivp(rhs, (0.0, T), np.concatenate([q0, p0] + [np.eye(8).ravel()] * 3), **base.RT).y[:, -1]
    qT, pT = sT[:3], sT[3:6]
    print(f"  background after one period: x {q0[0]:.4f} -> {qT[0]:.4f} (drift {(qT[0] - q0[0]) / q0[0]:+.2e}, eps model {4 * np.pi ** 2 * J0 / D0 ** 3:.2e}); y {q0[1]:.4f} -> {qT[1]:.4f}")
    for k, b in enumerate(base.BLOCKS):
        idx = base.block_idx(b)
        M = sT[6 + 64 * k:6 + 64 * (k + 1)].reshape(8, 8)
        G0 = np.array(Gf(*q0, *p0, CJ), dtype=float)[:, idx]
        GT = np.array(Gf(*qT, *pT, CJ), dtype=float)[:, idx]
        G0 = G0[np.abs(G0).sum(1) > 0]
        GT = GT[np.abs(GT).sum(1) > 0]
        mu, V = np.linalg.eig(M)
        order = np.argsort(-np.abs(mu))
        gv = [np.linalg.norm(G0 @ V[:, j]) / np.linalg.norm(V[:, j]) if len(G0) else 0.0 for j in order]
        # Gauss-clean subspace growth: restrict M to ker G0, project onto ker GT
        def kerbasis(G):
            if len(G) == 0:
                return np.eye(8)
            _, _, Vt = np.linalg.svd(G)
            return Vt[len(G):].T
        B0, BT = kerbasis(G0), kerbasis(GT)
        Mr = BT.T @ M @ B0
        leak = np.linalg.norm(GT @ M @ B0) / np.linalg.norm(M @ B0) if len(GT) else 0.0
        sv = np.linalg.svd(Mr, compute_uv=False)
        # frozen comparison: same restricted growth on the frozen background
        _, Mfz, _ = CH.monodromies(D0, J0, PSI, CJ)
        Mf = Mfz[b]
        Gf0 = np.array(Gf(*q0, *p0, CJ), dtype=float)[:, idx]
        Gf0 = Gf0[np.abs(Gf0).sum(1) > 0]
        Bf = kerbasis(Gf0)
        svf = np.linalg.svd(Bf.T @ Mf @ Bf, compute_uv=False)
        leakf = np.linalg.norm(Gf0 @ Mf @ Bf) / np.linalg.norm(Mf @ Bf) if len(Gf0) else 0.0
        print(f"  block {b}: exact-background transfer matrix, det-1 = {np.linalg.det(M) - 1:.1e}; eigenvalues by |.|: {np.round(mu[order], 4)}")
        print(f"           Gauss residual of each eigenvector at t=0: {np.round(gv, 3)}")
        print(f"           Gauss-conservation of the flow (leak of ker G(0) off ker G(T)): {leak:.1e}; frozen background: {leakf:.1e}")
        print(f"           singular values of the Gauss-clean restriction over one period: exact {np.round(sv, 4)};  frozen {np.round(svf, 4)}")


if __name__ == '__main__':
    main()
