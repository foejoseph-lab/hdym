"""dr_offdiag_krein.py -- B005: Krein signatures of the off-diagonal ghost pairs on the DR orbit,
computed on the half-period map N (M = N^2, see dr_offdiag_halfperiod.py).

What is done, per block and per E on the B001 grid (400 log-spaced points in [0.05, 1]):
  * Phi_{T/2} by integration; N = St Phi_{T/2}; M = N^2 (checked against Phi_T to 1e-11 in
    dr_offdiag_halfperiod.py, not repeated here).
  * Omega: the Ostrogradsky symplectic form of the block, pulled back to the state coordinates
    (u, u', u'', u''') at t = 0 from the exported quadratic Lagrangian L2 (ym/rung1_system.pkl).
    Check printed: |N^T Omega N - Omega| / |Omega|.
  * Trace formula vs drop rule for the Broucke invariants of M:
    A = tr M - 4,  B = (A^2 - tr M^2 + 4) / 2  (the gauge sector contributes exactly 4 to each
    trace on an exact solution) against A, B from the four multipliers not nearest 1.
  * Deflation: E1 = generalised eigenspace of M at 1 (4-dim, the gauge sector), W = its
    Omega-orthogonal complement (4-dim, N-invariant).  N_W = N restricted to W, Omega_W likewise.
  * Krein sign of each unit-circle multiplier of N_W:  kappa = sign( i v^* Omega_W v ).
    For each of the two pairs the Krein-signed angle psi+ (argument of the member with kappa > 0)
    is reported; when a pair is off the circle its modulus is reported instead and the sign is
    undefined.

Output: results/dr_offdiag_krein.txt (log) and results/dr_offdiag_krein.npy, rows
  E, T, sympl_err, A_tr, B_tr, A_drop, B_drop, then per block (12, 13, 23):
  |mu1|, psi1+, |mu2|, psi2+, kappa_defined(0/1), gauss_resid_of_W
where |mu| is for N (so M's modulus is its square) and psi is N's angle in (-pi, pi].

Usage (from islands/):  python code/dr_offdiag_krein.py [--nE 400] [--procs 4]
Added 2026-10-19 (review).
"""
import argparse
import os
import pickle
import sys
import time
from multiprocessing import Pool

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dr_offdiag_floquet as base  # noqa: E402
from dr_offdiag_halfperiod import S_BG, sign_matrix  # noqa: E402

PKL = base.PKL


# ---------------------------------------------------------------- Omega from L2 (Ostrogradsky)
def build_omega():
    """Returns omega_f(block_name)(q, p) -> 8x8 Omega in the block's state coordinates
    (u0_a, u0_b, u1_a, u1_b, u2_a, u2_b, u3_a, u3_b), evaluated at the background (q, p)."""
    d = pickle.load(open(PKL, 'rb'))
    L2, state, X, P, cj = d['L2'], d['state'], d['X'], d['P'], d['cj']
    L2 = L2.subs(cj, 1)
    U = [[sp.Symbol(f'u{k}_{j}') for j in range(6)] for k in range(5)]
    Vs = (X[0] ** 2 * X[1] ** 2 + X[1] ** 2 * X[2] ** 2 + X[2] ** 2 * X[0] ** 2) / 2
    qdd = [-Vs.diff(v) for v in X]

    def ddt(expr):
        """time derivative along the ghost-free background, u_k -> u_{k+1}."""
        out = sum(expr.diff(X[i]) * P[i] + expr.diff(P[i]) * qdd[i] for i in range(3))
        for k in range(4):
            for j in range(6):
                out += expr.diff(U[k][j]) * U[k + 1][j]
        return sp.expand(out)

    # canonical momenta: p1_j = dL/du1_j - d/dt dL/du2_j,  p2_j = dL/du2_j   (L2 depends on u0,u1,u2)
    assert all(L2.diff(U[3][j]) == 0 and L2.diff(U[4][j]) == 0 for j in range(6))
    p1 = [sp.expand(L2.diff(U[1][j]) - ddt(L2.diff(U[2][j]))) for j in range(6)]
    p2 = [sp.expand(L2.diff(U[2][j])) for j in range(6)]
    funcs = {}
    for name, comps in base.BLOCKS.items():
        idx = base.block_idx(name)
        xs = [state[i] for i in idx]                       # the block's 8 state symbols
        canon = ([U[0][j] for j in comps] + [U[1][j] for j in comps]
                 + [p1[j] for j in comps] + [p2[j] for j in comps])
        C = sp.Matrix(8, 8, lambda r, c: sp.expand(canon[r]).coeff(xs[c]))
        # canonical form on (q1, q2, p1, p2): omega = dq1^dp1 + dq2^dp2
        J = sp.zeros(8, 8)
        for r in range(4):
            J[r, r + 4] = 1
            J[r + 4, r] = -1
        Om = sp.simplify(C.T * J * C)
        # the block's Omega must not involve other components' symbols
        assert not (Om.free_symbols - set(X) - set(P)), Om.free_symbols
        funcs[name] = sp.lambdify(list(X) + list(P), Om, 'numpy')
    return funcs


# ---------------------------------------------------------------- per-E work
_Kf = _Of = _St = None


def _init():
    global _Kf, _Of, _St
    _Kf = base.load_K()
    _Of = build_omega()
    _St = {b: np.diag(sign_matrix(_Kf, b)[1]) for b in base.BLOCKS}


def half_map(E):
    q0, p0, T = base.orbit(E)
    nb = len(base.BLOCKS)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(_Kf(*q, *p, 1.0), dtype=float)
        out = [p, -base.gradV(q)]
        off = 6
        for b in base.BLOCKS:
            idx = base.block_idx(b)
            Mb = s[off:off + 64].reshape(8, 8)
            out.append((K[np.ix_(idx, idx)] @ Mb).ravel())
            off += 64
        return np.concatenate(out)

    s0 = np.concatenate([q0, p0] + [np.eye(8).ravel()] * nb)
    sh = solve_ivp(rhs, (0.0, T / 2), s0, **base.RT).y[:, -1]
    Ph = {}
    off = 6
    for b in base.BLOCKS:
        Ph[b] = sh[off:off + 64].reshape(8, 8)
        off += 64
    return T, q0, p0, Ph


def deflate(M, N, Om):
    """W = Omega-orthogonal complement of the generalised eigenspace of M at 1."""
    I = np.eye(8)
    A2 = (M - I) @ (M - I)
    _, s, vt = np.linalg.svd(A2)
    E1 = vt[4:].T                                  # 4 smallest singular directions = ker (M-1)^2
    # W = { x : E1^T Omega x = 0 }
    _, _, vt2 = np.linalg.svd(E1.T @ Om)
    W = vt2[4:].T                                  # 8x4 basis of W
    # restrict: N W = W N_W  ->  N_W = (W^T W)^{-1} W^T N W ; W is orthonormal
    N_W = W.T @ N @ W
    inv_err = np.linalg.norm(N @ W - W @ N_W) / np.linalg.norm(N @ W)   # invariance of W
    Om_W = W.T @ Om @ W
    return N_W, Om_W, inv_err, s[4:]


def krein_pairs(N_W, Om_W):
    """Returns two pairs: (modulus of N multiplier, Krein-signed angle psi+ or nan, sign-defined)."""
    mu, V = np.linalg.eig(N_W)
    kappa = np.array([np.imag(np.vdot(v, Om_W @ v)) for v in V.T])   # i v^* Om v is real: use Im(v^* Om v)
    on = np.abs(np.abs(mu) - 1) < 1e-6
    out = []
    used = np.zeros(4, bool)
    for i in np.argsort(-np.abs(mu)):
        if used[i]:
            continue
        # partner: 1/mu (on the circle this is conj(mu); for a real pair the reciprocal)
        j = np.argmin(np.where(used | (np.arange(4) == i), np.inf, np.abs(mu - 1 / mu[i])))
        used[i] = used[j] = True
        if on[i] and on[j] and abs(mu[i].imag) > 1e-9:
            pos = i if kappa[i] > 0 else j
            out.append((abs(mu[i]), np.angle(mu[pos]), 1.0, abs(kappa[i])))
        else:
            out.append((max(abs(mu[i]), abs(mu[j])), np.angle(mu[i]), 0.0, 0.0))
    return out


def one(E):
    T, q0, p0, Ph = half_map(E)
    row = [E, T]
    log = []
    for b in base.BLOCKS:
        N = _St[b] @ Ph[b]
        M = N @ N
        Om = np.array(_Of[b](*q0, *p0), dtype=float)
        se = np.linalg.norm(N.T @ Om @ N - Om) / np.linalg.norm(Om)
        A_tr = np.trace(M) - 4
        B_tr = (A_tr ** 2 - np.trace(M @ M) + 4) / 2
        g = base.ghost_multipliers(M)
        A_d = np.sum(g).real
        B_d = ((np.sum(g) ** 2 - np.sum(g ** 2)) / 2).real
        N_W, Om_W, inv_err, sv = deflate(M, N, Om)
        pairs = krein_pairs(N_W, Om_W)
        row += [se, A_tr, B_tr, A_d, B_d, inv_err]
        for (m, psi, defined, kap) in pairs:
            row += [m, psi, defined]
        log.append(f"blk{b} sympl {se:.1e} inv {inv_err:.1e} |A_tr-A_drop| {abs(A_tr - A_d):.1e} "
                   f"|B_tr-B_drop| {abs(B_tr - B_d):.1e} pairs "
                   + " ".join(f"({m:.5f},{psi / np.pi:+.4f}pi,{'K' if d_ else '-'})" for m, psi, d_, _ in pairs))
    return row, f"E={E:.4f} T={T:.4f} " + " | ".join(log)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nE', type=int, default=400)
    ap.add_argument('--procs', type=int, default=4)
    ap.add_argument('--E', type=float, nargs='*')
    a = ap.parse_args()
    t0 = time.time()
    grid = np.array(a.E) if a.E else np.logspace(np.log10(0.05), 0.0, a.nE)
    out_txt = os.path.join(HERE, '..', 'results', 'dr_offdiag_krein.txt')
    out_npy = os.path.join(HERE, '..', 'results', 'dr_offdiag_krein.npy')
    with Pool(a.procs, initializer=_init) as pool:
        res = pool.map(one, grid)
    rows = np.array([r for r, _ in res])
    with open(out_txt, 'w') as fh:
        fh.write(__doc__ + "\n")
        fh.write("columns per block (after E, T): sympl_err, A_tr, B_tr, A_drop, B_drop, W_invariance, "
                 "then pair1 (|mu_N|, psi+, Kdefined), pair2 (same)\n")
        for _, line in res:
            fh.write(line + "\n")
        se = rows[:, [2, 2 + 12, 2 + 24]].max()
        inv = rows[:, [7, 7 + 12, 7 + 24]].max()
        dA = max(np.abs(rows[:, 3 + 12 * k] - rows[:, 5 + 12 * k]).max() for k in range(3))
        dB = max(np.abs(rows[:, 4 + 12 * k] - rows[:, 6 + 12 * k]).max() for k in range(3))
        fh.write(f"\nmax over grid: symplectic error {se:.1e}, W-invariance error {inv:.1e}, "
                 f"|A_tr - A_drop| {dA:.1e}, |B_tr - B_drop| {dB:.1e}\n({time.time() - t0:.0f} s)\n")
    np.save(out_npy, rows)
    print(open(out_txt).read()[-600:])


if __name__ == "__main__":
    main()
