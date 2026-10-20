"""dr_offdiag_halfperiod.py -- is the off-diagonal monodromy on the DR orbit the square of a
symplectic half-period map?  (Review of the 19 Oct brief, step 1.)

The DR orbit returns to y = 0 at t = T/2 with (x, p_x) restored and p_y reversed
(islands/README.md, lemma RET), i.e.  q(t + T/2) = S q(t),  S = diag(1, -1, 1) on (x, y, z).
S is a symmetry of V, so the linearised off-diagonal system inherits it: for each 8-dim block
there is a constant sign matrix St with  K(S q, S p) = St K(q, p) St,  and then the monodromy
over the full period is  M = N^2  with  N = St Phi_{T/2}.  Consequences for the brief's claims:
a multiplier -1 of M is +-i of N, not a fold, so a single pair of M cannot leave the circle at -1
(no REAL- is forced), and a = mu + 1/mu touches -2 tangentially.

Printed per E: the half-period closure |q(T/2) - S q(0)|, the sign pattern St found from K,
|N^2 - M| / |M|, and the multipliers of N and M.  Also, at E = 1/2, the Jordan structure of
the gauge sector (ranks of (M - 1)^k) so that the 'four nearest 1' rule is checked against fact.

Usage (from islands/):  python code/dr_offdiag_halfperiod.py [E ...]
Added 2026-10-19 (review).
"""
import itertools
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dr_offdiag_floquet as base  # noqa: E402

S_BG = np.array([1.0, -1.0, 1.0])          # y -> -y


def sign_matrix(Kf, name, cj=1.0, trials=5, seed=0):
    """Find the diagonal sign matrix St (one sign per component, same at all four derivative
    levels) with K(Sq, Sp) = St K(q, p) St on the block, by testing random backgrounds."""
    rng = np.random.default_rng(seed)
    idx = base.block_idx(name)
    best = None
    for signs in itertools.product((1.0, -1.0), repeat=2):
        s = np.array([signs[0], signs[1]] * 4)
        St = np.diag(s)
        err = 0.0
        for _ in range(trials):
            q = rng.normal(size=3); p = rng.normal(size=3)
            K1 = np.array(Kf(*q, *p, cj), dtype=float)[np.ix_(idx, idx)]
            K2 = np.array(Kf(*(S_BG * q), *(S_BG * p), cj), dtype=float)[np.ix_(idx, idx)]
            err = max(err, np.max(np.abs(K2 - St @ K1 @ St)) / max(1e-300, np.max(np.abs(K1))))
        if best is None or err < best[0]:
            best = (err, s)
    return best


def half_maps(Kf, E, cj=1.0):
    """Phi_{T/2} and Phi_T for the three blocks, one integration each (checks closure both ways)."""
    q0, p0, T = base.orbit(E)
    nb = len(base.BLOCKS)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, cj), dtype=float)
        out = [p, -base.gradV(q)]
        off = 6
        for b in base.BLOCKS:
            idx = base.block_idx(b)
            Mb = s[off:off + 64].reshape(8, 8)
            out.append((K[np.ix_(idx, idx)] @ Mb).ravel())
            off += 64
        return np.concatenate(out)

    s0 = np.concatenate([q0, p0] + [np.eye(8).ravel()] * nb)
    sol = solve_ivp(rhs, (0.0, T), s0, t_eval=[T / 2, T], **base.RT)
    sh, sT = sol.y[:, 0], sol.y[:, 1]
    half_clo = np.max(np.abs(np.concatenate([sh[:3] - S_BG * q0, sh[3:6] - S_BG * p0])))
    full_clo = np.max(np.abs(sT[:6] - s0[:6]))
    Ph, Pf = {}, {}
    off = 6
    for b in base.BLOCKS:
        Ph[b] = sh[off:off + 64].reshape(8, 8)
        Pf[b] = sT[off:off + 64].reshape(8, 8)
        off += 64
    return T, half_clo, full_clo, Ph, Pf


def fmt(mu):
    mu = mu[np.argsort(-np.abs(mu))]
    return "  ".join(f"{m.real:+.5f}{m.imag:+.5f}i" for m in mu)


def main():
    Es = [float(a) for a in sys.argv[1:]] or [0.5, 0.3, 0.412, 0.0761, 0.7133]
    Kf = base.load_K()
    print("Symmetry of K under y -> -y, per block (sign pattern on (u_first, u_second), relative error):")
    St = {}
    for b in base.BLOCKS:
        err, s = sign_matrix(Kf, b)
        St[b] = np.diag(s)
        print(f"  block {b}: signs {s[:2]}  err {err:.1e}")
    for E in Es:
        T, hc, fc, Ph, Pf = half_maps(Kf, E)
        print(f"\nE = {E}:  T = {T:.6f}   |q(T/2) - S q(0)| = {hc:.1e}   closure(T) = {fc:.1e}")
        for b in base.BLOCKS:
            N = St[b] @ Ph[b]
            M = Pf[b]
            rel = np.linalg.norm(N @ N - M) / np.linalg.norm(M)
            print(f"  block {b}:  |N^2 - M|/|M| = {rel:.1e}   det N - 1 = {np.linalg.det(N) - 1:+.1e}")
            print(f"     mu(N): {fmt(np.linalg.eigvals(N))}")
            print(f"     mu(M): {fmt(np.linalg.eigvals(M))}")
            if E == Es[0]:
                I = np.eye(8)
                ranks = [np.linalg.matrix_rank(np.linalg.matrix_power(M - I, k), tol=1e-6) for k in (1, 2, 3, 4)]
                print(f"     rank (M - 1)^k, k=1..4 (tol 1e-6): {ranks}   -> multiplicity of 1 = {8 - ranks[-1]}, "
                      f"largest Jordan block = {next(k for k in range(1, 5) if ranks[k - 1] == ranks[-1])}")
                print(f"     tr M - 4 = {np.trace(M) - 4:+.8f}  vs  A from drop-rule = "
                      f"{np.sum(base.ghost_multipliers(M)).real:+.8f}")


if __name__ == "__main__":
    main()
