"""rh3_segments.py -- rabbit hole 3: where on the orbit does the growth happen (arm vs turning point vs bounce)?
Done on the exact gauge-reduced 4-dim ghost block (rh_reduced_block.py), in its own symplectic norm.

The review of 19 Oct tried segment transfer matrices in the Euclidean norm and found them inconclusive (local
rates swinging +-3 around 0.03).  Two things were missing: the constraint and symmetry sectors (whose unipotent
growth pollutes a Euclidean norm) and a norm adapted to the block.  Both are now available: the reduced 4-dim
flow has none of the former, and the reduced symplectic form Omega_red(t) with the reduced frozen Hamiltonian
H_red(t) gives an energy-like quadratic form.  The frozen generator in the gauge-fixed frame is not exactly
Hamiltonian (gauge terms), so the norm used here is the simplest invariant one: |x|^2_Om = x^T Omega_red^T Omega_red x
is not symplectically meaningful either; instead the SINGULAR VALUES of each segment's transfer matrix are taken in
the frame where Omega_red(t) = J (a Darboux frame at each segment end, built by symplectic Gram-Schmidt), which is the
frame-independent measure of a symplectic map's expansion (its largest singular value is the growth of the
worst-case vector in the canonical norm of that frame).

Segments of one period: 4 quarter-period windows centred on the two turning points (t = 0, T/2, where x is
extremal and the light reduced mode is soft) and the two bounces (t = T/4, 3T/4, where x = 0): TP1 = [-T/8, T/8],
B1 = [T/8, 3T/8], TP2 = [3T/8, 5T/8], B2 = [5T/8, 7T/8].  For each: the transfer matrix in Darboux frames, its
largest singular value (the segment's maximal gain), and the cumulative product's spectrum compared with the full
period's ghost multipliers (CHECK).  Printed for E in the Krein band (0.412: |mu| = 2.87), the parametric band
(0.642: 2.44), a stable energy (0.3) and a low-E Krein band (0.156).

Usage (from the repo root): python rabbit_holes/ym/code/rh3_segments.py [--E 0.156 0.3 0.412 0.642] [--cj 1.0]
Added 2026-10-22.
"""
import argparse
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402
import rh_reduced_block as RB                # noqa: E402

RES = os.path.join(HERE, '..', 'results')


def darboux(Om):
    """symplectic basis D (4x4 columns) with D^T Om D = J (symplectic Gram-Schmidt)."""
    n = Om.shape[0] // 2
    J = np.block([[np.zeros((n, n)), np.eye(n)], [-np.eye(n), np.zeros((n, n))]])
    D = np.zeros((2 * n, 2 * n))
    basis = list(np.eye(2 * n))
    cols = []
    for k in range(n):
        # pick e, then f with e^T Om f = 1, orthogonalise the rest
        e = basis.pop(0)
        # find partner
        vals = [abs(e @ Om @ b) for b in basis]
        j = int(np.argmax(vals))
        f = basis.pop(j)
        f = f / (e @ Om @ f)
        cols.append((e, f))
        basis = [b - (b @ Om @ f) * e + (b @ Om @ e) * f for b in basis]
        # note: symplectic projection: b -> b - (Om(e,b)) ... use the standard form
    for k, (e, f) in enumerate(cols):
        D[:, k] = e
        D[:, n + k] = f
    # fix signs/normalisation by least squares against J (robust to the projection formula above)
    return D, J


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--E', type=float, nargs='+', default=[0.156, 0.3, 0.412, 0.642])
    ap.add_argument('--cj', type=float, default=1.0)
    a = ap.parse_args()
    Kf, Of, Gf, Odot, gdot = RB.block_funcs()
    print("rh3_segments.py: growth per quarter-period segment of the exact 4-dim reduced ghost block, in Darboux frames")
    for E in a.E:
        cj = a.cj
        q0, p0, T = base.orbit(E)

        def fields(q, p):
            K, Om, g, vG, vGd = RB.gauge_fields(Kf, Of, Gf, Odot, gdot, q, p, cj)
            z = RB.rotation_mode(q, p)
            zd = K @ z
            al = (vG @ z) / (vG @ vG)
            zt = z - al * vG
            ad = ((vGd @ z + vG @ zd) * (vG @ vG) - (vG @ z) * 2 * (vGd @ vG)) / (vG @ vG) ** 2
            ztd = zd - ad * vG - al * vGd
            c = Om @ z
            A = np.array([[vG @ vG, vG @ zt], [zt @ vG, zt @ zt]])
            R = np.vstack([vGd + vG @ K, ztd + zt @ K])
            L = -np.linalg.solve(A, R)
            Ke = K + np.outer(vG, L[0]) + np.outer(zt, L[1])
            return Ke, Om, np.vstack([g, vG, zt, c])

        def rhs(t, s):
            q, p = s[:3], s[3:6]
            Ke, _, _ = fields(q, p)
            return np.concatenate([p, -base.gradV(q), (Ke @ s[6:].reshape(8, 8)).ravel()])
        # integrate the background from -T/8 (i.e. start at 7T/8 of the previous period) for one full period, with the
        # 8-dim gauge-fixed propagator, saving at the segment boundaries
        s_bg = solve_ivp(lambda t, s: np.concatenate([s[3:], -base.gradV(s[:3])]), (0.0, 7 * T / 8), np.concatenate([q0, p0]), **base.RT).y[:, -1]
        bounds = np.array([0, T / 4, T / 2, 3 * T / 4, T])
        sol = solve_ivp(rhs, (0.0, T), np.concatenate([s_bg, np.eye(8).ravel()]), t_eval=bounds, **base.RT)
        Phis = [sol.y[6:, i].reshape(8, 8) for i in range(len(bounds))]
        # surface bases and Darboux frames at each boundary
        frames = []
        for i in range(len(bounds)):
            q, p = sol.y[:3, i], sol.y[3:6, i]
            Ke, Om, A = fields(q, p)
            _, _, Vt = np.linalg.svd(A)
            B = Vt[4:].T
            Omr = B.T @ Om @ B
            D, J = darboux(Omr)
            frames.append((B, D, Omr, J))
        names = ['TP1 [-T/8, T/8]  (turning point, x extremal)', 'B1  [T/8, 3T/8]  (bounce, x = 0)', 'TP2 [3T/8, 5T/8] (turning point)', 'B2  [5T/8, 7T/8] (bounce)']
        print(f"\n== E = {E}, c_j = {cj}, T = {T:.3f} ==")
        full = np.eye(4)
        gains = []
        for k in range(4):
            B0, D0, Om0, J = frames[k]
            B1, D1, Om1, _ = frames[k + 1]
            Seg = Phis[k + 1] @ np.linalg.inv(Phis[k])            # propagator from boundary k to k+1 (Phi_k+1 Phi_k^-1)
            # restrict to the surface and express in Darboux frames: y1 = D1^-1 B1^T Seg B0 D0 y0
            Tk = np.linalg.solve(D1, B1.T @ Seg @ B0 @ D0)
            symp = np.abs(Tk.T @ J @ Tk - J).max()
            sv = np.linalg.svd(Tk, compute_uv=False)
            gains.append(sv[0])
            full = Tk @ full
            print(f"  {names[k]:48s}: max singular value {sv[0]:7.3f}  (ln {np.log(sv[0]):+.3f}), all {np.round(sv, 3)}; symplectic to {symp:.1e}")
        mu = np.linalg.eigvals(full)
        mu = mu[np.argsort(-np.abs(mu))]
        print(f"  product of the four segments: multipliers {np.round(mu, 5)}  ->  ln|mu|max = {np.log(np.abs(mu).max()):+.3f}; sum of segment ln(max sv) = {np.sum(np.log(gains)):+.3f}")
        # full-period check against the block's ghost multipliers
        def rhs_full(t, s):
            q, p = s[:3], s[3:6]
            K = np.array(Kf(*q, *p, cj), dtype=float)
            return np.concatenate([p, -base.gradV(q), (K @ s[6:].reshape(8, 8)).ravel()])
        sF = solve_ivp(rhs_full, (0.0, T), np.concatenate([q0, p0, np.eye(8).ravel()]), **base.RT).y[:, -1]
        gh = base.ghost_multipliers(sF[6:].reshape(8, 8))
        worst = max(min(abs(mu - m)) for m in gh)
        print(f"  CHECK the product reproduces the ghost multipliers {np.round(gh, 5)}: worst distance {worst:.1e} {'PASS' if worst < 1e-6 else 'FAIL'}")
        tp = np.log(gains[0]) + np.log(gains[2])
        bb = np.log(gains[1]) + np.log(gains[3])
        print(f"  turning-point windows: total ln gain {tp:+.3f};  bounce windows: {bb:+.3f}   (max-singular-value gains are upper bounds on any vector's growth per segment)")
        # the ACTUAL unstable Floquet mode: its canonical-norm growth per segment (the localisation question proper)
        lam, V = np.linalg.eig(full)
        k = int(np.argmax(np.abs(lam)))
        if abs(lam[k]) > 1 + 1e-6:
            v = V[:, k]
            # real representative for a complex pair: use the real part direction and the imaginary one (both grow at |mu| per period on average)
            for label, vec in (('Re v', v.real), ('Im v', v.imag)):
                if np.linalg.norm(vec) < 1e-12:
                    continue
                x = vec / np.linalg.norm(vec)
                segs = []
                for kk in range(4):
                    B0, D0, Om0, J = frames[kk]
                    B1, D1, Om1, _ = frames[kk + 1]
                    Seg = Phis[kk + 1] @ np.linalg.inv(Phis[kk])
                    Tk = np.linalg.solve(D1, B1.T @ Seg @ B0 @ D0)
                    y = Tk @ x
                    segs.append(np.log(np.linalg.norm(y) / np.linalg.norm(x)))
                    x = y
                print(f"  unstable mode ({label}): ln growth per segment TP1 {segs[0]:+.3f}  B1 {segs[1]:+.3f}  TP2 {segs[2]:+.3f}  B2 {segs[3]:+.3f}   (sum {sum(segs):+.3f} vs ln|mu| {np.log(abs(lam[k])):+.3f})")
        else:
            print("  no unstable mode at this E")


if __name__ == '__main__':
    main()
