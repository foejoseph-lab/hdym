"""rh12_direction.py -- rabbit hole 12: does a Krein quadruplet's runaway ROTATE in the off-diagonal state space while a
real pair's runs along a fixed direction?  From the state snapshots the kernel records at the crossing events
(islands/results/esc_E_snap_krein.npz, esc_E_snap_param.npz: section E, block-12 seeds, snap (nE, N, 5, 36)).

Per trajectory: the block-12 off-diagonal state (8 numbers: u12, u21 and their three derivatives, from the 36-component
state) at the two snapshots |w| = 0.1 (index 4) and |w| = 1 (index 0), as unit vectors d_a, d_b; the angle
theta_nl = arccos |d_a . d_b| between them (sign-free: a real pair may flip sign).  The linear prediction: the same
8-dim state propagated with the block's linearised flow from t_a to t_b (the background integrated from the DR point
to t_a, then the block's K along the orbit), giving d_b^lin and theta_lin = arccos |d_a . d_b^lin|.  A rotating
runaway: theta_nl tracks theta_lin, which for a quadruplet is set by the Krein angle arg(mu) per period times the
elapsed periods (mod the periodic part); a directional one: both near 0 (or the periodic part's own value).  Also
reported: the Euclidean angle between d_b and the unstable Floquet eigenvector(s) at the phase of t_b (the mode the
nonlinear state has converged to), and the periodic-part contribution (theta_lin for a vector started ON the mode).

Caveat: Euclidean angles in the (u, u', u'', u''') coordinates are not canonical; the comparison nonlinear-vs-linear in
the same coordinates is what carries the information, not the absolute angles.

Usage (from the repo root): python rabbit_holes/ym/code/rh12_direction.py [--files esc_E_snap_krein esc_E_snap_param]
Added 2026-10-22; waits on the two runs of hole12_queue.md.
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

ISL = os.path.join(ROOT, 'islands', 'results')
RES = os.path.join(HERE, '..', 'results')
OD_M = [1, 2, 3, 5, 6, 7]                    # full-state component indices m = 3i+a of the off-diagonal entries (rung1 ulist order)
BLOCK12_J = [0, 2]                           # ulist positions of (u12, u21)
IDX = base.block_idx('12')


def block12_state(s36):
    """(u12, u21, u12', u21', u12'', u21'', u12''', u21''') from the 36-component full state (s[9k + m], k = derivative order)."""
    return np.array([s36[9 * k + OD_M[j]] for k in range(4) for j in BLOCK12_J])


def linear_propagate(Kf, E, ta, tb, v8):
    """propagate the block-12 vector v8 with the linearised flow from t = ta to t = tb along the DR orbit."""
    q0, p0, T = base.orbit(E)
    bg = solve_ivp(lambda t, s: np.concatenate([s[3:], -base.gradV(s[:3])]), (0.0, ta), np.concatenate([q0, p0]), **base.RT).y[:, -1] if ta > 0 else np.concatenate([q0, p0])

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, 1.0), dtype=float)[np.ix_(IDX, IDX)]
        return np.concatenate([p, -base.gradV(q), K @ s[6:]])
    out = solve_ivp(rhs, (ta, tb), np.concatenate([bg, v8]), **base.RT).y[:, -1]
    return out[6:], T


def unstable_modes(Kf, E):
    _, _, Ms, _ = base.monodromies(Kf, E, 1.0)
    mu, V = np.linalg.eig(Ms['12'])
    keep = np.argsort(np.abs(mu - 1.0))[4:]
    unst = [k for k in keep if abs(mu[k]) > 1 + 1e-6]
    return mu[unst], V[:, unst]


def angle(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    return np.degrees(np.arccos(min(1.0, abs(a @ b))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--files', nargs='+', default=['esc_E_snap_krein', 'esc_E_snap_param'])
    ap.add_argument('--nmax', type=int, default=64, help='trajectories per energy to propagate linearly (cost ~1 s each)')
    a = ap.parse_args()
    Kf = base.load_K()
    print("rh12_direction.py: direction of the block-12 off-diagonal state between the |w| = 0.1 and |w| = 1 crossings, nonlinear vs linear")
    for tag in a.files:
        path = os.path.join(ISL, tag + '.npz')
        if not os.path.exists(path):
            print(f"\n{tag}: not found ({path})")
            continue
        d = np.load(path, allow_pickle=True)
        E, snap, tc, st, te = d['E'], d['snap'], d['tcross'], d['status'], d['tesc']
        print(f"\n== {tag}: {len(E)} energies x {snap.shape[1]} seeds ==")
        for i in range(len(E)):
            mu, V = unstable_modes(Kf, E[i])
            kind = 'Krein quadruplet' if len(mu) == 2 and abs(mu[0].imag) > 1e-6 else ('real pair' if len(mu) == 1 else ('stable' if len(mu) == 0 else f'{len(mu)} unstable'))
            _, _, T = base.orbit(E[i])
            rows = []
            for n in range(min(snap.shape[1], a.nmax)):
                ta, tb = tc[i, n, 4], tc[i, n, 0]
                if ta <= 0 or tb <= ta:
                    continue
                sa, sb = snap[i, n, 4].astype(float), snap[i, n, 0].astype(float)
                if not (np.isfinite(sa).all() and np.isfinite(sb).all()):
                    continue
                va, vb = block12_state(sa), block12_state(sb)
                if np.linalg.norm(va) < 1e-12 or np.linalg.norm(vb) < 1e-12:
                    continue
                vb_lin, _ = linear_propagate(Kf, E[i], ta, tb, va)
                th_nl = angle(va, vb)
                th_lin = angle(va, vb_lin)
                th_nl_vs_lin = angle(vb, vb_lin)
                # angle of the nonlinear end state to the unstable subspace at phase t_b (propagate the eigenvectors to t_b mod T)
                phase = tb % T
                sub = []
                for k in range(V.shape[1]):
                    vk, _ = linear_propagate(Kf, E[i], 0.0, phase, V[:, k].real) if phase > 0 else (V[:, k].real, T)
                    sub.append(vk)
                    if abs(mu[k].imag) > 1e-6:
                        vki, _ = linear_propagate(Kf, E[i], 0.0, phase, V[:, k].imag) if phase > 0 else (V[:, k].imag, T)
                        sub.append(vki)
                if len(sub) == 0:                                 # no unstable mode at this E (smoke energies): subspace angle undefined
                    rows.append((th_nl, th_lin, th_nl_vs_lin, np.nan, (tb - ta) / T))
                    continue
                S = np.array(sub).T
                Q, _ = np.linalg.qr(S)
                proj = np.linalg.norm(Q.T @ (vb / np.linalg.norm(vb)))
                rows.append((th_nl, th_lin, th_nl_vs_lin, np.degrees(np.arccos(min(1.0, proj))), (tb - ta) / T))
            if not rows:
                print(f"  E={E[i]:.4f} ({kind}): no usable trajectories")
                continue
            R = np.array(rows)
            mudesc = f"|mu| = {np.abs(mu).max():.3f}, arg/pi = {np.round(np.angle(mu) / np.pi, 3)}" if len(mu) else "no unstable block-12 mode"
            print(f"  E={E[i]:.4f} ({kind}, {mudesc}): {len(R)} trajectories; elapsed between crossings median {np.median(R[:, 4]):.1f} periods")
            print(f"     angle(d_0.1, d_1) nonlinear: median {np.median(R[:, 0]):5.1f} deg [p10 {np.percentile(R[:, 0], 10):5.1f}, p90 {np.percentile(R[:, 0], 90):5.1f}]   linear prediction: median {np.median(R[:, 1]):5.1f} [{np.percentile(R[:, 1], 10):5.1f}, {np.percentile(R[:, 1], 90):5.1f}]")
            print(f"     angle(nonlinear end, linear end): median {np.median(R[:, 2]):5.1f} deg [p10 {np.percentile(R[:, 2], 10):5.1f}, p90 {np.percentile(R[:, 2], 90):5.1f}]   angle(nonlinear end, unstable subspace at that phase): median {np.nanmedian(R[:, 3]):5.1f} [{np.nanpercentile(R[:, 3], 10):5.1f}, {np.nanpercentile(R[:, 3], 90):5.1f}]")
            np.savez(os.path.join(RES, f'rh12_{tag}_E{E[i]:.4f}.npz'), rows=R, mu=mu)
    print("\nreading: 'nonlinear end vs linear end' small => the runaway follows the linear mode's direction (rotating if the mode is a quadruplet, which the")
    print("linear angle itself shows); large => the nonlinear state has left the unstable subspace by the time |w| = 1 (the diagonal ghost has taken over).")


if __name__ == '__main__':
    main()
