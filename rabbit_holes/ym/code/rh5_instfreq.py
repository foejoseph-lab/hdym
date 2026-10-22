"""rh5_instfreq.py -- the reviewer's sweep-vs-beat test: the instantaneous frequency of each ghost pair's Floquet
solution along one period, omega_inst(t) = sum_c Im(xdot_c conj x_c) / sum_c |x_c|^2 over the field components
(u12, u21, u12', u21'), i.e. the energy-weighted phase rate.  A swept oscillator (near 1 for part of the orbit, near 0
for the rest: WKB zeroth order) shows omega_inst sweeping once per period; a true superposition of two steady lines
shows a beat at their difference instead.  Printed against the background (x, y) and |p_y| for orientation.
Usage: python rabbit_holes/ym/code/rh5_instfreq.py [--E 0.1 0.3 0.6] > rabbit_holes/ym/results/rh5_instfreq.txt   (writes .png)
Added 2026-10-21 (late).
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

RES = os.path.join(HERE, '..', 'results')
IDX = base.block_idx('12')


def run(Kf, E, nsamp=2000):
    q0, p0, T = base.orbit(E)

    def rhsM(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, 1.0), dtype=float)[np.ix_(IDX, IDX)]
        return np.concatenate([p, -base.gradV(q), (K @ s[6:].reshape(8, 8)).ravel()])
    sT = solve_ivp(rhsM, (0.0, T), np.concatenate([q0, p0, np.eye(8).ravel()]), **base.RT).y[:, -1]
    M = sT[6:].reshape(8, 8)
    mu, V = np.linalg.eig(M)
    keep = np.argsort(np.abs(mu - 1.0))[4:]
    g, Vg = mu[keep], V[:, keep]
    sel = [k for k in range(4) if g[k].imag > 1e-9] if (g.imag > 1e-9).sum() == 2 else list(np.argsort(-np.abs(g))[:2])
    sel.sort(key=lambda k: abs(np.angle(g[k])))
    nv = len(sel)
    Vs = Vg[:, sel]

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, 1.0), dtype=float)[np.ix_(IDX, IDX)]
        X = s[6:].reshape(2 * nv, 8)
        return np.concatenate([p, -base.gradV(q), (X @ K.T).ravel()])
    s0 = np.concatenate([q0, p0, np.concatenate([Vs.real.T, Vs.imag.T]).ravel()])
    te = np.linspace(0, T, nsamp)
    sol = solve_ivp(rhs, (0.0, T), s0, t_eval=te, **base.RT)
    q = sol.y[:3]
    p = sol.y[3:6]
    X = sol.y[6:].reshape(2 * nv, 8, nsamp)
    out = []
    for k in range(nv):
        x = (X[k] + 1j * X[nv + k])[:4]                           # field components
        xd = (X[k] + 1j * X[nv + k])[2:6]                         # their time derivatives are the next two entries... (u', u'') for (u, u')
        w = (np.imag(xd * np.conj(x)).sum(axis=0)) / (np.abs(x) ** 2).sum(axis=0)
        out.append((g[sel[k]], w))
    return T, te, q, p, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--E', type=float, nargs='+', default=[0.1, 0.3, 0.6])
    a = ap.parse_args()
    Kf = base.load_K()
    print("rh5_instfreq.py: instantaneous frequency of the ghost Floquet solutions along one period (energy-weighted phase rate of the field components)")
    allr = []
    for E in a.E:
        T, te, q, p, out = run(Kf, E)
        allr.append((E, T, te, q, p, out))
        print(f"\nE={E}, T={T:.2f}: columns t/T | x | y | pair 1: omega_inst | pair 2: omega_inst")
        for i in range(0, len(te), len(te) // 24):
            print(f"  {te[i] / T:5.3f} | {q[0, i]:+6.3f} | {q[1, i]:+6.3f} | {out[0][1][i]:+7.3f} | {out[1][1][i]:+7.3f}")
        for k, (m, w) in enumerate(out):
            print(f"  pair {k + 1}: mu-angle {np.angle(m) / np.pi:+.3f} pi; omega_inst: mean {w.mean():+.3f}, min {w.min():+.3f}, max {w.max():+.3f}, fraction of the period with |omega_inst| < 0.3: {(np.abs(w) < 0.3).mean():.2f}, with |omega_inst| > 0.7: {(np.abs(w) > 0.7).mean():.2f}; mean over period / (arg mu / T) = {w.mean() / (np.angle(m) / T):.2f}")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(len(a.E), 1, figsize=(9, 2.6 * len(a.E)), sharex=True)
        for row, (E, T, te, q, p, out) in zip(np.atleast_1d(ax), allr):
            row.plot(te / T, out[0][1], 'C0', lw=1, label='pair 1 omega_inst')
            row.plot(te / T, out[1][1], 'C1', lw=1, label='pair 2 omega_inst')
            row.plot(te / T, q[1] / np.abs(q).max(), 'k--', lw=0.7, label='y(t) (scaled)')
            row.plot(te / T, q[0] / np.abs(q).max(), 'k:', lw=0.7, label='x(t) (scaled)')
            row.axhline(1, color='grey', lw=0.5)
            row.axhline(0, color='grey', lw=0.5)
            row.set_ylabel(f'E={E}')
            row.set_ylim(-2.5, 2.5)
        np.atleast_1d(ax)[0].legend(fontsize=7, ncol=4)
        np.atleast_1d(ax)[-1].set_xlabel('t / T')
        fig.suptitle('instantaneous frequency of the block-12 ghost solutions along the DR orbit', fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, 'rh5_instfreq.png'), dpi=130)
        print("figure: results/rh5_instfreq.png")
    except Exception as e:                                        # noqa: BLE001
        print(f"(no figure: {e})")


if __name__ == '__main__':
    main()
