"""rh2_rotation.py -- rabbit hole 2, an ATTEMPT that did not give a canonical quantity (kept for the record; see
notes/RH2_ladder.md "What did not work"): the u12 component passes within 1e-3 of zero along the orbit, so the
winding count below is the mod-2pi angle plus an arbitrary integer.  Not run at size; no committed output.

Original intent -- the missing quantity: the ROTATION of each block-12 ghost pair per half period,
as an integer-valued-plus-fraction phase, so that the tongue ladder can be read as a formula.

B005 gives the angle of each eigenvalue of N = S~ Phi_T/2 only mod 2 pi; unwrapping it in E is ambiguous across a
band (rh2_ladder.py).  Here the phase is measured along the ORBIT instead: for an eigenvector v of N with
eigenvalue lambda, the component u12(t) of the solution Phi_t v (S~ acts as +1 on u12, so u12(T/2) = lambda u12(0))
is a complex scalar whose continuously unwrapped argument over [0, T/2] is theta with e^{i theta} = lambda.  Its
integer part is the winding of u12 around zero, which is canonical as long as the periodic part of the Floquet
solution does not itself encircle zero (checked: the fraction agrees with arg lambda, and theta(E) is continuous).
Off the circle (inside a band) the same construction gives the phase of the growing solution.

Output: results/rh2_rotation.txt, results/rh2_rotation.npz (E, T, the four theta_k/pi, |lambda_k|), results/rh2_rotation.png
(theta/pi of the two pairs vs E, the bands shaded; omega_k = 2 theta_k / T, the effective frequency per period).
Cost: one half-period monodromy plus four dense propagations per E, ~2.5 s each; default nE = 200: ~8 min.
Usage (from the repo root): python rabbit_holes/ym/code/rh2_rotation.py [--nE 200] > rabbit_holes/ym/results/rh2_rotation.txt
Float, DOP853 rtol 1e-12.  Added 2026-10-21 (rabbit-hole session).
"""
import argparse
import os
import re
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402
import krein_b005 as kb                      # noqa: E402

RES = os.path.join(HERE, '..', 'results')
B001 = os.path.join(ROOT, 'islands', 'results', 'dr_offdiag_floquet_typed.txt')
IDX = base.block_idx('12')
STILDE = np.diag(np.array(kb.STILDE['12'] * 4, dtype=float))


def b001_intervals():
    txt = open(B001).read().split('== Unstable intervals')[1]
    out, cur = [], None
    for line in txt.splitlines():
        m = re.match(r'\s*interval \d+: E in \[([\d.]+), ([\d.]+)\]', line)
        if m:
            cur = (float(m.group(1)), float(m.group(2)))
            continue
        m = re.match(r'\s*blk12: max \|mu\| = ([\d.]+) at E=([\d.]+); type sequence (.*?);', line)
        if m and cur:
            out.append((cur[0], cur[1], 'CPLX' if 'CPLX' in m.group(3) else 'REAL', float(m.group(1))))
    return out


def phases(Kf, E, ndense=3000):
    q0, p0, T = base.orbit(E)
    # half-period map
    Ph, M, _, _, T = kb.half_and_full(Kf, E, '12')
    N = STILDE @ Ph
    lam, V = np.linalg.eig(N)
    keep = kb.ghost_four(lam)
    lam4, V4 = lam[keep], V[:, keep]
    # propagate the four eigenvectors (real and imaginary parts) with the block ODE, dense in t
    nv = 4

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, 1.0), dtype=float)[np.ix_(IDX, IDX)]
        X = s[6:].reshape(2 * nv, 8)
        return np.concatenate([p, -base.gradV(q), (X @ K.T).ravel()])
    s0 = np.concatenate([q0, p0, np.concatenate([V4.real.T, V4.imag.T]).ravel()])
    te = np.linspace(0, T / 2, ndense)
    sol = solve_ivp(rhs, (0.0, T / 2), s0, t_eval=te, **base.RT)
    X = sol.y[6:].reshape(2 * nv, 8, ndense)
    out = []
    for k in range(nv):
        u12 = X[k, 0] + 1j * X[nv + k, 0]                  # the u12 component (index 0 of the block state)
        ph = np.unwrap(np.angle(u12))
        theta = ph[-1] - ph[0]
        # consistency: e^{i theta} vs lambda (phase), and |u12(T/2)/u12(0)| vs |lambda|
        resid = abs(np.exp(1j * theta) - lam4[k] / abs(lam4[k]))
        gain = abs(u12[-1] / u12[0]) / abs(lam4[k])
        minabs = np.abs(u12).min() / np.abs(u12).max()
        out.append((theta / np.pi, abs(lam4[k]), np.angle(lam4[k]) / np.pi, resid, gain, minabs))
    return T, lam4, np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nE', type=int, default=200)
    ap.add_argument('--emin', type=float, default=0.05)
    ap.add_argument('--emax', type=float, default=1.0)
    a = ap.parse_args()
    t0 = time.time()
    print(f"rh2_rotation.py  python {sys.version.split()[0]}  numpy {np.__version__}; nE={a.nE} log-spaced in [{a.emin}, {a.emax}], RT={base.RT}")
    print("columns: E | T | for each of the four ghost eigenvectors of N, sorted by theta: theta/pi (unwrapped phase of u12 over T/2), |lambda|, arg lambda/pi | worst phase residual | worst |u12| gain/|lambda| | min |u12| along the half period (relative)")
    Kf = base.load_K()
    grid = np.logspace(np.log10(a.emin), np.log10(a.emax), a.nE)
    TH, LA, TT = [], [], []
    worst = dict(resid=0.0, gain=0.0)
    minabs_all = 1.0
    for E in grid:
        T, lam4, out = phases(Kf, E)
        order = np.argsort(out[:, 0])
        out = out[order]
        TH.append(out[:, 0])
        LA.append(out[:, 1])
        TT.append(T)
        worst['resid'] = max(worst['resid'], out[:, 3].max())
        worst['gain'] = max(worst['gain'], abs(out[:, 4] - 1).max())
        minabs_all = min(minabs_all, out[:, 5].min())
        cols = '  '.join(f"{o[0]:+8.4f} {o[1]:.4f} {o[2]:+.3f}" for o in out)
        print(f"E={E:.5f} T={T:8.4f} | {cols} | {out[:, 3].max():.1e} | {abs(out[:, 4] - 1).max():.1e} | {out[:, 5].min():.2e}")
    TH, LA, TT = np.array(TH), np.array(LA), np.array(TT)
    print(f"CHECK e^(i theta) = lambda/|lambda| for every eigenvector, worst {worst['resid']:.1e} (tol 1e-6) {'PASS' if worst['resid'] < 1e-6 else 'FAIL'}")
    print(f"CHECK |u12(T/2)| / |u12(0)| = |lambda|, worst relative {worst['gain']:.1e} (tol 1e-6) {'PASS' if worst['gain'] < 1e-6 else 'FAIL'}")
    print(f"INFO  smallest |u12| along a half period relative to its max: {minabs_all:.2e} (the unwrapping is safe when this is not ~0)")
    # continuity of theta_k(E): largest jump between neighbours
    jumps = np.abs(np.diff(TH, axis=0)).max(axis=1)
    print(f"INFO  largest jump of a sorted theta/pi between neighbouring energies: {jumps.max():.3f} at E = {grid[np.argmax(jumps) + 1]:.4f} (a jump of ~1 would be a re-sorting across a band)")
    np.savez(os.path.join(RES, 'rh2_rotation.npz'), E=grid, T=TT, theta=TH, mod=LA)

    # the ladder read off: theta_k(E)/pi of the two positive members (the two largest theta), the band edges
    print("\n== the ladder: the two positive-rotation phases at each block-12 band edge (B001) ==")
    ivs = b001_intervals()
    for lo, hi, ty, mu in ivs:
        i0 = max(np.searchsorted(grid, lo) - 1, 0)
        i1 = min(np.searchsorted(grid, hi) + 1, len(grid) - 1)
        print(f"  [{lo:.4f}, {hi:.4f}] {ty}: at E={grid[i0]:.4f} theta/pi = {np.round(TH[i0], 3)}; at E={grid[i1]:.4f} theta/pi = {np.round(TH[i1], 3)}; T/2pi = {TT[i0] / 2 / np.pi:.2f} -> {TT[i1] / 2 / np.pi:.2f}")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
        for lo, hi, ty, mu in ivs:
            for axx in ax:
                axx.axvspan(lo, hi, color=('tab:red' if ty == 'CPLX' else 'tab:blue'), alpha=0.15, lw=0)
        for k in range(4):
            ax[0].plot(grid, TH[:, k], '.', ms=3, label=f'theta_{k}')
        for n in range(-12, 13):
            ax[0].axhline(n, color='k', lw=0.3, alpha=0.4)
            ax[0].axhline(n + 0.5, color='r', lw=0.3, alpha=0.3)
        ax[0].plot(grid, TT / 2 / np.pi, 'k--', lw=1, label='T/(2 pi): omega = 1')
        ax[0].plot(grid, -TT / 2 / np.pi, 'k--', lw=1)
        ax[0].set_ylabel('theta/pi over a half period')
        ax[0].legend(fontsize=7, ncol=3)
        for k in range(4):
            ax[1].plot(grid, 2 * np.pi * TH[:, k] / TT, '.', ms=3, label=f'omega_{k} = 2 theta_{k} / T')
        ax[1].axhline(1.0, color='k', ls='--', lw=1, label='Lee-Wick mass')
        ax[1].set_ylabel('effective frequency')
        ax[1].set_xlabel('E')
        ax[1].set_xscale('log')
        ax[1].legend(fontsize=7, ncol=2)
        ax[0].set_title('rabbit hole 2: rotation of the block-12 ghost solutions per half period (grey lines: N = +-1, red: N = +-i); bands: red Krein, blue parametric', fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, 'rh2_rotation.png'), dpi=130)
        print("figure: results/rh2_rotation.png")
    except Exception as e:                                            # noqa: BLE001
        print(f"(no figure: {e})")
    print(f"DONE {time.time() - t0:.0f}")


if __name__ == '__main__':
    main()
