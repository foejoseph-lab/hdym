"""plot_cantori.py -- figures for cantori_residence.py (rabbit hole 9, last section).  No integration; reads the npz.

    python islands/code/plot_cantori.py
writes islands/results/cantori_survivors.png (start points still near the island at t = 300, 3,000, 30,000 and
never left by tmax, in the island's round coordinates, coloured by exit time) and islands/results/cantori_tail.png
(survival of the exiting starts on log-log axes with the bulk's exponential, and the section points of the longest-
lived orbit).  Round coordinates: u = sqrt(w) * (U^T dz), with S = U diag(w) U^T the invariant form, so rho = |u|.
Exit times are for the smallest rho_out (1.5e-3).  Added 2026-10-29.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
BLUE = LinearSegmentedColormap.from_list('blue', ['#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
ORANGE = '#eb6834'
EDGE_IN, WALL = 2.15e-4, 3.06e-4          # last radius intact in cantori_edge.txt; trace's 10th-percentile rho

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK2, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.facecolor': SURF, 'figure.facecolor': SURF, 'savefig.facecolor': SURF})


def round_coords(y, vy, z, S):
    w, U = np.linalg.eigh(S)
    dz = np.stack([y - z[0], vy - z[1]], 1)
    v = dz @ U
    return v[:, 0] * np.sqrt(w[0]), v[:, 1] * np.sqrt(w[1])


def circles(ax, lim):
    th = np.linspace(0, 2 * np.pi, 400)
    for r, ls, lab in ((EDGE_IN, ':', 'last intact edge-scan radius'), (WALL, '--', 'ρ = 3.06×10⁻⁴ (10 % of the long-lived orbits\' points lie inside)')):
        ax.plot(r * np.cos(th) * 1e4, r * np.sin(th) * 1e4, ls, color=INK2, lw=0.9, label=lab)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
    ax.set_xlabel('u₁ (×10⁻⁴)'); ax.set_ylabel('u₂ (×10⁻⁴)')


def main():
    Z = np.load(os.path.join(RES, 'cantori_survive.npz'))
    f = Z['first'][:, 0]; tmax = float(Z['tmax']); z, S = Z['z'], Z['S']
    u1, u2 = round_coords(Z['y0'], Z['vy0'], z, S)
    lim = float(Z['rb']) * 1e4 * 1.05
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.6), constrained_layout=True)
    norm = LogNorm(vmin=30, vmax=tmax)
    for ax, T, title in zip(axs, (300, 3000, 30000, None),
                            ('still near at t = 300', 't = 3,000', 't = 30,000', f'never left by t = {tmax:,.0f}')):
        m = (f > T) if T is not None else np.isinf(f)
        ex = m & np.isfinite(f)
        reg = m & np.isinf(f)
        ax.scatter(u1[reg] * 1e4, u2[reg] * 1e4, s=0.3, color='#c3c2b7', lw=0, rasterized=True)
        sc = ax.scatter(u1[ex] * 1e4, u2[ex] * 1e4, c=f[ex], cmap=BLUE, norm=norm, s=3 if ex.sum() < 2000 else 0.6,
                        lw=0, rasterized=True)
        circles(ax, lim)
        ax.set_title(f'{title}\n{int(ex.sum()):,} exit later, {int(reg.sum()):,} never (grey)', color=INK, fontsize=9)
    axs[0].legend(loc='lower left', fontsize=7, frameon=False)
    cb = fig.colorbar(sc, ax=axs, shrink=0.85, pad=0.01); cb.set_label('exit time (ρ > 1.5×10⁻³)', color=INK2)
    fig.suptitle('DR island, planar ghost-free Yang–Mills, E = 1/2: where the starts that stay come from '
                 f'({len(f):,} starts, ρ₀ ∈ [{float(Z["ra"]):g}, {float(Z["rb"]):g}])', color=INK, fontsize=10)
    fig.savefig(os.path.join(RES, 'cantori_survivors.png'), dpi=160)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    ts = np.geomspace(30, tmax, 200)
    for j, col, lab in ((0, '#2a78d6', 'ρ_out = 1.5×10⁻³'), (len(Z['rho_out']) - 1, ORANGE, 'ρ_out = 10⁻¹')):
        fj = Z['first'][:, j]; fe = fj[np.isfinite(fj)]
        N = np.array([np.sum(fe > t) for t in ts])
        a1.loglog(ts[N > 0], N[N > 0], color=col, lw=2, label=lab)
    fe = f[np.isfinite(f)]
    n300, n500 = np.sum(fe > 300), np.sum(fe > 500)
    rate = np.log(n300 / n500) / 200
    tb = np.geomspace(300, 3000, 50)
    a1.loglog(tb, n300 * np.exp(-rate * (tb - 300)), '--', color=INK2, lw=1.2,
              label=f'bulk exponential (rate {rate:.4f}, from t = 300–500)')
    a1.set_ylim(0.8, 3e5); a1.set_xlabel('t'); a1.set_ylabel('exiting starts with exit time > t')
    a1.grid(True, which='major', color=GRID, lw=0.6); a1.legend(frameon=False, fontsize=8)
    a1.set_title('exit-time tail (the 28,164 starts that never exit are not counted)', color=INK, fontsize=9)
    T = np.load(os.path.join(RES, 'cantori_trace.npz'))
    P = T['traj0']
    sgn = np.where(np.abs(P[:, 1] - z[0]) < np.abs(P[:, 1] + z[0]), 1.0, -1.0)      # fold -c onto c
    v1, v2 = round_coords(sgn * P[:, 1], sgn * P[:, 2], z, S)
    sc = a2.scatter(v1 * 1e4, v2 * 1e4, c=P[:, 0], cmap=BLUE, s=0.8, lw=0, rasterized=True)
    circles(a2, 12)
    cb = fig.colorbar(sc, ax=a2, shrink=0.85, pad=0.01); cb.set_label('t', color=INK2)
    a2.set_title(f'longest-lived orbit: {len(P):,} section points to exit at t = {P[-1, 0]:,.0f}', color=INK, fontsize=9)
    a2.legend(loc='lower left', fontsize=7, frameon=False)
    fig.savefig(os.path.join(RES, 'cantori_tail.png'), dpi=160)
    print(f'bulk rate {rate:.5f}; exiting starts past t = 1000: {np.sum(fe > 1000)}, bulk exponential predicts '
          f'{n300 * np.exp(-rate * 700):.1f}; past 2000: {np.sum(fe > 2000)} vs {n300 * np.exp(-rate * 1700):.2g}')


if __name__ == '__main__':
    main()
