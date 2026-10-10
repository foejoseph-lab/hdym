"""make_fig_landscape.py -- Fig: where the ghost grows, in the (depth D, transverse action J) plane
of the channel description, at one energy, for three ghosts:
  (a) colour-diagonal ghost (hdym_gpu.py --ghost diag)           -> results_gpu/landscape_E0.1_diag.npz
  (b) off-diagonal block, (D^mu F_mu0)^2 term off (--cj 0)       -> results_rung1/landscape_E0.1_cj0.npz
  (c) off-diagonal block, the theory (--cj 1)                    -> results_rung1/landscape_E0.1_cj1.npz
Colour: share of the total log-growth per unit area of the (D, J) plane (growth density).  Thin
grey: time density on the same bins (where the trajectory *is*).  Lines: the turning-point
relation D_max J = 1 (dashed) and the resonant depths D = E^{-1/4} (diagonal: K = y^2 chirps at
2|x|) and D = 2 E^{-1/4} (off-diagonal: xy chirps at |x|), dotted.
Usage: python make_fig_landscape.py --E 0.1 [--files a.npz b.npz c.npz]
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 9})
ap = argparse.ArgumentParser()
ap.add_argument('--E', type=float, default=0.1)
ap.add_argument('--files', nargs=3, default=None)
ap.add_argument('--out', default='fig_landscape')
a = ap.parse_args()
E = a.E
files = a.files or [f'results_gpu/landscape_E{E:.4g}_diag.npz',
                    f'results_rung1/landscape_E{E:.4g}_cj0.npz',
                    f'results_rung1/landscape_E{E:.4g}_cj1.npz']
labels = ['(a) diagonal ghost', r'(b) off-diagonal, $(D^\mu F_{\mu0})^2$ off', '(c) off-diagonal, the theory']
fig, ax = plt.subplots(1, 3, figsize=(9.6, 3.2), sharey=True)
for k, (f, lab) in enumerate(zip(files, labels)):
    d = np.load(f, allow_pickle=True)
    hg, ht = d['hist_g'], d['hist_t']
    De, Je = d['D_edges'], d['J_edges']
    dD, dJ = De[1] - De[0], Je[1] - Je[0]
    tot = hg.sum()
    dens = hg / tot / (dD * dJ)                       # share of growth per unit (D, J) area
    tdens = ht / ht.sum() / (dD * dJ)
    Dc, Jc = 0.5 * (De[1:] + De[:-1]), 0.5 * (Je[1:] + Je[:-1])
    vmax = np.percentile(dens[dens > 0], 99.5) if np.any(dens > 0) else 1.0
    im = ax[k].pcolormesh(Jc, Dc, np.clip(dens, 0, None), cmap='magma', vmin=0, vmax=vmax, shading='nearest')
    ax[k].contour(Jc, Dc, tdens, levels=4, colors='w', linewidths=0.5, alpha=0.6)
    Js = np.linspace(0.13, 1, 50)
    ax[k].plot(Js, 1 / Js, 'w--', lw=0.8)
    ax[k].axhline(E ** -0.25, color='c', ls=':', lw=0.8)
    ax[k].axhline(2 * E ** -0.25, color='y', ls=':', lw=0.8)
    g = float(d['lam']) if 'lam' in d.files else float(d['g1'])
    share_ch = hg[Dc > 3].sum() / tot
    ax[k].set_title(f'{lab}\n' + rf'$\lambda={g:.2e}$, {100 * share_ch:.0f}% of growth at $D>3$', fontsize=8)
    ax[k].set_xlabel(r'transverse action $J\,E^{-3/4}$')
    ax[k].set_ylim(0, 8)
    ax[k].set_xlim(0, 1)
ax[0].set_ylabel(r'depth $D=\max|q|\,E^{-1/4}$')
cb = plt.colorbar(im, ax=ax.ravel().tolist(), fraction=0.025, pad=0.02)
cb.set_label('growth density (share of log-growth per unit area)')
plt.savefig(a.out + '.pdf', bbox_inches='tight')
plt.savefig(a.out + '.png', dpi=150, bbox_inches='tight')
print('wrote', a.out + '.pdf/.png')
