"""offdiag_escape_nodal.py -- the orient maps of offdiag_escape.py against the linear prediction of their structure.

A seed direction s (unit vector in the span of the map's three basis directions e1, e2, e3) grows along the
unstable Floquet mode k of a block at the rate ln|mu_k|/T with amplitude |w_k . s|, w_k the left eigenvector of
the block's monodromy (rows of inv(V), V the right eigenvectors).  Where that amplitude vanishes the linear phase
must start from second-order leakage and the escape is delayed.  So the nodal set of the unstable subspace is
drawn on the (theta, phi) sphere: for a real parametric pair (one unstable direction) it is a great circle; for a
Krein quadruplet (two unstable directions, mu and conj mu) it is two antipodal points.  This script computes the
unstable left eigenvectors from dr_offdiag_floquet.monodromies at the map's E, projects them on the stored basis,
overlays the predicted nodal set on the smoothed log10 t_esc map (5x5 boxcar mean, scipy.ndimage.uniform_filter,
reflecting edges; non-escaped cells filled with the map mean before smoothing), and prints the correlation of log t_esc with
log |projection| over the map.  No fitting.

Usage: python islands/code/offdiag_escape_nodal.py esc_orient_030:23 esc_orient_0156:12   (tag:block, from islands/results/)
Writes islands/results/<tag>_nodal.png per map (the figures are not committed; the numbers are in islands/README.md:
corr -0.44 for the real pair at E = 0.30, -0.11 for the quadruplet at 0.156).  Added 2026-10-21.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dr_offdiag_floquet as F                                    # noqa: E402

OD = [1, 2, 3, 5, 6, 7]                                           # rung1 ulist order, m = 3i + a
BLOCK_M = {'12': [1, 3], '13': [2, 6], '23': [5, 7]}


def seed_to_block(v24, b):
    """the 24-dim seed (k-major, ulist-minor) restricted to a block's 8-dim state in block_idx order"""
    jj = [OD.index(m) for m in BLOCK_M[b]]
    return np.array([v24[6 * k + j] for k in range(4) for j in jj])


def nodal(tag, b, results, Kf):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.ndimage import uniform_filter
    d = np.load(os.path.join(results, tag + '.npz'), allow_pickle=True)
    te, th, ph, basis, E = d['tesc'], d['theta'], d['phi'], d['basis'], float(d['E'])
    T, clo, Ms, _ = F.monodromies(Kf, E, 1.0)
    mu, V = np.linalg.eig(Ms[b])
    W = np.linalg.inv(V)
    unst = np.where(np.abs(mu) > 1 + 1e-6)[0]
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    amp = np.zeros_like(TH)
    Eb = [seed_to_block(basis[i], b) for i in range(3)]
    for k in unst:
        c = [W[k] @ e for e in Eb]
        amp += np.abs(np.cos(TH) * c[0] + np.sin(TH) * np.cos(PH) * c[1] + np.sin(TH) * np.sin(PH) * c[2]) ** 2
    amp = np.sqrt(amp)
    L = np.where(te > 0, np.log10(np.where(te > 0, te, 1)), np.nan)
    Ls = uniform_filter(np.nan_to_num(L, nan=np.nanmean(L)), size=5)
    r = np.corrcoef(np.nan_to_num(L, nan=np.nanmean(L)).ravel(), np.log10(amp).ravel())[0, 1]
    kind = 'Krein quadruplet, 2 unstable directions -> nodal points' if len(unst) == 2 else 'real pair, 1 unstable direction -> nodal great circle'
    print(f"{tag}: E={E} block {b}: unstable multipliers {np.round(mu[unst], 4)}; {kind}; corr(log t_esc, log |proj|) = {r:+.3f}")
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.4), facecolor='#fcfcfb', gridspec_kw=dict(width_ratios=[1.25, 1]))
    ax = axs[0]
    im = ax.imshow(Ls, origin='lower', aspect='auto', extent=[0, 2, 0, 1], cmap='viridis')
    plt.colorbar(im, ax=ax, label='log10 t_esc, 5×5 boxcar mean (non-escaped cells = map mean)')
    ax.contour(PH / np.pi, TH / np.pi, np.log10(amp / amp.max()), levels=[-2.5, -2, -1.5], colors='white', linewidths=[1.6, 1.0, 0.6])
    ax.set_title(f'{tag}: E = {E}, block {b}\n{kind}\nwhite: predicted nodal set from the linear monodromy; corr(log t, log |proj|) = {r:+.2f}', fontsize=9)
    ax.set_xlabel('φ/π')
    ax.set_ylabel('θ/π')
    # escape time against the projection: the lower envelope is the linear phase, the scatter above it the loiter
    tmax = float(d['tmax'])
    ax = axs[1]
    ax.scatter(np.log10(amp).ravel(), np.where(te > 0, te, tmax).ravel(), s=1.5, alpha=0.25, color='#2a78d6', lw=0)
    ax.axhline(tmax, color='#52514e', lw=0.8, ls='--')
    ax.set_xlabel('log10 |projection of the seed on the unstable subspace| (0 = largest)')
    ax.set_ylabel('t_esc  (at the dashed line: not escaped by tmax)')
    ax.set_title('lower envelope = linear phase set by the projection; scatter above it = the loiter', fontsize=9)
    ax.set_facecolor('#fcfcfb')
    fig.savefig(os.path.join(results, tag + '_nodal.png'), dpi=140, bbox_inches='tight')


if __name__ == '__main__':
    results = os.path.join(HERE, '..', 'results')
    Kf = F.load_K()
    for arg in sys.argv[1:]:
        tag, b = arg.split(':')
        nodal(tag, b, results, Kf)
