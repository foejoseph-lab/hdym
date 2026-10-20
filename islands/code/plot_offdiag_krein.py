"""plot_offdiag_krein.py -- figure for B005: the Krein-signed angles psi+ of the two ghost-sector
pairs of the half-period map N (M = N^2) on the DR orbit, block 12, tracked by continuity.

Top: psi+/pi (angle of the positive-Krein member of each pair) unwrapped vs E; a pair off the
circle is drawn hollow at 0 or +-1 (real, sign undefined) or, in a Krein band, at +-psi with both
pairs hollow.  Middle: (psi1+ + psi2+)/pi mod 2 -- a true (dangerous) collision needs this = 0;
(psi1+ - psi2+)/pi mod 2 -- a safe pass-through needs this = 0.  Bottom: |mu_N| of each pair.
Reads results/dr_offdiag_krein.npy.  Usage (from islands/): python code/plot_offdiag_krein.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
r = np.load(os.path.join(HERE, '..', 'results', 'dr_offdiag_krein.npy'))
E = r[:, 0]
BLOCK = sys.argv[1] if len(sys.argv) > 1 else '12'
k = {'12': 0, '13': 1, '23': 2}[BLOCK]
o = 2 + 12 * k
pairs = np.stack([r[:, [o + 6, o + 7, o + 8]], r[:, [o + 9, o + 10, o + 11]]], axis=1)   # (nE, 2, [m, psi, def])

# ---- track by continuity of the angle (mod 2pi), unwrap
n = len(E)
tr = np.empty((n, 2, 3))
tr[0] = pairs[0]
for i in range(1, n):
    a, b = pairs[i]
    prev = tr[i - 1]
    def dist(p, q):
        return abs((p[1] - q[1] + np.pi) % (2 * np.pi) - np.pi)
    keep = dist(a, prev[0]) + dist(b, prev[1])
    swap = dist(b, prev[0]) + dist(a, prev[1])
    tr[i] = pairs[i] if keep <= swap else pairs[i][::-1]
psi = np.unwrap(tr[:, :, 1], axis=0)
m = tr[:, :, 0]
defined = tr[:, :, 2] > 0.5
krein = (~defined[:, 0]) & (~defined[:, 1]) & (np.abs(tr[:, 0, 1] + tr[:, 1, 1]) < 1e-6) & (m[:, 0] > 1 + 1e-6)

SURF = '#fcfcfb'
fig, axs = plt.subplots(3, 1, figsize=(11, 10), sharex=True, facecolor=SURF, gridspec_kw=dict(hspace=0.12, height_ratios=[2, 1.4, 1]))
cols = ('#2a78d6', '#eb6834')
for ax in axs:
    ax.set_facecolor(SURF)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.grid(True, color='#eeede9', lw=0.7)
    ax.set_axisbelow(True)
    ax.fill_between(E, -10, 10, where=krein, color='#e34948', alpha=0.12, lw=0)
ax = axs[0]
for j, c in enumerate(cols):
    on = defined[:, j]
    ax.scatter(E[on], psi[on, j] / np.pi, s=9, color=c, lw=0)
    ax.scatter(E[~on], psi[~on, j] / np.pi, s=12 + 400 * np.abs(np.log(m[~on, j])), facecolor='none', edgecolor=c, lw=0.9)
ax.set_ylabel('ψ⁺ / π   (angle of the positive-Krein member of each pair, N)')
ax.set_ylim(np.floor(psi.min() / np.pi) - 0.1, np.ceil(psi.max() / np.pi) + 0.1)
for yv in np.arange(np.floor(psi.min() / np.pi), np.ceil(psi.max() / np.pi) + 1):
    ax.axhline(yv, color='#52514e', lw=0.5, ls=':' if yv % 2 else '-', alpha=0.6)
ax.set_title(f'Block {BLOCK}, half-period map N (M = N²): Krein-signed angles of the two ghost-sector pairs vs E.  '
             'Solid lines: N at +1 (M at +1); dotted: N at −1 (M at +1); ±½ between them: N at ±i (M at −1).', fontsize=9.5, loc='left')
ax = axs[1]
ssum = ((psi[:, 0] + psi[:, 1]) / np.pi + 1) % 2 - 1
sdif = ((psi[:, 0] - psi[:, 1]) / np.pi + 1) % 2 - 1
ax.scatter(E, ssum, s=8, color='#e34948', lw=0, label='(ψ₁⁺ + ψ₂⁺)/π mod 2   — zero = true collision (opposite signs meet)')
ax.scatter(E, sdif, s=8, color='#1baf7a', lw=0, label='(ψ₁⁺ − ψ₂⁺)/π mod 2   — zero = safe meeting (same signs meet)')
ax.axhline(0, color='#52514e', lw=0.6)
ax.axhline(1, color='#52514e', lw=0.6, ls=':'); ax.axhline(-1, color='#52514e', lw=0.6, ls=':')
ax.set_ylim(-1.05, 1.05); ax.legend(frameon=False, fontsize=8, loc='lower left')
ax.set_ylabel('sum / difference')
ax = axs[2]
for j, c in enumerate(cols):
    ax.plot(E, m[:, j], color=c, lw=1.2)
ax.set_ylabel('|μ_N|  (|μ_M| = square)'); ax.set_xlabel('E'); ax.set_ylim(0.95, 1.05 * m.max())
ax.set_xscale('log')
ax.set_xticks([0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1]); ax.set_xticklabels(['0.05', '0.07', '0.1', '0.15', '0.2', '0.3', '0.4', '0.5', '0.7', '1'])
fig.text(0.01, 0.005, 'Red shading: Krein quadruplet (both pairs off the circle at ±ψ). Hollow markers: pair off the circle, size ∝ |ln|μ||; '
         'inside a tongue or band the sign is undefined and the label is carried by continuity of the angle.', fontsize=8, color='#52514e')
out = os.path.join(HERE, '..', 'results', f'offdiag_krein_b{BLOCK}.png')
fig.savefig(out, dpi=150, bbox_inches='tight')
print(out)
# ---- text summary of the meetings
print("sum ~ 0 (true collisions) at E:", np.round(E[np.abs(ssum) < 0.02], 4))
print("diff ~ 0 (safe meetings) at E:", np.round(E[np.abs(sdif) < 0.02], 4))
print("sum ~ +-1 (fake M-degeneracies, N antipodal) at E:", np.round(E[np.abs(np.abs(ssum) - 1) < 0.02], 4))
