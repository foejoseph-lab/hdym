"""rh1_phase.py -- rabbit hole 1: the phase-dependent loiter at a Krein band, redone from the committed orient maps,
with the control the list asserts but never ran (a parametric band, where the phase is only a sign).

STATUS (notes/RH1_phase.md; README hole 1): CLOSED, NEGATIVE.  On the maps this script finds a 4-7 % swing at ~4 sigma
against the tile-permuting null (a cell-wise null overstated it), nothing for the parametric map or a placebo phase --
a signature on one map, not a measurement.  The ensemble built to test it (12 phases x 256 independent cells,
rh_overnight_2026-10-22.py) shows no phase effect once |c| is controlled: residual z = +0.7 to +1.3; the raw 96 % vs
90 % escape by phase was |c| varying 0.56 -> 1.21 with alpha (the section's seed normalisation).  Do not cite a Krein
seed-phase effect.

Input: islands/results/esc_orient_0156.npz (E = 0.156, block 12 Krein quadruplet, 128^2 cells, tmax 3000) and
esc_orient_030.npz (E = 0.30, block 23 real pair).  No integration beyond one monodromy per map (~2 s).

For a Krein quadruplet the unstable subspace is spanned by the eigenvectors of mu and conj(mu); the seed's
projection on it is one complex number c(s) = w . s (w the left eigenvector of mu).  |c| sets the linear phase
(the nodal structure of offdiag_escape_nodal.py); arg c is the PHASE of the growing mode -- it decides at which
point of the orbit the mode reaches order one -- and the seed sets it freely.  For a real pair c is real and only
its sign is free, and s -> -s is a symmetry of the linearised problem, so no phase dependence can exist at linear
order.  The hole's claim: at E = 0.156 the non-escape fraction and the median escape time depend on arg c.

What is computed, for each map:
  * c(s) per cell; log t_esc regressed on log |c| (the linear-phase predictor), residual r = log t_esc - fit;
  * binning of the escape fraction, the median t_esc and the median residual by arg c (12 bins over [0, 2 pi) and
    6 over [0, pi), the latter because s -> -s maps arg c -> arg c + pi);
  * a permutation null: 200 random re-labellings of whole 8 x 8 tiles of cells (see binned()) give the spread of each binned statistic under
    "no phase effect"; the observed range across bins is compared with it (a z-score on the across-bin range);
  * for the control map (real pair) the same with "phase" = sign c (2 bins), and with the phase of the projection
    on the complex STABLE pair of the same block as a placebo (it should do nothing).
Non-escaped cells are those with status = -1 (tesc = tmax, censored).

Usage (from the repo root): python rabbit_holes/ym/code/rh1_phase.py > rabbit_holes/ym/results/rh1_phase.txt
Writes results/rh1_phase.png.  Added 2026-10-21 (rabbit-hole session).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as F                                    # noqa: E402
from offdiag_escape_nodal import seed_to_block                    # noqa: E402

RES = os.path.join(HERE, '..', 'results')
ISL = os.path.join(ROOT, 'islands', 'results')


def load(tag, b, Kf):
    d = np.load(os.path.join(ISL, tag + '.npz'), allow_pickle=True)
    te, st, th, ph, basis, E, tmax = d['tesc'], d['status'], d['theta'], d['phi'], d['basis'], float(d['E']), float(d['tmax'])
    T, clo, Ms, _ = F.monodromies(Kf, E, 1.0)
    mu, V = np.linalg.eig(Ms[b])
    W = np.linalg.inv(V)
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    Eb = [seed_to_block(basis[i], b) for i in range(3)]

    def proj(k):
        c = [W[k] @ e for e in Eb]
        return np.cos(TH) * c[0] + np.sin(TH) * np.cos(PH) * c[1] + np.sin(TH) * np.sin(PH) * c[2]
    return dict(te=te, st=st, E=E, T=T, tmax=tmax, mu=mu, proj=proj, n=te.shape[0])


def tile_perm(n, tile, rng):
    """a permutation of the n x n cells that moves whole tile x tile blocks (keeps the map's local correlation)"""
    nt = n // tile
    order = rng.permutation(nt * nt)
    idx = np.arange(n * n).reshape(n, n)
    out = np.empty_like(idx)
    for a in range(nt * nt):
        i, j = divmod(a, nt)
        bi, bj = divmod(order[a], nt)
        out[i * tile:(i + 1) * tile, j * tile:(j + 1) * tile] = idx[bi * tile:(bi + 1) * tile, bj * tile:(bj + 1) * tile]
    return out.ravel()


def binned(phase, nb, esc, logt, resid, rng, nperm=200, n=128, tile=8):
    """per-bin escape fraction, median log t (escaped only), median residual; and permutation nulls on the
    across-bin RANGE of each.  The null permutes 8 x 8 TILES of the map, not cells, so that the spatial
    correlation of the loiter texture (lag-1 corr 0.09-0.28) is kept: phase bins are wedges of the sphere and a
    cell-wise null would understate the range a correlated texture produces by chance."""
    edges = np.linspace(phase.min(), phase.max() + 1e-12, nb + 1)
    idx = np.clip(np.digitize(phase, edges) - 1, 0, nb - 1)

    def stats(esc, logt, resid):
        f = np.array([esc[idx == i].mean() for i in range(nb)])
        m = np.array([np.median(logt[(idx == i) & esc]) if (esc & (idx == i)).sum() > 5 else np.nan for i in range(nb)])
        r = np.array([np.median(resid[(idx == i) & esc]) if (esc & (idx == i)).sum() > 5 else np.nan for i in range(nb)])
        return f, m, r
    f, m, r = stats(esc, logt, resid)
    obs = np.array([f.max() - f.min(), np.nanmax(m) - np.nanmin(m), np.nanmax(r) - np.nanmin(r)])
    null = []
    for _ in range(nperm):
        p = tile_perm(n, tile, rng)
        fp, mp, rp = stats(esc[p], logt[p], resid[p])
        null.append([fp.max() - fp.min(), np.nanmax(mp) - np.nanmin(mp), np.nanmax(rp) - np.nanmin(rp)])
    null = np.array(null)
    z = (obs - null.mean(0)) / null.std(0)
    return edges, f, m, r, obs, null.mean(0), null.std(0), z


def analyse(name, m, phase_fn, nb_list, rng, placebo=None):
    te, st = m['te'].ravel(), m['st'].ravel()
    esc = st >= 0
    logt = np.where(esc, np.log10(np.where(te > 0, te, 1.0)), np.nan)
    c = m['c'].ravel()
    amp = np.abs(c)
    # linear-phase predictor: log t vs log |c| on escaped cells
    A = np.vstack([np.log10(amp[esc]), np.ones(esc.sum())]).T
    coef = np.linalg.lstsq(A, logt[esc], rcond=None)[0]
    pred = coef[0] * np.log10(amp) + coef[1]
    resid = logt - pred
    rr = np.corrcoef(np.log10(amp[esc]), logt[esc])[0, 1]
    print(f"\n== {name}: E = {m['E']}, T = {m['T']:.3f}; {esc.sum()} of {len(esc)} cells escape by tmax ({100 * esc.mean():.1f} %), median t_esc {np.median(te[esc]):.0f} ==")
    print(f"  linear predictor: corr(log t, log |c|) = {rr:+.3f}, slope {coef[0]:+.3f} dex/dex; residual sd {np.nanstd(resid[esc]):.3f} dex")
    out = {}
    for label, nb, ph in phase_fn(c):
        edges, f, md, rs, obs, nm, ns, z = binned(ph, nb, esc, logt, resid, rng)
        print(f"  -- {label}, {nb} bins --")
        print(f"     escape fraction per bin: {np.round(100 * f, 1)} %   range {100 * obs[0]:.1f} pts vs null {100 * nm[0]:.1f} +- {100 * ns[0]:.1f}  (z = {z[0]:+.1f})")
        print(f"     median t_esc per bin:    {np.round(10 ** md, 0)}   range in log10: {obs[1]:.3f} vs null {nm[1]:.3f} +- {ns[1]:.3f}  (z = {z[1]:+.1f})")
        print(f"     median residual per bin (after |c|): {np.round(rs, 3)}   range {obs[2]:.3f} vs null {nm[2]:.3f} +- {ns[2]:.3f}  (z = {z[2]:+.1f})")
        out[label] = (edges, f, md, rs, z)
    return out


def main():
    rng = np.random.default_rng(1)
    Kf = F.load_K()
    print("rh1_phase.py: phase of the seed's projection on the unstable Floquet subspace vs escape, from the committed orient maps")

    # Krein map
    m = load('esc_orient_0156', '12', Kf)
    unst = np.where(np.abs(m['mu']) > 1 + 1e-6)[0]
    k = unst[np.argmax(m['mu'][unst].imag)]                         # the member with Im mu > 0
    m['k'] = k
    m['c'] = m['proj'](k)
    print(f"Krein map: unstable multipliers {np.round(m['mu'][unst], 4)}; using the eigenvalue {np.round(m['mu'][k], 4)} (its conjugate carries the same information)")

    def krein_phases(c):
        a = np.angle(c)
        return [('arg c over [0, 2pi)', 12, np.mod(a, 2 * np.pi)), ('arg c mod pi (s -> -s identified)', 6, np.mod(a, np.pi))]
    outK = analyse('Krein quadruplet, block 12', m, krein_phases, None, rng)

    # control: parametric map
    mc = load('esc_orient_030', '23', Kf)
    unst = np.where(np.abs(mc['mu']) > 1 + 1e-6)[0]
    kc = unst[0]
    mc['k'] = kc
    cc = mc['proj'](kc)
    print(f"\ncontrol map: unstable multipliers {np.round(mc['mu'][unst], 4)}; projection is real to {np.abs(cc.imag).max() / np.abs(cc).max():.1e}")
    mc['c'] = cc.real + 0j
    # placebo phase: the projection on a stable complex pair of the same block (farthest from 1 among |mu| = 1)
    stab = [i for i in range(len(mc['mu'])) if abs(abs(mc['mu'][i]) - 1) < 1e-6 and mc['mu'][i].imag > 1e-3]
    kp = max(stab, key=lambda i: abs(mc['mu'][i] - 1))
    cp = mc['proj'](kp)
    print(f"  placebo: phase of the projection on the stable pair {np.round(mc['mu'][kp], 4)}")

    def para_phases(c):
        return [('sign c', 2, np.where(c.real >= 0, 1.0, 0.0)), ('placebo: arg of the projection on a stable pair, mod pi', 6, np.mod(np.angle(cp.ravel()), np.pi))]
    outP = analyse('parametric real pair, block 23 (control)', mc, para_phases, None, rng)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(14, 4))
        e, f, md, rs, z = outK['arg c over [0, 2pi)']
        x = 0.5 * (e[1:] + e[:-1]) / np.pi
        ax[0].bar(x, 100 * (1 - f), width=0.95 * (e[1] - e[0]) / np.pi, color='tab:red', alpha=0.7)
        ax[0].set_xlabel('arg c / pi (Krein, E = 0.156)')
        ax[0].set_ylabel('non-escape by 3000 (%)')
        ax[0].set_title(f'non-escape fraction vs phase (z = {z[0]:+.1f})', fontsize=9)
        ax[1].plot(x, 10 ** md, 'o-', color='tab:red', label='median t_esc')
        ax[1].set_xlabel('arg c / pi')
        ax[1].set_ylabel('median t_esc')
        ax[1].set_title(f'median escape time vs phase (z = {z[1]:+.1f})', fontsize=9)
        ax[2].plot(x, rs, 'o-', color='tab:red', label='Krein map')
        e2, f2, md2, rs2, z2 = outP['placebo: arg of the projection on a stable pair, mod pi']
        ax[2].plot(0.5 * (e2[1:] + e2[:-1]) / np.pi, rs2, 's--', color='tab:blue', label='control map, placebo phase')
        ax[2].axhline(0, color='k', lw=0.5)
        ax[2].set_xlabel('phase / pi')
        ax[2].set_ylabel('median residual of log10 t_esc after |c|')
        ax[2].set_title(f'phase effect after the linear predictor (z = {z[2]:+.1f} vs {z2[2]:+.1f})', fontsize=9)
        ax[2].legend(fontsize=7)
        fig.suptitle('rabbit hole 1: phase-dependent escape at a Krein band, from the committed orient maps', fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, 'rh1_phase.png'), dpi=130)
        print("figure: results/rh1_phase.png")
    except Exception as ex:                                       # noqa: BLE001
        print(f"(no figure: {ex})")


if __name__ == '__main__':
    main()
