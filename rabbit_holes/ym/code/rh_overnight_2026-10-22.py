"""rh_overnight_2026-10-22.py -- first reading of the three overnight runs of 21-22 Oct (holes 1, 4, 12).

Inputs (islands/results/, from Joseph's machine, all checks PASS):
  esc_phase_0156.npz     --section phase, E = 0.156, block 12: 12 phases x 256 cells, |c| ~ fixed, arg c swept; chk + tcross
  esc_orient_030_n64.npz --section orient, E = 0.30, 64^2 cells, chk kept, tcross
  esc_E_cross.npz        --section E, 20 energies x 64 seeds, tcross (first times |w| >= 1, sector norms >= 1, |w| >= 0.1)

Hole 1 (phase ensemble).  Cells are independent draws (random admixture and random diagonal seed per cell), so a plain
permutation null is valid here (no spatial correlation, unlike the orient maps).  The recorded c = w.s per cell has
|c| varying with alpha (0.56-1.2: the direction Re(e^{i alpha} v) has alpha-dependent norm before normalisation --
a design flaw of the section, recorded), so the analysis regresses log t_esc on log |c| and bins the RESIDUAL by the
actual arg c, exactly as rh1_phase.py did for the maps; the escape fraction is binned the same way.  Prints z-scores
on the across-bin range against 500 permutations, for arg c over [0, 2 pi) (12 bins) and mod pi (6 bins).
Hole 4 (crossing times).  Loiter L = t_esc - t_cross(|w| = 1), now resolved to the integrator step: per-energy
quartiles in time units and periods; the linear-phase check t_cross(1) - t_cross(0.1) against ln(10)/rate_diag where
the diagonal ghost is Floquet-unstable; survival of L in periods pooled by E range, Pareto vs exponential as in
rh4_loiter.py but without the clamping.
Hole 12.  chk holds sector NORMS, not directions, so "a quadruplet's runaway rotates while a real pair's runs along a
direction" is not testable from these files either; what is testable is which sector crosses 1 first (tcross) in the
Krein ensemble against the parametric map, and the ordering of the crossings.  Reported; the directional question
needed state snapshots -- added to the kernel later on 22 Oct and answered by rh12_direction.py (hole 12 closed:
not observable nonlinearly, the diagonal ghost has taken over; notes/RH12_direction.md).

Outcomes (README holes 1, 4): hole 1 closed negative (z = +0.7 to +1.3 after controlling |c|; notes/RH1_phase.md);
hole 4 a measurement -- the loiter scales with the period, pooled median 0.58 periods, stretched tail, the linear
phase the clock (notes/RH4_loiter.md); "E-independent in absolute time" withdrawn.

Usage (from the repo root): python rabbit_holes/ym/code/rh_overnight_2026-10-22.py > rabbit_holes/ym/results/rh_overnight_2026-10-22.txt
Added 2026-10-22.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
ISL = os.path.join(ROOT, 'islands', 'results')
RES = os.path.join(HERE, '..', 'results')
SECTOR_NAMES = ['diagonal', 'block 12', 'block 13', 'block 23']


def binned_z(phase, nb, esc, logt, resid, rng, nperm=500):
    edges = np.linspace(phase.min(), phase.max() + 1e-12, nb + 1)
    idx = np.clip(np.digitize(phase, edges) - 1, 0, nb - 1)

    def stats(esc, logt, resid):
        f = np.array([esc[idx == i].mean() for i in range(nb)])
        m = np.array([np.median(logt[(idx == i) & esc]) for i in range(nb)])
        r = np.array([np.median(resid[(idx == i) & esc]) for i in range(nb)])
        return f, m, r
    f, m, r = stats(esc, logt, resid)
    obs = np.array([f.max() - f.min(), m.max() - m.min(), r.max() - r.min()])
    null = []
    for _ in range(nperm):
        p = rng.permutation(len(esc))
        fp, mp, rp = stats(esc[p], logt[p], resid[p])
        null.append([fp.max() - fp.min(), mp.max() - mp.min(), rp.max() - rp.min()])
    null = np.array(null)
    return f, m, r, obs, null.mean(0), null.std(0), (obs - null.mean(0)) / null.std(0)


def fit_exp(x, cens):
    d = (~cens).sum()
    tau = x.sum() / max(d, 1)
    return tau, -d * np.log(tau) - x.sum() / tau


def fit_pareto(x, cens):
    best = None
    for x0 in np.exp(np.linspace(np.log(max(x[~cens].min(), 1e-3) + 1e-3), np.log(x.max() * 5), 60)):
        s = np.log1p(x / x0)
        d = (~cens).sum()
        g = d / max(s.sum(), 1e-12)
        ll = d * np.log(g / x0) - (g + 1) * s[~cens].sum() - g * s[cens].sum()
        if best is None or ll > best[2]:
            best = (x0, g, ll)
    return best


def main():
    rng = np.random.default_rng(1)
    print("rh_overnight_2026-10-22.py")

    # ---------------------------------------------------------------- hole 1
    d = np.load(os.path.join(ISL, 'esc_phase_0156.npz'), allow_pickle=True)
    te, st, c, al = d['tesc'], d['status'], d['cproj'], d['alpha']
    tmax = float(d['tmax'])
    esc = (st >= 0).ravel()
    logt = np.where(esc, np.log10(np.where(te.ravel() > 0, te.ravel(), 1.0)), np.nan)
    amp = np.abs(c).ravel()
    arg = np.angle(c).ravel()
    A = np.vstack([np.log10(amp[esc]), np.ones(esc.sum())]).T
    coef = np.linalg.lstsq(A, logt[esc], rcond=None)[0]
    resid = logt - (coef[0] * np.log10(amp) + coef[1])
    print(f"\n== hole 1: phase ensemble at E = {float(d['E'])}, block 12 Krein; {te.size} cells, {esc.sum()} escape by {tmax:.0f} ({100 * esc.mean():.1f} %), median t_esc {np.median(te.ravel()[esc]):.0f} ==")
    print(f"  |c| per cell: {amp.min():.2f}-{amp.max():.2f} (design flaw: varies with alpha); corr(log t, log |c|) = {np.corrcoef(np.log10(amp[esc]), logt[esc])[0, 1]:+.3f}, slope {coef[0]:+.3f}, residual sd {np.nanstd(resid[esc]):.3f} dex")
    print("  per nominal alpha: escape %, median t_esc, mean |c|:")
    for i in range(len(al)):
        e = st[i] >= 0
        print(f"    alpha/pi = {al[i] / np.pi:.3f}: {100 * e.mean():5.1f} %  {np.median(te[i][e]):7.1f}  |c| {np.abs(c[i]).mean():.2f}")
    for label, nb, ph in (('arg c over [0, 2pi)', 12, np.mod(arg, 2 * np.pi)), ('arg c mod pi', 6, np.mod(arg, np.pi))):
        f, m, r, obs, nm, ns, z = binned_z(ph, nb, esc, logt, resid, rng)
        print(f"  -- {label}, {nb} bins, permutation null (cells independent) --")
        print(f"     escape fraction by bin {np.round(100 * f, 1)} %: range {100 * obs[0]:.1f} pts vs null {100 * nm[0]:.1f} +- {100 * ns[0]:.1f}  z = {z[0]:+.1f}")
        print(f"     median t_esc by bin {np.round(10 ** m, 0)}: range {obs[1]:.3f} dex vs null {nm[1]:.3f} +- {ns[1]:.3f}  z = {z[1]:+.1f}")
        print(f"     median residual by bin {np.round(r, 3)}: range {obs[2]:.3f} vs null {nm[2]:.3f} +- {ns[2]:.3f}  z = {z[2]:+.1f}")

    # ---------------------------------------------------------------- hole 4
    print("\n== hole 4: loiter from the crossing times (section E, 20 energies x 64) ==")
    e = np.load(os.path.join(ISL, 'esc_E_cross.npz'), allow_pickle=True)
    E, T, teE, stE, tc, rate = e['E'], e['T'], e['tesc'], e['status'], e['tcross'], e['rate_diag']
    print("  E | T | escaped | reached |w|=1 (censored) | growth 0.1 -> 1: median t_cross(1)-t_cross(0.1) vs ln(10)/rate | loiter L = t_esc - t_cross(1): median [p25, p75] p90, in periods | sector first to cross 1 (counts)")
    allL, allLp, allc, allE = [], [], [], []
    for i in range(len(E)):
        t1, t01 = tc[i, :, 0], tc[i, :, 4]
        reached = t1 > 0
        escd = stE[i] >= 0
        n = int(reached.sum())
        if n == 0:
            print(f"  {E[i]:.4f} {T[i]:6.2f} {int(escd.sum()):2d}/64  none reached |w| = 1")
            continue
        g = (t1 - t01)[reached & (t01 > 0)]
        lin = np.log(10) / rate[i] if rate[i] > 1e-4 else np.nan
        L = np.where(escd, teE[i], tmax) - t1
        L = L[reached]
        cens = ~escd[reached]
        Lp = L / T[i]
        unc = L[~cens]
        q = np.percentile(unc, [25, 50, 75, 90]) if len(unc) else [np.nan] * 4
        # which sector crosses 1 first
        first = []
        for jcell in np.where(reached)[0]:
            tt = tc[i, jcell].copy()
            tt[4] = -1
            tt = np.where(tt > 0, tt, np.inf)
            first.append(int(np.argmin(tt)))
        cnt = {SECTOR_NAMES[s]: first.count(s) for s in range(4) if first.count(s)}
        allL += list(L)
        allLp += list(Lp)
        allc += list(cens)
        allE += [E[i]] * len(L)
        print(f"  {E[i]:.4f} {T[i]:6.2f} {int(escd.sum()):2d}/64  {n:2d} ({int(cens.sum()):2d})  growth {np.median(g) if len(g) else np.nan:6.1f} vs {lin:6.1f}  L {q[1]:6.1f} [{q[0]:5.1f},{q[2]:6.1f}] p90 {q[3]:6.1f}  {q[1] / T[i]:5.2f} T  {cnt}")
    allL, allLp, allc, allE = map(np.array, (allL, allLp, allc, allE))
    print(f"  pooled: {len(allL)} loiters ({int(allc.sum())} censored); uncensored median {np.median(allL[~allc]):.1f} time units, {np.median(allLp[~allc]):.2f} periods; min {allL[~allc].min():.1f}")
    from scipy.stats import spearmanr
    medE = [(E[i], np.median((np.where(stE[i] >= 0, teE[i], tmax) - tc[i, :, 0])[(tc[i, :, 0] > 0) & (stE[i] >= 0)])) for i in range(len(E)) if ((tc[i, :, 0] > 0) & (stE[i] >= 0)).sum() >= 20]
    if len(medE) > 4:
        Em = np.array([m[0] for m in medE])
        Lm = np.array([m[1] for m in medE])
        Tm = np.interp(Em, E, T)
        print(f"  scaling over {len(medE)} energies with >= 20 escaped loiters: Spearman(E, median L) = {spearmanr(Em, Lm)[0]:+.2f} (time units), {spearmanr(Em, Lm / Tm)[0]:+.2f} (periods)")
    print("  survival of L in periods, no clamping (fits on L > 0.02 periods):")
    for lo, hi in [(0.05, 0.2), (0.2, 0.35), (0.35, 0.6), (0.6, 1.001)]:
        m = (allE >= lo) & (allE < hi) & (allLp > 0.02)
        if (m & ~allc).sum() < 20:
            continue
        x, cc = allLp[m], allc[m]
        tau, lle = fit_exp(x, cc)
        x0, g, llp = fit_pareto(x, cc)
        print(f"    E in [{lo:.2f}, {hi:.2f}): n={m.sum():4d} (cens {int(cc.sum())}) median {np.median(x[~cc]):.2f} T; exponential tau={tau:.2f} ll={lle:.1f} | Pareto x0={x0:.2f} gamma={g:.2f} ll={llp:.1f} (delta ll {llp - lle:+.1f}); p90/median {np.percentile(x[~cc], 90) / np.median(x[~cc]):.2f} (exp 3.32)")

    # the Krein ensemble's own loiter, and the parametric map's
    print("\n== hole 4 bis: loiter at the two bands (phase ensemble = Krein, E = 0.156; n64 map = parametric, E = 0.30) ==")
    for name, tes, sts, tcs, Tper in (('Krein E=0.156', te, st, d['tcross'], float(d['T'])), ('parametric E=0.30', np.load(os.path.join(ISL, 'esc_orient_030_n64.npz'), allow_pickle=True)['tesc'], None, None, None)):
        if sts is None:
            o = np.load(os.path.join(ISL, 'esc_orient_030_n64.npz'), allow_pickle=True)
            tes, sts, tcs, Tper = o['tesc'], o['status'], o['tcross'], float(o['T'])
        t1 = tcs[..., 0].ravel()
        es = (sts >= 0).ravel()
        ok = (t1 > 0) & es
        L = tes.ravel()[ok] - t1[ok]
        print(f"  {name}: {ok.sum()} escaped cells with |w| = 1 crossed; loiter median {np.median(L):.1f} ({np.median(L) / Tper:.2f} T), p10 {np.percentile(L, 10):.1f}, p90 {np.percentile(L, 90):.1f} ({np.percentile(L, 90) / Tper:.2f} T); t_cross(1) median {np.median(t1[ok]):.0f}")
        # first sector to cross 1
        first = []
        for row in tcs.reshape(-1, 5)[ok]:
            tt = row.copy()
            tt[4] = -1
            tt = np.where(tt > 0, tt, np.inf)
            first.append(int(np.argmin(tt)))
        print(f"     first sector to reach norm 1: {{{', '.join(f'{SECTOR_NAMES[s]}: {first.count(s)}' for s in range(4) if first.count(s))}}}")

    print("\n== hole 12 ==")
    print("  chk stores sector norms only; the rotating-vs-directional question needs state snapshots and is not answerable from these files.")


if __name__ == '__main__':
    main()
