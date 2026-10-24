"""rh4_loiter.py -- rabbit hole 4: the statistics of the nonlinear loiter, from the section-E checkpoints.

STATUS (notes/RH4_loiter.md; README hole 4): the checkpoint reading here (50-unit resolution) is SUPERSEDED by the
kernel's crossing times (rh_overnight_2026-10-22.py, rh4_lowE_long.py).  The answer to question 3 below is: the loiter
scales with the period (pooled median 0.58 periods = 18 time units; Spearman of the per-energy median with E: -0.05 in
periods, -0.46 in time units).  "E-independent in absolute time" (the hole as posed, from these checkpoints) is
WITHDRAWN -- the opposite holds.  Tail: stretched (p90/median 4.3-5.7 vs 3.3 exponential), as this script also found.

Input: islands/results/esc_E.npz (100 energies x 64 seeds; chk[E, seed, k, :] = four sector norms, diagonal-ghost
amplitude |w| and the Ostrogradsky energy at t_k = 50 (k+1), NaN after escape).  No integration.

Definitions.  t_1 = the time the diagonal ghost first reaches amplitude 1, interpolated log-linearly between the
last checkpoint below 1 and the first at or above it (exponential growth is the linear phase; the interpolation
error is a fraction of the 50-unit checkpoint spacing, estimated from the local growth rate).  The loiter is
L = t_esc - t_1, in time units and in orbit periods T(E).  A trajectory that reaches |w| >= 1 and does not escape by
tmax = 3000 is a right-censored loiter (L >= tmax - t_1); a trajectory that never reaches 1 contributes nothing.

Questions answered (results/rh4_loiter.txt, figure results/rh4_loiter.png):
  1. per energy: number of loiters, censoring, quartiles of L (time units and periods), the linear clock t_1
     against ln(1/delta)/rate_diag where the diagonal ghost is Floquet-unstable on the orbit.
  2. the loiter distribution pooled over E (Kaplan-Meier survival S(L) in periods, with censoring): is it
     exponential (memoryless escape from a repeller: S = exp(-L/tau)) or does it carry a power-law tail (KAM
     stickiness: S ~ L^-gamma at large L)?  Both fits are reported with their log-likelihoods; the tail exponent
     is the Hill estimator over the top decade.
  3. does L scale with T(E) (loiter = a number of orbits) or is it E-independent in absolute time?  The rank
     correlation of the median loiter with E in each unit, over the energies with >= 20 loiters.
  4. which sector is largest at escape against the loiter length (does a long loiter end differently?).

Usage (from the repo root): python rabbit_holes/ym/code/rh4_loiter.py > rabbit_holes/ym/results/rh4_loiter.txt
Added 2026-10-21 (rabbit-hole session).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
RES = os.path.join(HERE, '..', 'results')
NPZ = os.path.join(ROOT, 'islands', 'results', 'esc_E.npz')
SECTOR_NAMES = ['diagonal', 'block 12', 'block 13', 'block 23']


def first_crossing(times, amp, level=1.0, escaped=False, seed=1e-4):
    """log-linear interpolated time at which amp first reaches level; None if never.  Returns (t1, k) with k the
    index of the first checkpoint at or above level.  If the trajectory escaped with every recorded amplitude
    below level (the crossing and the blow-up fell inside one checkpoint interval), the crossing is EXTRAPOLATED
    from the last two recorded amplitudes with their local exponential rate (k = -1 marks it); with fewer than
    two finite checkpoints the rate from the seed is used."""
    ok = np.isfinite(amp)
    above = np.where(ok & (amp >= level))[0]
    if len(above) == 0:
        if not escaped:
            return None, None
        fin = np.where(ok)[0]
        if len(fin) == 0:
            return None, None                                       # escaped before the first checkpoint: no growth curve
        kl = fin[-1]
        if len(fin) >= 2 and amp[fin[-2]] > 0 and amp[kl] > amp[fin[-2]]:
            r = np.log(amp[kl] / amp[fin[-2]]) / (times[kl] - times[fin[-2]])
        elif amp[kl] > seed:
            r = np.log(amp[kl] / seed) / times[kl]
        else:
            return None, None
        return times[kl] + np.log(level / amp[kl]) / r, -1
    k = above[0]
    if k == 0:
        return times[0] * np.log(level / 1e-4) / np.log(amp[0] / 1e-4), 0     # from the seed amplitude
    a0, a1 = amp[k - 1], amp[k]
    if not np.isfinite(a0) or a0 <= 0:
        return times[k], k
    f = np.log(level / a0) / np.log(a1 / a0)
    return times[k - 1] + f * (times[k] - times[k - 1]), k


def kaplan_meier(L, cens):
    """survival S(t) at the sorted event times; cens[i] = True if the i-th value is a lower bound."""
    order = np.argsort(L)
    L, cens = L[order], cens[order]
    n = len(L)
    S = 1.0
    ts, Ss = [0.0], [1.0]
    at_risk = n
    i = 0
    while i < n:
        t = L[i]
        d = c = 0
        while i < n and L[i] == t:
            if cens[i]:
                c += 1
            else:
                d += 1
            i += 1
        if d > 0:
            S *= 1 - d / at_risk
            ts.append(t)
            Ss.append(S)
        at_risk -= d + c
    return np.array(ts), np.array(Ss)


def fit_exponential(L, cens):
    """MLE of tau for exp(-L/tau) with right censoring; log-likelihood."""
    d = (~cens).sum()
    tau = L.sum() / max(d, 1)
    ll = -d * np.log(tau) - L.sum() / tau
    return tau, ll


def fit_exp_power(L, cens, grid_t0=None):
    """S(L) = exp(-L/tau) for L < t0, continuous into (L/t0)^-g exp(-t0/tau) beyond (a crude but standard two-piece
    model: Poisson escape then a sticky tail).  Profile the break t0 on a grid, MLE in tau and g at each t0."""
    L = np.asarray(L)
    if grid_t0 is None:
        grid_t0 = np.quantile(L[~cens], np.linspace(0.3, 0.95, 14))
    best = None
    for t0 in grid_t0:
        head = L < t0
        tail = ~head
        # tau from the head (events below t0 plus everything beyond t0 as censored at t0)
        d_h = (head & ~cens).sum()
        T_h = L[head].sum() + t0 * tail.sum()
        tau = T_h / max(d_h, 1)
        # g from the tail: density g/t0 (L/t0)^-(g+1) for events, survival (L/t0)^-g for censored
        ev = tail & ~cens
        ce = tail & cens
        if ev.sum() < 3:
            continue
        s = np.log(L[tail] / t0).sum()
        g = ev.sum() / max(s, 1e-12)
        ll = (-d_h * np.log(tau) - T_h / tau
              + ev.sum() * np.log(g / t0) - (g + 1) * np.log(L[ev] / t0).sum() - g * np.log(L[ce] / t0).sum())
        if best is None or ll > best[0]:
            best = (ll, t0, tau, g)
    return best


def fit_pareto(x, cens):
    """S(x) = (1 + x/x0)^-gamma (Lomax / shifted Pareto), MLE in gamma at each x0 on a grid, with right censoring."""
    x = np.asarray(x)
    best = None
    for x0 in np.exp(np.linspace(np.log(max(x[~cens].min(), 1e-3) + 1e-3), np.log(x.max() * 5), 60)):
        s = np.log1p(x / x0)
        d = (~cens).sum()
        g = d / max(s.sum(), 1e-12)
        ll = d * np.log(g / x0) - (g + 1) * s[~cens].sum() - g * s[cens].sum()
        if best is None or ll > best[2]:
            best = (x0, g, ll)
    return best


def hill(L, top=0.1):
    """Hill tail-index estimator over the top fraction of the (uncensored) sample: S ~ L^-gamma."""
    x = np.sort(L)[::-1]
    k = max(int(top * len(x)), 5)
    return 1.0 / np.mean(np.log(x[:k] / x[k]))


def main():
    d = np.load(NPZ)
    E, T, tesc, sector, status, chk, rate = d['E'], d['T'], d['tesc'], d['sector'], d['status'], d['chk'], d['rate_diag']
    nchk, tmax, delta = int(d['nchk']), float(d['tmax']), float(d['delta'])
    times = tmax * np.arange(1, nchk + 1) / nchk
    print(f"rh4_loiter.py: esc_E.npz, {len(E)} energies x {tesc.shape[1]} seeds, tmax {tmax}, checkpoints every {times[0]:.0f}, seed amplitude {delta}")
    rows = []
    allL, allLp, allcens, allE, allsec, allt1 = [], [], [], [], [], []
    n_clamped = 0
    print("\n== per energy ==")
    print("columns: E | T | escaped/64 | reached |w|>=1: n (censored; crossings extrapolated from the last checkpoint) | t_1 median [p10,p90] | ln(1/delta)/rate (diag Floquet) | loiter L: median [p25,p75] p90 time units | in periods | sector largest at escape")
    for i in range(len(E)):
        Ls, Lp, cs, secs, t1s = [], [], [], [], []
        n_extrap = 0
        for j in range(tesc.shape[1]):
            amp = chk[i, j, :, 4]
            t1, k = first_crossing(times, amp, escaped=status[i, j] >= 0)
            if t1 is None:
                continue
            if k == -1:
                n_extrap += 1
            if status[i, j] == -1:                      # no escape by tmax: censored
                L = tmax - t1
                c = True
            else:
                L = tesc[i, j] - t1
                c = False
            if L < 0:                                   # extrapolated crossing later than the escape: unresolved short loiter
                L = 0.0
                n_clamped += 1
            Ls.append(L)
            Lp.append(L / T[i])
            cs.append(c)
            secs.append(sector[i, j])
            t1s.append(t1)
        n = len(Ls)
        esc = int((status[i] >= 0).sum())
        if n == 0:
            print(f"E={E[i]:.5f} T={T[i]:7.3f} | {esc:2d}/64 | none reached |w| = 1")
            continue
        Ls, Lp, cs, secs, t1s = map(np.array, (Ls, Lp, cs, secs, t1s))
        allL += list(Ls)
        allLp += list(Lp)
        allcens += list(cs)
        allE += [E[i]] * n
        allsec += list(secs)
        allt1 += list(t1s)
        lin = np.log(1 / delta) / rate[i] if rate[i] > 1e-4 else np.nan
        unc = Ls[~cs]
        q = np.percentile(unc, [25, 50, 75, 90]) if len(unc) else [np.nan] * 4
        qp = np.percentile(Lp[~cs], [50]) if len(unc) else [np.nan]
        secn = {SECTOR_NAMES[s]: int((secs[~cs] == s).sum()) for s in range(4) if (secs[~cs] == s).sum() > 0}
        rows.append((E[i], T[i], n, int(cs.sum()), np.median(t1s), lin, q[1], q[0], q[2], q[3], qp[0], n_extrap))
        print(f"E={E[i]:.5f} T={T[i]:7.3f} | {esc:2d}/64 | {n:2d} ({int(cs.sum()):2d}; {n_extrap:2d} extrapolated) | t1 {np.median(t1s):7.1f} [{np.percentile(t1s, 10):6.1f},{np.percentile(t1s, 90):6.1f}] | {lin:7.1f} | "
              f"L {q[1]:6.1f} [{q[0]:5.1f},{q[2]:6.1f}] p90 {q[3]:6.1f} | {qp[0]:5.2f} T | {secn}")
    allL, allLp, allcens, allE, allsec, allt1 = map(np.array, (allL, allLp, allcens, allE, allsec, allt1))
    print(f"\nloiters pooled: {len(allL)} trajectories reached |w| = 1 ({int(allcens.sum())} censored at tmax); escaped {int((~allcens).sum())}")
    print(f"INFO  loiters clamped to zero (extrapolated crossing later than the escape): {int(n_clamped)} of {int((~allcens).sum())} escaped; "
          f"the checkpoint spacing (50) and the extrapolation set the resolution of the shortest loiters")

    # 3. absolute time vs periods
    from scipy.stats import spearmanr
    good = [r for r in rows if r[2] - r[3] >= 20]
    Eg = np.array([r[0] for r in good])
    medL = np.array([r[6] for r in good])
    medLp = np.array([r[10] for r in good])
    print("\n== 3. does the loiter scale with the period? (energies with >= 20 escaped loiters) ==")
    print(f"  {len(good)} energies; Spearman(E, median L in time units) = {spearmanr(Eg, medL)[0]:+.3f} (p {spearmanr(Eg, medL)[1]:.2g});"
          f" Spearman(E, median L in periods) = {spearmanr(Eg, medLp)[0]:+.3f} (p {spearmanr(Eg, medLp)[1]:.2g})")
    print(f"  median L over those energies: {np.median(medL):.1f} time units (range {medL.min():.1f}-{medL.max():.1f}); {np.median(medLp):.2f} periods (range {medLp.min():.2f}-{medLp.max():.2f})")

    # 2. pooled survival, per unit, in E bins -- on the RESOLVED part of the distribution only: the checkpoint
    #    spacing (50 time units, 0.9-1.9 periods) is the resolution of the shortest loiters (most escapes cross
    #    |w| = 1 and blow up inside one interval and are extrapolated or clamped), so the fits are made on the
    #    conditional distribution beyond L_min (= 50 time units, or 2 periods), where every loiter is resolved.
    print("\n== 2. the loiter distribution beyond the resolution limit: exponential vs power-law tail (Kaplan-Meier with censoring) ==")
    print("   (fits on the excess L - L_min conditional on L >= L_min; exponential S = exp(-x/tau), Pareto S = (1 + x/x0)^-gamma with x0 profiled;")
    print("    'frac' = fraction of the trajectories that reached |w| = 1 whose loiter exceeds L_min)")
    bins = [(0.05, 0.2), (0.2, 0.35), (0.35, 0.6), (0.6, 1.001)]
    km = {}
    for unit, Lu, Lmin in (('time units', allL, 50.0), ('periods', allLp, 2.0)):
        print(f"  -- loiter in {unit}, L_min = {Lmin} --")
        for lo, hi in bins + [(0.05, 1.001)]:
            m = (allE >= lo) & (allE < hi)
            tot = int(m.sum())
            m = m & (Lu >= Lmin)
            if (m & ~allcens).sum() < 30:
                continue
            x, c = Lu[m] - Lmin, allcens[m]
            tau, ll_e = fit_exponential(x, c)
            fp = fit_pareto(x, c)
            ts, Ss = kaplan_meier(Lu[m], c)
            km[(unit, lo, hi)] = (ts, Ss)
            tag = f"E in [{lo:.2f}, {hi:.2f})" if (lo > 0.05 or hi < 1.001) else "all E"
            med = np.median(x[~c]) + Lmin
            h = hill(x[~c] + Lmin, 0.2)
            x0, g, ll_p = fp
            dll = ll_p - ll_e
            print(f"  {tag:18s} n={len(x):4d} of {tot:4d} (frac {len(x) / tot:4.2f}; cens {int(c.sum()):3d}) median L {med:7.2f} | exponential tau={tau:7.2f} ll={ll_e:9.1f}"
                  f" | Pareto x0={x0:7.2f} gamma={g:5.2f} ll={ll_p:9.1f} (1 extra parameter; delta ll = {dll:+6.1f}, AIC prefers {'power law' if dll > 1 else 'exponential'})"
                  f" | Hill gamma (top 20 %) = {h:.2f} | p90/median of the excess = {np.percentile(x[~c], 90) / max(np.median(x[~c]), 1e-9):5.2f} (exp 3.32)")

    # 4. sector at escape vs loiter length
    print("\n== 4. the sector largest at escape against the loiter length (uncensored, pooled; L in periods) ==")
    unc = ~allcens
    for lo, hi, name in [(-1, 1e-9, 'unresolved (clamped to 0)'), (1e-9, 1.0, '(0, 1] period'), (1.0, 3.0, '(1, 3] periods'), (3.0, np.inf, '> 3 periods')]:
        m = unc & (allLp > lo) & (allLp <= hi)
        cnt = {SECTOR_NAMES[s]: int((allsec[m] == s).sum()) for s in range(4)}
        offd = 100.0 * (allsec[m] != 0).sum() / max(m.sum(), 1)
        print(f"  loiter {name:26s}: n={int(m.sum()):4d}  {cnt}  off-diagonal-led {offd:4.1f} %")

    np.savez(os.path.join(RES, 'rh4_loiter.npz'), L=allL, Lp=allLp, cens=allcens, E=allE, sector=allsec, t1=allt1,
             rows=np.array(rows))
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
        r = np.array(rows)
        ax[0].errorbar(r[:, 0], r[:, 6], yerr=[r[:, 6] - r[:, 7], r[:, 8] - r[:, 6]], fmt='o', ms=3, lw=0.8, label='median L [p25, p75], time units')
        ax[0].plot(r[:, 0], r[:, 9], 'r.', ms=4, label='p90')
        ax[0].plot(r[:, 0], r[:, 1], 'k--', lw=0.8, label='one period T(E)')
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        ax[0].set_xlabel('E')
        ax[0].set_ylabel('loiter after |w| = 1')
        ax[0].legend(fontsize=7)
        for (unit, lo, hi), (ts, Ss) in km.items():
            if unit != 'periods':
                continue
            lab = 'all E' if (lo <= 0.05 and hi > 1) else f'E in [{lo:.2f},{hi:.2f})'
            ax[1].step(ts, Ss, where='post', lw=1.2 if hi > 1 else 0.8, label=lab)
        ax[1].set_yscale('log')
        ax[1].set_xlabel('loiter in periods')
        ax[1].set_ylabel('survival S(L)  (Kaplan-Meier)')
        ax[1].legend(fontsize=7)
        ax[1].set_title('exponential = straight line', fontsize=8)
        for (unit, lo, hi), (ts, Ss) in km.items():
            if unit != 'periods':
                continue
            lab = 'all E' if (lo <= 0.05 and hi > 1) else f'E in [{lo:.2f},{hi:.2f})'
            ax[2].loglog(ts[1:], Ss[1:], drawstyle='steps-post', lw=1.2 if hi > 1 else 0.8, label=lab)
        ax[2].set_xlabel('loiter in periods')
        ax[2].set_ylabel('S(L)')
        ax[2].set_title('power law = straight line', fontsize=8)
        ax[2].legend(fontsize=7)
        fig.suptitle('rabbit hole 4: the nonlinear loiter between |w| = 1 and the blow-up (esc_E.npz, 6400 trajectories)', fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, 'rh4_loiter.png'), dpi=130)
        print("figure: results/rh4_loiter.png")
    except Exception as e:                                            # noqa: BLE001
        print(f"(no figure: {e})")


if __name__ == '__main__':
    main()
