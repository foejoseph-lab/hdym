"""rh4_lowE_long.py -- the low-E loiter at tmax = 12000 (islands/results/esc_E_lowE_long.npz, 22 Oct day run: 20 energies
in [0.05, 0.13] x 64 seeds, crossing times).  Loiter L = t_esc - t(|w| = 1); Kaplan-Meier survival with censoring at
tmax; exponential time constant tau = sum L / #escaped (censored MLE).
Usage: python rabbit_holes/ym/code/rh4_lowE_long.py > rabbit_holes/ym/results/rh4_lowE_long.txt
Added 2026-10-22.
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))


def km(L, c):
    order = np.argsort(L)
    Ls, cs = L[order], c[order]
    S, out, at_risk = 1.0, [], len(L)
    for j in range(len(Ls)):
        if not cs[j]:
            S *= 1 - 1 / at_risk
        at_risk -= 1
        out.append((Ls[j], S))
    out = np.array(out)
    return lambda x: out[out[:, 0] <= x, 1][-1] if (out[:, 0] <= x).any() else 1.0


def main():
    d = np.load(os.path.join(ROOT, 'islands', 'results', 'esc_E_lowE_long.npz'))
    E, T, te, st, tc, tmax = d['E'], d['T'], d['tesc'], d['status'], d['tcross'], float(d['tmax'])
    print("rh4_lowE_long.py: esc_E_lowE_long.npz, tmax 12000")
    print(" E      T     reached|w|=1 escaped cens | t_cross(1) med | loiter of escaped: median p25 p75 p90 [med/T] | censored loiters >= min/median (periods) | tau (periods)")
    for i in range(len(E)):
        t1 = tc[i, :, 0]
        r, e = t1 > 0, st[i] >= 0
        if r.sum() == 0:
            print(f"{E[i]:.4f} {T[i]:6.2f}    0")
            continue
        L = (np.where(e, te[i], tmax) - t1)[r]
        c = ~e[r]
        unc, cl = L[~c], L[c]
        q = np.percentile(unc, [50, 25, 75, 90]) if len(unc) else [np.nan] * 4
        tau = L.sum() / max((~c).sum(), 1) / T[i] if (~c).sum() else np.nan
        print(f"{E[i]:.4f} {T[i]:6.2f}   {r.sum():2d}       {e[r].sum():2d}     {c.sum():2d} | {np.median(t1[r]):7.0f} | {q[0]:7.0f} {q[1]:6.0f} {q[2]:6.0f} {q[3]:6.0f} [{q[0] / T[i]:5.1f}] | {(cl.min() / T[i] if len(cl) else np.nan):5.1f} {(np.median(cl) / T[i] if len(cl) else np.nan):5.1f} | {tau:6.1f}")
    print("\nsurvival of the loiter at the three off-diagonal-led energies (Kaplan-Meier vs exponential with the censored-MLE tau):")
    for Et in (0.0581, 0.1063, 0.1118):
        i = int(np.argmin(abs(E - Et)))
        t1 = tc[i, :, 0]
        r, e = t1 > 0, st[i] >= 0
        L = (np.where(e, te[i], tmax) - t1)[r]
        c = ~e[r]
        tau = L.sum() / max((~c).sum(), 1)
        S = km(L, c)
        xs = [500, 1000, 2000, 4000, 6000, 8000]
        print(f"  E={E[i]:.4f}: n={r.sum()}, censored {c.sum()}, tau = {tau:.0f} = {tau / T[i]:.0f} periods; S(L) at L = {xs}: KM {[round(S(x), 2) for x in xs]}  exp {[round(float(np.exp(-x / tau)), 2) for x in xs]}")


if __name__ == '__main__':
    main()
