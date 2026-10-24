"""rh_krein_width_cj.py -- open item (iii) of RH_reduced_block.md: the width and strength of block 12's low Krein band
(E in [0.151, 0.163] at c_j = 1) across c_j = 0.25, 0.5, 0.75, 1, on a fine grid, as the one first-order test.

The two-mode resonance normal form (RH27_averaging.md, part B) makes a sum resonance of two signed-frequency modes
unstable in a band of half-width 2|g| in detuning and maximal rate |g|, with g the resonant matrix element of the
coupling; the coupling of the ghost to the background scales with c_j, so first-order theory predicts the band's
width in detuning and its peak rate both proportional to c_j (and its centre fixed by the unperturbed frequencies,
which depend on c_j through the frozen spectrum).  Measured here, per c_j: the band edges, the width in ln E and in
the ladder's detuning (orbit period, RH2_ladder.md: T_c = 4.057 (n + 0.210), so detuning in units of a rung is
dT/4.057), the peak |mu| and the peak rate ln|mu|/T, and the ratios against c_j = 1.  The committed coarse scan
rh27_cj_scan.txt (120 energies in [0.04, 1]) already has the band at every c_j except 0.5 on its 0.027-in-ln-E grid.

Usage (from the repo root): python rabbit_holes/ym/code/rh_krein_width_cj.py [--cj 0.25 0.5 0.75 1.0] [--Emin 0.12 --Emax 0.19 --nE 71]
Cost: one 8-dim block monodromy per (c_j, E), ~0.6 s each.  Writes results/rh_krein_width_cj.{txt,npz}.  Added 2026-10-24.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402
import rh27_mu_scan as MS                    # noqa: E402

RES = os.path.join(HERE, '..', 'results')
RUNG = 4.057


def bands(Es, types, rates):
    out, cur = [], None
    for E, ty, r in zip(Es, types, rates):
        if ty != 'stable':
            if cur is None or cur[2] != ty:
                if cur is not None:
                    out.append(cur)
                cur = [E, E, ty, r, E]
            else:
                cur[1] = E
                if r > cur[3]:
                    cur[3], cur[4] = r, E
        else:
            if cur is not None:
                out.append(cur); cur = None
    if cur is not None:
        out.append(cur)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cj', type=float, nargs='+', default=[0.25, 0.5, 0.75, 1.0])
    ap.add_argument('--Emin', type=float, default=0.12)
    ap.add_argument('--Emax', type=float, default=0.19)
    ap.add_argument('--nE', type=int, default=71)
    a = ap.parse_args()
    Kf = base.load_K()
    Es = np.exp(np.linspace(np.log(a.Emin), np.log(a.Emax), a.nE))
    print(f"rh_krein_width_cj.py: block 12 on the DR orbit, c_j in {a.cj}, {a.nE} energies in [{a.Emin}, {a.Emax}] (ln-E step {np.log(a.Emax / a.Emin) / (a.nE - 1):.4f})")
    res = {}
    for cj in a.cj:
        types, rates, Ts = [], [], []
        for E in Es:
            T, M = MS.monodromy_mu(Kf, E, cj, 'cj')
            ty, mu = MS.classify(M)
            types.append(ty); rates.append(np.log(mu) / T); Ts.append(T)
        res[cj] = (np.array(types), np.array(rates), np.array(Ts))
        print(f"\n== c_j = {cj} ==")
        for lo, hi, ty, r, Ep in bands(Es, types, rates):
            i0, i1 = np.searchsorted(Es, lo), np.searchsorted(Es, hi)
            print(f"  [{lo:.4f}, {hi:.4f}] {ty:14s} width {np.log(hi / lo):.4f} in ln E, {abs(Ts[i0] - Ts[i1]) / RUNG:.4f} rung in period; peak rate {r:.5f} per unit time (|mu| {np.exp(r * Ts[np.searchsorted(Es, Ep)]):.4f}) at E = {Ep:.4f}")
    # the low Krein band: the Krein band nearest E = 0.157 at each c_j
    print("\nThe low Krein band (nearest 0.157) across c_j, and the first-order ratios against c_j = 1 (width and peak rate both ~ c_j at first order):")
    print("  c_j  |  [lo, hi]           width/lnE  width/rung | peak rate   |mu|   centre | width ratio  rate ratio  (c_j ratio)")
    ref = None
    rows = {}
    for cj in a.cj:
        types, rates, Ts = res[cj]
        kb = [b for b in bands(Es, types, rates) if b[2] == 'Krein']
        if not kb:
            print(f"  {cj:4.2f} |  no Krein band in the window"); continue
        b = min(kb, key=lambda b: abs(np.sqrt(b[0] * b[1]) - 0.157))
        lo, hi, ty, r, Ep = b
        i0, i1 = np.searchsorted(Es, lo), np.searchsorted(Es, hi)
        w, wr, c = np.log(hi / lo), abs(Ts[i0] - Ts[i1]) / RUNG, np.sqrt(lo * hi)
        rows[cj] = (lo, hi, w, wr, r, c)
    if 1.0 in rows:
        ref = rows[1.0]
    for cj, (lo, hi, w, wr, r, c) in rows.items():
        rr = f"{w / ref[2]:8.2f}    {r / ref[4]:8.2f}    ({cj:.2f})" if ref else ''
        print(f"  {cj:4.2f} |  [{lo:.4f}, {hi:.4f}]   {w:.4f}     {wr:.4f}   | {r:.5f}   {np.exp(r * np.mean([t for t in res[cj][2]])):.3f}  {c:.4f} | {rr}")
    np.savez(os.path.join(RES, 'rh_krein_width_cj.npz'), E=Es, **{f'cj{cj:g}_types': res[cj][0] for cj in res}, **{f'cj{cj:g}_rates': res[cj][1] for cj in res}, **{f'cj{cj:g}_T': res[cj][2] for cj in res})


if __name__ == '__main__':
    main()
