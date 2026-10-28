"""rh27_mu_scan.py -- the reviewer's test 4, in two versions: --dial mu (INVALID, kept for the record: scaling the
component cross-coupling breaks the gauge structure and the instabilities it produces violate the linearised
Gauss law at the 0.9 level -- the F118 trap of B002) and --dial cj (valid: c_j is the ghost's coupling to the
background in the Lagrangian, the block at every c_j is the linearisation of a gauge-invariant theory; B001 did
c_j = 0 and 1).

Original description of the mu version (notes/RH27_averaging.md, item 7): scale the COMPONENT cross-coupling of
block 12 by mu and scan the Floquet bands in E.  (The "Krein bands are the odd rungs" prediction it was meant to test
was itself withdrawn the same evening, RH27 section C; the c_j scan's result -- c_j = 1 the least unstable member
except c_j = 0, c_j = 0.25 a double +1 fold with |mu| up to 34 -- is RH27 item 8, results/rh27_cj_scan.txt.)

In the (u12, u21) basis the half-period symmetry makes the self blocks of K even-harmonic and the cross blocks
odd-harmonic (rh7_selection_rule.py, step 1).  K_mu = K with the u12-u21 cross blocks multiplied by mu keeps the
symmetry for every mu (S~ acts within each component), so the Lemma holds at every mu: a single pair never folds
at -1.  At mu = 0 the two components are independent T/2-periodic systems, each with one ghost pair -- the only
instabilities possible are +1 folds of either component (the even rungs).  Prediction under "Krein bands are the
odd rungs dressed by the cross-coupling": at mu = 0 the Krein bands are closed and the parametric bands remain near
their places; as mu grows the Krein bands open at the odd rungs with width growing with mu, while the parametric
bands' widths change little.  Counter-prediction (the Krein bands are a combination resonance of two modes at
their own frequencies): their positions move with mu.

Usage (from the repo root): python rabbit_holes/ym/code/rh27_mu_scan.py [--dial mu|cj] [--mus ...] [--nE 150] [--tag ...]
(defaults: --dial cj --tag rh27_cj_scan --mus 0 0.3 0.6 1.  The committed logs, by their headers: rh27_mu_scan.txt
--dial mu --mus 0 0.3 0.6 1 --nE 120 --tag rh27_mu_scan; rh27_cj_scan.txt --dial cj --mus 0.25 0.5 0.75 1 --nE 120.)
Cost: one block monodromy per (mu, E), ~1.5 s each: four values of mu at 150 energies ~15 min.
Writes results/<tag>.txt (this print), .npz, .png.  (The final "figure: results/rh27_mu_scan.png" print is wrong for
any other tag: the file is <tag>.png.)  Added 2026-10-21 (late).
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

RES = os.path.join(HERE, '..', 'results')
B001 = os.path.join(ROOT, 'islands', 'results', 'dr_offdiag_floquet_typed.txt')
IDX = base.block_idx('12')
# block state order (u12, u21, u12', u21', ...): component a = even positions, b = odd positions
CA = [0, 2, 4, 6]
CB = [1, 3, 5, 7]


def monodromy_mu(Kf, E, mu, dial='mu'):
    q0, p0, T = base.orbit(E)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        if dial == 'cj':
            K = np.array(Kf(*q, *p, mu), dtype=float)[np.ix_(IDX, IDX)]
        else:
            K = np.array(Kf(*q, *p, 1.0), dtype=float)[np.ix_(IDX, IDX)]
            K[np.ix_(CA, CB)] *= mu
            K[np.ix_(CB, CA)] *= mu
        return np.concatenate([p, -base.gradV(q), (K @ s[6:].reshape(8, 8)).ravel()])
    sT = solve_ivp(rhs, (0.0, T), np.concatenate([q0, p0, np.eye(8).ravel()]), **base.RT).y[:, -1]
    return T, sT[6:].reshape(8, 8)


def classify(M):
    g = base.ghost_multipliers(M)
    off = np.abs(np.abs(g) - 1) > 1e-6
    if not off.any():
        return 'stable', 1.0
    if off.sum() == 4 and np.abs(g[off].imag).max() > 1e-6 * np.abs(g[off]).max():
        return 'Krein', np.abs(g).max()
    m = g[off][np.argmax(np.abs(g[off]))]
    tag = 'fold+1' if m.real > 0 else 'fold-1'
    if off.sum() == 4:
        tag = 'double ' + tag                                 # both pairs real at the same sign: two folds, not a collision
    return tag, np.abs(g).max()


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
            out.append((cur[0], cur[1], 'CPLX' if 'CPLX' in m.group(3) else 'REAL'))
    return out


def bands(Es, types):
    out = []
    i = 0
    while i < len(Es):
        if types[i] != 'stable':
            j = i
            while j + 1 < len(Es) and types[j + 1] == types[i]:
                j += 1
            out.append((Es[i], Es[j], types[i]))
            i = j + 1
        else:
            i += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mus', type=float, nargs='+', default=[0.0, 0.3, 0.6, 1.0])
    ap.add_argument('--nE', type=int, default=150)
    ap.add_argument('--emin', type=float, default=0.04)
    ap.add_argument('--emax', type=float, default=1.0)
    ap.add_argument('--dial', choices=['mu', 'cj'], default='cj')
    ap.add_argument('--tag', default='rh27_cj_scan')
    a = ap.parse_args()
    t0 = time.time()
    print(f"rh27_mu_scan.py: block 12, dial = {a.dial} ({'component cross-coupling scale, GAUGE-BREAKING' if a.dial == 'mu' else 'c_j, the ghost-background coupling of the Lagrangian'}) in {a.mus}; nE={a.nE} in [{a.emin}, {a.emax}]; RT={base.RT}")
    Kf = base.load_K()
    grid = np.logspace(np.log10(a.emin), np.log10(a.emax), a.nE)
    allres = {}
    for mu in a.mus:
        types, amps = [], []
        for E in grid:
            T, M = monodromy_mu(Kf, E, mu, a.dial)
            ty, am = classify(M)
            types.append(ty)
            amps.append(am)
        allres[mu] = (types, amps)
        print(f"\n== {a.dial} = {mu} ==")
        nminus = sum(1 for ty in types if ty == 'fold-1')        # a 'double fold-1' is two pairs: allowed by the Lemma
        print(f"  CHECK no single-pair fold at -1 (Lemma, any mu): {nminus} found {'PASS' if nminus == 0 else 'FAIL'}")
        for lo, hi, ty in bands(grid, types):
            print(f"  [{lo:.4f}, {hi:.4f}] {ty:7s} max |mu| {max(am for E_, am, t_ in zip(grid, amps, types) if lo <= E_ <= hi):.4f}  width {np.log(hi / lo):.3f} in log E")
    print("\n== measured block 12 (B001, mu = 1) ==")
    for lo, hi, ty in b001_intervals():
        print(f"  [{lo:.4f}, {hi:.4f}] {ty}")
    np.savez(os.path.join(RES, a.tag + '.npz'), E=grid, mus=np.array(a.mus),
             types=np.array([allres[m][0] for m in a.mus]), amps=np.array([allres[m][1] for m in a.mus]))
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(len(a.mus), 1, figsize=(9, 1.6 * len(a.mus) + 1), sharex=True)
        for k, mu in enumerate(a.mus):
            types, amps = allres[mu]
            col = np.where(np.array(types) == 'Krein', 'tab:red', np.where(np.array(types) == 'fold+1', 'tab:blue', 'lightgrey'))
            ax[k].scatter(grid, amps, c=col, s=8)
            ax[k].set_ylabel(f'{a.dial}={mu}\n|mu_max|', fontsize=8)
            for lo, hi, ty in b001_intervals():
                ax[k].axvspan(lo, hi, color=('tab:red' if ty == 'CPLX' else 'tab:blue'), alpha=0.1, lw=0)
        ax[-1].set_xscale('log')
        ax[-1].set_xlabel('E')
        ax[0].set_title(f'block 12, dial {a.dial} (red Krein, blue +1 fold); shading: measured bands at c_j = 1', fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, a.tag + '.png'), dpi=130)
        print("figure: results/rh27_mu_scan.png")
    except Exception as e:                                            # noqa: BLE001
        print(f"(no figure: {e})")
    print(f"DONE {time.time() - t0:.0f}")


if __name__ == '__main__':
    main()
