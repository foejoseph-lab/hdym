"""dr_offdiag_floquet_typed.py -- B001: the off-diagonal ghost on the DR orbit, with the TYPE of
every unstable multiplier.

Same system, orbit and monodromy as dr_offdiag_floquet.py (imported from it, so there is one rhs);
this file adds what the bridge needs (WORKORDER_bridge_2026-10-18_v2.md, B001):

  * a type column per block.  For the largest ghost multiplier off the unit circle:
      REAL+   mu real, mu > 1          parametric (period-1 tongue, multipliers leave the circle at +1)
      REAL-   mu real, mu < -1         parametric (period-doubling tongue, leave at -1)
      CPLX    |Im mu| > 1e-6 |mu|      Krein collision: a complex quadruplet (mu, conj mu, 1/mu, 1/conj mu);
                                       two unit-circle pairs of opposite Krein signature met and left
                                       the circle (the ghost supplies the negative signature)
    together with |mu| and arg(mu)/pi.
  * the count of ALL 8 multipliers off the circle (|ln|mu|| > OFF_TOL), not just the 4 that
    ghost_multipliers() keeps: wherever that count changes between consecutive grid points, all 8
    multipliers at both points are printed, so a collision at arg -> 0 or pi cannot be masked by the
    "drop the 4 nearest 1" rule.  A multiplier off the circle that is NOT among the ghost 4 is flagged
    MASKED (never seen so far; the flag is there so it cannot go unnoticed).
  * the linearised Gauss-law violation |G v| / |v| of the leading unstable eigenvector, G from the
    pickle (row matched to the block by its non-zero columns); a growing mode that violates Gauss is
    not physical (cf. the in-channel frozen-coefficient trap, F118).
  * the unstable intervals of the fine scan, each typed, and whether the type changes inside it.

Checks printed first (CHECK ... PASS/FAIL): background closure, det M = 1, (mu, 1/mu) pairing,
multiplier 1 present, and the c_j = 0 reduction to the scalar Hill equations.  Ends with DONE <s>.

Float only (DOP853, rtol 1e-12); nothing here is certified.

Usage (Spyder, from the repo root, kernel's own interpreter):
    import sys; %cd C:\\GitHub\\hdym
    !"{sys.executable}" -B -u islands\\code\\dr_offdiag_floquet_typed.py --ncores 20 > islands\\results\\dr_offdiag_floquet_typed.txt 2>&1
Options: --grid full (400 log-spaced points in [0.05, 1], the committed scan; default) | coarse (the 16
energies of the first table only); --ncores N (energies are independent; default 1); --emin/--emax/--npts
override the fine grid.

Run parameters and why: RT = DOP853 rtol 1e-12 atol 1e-13 (as dr_offdiag_floquet.py; its tolerance
check at E = 0.3 showed the ghost multipliers stable to 1e-6 between rtol 1e-11 and 1e-12).
UNST_TOL = 1e-6 on |mu| - 1 (as before: the verdict threshold).  OFF_TOL = 1e-4 on |ln|mu|| for the
all-8 count: the defective multiplier 1 splits by ~1e-5 in floating point and must not count as off
the circle; det M - 1 is ~1e-10.  IM_TOL = 1e-6 relative for real vs complex: an off-circle real pair
has Im mu at rounding level (~1e-12).
Added 2026-10-19.
"""
import argparse
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dr_offdiag_floquet as base                    # noqa: E402  (one rhs: orbit, monodromies, load_K)

UNST_TOL, OFF_TOL, IM_TOL = 1e-6, 1e-4, 1e-6
COARSE = [1.0, 0.75, 0.6, 0.5, 0.45, 0.4, 0.35, 0.315, 0.3, 0.28, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05]

_Kf = None
_Gf = None


def load_G():
    """The 3 x 24 linearised Gauss constraint, one row per block, as a numpy function; returns
    {block: (row function)} by matching each row's non-zero columns to a block's index set."""
    import pickle
    d = pickle.load(open(base.PKL, 'rb'))
    X, P, cj = d['X'], d['P'], d['cj']
    G = d['G']
    out = {}
    for r in range(G.rows):
        cols = {c for c in range(G.cols) if G[r, c] != 0}
        for b in base.BLOCKS:
            if cols <= set(base.block_idx(b)):
                out[b] = (r, sp.lambdify(list(X) + list(P) + [cj], G.row(r), 'numpy'))
    assert set(out) == set(base.BLOCKS), f"Gauss rows do not match blocks: {list(out)}"
    return out


def _init():
    global _Kf, _Gf
    _Kf = base.load_K()
    _Gf = load_G()


def classify(mu, tol=UNST_TOL):
    """mu: multipliers (any number).  Returns a list of (mu, type) for those with |mu| > 1 + tol,
    type in REAL+/REAL-/CPLX, sorted by |mu| descending.  Call it on the ghost 4 with UNST_TOL (the
    verdict, as in dr_offdiag_floquet.py); on all 8 only with OFF_TOL, or the split defective 1
    (1.00002) is counted as unstable."""
    out = []
    for m in mu:
        if abs(m) > 1.0 + tol:
            if abs(m.imag) > IM_TOL * abs(m):
                t = 'CPLX'
            else:
                t = 'REAL+' if m.real > 0 else 'REAL-'
            out.append((m, t))
    return sorted(out, key=lambda z: -abs(z[0]))


def gauss_violation(b, q, p, cj, v):
    r, gf = _Gf[b]
    row = np.array(gf(*q, *p, cj), dtype=complex).ravel()
    idx = base.block_idx(b)
    return abs(row[idx] @ v) / np.linalg.norm(v)


def one_energy(E, cj=1.0):
    """Everything about one energy: per block the 8 multipliers, the ghost 4, the typed unstable
    ones, the Gauss violation of the leading one, and the checks."""
    T, clo, Ms, Ys = base.monodromies(_Kf, E, cj)
    q0, p0, _ = base.orbit(E)
    res = dict(E=E, T=T, closure=clo, blocks={})
    for b in base.BLOCKS:
        M = Ms[b]
        mu, vec = np.linalg.eig(M)
        order = np.argsort(-np.abs(mu))
        mu, vec = mu[order], vec[:, order]
        g = base.ghost_multipliers(M)
        off8 = [m for m in mu if abs(np.log(abs(m))) > OFF_TOL]
        unst = classify(g, UNST_TOL)                                   # verdict + type: the ghost 4
        masked = any(min(abs(m - gg) for gg in g) > 1e-8 * abs(m) for m in off8)   # any off-circle one not in the ghost 4?
        gv = None
        if unst:
            k = int(np.argmin(np.abs(mu - unst[0][0])))
            gv = gauss_violation(b, q0, p0, cj, vec[:, k])
        res['blocks'][b] = dict(
            mu=mu, ghost=g, n_off=len(off8), unst=unst, masked=masked, gauss=gv,
            det=abs(np.linalg.det(M) - 1.0),
            pair=abs(np.prod(np.abs(g)) - 1.0),
            one=min(abs(mu - 1.0)),
            ym=base.summary(Ys[b])[1][0])
    return res


def fmt_mu(m):
    return f"{m.real:+.5f}{m.imag:+.5f}i"


def type_str(bk):
    """Short type string for a block: '-' if stable, else type of the largest, with |mu| and arg/pi,
    and '(x2)' if two independent unstable pairs/quadruplets are present."""
    if not bk['unst']:
        return '-'
    m, t = bk['unst'][0]
    n_unst = len(bk['unst'])                     # multipliers with |mu| > 1: 1 per real pair, 2 per quadruplet
    extra = ''
    if t == 'CPLX' and n_unst > 2:
        extra = ' (+more)'
    if t.startswith('REAL') and n_unst > 1:
        extra = f' (x{n_unst})'
    return f"{t} |mu|={abs(m):.4f} arg/pi={np.angle(m) / np.pi:+.4f}{extra}"


def print_row(r, prev=None):
    cols = []
    for b in base.BLOCKS:
        bk = r['blocks'][b]
        mx = abs(bk['ghost'][0])
        cols.append(f"{mx:8.4f} {type_str(bk):<32s} n8={bk['n_off']} "
                    + (f"Gv={bk['gauss']:.1e}" if bk['gauss'] is not None else "Gv=   -   ")
                    + (" MASKED" if bk['masked'] else ""))
    ym = max(r['blocks'][b]['ym'] for b in base.BLOCKS)
    det = max(r['blocks'][b]['det'] for b in base.BLOCKS)
    pair = max(r['blocks'][b]['pair'] for b in base.BLOCKS)
    one = max(r['blocks'][b]['one'] for b in base.BLOCKS)
    print(f"E={r['E']:.5f} T={r['T']:8.4f} clo={r['closure']:.1e} det-1={det:.1e} pair={pair:.1e} "
          f"|mu-1|min={one:.1e} YM={ym:.4f}", flush=True)
    for b, c in zip(base.BLOCKS, cols):
        print(f"      blk{b} {c}", flush=True)
    # all 8 wherever the off-circle count changed from the previous grid point
    if prev is not None:
        for b in base.BLOCKS:
            if prev['blocks'][b]['n_off'] != r['blocks'][b]['n_off']:
                print(f"      n8 changed in blk{b} ({prev['blocks'][b]['n_off']} -> {r['blocks'][b]['n_off']}); "
                      f"all 8 at E={prev['E']:.5f}: " + "  ".join(fmt_mu(m) for m in prev['blocks'][b]['mu']), flush=True)
                print(f"      {'':>44s} all 8 at E={r['E']:.5f}: " + "  ".join(fmt_mu(m) for m in r['blocks'][b]['mu']),
                      flush=True)


def checks(res_list, label):
    clo = max(r['closure'] for r in res_list)
    det = max(r['blocks'][b]['det'] for r in res_list for b in base.BLOCKS)
    pair = max(r['blocks'][b]['pair'] for r in res_list for b in base.BLOCKS)
    one = max(r['blocks'][b]['one'] for r in res_list for b in base.BLOCKS)
    masked = any(r['blocks'][b]['masked'] for r in res_list for b in base.BLOCKS)
    ok = True
    for name, val, tol in (('closure', clo, 1e-9), ('det-1', det, 1e-8), ('pairing', pair, 1e-6), ('|mu-1|min', one, 1e-3)):
        print(f"CHECK {label} max {name} = {val:.1e} (tol {tol:.0e}) {'PASS' if val < tol else 'FAIL'}")
        ok &= val < tol
    print(f"CHECK {label} no unstable multiplier masked by the ghost-4 rule: {'PASS' if not masked else 'FAIL'}")
    return ok and not masked


def hill_check(E):
    """c_j = 0: the block's multipliers must contain those of the two scalar Hill equations."""
    T, clo, Ms, _ = base.monodromies(_Kf, E, 0.0)
    hill = sorted([m for s in (+1, -1) for m in np.linalg.eigvals(base.hill_monodromy(E, s))], key=lambda z: -abs(z))
    mu = np.linalg.eigvals(Ms['12'])
    worst = max(min(abs(h - m) for m in mu) / max(abs(h), 1.0) for h in hill)
    print(f"CHECK c_j=0 Hill reduction at E={E}: worst relative mismatch {worst:.1e} (tol 1e-6) "
          f"{'PASS' if worst < 1e-6 else 'FAIL'}")
    return worst < 1e-6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--grid', choices=['full', 'coarse'], default='full')
    ap.add_argument('--ncores', type=int, default=1)
    ap.add_argument('--emin', type=float, default=0.05)
    ap.add_argument('--emax', type=float, default=1.0)
    ap.add_argument('--npts', type=int, default=400)
    a = ap.parse_args()
    t0 = time.time()
    print(f"dr_offdiag_floquet_typed.py  python {sys.version.split()[0]}  numpy {np.__version__}  "
          f"sympy {sp.__version__}  ({sys.executable})")
    print(f"run parameters: RT={base.RT}  UNST_TOL={UNST_TOL}  OFF_TOL={OFF_TOL}  IM_TOL={IM_TOL}  "
          f"grid={a.grid} [{a.emin}, {a.emax}] x {a.npts} log-spaced  ncores={a.ncores}  c_j=1 (the theory)")
    print("orbit: DR, X0=%r PX0=%r T_HALF=%r (E=1/2, scaled)" % (base.X0, base.PX0, base.T_HALF))
    _init()

    ok = True
    for E in (0.5, 0.3):
        ok &= hill_check(E)

    print("\n== Coarse table (the 16 energies of the 2026-10-18 session checks, for comparison) ==")
    coarse = []
    prev = None
    for E in COARSE:
        r = one_energy(E)
        print_row(r, prev)
        coarse.append(r)
        prev = r
    ok &= checks(coarse, 'coarse')
    if a.grid == 'coarse':
        print(f"\n{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}")
        print(f"DONE {time.time() - t0:.0f}")
        return

    print(f"\n== Fine scan, E in [{a.emin}, {a.emax}], {a.npts} log-spaced points, typed ==")
    grid = np.logspace(np.log10(a.emin), np.log10(a.emax), a.npts)
    results = []
    prev = None
    if a.ncores > 1:
        with Pool(a.ncores, initializer=_init) as pool:
            for r in pool.imap(one_energy, grid):
                print_row(r, prev)
                results.append(r)
                prev = r
    else:
        for E in grid:
            r = one_energy(E)
            print_row(r, prev)
            results.append(r)
            prev = r
    ok &= checks(results, 'fine')

    print("\n== Unstable intervals of the fine scan (any block), typed ==")
    flags = np.array([any(r['blocks'][b]['unst'] for b in base.BLOCKS) for r in results])
    k = 0
    n_int = 0
    while k < len(grid):
        if flags[k]:
            j = k
            while j + 1 < len(grid) and flags[j + 1]:
                j += 1
            n_int += 1
            print(f"  interval {n_int}: E in [{grid[k]:.4f}, {grid[j]:.4f}]  ({j - k + 1} points)")
            for b in base.BLOCKS:
                types = [r['blocks'][b]['unst'][0][1] if r['blocks'][b]['unst'] else '-' for r in results[k:j + 1]]
                if all(t == '-' for t in types):
                    continue
                mxs = [abs(r['blocks'][b]['unst'][0][0]) if r['blocks'][b]['unst'] else 1.0 for r in results[k:j + 1]]
                im = int(np.argmax(mxs))
                runs = []
                for t in types:
                    if not runs or runs[-1][0] != t:
                        runs.append([t, 1])
                    else:
                        runs[-1][1] += 1
                seq = " -> ".join(f"{t}x{n}" for t, n in runs)
                changes = len({t for t in types if t != '-'}) > 1
                print(f"      blk{b}: max |mu| = {mxs[im]:.4f} at E={grid[k + im]:.4f}; type sequence {seq}; "
                      f"type changes inside: {'YES' if changes else 'no'}")
            k = j + 1
        else:
            k += 1
    print(f"  unstable fraction of the grid: {flags.mean():.3f}")
    n_c = sum(1 for r in results for b in base.BLOCKS if r['blocks'][b]['unst'] and r['blocks'][b]['unst'][0][1] == 'CPLX')
    n_r = sum(1 for r in results for b in base.BLOCKS if r['blocks'][b]['unst'] and r['blocks'][b]['unst'][0][1] != 'CPLX')
    print(f"  (block, E) pairs unstable: {n_c} Krein (CPLX), {n_r} parametric (REAL)")
    print(f"\n{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}")
    print(f"DONE {time.time() - t0:.0f}")


if __name__ == "__main__":
    main()
