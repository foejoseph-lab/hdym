"""channel_offdiag_monodromy.py -- B002: the off-diagonal ghost inside a channel, monodromy over one
transverse period at fixed depth, typed.  Open item (iv)'s proper tool (WORKORDER_bridge_2026-10-18_v2.md).

Background (the channel model, not an exact orbit): the deep coordinate is frozen at x = D, px = 0; the two
transverse coordinates oscillate at the channel's frequency D with action J each,
    y = A cos(D t),   z = A cos(D t + psi),   A = sqrt(2J/D),   py = ydot, pz = zdot,
over one transverse period T = 2 pi / D.  Frozen-depth error: with px = 0 the background sits at the deep
turning point, where really x'' = -x(y^2+z^2) ~ -2J, so over one transverse period the depth changes by
(2J/2)(2 pi/D)^2 / D = 4 pi^2 J / D^3 =: eps relative to D (the dropped z^2 y, y^2 z cross terms are smaller,
2J/D^3).  eps is printed per row; rows with eps > 0.1 are outside the model and marked OUTSIDE-MODEL.  Away
from the turning point the deep coordinate moves at speed ~ sqrt(2E) ~ 1 and the adiabaticity of the transverse
oscillation is 1/D^2; the model is a fixed-depth idealisation either way.  The 24-dim linearisation K(q, p, c_j) is the one hdym_rung1.py exports (ym/rung1_system.pkl);
the three 8-dim blocks (u12,u21), (u13,u31), (u23,u32) are integrated as in dr_offdiag_floquet.py (its
BLOCKS/block_idx/ym_mass/gradV are imported: one definition).  Blocks 12 and 13 contain the deep coordinate;
block 23 does not (F118: it is expected stable with frequencies ~ D).

On a background that is not an exact solution the gauge zero mode is lost at O(eps) and its defective
multiplier 1 (a Jordan pair) splits by ~sqrt(eps), so dr_offdiag_floquet.py's "drop the 4 nearest 1" rule has
no meaning here: the verdict and type come from ALL 8 multipliers with the off-circle tolerance OFF_TOL, and
the Gauss-law column tells a physical growing mode from growth in the broken constraint sector.  The
(mu, 1/mu) pairing of all 8 is checked (the linearised flow is symplectic on any background: 1e-14 measured).

Per (c_j, psi, D, J) and block: the largest multiplier off the circle, its TYPE
    REAL+/REAL-  real pair, parametric;   CPLX  complex quadruplet, Krein collision
(definitions as in dr_offdiag_floquet_typed.py), |mu|, arg/pi, the count n8 of all 8 multipliers off the
circle (n8 = 2 one real pair, 4 one quadruplet or two real pairs, 6 both), the Gauss-law violation |G v|/|v| of the leading unstable eigenvector (F118: in-channel frozen
coefficients gave Gauss-violating growth, |Gv|/|v| = 1-100; a monodromy over the real transverse motion should
not), min |mu - 1| (the gauge zero mode; on a non-exact background it need not be exactly 1, so it is a column,
not a check), and the pure-YM off-diagonal block's max |mu| on the same background.

Checks (CHECK ... PASS/FAIL): det M = 1 and all-8 (mu, 1/mu) pairing for every row; the c_j = 0 reduction of
blocks 12 and 13 to the scalar Hill equation w'' = -(1 + q_c^2 -+ q_a q_b) w (q_c the third coordinate) at one
small-eps (D, J), with tolerance eps: the exported K has the background equations of motion substituted
(q'' -> -grad V), so on the frozen background it differs from the direct linearisation by O(eps) (measured:
mismatch 2e-2 at eps 0.25, 1e-5 at eps 3e-3, 1e-6 at eps 1e-3).  Block 23 is reported as INFO only: at c_j = 0
its frequencies are ~ D, the 1:1 resonance with the drive, so its multipliers sit in the sqrt(eps)-split
cluster near 1 and the comparison is degenerate there.
Ends with a type map over the (D, J) grid per (psi, block) and DONE <s>.

Cap (work order): one scan in (D, J); if the real/complex split is not clean, record what was seen and stop.
Found (results/channel_offdiag_monodromy.txt, 19 Oct; islands/README.md): Floquet-stable in every block within the
model (eps <= 0.1, 188 of 576 cases); the one in-model "weak Krein cell" (D = 3.89, J = 0.0315, |mu| = 1.002, rate
0.001) is the gauge Jordan pair split by sqrt(eps) and labelled Krein by the all-8 rule, not a ghost instability
(rabbit_holes/ym/notes/RH11_weakcell.md); every strong instability lies outside the model and most of it violates
Gauss.  The DR orbit itself has eps ~ 0.2, outside this model.

Usage (Spyder, from the repo root, kernel's own interpreter):
    import sys; %cd C:\\GitHub\\hdym
    !"{sys.executable}" -B -u islands\\code\\channel_offdiag_monodromy.py --ncores 20 > islands\\results\\channel_offdiag_monodromy.txt 2>&1

Run parameters and why: D in [0.5, 6] (24 log points): the DR orbit's deep turning point is X0 (2E)^{1/4} =
1.8 ... 3.7 over E = 0.05 ... 1, so the paper's channels live in [1, 4]; the ends show the trend.  J in
[0.005, 0.5] (6 log points): J = E_perp / D; E_perp of order the DR orbit's transverse energy and below.
psi in {0, pi/2}: the two symmetric relative phases (y = z, and y^2 + z^2 = A^2 constant); the real channel has
any phase.  c_j = 1 (the theory) and 0 (the reduction check and the ghost-free comparison).  Tolerances as
dr_offdiag_floquet_typed.py (RT, UNST_TOL, OFF_TOL, IM_TOL).  Units g = M = 1.
Added 2026-10-19.
"""
import argparse
import os
import sys
import time
from itertools import product
from multiprocessing import Pool

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dr_offdiag_floquet as base                               # noqa: E402
from dr_offdiag_floquet_typed import (UNST_TOL, OFF_TOL, IM_TOL, classify, load_G)   # noqa: E402

THIRD = {'12': 2, '13': 1, '23': 0}            # index of the coordinate not in the block
PAIR = {'12': (0, 1), '13': (0, 2), '23': (1, 2)}
_Kf = None
_Gf = None


def _init():
    global _Kf, _Gf
    _Kf = base.load_K()
    _Gf = load_G()


def background(D, J, psi):
    A = np.sqrt(2.0 * J / D)

    def qp(t):
        q = np.array([D, A * np.cos(D * t), A * np.cos(D * t + psi)])
        p = np.array([0.0, -A * D * np.sin(D * t), -A * D * np.sin(D * t + psi)])
        return q, p
    return qp, 2.0 * np.pi / D


def monodromies(D, J, psi, cj):
    qp, T = background(D, J, psi)

    def rhs(t, s):
        q, p = qp(t)
        K = np.array(_Kf(*q, *p, cj), dtype=float)
        out = []
        off = 0
        for b in base.BLOCKS:
            idx = base.block_idx(b)
            M = s[off:off + 64].reshape(8, 8)
            out.append((K[np.ix_(idx, idx)] @ M).ravel())
            off += 64
        for b in base.BLOCKS:
            M = s[off:off + 16].reshape(4, 4)
            Am = np.zeros((4, 4))
            Am[0, 2] = Am[1, 3] = 1.0
            Am[2:, :2] = -base.ym_mass(q, b)
            out.append((Am @ M).ravel())
            off += 16
        return np.concatenate(out)

    nb = len(base.BLOCKS)
    s0 = np.concatenate([np.eye(8).ravel()] * nb + [np.eye(4).ravel()] * nb)
    sT = solve_ivp(rhs, (0.0, T), s0, **base.RT).y[:, -1]
    Ms, Ys = {}, {}
    off = 0
    for b in base.BLOCKS:
        Ms[b] = sT[off:off + 64].reshape(8, 8)
        off += 64
    for b in base.BLOCKS:
        Ys[b] = sT[off:off + 16].reshape(4, 4)
        off += 16
    return T, Ms, Ys


def hill_monodromy(D, J, psi, b):
    """c_j = 0 reduction: w'' = -(1 + q_c^2 - sign q_a q_b) w, sign = +-1, 2 x 2 monodromies, on the same background."""
    qp, T = background(D, J, psi)
    ia, ib = PAIR[b]
    ic = THIRD[b]
    out = []
    for sign in (+1, -1):
        def rhs(t, s):
            q, _ = qp(t)
            om2 = 1.0 + q[ic] ** 2 - sign * q[ia] * q[ib]
            return [s[1], -om2 * s[0], s[3], -om2 * s[2]]
        sT = solve_ivp(rhs, (0.0, T), [1.0, 0.0, 0.0, 1.0], **base.RT).y[:, -1]
        out.extend(np.linalg.eigvals(np.array([[sT[0], sT[2]], [sT[1], sT[3]]])))
    return out


def one_case(args):
    cj, psi, D, J = args
    T, Ms, Ys = monodromies(D, J, psi, cj)
    q0, p0 = background(D, J, psi)[0](0.0)
    res = dict(cj=cj, psi=psi, D=D, J=J, T=T, eps=4.0 * np.pi ** 2 * J / D ** 3, blocks={})
    for b in base.BLOCKS:
        M = Ms[b]
        mu, vec = np.linalg.eig(M)
        order = np.argsort(-np.abs(mu))
        mu, vec = mu[order], vec[:, order]
        off8 = [m for m in mu if abs(np.log(abs(m))) > OFF_TOL]
        unst = classify(mu, np.expm1(OFF_TOL))                 # all 8, |mu| > exp(OFF_TOL)
        pair = max(min(abs(m * m2 - 1.0) for m2 in mu) for m in mu)   # every mu has its 1/mu among the 8
        gv = None
        if unst:
            k = int(np.argmin(np.abs(mu - unst[0][0])))
            r, gf = _Gf[b]
            row = np.array(gf(*q0, *p0, cj), dtype=complex).ravel()[base.block_idx(b)]
            gv = abs(row @ vec[:, k]) / np.linalg.norm(vec[:, k])
        res['blocks'][b] = dict(mu=mu, n_off=len(off8), unst=unst, gauss=gv,
                                det=abs(np.linalg.det(M) - 1.0), pair=pair,
                                one=min(abs(mu - 1.0)), ym=base.summary(Ys[b])[1][0])
    return res


def type_code(bk):
    """One character for the map: '.' stable, R/C by type, lower case if the leading eigenvector violates Gauss."""
    if not bk['unst']:
        return '.'
    c = 'C' if bk['unst'][0][1] == 'CPLX' else 'R'
    return c.lower() if (bk['gauss'] is not None and bk['gauss'] > 0.1) else c


def print_row(r):
    head = (f"cj={r['cj']:.0f} psi/pi={r['psi'] / np.pi:.2f} D={r['D']:7.4f} J={r['J']:7.4f} T={r['T']:7.4f} "
            f"eps={r['eps']:.2e}{' OUTSIDE-MODEL' if r['eps'] > 0.1 else ''}")
    print(head, flush=True)
    for b in base.BLOCKS:
        bk = r['blocks'][b]
        if bk['unst']:
            m, t = bk['unst'][0]
            ts = f"{t} |mu|={abs(m):.4f} arg/pi={np.angle(m) / np.pi:+.4f} rate={np.log(abs(m)) / r['T']:.4f}"
        else:
            ts = '-'
        print(f"      blk{b} {abs(bk['mu'][0]):8.4f} {ts:<48s} n8={bk['n_off']} "
              + (f"Gv={bk['gauss']:.1e} " if bk['gauss'] is not None else "Gv=   -    ")
              + f"|mu-1|min={bk['one']:.1e} det-1={bk['det']:.1e} pair={bk['pair']:.1e} YM={bk['ym']:.4f}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ncores', type=int, default=1)
    ap.add_argument('--nD', type=int, default=24)
    ap.add_argument('--nJ', type=int, default=6)
    ap.add_argument('--Dmin', type=float, default=0.5)
    ap.add_argument('--Dmax', type=float, default=6.0)
    ap.add_argument('--Jmin', type=float, default=0.005)
    ap.add_argument('--Jmax', type=float, default=0.5)
    a = ap.parse_args()
    t0 = time.time()
    print(f"channel_offdiag_monodromy.py  python {sys.version.split()[0]}  numpy {np.__version__}  "
          f"sympy {sp.__version__}  ({sys.executable})")
    Ds = np.logspace(np.log10(a.Dmin), np.log10(a.Dmax), a.nD)
    Js = np.logspace(np.log10(a.Jmin), np.log10(a.Jmax), a.nJ)
    psis = [0.0, np.pi / 2]
    print(f"run parameters: RT={base.RT} UNST_TOL={UNST_TOL} OFF_TOL={OFF_TOL} IM_TOL={IM_TOL}  "
          f"D in [{a.Dmin}, {a.Dmax}] x {a.nD} log  J in [{a.Jmin}, {a.Jmax}] x {a.nJ} log  psi in {{0, pi/2}}  "
          f"c_j in {{1, 0}}  ncores={a.ncores}")
    _init()

    ok = True
    D0, J0, psi0 = 4.0, 0.005, np.pi / 2
    r0 = one_case((0.0, psi0, D0, J0))
    eps0 = r0['eps']
    for b in base.BLOCKS:
        hill = hill_monodromy(D0, J0, psi0, b)
        mu = r0['blocks'][b]['mu']
        worst = max(min(abs(h - m) for m in mu) / max(abs(h), 1.0) for h in hill)
        if b == '23':
            print(f"INFO  c_j=0 Hill reduction, block 23 at D={D0} J={J0} psi=pi/2: mismatch {worst:.1e} "
                  f"(degenerate: 1:1 resonance with the drive, cluster near 1 split by sqrt(eps) = {np.sqrt(eps0):.1e})")
            continue
        print(f"CHECK c_j=0 Hill reduction, block {b} at D={D0} J={J0} psi=pi/2 (eps={eps0:.1e}): worst relative "
              f"mismatch {worst:.1e} (tol eps) {'PASS' if worst < eps0 else 'FAIL'}")
        ok &= worst < eps0

    cases = [(cj, psi, D, J) for cj in (1.0, 0.0) for psi in psis for J in Js for D in Ds]
    print(f"\n== Scan: {len(cases)} cases (order: c_j, psi, J, D) ==")
    results = []
    if a.ncores > 1:
        with Pool(a.ncores, initializer=_init) as pool:
            for r in pool.imap(one_case, cases):
                print_row(r)
                results.append(r)
    else:
        for c in cases:
            r = one_case(c)
            print_row(r)
            results.append(r)

    inm = [r for r in results if r['eps'] <= 0.1]
    det = max(r['blocks'][b]['det'] for r in inm for b in base.BLOCKS)
    pair = max(r['blocks'][b]['pair'] for r in inm for b in base.BLOCKS)
    print(f"\nCHECK in-model ({len(inm)} of {len(results)} cases) max det-1 = {det:.1e} (tol 1e-8) {'PASS' if det < 1e-8 else 'FAIL'}")
    print(f"CHECK in-model max all-8 (mu, 1/mu) pairing = {pair:.1e} (tol 1e-6) {'PASS' if pair < 1e-6 else 'FAIL'}")
    ok &= det < 1e-8 and pair < 1e-6
    print(f"INFO  whole scan max det-1 = {max(r['blocks'][b]['det'] for r in results for b in base.BLOCKS):.1e}, "
          f"max pairing = {max(r['blocks'][b]['pair'] for r in results for b in base.BLOCKS):.1e} "
          f"(outside-model cells reach |mu| ~ 1e5, where the float determinant is meaningless)")
    print(f"INFO  scan max |mu-1|min (gauge zero mode on the frozen-depth background) = "
          f"{max(r['blocks'][b]['one'] for r in results for b in base.BLOCKS):.1e}")

    print("\n== Type map over (D, J): one character per D (left = Dmin), one line per J; "
          "'.' stable, R real pair (parametric), C complex quadruplet (Krein); lower case = leading eigenvector "
          "violates Gauss (|Gv|/|v| > 0.1); '!' after the J if any cell has eps = 4 pi^2 J/D^3 > 0.1 ==")
    for cj in (1.0, 0.0):
        for psi in psis:
            for b in base.BLOCKS:
                print(f"  c_j={cj:.0f} psi/pi={psi / np.pi:.2f} block {b}:")
                for J in Js:
                    row = [r for r in results if r['cj'] == cj and r['psi'] == psi and r['J'] == J]
                    row.sort(key=lambda r: r['D'])
                    s = "".join(type_code(r['blocks'][b]) for r in row)
                    bang = '!' if any(r['eps'] > 0.1 for r in row) else ' '
                    print(f"     J={J:7.4f}{bang} {s}")
    print(f"     D grid: " + " ".join(f"{D:.2f}" for D in Ds))

    print("\n== Split summary, c_j = 1, within the model (eps <= 0.1) ==")
    for psi in psis:
        for b in base.BLOCKS:
            cells = [r for r in results if r['cj'] == 1.0 and r['psi'] == psi and r['eps'] <= 0.1]
            R = [(r['D'], r['J']) for r in cells if r['blocks'][b]['unst'] and r['blocks'][b]['unst'][0][1] != 'CPLX']
            C = [(r['D'], r['J']) for r in cells if r['blocks'][b]['unst'] and r['blocks'][b]['unst'][0][1] == 'CPLX']
            G = [(r['D'], r['J']) for r in cells if r['blocks'][b]['gauss'] is not None and r['blocks'][b]['gauss'] > 0.1]

            def rng(L):
                if not L:
                    return "none"
                return f"{len(L)} cells, D in [{min(d for d, _ in L):.3f}, {max(d for d, _ in L):.3f}], J in [{min(j for _, j in L):.4f}, {max(j for _, j in L):.4f}]"
            print(f"  psi/pi={psi / np.pi:.2f} block {b}: {len(cells)} cells; REAL {rng(R)}; CPLX {rng(C)}; Gauss-violating {rng(G)}")
    print(f"\n{'ALL CHECKS PASS' if ok else 'SOME CHECK FAILED'}")
    print(f"DONE {time.time() - t0:.0f}")


if __name__ == "__main__":
    main()
