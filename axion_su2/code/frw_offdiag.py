"""frw_offdiag.py -- the off-diagonal Lee-Wick blocks on the isotropic orbit with M -> M a(eta).

Uses rung1_frw_system.pkl (rung1_frw.py): the 24-dim first-order system x' = K(q, q', c_j; a, a', a'') x
of the six off-diagonal components on the ghost-free background, which in conformal time is the flat
isotropic orbit x = y = z = s, s = cn(sqrt2 eta | 1/2) (comoving amplitude 1; any amplitude rescales to
this), period T0 = 4 K(1/2) / sqrt2.  Units M = 1 (the fixed Lee-Wick parameter).  The three 8-dim
blocks (u12,u21), (u13,u31), (u23,u32) decouple and are equal by isotropy (checked).

Part A, frozen a (a' = a'' = 0): the Floquet multipliers of a block over one period.  A frozen a
equals flat space with M -> a, which by the scaling q -> q/M, t -> M t equals the committed flat
system (ym/rung1_system.pkl, iso_floquet.py) at E = (3/2)/a^4 with the same multipliers; checked at
spot values against iso_floquet.monodromies.  The frozen scan is then run over a in [0.25, 5], i.e.
E in [0.0024, 384], which extends the stage-1 scan (E <= 10) down to a = 0.25.

Part B, radiation era a = eta/eta_s (a' = 1/eta_s, a'' = 0): the transfer matrix Phi of a block from
a0 to a1 and ln sigma_max(Phi), for eta_s spanning H/M = a'/a^2 = 1/(eta_s a^2) from ~1 to ~1e-3 at
a = 1.  With every frozen a stable (Part A, c_j = 1) the question is whether the time dependence
itself (the a', a'' terms) drives the off-diagonal ghost: if not, ln sigma_max stays bounded as the
crossing slows (eta_s grows); if a frozen band exists (c_j = 0, the diagnostic, unstable for
E >~ 1.5 i.e. a <~ 1) it grows in proportion to the conformal time spent in the band.  The slowest
crossing is the adiabatic reference; no WKB normalisation is attempted for the 4th-order block, so
ln sigma_max carries an O(1) adiabatic part from the change of scale of the derivative components.
Also checked: the linearised Gauss charge is conserved along the FRW evolution (G(eta1) Phi v = 0
for v in ker G(eta0)).

Float only (DOP853, rtol 1e-11); nothing here is certified.
Usage (from the repository root, after rung1_frw.py; a few minutes):
    python axion_su2/code/frw_offdiag.py > axion_su2/results/frw_offdiag.txt
Added 2026-10-14.
"""
import os
import pickle
import sys
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.special import ellipk

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import iso_floquet as iso                                   # noqa: E402  (committed flat system, for the check)

PKL = os.path.join(HERE, '..', 'rung1_frw_system.pkl')
BLOCKS = iso.BLOCKS
RT = dict(method='DOP853', rtol=1e-11, atol=1e-13)
T0 = 4.0 * ellipk(0.5) / np.sqrt(2.0)


def load_frw():
    d = pickle.load(open(PKL, 'rb'))
    args = list(d['X']) + list(d['P']) + [d['cj']] + list(d['asyms'])
    Kb, Gb = {}, {}
    for name in BLOCKS:
        idx = iso.block_idx(name)
        Kb[name] = sp.lambdify(args, d['K'].extract(idx, idx), 'numpy')
        Gb[name] = sp.lambdify(args, d['G'].extract([0, 1, 2], idx), 'numpy')
    return Kb, Gb


def bg_rhs(s, sd):
    return sd, -2.0 * s ** 3


def transfer(Kf, cj, afun, e0, e1, s0=(1.0, 0.0)):
    """Phi(e0 -> e1) of one 8-dim block along the isotropic orbit; afun(eta) -> (a, a', a'')."""
    def rhs(e, y):
        s, sd = y[0], y[1]
        a, a1, a2 = afun(e)
        K = np.array(Kf(s, s, s, sd, sd, sd, cj, a, a1, a2), dtype=float)
        M = y[2:].reshape(8, 8)
        return np.concatenate([bg_rhs(s, sd), (K @ M).ravel()])
    y0 = np.concatenate([s0, np.eye(8).ravel()])
    sol = solve_ivp(rhs, (e0, e1), y0, **RT)
    y = sol.y[:, -1]
    return y[2:].reshape(8, 8), y[:2]


def phase_state(ph):
    if ph == 0.0:
        return (1.0, 0.0)
    y = solve_ivp(lambda t, y: bg_rhs(*y), (0.0, ph * T0), [1.0, 0.0], **RT).y[:, -1]
    return (y[0], y[1])


def main():
    t0 = time.time()
    print(f"frw_offdiag.py  python {sys.version.split()[0]}  numpy {np.__version__}  ({sys.executable})")
    print(f"isotropic orbit, comoving amplitude 1, period T0 = {T0:.4f}; M = 1")
    Kb, Gb = load_frw()
    Kflat = iso.load_K()

    # ---------------------------------------------------------------- Part A: frozen a
    print("\nPart A: frozen a.  Check against the committed flat system at E = (3/2)/a^4 (largest ghost multiplier)")
    print(f"{'a':>6} {'E=1.5/a^4':>10} {'FRW frozen':>11} {'iso_floquet':>11} {'diff':>8} {'spread/3 blocks':>15}")
    for a in (0.7, 1.0, 1.5, 2.0, 3.0):
        E = 1.5 / a ** 4
        g = []
        for name in BLOCKS:
            Phi, _ = transfer(Kb[name], 1.0, lambda e: (a, 0.0, 0.0), 0.0, T0)
            g.append(abs(iso.ghost_multipliers(Phi)[0]))
        _, _, Ms = iso.monodromies(Kflat, E, 1.0)
        ref = abs(iso.ghost_multipliers(Ms['Ogh12'])[0])
        print(f"{a:6.2f} {E:10.4f} {g[0]:11.6f} {ref:11.6f} {g[0] - ref:8.1e} {max(g) - min(g):15.1e}", flush=True)

    for cj in (1.0, 0.0):
        grid = np.logspace(np.log10(0.25), np.log10(5.0), 120)
        flags, rates = [], []
        for a in grid:
            Phi, _ = transfer(Kb['12'], cj, lambda e: (a, 0.0, 0.0), 0.0, T0)
            mu = abs(iso.ghost_multipliers(Phi)[0])
            flags.append(mu > 1.0 + 1e-5)
            rates.append(np.log(mu) / T0)
        iv = iso.intervals(grid, np.array(flags))
        print(f"\n  frozen scan a in [0.25, 5] (120 pts, log), c_j = {cj:g}: off-diagonal ghost unstable on "
              f"{len(iv)} interval(s); fraction {np.mean(flags):.2f}")
        for lo, hi, n in iv:
            sel = (grid >= lo) & (grid <= hi)
            print(f"     a in [{lo:.4f}, {hi:.4f}]  (E in [{1.5 / hi ** 4:.4g}, {1.5 / lo ** 4:.4g}], {n} pts)  "
                  f"max rate {np.max(np.array(rates)[sel]):.4g} per unit eta", flush=True)

    # ---------------------------------------------------------------- Part B: radiation-era crossings
    print("\nPart B: radiation era a = eta/eta_s.  ln sigma_max of the block transfer matrix from a0 to a1,")
    print("        two background phases (0, T0/4); Gauss-charge conservation |G1 Phi v|, v in ker G0 (max).")
    for (a0, a1) in ((0.3, 3.0), (1.4, 2.3)):
        for cj in (1.0, 0.0):
            print(f"\n  a from {a0} to {a1}, c_j = {cj:g}")
            print(f"{'eta_s':>8} {'H/M at a=1':>11} {'eta span':>9} | {'ln s_max (ph 0)':>15} {'(ph 1/4)':>9} | {'Gauss':>8}")
            for eta_s in (1, 3, 10, 30, 100, 300, 1000):
                afun = lambda e, es=eta_s: (e / es, 1.0 / es, 0.0)
                e0, e1 = a0 * eta_s, a1 * eta_s
                vals = []
                for ph in (0.0, 0.25):
                    s0 = phase_state(ph)
                    Phi, sb = transfer(Kb['12'], cj, afun, e0, e1, s0)
                    vals.append(np.log(np.linalg.svd(Phi, compute_uv=False)[0]))
                    if ph == 0.0:
                        G0 = np.array(Gb['12'](*([s0[0]] * 3), *([s0[1]] * 3), cj, *afun(e0)), dtype=float)
                        G1 = np.array(Gb['12'](*([sb[0]] * 3), *([sb[1]] * 3), cj, *afun(e1)), dtype=float)
                        _, _, Vt = np.linalg.svd(G0)
                        ker = Vt[np.linalg.matrix_rank(G0):].T
                        gauss = np.max(np.abs(G1 @ Phi @ ker)) / max(1.0, np.linalg.norm(Phi, 2))
                print(f"{eta_s:8d} {1.0 / eta_s:11.2e} {e1 - e0:9.1f} | {vals[0]:15.3f} {vals[1]:9.3f} | {gauss:8.1e}", flush=True)
    print(f"\n({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
