"""rh27_averaging.py -- rabbit holes 2 + 7: first-order Floquet (averaging) theory for block 12 on the DR orbit,
with the pen-and-paper parts done in sympy.  The question: does the averaged canonical system, with the first
Fourier harmonics of the background treated as a perturbation, PREDICT the ladder of band centres
T_c = 4.057 (n + 0.210) (notes/RH2_ladder.md) and the parametric/Krein alternation (notes/RH7_selection_rule.md)?

STATUS (notes/RH27_averaging.md, section C; README hole 2+7).  Parts A and B are exact and stand.  Part C's reading
"omega_0 T = n pi of the averaged ghost mode predicts the bands, odd n Krein" was made and WITHDRAWN the same evening
(21 Oct): the half-integer rungs (n + 1/2) pi also land within 8 % of a band 7 times in 9 because bands plus an 8 % halo
cover 79 % of the log-E range; the parity is a theorem given existence (RH7 lemma), not a prediction; rung 11 fails it
(E = 1.049 inside the parametric band [1.014, 1.271]).  What survives is slope plus count.  The averaged 8-dim block
is the wrong object (its hyperbolic sector is constraint contamination; the mode that tracked a pair was the dirty
one); the right zeroth order is the gauge-reduced 4-dim block (rh_reduced_block.py), and even there first-order
averaging is not Floquet-accurate.  The output of this script is kept as the record of that test.

Three parts.

A. CANONICAL FORM (sympy).  L2(u, u', u''; q, p) from ym/rung1_system.pkl, c_j = 1, has u''^2 coefficient -1/2.
   Ostrogradsky: Q1 = u, Q2 = u', P2 = dL2/du'', P1 = dL2/du' - d/dt P2, H = P1 Q2 + P2 u'' - L2 with u'' eliminated
   through P2.  H is a quadratic form z^T Hm(q, p) z / 2 in z = (Q1, Q2, P1, P2) with the canonical J, and the state
   map z = Jx(q, p) x from x = (u, u', u'', u''') is the one krein_b005.build_omega uses (Omega = Jx^T J_can Jx).
   CHECK (numerical, along the orbit): J Hm = (d/dt Jx + Jx K) Jx^{-1} -- the canonical equations reproduce K.
B. NORMAL FORM (sympy).  Two modes with signed frequencies nu_a = sigma_a omega_a, nu_b = sigma_b omega_b
   (sigma = Krein sign = sign of the mode's energy) and a resonant coupling 2 Re(g e^{i n Omega t} zeta_a zeta_b),
   zeta = (Q + iP)/sqrt2.  In the frame rotating at n Omega / 2 the Hamiltonian is constant; sympy gives the Floquet
   exponent lambda(detuning) exactly.  Result (printed): with Delta = nu_a + nu_b - n Omega,
   lambda^2 = |g|^2 - Delta^2/4 for the SUM resonance of signed frequencies (unstable inside |Delta| < 2|g|), while
   the difference resonance nu_a - nu_b = n Omega (coupling zeta_a conj zeta_b) is stable.  So: parametric band of
   mode a at 2 nu_a = n Omega, combination band of a, b at nu_a + nu_b = n Omega -- with SIGNED frequencies, which is
   where the Krein signature enters: two positive-energy modes combine at omega_a + omega_b, a ghost (sigma = -1)
   with a normal mode at omega_normal - omega_ghost.
C. SCAN (numpy).  Along the DR orbit at each E: sample Hm(t), Fourier-decompose, average -> H0; the eigenmodes of
   J H0 (zero modes dropped) give nu_a(E) and sigma_a(E); the matrix elements g_ab^(n) = <a| H_n |b> in the mode
   basis give the couplings.  Predicted bands: the energies where nu_a + nu_b = n Omega(E) (a <= b), their type
   (a = b: parametric; a != b: combination, Krein if sigma_a sigma_b = -1 ... as the normal form says), and their
   first-order half-width 2|g| in detuning.  Compared with B001's measured block-12 intervals and with the ladder's
   0.774 = pi/4.057 (plotted as a reference line; NOT a frequency -- the ladder is branch-independent, omega -> omega
   + k Omega shifts omega T by 2 pi k; the physical line is at 0.85-1.2 with a quasi-static companion, rh5_branch.py).

Cost: ~1.5 s per E (orbit + 1024 samples of the symbolic Hm); default nE = 120: ~3 min, no sized run.
Usage (from the repo root): python rabbit_holes/ym/code/rh27_averaging.py [--nE 120] > rabbit_holes/ym/results/rh27_averaging.txt
Writes results/rh27_averaging.npz and results/rh27_averaging.png.  Float where it is numeric.  Added 2026-10-21.
"""
import argparse
import os
import pickle
import re
import sys
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402

RES = os.path.join(HERE, '..', 'results')
B001 = os.path.join(ROOT, 'islands', 'results', 'dr_offdiag_floquet_typed.txt')
BLOCK_J = [0, 2]                              # ulist indices of (u12, u21)


# ---------------------------------------------------------------- A. canonical form from L2 (sympy)
def canonical_form():
    d = pickle.load(open(base.PKL, 'rb'))
    L2, X, P, cj = d['L2'], d['X'], d['P'], d['cj']
    L2 = L2.subs(cj, 1)
    U = [[sp.Symbol(f'u{k}_{j}') for j in range(6)] for k in range(4)]
    # restrict to block 12: set the other components (and their derivatives) to zero
    others = {U[k][j]: 0 for k in range(4) for j in range(6) if j not in BLOCK_J}
    L = sp.expand(L2.subs(others))
    Vs = (X[0] ** 2 * X[1] ** 2 + X[1] ** 2 * X[2] ** 2 + X[2] ** 2 * X[0] ** 2) / 2
    bg_rate = {X[i]: P[i] for i in range(3)}
    bg_rate.update({P[i]: -Vs.diff(X[i]) for i in range(3)})

    def ddt(expr):
        out = 0
        for k in range(3):
            for j in BLOCK_J:
                out += expr.diff(U[k][j]) * U[k + 1][j]
        for s_, r in bg_rate.items():
            out += expr.diff(s_) * r
        return sp.expand(out)
    x = [U[k][j] for k in range(4) for j in BLOCK_J]                         # (u12, u21, u12', u21', u12'', u21'', u12''', u21''')
    Q1 = [U[0][j] for j in BLOCK_J]
    Q2 = [U[1][j] for j in BLOCK_J]
    P2 = [sp.expand(L.diff(U[2][j])) for j in BLOCK_J]
    P1 = [sp.expand(L.diff(U[1][j]) - ddt(P2[jj])) for jj, j in enumerate(BLOCK_J)]
    canon = Q1 + Q2 + P1 + P2
    Jx = sp.Matrix(8, 8, lambda r, c: canon[r].coeff(x[c]))
    assert all(sp.expand(canon[r] - sum(Jx[r, c] * x[c] for c in range(8))) == 0 for r in range(8)), "canonical coordinates not linear in x"
    # H = P1 Q2 + P2 u'' - L with u'' eliminated: P2 = dL/du'' is linear in u'' with coefficient -1 (u''^2 term -1/2)
    z = sp.symbols('Q1a Q1b Q2a Q2b P1a P1b P2a P2b')
    udd = sp.symbols('udd_a udd_b')
    subs_u = {U[0][BLOCK_J[0]]: z[0], U[0][BLOCK_J[1]]: z[1], U[1][BLOCK_J[0]]: z[2], U[1][BLOCK_J[1]]: z[3],
              U[2][BLOCK_J[0]]: udd[0], U[2][BLOCK_J[1]]: udd[1]}
    P2z = [sp.expand(P2[jj].subs(subs_u)) for jj in range(2)]
    sol = sp.solve([sp.Eq(P2z[0], z[6]), sp.Eq(P2z[1], z[7])], udd, dict=True)[0]
    Lz = L.subs(subs_u).subs(sol)
    H = sp.expand(z[4] * z[2] + z[5] * z[3] + z[6] * sol[udd[0]] + z[7] * sol[udd[1]] - Lz)
    # the u''' entries must have dropped out (they do not appear in L2); H is quadratic in z
    Hm = sp.Matrix(8, 8, lambda r, c: sp.diff(H, z[r], z[c]))
    assert sp.expand(H - sum(Hm[r, c] * z[r] * z[c] for r in range(8) for c in range(8)) / 2) == 0, "H not quadratic"
    args = list(X) + list(P)
    Hf = sp.lambdify(args, Hm, 'numpy')
    Jf = sp.lambdify(args, Jx, 'numpy')
    Jdotf = sp.lambdify(args, Jx.applyfunc(ddt), 'numpy')
    return Hf, Jf, Jdotf, Hm, Jx


JCAN = np.block([[np.zeros((4, 4)), np.eye(4)], [-np.eye(4), np.zeros((4, 4))]])   # (Q, P) ordering


# ---------------------------------------------------------------- B. the two-mode normal form (sympy)
def normal_form():
    nu_a, nu_b, n, Om, dlt = sp.symbols('nu_a nu_b n Omega Delta', real=True)
    gr, gi = sp.symbols('g_r g_i', real=True)
    Qa, Pa, Qb, Pb = sp.symbols('Q_a P_a Q_b P_b', real=True)
    t = sp.symbols('t', real=True)
    za = (Qa + sp.I * Pa) / sp.sqrt(2)
    zb = (Qb + sp.I * Pb) / sp.sqrt(2)
    g = gr + sp.I * gi
    # H = nu_a |z_a|^2 + nu_b |z_b|^2 + 2 Re(g e^{i n Omega t} z_a z_b): the sum-resonant coupling
    Hres = nu_a * (Qa ** 2 + Pa ** 2) / 2 + nu_b * (Qb ** 2 + Pb ** 2) / 2 + 2 * sp.re(sp.expand(g * sp.exp(sp.I * n * Om * t) * za * zb))
    # rotating frame: each mode rotated by n Omega t / 2 removes the time dependence and shifts nu -> nu - n Omega / 2
    Hrot = Hres.subs(t, 0).subs({nu_a: nu_a - n * Om / 2, nu_b: nu_b - n * Om / 2})
    zz = [Qa, Qb, Pa, Pb]
    Hm = sp.Matrix(4, 4, lambda r, c: sp.diff(Hrot, zz[r], zz[c]))
    J = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
    A = J * Hm
    lam = sp.symbols('lambda')
    charpoly = sp.factor((A - lam * sp.eye(4)).det())
    # at exact resonance nu_a + nu_b = n Omega: the two shifted frequencies are opposite, +-Delta/2 in general
    dd = sp.symbols('d', real=True)                                   # d = nu_a - n Omega / 2, the frame-shifted frequency of a
    cp = sp.factor(charpoly.subs(nu_b, n * Om - nu_a + dlt).subs(nu_a, dd + n * Om / 2))
    sols = sp.solve(sp.Eq(cp, 0), lam)
    lam2 = sorted({sp.simplify(s_) for s_ in sols}, key=str)
    # the difference-resonant coupling 2 Re(g e^{i n Omega t} z_a conj z_b) for comparison
    Hdiff = nu_a * (Qa ** 2 + Pa ** 2) / 2 + nu_b * (Qb ** 2 + Pb ** 2) / 2 + 2 * sp.re(sp.expand(g * sp.exp(sp.I * n * Om * t) * za * sp.conjugate(zb)))
    Hrot2 = Hdiff.subs(t, 0).subs({nu_a: nu_a - n * Om / 2, nu_b: nu_b + n * Om / 2})
    Hm2 = sp.Matrix(4, 4, lambda r, c: sp.diff(Hrot2, zz[r], zz[c]))
    cp2 = sp.factor(((J * Hm2) - lam * sp.eye(4)).det().subs(nu_b, nu_a - n * Om + dlt).subs(nu_a, dd + n * Om / 2))
    sols2 = sp.solve(sp.Eq(cp2, 0), lam)
    lam2d = sorted({sp.simplify(s_) for s_ in sols2}, key=str)
    # the symmetric case d = 0 (both frame-shifted frequencies equal and opposite): the textbook forms
    lam2_0 = sorted({sp.simplify(s_.subs(dd, 0)) for s_ in sols}, key=str)
    lam2d_0 = sorted({sp.simplify(s_.subs(dd, 0)) for s_ in sols2}, key=str)
    return lam2, lam2d, charpoly, lam2_0, lam2d_0


# ---------------------------------------------------------------- C. the scan
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
            out.append((cur[0], cur[1], 'CPLX' if 'CPLX' in m.group(3) else 'REAL', float(m.group(1))))
    return out


def modes_at(E, Hf, Jf, Jdotf, Kf, nsamp=1024, nharm=24, check=False):
    q0, p0, T = base.orbit(E)
    sol = solve_ivp(lambda t, s: np.concatenate([s[3:], -base.gradV(s[:3])]), (0.0, T), np.concatenate([q0, p0]),
                    t_eval=np.linspace(0, T, nsamp, endpoint=False), **base.RT)
    Hs = np.array([np.array(Hf(*sol.y[:3, i], *sol.y[3:, i]), dtype=float) for i in range(nsamp)])
    resid = 0.0
    if check:
        idx = base.block_idx('12')
        for i in range(0, nsamp, nsamp // 8):
            q, p = sol.y[:3, i], sol.y[3:, i]
            Jx = np.array(Jf(*q, *p), dtype=float)
            Jd = np.array(Jdotf(*q, *p), dtype=float)
            K = np.array(Kf(*q, *p, 1.0), dtype=float)[np.ix_(idx, idx)]
            lhs = JCAN @ Hs[i]
            rhs = (Jd + Jx @ K) @ np.linalg.inv(Jx)
            resid = max(resid, np.abs(lhs - rhs).max() / max(np.abs(lhs).max(), 1e-300))
    F = np.fft.fft(Hs, axis=0) / nsamp                         # F[n] = (1/T) int H e^{-i n Omega t} dt
    H0 = F[0].real
    Hn = [F[n] for n in range(1, nharm + 1)]
    return T, H0, Hn, Hs, resid


def eigenmodes(H0, wmin=0.05):
    A = JCAN @ H0
    lam, V = np.linalg.eig(A)
    modes = []
    for i in range(len(lam)):
        if lam[i].imag > wmin and abs(lam[i].real) < 1e-6:
            e = V[:, i]
            en = (np.conj(e) @ H0 @ e).real
            sig = np.sign(en)
            # normalise so that the energy form is sigma * omega (|zeta|^2 = 1 convention: e* (i J^-1) e = ... ); here only
            # the sign and the frequency are used, and the couplings are normalised by the form below
            modes.append(dict(omega=lam[i].imag, sigma=sig, e=e / np.sqrt(abs(en) / lam[i].imag)))
    unstable0 = [l for l in lam if abs(l.real) > 1e-6]
    return sorted(modes, key=lambda m: m['omega']), unstable0


def coupling(ma, mb, Hn):
    """|<a| H_n |b>| with the modes normalised to unit energy/omega: the first-order half-width in detuning is 2|g|."""
    return abs(np.conj(ma['e']) @ Hn @ np.conj(mb['e'])) if ma is not mb else abs(np.conj(ma['e']) @ Hn @ np.conj(ma['e'])) / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nE', type=int, default=120)
    ap.add_argument('--emin', type=float, default=0.04)
    ap.add_argument('--emax', type=float, default=1.0)
    ap.add_argument('--nharm', type=int, default=24)
    a = ap.parse_args()
    t0 = time.time()
    print(f"rh27_averaging.py  python {sys.version.split()[0]}  numpy {np.__version__}  sympy {sp.__version__}")
    print("\n== A. canonical form of block 12 from L2 (sympy) ==")
    Hf, Jf, Jdotf, Hm_sym, Jx_sym = canonical_form()
    print(f"  H quadratic in (Q1, Q2, P1, P2) of (u12, u21): {Hm_sym.shape}, {sum(1 for r in range(8) for c in range(8) if Hm_sym[r, c] != 0)} nonzero entries; Jx: {sum(1 for r in range(8) for c in range(8) if Jx_sym[r, c] != 0)} nonzero")
    Kf = base.load_K()
    for E in (0.1, 0.3, 0.6):
        T, H0, Hn, Hs, resid = modes_at(E, Hf, Jf, Jdotf, Kf, check=True)
        print(f"  CHECK E={E}: J Hm = (dJx/dt + Jx K) Jx^-1 along the orbit to {resid:.1e} (tol 1e-8) {'PASS' if resid < 1e-8 else 'FAIL'}; Hm symmetric to {np.abs(Hs - Hs.transpose(0, 2, 1)).max():.1e}")

    print("\n== B. the two-mode resonance normal form (sympy) ==")
    lam2, lam2d, cp, lam2_0, lam2d_0 = normal_form()
    print("  (d = nu_a - n Omega/2, the frame-shifted frequency of mode a; the exponents are given in the rotating frame)")
    print("  sum resonance nu_a + nu_b = n Omega + Delta, coupling 2 Re(g e^{i n Omega t} zeta_a zeta_b): lambda in")
    for s_ in lam2:
        print("     ", s_)
    print("  difference resonance nu_a - nu_b = n Omega + Delta, coupling 2 Re(g e^{i n Omega t} zeta_a conj(zeta_b)): lambda in")
    for s_ in lam2d:
        print("     ", s_)
    print("  at d = 0 (symmetric detuning): sum resonance lambda in", [str(s_) for s_ in lam2_0], "; difference resonance lambda in", [str(s_) for s_ in lam2d_0])
    print("  reading: with SIGNED frequencies nu = sigma omega (sigma the Krein/energy sign), instability needs |Delta| < 2|g| at a")
    print("  sum resonance; difference resonances are stable.  A ghost mode (sigma = -1) therefore combines with a normal one at")
    print("  omega_normal - omega_ghost = n Omega, and a single mode is parametrically unstable at 2 nu = n Omega.")

    print("\n== C. the averaged modes and the predicted resonances along E ==")
    print("columns: E | T | Omega | modes of the averaged system: nu = sigma*omega (sigma, omega) ... | modes dropped (|omega| < 0.05, or unstable in H0)")
    grid = np.logspace(np.log10(a.emin), np.log10(a.emax), a.nE)
    rows = []
    for E in grid:
        T, H0, Hn, Hs, _ = modes_at(E, Hf, Jf, Jdotf, Kf, nharm=a.nharm)
        Om = 2 * np.pi / T
        modes, unst0 = eigenmodes(H0)
        rows.append(dict(E=E, T=T, Om=Om, modes=modes, Hn=Hn, unst0=unst0))
        desc = '  '.join(f"{m['sigma'] * m['omega']:+.4f} ({'+' if m['sigma'] > 0 else '-'}, {m['omega']:.4f})" for m in modes)
        print(f"E={E:.5f} T={T:8.4f} Om={Om:.4f} | {desc} | dropped {8 - 2 * len(modes)} ({len(unst0)} real/complex in H0)")

    # resonances: for each ordered pair (a <= b) and n, find sign changes of nu_a + nu_b - n Omega along E
    print("\n== predicted first-order resonances (sign changes along E of nu_a + nu_b - n Omega, n of either sign; modes matched by index after sorting by omega) ==")
    print("columns: E_res | n | modes (a, b) with sigma | type | coupling |g| (half-width in detuning 2|g|) | nearest measured block-12 band (B001)")
    ivs = b001_intervals()
    nm = min(len(r['modes']) for r in rows)
    preds = []
    for ia in range(nm):
        for ib in range(ia, nm):
            f = np.array([r['modes'][ia]['sigma'] * r['modes'][ia]['omega'] + r['modes'][ib]['sigma'] * r['modes'][ib]['omega'] for r in rows])
            Oms = np.array([r['Om'] for r in rows])
            for n in [k for k in range(-a.nharm, a.nharm + 1) if k != 0]:
                h = f - n * Oms
                for i in range(len(rows) - 1):
                    if h[i] * h[i + 1] < 0:
                        w = abs(h[i]) / (abs(h[i]) + abs(h[i + 1]))
                        Er = np.exp((1 - w) * np.log(rows[i]['E']) + w * np.log(rows[i + 1]['E']))
                        ma, mb = rows[i]['modes'][ia], rows[i]['modes'][ib]
                        Hn_ = rows[i]['Hn'][abs(n) - 1]
                        g = coupling(ma, mb, Hn_ if n > 0 else np.conj(Hn_))
                        typ = 'parametric' if ia == ib else ('combination, opposite sigma' if ma['sigma'] != mb['sigma'] else 'combination, same sigma')
                        best = min(ivs, key=lambda iv: abs(np.log(Er / np.sqrt(iv[0] * iv[1]))))
                        inside = best[0] <= Er <= best[1]
                        preds.append((Er, n, ia, ib, typ, g, best, inside))
    preds.sort()
    for Er, n, ia, ib, typ, g, best, inside in preds:
        print(f"  E={Er:.4f} | n={n:2d} | ({ia}, {ib}) sigma ({'+' if rows[0]['modes'][ia]['sigma'] > 0 else '-'},{'+' if rows[0]['modes'][ib]['sigma'] > 0 else '-'}) | {typ:28s} | |g| = {g:.4f} | [{best[0]:.4f}, {best[1]:.4f}] {best[2]} {'<- inside' if inside else ''}")
    print(f"  {sum(1 for p in preds if p[7])} of {len(preds)} predicted resonances fall inside a measured block-12 band; {len(ivs)} measured bands")

    np.savez(os.path.join(RES, 'rh27_averaging.npz'), E=grid, T=np.array([r['T'] for r in rows]),
             nu=np.array([[m['sigma'] * m['omega'] for m in r['modes'][:nm]] for r in rows]),
             preds=np.array([(p[0], p[1], p[2], p[3], p[5]) for p in preds]))
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(9, 5))
        for lo, hi, ty, mu in ivs:
            ax.axvspan(lo, hi, color=('tab:red' if ty == 'CPLX' else 'tab:blue'), alpha=0.15, lw=0)
        for ia in range(nm):
            ax.plot(grid, [abs(r['modes'][ia]['omega']) for r in rows], '.', ms=3, label=f'averaged mode {ia} (|omega|)')
        ax.plot(grid, 0.774 * np.ones_like(grid), 'k--', lw=1, label='ladder: 0.774')
        ax.set_xscale('log')
        ax.set_xlabel('E')
        ax.set_ylabel('frequency')
        ax.legend(fontsize=7)
        ax.set_title('rabbit holes 2+7: frequencies of the averaged block-12 system vs the measured ladder frequency; bands shaded', fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, 'rh27_averaging.png'), dpi=130)
        print("figure: results/rh27_averaging.png")
    except Exception as e:                                            # noqa: BLE001
        print(f"(no figure: {e})")
    print(f"DONE {time.time() - t0:.0f}")


if __name__ == '__main__':
    main()
