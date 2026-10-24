"""rh_reduced_block.py -- the gauge-reduced block 12 (the reviewer's construction, numerical first cut), verified on
the Floquet side before anything is averaged.

Block 12 in the state x = (u12, u21, u12', u21', u12'', u21'', u12''', u21''') (8-dim) carries ONE linearised Gauss
row G(t) x = 0 (rung1_system.pkl: the third component of G touches this block only), first class.  With Omega(t) the
Ostrogradsky symplectic form (krein_b005.build_omega) the gauge direction is the Hamiltonian vector field of G,
v_G = Omega^{-1} g(t), g the row of G.  Constant gauge fixing (a fixed functional f.x = 0 with f.v_G != 0 along the
orbit; reduced_K / run / admissibility below) is NOT admissible on the DR orbit: v_G(t) sweeps the whole 8-space over
a period and the best fixed f has min |f.v_G| ~ 1e-3 somewhere, where the gauge-fixed flow blows up (a multiplier
4e8 in the first attempt).  Used instead (run_moving, run_4dim): the moving gauge f(t) = v_G(t) with the connection
term from the exact orbit derivatives of Omega and g, x' = K x + lambda(t) v_G, lambda = -(v_G'.x + v_G.K x)/|v_G|^2
(the added term is along the null direction of the reduced bracket and changes no physics); off-surface components
then conserved to 1e-11, surface invariant to 1e-14.  Omega keeps c_j symbolic (the c_j = 1 form gave a wrong gauge
direction at c_j != 1).

CHECK A: the ghost multipliers of the 6-dim reduced monodromy equal the four ghost multipliers of the full block (B001's
drop-4 rule) -- this pins the construction, the admissibility and any sign error, as the reviewer prescribed.
CHECK B: the remaining 2 of the 6 reduced dims are the unipotent pair at 1 (the global colour-rotation zero mode
z = (-y, x, ...) of the planar background and its conserved charge Omega z); the second reduction by that symmetry
(run_4dim) gives the 4-dim ghost block, and CHECK A4: its spectrum IS the four ghost multipliers, to <= 8e-11 at
every E and c_j in {0.25, 0.5, 1}, including the c_j = 0.25, E = 0.6 double +1 fold (|mu| = 34, 16.5).
Then, on the 4-dim block: the reduced symplectic form (nondegenerate; monodromy symplectic for it to 3e-13), the
frozen spectrum of the restricted generator along the orbit (the gauge terms make it not exactly Hamiltonian, real
parts <= 0.1 where elliptic; read the imaginary parts as instantaneous frequencies: a light mode soft at the TURNING
POINTS, 0.13-0.7 at t = 0, hyperbolic on the arms around them, and a heavy mode at 1.0-1.7), and the averaged
generator in a parallel-transported frame -- NOT Floquet-accurate (omega T mod 2 pi off by up to a half-turn); do not
quote its frequencies as predictions.  The low-c_j static instability is not negative mean stiffness (averaged
generator purely oscillatory at c_j = 0.25 while the exact multipliers reach 34).

Usage (from the repo root): python rabbit_holes/ym/code/rh_reduced_block.py [--E 0.1 0.3 0.6] [--cj 1.0]
Cost: seconds per energy.  Added 2026-10-22 (numerical; the Ostrogradsky form, the Gauss row and their orbit
derivatives are exact sympy objects, the reduction and integration float).  Logs results/rh_reduced_block_cj*.txt,
spectra results/rh_reduced_block_E*_cj*.npz.  Status (notes/RH_reduced_block.md, README "Open"): item (i), the
hyperbolic fraction per (E, c_j), done 23 Oct in rh_reduced_block_hyp.py; (ii)-(iv) closed 24 Oct in
rh_reduced_block_transfer.py and rh_krein_width_cj.py (notes/RH_reduced_block_transfer.md): the arm/bounce split is a
frame convention (a Darboux frame is defined up to Sp(4)), the Krein width is not first order in c_j, and a symplectic
frame with an exactly Hamiltonian generator exists -- in which the frozen spectrum printed here changes (the frame
here is orthonormal, not symplectic; quote hyperbolic fractions with the frame named).
"""
import argparse
import os
import pickle
import sys

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402
import krein_b005 as kb                      # noqa: E402

IDX = base.block_idx('12')


def block_funcs():
    """8x8 lambdified K, Omega and the 8-vector g for block 12 only (the full 24x24 Omega per RHS call was the bottleneck)."""
    d = pickle.load(open(base.PKL, 'rb'))
    args = list(d['X']) + list(d['P'])
    Kf8 = sp.lambdify(args + [d['cj']], d['K'][IDX, IDX], 'numpy')
    # Omega: rebuild via krein_b005.build_omega's construction but keep only the block (lambdify the 8x8 slice)
    L2, X, P, cj = d['L2'], d['X'], d['P'], d['cj']            # c_j kept symbolic: Omega depends on it (the c_j = 1 Omega gave a wrong gauge direction at c_j != 1)
    U = [[sp.Symbol(f'u{k}_{j}') for j in range(6)] for k in range(4)]
    state = [U[k][j] for k in range(4) for j in range(6)]
    Vs = (X[0] ** 2 * X[1] ** 2 + X[1] ** 2 * X[2] ** 2 + X[2] ** 2 * X[0] ** 2) / 2
    bg_rate = {X[i]: P[i] for i in range(3)}
    bg_rate.update({P[i]: -Vs.diff(X[i]) for i in range(3)})

    def ddt(expr):
        out = 0
        for k in range(3):
            for j in range(6):
                out += expr.diff(U[k][j]) * U[k + 1][j]
        for s_, r in bg_rate.items():
            out += expr.diff(s_) * r
        return sp.expand(out)
    Q1 = [U[0][j] for j in range(6)]
    Q2 = [U[1][j] for j in range(6)]
    P2 = [sp.expand(L2.diff(U[2][j])) for j in range(6)]
    P1 = [sp.expand(L2.diff(U[1][j]) - ddt(P2[j])) for j in range(6)]
    canon = Q1 + Q2 + P1 + P2
    Jm = sp.Matrix(24, 24, lambda r, c: canon[r].coeff(state[c]))
    Jcan = sp.zeros(24, 24)
    for i in range(12):
        Jcan[i, 12 + i] = 1
        Jcan[12 + i, i] = -1
    Omega = (Jm.T * Jcan * Jm).applyfunc(sp.expand)
    Om8 = Omega[IDX, IDX]
    Of8 = sp.lambdify(args + [cj], Om8, 'numpy')
    G = d['G'][:, IDX]
    rows = [r for r in range(G.rows) if any(G[r, c] != 0 for c in range(8))]
    assert len(rows) == 1, rows
    g8 = G[rows[0], :]
    Gf8 = sp.lambdify(args + [cj], g8, 'numpy')
    # orbit derivatives (d/dt along the ghost-free background) of Omega and g, exact in sympy
    def dorb(expr):
        out = 0
        for s_, r in bg_rate.items():
            out += expr.diff(s_) * r
        return sp.expand(out)
    Odot8 = sp.lambdify(args + [cj], Om8.applyfunc(dorb), 'numpy')
    gdot8 = sp.lambdify(args + [cj], g8.applyfunc(dorb), 'numpy')
    return Kf8, Of8, Gf8, Odot8, gdot8


def load_G():
    d = pickle.load(open(base.PKL, 'rb'))
    G = d['G'][:, IDX]
    rows = [r for r in range(G.rows) if any(G[r, c] != 0 for c in range(8))]
    assert len(rows) == 1, rows
    return sp.lambdify(list(d['X']) + list(d['P']) + [d['cj']], G[rows[0], :], 'numpy')


def reduced_K(K, Om, g, f):
    """6x6 reduced generator on the surface {g.x = 0, f.x = 0}, in an orthonormal basis B (8x6) of that surface.
    x = B y;  x' = K x + lam v_G with f.(K x + lam v_G) = 0 -> lam = -(f.K x)/(f.v_G);  y' = B^T x'.
    (g.x' = 0 automatically on g.x = 0 for a first-class constraint: checked separately.)"""
    vG = np.linalg.solve(Om, g)
    A = np.vstack([g, f])
    # basis of the null space of A
    _, s, Vt = np.linalg.svd(A)
    B = Vt[2:].T                                              # 8 x 6
    fv = f @ vG
    Kg = K - np.outer(vG, f @ K) / fv
    return B.T @ Kg @ B, B, vG, fv


def admissibility(E, cj, Kf, Of, Gf, f, nsamp=128):
    q0, p0, T = base.orbit(E)
    te = np.linspace(0, T, nsamp, endpoint=False)
    sol = solve_ivp(lambda t, s: np.concatenate([s[3:], -base.gradV(s[:3])]), (0.0, T), np.concatenate([q0, p0]), t_eval=te, **base.RT)
    adm = np.inf
    for i in range(nsamp):
        Om = np.array(Of(*sol.y[:3, i], *sol.y[3:, i], cj), dtype=float)
        g = np.array(Gf(*sol.y[:3, i], *sol.y[3:, i], cj), dtype=float).ravel()
        vG = np.linalg.solve(Om, g)
        adm = min(adm, abs(f @ vG) / (np.linalg.norm(f) * np.linalg.norm(vG)))
    return adm


def run(E, cj, Kf, Of, Gf, f, nsamp=512):
    q0, p0, T = base.orbit(E)

    def bg(t, s):
        return np.concatenate([s[3:], -base.gradV(s[:3])])
    te = np.linspace(0, T, nsamp, endpoint=False)
    sol = solve_ivp(bg, (0.0, T), np.concatenate([q0, p0]), t_eval=te, **base.RT)
    Q, Pm = sol.y[:3].T, sol.y[3:].T
    # admissibility and first-class check along the orbit
    adm, fc = np.inf, 0.0
    Kred_t = []
    for i in range(nsamp):
        K = np.array(Kf(*Q[i], *Pm[i], cj), dtype=float)
        Om = np.array(Of(*Q[i], *Pm[i], cj), dtype=float)
        g = np.array(Gf(*Q[i], *Pm[i], cj), dtype=float).ravel()
        Kr, B, vG, fv = reduced_K(K, Om, g, f)
        adm = min(adm, abs(fv) / (np.linalg.norm(f) * np.linalg.norm(vG)))
        Kred_t.append(Kr)
    # first-class: d/dt (g.x) = (g' + g K).x must vanish on g.x = 0 -> (g' + gK) parallel to g.  Check numerically at a few points
    # (g' by finite difference along the orbit)
    for i in (0, nsamp // 4, nsamp // 2):
        j = (i + 1) % nsamp
        K = np.array(Kf(*Q[i], *Pm[i], cj), dtype=float)
        g = np.array(Gf(*Q[i], *Pm[i], cj), dtype=float).ravel()
        g2 = np.array(Gf(*Q[j], *Pm[j], cj), dtype=float).ravel()
        gd = (g2 - g) / (T / nsamp) + g @ K
        fc = max(fc, np.linalg.norm(gd - (gd @ g) / (g @ g) * g) / np.linalg.norm(gd))

    # reduced monodromy: integrate y' = K_red(t) y with K_red interpolated... better: integrate the full x' = K x + lam v_G
    # and read off in the basis; simplest exact route: integrate the 8-dim gauge-fixed generator Kg(t) directly
    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, cj), dtype=float)
        Om = np.array(Of(*q, *p, cj), dtype=float)
        g = np.array(Gf(*q, *p, cj), dtype=float).ravel()
        vG = np.linalg.solve(Om, g)
        Kg = K - np.outer(vG, f @ K) / (f @ vG)
        return np.concatenate([p, -base.gradV(q), (Kg @ s[6:].reshape(8, 8)).ravel()])
    sT = solve_ivp(rhs, (0.0, T), np.concatenate([q0, p0, np.eye(8).ravel()]), **base.RT).y[:, -1]
    Mg = sT[6:].reshape(8, 8)
    # restrict to the surface at t = 0 (the surface is time-dependent through g(t), but g(T) = g(0))
    K0 = np.array(Kf(*q0, *p0, cj), dtype=float)
    Om0 = np.array(Of(*q0, *p0, cj), dtype=float)
    g0 = np.array(Gf(*q0, *p0, cj), dtype=float).ravel()
    _, B, _, _ = reduced_K(K0, Om0, g0, f)
    Mred = B.T @ Mg @ B
    leak = np.linalg.norm(Mg @ B - B @ Mred) / np.linalg.norm(Mg @ B)      # does the surface map to itself?
    # full block monodromy for comparison (8-dim, same block functions)
    def rhs_full(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, cj), dtype=float)
        return np.concatenate([p, -base.gradV(q), (K @ s[6:].reshape(8, 8)).ravel()])
    sF = solve_ivp(rhs_full, (0.0, T), np.concatenate([q0, p0, np.eye(8).ravel()]), **base.RT).y[:, -1]
    ghost_full = base.ghost_multipliers(sF[6:].reshape(8, 8))
    mu_red = np.linalg.eigvals(Mred)
    # reduced symplectic form and Hamiltonian along the orbit
    Om_red0 = B.T @ Om0 @ B
    return dict(T=T, adm=adm, fc=fc, leak=leak, mu_red=mu_red, ghost_full=ghost_full, B=B, Om_red0=Om_red0, Kred_t=np.array(Kred_t), te=te, Q=Q, Pm=Pm)


def gauge_fields(Kf, Of, Gf, Odot, gdot, q, p, cj):
    K = np.array(Kf(*q, *p, cj), dtype=float)
    Om = np.array(Of(*q, *p, cj), dtype=float)
    g = np.array(Gf(*q, *p, cj), dtype=float).ravel()
    Od = np.array(Odot(*q, *p, cj), dtype=float)
    gd = np.array(gdot(*q, *p, cj), dtype=float).ravel()
    vG = np.linalg.solve(Om, g)
    vGd = np.linalg.solve(Om, gd - Od @ vG)                   # d/dt of Omega^-1 g from Omega v = g
    return K, Om, g, vG, vGd


def run_moving(E, cj, Kf, Of, Gf, Odot, gdot):
    """moving gauge fixing f(t) = v_G(t): the surface S(t) = {g.x = 0, v_G.x = 0}; x' = K x + lam v_G with
    lam = -(v_G' . x + v_G . K x)/|v_G|^2 keeps v_G.x = 0 (g.x = 0 is automatic).  Off-surface components g.x and
    v_G.x are then conserved, so the 8-dim gauge-fixed flow has no spurious growth and its restriction to S(0) is the
    reduced 6-dim monodromy."""
    q0, p0, T = base.orbit(E)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K, Om, g, vG, vGd = gauge_fields(Kf, Of, Gf, Odot, gdot, q, p, cj)
        X = s[6:].reshape(8, 8)
        lam = -(vGd @ X + vG @ (K @ X)) / (vG @ vG)             # row vector over the 8 columns
        return np.concatenate([p, -base.gradV(q), (K @ X + np.outer(vG, lam)).ravel()])
    sT = solve_ivp(rhs, (0.0, T), np.concatenate([q0, p0, np.eye(8).ravel()]), **base.RT)
    Mg = sT.y[6:, -1].reshape(8, 8)
    K0, Om0, g0, vG0, _ = gauge_fields(Kf, Of, Gf, Odot, gdot, q0, p0, cj)
    A = np.vstack([g0, vG0])
    _, _, Vt = np.linalg.svd(A)
    B = Vt[2:].T
    Mred = B.T @ Mg @ B
    leak = np.linalg.norm(Mg @ B - B @ Mred) / np.linalg.norm(Mg @ B)
    # conservation of the off-surface components: g.Mg x = g.x and vG.Mg x = vG.x for all x
    cons = max(np.linalg.norm(g0 @ Mg - g0) / np.linalg.norm(g0), np.linalg.norm(vG0 @ Mg - vG0) / np.linalg.norm(vG0))

    def rhs_full(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, cj), dtype=float)
        return np.concatenate([p, -base.gradV(q), (K @ s[6:].reshape(8, 8)).ravel()])
    sF = solve_ivp(rhs_full, (0.0, T), np.concatenate([q0, p0, np.eye(8).ravel()]), **base.RT).y[:, -1]
    ghost_full = base.ghost_multipliers(sF[6:].reshape(8, 8))
    return dict(T=T, leak=leak, cons=cons, mu_red=np.linalg.eigvals(Mred), ghost_full=ghost_full, B=B, Om_red0=B.T @ Om0 @ B)


def rotation_mode(q, p):
    """the global colour-rotation zero mode (about colour axis 3) applied to the diagonal background: (u12, u21) = (-y, x)
    and its time derivatives along the ghost-free orbit -- an exact solution of x' = K x, Gauss-clean."""
    x, y, z = q
    gV = base.gradV(q)
    H = np.array([[y * y + z * z, 2 * x * y, 2 * x * z], [2 * x * y, z * z + x * x, 2 * y * z], [2 * x * z, 2 * y * z, x * x + y * y]])
    q2 = -gV
    q3 = -(H @ p)
    return np.array([-y, x, -p[1], p[0], -q2[1], q2[0], -q3[1], q3[0]])


def run_4dim(E, cj, Kf, Of, Gf, Odot, gdot, nsamp=400):
    """second reduction, by the colour-rotation symmetry: S2(t) = {g.x = 0, v_G.x = 0, z.x = 0, c.x = 0} with z(t) the
    rotation zero mode made orthogonal to v_G and c(t) = Omega z (its conserved charge), dynamics x' = K x + lam v_G +
    lam2 z keeping both gauge conditions.  Returns the 4-dim reduced monodromy on S2(0), the frozen spectrum of the
    restricted generator along the orbit (instantaneous eigenfrequencies), and the time-average of the restricted
    generator in a parallel-transported orthonormal frame of S2(t) (the averaged block; frame connection by finite
    differences, reported)."""
    q0, p0, T = base.orbit(E)

    def fields(q, p):
        K, Om, g, vG, vGd = gauge_fields(Kf, Of, Gf, Odot, gdot, q, p, cj)
        z = rotation_mode(q, p)
        zd = K @ z                                               # exact: z is a solution
        a = (vG @ z) / (vG @ vG)
        zt = z - a * vG                                          # on the surface v_G.x = 0
        # d/dt of zt: zd - a' vG - a vGd, with a' from the quotient rule
        ad = ((vGd @ z + vG @ zd) * (vG @ vG) - (vG @ z) * 2 * (vGd @ vG)) / (vG @ vG) ** 2
        ztd = zd - ad * vG - a * vGd
        c = Om @ z
        return K, Om, g, vG, vGd, zt, ztd, c

    def Keff(K, vG, vGd, zt, ztd):
        # x' = K x + vG lam(x) + zt lam2(x);  lam keeps vG.x = 0, lam2 keeps zt.x = 0 (solve the 2x2 system)
        A = np.array([[vG @ vG, vG @ zt], [zt @ vG, zt @ zt]])
        R = np.vstack([vGd + vG @ K, ztd + zt @ K])              # rows: d/dt of (vG.x), (zt.x) before correction
        L = -np.linalg.solve(A, R)                               # 2 x 8: lam rows
        return K + np.outer(vG, L[0]) + np.outer(zt, L[1])

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K, Om, g, vG, vGd, zt, ztd, c = fields(q, p)
        X = s[6:].reshape(8, 8)
        return np.concatenate([p, -base.gradV(q), (Keff(K, vG, vGd, zt, ztd) @ X).ravel()])
    sol = solve_ivp(rhs, (0.0, T), np.concatenate([q0, p0, np.eye(8).ravel()]), t_eval=np.linspace(0, T, nsamp), **base.RT)
    Mg = sol.y[6:, -1].reshape(8, 8)
    K0, Om0, g0, vG0, _, zt0, _, c0 = fields(q0, p0)
    A0 = np.vstack([g0, vG0, zt0, c0])
    _, sv, Vt = np.linalg.svd(A0)
    B0 = Vt[4:].T                                                # 8 x 4
    Mred = B0.T @ Mg @ B0
    leak = np.linalg.norm(Mg @ B0 - B0 @ Mred) / np.linalg.norm(Mg @ B0)
    cons = max(np.linalg.norm(w @ Mg - w) / np.linalg.norm(w) for w in (g0, vG0, c0))
    # frozen spectrum of the restricted generator along the orbit, and the averaged generator in a transported frame
    inst, Kbar, nfr = [], np.zeros((4, 4)), 0
    Bprev = None
    for i in range(nsamp):
        q, p = sol.y[:3, i], sol.y[3:6, i]
        K, Om, g, vG, vGd, zt, ztd, c = fields(q, p)
        Ke = Keff(K, vG, vGd, zt, ztd)
        A = np.vstack([g, vG, zt, c])
        _, _, Vt = np.linalg.svd(A)
        B = Vt[4:].T
        if Bprev is not None:                                    # align the frame with the previous one (parallel transport by polar decomposition)
            U, _, Wt = np.linalg.svd(Bprev.T @ B)
            B = B @ (U @ Wt).T
        Kr = B.T @ Ke @ B
        lam = np.linalg.eigvals(Kr)
        inst.append(np.sort_complex(lam))
        Omr = B.T @ Om @ B
        if i < nsamp - 1:
            Kbar += Kr
            nfr += 1
        Bprev = B
    Kbar /= nfr
    return dict(T=T, leak=leak, cons=cons, mu_red=np.linalg.eigvals(Mred), inst=np.array(inst), te=sol.t, Kbar=Kbar,
                avg_eig=np.linalg.eigvals(Kbar), B0=B0, Om_red0=B0.T @ Om0 @ B0, Mred=Mred)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--E', type=float, nargs='+', default=[0.1, 0.3, 0.6])
    ap.add_argument('--cj', type=float, default=1.0)
    a = ap.parse_args()
    Kf, Of, Gf, Odot, gdot = block_funcs()
    # gauge-fixing functional: try the unit vectors of the state and keep the one with the best admissibility
    print("rh_reduced_block.py: gauge-reduced block 12, numerical first cut.  Constant gauge fixing is NOT admissible on the DR orbit")
    print("  (the unit gauge direction v_G(t) sweeps the whole 8-space: singular values of its orbit samples 0.74/0.39/.../0.09 at E=0.1,")
    print("  and the best fixed f has min |f.v_G| ~ 1e-3), so the moving gauge f(t) = v_G(t) is used, with the connection term from the")
    print("  exact orbit derivatives of Omega and g (sympy).")
    for E in a.E:
        r = run_moving(E, a.cj, Kf, Of, Gf, Odot, gdot)
        print(f"\n== E = {E}, c_j = {a.cj}, T = {r['T']:.3f}: first reduction (Gauss + gauge), 6-dim ==")
        print(f"  CHECK off-surface components conserved over one period: {r['cons']:.1e} (tol 1e-7) {'PASS' if r['cons'] < 1e-7 else 'FAIL'}; surface invariant: leak {r['leak']:.1e} {'PASS' if r['leak'] < 1e-7 else 'FAIL'}")
        mr = r['mu_red'][np.argsort(-np.abs(r['mu_red']))]
        worst = max(min(abs(mr - m)) for m in r['ghost_full'])
        print(f"  CHECK A: the four ghost multipliers of the full block are reduced multipliers, worst distance {worst:.1e} {'PASS' if worst < 1e-6 else 'FAIL'}; the other two: {np.round([m for m in mr if min(abs(m - r['ghost_full'])) > 1e-6], 6)}")
        r4 = run_4dim(E, a.cj, Kf, Of, Gf, Odot, gdot)
        print(f"  -- second reduction (colour rotation + its charge), 4-dim --")
        print(f"  CHECK conserved off-surface functionals {r4['cons']:.1e}, surface invariant leak {r4['leak']:.1e} (tol 1e-7) {'PASS' if max(r4['cons'], r4['leak']) < 1e-7 else 'FAIL'}")
        m4 = r4['mu_red'][np.argsort(-np.abs(r4['mu_red']))]
        worst4 = max(min(abs(m4 - m)) for m in r['ghost_full'])
        print(f"  CHECK A4: 4-dim reduced multipliers {np.round(m4, 6)} = the ghost multipliers to {worst4:.1e} {'PASS' if worst4 < 1e-6 else 'FAIL'}")
        Omr = r4['Om_red0']
        print(f"  reduced symplectic form (4x4) at t=0: antisymmetric to {np.abs(Omr + Omr.T).max():.1e}, rank {np.linalg.matrix_rank(Omr)}; Mred^T Om Mred = Om to {np.abs(r4['Mred'].T @ Omr @ r4['Mred'] - Omr).max() / np.abs(Omr).max():.1e}")
        inst = r4['inst']
        print(f"  instantaneous (frozen) eigenvalues of the restricted 4-dim generator along the orbit, t/T = 0, 1/8, ..., 7/8:")
        for i in range(0, len(r4['te']), len(r4['te']) // 8):
            print(f"     t/T={r4['te'][i] / r4['T']:.3f}: {np.round(inst[i], 4)}")
        print(f"  fraction of the period with a real (hyperbolic) frozen eigenvalue: {np.mean([np.any(np.abs(l.real) > 1e-6) for l in inst]):.2f}")
        print(f"  averaged generator (transported frame) eigenvalues: {np.round(np.sort_complex(r4['avg_eig']), 4)}   (frequencies = imaginary parts; compare the measured lines 0.85-1.2 and the quasi-static one)")
        np.savez(os.path.join(HERE, '..', 'results', f'rh_reduced_block_E{E:g}_cj{a.cj:g}.npz'), te=r4['te'], inst=inst, Kbar=r4['Kbar'], mu_red=m4, Om_red0=Omr)


if __name__ == '__main__':
    main()
