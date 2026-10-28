"""rh_reduced_block_transfer.py -- open items (ii) and (iv) of RH_reduced_block.md on the exact gauge-reduced 4-dim
block 12: the transfer matrix across one hyperbolic arm and across one bounce, how much of the arm's expansion the
bounce returns, the sector mixing per half period against Delta(E), and a symplectic (Darboux) frame along the orbit
in which the frozen generator is exactly Hamiltonian.

Segments.  The light reduced mode is hyperbolic on the arms around the turning points (t = 0, T/2) for |t/T| <= 0.13
at E >= 0.3 (results/rh_reduced_block_hyp.txt) and never at the bounce (t = T/4).  So one period is cut into
ARM1 = [-aT, aT], BOUNCE1 = [aT, (1/2 - a)T], ARM2 = [(1/2 - a)T, (1/2 + a)T], BOUNCE2 = [(1/2 + a)T, (1 - a)T],
a = 0.13 by default (--a; the quarter-period cut of rh3_segments.py is a = 0.125).  The moving-gauge 8-dim propagator
of rh_reduced_block.run_4dim is integrated from -aT over one period and restricted to the surface S(t) at every
boundary; each segment map is expressed between Darboux frames of Omega_red (symplectic Gram-Schmidt, as in
rh3_segments.py), where its singular values are the expansion of a symplectic map in that frame (frame-DEPENDENT: see the REFRAME check
and notes/RH_reduced_block_transfer.md).

Per (E, c_j), printed:
  * the arm's and the bounce's largest singular value (ln gain), the gain of the composed half period ARM then BOUNCE,
    and the RETURN: the arm's most-stretched output direction pushed through the bounce -- the fraction of the arm's
    ln-stretch that survives (negative = the bounce undoes it);
  * the half-period map N = St Phi_{T/2} on S(-aT) (St = diag(+1, -1) x 4 the reflection of
    islands/code/dr_offdiag_halfperiod.py): CHECK N^2 has the ghost multipliers; its eigenangles and
    Delta = (phi_A + phi_B) mod 2pi (CHECK against rh5_splitting's convention: Delta = 0 means +-sqrt(ab) antipodal pairs);
  * the SECTOR MIXING of the arm, the bounce and the half period: with P+- = (1 +- St)/2 the reflection projectors
    and W+- = B^T P+- B their restrictions to the surface, the sector-diagonal part of a surface map A is
    A_d = W+ A W+ + W- A W- and the mixing is |A - A_d|_F / |A|_F (Euclidean, as the sector fractions f+ of
    RH6_polarisation.md); RH27_averaging.md part E conjectured Delta(E) to be this mixing of Phi_{T/2}.
  * item (iv): a continuous Darboux frame D(t) of Omega_red(t) (symplectic Gram-Schmidt seeded by the previous
    frame on a fine grid, on top of the polar-transported orthonormal frame of run_4dim) and the generator in it,
    K_z = D^-1 (K_red D - dD/dt): CHECK J K_z is symmetric (exactly Hamiltonian up to the finite-difference dD/dt);
    its frozen spectrum's hyperbolic fraction against the orthonormal frame's.  A symplectic frame is not unique, so
    the frozen spectrum in it is one choice among many; the segment singular values above are not frame-free either (REFRAME); only N's spectrum is.

Usage (from the repo root):
  python rabbit_holes/ym/code/rh_reduced_block_transfer.py [--E 0.1 0.156 0.3 0.412 0.6 0.642] [--cj 1.0] [--a 0.13]
Cost: ~10 s per (E, c_j).  Writes results/rh_reduced_block_transfer_cj<cj>_a<a>.npz (segment data and the two frozen
spectra at 2000 samples per period).  Added 2026-10-24.
"""
import argparse
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402
import rh_reduced_block as RB                # noqa: E402

RES = os.path.join(HERE, '..', 'results')
ST = np.diag([1.0, -1.0] * 4)                # reflection y -> -y on block 12: u12 even, u21 odd, all derivative levels
PP, PM = (np.eye(8) + ST) / 2, (np.eye(8) - ST) / 2
J4 = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])


def darboux(Om, seed=None):
    """D with D^T Om D = J4 by symplectic Gram-Schmidt; seeded by the columns of a previous frame for continuity."""
    basis = list((seed if seed is not None else np.eye(4)).T)
    cols = []
    for k in range(2):
        e = basis.pop(0)
        vals = [abs(e @ Om @ b) for b in basis]
        f = basis.pop(int(np.argmax(vals)))
        f = f / (e @ Om @ f)
        cols.append((e, f))
        basis = [b - (b @ Om @ f) * e + (b @ Om @ e) * f for b in basis]
    D = np.zeros((4, 4))
    for k, (e, f) in enumerate(cols):
        D[:, k], D[:, 2 + k] = e, f
    return D


def make_fields(Kf, Of, Gf, Odot, gdot, cj):
    def fields(q, p):
        K, Om, g, vG, vGd = RB.gauge_fields(Kf, Of, Gf, Odot, gdot, q, p, cj)
        z = RB.rotation_mode(q, p)
        zd = K @ z
        al = (vG @ z) / (vG @ vG)
        zt = z - al * vG
        ad = ((vGd @ z + vG @ zd) * (vG @ vG) - (vG @ z) * 2 * (vGd @ vG)) / (vG @ vG) ** 2
        ztd = zd - ad * vG - al * vGd
        c = Om @ z
        A = np.array([[vG @ vG, vG @ zt], [zt @ vG, zt @ zt]])
        R = np.vstack([vGd + vG @ K, ztd + zt @ K])
        L = -np.linalg.solve(A, R)
        Ke = K + np.outer(vG, L[0]) + np.outer(zt, L[1])
        return Ke, Om, np.vstack([g, vG, zt, c])
    return fields


def surface(fields, q, p):
    Ke, Om, A = fields(q, p)
    _, _, Vt = np.linalg.svd(A)
    B = Vt[4:].T
    return Ke, Om, B


def mixing(A, B0, B1):
    """sector mixing of a surface map A: S(t0) -> S(t1) (8x8 propagator restricted through the orthonormal bases)."""
    Ar = B1.T @ A @ B0
    W0p, W0m = B0.T @ PP @ B0, B0.T @ PM @ B0
    W1p, W1m = B1.T @ PP @ B1, B1.T @ PM @ B1
    Ad = W1p @ Ar @ W0p + W1m @ Ar @ W0m
    return np.linalg.norm(Ar - Ad) / np.linalg.norm(Ar)


def run(E, cj, a, Kf, Of, Gf, Odot, gdot, nfine=2000):
    fields = make_fields(Kf, Of, Gf, Odot, gdot, cj)
    q0, p0, T = base.orbit(E)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        Ke, _, _ = fields(q, p)
        return np.concatenate([p, -base.gradV(q), (Ke @ s[6:].reshape(8, 8)).ravel()])
    # background to t0 = -aT, then one period of the gauge-fixed propagator with the segment boundaries saved
    s_bg = solve_ivp(lambda t, s: np.concatenate([s[3:], -base.gradV(s[:3])]), (0.0, (1 - a) * T), np.concatenate([q0, p0]), **base.RT).y[:, -1]
    bounds = np.array([0.0, 2 * a, 0.5, 0.5 + 2 * a, 1.0]) * T
    sol = solve_ivp(rhs, (0.0, T), np.concatenate([s_bg, np.eye(8).ravel()]), t_eval=bounds, **base.RT)
    Phis = [sol.y[6:, i].reshape(8, 8) for i in range(len(bounds))]
    fr = []
    for i in range(len(bounds)):
        Ke, Om, B = surface(fields, sol.y[:3, i], sol.y[3:6, i])
        Omr = B.T @ Om @ B
        fr.append((B, darboux(Omr), Omr))
    out = dict(T=T)
    names = ['ARM1', 'BOUNCE1', 'ARM2', 'BOUNCE2']
    segs, Tk_list, sv_list = [], [], []
    for k in range(4):
        B0, D0, _ = fr[k]
        B1, D1, _ = fr[k + 1]
        Seg = Phis[k + 1] @ np.linalg.inv(Phis[k])
        Tk = np.linalg.solve(D1, B1.T @ Seg @ B0 @ D0)
        symp = np.abs(Tk.T @ J4 @ Tk - J4).max()
        U, sv, Vt = np.linalg.svd(Tk)
        segs.append((names[k], sv, symp, mixing(Seg, B0, B1)))
        Tk_list.append(Tk)
        sv_list.append(sv)
    # the return: the arm's most-stretched output (left singular vector) through the bounce
    U1, s1, _ = np.linalg.svd(Tk_list[0])
    y = Tk_list[1] @ U1[:, 0]
    ret = np.log(np.linalg.norm(y)) / np.log(s1[0]) if s1[0] > 1 + 1e-9 else np.nan
    # REFRAME test: a Darboux frame is defined up to Sp(4), and singular values are invariant only under its unitary
    # subgroup, so the segment gains and the return depend on the frame chosen at each boundary.  Three random
    # symplectic reframes S_k = exp(J H_k) (H_k symmetric, entries ~ N(0, 1/2)) at the three boundaries of the half
    # period: the segment singular values and the return change, the spectrum of N does not.
    rng = np.random.default_rng(0)
    reframe = []
    for trial in range(3):
        Ss = []
        for _ in range(3):
            Aq = rng.normal(size=(4, 4)); Hs = (Aq + Aq.T) / 2
            Ss.append(expm(J4 @ Hs))
        Ta = np.linalg.solve(Ss[1], Tk_list[0] @ Ss[0])
        Tb = np.linalg.solve(Ss[2], Tk_list[1] @ Ss[1])
        sa, sb = np.linalg.svd(Ta, compute_uv=False)[0], np.linalg.svd(Tb, compute_uv=False)[0]
        Ua, s1a, _ = np.linalg.svd(Ta)
        reta = np.log(np.linalg.norm(Tb @ Ua[:, 0])) / np.log(s1a[0])
        reframe.append((np.log(sa), np.log(sb), reta))
    # and the arm's most-stretched direction's sector content before and after
    half = Tk_list[1] @ Tk_list[0]
    sv_half = np.linalg.svd(half, compute_uv=False)
    full = Tk_list[3] @ Tk_list[2] @ half
    mu = np.linalg.eigvals(full)
    # full-period check against the block's ghost multipliers
    def rhs_full(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, cj), dtype=float)
        return np.concatenate([p, -base.gradV(q), (K @ s[6:].reshape(8, 8)).ravel()])
    sF = solve_ivp(rhs_full, (0.0, T), np.concatenate([q0, p0, np.eye(8).ravel()]), **base.RT).y[:, -1]
    gh = base.ghost_multipliers(sF[6:].reshape(8, 8))
    worst = max(min(abs(mu - m)) for m in gh)
    # the half-period map N = St Phi_{T/2} on S(t0): St maps S(t0 + T/2) onto S(t0) (leak checked)
    B0, D0, Om0 = fr[0]
    B2, D2, Om2 = fr[2]
    Phalf = Phis[2]                                            # Phi(t0 + T/2, t0)
    Nfull = ST @ Phalf
    Nred = B0.T @ Nfull @ B0
    leakN = np.linalg.norm(Nfull @ B0 - B0 @ Nred) / np.linalg.norm(Nfull @ B0)
    N2 = Nred @ Nred
    muN2 = np.linalg.eigvals(N2)
    worstN = max(min(abs(muN2 - m)) for m in gh)
    lamN = np.linalg.eigvals(np.linalg.solve(D0, Nred @ D0))
    # Delta in rh5_splitting's convention: psi the unsigned half-period angles of the two pairs, phi = 2|psi| the
    # full-period angles, Delta = (phi_A + phi_B) mod 2pi folded to (-pi, pi]; 0 = pairs degenerate as sets, +-pi antipodal
    on_circle = np.all(np.abs(np.abs(lamN) - 1) < 1e-6)
    ang = np.sort(np.abs(np.angle(lamN)))                      # four values, two coincident pairs
    phiA, phiB = ang[0], ang[2]
    Delta = np.mod(2 * (phiA + phiB) + np.pi, 2 * np.pi) - np.pi
    mix_half = mixing(Phalf, B0, B2)
    mix_N = mixing(Nfull, B0, B0)
    out.update(reframe=reframe, segs=segs, ret=ret, sv_half=sv_half, mu=mu, gh=gh, worst=worst, leakN=leakN, worstN=worstN,
               lamN=lamN, phiA=phiA, phiB=phiB, Delta=Delta, on_circle=on_circle, mix_half=mix_half, mix_N=mix_N, absN=np.abs(lamN))

    # ---- item (iv): a continuous symplectic frame on a fine grid over one period from t = 0 ----
    # Frame: in the polar-transported orthonormal basis B(t) of run_4dim, Omega_red(t) = R Lambda R^T with two
    # invariant 2-planes of strengths w_1 >= w_2 (the eigen-planes of -Omega_red^2); e_k the unit vector of plane k
    # rotated as little as possible from the previous sample, f_k = -Omega_red e_k / w_k^2, D = [e_1, e_2, f_1, f_2]
    # gives D^T Omega_red D = J4.  The only freedom left is a rotation within each plane, fixed by minimal rotation.
    sol2 = solve_ivp(lambda t, s: np.concatenate([s[3:], -base.gradV(s[:3])]), (0.0, T), np.concatenate([q0, p0]),
                     t_eval=np.linspace(0, T, nfine + 1), **base.RT)
    Bprev, eprev = None, None
    Ds, Kr_list, wmin = [], [], np.inf
    for i in range(nfine + 1):
        q, p = sol2.y[:3, i], sol2.y[3:6, i]
        Ke, Om, B = surface(fields, q, p)
        if Bprev is not None:
            U, _, Wt = np.linalg.svd(Bprev.T @ B)
            B = B @ (U @ Wt).T
        Omr = B.T @ Om @ B
        w2, R = np.linalg.eigh(-Omr @ Omr)                     # eigenvalues w_k^2 (each twice), planes = eigenvectors
        order = np.argsort(-w2)
        w2, R = w2[order], R[:, order]
        wmin = min(wmin, np.sqrt(w2[2]) / np.sqrt(w2[0]))       # plane separation (1 = degenerate planes)
        planes = [R[:, :2], R[:, 2:]]
        es = []
        for k, Pk in enumerate(planes):
            if eprev is None:
                e = Pk[:, 0]
            else:
                e = Pk @ (Pk.T @ eprev[k])
                e = e / np.linalg.norm(e)
            es.append(e)
        wk = np.sqrt(w2[[0, 2]])
        fs = [-Omr @ es[k] / wk[k] ** 2 for k in range(2)]
        D = np.column_stack([es[0], es[1], fs[0], fs[1]])
        Ds.append(D); Kr_list.append(B.T @ Ke @ B)
        Bprev, eprev = B, es
    h = T / nfine
    Ds = np.array(Ds)
    darb = max(np.abs(Ds[i].T @ (np.linalg.inv(Ds[i]).T @ J4 @ np.linalg.inv(Ds[i])) @ Ds[i] - J4).max() for i in (0, nfine // 2))
    hyp_orth, hyp_symp, asym, inst_z, inst_o, conn = [], [], [], [], [], []
    for i in range(1, nfine):
        Dd = (Ds[i + 1] - Ds[i - 1]) / (2 * h)
        Kz = np.linalg.solve(Ds[i], Kr_list[i] @ Ds[i] - Dd)
        Hm = J4 @ Kz                                           # J K_z symmetric iff K_z is Hamiltonian
        asym.append(np.abs(Hm - Hm.T).max() / np.abs(Hm).max())
        conn.append(np.linalg.norm(Dd) / np.linalg.norm(Kr_list[i] @ Ds[i]))
        lz = np.linalg.eigvals(Kz)
        lo = np.linalg.eigvals(Kr_list[i])
        inst_z.append(np.sort_complex(lz)); inst_o.append(np.sort_complex(lo))
        hyp_symp.append(np.any((np.abs(lz.imag) < 1e-3) & (np.abs(lz.real) > 0.1)))
        hyp_orth.append(np.any((np.abs(lo.imag) < 1e-3) & (np.abs(lo.real) > 0.1)))
    inst_z, inst_o = np.array(inst_z), np.array(inst_o)
    # non-Hamiltonian fuzz: the largest |Re| among eigenvalues that are not real (a Hamiltonian 4x4 has either
    # purely imaginary pairs, real pairs, or a quadruplet; a complex pair with a small real part is frame fuzz)
    def fuzz(inst):
        v = []
        for l in inst:
            cplx = l[np.abs(l.imag) > 1e-3]
            if len(cplx) == 4 and np.abs(cplx.real).max() > 1e-3:
                # quadruplet test: real parts come in +-pairs
                rs = np.sort(cplx.real)
                if abs(rs[0] + rs[3]) < 3e-2 * np.abs(rs).max() and abs(rs[1] + rs[2]) < 3e-2 * np.abs(rs).max():
                    continue
            v.append(np.abs(cplx.real).max() if len(cplx) else 0.0)
        return max(v) if v else 0.0
    out.update(asym=np.median(asym), asym_max=np.max(asym), conn=np.median(conn), wmin=wmin, darb=darb,
               hyp_orth=np.mean(hyp_orth), hyp_symp=np.mean(hyp_symp), fuzz_orth=fuzz(inst_o), fuzz_symp=fuzz(inst_z),
               te=sol2.t[1:nfine], inst_z=inst_z, inst_o=inst_o)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--E', type=float, nargs='+', default=[0.1, 0.156, 0.3, 0.412, 0.6, 0.642])
    ap.add_argument('--cj', type=float, default=1.0)
    ap.add_argument('--a', type=float, default=0.13)
    args = ap.parse_args()
    Kf, Of, Gf, Odot, gdot = RB.block_funcs()
    print(f"rh_reduced_block_transfer.py: arm / bounce transfer matrices of the 4-dim reduced block 12, a = {args.a} "
          f"(ARM = [-aT, aT] around the turning point, BOUNCE = [aT, (1/2 - a)T]); c_j = {args.cj}")
    rows = {}
    for E in args.E:
        r = run(E, args.cj, args.a, Kf, Of, Gf, Odot, gdot)
        rows[E] = r
        print(f"\n== E = {E}, c_j = {args.cj}, T = {r['T']:.3f} ==")
        for name, sv, symp, mx in r['segs']:
            print(f"  {name:8s}: ln max sv {np.log(sv[0]):+.3f}  all {np.round(sv, 3)}  symplectic to {symp:.1e}  sector mixing {mx:.3f}")
        print(f"  CHECK product of the four segments = ghost multipliers {np.round(r['gh'], 5)}: worst {r['worst']:.1e} {'PASS' if r['worst'] < 1e-6 else 'FAIL'}")
        print(f"  half period ARM then BOUNCE: ln max sv {np.log(r['sv_half'][0]):+.3f} against ARM + BOUNCE {np.log(r['segs'][0][1][0]) + np.log(r['segs'][1][1][0]):+.3f}")
        print(f"  RETURN: the arm's most-stretched output through the bounce keeps {r['ret']:+.2f} of its ln-stretch "
              f"({'the bounce undoes it' if r['ret'] < 0 else 'the bounce keeps part of it'})")
        print("  REFRAME (three random symplectic frames at the boundaries): arm ln sv / bounce ln sv / return = " +
              "; ".join(f"{x:+.2f} / {y:+.2f} / {z:+.2f}" for x, y, z in r['reframe']) + "  -- frame-dependent (the spectrum of N below is not)")
        print(f"  N = St Phi_T/2 on S(-aT): St maps the surface to itself to {r['leakN']:.1e}; CHECK N^2 = ghost multipliers to {r['worstN']:.1e} "
              f"{'PASS' if r['worstN'] < 1e-6 else 'FAIL'}; |lambda_N| {np.round(r['absN'], 4)}")
        print(f"  N half-period angles |psi| {r['phiA'] / np.pi:.3f} pi, {r['phiB'] / np.pi:.3f} pi -> Delta = {r['Delta'] / np.pi:+.3f} pi "
              f"({'on circle' if r['on_circle'] else 'OFF circle, Delta not meaningful'}; rh5_splitting convention); sector mixing of Phi_T/2 {r['mix_half']:.3f}, of N {r['mix_N']:.3f}")
        print(f"  (iv) symplectic frame: D^T Om D = J to {r['darb']:.1e}; plane separation min w2/w1 {r['wmin']:.3f}; connection |dD/dt|/|K D| median {r['conn']:.2f}; "
              f"J K_z symmetric to {r['asym']:.1e} (median; max {r['asym_max']:.1e})")
        print(f"       hyperbolic fraction: orthonormal frame {r['hyp_orth']:.2f} -> symplectic frame {r['hyp_symp']:.2f}; non-Hamiltonian fuzz max |Re| of complex pairs {r['fuzz_orth']:.3f} -> {r['fuzz_symp']:.3f}")
    print("\nSummary (E | arm ln sv | bounce ln sv | half ln sv | return | Delta/pi (* = off circle) | mix Phi_T/2 | mix arm | mix bounce | hyp frac orth -> symp)")
    for E, r in rows.items():
        s = r['segs']
        print(f"  {E:5.3f} | {np.log(s[0][1][0]):+.3f} | {np.log(s[1][1][0]):+.3f} | {np.log(r['sv_half'][0]):+.3f} | {r['ret']:+.2f} | {r['Delta'] / np.pi:+.3f}{'' if r['on_circle'] else '*'} | {r['mix_half']:.3f} | {s[0][3]:.3f} | {s[1][3]:.3f} | {r['hyp_orth']:.2f} -> {r['hyp_symp']:.2f}")
    np.savez(os.path.join(RES, f'rh_reduced_block_transfer_cj{args.cj:g}_a{args.a:g}.npz'),
             E=np.array(list(rows)), **{f'E{E:g}_{k}': v for E, r in rows.items() for k, v in r.items() if isinstance(v, (np.ndarray, float, int))})


if __name__ == '__main__':
    main()
