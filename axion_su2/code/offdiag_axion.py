"""offdiag_axion.py -- the off-diagonal Lee-Wick blocks on the isotropic axion-SU(2) background with the drive.

Uses rung1_axion_system.pkl (rung1_axion.py): the 24-dim first-order system of the six off-diagonal
components on the diagonal background with L_axion = kappa(t) det A, the background exact (q, q', q'', q'''
free, q'''' from the diagonal 4th-order equations).  Three 8-dim blocks (u12,u21), (u13,u31), (u23,u32), equal
by isotropy on x = y = z.  Units g = M = 1; c_j = 1 (the theory).  The drive adds to each block exactly the
potential -kappa z u12 u21 (rung1_axion.py), nothing else.

Prediction (confirmed in Part A exactly).  On the isotropic background the k = 0 fluctuations decompose under the diagonal SO(3) into
1 + 3 + 5 (Khan, Pirzada, Mustafa, arXiv:2604.07044): the two diagonal anisotropic modes of axion_aniso.py and
the three symmetric off-diagonal combinations (u12 + u21)/sqrt2 etc. form the traceless-symmetric quintet,
the antisymmetric combinations the triplet, which at m = 0 is the gauge direction u = theta x A and its Gauss
partner.  So the physical off-diagonal pair should obey exactly the anisotropic equation
eps''''/M^2 + eps'' + kappa Q eps = 0, degenerate with it, and whatever holds for the anisotropic sector carries
over unchanged (on 15 Oct that was the sign rule of axion_twins.py, withdrawn the same evening: on the
self-consistent background the sector is marginal at second order, axion_dressed.py A1); the self-coupling's only
role in this sector is the Yang-Mills mass 2 Q^2 of the triplet, which the Gauss law removes.

Parts:
  A  algebra at random background points on x = y = z: the symmetric subspace of the (u12,u21) block is
     invariant under K and K restricted to it is the companion matrix of eps'''' = -eps'' - kappa Q eps
     (PASS/FAIL); the three blocks agree; the gauge zero mode is in the kernel of G.
  B  direct two-frequency integration on the exact 4th-order background (ghost included), phi = 1, m_phi/gQ = 100,
     6 Q periods: the 8-dim (u12,u21) block, its symmetric restriction, the anisotropic 4-dim sector alongside, the
     block on the Gauss surface modulo the gauge direction, the gauge zero mode's norm and the Gauss drift;
     omega = 0.8, 2, 3 and (u13,u31) at 2.  The first version of this part used the M -> infinity background and
     found the gauge direction growing exponentially at omega = 0.8: the gauge zero mode is exact only on a true
     solution, and that background is not one at O(1) when omega ~ 1 (the ghost's response to the drive).  That
     is what led to axion_dressed.py, which holds the corrected result for both sectors.

Float only (DOP853, rtol 1e-11); nothing here is certified.  Usage (from the repository root, after
rung1_axion.py; ~2 min on the cloud CPU):
    python axion_su2/code/offdiag_axion.py > axion_su2/results/offdiag_axion.txt
Added 2026-10-15.
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
import iso_floquet as iso                                   # noqa: E402

PKL = os.path.join(HERE, '..', 'rung1_axion_system.pkl')
RT = dict(method='DOP853', rtol=1e-11, atol=1e-13)
T0 = time.time()


def stamp(msg=''):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def load():
    d = pickle.load(open(PKL, 'rb'))
    args = list(d['X']) + list(d['P']) + list(d['X2']) + list(d['X3']) + [d['cj'], d['ksyms'][0]]
    Kb, Gb, Zb = {}, {}, {}
    for name in iso.BLOCKS:
        idx = iso.block_idx(name)
        Kb[name] = sp.lambdify(args, d['K'].extract(idx, idx), 'numpy')
        Gb[name] = sp.lambdify(args, d['G'].extract([0, 1, 2], idx), 'numpy')
        Zb[name] = sp.lambdify(args, d['Zm'].extract(idx, [0, 1, 2]), 'numpy')
    q4 = sp.lambdify(args, d['q4'], 'numpy')
    return Kb, Gb, Zb, q4


def bgargs(Q, Q1, Q2, Q3, kap, cj=1.0):
    return (Q, Q, Q, Q1, Q1, Q1, Q2, Q2, Q2, Q3, Q3, Q3, cj, kap)


SYM = np.zeros((8, 4))            # symmetric subspace of the (u_ab, u_ba) block: u_ab = u_ba at each derivative order
for k in range(4):
    SYM[2 * k, k] = SYM[2 * k + 1, k] = 1 / np.sqrt(2)


def companion(c):
    C = np.zeros((4, 4))
    C[0, 1] = C[1, 2] = C[2, 3] = 1.0
    C[3, 2] = -1.0
    C[3, 0] = -c
    return C


# ------------------------------------------------------------------ Part A
def part_A(Kb, Gb, Zb, q4):
    print("\nPart A: algebra at random isotropic background points (x = y = z = Q, kappa given)\n")
    rng = np.random.default_rng(1)
    worst_inv = worst_comp = worst_iso = worst_z = 0.0
    for _ in range(20):
        Q, Q1, Q2, Q3, kap = rng.normal(size=5) * [0.5, 0.3, 0.3, 0.3, 0.4]
        a = bgargs(Q, Q1, Q2, Q3, kap)
        K8 = np.array(Kb['12'](*a), dtype=float)
        KS = K8 @ SYM
        proj = SYM @ (SYM.T @ KS)
        worst_inv = max(worst_inv, np.max(np.abs(KS - proj)))
        worst_comp = max(worst_comp, np.max(np.abs(SYM.T @ KS - companion(kap * Q))))
        for other in ('13', '23'):
            worst_iso = max(worst_iso, np.max(np.abs(np.array(Kb[other](*a), dtype=float) - K8)))
        G8 = np.array(Gb['12'](*a), dtype=float)
        Z8 = np.array(Zb['12'](*a), dtype=float)
        worst_z = max(worst_z, np.max(np.abs(G8 @ Z8)))
    print(f"    symmetric subspace invariant under K: max residual {worst_inv:.1e}  {'PASS' if worst_inv < 1e-10 else 'FAIL'}")
    print(f"    K on it = companion of eps'''' = -eps'' - kappa Q eps: max difference {worst_comp:.1e}  {'PASS' if worst_comp < 1e-10 else 'FAIL'}")
    print(f"    the three blocks agree: {worst_iso:.1e}  {'PASS' if worst_iso < 1e-10 else 'FAIL'}")
    print(f"    gauge zero mode in ker G: {worst_z:.1e}  {'PASS' if worst_z < 1e-10 else 'FAIL'}")
    print("    So on the isotropic background the symmetric off-diagonal pair IS the anisotropic sector of")
    print("    axion_aniso.py (the quintet is one representation); the antisymmetric pair is the gauge direction.")


# ------------------------------------------------------------------ backgrounds
def drive(omega, phi_amp):
    kappa = lambda t: -omega * phi_amp * np.sin(omega * t)
    kappa1 = lambda t: -omega ** 2 * phi_amp * np.cos(omega * t)
    return kappa, kappa1


def minf_state(Q, Q1, t, kappa, kappa1):
    """(q, q', q'', q''') of the M -> infinity isotropic background Q'' = -2Q^3 + kappa Q^2."""
    Q2 = -2.0 * Q ** 3 + kappa(t) * Q ** 2
    Q3 = -6.0 * Q ** 2 * Q1 + kappa1(t) * Q ** 2 + 2.0 * kappa(t) * Q * Q1
    return Q, Q1, Q2, Q3


# ------------------------------------------------------------------ Part B
def part_B(Kb, Gb, Zb, q4):
    print("\nPart B: direct two-frequency integration on the exact 4th-order background (ghost included, dressed data),")
    print("        m_phi/gQ = 100, phi = 1 (below the isotropic threshold of axion_dressed.py at every omega here), 6 Q periods;")
    print("        per cycle: ln sigma_max of the 8-dim (u12,u21) block, of its symmetric restriction, of the anisotropic")
    print("        4-dim sector alongside, of the block on the Gauss surface modulo the gauge direction; |Phi Z0| (the gauge")
    print("        zero mode, exact: must stay O(1)); Gauss drift |G(t) Phi v|, v in ker G(0).\n")
    K = ellipk(0.5)
    phi_amp, ratio = 1.0, 100
    for name, omega in (('12', 0.8), ('12', 2.0), ('13', 2.0), ('12', 3.0)):
        Q_amp = omega / ratio
        phi_eff = phi_amp / (1 - omega ** 2)
        lam = 2 + phi_eff ** 2
        X_amp = Q_amp * (1 + 2 * phi_eff ** 2 / lam) ** 0.25
        T_Q = 4 * K / (np.sqrt(lam) * X_amp)
        kappa, kappa1 = drive(omega, phi_amp)

        def bg4(t, b):
            Q, Q1, Q2, Q3 = b
            return np.array([Q1, Q2, Q3, q4(*bgargs(Q, Q1, Q2, Q3, kappa(t)))[0]])
        b0 = np.array([Q_amp, 0.0, -2 * Q_amp ** 3, kappa1(0.0) * Q_amp ** 2 / (1 - omega ** 2)])
        n_cyc = 6

        def rhs(t, y):
            b = y[:4]
            a = bgargs(*b, kappa(t))
            K8 = np.array(Kb[name](*a), dtype=float)
            C4 = companion(kappa(t) * b[0])
            return np.concatenate([bg4(t, b), (K8 @ y[4:68].reshape(8, 8)).ravel(), (C4 @ y[68:].reshape(4, 4)).ravel()])
        y0 = np.concatenate([b0, np.eye(8).ravel(), np.eye(4).ravel()])
        t_eval = np.arange(1, n_cyc + 1) * T_Q
        sol = solve_ivp(rhs, (0.0, n_cyc * T_Q), y0, t_eval=t_eval, **RT)
        a0 = bgargs(*b0, kappa(0.0))
        G0 = np.array(Gb[name](*a0), dtype=float)
        Z0 = np.array(Zb[name](*a0), dtype=float)
        _, _, Vt = np.linalg.svd(G0)
        ker = Vt[np.linalg.matrix_rank(G0):].T
        zc = Z0[:, [np.argmax(np.linalg.norm(Z0, axis=0))]]
        zc = zc / np.linalg.norm(zc)
        km = np.linalg.qr((np.eye(8) - zc @ zc.T) @ ker)[0][:, :ker.shape[1] - 1]
        print(f"  block {name}, omega = {omega}, T_Q = {T_Q:.1f}:")
        print(f"    {'cycle':>5} {'8-d':>7} {'sym':>7} {'aniso':>7} {'constr':>7} {'|Phi Z0|':>9} {'Gauss':>9}")
        for k in range(n_cyc):
            P8 = sol.y[4:68, k].reshape(8, 8)
            P4 = sol.y[68:, k].reshape(4, 4)
            bk = sol.y[:4, k]
            G1 = np.array(Gb[name](*bgargs(*bk, kappa(t_eval[k]))), dtype=float)
            print(f"    {k + 1:5d} {np.log(np.linalg.svd(P8, compute_uv=False)[0]):7.2f} "
                  f"{np.log(np.linalg.svd(SYM.T @ P8 @ SYM, compute_uv=False)[0]):7.2f} "
                  f"{np.log(np.linalg.svd(P4, compute_uv=False)[0]):7.2f} "
                  f"{np.log(np.linalg.svd(P8 @ km, compute_uv=False)[0]):7.2f} "
                  f"{np.linalg.norm(P8 @ zc):9.3f} {np.max(np.abs(G1 @ P8 @ ker)):9.1e}", flush=True)
    print("\n  Reading: 'sym' equals 'aniso' (the quintet identity of Part A, now on a true solution); the 8-d and constrained")
    print("  values saturate after the O(1) transient, the gauge mode stays O(1), Gauss is conserved to roundoff: on the")
    print("  self-consistent background the off-diagonal sector is as marginal as the diagonal anisotropic one, on both")
    print("  sides of M = m_phi, below the isotropic threshold.  (On the M -> infinity background the gauge direction grew")
    print("  exponentially at omega = 0.8 -- the artefact that exposed the inconsistency; see axion_dressed.py.)")


if __name__ == "__main__":
    print(f"offdiag_axion.py  python {sys.version.split()[0]}  numpy {np.__version__}  ({sys.executable})")
    Kb, Gb, Zb, q4 = load()
    part_A(Kb, Gb, Zb, q4)
    stamp("Part A done")
    part_B(Kb, Gb, Zb, q4)
    stamp("Part B done")
    print(f"\n({time.time() - T0:.0f} s)")
