#!/usr/bin/env python3
"""
probe_derive.py -- plan item 5 (16 Oct): representative real particles coupled to the Lee-Wick
ghost on the homogeneous SU(2) backgrounds, derived symbolically with PASS/FAIL checks, and the
shared right-hand side probe_core.py checked against the derivation.  Drive-agnostic: the axion's
kappa(t) x y z and the dilaton's Theta(t) (column A of dilaton_su2/code/dilaton_derive.py) are both
carried, so the same probe sector serves axion_su2/ and dilaton_su2/.

Why these probes.  The Yang-Mills paper could couple nothing to the ghost; these backgrounds can:
anything with SU(2) charge couples to the full gauge field q = p - w, ghost included, through its
covariant derivative, with no coupling constant to invent (the charge is the representation, the
strength g, the mass its own).  Boson: a fundamental scalar doublet (the Higgs is one).  Fermion:
a Dirac doublet (a quark or lepton doublet), classical (c-number) spinor, kept in phase space.
Light: a photon is not an SU(2) particle and its homogeneous mode is a free constant field; it
enters this space only through a carrier charged under both groups -- not done here.

Theory (metric eta = diag(+1,-1,-1,-1); F, J as in ym/hdym_derive.py; T^a = sigma^a/2):

    L = -1/4 F.F + 1/(2M^2) J.J  [+ axion, + dilaton as in the two derive scripts]
        + (D_mu Phi)^dagger (D^mu Phi) - ms^2 Phi^dagger Phi                 D_mu Phi = d_mu Phi - i g A^a_mu T^a Phi
        + psibar (i gamma^mu D_mu - mD) psi                                  D_mu psi likewise

Checks (each prints PASS/FAIL):

  CHECK 1  Scalar doublet, all 12 gauge functions and A_0 kept.  (a) With A_0 = 0 the scalar
           Lagrangian depends on A only through sum_{a,i} (A^a_i)^2 (the epsilon term of
           T^a T^b cancels against the symmetric A^a_i A^b_i), so its current is diagonal on the
           diagonal ansatz and the truncation survives for ANY Phi.  (b) The A_0 equation gains the
           colour charge rho^a = i g (Phi^dagger T^a Phidot - h.c.); on the ansatz the gauge part
           vanishes (hdym_derive.py CHECK 2), so the Gauss law is rho^a = 0 for all a, which for a
           doublet means Phidot parallel to Phi with a real coefficient (the overall phase carries
           SU(2) charge on a doublet): Phi = f(t) Phi_0, f real, Phi_0 a constant unit doublet.
           (c) On that reduction L_Phi = fdot^2 - (ms^2 + (g^2/4)|q|^2) f^2 for every Phi_0, and the
           off-diagonal gauge equations vanish.
  CHECK 2  Adjoint triplet (noted, not in the core): the Gauss law forces a fixed direction n, the
           truncation survives only for n along an axis, and the probe is then anisotropic
           (mass g^2 (x^2 + y^2) f^2 for n = z-hat).  Reported for the record.
  CHECK 3  Dirac doublet.  H_D = mD beta - g A^a_i alpha^i T^a is Hermitian.  The diagonal SO(3)
           (spin + colour, the background's symmetry) has a two-dimensional singlet sector of the
           8-dimensional spinor; on the isotropic background H_D maps it to itself (consistent
           truncation), there rho^a = psi^dagger T^a psi = 0 and the current J^a_i = psi^dagger
           alpha^i T^a psi is j delta^a_i (the Gauss law holds and the diagonal truncation of the
           gauge field survives); the reduced Hamiltonian is h = [[mD, c g Q],[c g Q, -mD]] with c
           found here (and written into probe_core.C_FERMION), j = -(c/3) zeta^dagger sigma_x zeta.
           Each alpha^i T^i preserves the singlet separately, so the sector is invariant on EVERY
           diagonal background x, y, z with h as above and Q = (x + y + z)/3 (the fermion sees only
           the trace); a quintet (traceless symmetric) perturbation of A^a_i annihilates the singlet,
           so the fermion does not couple to the anisotropic sector at linear order.  The singlet spinor's current is j delta^a_i
           for ANY singlet state (so it never sources anisotropy).
  CHECK 4  The full (q, w, f, zeta) system from L' (the (q, w) auxiliary form, both drives, sigma =
           +1 ghost / -1 twin): the Euler-Lagrange equations are derived with sympy and
           probe_core.rhs is checked against them at 30 random points for each of the 2 x 2 x 2 x 2
           switch combinations (ghost/twin, axion drive on/off, dilaton on/off, each probe on/off;
           the twin only with Theta = 1), to 1e-10 relative.  probe_core.energy is checked to be
           the Legendre transform of L' symbolically and to be conserved numerically (constant
           drives) over a short integration.
  CHECK 5  Linearisation on x = y = z with the probes present, in the (q, w) variables (M, g kept):
           the anisotropic sector is
               eps''''/M^2 + Theta eps'' + Theta-dot eps' + (kappa Q + (g^2/2) f^2) eps = 0
           (the scalar adds a direct, isotropic effective mass (g^2/2) f^2 to the formerly free
           directions -- a second Hill drive at twice the scalar's frequency; the fermion adds
           nothing there, CHECK 3), and the isotropic tangent and the background's own equations
           with the probe forces are printed.

Usage (from the repository root, seconds):
    python axion_su2/code/probe_derive.py > axion_su2/results/probe_derive.txt
Added 2026-10-16.
"""
import os
import sys
import time

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import probe_core as pc  # noqa: E402

T0 = time.time()


def stamp(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def report(name, ok):
    print(f"    {name}: {'PASS' if ok else 'FAIL'}", flush=True)
    return bool(ok)


t = sp.Symbol('t', real=True)
g, M, ms, mD = sp.symbols('g M m_s m_D', positive=True)
eps = sp.LeviCivita
eta = sp.diag(1, -1, -1, -1)
I = sp.I
sig = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
Tc = [s / 2 for s in sig]
I2 = sp.eye(2)

# ---------------------------------------------------------------- gauge field
A = [[sp.Function(f'A{mu}{c}')(t) for c in range(3)] for mu in range(4)]
x, y, z = [sp.Function(n)(t) for n in 'xyz']
q = [x, y, z]
ans = {}
for c in range(3):
    ans[A[0][c]] = sp.Integer(0)
    for i in range(3):
        ans[A[i + 1][c]] = q[i] if i == c else sp.Integer(0)


def on_ansatz(e):
    return sp.expand(e.subs(ans).doit())


Amat = lambda mu: sum((A[mu][a] * Tc[a] for a in range(3)), sp.zeros(2))   # A_mu = A^a_mu T^a

# ---------------------------------------------------------------- CHECK 1: scalar doublet
stamp("CHECK 1: scalar doublet")
ph = [sp.Function(f'phi{k}', real=True)(t) for k in range(4)]          # Phi = (ph0 + i ph1, ph2 + i ph3)
Phi = sp.Matrix([ph[0] + I * ph[1], ph[2] + I * ph[3]])
def cj(Mx):
    """conjugate for matrices of real functions (sympy leaves conjugate(Derivative(.)) unevaluated)"""
    return Mx.applyfunc(lambda e: e.subs(I, -I)).T


Phid = cj(Phi)


def Dmu(mu, X):
    return (X.diff(t) if mu == 0 else sp.zeros(2, 1)) - I * g * Amat(mu) * X


LPhi = sp.Integer(0)
for mu in range(4):
    DX = Dmu(mu, X=Phi)
    LPhi += eta[mu, mu] * (cj(DX) * DX)[0, 0]
LPhi = sp.expand(LPhi - ms**2 * (Phid * Phi)[0, 0])
LPhi_A0 = sp.expand(LPhi.subs({A[0][c]: 0 for c in range(3)}))
S2 = sum(A[i + 1][a]**2 for i in range(3) for a in range(3))
target = sp.expand(sum(d_.diff(t)**2 for d_ in ph) - (ms**2 + g**2 / 4 * S2) * sum(d_**2 for d_ in ph))
ok1 = report("(a) with A_0 = 0: L_Phi = |Phidot|^2 - (ms^2 + (g^2/4) sum (A^a_i)^2) |Phi|^2", sp.expand(LPhi_A0 - target) == 0)
# (b) colour charge: the A_0 equation of L_Phi
rho = [sp.expand(sp.euler_equations(LPhi, [A[0][c]], t)[0].lhs.subs({A[0][k]: 0 for k in range(3)})) for c in range(3)]
rho_pred = [sp.expand(I * g * ((Phid * Tc[c] * Phi.diff(t))[0, 0] - (cj(Phi.diff(t)) * Tc[c] * Phi)[0, 0]))
            for c in range(3)]
ok1 &= report("(b) the A_0 equation of L_Phi is the colour charge i g (Phi^+ T^a Phidot - h.c.)",
              all(sp.expand(rho[c] - rho_pred[c]) == 0 for c in range(3)))
# rho = 0 for all a  <=>  Phidot = (real) Phi.  Show: with Phi = f(t) Phi_0 (f real, Phi_0 constant) rho = 0 for
# every Phi_0, and with a phase e^{i alpha(t)} on top rho^a = -g alphadot Phi_0^+ sigma^a Phi_0 f^2 != 0.
f_ = sp.Function('f', real=True)(t)
al = sp.Function('alpha', real=True)(t)
c0 = sp.symbols('c0:4', real=True)
Phi0 = sp.Matrix([c0[0] + I * c0[1], c0[2] + I * c0[3]])
sub_f = {ph[0]: f_ * c0[0], ph[1]: f_ * c0[1], ph[2]: f_ * c0[2], ph[3]: f_ * c0[3]}
ok1 &= report("(b) Phi = f(t) Phi_0, f real: rho^a = 0 for every constant Phi_0",
              all(sp.expand(rho[c].subs(sub_f).doit()) == 0 for c in range(3)))
Phi_al = sp.exp(I * al) * f_ * Phi0
sub_al = {ph[0]: sp.re(Phi_al[0]), ph[1]: sp.im(Phi_al[0]), ph[2]: sp.re(Phi_al[1]), ph[3]: sp.im(Phi_al[1])}
rho_al = [sp.simplify(rho[c].subs(sub_al).doit()) for c in range(3)]
pred_al = [sp.simplify(-2 * g * al.diff(t) * f_**2 * (cj(Phi0) * Tc[c] * Phi0)[0, 0]) for c in range(3)]
ok1 &= report("(b) an overall phase alpha(t) carries SU(2) charge: rho^a = -g alphadot f^2 Phi_0^+ sigma^a Phi_0",
              all(sp.simplify(rho_al[c] - pred_al[c]) == 0 for c in range(3)))
# (c) on the ansatz with Phi = f Phi_0, |Phi_0| = 1
Lred = sp.expand(on_ansatz(LPhi).subs(sub_f).doit())
norm = sum(ci**2 for ci in c0)
Lred_target = sp.expand(norm * (f_.diff(t)**2 - (ms**2 + g**2 / 4 * (x**2 + y**2 + z**2)) * f_**2))
ok1 &= report("(c) on the ansatz: L_Phi = |Phi_0|^2 [fdot^2 - (ms^2 + (g^2/4)(x^2 + y^2 + z^2)) f^2]",
              sp.expand(Lred - Lred_target) == 0)
funcs_i = [A[mu][c] for mu in range(1, 4) for c in range(3)]
EOM_Phi = sp.euler_equations(LPhi_A0, funcs_i, t)
ok1 &= report("(c) the off-diagonal gauge equations of L_Phi vanish on the ansatz, for any Phi",
              all(on_ansatz(EOM_Phi[3 * i + c].lhs) == 0 for i in range(3) for c in range(3) if i != c))
diag_force = [on_ansatz(EOM_Phi[3 * i + i].lhs).subs(sub_f).doit() for i in range(3)]
ok1 &= report("(c) its force on q_i is -(g^2/2) q_i |Phi_0|^2 f^2",
              all(sp.expand(diag_force[i] + g**2 / 2 * q[i] * norm * f_**2) == 0 for i in range(3)))

# ---------------------------------------------------------------- CHECK 2: adjoint triplet (record only)
stamp("CHECK 2: adjoint triplet (for the record)")
ch = [sp.Function(f'chi{a}')(t) for a in range(3)]
Lchi = sp.Integer(0)
for mu in range(4):
    Dchi = [(ch[a].diff(t) if mu == 0 else 0) + g * sum(eps(a, b, c) * A[mu][b] * ch[c] for b in range(3) for c in range(3))
            for a in range(3)]
    Lchi += eta[mu, mu] * sum(dc**2 for dc in Dchi) / 2
Lchi = sp.expand(Lchi - ms**2 / 2 * sum(c_**2 for c_ in ch))
rho3 = [sp.expand(sp.euler_equations(Lchi, [A[0][c]], t)[0].lhs.subs({A[0][k]: 0 for k in range(3)})) for c in range(3)]
nvec = sp.symbols('n0:3', real=True)
sub_n = {ch[a]: f_ * nvec[a] for a in range(3)}
ok2 = report("triplet: chi = f(t) n (fixed direction) has zero colour charge",
             all(sp.expand(rho3[c].subs(sub_n).doit()) == 0 for c in range(3)))
EOM_chi = sp.euler_equations(sp.expand(Lchi.subs({A[0][c]: 0 for c in range(3)})), funcs_i, t)
offd = [sp.factor(on_ansatz(EOM_chi[3 * i + c].lhs).subs(sub_n).doit()) for i in range(3) for c in range(3) if i != c]
print("    triplet: off-diagonal gauge equations on the ansatz:", offd)
ok2 &= report("triplet: they vanish for n along an axis (n = z-hat) and not otherwise",
              all(sp.expand(e.subs({nvec[0]: 0, nvec[1]: 0, nvec[2]: 1})) == 0 for e in offd)
              and any(sp.expand(e.subs({nvec[0]: 1, nvec[1]: 1, nvec[2]: 0})) != 0 for e in offd))
Lchi_red = sp.expand(on_ansatz(Lchi).subs(sub_n).subs({nvec[0]: 0, nvec[1]: 0, nvec[2]: 1}).doit())
ok2 &= report("triplet, n = z-hat: L = fdot^2/2 - (ms^2/2 + (g^2/2)(x^2 + y^2)) f^2  (anisotropic probe)",
              sp.expand(Lchi_red - (f_.diff(t)**2 / 2 - (ms**2 / 2 + g**2 / 2 * (x**2 + y**2)) * f_**2)) == 0)
print("    not in probe_core.py: it breaks the isotropy of the background by construction.")

# ---------------------------------------------------------------- CHECK 3: Dirac doublet
stamp("CHECK 3: Dirac doublet, the diagonal-SO(3) singlet sector")
Z2 = sp.zeros(2)
beta = sp.diag(1, 1, -1, -1)
alpha = [sp.Matrix(sp.BlockMatrix([[Z2, s], [s, Z2]])) for s in sig]
Sspin = [sp.Matrix(sp.BlockMatrix([[s / 2, Z2], [Z2, s / 2]])) for s in sig]


def kron(a, b):
    return sp.Matrix(sp.kronecker_product(a, b))


I4 = sp.eye(4)
HD_of = lambda Aij: mD * kron(beta, I2) - g * sum((Aij[i][a] * kron(alpha[i], Tc[a]) for i in range(3) for a in range(3)), sp.zeros(8))
Qs = sp.Symbol('Q', real=True)
Aiso = [[Qs if i == a else 0 for a in range(3)] for i in range(3)]
xs, ys, zs = sp.symbols('x y z', real=True)
Adiag = [[[xs, ys, zs][i] if i == a else 0 for a in range(3)] for i in range(3)]
Agen = [[sp.Symbol(f'a{i}{a}', real=True) for a in range(3)] for i in range(3)]
Hgen = HD_of(Agen)
ok3 = report("H_D = mD beta - g A^a_i alpha^i T^a is Hermitian (12 gauge components free)",
             sp.simplify(Hgen - Hgen.conjugate().T) == sp.zeros(8))
Jtot = [kron(Sspin[i], I2) + kron(I4, Tc[i]) for i in range(3)]
J2 = sum((Ji * Ji for Ji in Jtot), sp.zeros(8))
sing = J2.nullspace()
ok3 &= report("the diagonal SO(3) singlet sector of the 8-spinor is two-dimensional", len(sing) == 2)
# orthonormal basis, ordered (upper block, lower block)
Sb = sp.Matrix.hstack(*[v / sp.sqrt((v.conjugate().T * v)[0, 0]) for v in sing])
if Sb[0, 0] == 0 and Sb[1, 0] == 0 and Sb[2, 0] == 0 and Sb[3, 0] == 0:
    Sb = sp.Matrix.hstack(Sb[:, 1], Sb[:, 0])
print("    singlet basis (columns):", [list(Sb[:, k]) for k in range(2)])
Proj = Sb * Sb.conjugate().T
Hiso = HD_of(Aiso)
leak_iso = sp.simplify((sp.eye(8) - Proj) * Hiso * Sb)
ok3 &= report("isotropic background: H_D maps the singlet sector to itself", leak_iso == sp.zeros(8, 2))
h = sp.simplify(Sb.conjugate().T * Hiso * Sb)
print(f"    reduced h(Q) = {h.tolist()}")
cval = sp.simplify(h[0, 1] / (g * Qs))
ok3 &= report(f"h = [[mD, c g Q],[c g Q, -mD]] with c = {cval}",
              sp.simplify(h - sp.Matrix([[mD, cval * g * Qs], [cval * g * Qs, -mD]])) == sp.zeros(2) and cval.is_number)
ok3 &= report(f"probe_core.C_FERMION = {pc.C_FERMION} agrees", abs(float(cval) - pc.C_FERMION) < 1e-15)
zr = sp.symbols('z0r z0i z1r z1i', real=True)
zeta = sp.Matrix([zr[0] + I * zr[1], zr[2] + I * zr[3]])
psi = Sb * zeta
rho_f = [sp.simplify((psi.conjugate().T * kron(I4, Tc[a]) * psi)[0, 0]) for a in range(3)]
ok3 &= report("any singlet state has rho^a = psi^+ T^a psi = 0", all(r == 0 for r in rho_f))
Jcur = [[sp.simplify((psi.conjugate().T * kron(alpha[i], Tc[a]) * psi)[0, 0]) for a in range(3)] for i in range(3)]
jpred = sp.simplify(-cval / 3 * (zeta.conjugate().T * sig[0] * zeta)[0, 0])
ok3 &= report("any singlet state has J^a_i = j delta^a_i with j = -(c/3) zeta^+ sigma_x zeta",
              all(sp.simplify(Jcur[i][a] - (jpred if i == a else 0)) == 0 for i in range(3) for a in range(3)))
Hdiag = HD_of(Adiag)
leak = sp.simplify((sp.eye(8) - Proj) * Hdiag * Sb)
ok3 &= report("every diagonal background x, y, z: H_D maps the singlet to itself (each alpha^i T^i does separately)",
              leak == sp.zeros(8, 2))
hdiag = sp.simplify(Sb.conjugate().T * Hdiag * Sb)
ok3 &= report("there h = [[mD, c g Q],[c g Q, -mD]] with Q = (x + y + z)/3, the same c",
              sp.simplify(hdiag - sp.Matrix([[mD, cval * g * (xs + ys + zs) / 3], [cval * g * (xs + ys + zs) / 3, -mD]])) == sp.zeros(2))
# quintet perturbation: traceless symmetric delta A^a_i annihilates the singlet
dA = sp.Matrix(3, 3, lambda i, a: sp.Symbol(f'd{i}{a}', real=True))
dA = (dA + dA.T) / 2
dA = dA - sp.eye(3) * dA.trace() / 3
dH = -g * sum((dA[i, a] * kron(alpha[i], Tc[a]) for i in range(3) for a in range(3)), sp.zeros(8))
ok3 &= report("a traceless symmetric (quintet) perturbation of A^a_i annihilates the singlet: no linear coupling",
              sp.simplify(dH * Sb) == sp.zeros(8, 2))

# ---------------------------------------------------------------- CHECK 4: the full system and probe_core
stamp("CHECK 4: the full (q, w, f, zeta) system against probe_core.rhs")
wv = [sp.Function(f'w{n}')(t) for n in 'xyz']
Th = sp.Function('Theta', positive=True)(t)
kt = sp.Function('kappa_t')(t)
V = g**2 / 2 * (x**2 * y**2 + y**2 * z**2 + z**2 * x**2)
gradV = [V.diff(v) for v in q]
zf = [sp.Function(f'zeta{k}', real=True)(t) for k in range(4)]
zc = sp.Matrix([zf[0] + I * zf[1], zf[2] + I * zf[3]])
Qf = (x + y + z) / 3
hQ = sp.Matrix([[mD, cval * g * Qf], [cval * g * Qf, -mD]])
# L_D on the singlet: i zeta^+ zetadot - zeta^+ h zeta  (its EL equation is i zetadot = h zeta); the gauge force is
# dL_D/dq_i = g j  (CHECK 3), which we add to L' through the term g j Q * 3 ... i.e. L_D depends on q only via h(Q).
# L_D = Re[i zeta^+ zetadot - zeta^+ h zeta] in real components (its EL equation is i zetadot = h zeta)
LD = sp.expand(zf[1] * zf[0].diff(t) - zf[0] * zf[1].diff(t) + zf[3] * zf[2].diff(t) - zf[2] * zf[3].diff(t)
               - mD * (zf[0]**2 + zf[1]**2 - zf[2]**2 - zf[3]**2) - 2 * cval * g * Qf * (zf[0] * zf[2] + zf[1] * zf[3]))
results = []
for sigma in (+1, -1):
    Lp = (Th * (sp.Rational(1, 2) * sum(v.diff(t)**2 for v in q) - V)
          + sum(wv[i].diff(t) * q[i].diff(t) for i in range(3))
          - sum(wv[i] * gradV[i] for i in range(3)) + sigma * M**2 / 2 * sum(wi**2 for wi in wv)
          + kt * x * y * z - g**2 / 4 * (x**2 + y**2 + z**2) * f_**2 + f_.diff(t)**2 - ms**2 * f_**2 + LD)
    if sigma == -1:
        Lp = Lp.subs(Th, 1) + sum(wi.diff(t)**2 for wi in wv)
    EL = [sp.expand(e.lhs) for e in sp.euler_equations(Lp, q + wv + [f_], t)]
    acc = [q[i].diff(t, 2) for i in range(3)] + [wv[i].diff(t, 2) for i in range(3)] + [f_.diff(t, 2)]
    sol = sp.solve(EL, acc, dict=True)[0]
    results.append((sigma, Lp, sol))
    # energy (Legendre transform): sum vdot dL/dvdot - L, in terms of the velocities
    vel = [q[i].diff(t) for i in range(3)] + [wv[i].diff(t) for i in range(3)] + [f_.diff(t)] + [zz.diff(t) for zz in zf]
    E_sym = sp.expand(sum(vv * Lp.diff(vv) for vv in vel) - Lp)
    results[-1] = (sigma, Lp, sol, E_sym)
print("    Euler-Lagrange equations solved for the accelerations in both sigma", flush=True)

rng = np.random.default_rng(16)
ok4 = True
Thn = lambda t_: 1.0 + 0.3 * np.cos(1.7 * t_)
Thdn = lambda t_: -0.3 * 1.7 * np.sin(1.7 * t_)
kan = lambda t_: 0.4 * np.sin(2.3 * t_)
syms_q = [sp.Symbol(f'q{i}') for i in range(3)]
syms_v = [sp.Symbol(f'v{i}') for i in range(3)]
syms_w = [sp.Symbol(f'w{i}') for i in range(3)]
syms_u = [sp.Symbol(f'u{i}') for i in range(3)]
sf, sfd = sp.symbols('f fd')
sz = sp.symbols('z0r z0i z1r z1i')
sTh, sThd, ska = sp.symbols('Th Thd ka')
plain = {}
for i in range(3):
    plain[q[i].diff(t)] = syms_v[i]; plain[q[i]] = syms_q[i]
    plain[wv[i].diff(t)] = syms_u[i]; plain[wv[i]] = syms_w[i]
plain[f_.diff(t)] = sfd; plain[f_] = sf
for k in range(4):
    plain[zf[k]] = sz[k]
plain[Th.diff(t)] = sThd; plain[Th] = sTh; plain[kt] = ska
argl = syms_q + syms_v + syms_w + syms_u + [sf, sfd] + list(sz) + [sTh, sThd, ska, g, M, ms, mD]


def to_num(expr):
    e = expr.subs({zz.diff(t): 0 for zz in zf})   # zeta velocities appear only in L_D's first-order term
    # derivatives before the functions they are derivatives of (a dict subs would zero Theta-dot)
    for k_, v_ in sorted(plain.items(), key=lambda kv: -int(isinstance(kv[0], sp.Derivative))):
        e = e.subs(k_, v_)
    return sp.lambdify(argl, e, 'numpy')


maxerr = 0.0
for sigma, Lp, sol, E_sym in results:
    fq = [to_num(sol[q[i].diff(t, 2)]) for i in range(3)]
    fw = [to_num(sol[wv[i].diff(t, 2)]) for i in range(3)]
    ff = to_num(sol[f_.diff(t, 2)])
    fE = to_num(E_sym)
    for ax in (False, True):
        for dil in ((False,) if sigma == -1 else (False, True)):
            for bos in (False, True):
                for fer in (False, True):
                    P = pc.Params(g=0.9, M=1.3, ms=0.7, mD=0.6, sigma=sigma, boson=bos, fermion=fer,
                                  kappa=kan if ax else None, Theta=Thn if dil else None, Thetadot=Thdn if dil else None)
                    err_c = 0.0
                    for _ in range(30):
                        s = rng.normal(size=18)
                        if not bos:
                            s[12:14] = 0
                        if not fer:
                            s[14:18] = 0
                        tt = rng.uniform(0, 10)
                        vals = list(s[0:12]) + [s[12], s[13]] + [s[14], s[16], s[15], s[17]] + [P.Theta(tt), P.Thetadot(tt), P.kappa(tt), 0.9, 1.3, 0.7, 0.6]   # zeta = (z0r, z0i, z1r, z1i)
                        num = pc.rhs(tt, s, P)
                        ref_q = np.array([fq[i](*vals) for i in range(3)])
                        ref_w = np.array([fw[i](*vals) for i in range(3)])
                        ref_f = ff(*vals) if bos else 0.0
                        zeta_n = s[14:16] + 1j * s[16:18]
                        cq = float(cval) * 0.9 * (s[0] + s[1] + s[2]) / 3
                        hn = np.array([[0.6, cq], [cq, -0.6]])
                        ref_z = -1j * (hn @ zeta_n) if fer else np.zeros(2)
                        ref = np.concatenate([s[3:6], ref_q, s[9:12], ref_w, [s[13] if bos else 0, ref_f], ref_z.real, ref_z.imag])
                        scale = max(1.0, np.max(np.abs(ref)))
                        err_c = max(err_c, np.max(np.abs(num - ref)) / scale)
                        Eref = fE(*vals)
                        Enum = pc.energy(tt, s, P)
                        err_c = max(err_c, abs(Enum - Eref) / max(1.0, abs(Eref)))
                    print(f"      sigma={sigma:+d} axion={int(ax)} dilaton={int(dil)} boson={int(bos)} fermion={int(fer)}: max rel err {err_c:.1e}")
                    maxerr = max(maxerr, err_c)
ok4 &= report(f"probe_core.rhs and .energy agree with the Euler-Lagrange equations and the Legendre transform "
              f"at random points, all switch combinations (max relative error {maxerr:.1e})", maxerr < 1e-10)
# energy conservation, constant drives, both sigma, both probes on
from scipy.integrate import solve_ivp  # noqa: E402
drift = 0.0
for sigma in (+1, -1):
    P = pc.Params(sigma=sigma, boson=True, fermion=True, ms=0.8, mD=0.5)
    s0 = np.zeros(18)
    s0[0:3] = [0.5, 0.45, 0.55]; s0[3:6] = [0.1, -0.05, 0.02]; s0[6:9] = 0.01 * rng.normal(size=3)
    s0[12], s0[13] = 0.3, 0.0; s0[14:18] = [0.8, 0.1, 0.3, -0.2]
    sol_ = solve_ivp(lambda tt, ss: pc.rhs(tt, ss, P), (0, 20.0), s0, method='DOP853', rtol=1e-11, atol=1e-13, dense_output=False,
                     t_eval=np.linspace(0, 20, 21))
    E = np.array([pc.energy(tt, ss, P) for tt, ss in zip(sol_.t, sol_.y.T)])
    nrm = np.array([np.sum(ss[14:18]**2) for ss in sol_.y.T])
    drift = max(drift, np.max(np.abs(E - E[0])) / max(1.0, abs(E[0])), np.max(np.abs(nrm - nrm[0])))
ok4 &= report(f"energy and the spinor norm are conserved over t = 20 with both probes on, ghost and twin (drift {drift:.1e})", drift < 1e-8)

# ---------------------------------------------------------------- CHECK 5: linearisation on x = y = z with the probes
stamp("CHECK 5: linearisation on x = y = z with the probes (sigma = +1)")
sigma, Lp, sol, E_sym = results[0]
EL = [sp.expand(e.lhs) for e in sp.euler_equations(Lp, q + wv, t)]
Qb = sp.Function('Q')(t)
w0 = sp.Function('w0')(t)
eps_ = sp.Function('epsilon')(t)
zet = sp.Function('zeta_w')(t)
lam = sp.Symbol('lambda_')


def project(vec):
    sub = {}
    for i in range(3):
        sub[q[i]] = Qb + lam * eps_ * vec[i]
        sub[wv[i]] = w0 + lam * zet * vec[i]
    lin = [sp.expand(e.subs(sub).doit()).coeff(lam, 1) for e in EL]
    qeq = weq = None
    for i in range(3):
        if vec[i] != 0:
            cq = sp.expand(lin[i] / vec[i]); cw = sp.expand(lin[3 + i] / vec[i])
            if qeq is None:
                qeq, weq = cq, cw
            else:
                assert sp.expand(cq - qeq) == 0 and sp.expand(cw - weq) == 0
        else:
            assert sp.expand(lin[i]) == 0 and sp.expand(lin[3 + i]) == 0
    return qeq, weq


# the fermion and scalar are held at their background values (f(t), zeta(t) given): the fermion does not respond
# to the quintet (CHECK 3); the scalar responds to |q|^2, which is second order along an anisotropic direction.
an_q, an_w = project([1, -1, 0])
zsol = sp.solve(an_w, zet)[0]
an4 = sp.expand(an_q.subs(zet, zsol).doit())
ref = eps_.diff(t, 4) / M**2 + (Th * eps_.diff(t)).diff(t) + (kt * Qb + g**2 / 2 * f_**2) * eps_
ok5 = report("anisotropic sector:  eps''''/M^2 + (Theta eps')' + (kappa Q + (g^2/2) f^2) eps = 0  (Q, w0, zeta drop out)",
             sp.expand(an4 + ref) == 0)
sub_bg = {}
for i in range(3):
    sub_bg[q[i]] = Qb; sub_bg[wv[i]] = w0
bg = [sp.expand(e.subs(sub_bg).doit()) for e in EL]
wsol = sp.solve(bg[3], w0)[0]
iso4 = sp.expand(bg[0].subs(w0, wsol).doit())
print(f"    background (x = y = z) 4th-order equation with the probes (= 0):\n      {sp.collect(iso4, [Th, Th.diff(t), kt, f_, zf[0]])}")
lim = sp.expand(sp.limit(iso4, M, sp.oo))
jsym = sp.expand(-cval / 3 * 2 * (zf[0] * zf[2] + zf[1] * zf[3]))
ok5 &= report("background, M -> infinity: Theta (Q'' + 2 g^2 Q^3) + Theta-dot Q' = kappa Q^2 - (g^2/2) f^2 Q + g j",
              sp.expand(lim + Th * (Qb.diff(t, 2) + 2 * g**2 * Qb**3) + Th.diff(t) * Qb.diff(t) - kt * Qb**2 + g**2 / 2 * f_**2 * Qb - g * jsym) == 0)
is_q, is_w = project([1, 1, 1])
zsol_i = sp.solve(is_w, zet)[0]
is4 = sp.expand(is_q.subs(zet, zsol_i).doit())
print(f"    isotropic tangent (f, zeta frozen; = 0):\n      {sp.collect(is4, [eps_.diff(t, 4), eps_.diff(t, 2), eps_.diff(t), eps_])}")
print("    (with the probes responding, the isotropic sector is the coupled (Q, w, f, zeta) tangent: run, not printed.)")

# ---------------------------------------------------------------- summary
allok = all([ok1, ok2, ok3, ok4, ok5])
stamp(f"ALL CHECKS {'PASSED' if allok else 'NOT ALL PASSED'}")
print("""
Summary (units g = M = 1):
  boson   one real field f (the modulus of a scalar doublet; its phase is fixed by the Gauss law), effective
          mass^2 = ms^2 + (g^2/4)|q|^2, force -(g^2/2) q f^2 on the gauge field (ghost included through q = p - w);
          it gives the anisotropic directions the effective mass (g^2/2) f^2 and the off-diagonal ones the same.
  fermion a two-level classical Dirac doublet in the diagonal-SO(3) singlet, i zeta' = h(Q) zeta,
          h = [[mD, cgQ],[cgQ, -mD]], c = 3/2, force g j with j = -(c/3) zeta^+ sigma_x zeta on each gauge direction;
          exact on every diagonal background (it sees only x + y + z), blind to the anisotropic sector at linear order.
  twin    sigma = -1: a healthy heavy sector in the (p, w) variables, Theta = 1 only; a control, not a model.
  probe_core.py carries all of it; this script is its test.
""")
