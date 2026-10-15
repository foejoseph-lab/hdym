#!/usr/bin/env python3
"""
axion_derive.py -- stage 3, step 1: what the axion coupling phi F Ftilde adds to the homogeneous
SU(2) Lee-Wick system, derived symbolically with PASS/FAIL checks (no numbers from the literature
are trusted; the one convention taken from Fujita, Murai, Nakayama, Yin, arXiv:2312.06889v2, is
checked against their eq. (13), see CHECK 2).

Theory (metric eta = diag(+1,-1,-1,-1), same construction as ym/hdym_derive.py):

    L = -1/4 F^a_{mu nu} F^{a mu nu} + 1/(2 M^2) (D^mu F_{mu nu})^a (D_rho F^{rho nu})^a
        + (1/(4 f)) phi(t) F^a_{mu nu} Ftilde^{a mu nu},
    Ftilde^{a mu nu} = eps^{mu nu rho sigma} F^a_{rho sigma} / (2 sqrt(-det g)),

phi(t) the homogeneous axion, treated as a given function of time (its own equation, their (12),
is not used here).  F is the same object as in ym/hdym_derive.py (F = dA - dA + g eps A A), so the
Lee-Wick term, built from F alone, is untouched by the axion; only the new term is examined.

Checks (each prints PASS/FAIL):

  CHECK 1  Homogeneous reduction of the axion term (all 12 functions A^a_mu(t), A_0 kept):
             (1/4f) phi F.Ftilde  =  s (g/f) phi d/dt det(A^a_i),   s = +1 or -1
           (which sign is reported: it is the choice of eps^{0123}), i.e. the term is
             L_axion = s (g/f) phidot det A   up to a total time derivative,
           and it does not contain A_0 at all (so the Gauss law is unchanged).
  CHECK 2  The sign s that reproduces the isotropic equation of Fujita et al. (13) at H = 0,
             Qddot + 2 g^2 Q^3 = (g/f) Q^2 phidot,
           is recorded; with it L_axion = (g/f) phidot det A.  (This fixes the convention; the
           physics after inflation, where phidot oscillates, does not depend on it.)
  CHECK 3  The diagonal ansatz A^a_i = delta^a_i q_i, A_0 = 0 remains a consistent truncation:
           the off-diagonal Euler-Lagrange equations of L_axion vanish on it, and on it
           L_axion = (g/f) phidot x y z.
  CHECK 4  Auxiliary-field form: adding (g/f) phidot x y z to L' of hdym_derive.py CHECK 4
           reproduces the 4th-order equations of L_red + (g/f) phidot x y z.  In the reduced
           system of ym/ the axion is therefore the time-dependent potential  -kappa(t) x y z,
           kappa = (g/f) phidot, in the field q = p - w, and nothing else.
  CHECK 5  Scaling t -> tau/M, q -> (M/g) u (hdym_derive.py CHECK 5):
             kappa x y z = (M^4/g^2) (kappa/(g M)) u_x u_y u_z,
           so in units g = M = 1 the one new quantity is kappa/(g M) = phidot/(f M) = 2 xi H/M.
  CHECK 6  Linearisation on the isotropic orbit x = y = z in the (q, w) variables of L' (CHECK 4),
           M and g kept.  The drive sources the ghost: on shell w = (qddot + grad V)/M^2 and the
           driven background has qddot + grad V = kappa grad(xyz), so w = W = 0 is not invariant
           and the isotropic background is (Q, Q, Q; w0, w0, w0).  Checks: the isotropic
           truncation is consistent; its 4th-order equation tends to Fujita (13) at H = 0 as
           M -> infinity; the two anisotropic directions (the permutation-odd ones, degenerate)
           obey, with w0 dropping out exactly,
               eps''''/M^2 + eps'' + kappa(t) Q(t) eps = 0,
           a Pais-Uhlenbeck pair with frozen frequencies omega^2 = (M^2/2)(1 +- sqrt(1 - 4 kappa Q/M^2)):
           Yang-Mills branch ~ kappa Q (tachyonic for kappa Q < 0), ghost branch ~ M^2 - kappa Q,
           complex pair for kappa Q > M^2/4.  Without the axion this sector is exactly free
           (axion_su2/README.md, stage 1).  The isotropic tangent equation is printed; as
           M -> infinity it is eps'' + (6 g^2 Q^2 - 2 kappa Q) eps = 0.
  CHECK 7  FRW, ds^2 = a(eta)^2 (d eta^2 - dx^2): sqrt(-g) (1/4f) phi F.Ftilde is independent of
           a(eta) (the Chern-Simons density is conformal), so in conformal time the term is
           (g/f) (d phi/d eta) det A with comoving fields, and stage 2 (M -> M a(eta)) carries over.

Usage (from the repository root, seconds):
    python axion_su2/code/axion_derive.py > axion_su2/results/axion_derive.txt
Added 2026-10-14.
"""
import time

import sympy as sp

T0 = time.time()


def stamp(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def report(name, ok):
    print(f"    {name}: {'PASS' if ok else 'FAIL'}", flush=True)
    return bool(ok)


t = sp.Symbol('t', real=True)
g, M, f = sp.symbols('g M f', positive=True)
eps = sp.LeviCivita
eta = sp.diag(1, -1, -1, -1)
phi = sp.Function('phi')(t)

# ---------------------------------------------------------------- the homogeneous field strength
A = [[sp.Function(f'A{mu}{c}')(t) for c in range(3)] for mu in range(4)]


def d(mu, expr):
    return expr.diff(t) if mu == 0 else sp.Integer(0)


F = [[[None] * 3 for _ in range(4)] for _ in range(4)]
for mu in range(4):
    for nu in range(4):
        for c in range(3):
            F[mu][nu][c] = (d(mu, A[nu][c]) - d(nu, A[mu][c])
                            + g * sum(eps(c, b, e) * A[mu][b] * A[nu][e]
                                      for b in range(3) for e in range(3)))


def FFdual(metric_scale):
    """sqrt(-g) * F^a_{mu nu} Ftilde^{a mu nu} with g_{mu nu} = metric_scale^2 * eta, eps^{0123} = +1
    (the Levi-Civita symbol).  Ftilde^{mu nu} = eps^{mu nu rho sigma} F_{rho sigma} / (2 sqrt(-g)), so the
    contraction F_{mu nu} Ftilde^{mu nu} needs no inverse metric; sqrt(-g) = metric_scale^4 is kept
    explicitly so that CHECK 7 sees it cancel rather than assumes it."""
    s = metric_scale
    sqrtg = s**4
    tot = sp.Integer(0)
    for mu in range(4):
        for nu in range(4):
            for rho in range(4):
                for sig in range(4):
                    e = eps(mu, nu, rho, sig)
                    if e == 0:
                        continue
                    for c in range(3):
                        tot += F[mu][nu][c] * e * F[rho][sig][c] / (2 * sqrtg)
    return sp.expand(sqrtg * tot)


Amat = sp.Matrix(3, 3, lambda c, i: A[i + 1][c])     # A^a_i, row = colour a, column = space i
detA = Amat.det()

# ---------------------------------------------------------------- CHECK 1
stamp("CHECK 1: homogeneous reduction of the axion term")
Lax_full = sp.expand(phi / (4 * f) * FFdual(sp.Integer(1)))
target = sp.expand(g / f * phi * detA.diff(t))
sgn = None
for s_try in (+1, -1):
    if sp.expand(Lax_full - s_try * target) == 0:
        sgn = s_try
ok1 = report(f"(1/4f) phi F.Ftilde = s (g/f) phi d(det A)/dt with s = {sgn} (eps^0123 = +1)", sgn is not None)
has_A0 = any(Lax_full.has(A[0][c]) for c in range(3))
ok1 &= report("the term contains no A_0 (Gauss law unchanged)", not has_A0)
# the Euler-Lagrange equations w.r.t. A_0 vanish identically (same statement, from the equations)
A0eqs = sp.euler_equations(Lax_full, [A[0][c] for c in range(3)], t)
ok1 &= report("A_0 Euler-Lagrange equations of the term vanish", all(sp.expand(e.lhs) == 0 for e in A0eqs))

# ---------------------------------------------------------------- CHECK 2
stamp("CHECK 2: the sign against Fujita et al. (13) at H = 0")
Q = sp.Function('Q')(t)
if sgn is None:
    raise SystemExit("CHECK 1 failed; the rest depends on it")
# L_axion as a Lagrangian: (1/4f) phi F.Ftilde = s (g/f) phi d(det A)/dt = -s (g/f) phidot det A + total derivative
Lax_lag = sp.expand(-sgn * g / f * phi.diff(t) * detA)
iso = {}
for c in range(3):
    iso[A[0][c]] = sp.Integer(0)
    for i in range(3):
        iso[A[i + 1][c]] = Q if i == c else sp.Integer(0)
LYM_iso = sp.expand(sum(-sp.Rational(1, 4) * eta[mu, mu] * eta[nu, nu] * F[mu][nu][c]**2
                        for mu in range(4) for nu in range(4) for c in range(3)).subs(iso).doit())
L_iso = sp.expand(LYM_iso + Lax_lag.subs(iso).doit())
eq_iso = sp.expand(sp.euler_equations(L_iso, [Q], t)[0].lhs)
# Fujita (13) at H = 0, multiplied by 3 (three equal components) and with their sign of the EL expression
ref13 = 3 * (Q.diff(t, 2) + 2 * g**2 * Q**3 - g / f * Q**2 * phi.diff(t))
ratio = sp.simplify(eq_iso / ref13)
# sympy's euler_equations return dL/dQ - d/dt dL/dQdot = -(Qddot + ...), so a match is ratio = -1;
# the opposite sign of the axion term is not proportional to (13) at all.
match = (ratio == -1)
print(f"    isotropic EL equation from eps^0123 = +1 equals -(3 x Fujita (13) at H = 0): {match}")
eps_conv = +1 if match else -1
print(f"    Fujita's (13) corresponds to L_axion = +(g/f) phidot det A, i.e. to eps^0123 = {eps_conv:+d}"
      f" in this construction (with eps^0123 = +1 the term came out as {-sgn:+d} (g/f) phidot det A).")
kappa = sp.Symbol('kappa', real=True)                # kappa(t) = (g/f) phidot(t)
Lax = sp.expand(g / f * phi.diff(t) * detA)          # the convention that reproduces (13), used from here on
eq_check = sp.expand(sp.euler_equations(sp.expand(LYM_iso + Lax.subs(iso).doit()), [Q], t)[0].lhs)
ok2 = report("with L_axion = (g/f) phidot det A the isotropic equation is Fujita (13) at H = 0 (x3)",
             sp.simplify(eq_check / ref13) == -1)

# ---------------------------------------------------------------- CHECK 3
stamp("CHECK 3: diagonal ansatz stays a consistent truncation")
x, y, z = [sp.Function(n)(t) for n in 'xyz']
q = [x, y, z]
ans = {}
for c in range(3):
    ans[A[0][c]] = sp.Integer(0)
    for i in range(3):
        ans[A[i + 1][c]] = q[i] if i == c else sp.Integer(0)
funcs = [A[mu][c] for mu in range(1, 4) for c in range(3)]     # A_0 is absent from Lax (CHECK 1)
ax_eqs = sp.euler_equations(Lax, funcs, t)
EOMax = {(mu, c): ax_eqs[3 * (mu - 1) + c].lhs for mu in range(1, 4) for c in range(3)}
ok3 = True
for i in range(3):
    for c in range(3):
        if i != c:
            ok3 &= report(f"off-diagonal axion EOM (i={i + 1},a={c + 1}) vanishes on the ansatz",
                          sp.expand(EOMax[(i + 1, c)].subs(ans).doit()) == 0)
Lax_red = sp.expand(Lax.subs(ans).doit())
ok3 &= report("L_axion on the ansatz = (g/f) phidot x y z", sp.expand(Lax_red - g / f * phi.diff(t) * x * y * z) == 0)

# ---------------------------------------------------------------- CHECK 4
stamp("CHECK 4: auxiliary-field form with the axion")
V = g**2 / 2 * (x**2 * y**2 + y**2 * z**2 + z**2 * x**2)
Vtot = V - kappa * x * y * z                          # kappa stands for (g/f) phidot, a given function of t
kt = sp.Function('kappa_t')(t)
Vtot_t = V - kt * x * y * z
gradVt = [Vtot_t.diff(v) for v in q]
# 4th-order closed form (hdym_derive.py CHECK 3) with the axion term in V only -- the Lee-Wick term keeps grad V
gradV = [V.diff(v) for v in q]
Lclosed = (sp.Rational(1, 2) * sum(v.diff(t)**2 for v in q) - Vtot_t
           - sum((q[i].diff(t, 2) + gradV[i])**2 for i in range(3)) / (2 * M**2))
red_eqs = sp.euler_equations(Lclosed, q, t)
wv = [sp.Function(f'w{n}')(t) for n in 'xyz']
Lp = (sp.Rational(1, 2) * sum(v.diff(t)**2 for v in q)
      + sum(wv[i].diff(t) * q[i].diff(t) for i in range(3))
      - Vtot_t - sum(wv[i] * gradV[i] for i in range(3))
      + M**2 / 2 * sum(wv[i]**2 for i in range(3)))
aux_eqs = sp.euler_equations(Lp, q + wv, t)
wsol = sp.solve([e.lhs for e in aux_eqs[3:]], wv, dict=True)[0]
ok4 = True
for i in range(3):
    lhs = sp.expand(aux_eqs[i].lhs.subs(wsol).doit())
    ref = sp.expand(red_eqs[i].lhs)
    r = sp.simplify(lhs / ref)
    ok4 &= report(f"aux q-equation {i} proportional to 4th-order EOM (ratio={r})",
                  r.free_symbols <= {g, M} and not r.has(t))
print("    so in the (p, w) system of ym/: U(p, w) = V_tot(q) + w.grad V(q) - 1/2|w|^2, q = p - w,")
print("    V_tot = V - kappa(t) x y z; the Lee-Wick term keeps grad V (not grad V_tot): CHECK 4 confirms")
print("    that this is the correct reduction of L_YM + L_LW + L_axion.")

# ---------------------------------------------------------------- CHECK 5
stamp("CHECK 5: scaling to g = M = 1")
tau = sp.Symbol('tau', real=True)
u = [sp.Function(f'u{n}')(tau) for n in 'xyz']
sub = {q[i]: (M / g) * u[i] for i in range(3)}
Lax_s = sp.expand((kappa * x * y * z).subs(sub))
ok5 = report("kappa x y z = (M^4/g^2) (kappa/(g M)) u_x u_y u_z",
             sp.expand(Lax_s - M**4 / g**2 * (kappa / (g * M)) * u[0] * u[1] * u[2]) == 0)
print("    units g = M = 1: time in 1/M, fields in M/g; the drive enters as kappa/(g M) = phidot/(f M)")
print("    (= 2 xi H/M with Fujita's xi = phidot/(2 f H)), and Fujita's (13) reads u'' + 2u^3 = (phidot/(fM)) u^2.")

# ---------------------------------------------------------------- CHECK 6
stamp("CHECK 6: linearisation on x = y = z in the (q, w) variables, M and g kept")
# The axion drive sources the ghost: on shell w = (qddot + grad V)/M^2 (hdym_derive.py CHECK 4), and the
# driven background has qddot + grad V = kappa grad(xyz) != 0.  So w = W = 0 is no longer invariant and the
# isotropic background is (Q, Q, Q; w0, w0, w0) with w0 != 0.  Linearise the Euler-Lagrange equations of
# L' (CHECK 4) about it with Q(t), w0(t) free functions (the background equations are printed, not imposed),
# and project on the isotropic and anisotropic directions.
Qf = sp.Function('Q')(t)
w0 = sp.Function('w0')(t)
eps_ = sp.Function('epsilon')(t)      # anisotropic / isotropic q-perturbation amplitude
zet = sp.Function('zeta')(t)          # w-perturbation amplitude
lam = sp.Symbol('lambda_')
EL = [sp.expand(e.lhs) for e in aux_eqs]           # 3 q-equations, 3 w-equations of L' with kappa_t
subs_bg = {}
for i in range(3):
    subs_bg[q[i]] = Qf
    subs_bg[wv[i]] = w0
bg_eqs = [sp.expand(e.subs(subs_bg).doit()) for e in EL]
ok6 = report("isotropic truncation consistent: the three q-equations coincide, the three w-equations coincide",
             all(sp.expand(bg_eqs[i] - bg_eqs[0]) == 0 for i in range(3))
             and all(sp.expand(bg_eqs[3 + i] - bg_eqs[3]) == 0 for i in range(3)))
wsol_iso = sp.solve(bg_eqs[3], w0)[0]
print(f"    background w-equation: w0 = {sp.factor(wsol_iso)}")
iso4 = sp.expand(bg_eqs[0].subs(w0, wsol_iso).doit())
print(f"    background 4th-order isotropic equation (= 0): {sp.collect(iso4, kt)}")
lim = sp.expand(sp.limit(iso4, M, sp.oo))
ok6 &= report("M -> infinity: background equation is -(Q'' + 2 g^2 Q^3 - kappa Q^2)",
              sp.expand(lim + (Qf.diff(t, 2) + 2 * g**2 * Qf**3 - kt * Qf**2)) == 0)


def project(vec):
    """linearised equations along the direction vec: q = Q + lam*eps*vec, w = w0 + lam*zeta*vec,
    O(lam) terms, then the scalar equation along vec (all components equal up to the vector)."""
    sub = {}
    for i in range(3):
        sub[q[i]] = Qf + lam * eps_ * vec[i]
        sub[wv[i]] = w0 + lam * zet * vec[i]
    lin = [sp.expand(e.subs(sub).doit()).coeff(lam, 1) for e in EL]
    qeq = None
    weq = None
    for i in range(3):
        if vec[i] != 0:
            cq = sp.expand(lin[i] / vec[i])
            cw = sp.expand(lin[3 + i] / vec[i])
            if qeq is None:
                qeq, weq = cq, cw
            else:
                assert sp.expand(cq - qeq) == 0 and sp.expand(cw - weq) == 0, "not proportional to vec"
        else:
            assert sp.expand(lin[i]) == 0 and sp.expand(lin[3 + i]) == 0
    return qeq, weq


an_q, an_w = project([1, -1, 0])
an_q2, an_w2 = project([1, 1, -2])
ok6 &= report("the two anisotropic directions give the same equations",
              sp.expand(an_q - an_q2) == 0 and sp.expand(an_w - an_w2) == 0)
zsol = sp.solve(an_w, zet)[0]
print(f"    anisotropic w-equation: zeta = {sp.factor(zsol)}")
an4 = sp.expand(an_q.subs(zet, zsol).doit())
print(f"    anisotropic 4th-order equation (= 0): {an4}")
ok6 &= report("anisotropic sector:  eps''''/M^2 + eps'' + kappa(t) Q(t) eps = 0  (w0 drops out)",
              sp.expand(an4 + (eps_.diff(t, 4) / M**2 + eps_.diff(t, 2) + kt * Qf * eps_)) == 0)
print("    frozen coefficients: omega^2 = (M^2/2)(1 +- sqrt(1 - 4 kappa Q / M^2)): the Yang-Mills branch")
print("    omega^2 ~ kappa Q is tachyonic for kappa Q < 0, the ghost branch omega^2 ~ M^2 - kappa Q, and the")
print("    two merge into a complex pair for kappa Q > M^2/4.  Without the axion the sector is exactly free.")
is_q, is_w = project([1, 1, 1])
zsol_i = sp.solve(is_w, zet)[0]
print(f"    isotropic w-equation: zeta = {sp.factor(zsol_i)}")
is4 = sp.expand(is_q.subs(zet, zsol_i).doit())
print(f"    isotropic 4th-order equation (= 0): {sp.collect(is4, [eps_.diff(t, 4), eps_.diff(t, 2), eps_])}")
is4_lim = sp.expand(sp.limit(is4.subs(w0, wsol_iso).doit(), M, sp.oo))   # on the background, w0 = O(1/M^2)
ok6 &= report("isotropic sector, M -> infinity:  -(eps'' + (6 g^2 Q^2 - 2 kappa Q) eps)",
              sp.expand(is4_lim + (eps_.diff(t, 2) + (6 * g**2 * Qf**2 - 2 * kt * Qf) * eps_)) == 0)

# ---------------------------------------------------------------- CHECK 7
stamp("CHECK 7: FRW in conformal time")
a = sp.Function('a')(t)
Lax_frw = sp.expand(phi / (4 * f) * FFdual(a))
ok7 = report("sqrt(-g) (1/4f) phi F.Ftilde on a(eta)^2 eta equals the flat expression (no a)",
             sp.expand(Lax_frw - Lax_full) == 0)

# ---------------------------------------------------------------- summary
allok = all([ok1, ok2, ok3, ok4, ok5, ok6, ok7])
stamp(f"ALL CHECKS {'PASSED' if allok else 'NOT ALL PASSED'}")
print("""
Summary (units g = M = 1, time in 1/M, fields in M/g):
  the axion adds to the 12-d system of ym/ exactly the time-dependent potential  -kappa(t) x y z
  in q = p - w, kappa = (g/f) phidot (in flat space) = (g/f) dphi/d eta (on FRW, conformal time,
  comoving fields), and nothing else: Lee-Wick term, Gauss law and the diagonal truncation are as
  before.  On x = y = z = Q(t), as M -> infinity:  Q'' = -2Q^3 + kappa Q^2 (Fujita (13), H = 0);
  the drive sources the ghost (w = kappa Q^2 / M^2 + O(1/M^4) on the isotropic background), and the
  anisotropic sector obeys eps''''/M^2 + eps'' + kappa Q eps = 0 (frozen: YM branch ~ kappa Q, ghost
  branch ~ M^2 - kappa Q, complex pair for kappa Q > M^2/4).
""")
