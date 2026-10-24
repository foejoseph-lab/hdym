#!/usr/bin/env python3
"""
dilaton_derive.py -- plan item 5b, step 1: what a gauge-kinetic (dilaton) coupling Theta(phi) F^2
adds to the homogeneous SU(2) Lee-Wick system, derived symbolically with PASS/FAIL checks, in two
columns side by side:

    column A  the dilaton multiplies F^2 only          (Adshead-Giblin-Weiner 1708.02944 eq. (1),
                                                        Obata-Soda 1602.06024 eq. (3): the cited form)
    column B  the dilaton multiplies the whole operator, Lee-Wick term included
                                                       (string-frame-natural; the control -- A primary, B
                                                        control, decision of 17 Oct, README)

Theory (metric eta = diag(+1,-1,-1,-1), same construction as ym/hdym_derive.py and
axion_su2/code/axion_derive.py):

    L_A = -Theta/4 F^a_{mu nu} F^{a mu nu} + 1/(2 M^2) (D^mu F_{mu nu})^a (D_rho F^{rho nu})^a
    L_B = Theta [ -1/4 F^a_{mu nu} F^{a mu nu} + 1/(2 M^2) (D^mu F_{mu nu})^a (D_rho F^{rho nu})^a ]

Theta(t) = exp(delta(t)), delta = phi/Lambda, the homogeneous dilaton treated as a given function of
time (Adshead et al. write exp(phi/M) and -Theta/2 Tr F^2; their M is the coupling scale, written
Lambda here because M is the Lee-Wick parameter; the trace normalisation only rescales Lambda).
Theta is kept exact (not expanded in delta): Adshead et al. have phi_end ~ 0.2 m_pl, so delta need
not be small.  This is the parity-even twin of the axion coupling: the axion adds a potential
-kappa(t) x y z that vanishes with the field; the dilaton multiplies the operator and does not.

Checks (each prints PASS/FAIL; "both columns" means two PASS lines):

  CHECK 1  Homogeneous reduction, all 12 functions A^a_mu(t), A_0 kept.  Column A: the A_0 equation
           of Theta L_YM is Theta times the Yang-Mills Gauss law (L_YM has no A_0-dot), so the Gauss law
           of the full theory is unchanged in content.  Column B: the A_0 equation picks up
           Theta-dot dL_LW/dA_0-dot, which is printed on the diagonal ansatz (it vanishes there) --
           off the ansatz the ghost sector's constraint carries the dilaton, which ym/hdym_rung1.py's
           generator would have to be regenerated with.
  CHECK 2  The diagonal ansatz A^a_i = delta^a_i q_i(t), A_0 = 0 stays a consistent truncation in both
           columns (off-diagonal equations and Gauss law vanish on it), and on it
             L_A = Theta (1/2|qdot|^2 - V) - |qddot + grad V|^2/(2M^2),
             L_B = Theta (1/2|qdot|^2 - V - |qddot + grad V|^2/(2M^2)),
           with the Euler-Lagrange equations of the closed forms equal to the restricted full ones
           (symmetric criticality, verified).
  CHECK 3  Auxiliary-field form (ym/hdym_derive.py CHECK 4, L' = 1/2|qdot|^2 + wdot.qdot - V - w.grad V
           + M^2/2 |w|^2):
             column A  L'_A = Theta (1/2|qdot|^2 - V) + wdot.qdot - w.grad V + M^2/2 |w|^2
             column B  L'_B = Theta L' + Theta-dot w.qdot
           reproduce the 4th-order equations.  In the (q, w) system of ym/ the dilaton of column A
           is the factor Theta on the normal sector's kinetic term and potential, and nothing else;
           on shell w = (qddot + grad V)/M^2 in both columns.
  CHECK 4  Scaling t -> tau/M, q -> (M/g) u (hdym_derive.py CHECK 5): Theta is dimensionless, so in
           units g = M = 1 the new quantities are delta(t) = phi/Lambda (amplitude delta_a = phi_amp/Lambda)
           and its frequency m_phi/M, and nothing else.
  CHECK 5  M -> infinity (both columns): in the time tau with d tau = dt/Theta the equations are
           Yang-Mills with g -> g Theta(t),  q_tautau = -Theta^2 grad V,  on x = y = z: Q_tautau = -2 g^2
           Theta^2 Q^3.  (Consequence, not a check: for a fast dilaton the leading averaged coupling is the
           tau-average <Theta^2>_tau = int e^delta dt / int e^-delta dt, which is positive, and exactly 1 when
           delta is symmetric about 0 (a cosine drive): the M -> infinity control never inverts the
           potential, and at leading order does not even renormalise it; its only loss of the isotropic
           solution is the slow-drive parametric resonance of Q itself.)
  CHECK 6  Linearisation on x = y = z in the (q, w) variables, M and g kept, both columns.  The
           dilaton sources the ghost: w = (qddot + grad V)/M^2 and the driven background has
           qddot + grad V = -(Theta-dot/Theta) qdot + ... != 0, so w = W = 0 is not invariant and the
           isotropic background is (Q, Q, Q; w0, w0, w0).  Checks: the isotropic truncation is
           consistent; the background's 4th-order equation tends to Theta (Q'' + 2 g^2 Q^3) +
           Theta-dot Q' = 0 as M -> infinity; the two anisotropic directions (permutation-odd,
           degenerate) obey, with Q and w0 dropping out exactly,
             column A   eps''''/M^2 + Theta eps'' + Theta-dot eps' = 0,
             column B   (Theta eps'')''/M^2 + (Theta eps')' = 0,
           i.e. with u = eps' the first integrals
             column A   u'' + M^2 Theta(t) u = C          (the ghost's effective mass M^2 Theta is modulated),
             column B   (Theta u')' + M^2 Theta u = C M^2   (the ghost's damping delta-dot is modulated;
                                                            u = exp(-delta/2) v gives v'' + (M^2 - delta''/2
                                                            - delta'^2/4) v = C M^2 exp(-delta/2)).
           So the anisotropic sector does not see the background at all -- there is nothing for the
           axion's Pi-cancellation to act on -- and its only instabilities are the ghost's parametric
           tongues, first order in delta_a and not suppressed by Q (printed: the Mathieu parameters to
           first order in delta_a; heuristics, the runs decide).  The isotropic tangent equation is
           printed; as M -> infinity it is Theta eps'' + Theta-dot eps' + 6 g^2 Theta Q^2 eps = 0.
  CHECK 7  FRW, ds^2 = a(eta)^2 (deta^2 - dx^2): sqrt(-g) Theta F.F is independent of a(eta) (Yang-Mills
           is conformal), so in conformal time the dilaton term is the flat one with comoving fields,
           and stage 2 of axion_su2/ (M -> M a(eta)) carries over to both columns.

Usage (from the repository root, about a minute):
    python dilaton_su2/code/dilaton_derive.py > dilaton_su2/results/dilaton_derive.txt
Added 2026-10-16.  The runs the closing summary calls "dilaton_dressed.py" are dilaton_iso.py (17 Oct): the
isotropic loss is the ghost's parametric tongue at m_phi = 2M/n, not a threshold and not a pump (README, Result).
"""
import time

import sympy as sp

T0 = time.time()


def stamp(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def report(name, ok):
    print(f"    {name}: {'PASS' if ok else 'FAIL'}", flush=True)
    return bool(ok)


COLS = ('A', 'B')
t = sp.Symbol('t', real=True)
g, M = sp.symbols('g M', positive=True)
eps = sp.LeviCivita
eta = sp.diag(1, -1, -1, -1)
Th = sp.Function('Theta', positive=True)(t)       # Theta(t) = exp(delta(t)), a given function

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


def Dcov(mu, X):
    return [d(mu, X[c]) + g * sum(eps(c, b, e) * A[mu][b] * X[e]
                                  for b in range(3) for e in range(3))
            for c in range(3)]


J = []
for nu in range(4):
    acc = [sp.Integer(0)] * 3
    for mu in range(4):
        col = [F[mu][nu][c] for c in range(3)]
        Dm = Dcov(mu, col)
        acc = [acc[c] + eta[mu, mu] * Dm[c] for c in range(3)]
    J.append(acc)


def FF(metric_scale):
    """sqrt(-g) F^a_{mu nu} F^{a mu nu} on g_{mu nu} = metric_scale^2 eta."""
    s = metric_scale
    return sp.expand(s**4 * sum((eta[mu, mu] * eta[nu, nu] / s**4) * F[mu][nu][c]**2
                                for mu in range(4) for nu in range(4) for c in range(3)))


LYM = sp.expand(-sp.Rational(1, 4) * FF(sp.Integer(1)))
LHD = sp.expand(sum(eta[nu, nu] * J[nu][c]**2 for nu in range(4) for c in range(3)) / (2 * M**2))
stamp(f"Lagrangians built: L_YM {len(LYM.args)} terms, L_LW {len(LHD.args)} terms")

funcs = [A[mu][c] for mu in range(4) for c in range(3)]
# Euler-Lagrange is linear in L: three calls on the pieces instead of two on the products
EL_YMth = [e.lhs for e in sp.euler_equations(sp.expand(Th * LYM), funcs, t)]
stamp("EL of Theta L_YM")
EL_LW = [e.lhs for e in sp.euler_equations(LHD, funcs, t)]
stamp("EL of L_LW")
EL_LWth = [e.lhs for e in sp.euler_equations(sp.expand(Th * LHD), funcs, t)]
stamp("EL of Theta L_LW")
EOM = {'A': {(mu, c): EL_YMth[3 * mu + c] + EL_LW[3 * mu + c] for mu in range(4) for c in range(3)},
       'B': {(mu, c): EL_YMth[3 * mu + c] + EL_LWth[3 * mu + c] for mu in range(4) for c in range(3)}}

x, y, z = [sp.Function(n)(t) for n in 'xyz']
q = [x, y, z]
ans = {}
for c in range(3):
    ans[A[0][c]] = sp.Integer(0)
    for i in range(3):
        ans[A[i + 1][c]] = q[i] if i == c else sp.Integer(0)


def on_ansatz(expr):
    return sp.expand(expr.subs(ans).doit())


# ---------------------------------------------------------------- CHECK 1
stamp("CHECK 1: the Gauss law with the dilaton")
ok1 = report("L_YM contains no A_0-dot (column A's A_0 equation is Theta x the Yang-Mills Gauss law + the LW one)",
             not any(LYM.has(A[0][c].diff(t)) for c in range(3)))
# column A: A_0 equation of Theta L_YM equals Theta times that of L_YM, identically (all 12 functions free)
EL_YM = [e.lhs for e in sp.euler_equations(LYM, [A[0][c] for c in range(3)], t)]
ok1 &= report("column A: A_0 equation of Theta L_YM = Theta x (A_0 equation of L_YM), identically",
              all(sp.expand(EL_YMth[c] - Th * EL_YM[c]) == 0 for c in range(3)))
# column B: the extra piece is -Theta-dot dL_LW/dA_0-dot (from -d/dt dL/dA_0-dot)
extra = [sp.expand(EL_LWth[c] - Th * EL_LW[c]) for c in range(3)]
pred = [sp.expand(-Th.diff(t) * LHD.diff(A[0][c].diff(t))) for c in range(3)]
ok1 &= report("column B: A_0 equation of Theta L_LW = Theta x (A_0 equation of L_LW) - Theta-dot dL_LW/dA_0-dot",
              all(sp.expand(extra[c] - pred[c]) == 0 for c in range(3)))
print("    column B, the extra constraint term dL_LW/dA_0-dot on the diagonal ansatz:",
      [on_ansatz(LHD.diff(A[0][c].diff(t))) for c in range(3)])
print("    (off the ansatz it is not zero: the constraint of the ghost sector carries Theta-dot in column B)")

# ---------------------------------------------------------------- CHECK 2
stamp("CHECK 2: diagonal ansatz stays a consistent truncation; closed forms")
V = g**2 / 2 * (x**2 * y**2 + y**2 * z**2 + z**2 * x**2)
gradV = [V.diff(v) for v in q]
LW_closed = -sum((q[i].diff(t, 2) + gradV[i])**2 for i in range(3)) / (2 * M**2)
YM_closed = sp.Rational(1, 2) * sum(v.diff(t)**2 for v in q) - V
Lclosed = {'A': sp.expand(Th * YM_closed + LW_closed), 'B': sp.expand(Th * (YM_closed + LW_closed))}
Lfull = {'A': sp.expand(Th * LYM + LHD), 'B': sp.expand(Th * (LYM + LHD))}
ok2 = True
red_eqs = {}
for col in COLS:
    E = EOM[col]
    okc = all(on_ansatz(E[(0, c)]) == 0 for c in range(3))
    ok2 &= report(f"column {col}: Gauss law vanishes on the ansatz", okc)
    okc = all(on_ansatz(E[(i + 1, c)]) == 0 for i in range(3) for c in range(3) if i != c)
    ok2 &= report(f"column {col}: off-diagonal equations vanish on the ansatz", okc)
    ok2 &= report(f"column {col}: L|ansatz == closed form", sp.expand(on_ansatz(Lfull[col]) - Lclosed[col]) == 0)
    red_eqs[col] = [sp.expand(e.lhs) for e in sp.euler_equations(Lclosed[col], q, t)]
    okc = all(sp.expand(red_eqs[col][i] - on_ansatz(E[(i + 1, i)])) == 0 for i in range(3))
    ok2 &= report(f"column {col}: reduced EL equations == restricted full equations", okc)

# ---------------------------------------------------------------- CHECK 3
stamp("CHECK 3: auxiliary-field forms")
wv = [sp.Function(f'w{n}')(t) for n in 'xyz']
Lp0 = (sp.Rational(1, 2) * sum(v.diff(t)**2 for v in q)
       + sum(wv[i].diff(t) * q[i].diff(t) for i in range(3))
       - V - sum(wv[i] * gradV[i] for i in range(3))
       + M**2 / 2 * sum(wv[i]**2 for i in range(3)))
Lp = {'A': sp.expand(Th * YM_closed + sum(wv[i].diff(t) * q[i].diff(t) for i in range(3))
                     - sum(wv[i] * gradV[i] for i in range(3)) + M**2 / 2 * sum(wv[i]**2 for i in range(3))),
      'B': sp.expand(Th * Lp0 + Th.diff(t) * sum(wv[i] * q[i].diff(t) for i in range(3)))}
ok3 = True
aux_eqs = {}
for col in COLS:
    aux_eqs[col] = [sp.expand(e.lhs) for e in sp.euler_equations(Lp[col], q + wv, t)]
    wsol = sp.solve(aux_eqs[col][3:], wv, dict=True)[0]
    ok3 &= report(f"column {col}: on shell w = (qddot + grad V)/M^2",
                  all(sp.expand(wsol[wv[i]] - (q[i].diff(t, 2) + gradV[i]) / M**2) == 0 for i in range(3)))
    for i in range(3):
        lhs = sp.expand(aux_eqs[col][i].subs(wsol).doit())
        r = sp.simplify(lhs / red_eqs[col][i])
        ok3 &= report(f"column {col}: aux q-equation {i} proportional to the 4th-order equation (ratio={r})",
                      r.free_symbols <= {g, M} and not r.has(t))
print("    column A in the (q, w) system of ym/: q-equation  Theta qddot + Theta-dot qdot + Theta grad V"
      " + wddot + Hess V w = 0,")
print("    w-equation  M^2 w = qddot + grad V; first-order-ready:  qddot = M^2 w - grad V,")
print("    wddot = -Theta-dot qdot - (Theta M^2 + Hess V) w.  On x = y = z (g = M = 1):")
print("      Q'' = w - 2Q^3,   w'' + (Theta + 6Q^2) w = -Theta-dot Q'      [compare axion: w'' + (1 + 6Q^2) w = kappa Q^2].")

# ---------------------------------------------------------------- CHECK 4
stamp("CHECK 4: scaling to g = M = 1")
tau = sp.Symbol('tau', real=True)
u = [sp.Function(f'u{n}')(tau) for n in 'xyz']
Thtau = sp.Function('Theta_tau', positive=True)(tau)
sub = {}
for i in range(3):
    sub[q[i].diff(t, 2)] = (M / g) * M**2 * u[i].diff(tau, 2)
    sub[q[i].diff(t)] = (M / g) * M * u[i].diff(tau)
    sub[q[i]] = (M / g) * u[i]
Vu = (u[0]**2 * u[1]**2 + u[1]**2 * u[2]**2 + u[2]**2 * u[0]**2) / 2
gradVu = [Vu.diff(v) for v in u]
ok4 = True
for col in COLS:
    Ls = sp.expand(Lclosed[col].subs(sub).subs(Th, Thtau))
    YMu = sp.Rational(1, 2) * sum(v.diff(tau)**2 for v in u) - Vu
    LWu = -sum((u[i].diff(tau, 2) + gradVu[i])**2 for i in range(3)) / 2
    target = Thtau * YMu + LWu if col == 'A' else Thtau * (YMu + LWu)
    ok4 &= report(f"column {col}: L = (M^4/g^2) L~(u, Theta) with g = M = 1 and Theta untouched",
                  sp.expand(Ls - M**4 / g**2 * target) == 0)
print("    units g = M = 1: time in 1/M, fields in M/g; the dilaton enters as delta(t) = phi/Lambda, amplitude")
print("    delta_a = phi_amp/Lambda and frequency omega = m_phi/M; no other new quantity.")

# ---------------------------------------------------------------- CHECK 5
stamp("CHECK 5: M -> infinity is Yang-Mills with g -> g Theta in the time d tau = dt / Theta")
LYMth = sp.expand(Th * YM_closed)                       # both columns at M -> infinity
eqs_inf = [sp.expand(e.lhs) for e in sp.euler_equations(LYMth, q, t)]
# q(t) = U(tau(t)) with d tau/dt = 1/Theta:  qdot = U'/Theta,  qddot = U''/Theta^2 - U' Theta-dot/Theta^2
Uf = [sp.Function(f'U{n}')(tau) for n in 'xyz']
ok5 = True
for i in range(3):
    e = eqs_inf[i]
    subs_t = {}
    for j in range(3):
        subs_t[q[j].diff(t, 2)] = Uf[j].diff(tau, 2) / Th**2 - Uf[j].diff(tau) * Th.diff(t) / Th**2
        subs_t[q[j].diff(t)] = Uf[j].diff(tau) / Th
        subs_t[q[j]] = Uf[j]
    e_tau = sp.expand(e.subs(subs_t))
    # expected: -(U'' + Theta^2 grad V(U)) / Theta  (EL sign: dL/dq - d/dt dL/dqdot)
    gV = [V.diff(v).subs({q[k]: Uf[k] for k in range(3)}) for v in q]
    ok5 &= report(f"q_{i}: equation in tau is  U'' = -Theta^2 grad V(U)",
                  sp.expand(e_tau + (Uf[i].diff(tau, 2) + Th**2 * gV[i]) / Th) == 0)
print("    on x = y = z:  Q_tautau = -2 g^2 Theta^2 Q^3.  For a fast dilaton the leading averaged coupling is the")
print("    tau-average <Theta^2>_tau = int e^delta dt / int e^-delta dt: positive, and exactly 1 for delta symmetric")
print("    about 0 (a cosine drive).  So the M -> infinity control never inverts the potential and at leading order")
print("    does not renormalise it; its only loss of the isotropic solution is the slow-drive parametric resonance")
print("    of Q itself, as in axion_su2/ (B').  Clock: dt = Theta d tau, <Theta>_tau = 1/I_0(delta_a) for the cosine.")

# ---------------------------------------------------------------- CHECK 6
stamp("CHECK 6: linearisation on x = y = z in the (q, w) variables, M and g kept")
Qf = sp.Function('Q')(t)
w0 = sp.Function('w0')(t)
eps_ = sp.Function('epsilon')(t)
zet = sp.Function('zeta')(t)
lam = sp.Symbol('lambda_')
ok6 = True
subs_bg = {}
for i in range(3):
    subs_bg[q[i]] = Qf
    subs_bg[wv[i]] = w0


def project(EL, vec):
    sub = {}
    for i in range(3):
        sub[q[i]] = Qf + lam * eps_ * vec[i]
        sub[wv[i]] = w0 + lam * zet * vec[i]
    lin = [sp.expand(e.subs(sub).doit()).coeff(lam, 1) for e in EL]
    qeq = weq = None
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


aniso4 = {}
for col in COLS:
    print(f"  -- column {col}")
    EL = aux_eqs[col]
    bg = [sp.expand(e.subs(subs_bg).doit()) for e in EL]
    ok6 &= report(f"column {col}: isotropic truncation consistent (q-equations coincide, w-equations coincide)",
                  all(sp.expand(bg[i] - bg[0]) == 0 for i in range(3))
                  and all(sp.expand(bg[3 + i] - bg[3]) == 0 for i in range(3)))
    wsol_iso = sp.solve(bg[3], w0)[0]
    print(f"    background w-equation: w0 = {sp.factor(wsol_iso)}")
    iso4 = sp.expand(bg[0].subs(w0, wsol_iso).doit())
    print(f"    background 4th-order isotropic equation (= 0): {sp.collect(iso4, [Th, Th.diff(t)])}")
    lim = sp.expand(sp.limit(iso4, M, sp.oo))
    ok6 &= report(f"column {col}: M -> infinity: background equation is -(Theta (Q'' + 2 g^2 Q^3) + Theta-dot Q')",
                  sp.expand(lim + Th * (Qf.diff(t, 2) + 2 * g**2 * Qf**3) + Th.diff(t) * Qf.diff(t)) == 0)
    an_q, an_w = project(EL, [1, -1, 0])
    an_q2, an_w2 = project(EL, [1, 1, -2])
    ok6 &= report(f"column {col}: the two anisotropic directions give the same equations",
                  sp.expand(an_q - an_q2) == 0 and sp.expand(an_w - an_w2) == 0)
    zsol = sp.solve(an_w, zet)[0]
    print(f"    anisotropic w-equation: zeta = {sp.factor(zsol)}")
    an4 = sp.expand(an_q.subs(zet, zsol).doit())
    aniso4[col] = an4
    print(f"    anisotropic 4th-order equation (= 0): {an4}")
    ok6 &= report(f"column {col}: anisotropic equation contains neither Q nor w0", not an4.has(Qf) and not an4.has(w0))
    if col == 'A':
        ref = eps_.diff(t, 4) / M**2 + (Th * eps_.diff(t)).diff(t)
        ok6 &= report("column A: anisotropic sector  eps''''/M^2 + (Theta eps')' = 0", sp.expand(an4 + ref) == 0)
        # first integral with u = eps':  u''/M^2 + Theta u = C
        uf = sp.Function('u')(t)
        first = uf.diff(t, 2) / M**2 + Th * uf
        ok6 &= report("column A: d/dt [u''/M^2 + Theta u] reproduces the equation (u = eps')",
                      sp.expand(first.diff(t).subs(uf, eps_.diff(t)).doit() + an4) == 0)
    else:
        ref = (Th * eps_.diff(t, 2)).diff(t, 2) / M**2 + (Th * eps_.diff(t)).diff(t)
        ok6 &= report("column B: anisotropic sector  (Theta eps'')''/M^2 + (Theta eps')' = 0", sp.expand(an4 + ref) == 0)
        uf = sp.Function('u')(t)
        first = (Th * uf.diff(t)).diff(t) / M**2 + Th * uf
        ok6 &= report("column B: d/dt [(Theta u')'/M^2 + Theta u] reproduces the equation (u = eps')",
                      sp.expand(first.diff(t).subs(uf, eps_.diff(t)).doit() + an4) == 0)
        # u = exp(-delta/2) v removes the damping: v'' + (M^2 - delta''/2 - delta'^2/4) v = C M^2 exp(-delta/2)
        dl = sp.Function('delta')(t)
        vf = sp.Function('v')(t)
        lhsB = ((sp.exp(dl) * (sp.exp(-dl / 2) * vf).diff(t)).diff(t) / M**2 + sp.exp(dl) * sp.exp(-dl / 2) * vf)
        lhsB = sp.expand(lhsB * M**2 * sp.exp(-dl / 2))
        ok6 &= report("column B: u = exp(-delta/2) v gives  v'' + (M^2 - delta''/2 - delta'^2/4) v = C M^2 exp(-delta/2)",
                      sp.simplify(lhsB - (vf.diff(t, 2) + (M**2 - dl.diff(t, 2) / 2 - dl.diff(t)**2 / 4) * vf)) == 0)
    is_q, is_w = project(EL, [1, 1, 1])
    zsol_i = sp.solve(is_w, zet)[0]
    print(f"    isotropic w-equation: zeta = {sp.factor(zsol_i)}")
    is4 = sp.expand(is_q.subs(zet, zsol_i).doit())
    print(f"    isotropic 4th-order equation (= 0): {sp.collect(is4, [eps_.diff(t, 4), eps_.diff(t, 3), eps_.diff(t, 2), eps_.diff(t), eps_])}")
    is4_lim = sp.expand(sp.limit(is4.subs(w0, wsol_iso).doit(), M, sp.oo))
    ok6 &= report(f"column {col}: isotropic sector, M -> infinity:  -(Theta eps'' + Theta-dot eps' + 6 g^2 Theta Q^2 eps)",
                  sp.expand(is4_lim + Th * eps_.diff(t, 2) + Th.diff(t) * eps_.diff(t) + 6 * g**2 * Th * Qf**2 * eps_) == 0)

print("  -- the anisotropic ghost as a Mathieu problem, first order in delta_a (heuristic; the runs decide).")
print("     delta = delta_a cos(omega t), z = omega t / 2, Mathieu  u_zz + (a - 2 q cos 2z) u = 0:")
print("       column A  u'' + M^2 e^delta u = C:          a = 4M^2/omega^2,  q = -2 M^2 delta_a / omega^2;")
print("                 first tongue omega = 2M: q = -delta_a/2, growth delta_a M / 4 per unit time;")
print("       column B  v'' + (M^2 - delta''/2 - ...) v:  a = 4M^2/omega^2 - delta_a^2/2,  q = -delta_a;")
print("                 first tongue omega = 2M: growth delta_a M / 2 per unit time (twice column A).")
print("     Both are first order in delta_a with no factor of Q (the axion's first tongue has c1/4, c1 = kappa Q).")
print("     Off the tongues (omega != 2M/n) the sector is stable with eps drifting linearly (the C term).")

# ---------------------------------------------------------------- CHECK 7
stamp("CHECK 7: FRW in conformal time")
a = sp.Function('a')(t)
ok7 = report("sqrt(-g) Theta F.F on a(eta)^2 eta equals the flat expression (no a)",
             sp.expand(Th * FF(a) - Th * FF(sp.Integer(1))) == 0)
print("    so on FRW the dilaton drive is Theta(eta) (-1/4 F.F)_flat with comoving fields in both columns, and the")
print("    Lee-Wick term is the flat one with M -> M a(eta) (axion_su2/code/frw_identity.py); column B multiplies")
print("    that by Theta(eta) too.")

# ---------------------------------------------------------------- summary
allok = all([ok1, ok2, ok3, ok4, ok5, ok6, ok7])
stamp(f"ALL CHECKS {'PASSED' if allok else 'NOT ALL PASSED'}")
print("""
Summary (units g = M = 1, time in 1/M, fields in M/g; Theta = exp(delta), delta = phi/Lambda prescribed):
  column A (F^2 only, the cited coupling): the dilaton is the factor Theta(t) on the normal sector's kinetic
    term and potential in the (q, w) system of ym/ and nothing else; Gauss law and the diagonal truncation
    unchanged; on x = y = z:  Q'' = w - 2Q^3,  w'' + (Theta + 6Q^2) w = -Theta-dot Q'.
  column B (whole operator): L'_B = Theta L' + Theta-dot w.qdot; the ghost sector's constraint carries
    Theta-dot off the diagonal ansatz (the off-diagonal generator would need regenerating).
  Both: M -> infinity is Yang-Mills with g -> g Theta in d tau = dt/Theta (never an inverted potential);
    the anisotropic sector is independent of the background -- A: u'' + Theta u = C, B: (Theta u')' + Theta u = C
    (u = eps') -- so the axion's Pi-cancellation has nothing to act on there, and the sector's only
    instabilities are the ghost's parametric tongues at m_phi = 2M/n, first order in delta_a and not
    suppressed by Q.  What is left to measure is the isotropic sector: the background's own response to
    the drive through the full propagator (the axion's threshold lambda_eff = 0 had its analogue here
    only if the second-order DC term is a potential term; the drive -Theta-dot Q' is a friction, so it
    may instead be a secular energy pump).  That is dilaton_dressed.py.
""")
