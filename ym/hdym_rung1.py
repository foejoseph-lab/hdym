#!/usr/bin/env python3
"""
hdym_rung1.py -- Rung 1: linearisation of the FULL homogeneous SU(2) Lee-Wick system about the
colour-diagonal background A^a_i = delta^a_i q_i(t), A_0 = 0 (units g = M = 1).

Builds (symbolically, then exports to numpy)
  * the quadratic Lagrangian L2 of the six off-diagonal components u^a_i (a != i),
    including the Gauss-law piece (J_0)^2 of the Lee-Wick term, which vanishes on the
    ansatz but is quadratic in u off it;
  * the linearised Gauss constraints (A_0 equations at O(u));
  * the six linear 4th-order Euler-Lagrange equations for u, as a first-order system
    xdot = K(q, qdot) x of dimension 24 on the ghost-free background (qddot = -grad V);
  * the three exact gauge zero modes  u^a_i = eps^{abi} theta^b q_i  (constant theta).

Writes rung1_system.pkl (lambdified K, constraint matrix, zero-mode matrix, L2) for the
numerical stage (hdym_rung1_lyap.py) and prints the structural checks.
"""
import pickle
import time

import sympy as sp

T0 = time.time()


def stamp(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def report(name, ok):
    print(f"    {name}: {'PASS' if ok else 'FAIL'}", flush=True)
    return ok


t = sp.Symbol('t', real=True)
eps = sp.LeviCivita
eta = sp.diag(1, -1, -1, -1)
g = M = 1

# ---------------------------------------------------------------- full homogeneous Lagrangian, 12 functions
stamp("building the full homogeneous Lagrangian (12 functions)")
A = [[sp.Function(f'A{mu}{c}')(t) for c in range(3)] for mu in range(4)]


def d(mu, expr):
    return expr.diff(t) if mu == 0 else sp.Integer(0)


F = [[[None] * 3 for _ in range(4)] for _ in range(4)]
for mu in range(4):
    for nu in range(4):
        for c in range(3):
            F[mu][nu][c] = (d(mu, A[nu][c]) - d(nu, A[mu][c])
                            + sum(eps(c, b, e) * A[mu][b] * A[nu][e] for b in range(3) for e in range(3)))


def Dcov(mu, X):
    return [d(mu, X[c]) + sum(eps(c, b, e) * A[mu][b] * X[e] for b in range(3) for e in range(3))
            for c in range(3)]


J = []
for nu in range(4):
    acc = [sp.Integer(0)] * 3
    for mu in range(4):
        Dm = Dcov(mu, [F[mu][nu][c] for c in range(3)])
        acc = [acc[c] + eta[mu, mu] * Dm[c] for c in range(3)]
    J.append(acc)
LYM = sum(-sp.Rational(1, 4) * eta[mu, mu] * eta[nu, nu] * F[mu][nu][c] ** 2
          for mu in range(4) for nu in range(4) for c in range(3))
cj = sp.Symbol('c_j')   # multiplies the (D^mu F_{mu 0})^2 = J_0^2 piece; c_j = 1 is the theory, c_j = 0 a diagnostic
LHD = (cj * sum(J[0][c] ** 2 for c in range(3)) - sum(J[nu][c] ** 2 for nu in range(1, 4) for c in range(3))) / 2
L = sp.expand(LYM + LHD)
stamp(f"   {len(L.args)} terms")

# ---------------------------------------------------------------- perturbation about the diagonal background
stamp("perturbation ansatz")
x, y, z = [sp.Function(n)(t) for n in 'xyz']
q = [x, y, z]
e_ = sp.Symbol('epsilon')
# off-diagonal components u[i][a], a != i   (i = spatial index 0..2, a = colour 0..2)
u = {(i, a): sp.Function(f'u{i + 1}{a + 1}')(t) for i in range(3) for a in range(3) if i != a}
ulist = [u[(i, a)] for i in range(3) for a in range(3) if i != a]
sub = {}
for c in range(3):
    sub[A[0][c]] = sp.Integer(0)
    for i in range(3):
        sub[A[i + 1][c]] = (q[i] if i == c else sp.Integer(0)) + (e_ * u[(i, c)] if i != c else 0)
# Gauss law = A_0 equations, needed before A_0 is set to zero
gauss_full = [sp.euler_equations(L, [A[0][c]], t)[0].lhs for c in range(3)]
stamp("   A_0 equations derived")
Lp = sp.expand(L.subs(sub).doit())
L0 = Lp.coeff(e_, 0)
L1 = Lp.coeff(e_, 1)
L2 = Lp.coeff(e_, 2)
stamp(f"   L expanded: L2 has {len(L2.args)} terms")
V = (x ** 2 * y ** 2 + y ** 2 * z ** 2 + z ** 2 * x ** 2) / 2
gradV = [V.diff(v) for v in q]
Lclosed = (sum(v.diff(t) ** 2 for v in q) / 2 - V - sum((q[i].diff(t, 2) + gradV[i]) ** 2 for i in range(3)) / 2)
ok = report("O(eps^0): diagonal reduced Lagrangian", sp.simplify(L0 - sp.expand(Lclosed)) == 0)
# O(eps^1) must vanish up to a total derivative because the ansatz is a consistent truncation:
# check via its Euler-Lagrange expressions in u, which are then identically zero
ok &= report("O(eps^1): vanishes (consistent truncation)",
             L1 == 0 or all(sp.simplify(sp.euler_equations(L1, [uu], t)[0].lhs) == 0 for uu in ulist))

# ---------------------------------------------------------------- Gauss constraints at O(u)
gauss_lin = [sp.expand(gc.subs(sub).doit()).coeff(e_, 1) for gc in gauss_full]
ok &= report("Gauss law vanishes at O(eps^0)",
             all(sp.expand(gc.subs(sub).doit()).coeff(e_, 0) == 0 for gc in gauss_full))
stamp("   linearised Gauss constraints built")

# ---------------------------------------------------------------- linear equations for u
stamp("Euler-Lagrange equations of L2")
eqs = [sp.expand(sp.euler_equations(L2, [uu], t)[0].lhs) for uu in ulist]
# highest derivative structure: coefficient of u'''' in each equation
D4 = sp.Matrix(6, 6, lambda r, c: eqs[r].coeff(ulist[c].diff(t, 4)))
stamp(f"   4th-derivative coefficient matrix: {D4}")
ok &= report("u'''' coefficients form -I (Lee-Wick kinetic term, sign as in the diagonal sector)",
             D4 == -sp.eye(6))

# ---------------------------------------------------------------- first-order system on the ghost-free background
stamp("first-order form on the ghost-free background (qddot = -grad V)")
# state x = (u, u', u'', u''') for the 6 components; background symbols
X = sp.symbols('x y z')
P = sp.symbols('px py pz')
bg = {}
qs = list(X)
Vs = (X[0] ** 2 * X[1] ** 2 + X[1] ** 2 * X[2] ** 2 + X[2] ** 2 * X[0] ** 2) / 2
gV = [Vs.diff(v) for v in X]
HV = sp.hessian(Vs, X)
qdd = [-gV[i] for i in range(3)]
qddd = [-sum(HV[i, j] * P[j] for j in range(3)) for i in range(3)]
qdddd = [sp.expand(-sum(HV[i, j].diff(X[k]) * P[k] * P[j] for j in range(3) for k in range(3))
                   - sum(HV[i, j] * qdd[j] for j in range(3))) for i in range(3)]
for i in range(3):
    bg[q[i].diff(t, 4)] = qdddd[i]
    bg[q[i].diff(t, 3)] = qddd[i]
    bg[q[i].diff(t, 2)] = qdd[i]
    bg[q[i].diff(t)] = P[i]
    bg[q[i]] = X[i]
U = [[sp.Symbol(f'u{k}_{j}') for j in range(6)] for k in range(5)]   # U[k][j] = d^k u_j / dt^k
usub = {}
for j, uu in enumerate(ulist):
    for k in range(4, 0, -1):
        usub[uu.diff(t, k)] = U[k][j]
    usub[uu] = U[0][j]


def to_syms(expr):
    e = expr
    for k in range(4, -1, -1):
        e = e.subs({uu.diff(t, k) if k else uu: U[k][j] for j, uu in enumerate(ulist)})
    for k in range(4, -1, -1):
        e = e.subs({q[i].diff(t, k) if k else q[i]: bg[q[i].diff(t, k) if k else q[i]] for i in range(3)})
    return sp.expand(e)


eqs_s = [to_syms(e) for e in eqs]
# solve for u'''' : eq = -u'''' + rest  ->  u'''' = rest
u4 = [sp.expand(eqs_s[j] + U[4][j]) for j in range(6)]
ok &= report("equations are linear in the state", all(sp.Poly(u4[j], *[s for row in U[:4] for s in row]).total_degree() == 1 for j in range(6)))
state = [s for row in U[:4] for s in row]          # 24
K = sp.zeros(24, 24)
for j in range(6):
    for k in range(3):
        K[k * 6 + j, (k + 1) * 6 + j] = 1
    for c, s in enumerate(state):
        K[18 + j, c] = u4[j].coeff(s)
gauss_s = [to_syms(gc) for gc in gauss_lin]
G = sp.Matrix(3, 24, lambda r, c: gauss_s[r].coeff(state[c]))
ok &= report("linearised Gauss constraints are linear in the state and depend on u up to u'''",
             all(sp.expand(gauss_s[r] - sum(G[r, c] * state[c] for c in range(24))) == 0 for r in range(3)))
# zero modes as state vectors: u = theta x A and its derivatives
Zm = sp.zeros(24, 3)
for b in range(3):
    for j, (i, a) in enumerate([(i, a) for i in range(3) for a in range(3) if i != a]):
        val = eps(a, b, i) * q[i]
        for k in range(4):
            Zm[k * 6 + j, b] = to_syms(val.diff(t, k)) if k else to_syms(val)
ok &= report("gauge zero modes satisfy the linearised Gauss law", sp.simplify(G * Zm) == sp.zeros(3, 3))
# the zero modes must be solutions: d/dt Zm along the ghost-free flow == K Zm
Zdot = sp.zeros(24, 3)
for r in range(24):
    for c in range(3):
        Zdot[r, c] = sum(Zm[r, c].diff(X[i]) * P[i] + Zm[r, c].diff(P[i]) * qdd[i] for i in range(3))
ok &= report("u = theta x A (constant theta) solves the linear equations on the ghost-free background",
             sp.simplify(sp.expand(K * Zm - Zdot)) == sp.zeros(24, 3))
stamp(f"ALL STRUCTURAL CHECKS {'PASSED' if ok else 'NOT ALL PASSED'}")

# ---------------------------------------------------------------- quadratic Lagrangian in readable form
L2s = to_syms(L2)
with open('rung1_system.pkl', 'wb') as fh:
    pickle.dump(dict(K=K, G=G, Zm=Zm, L2=L2s, state=state, X=X, P=P, cj=cj, ulist=[str(s) for s in ulist]), fh)
stamp("wrote rung1_system.pkl")
print("\nL2 (background symbols x,y,z, px,py,pz; u_k_j = d^k u_j/dt^k; u order:", [str(s) for s in ulist], ")")
print(sp.collect(L2s, [s for row in U for s in row]))
