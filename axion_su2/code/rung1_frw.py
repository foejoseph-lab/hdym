#!/usr/bin/env python3
"""
rung1_frw.py -- ym/hdym_rung1.py on spatially flat FRW in conformal time: the linearisation of the
full homogeneous SU(2) Lee-Wick system about the colour-diagonal background A^a_i = delta^a_i q_i(eta),
A_0 = 0, with the Lee-Wick term carrying 1/a(eta)^2 (M -> M a(eta), see frw_identity.py).  Units
g = M = 1, M the fixed Lee-Wick parameter; a(eta) is a free function.

Same construction as ym/hdym_rung1.py (which is the a = 1 case): quadratic Lagrangian L2 of the six
off-diagonal components, linearised Gauss constraints, the six linear 4th-order equations as a
first-order system x' = K(q, q', c_j; a, a', a'') x of dimension 24 on the ghost-free background
(q'' = -grad V: the background is pure Yang-Mills, hence conformal, hence unchanged), and the three
gauge zero modes u^a_i = eps^{abi} theta^b q_i.

Checks printed (in addition to those of hdym_rung1.py):
  * the u'''' coefficients are -I / a^2;
  * only a, a', a'' occur in K and in the Gauss matrix (no a''' or a'''');
  * K, G at a = 1, a' = a'' = 0 equal the committed ym/rung1_system.pkl exactly (symbolically);
  * the gauge zero modes solve the linear equations for arbitrary a(eta) (gauge symmetry is intact).

Writes rung1_frw_system.pkl next to this folder's README (lambdify over X + P + [c_j, a, a1, a2]).
Usage (from the repository root, ~15 s):
    python axion_su2/code/rung1_frw.py > axion_su2/results/rung1_frw.txt
Added 2026-10-14.
"""
import os
import pickle
import time

import sympy as sp

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'rung1_frw_system.pkl')
REF = os.path.join(HERE, '..', '..', 'ym', 'rung1_system.pkl')


def stamp(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def report(name, ok):
    print(f"    {name}: {'PASS' if ok else 'FAIL'}", flush=True)
    return ok


t = sp.Symbol('t', real=True)          # conformal time eta
eps = sp.LeviCivita
eta = sp.diag(1, -1, -1, -1)
a = sp.Function('a')(t)

# ---------------------------------------------------------------- full homogeneous Lagrangian, 12 functions
stamp("building the full homogeneous Lagrangian (12 functions), Lee-Wick term / a(eta)^2")
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
cj = sp.Symbol('c_j')
LHD = (cj * sum(J[0][c] ** 2 for c in range(3)) - sum(J[nu][c] ** 2 for nu in range(1, 4) for c in range(3))) / (2 * a ** 2)
L = sp.expand(LYM + LHD)
stamp(f"   {len(L.args)} terms")

# ---------------------------------------------------------------- perturbation about the diagonal background
stamp("perturbation ansatz")
x, y, z = [sp.Function(n)(t) for n in 'xyz']
q = [x, y, z]
e_ = sp.Symbol('epsilon')
u = {(i, aa): sp.Function(f'u{i + 1}{aa + 1}')(t) for i in range(3) for aa in range(3) if i != aa}
ulist = [u[(i, aa)] for i in range(3) for aa in range(3) if i != aa]
sub = {}
for c in range(3):
    sub[A[0][c]] = sp.Integer(0)
    for i in range(3):
        sub[A[i + 1][c]] = (q[i] if i == c else sp.Integer(0)) + (e_ * u[(i, c)] if i != c else 0)
gauss_full = [sp.euler_equations(L, [A[0][c]], t)[0].lhs for c in range(3)]
stamp("   A_0 equations derived")
Lp = sp.expand(L.subs(sub).doit())
L0 = Lp.coeff(e_, 0)
L1 = Lp.coeff(e_, 1)
L2 = Lp.coeff(e_, 2)
stamp(f"   L expanded: L2 has {len(L2.args)} terms")
V = (x ** 2 * y ** 2 + y ** 2 * z ** 2 + z ** 2 * x ** 2) / 2
gradV = [V.diff(v) for v in q]
Lclosed = (sum(v.diff(t) ** 2 for v in q) / 2 - V
           - sum((q[i].diff(t, 2) + gradV[i]) ** 2 for i in range(3)) / (2 * a ** 2))
ok = report("O(eps^0): diagonal reduced Lagrangian with M -> M a(eta)", sp.simplify(L0 - sp.expand(Lclosed)) == 0)
ok &= report("O(eps^1): vanishes (consistent truncation)",
             L1 == 0 or all(sp.simplify(sp.euler_equations(L1, [uu], t)[0].lhs) == 0 for uu in ulist))
gauss_lin = [sp.expand(gc.subs(sub).doit()).coeff(e_, 1) for gc in gauss_full]
ok &= report("Gauss law vanishes at O(eps^0)",
             all(sp.expand(gc.subs(sub).doit()).coeff(e_, 0) == 0 for gc in gauss_full))
stamp("   linearised Gauss constraints built")

# ---------------------------------------------------------------- linear equations for u
stamp("Euler-Lagrange equations of L2")
eqs = [sp.expand(sp.euler_equations(L2, [uu], t)[0].lhs) for uu in ulist]
D4 = sp.Matrix(6, 6, lambda r, c: eqs[r].coeff(ulist[c].diff(t, 4)))
ok &= report("u'''' coefficients form -I / a^2", sp.simplify(D4 + sp.eye(6) / a ** 2) == sp.zeros(6, 6))

# ---------------------------------------------------------------- first-order form on the ghost-free background
stamp("first-order form on the ghost-free (pure Yang-Mills, conformal) background")
X = sp.symbols('x y z')
P = sp.symbols('px py pz')
AS = sp.symbols('a a1 a2 a3 a4')        # a and its first four eta-derivatives
bg = {}
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
U = [[sp.Symbol(f'u{k}_{j}') for j in range(6)] for k in range(5)]


def to_syms(expr):
    e = expr
    for k in range(4, -1, -1):
        e = e.subs({uu.diff(t, k) if k else uu: U[k][j] for j, uu in enumerate(ulist)})
    for k in range(4, -1, -1):
        e = e.subs({q[i].diff(t, k) if k else q[i]: bg[q[i].diff(t, k) if k else q[i]] for i in range(3)})
    for k in range(4, -1, -1):
        e = e.subs({a.diff(t, k) if k else a: AS[k]})
    return sp.expand(e)


eqs_s = [to_syms(e) for e in eqs]
# eq = -(1/a^2) u'''' + rest  ->  u'''' = a^2 (eq + u''''/a^2)
u4 = [sp.expand(AS[0] ** 2 * eqs_s[j] + U[4][j]) for j in range(6)]
state = [s for row in U[:4] for s in row]
ok &= report("equations are linear in the state", all(sp.Poly(u4[j], *state).total_degree() == 1 for j in range(6)))
ok &= report("u'''' does not appear on the right-hand side", all(not u4[j].has(*U[4]) for j in range(6)))
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
ok &= report("only a, a', a'' occur in K and G (no a''', a'''')",
             not K.has(AS[3], AS[4]) and not G.has(AS[3], AS[4]))
has2 = K.has(AS[2])
print(f"    (a'' occurs in K: {has2}; a' occurs in K: {K.has(AS[1])}; a' occurs in G: {G.has(AS[1])})")

# zero modes u = theta x A, constant theta (gauge symmetry is intact on FRW)
Zm = sp.zeros(24, 3)
for b in range(3):
    for j, (i, aa) in enumerate([(i, aa) for i in range(3) for aa in range(3) if i != aa]):
        val = eps(aa, b, i) * q[i]
        for k in range(4):
            Zm[k * 6 + j, b] = to_syms(val.diff(t, k)) if k else to_syms(val)
ok &= report("gauge zero modes satisfy the linearised Gauss law for arbitrary a(eta)", sp.simplify(G * Zm) == sp.zeros(3, 3))
Zdot = sp.zeros(24, 3)
for r in range(24):
    for c in range(3):
        Zdot[r, c] = sum(Zm[r, c].diff(X[i]) * P[i] + Zm[r, c].diff(P[i]) * qdd[i] for i in range(3))
ok &= report("u = theta x A solves the linear equations for arbitrary a(eta)",
             sp.simplify(sp.expand(K * Zm - Zdot)) == sp.zeros(24, 3))

# ---------------------------------------------------------------- a = 1 reproduces the committed system
stamp("comparing with the committed ym/rung1_system.pkl at a = 1, a' = a'' = 0")
ref = pickle.load(open(REF, 'rb'))
flat = {AS[0]: 1, AS[1]: 0, AS[2]: 0}
ren = dict(zip(list(X) + list(P) + [cj], list(ref['X']) + list(ref['P']) + [ref['cj']]))
Kflat = sp.expand(K.subs(flat).subs(ren))
Gflat = sp.expand(G.subs(flat).subs(ren))
ok &= report("K at a = 1 equals the committed K", sp.simplify(Kflat - ref['K']) == sp.zeros(24, 24))
ok &= report("Gauss matrix at a = 1 equals the committed G", sp.simplify(Gflat - ref['G']) == sp.zeros(3, 24))
stamp(f"ALL STRUCTURAL CHECKS {'PASSED' if ok else 'NOT ALL PASSED'}")

L2s = to_syms(L2)
with open(OUT, 'wb') as fh:
    pickle.dump(dict(K=K, G=G, Zm=Zm, L2=L2s, state=state, X=X, P=P, cj=cj, asyms=AS[:3],
                     ulist=[str(s) for s in ulist]), fh)
stamp(f"wrote {os.path.relpath(OUT)}")
print("\nterms of K's u'''' rows that carry a' or a'' (the non-adiabatic part), by row:")
for j in range(6):
    row = sp.expand(u4[j])
    na = sum(term for term in row.args if term.has(AS[1], AS[2]))
    print(f"  u''''_{j} ({ulist[j]}):  {sp.collect(na, [AS[1], AS[2]])}")
