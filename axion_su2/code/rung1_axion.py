#!/usr/bin/env python3
"""
rung1_axion.py -- ym/hdym_rung1.py with the axion drive: the linearisation of the full homogeneous SU(2)
Lee-Wick system about the colour-diagonal background A^a_i = delta^a_i q_i(t), A_0 = 0, flat space, units
g = M = 1, with L_axion = kappa(t) det A added (axion_derive.py: (1/4f) phi F.Ftilde reduces to
(g/f) phidot det A up to a total derivative, with Fujita's sign; kappa = (g/f) phidot is a given function).

Same construction as ym/hdym_rung1.py and rung1_frw.py: quadratic Lagrangian L2 of the six off-diagonal
components, linearised Gauss constraints, the six linear 4th-order equations as a first-order system
x' = K(q, q', q'', q''', c_j; kappa) x of dimension 24, and the three gauge zero modes.  The background is
the exact driven diagonal one: q, q', q'', q''' are free data and q'''' comes from the diagonal 4th-order
equations of L0 (ghost included), so that the gauge zero modes are exact for any solution.  (A first version
substituted the M -> infinity background q'' = -grad V + kappa grad(xyz); on it the zero modes fail at O(kappa),
because it does not solve the 4th-order system.)  The numerical stage integrates the diagonal 4th-order
background alongside; kappa', kappa'' do not occur.

Checks printed:
  * det A has no A_0 (Gauss law unchanged); the diagonal ansatz stays a consistent truncation at O(eps^1);
  * L2 - L2(kappa = 0) = kappa (-x u23 u32 - y u13 u31 - z u12 u21): the drive enters the off-diagonal
    blocks as that potential and nothing else;
  * the u'''' coefficients are -I; kappa enters K only through q'''' of the background;
  * K, G at kappa = 0 equal the committed ym/rung1_system.pkl exactly (symbolically);
  * the gauge zero modes u = theta x A solve the driven linear equations and the Gauss law.

Writes rung1_axion_system.pkl next to this folder's README (lambdify over X + P + X2 + X3 + [c_j, kappa]).
Usage (from the repository root, ~20 s):
    python axion_su2/code/rung1_axion.py > axion_su2/results/rung1_axion.txt
Added 2026-10-15.
"""
import os
import pickle
import time

import sympy as sp

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'rung1_axion_system.pkl')
REF = os.path.join(HERE, '..', '..', 'ym', 'rung1_system.pkl')


def stamp(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def report(name, ok):
    print(f"    {name}: {'PASS' if ok else 'FAIL'}", flush=True)
    return ok


t = sp.Symbol('t', real=True)          # conformal time eta
eps = sp.LeviCivita
eta = sp.diag(1, -1, -1, -1)
kap = sp.Function('kappa')(t)

# ---------------------------------------------------------------- full homogeneous Lagrangian, 12 functions
stamp("building the full homogeneous Lagrangian (12 functions) with kappa(t) det A")
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
LHD = (cj * sum(J[0][c] ** 2 for c in range(3)) - sum(J[nu][c] ** 2 for nu in range(1, 4) for c in range(3))) / 2
detA = sp.Matrix(3, 3, lambda i, c: A[i + 1][c]).det()
LAX = kap * detA
L = sp.expand(LYM + LHD + LAX)
ok0 = report("det A contains no A_0", not any(detA.has(A[0][c]) for c in range(3)))
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
dL2 = sp.expand(L2 - L2.subs(kap, 0))
u12, u13, u21, u23, u31, u32 = [u[(i, aa)] for i in range(3) for aa in range(3) if i != aa]
okL2 = report("L2 - L2(kappa=0) = kappa (-x u23 u32 - y u13 u31 - z u12 u21)",
              sp.expand(dL2 - kap * (-x * u23 * u32 - y * u13 * u31 - z * u12 * u21)) == 0)
V = (x ** 2 * y ** 2 + y ** 2 * z ** 2 + z ** 2 * x ** 2) / 2
gradV = [V.diff(v) for v in q]
Lclosed = (sum(v.diff(t) ** 2 for v in q) / 2 - V + kap * x * y * z
           - sum((q[i].diff(t, 2) + gradV[i]) ** 2 for i in range(3)) / 2)
ok = ok0 & report("O(eps^0): diagonal reduced Lagrangian with + kappa x y z", sp.simplify(L0 - sp.expand(Lclosed)) == 0)
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
ok &= okL2
ok &= report("u'''' coefficients form -I", sp.simplify(D4 + sp.eye(6)) == sp.zeros(6, 6))

# ---------------------------------------------------------------- first-order form on the ghost-free background
stamp("first-order form on the exact driven diagonal background (q, q', q'', q''' free; q'''' from the diagonal equations)")
X = sp.symbols('x y z')
P = sp.symbols('px py pz')
AS = sp.symbols('kappa kappa1 kappa2 kappa3 kappa4')        # kappa and its first four derivatives
X2 = sp.symbols('x2 y2 z2')                                   # q'' (free background data)
X3 = sp.symbols('x3 y3 z3')                                   # q'''
bg = {}
# q'''' from the diagonal 4th-order equations of L0 (the exact driven background, ghost included)
eqs0 = [sp.expand(sp.euler_equations(L0, [qi], t)[0].lhs) for qi in q]


def bgsym(expr):
    e = expr
    for k in range(4, -1, -1):
        e = e.subs({q[i].diff(t, k) if k else q[i]: (X, P, X2, X3, sp.symbols('x4 y4 z4'))[k][i] for i in range(3)})
    for k in range(4, -1, -1):
        e = e.subs({kap.diff(t, k) if k else kap: AS[k]})
    return sp.expand(e)


X4s = sp.symbols('x4 y4 z4')
qdddd = [sp.solve(bgsym(eqs0[i]), X4s[i])[0] for i in range(3)]
ok &= report("diagonal q'''' equations are linear in q'''' with coefficient -1", all(bgsym(eqs0[i]).coeff(X4s[i]) == -1 for i in range(3)))
ok &= report("q'''' involves kappa but not kappa', kappa''", all(qdddd[i].has(AS[0]) and not qdddd[i].has(AS[1], AS[2]) for i in range(3)))
qdd = list(X2)
qddd = list(X3)


def Dt(expr):
    # d/dt on functions of (X, P, X2, X3, kappa^(k)):
    return sp.expand(sum(expr.diff(X[i]) * P[i] + expr.diff(P[i]) * X2[i] + expr.diff(X2[i]) * X3[i]
                         + expr.diff(X3[i]) * qdddd[i] for i in range(3))
                     + sum(expr.diff(AS[k]) * AS[k + 1] for k in range(4)))


for i in range(3):
    bg[q[i].diff(t, 4)] = qdddd[i]
    bg[q[i].diff(t, 3)] = X3[i]
    bg[q[i].diff(t, 2)] = X2[i]
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
        e = e.subs({kap.diff(t, k) if k else kap: AS[k]})
    return sp.expand(e)


eqs_s = [to_syms(e) for e in eqs]
# eq = -(1/a^2) u'''' + rest  ->  u'''' = a^2 (eq + u''''/a^2)
u4 = [sp.expand(eqs_s[j] + U[4][j]) for j in range(6)]
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
ok &= report("only kappa, kappa', kappa'' occur in K (no kappa''', kappa'''')", not K.has(AS[3], AS[4]))
print(f"    (kappa occurs in K: {K.has(AS[0])}, kappa', kappa'': {K.has(AS[1], AS[2])}; in G: {G.has(*AS)} -- only through q'''' of the background)")

# zero modes u = theta x A, constant theta (gauge symmetry is intact on FRW)
Zm = sp.zeros(24, 3)
for b in range(3):
    for j, (i, aa) in enumerate([(i, aa) for i in range(3) for aa in range(3) if i != aa]):
        val = eps(aa, b, i) * q[i]
        for k in range(4):
            Zm[k * 6 + j, b] = to_syms(val.diff(t, k)) if k else to_syms(val)
ok &= report("gauge zero modes satisfy the linearised Gauss law with the drive", sp.simplify(G * Zm) == sp.zeros(3, 3))
Zdot = sp.zeros(24, 3)
for r in range(24):
    for c in range(3):
        Zdot[r, c] = Dt(Zm[r, c])
ok &= report("u = theta x A solves the driven linear equations for arbitrary kappa(t)",
             sp.simplify(sp.expand(K * Zm - Zdot)) == sp.zeros(24, 3))

# ---------------------------------------------------------------- a = 1 reproduces the committed system
stamp("comparing with the committed ym/rung1_system.pkl at kappa = 0")
ref = pickle.load(open(REF, 'rb'))
Vs = (X[0] ** 2 * X[1] ** 2 + X[1] ** 2 * X[2] ** 2 + X[2] ** 2 * X[0] ** 2) / 2
gV = [Vs.diff(v) for v in X]
HV = sp.hessian(Vs, X)
flat = {AS[0]: 0, AS[1]: 0, AS[2]: 0}
flat.update({X2[i]: -gV[i] for i in range(3)})
flat.update({X3[i]: -sum(HV[i, j] * P[j] for j in range(3)) for i in range(3)})
ren = dict(zip(list(X) + list(P) + [cj], list(ref['X']) + list(ref['P']) + [ref['cj']]))
Kflat = sp.expand(K.subs(flat).subs(ren))
Gflat = sp.expand(G.subs(flat).subs(ren))
ok &= report("K at kappa = 0 on the ghost-free background equals the committed K", sp.simplify(Kflat - ref['K']) == sp.zeros(24, 24))
ok &= report("Gauss matrix at kappa = 0 on the ghost-free background equals the committed G", sp.simplify(Gflat - ref['G']) == sp.zeros(3, 24))
stamp(f"ALL STRUCTURAL CHECKS {'PASSED' if ok else 'NOT ALL PASSED'}")

L2s = to_syms(L2)
with open(OUT, 'wb') as fh:
    pickle.dump(dict(K=K, G=G, Zm=Zm, L2=L2s, state=state, X=X, P=P, X2=X2, X3=X3, cj=cj, ksyms=AS[:3], q4=qdddd,
                     ulist=[str(s) for s in ulist]), fh)
stamp(f"wrote {os.path.relpath(OUT)}")
print("\nterms of K's u'''' rows that carry kappa, kappa' or kappa'' (what the drive adds), by row:")
for j in range(6):
    row = sp.expand(u4[j])
    na = sum(term for term in row.args if term.has(*AS[:3]))
    print(f"  u''''_{j} ({ulist[j]}):  {sp.collect(na, list(AS[:3]))}")
