"""
bianchi9_derive.py -- diagonal Bianchi IX minisuperspace of quadratic gravity

    S = int d^4x sqrt(-g) [ R/2 - C_{abcd}C^{abcd} / (4 M^2) ]      (M_P-bar = 1)

Mass normalisation: Salvio, "Quadratic gravity" (2018) writes
    L = (M_P^2/2) R - W^2/(2 f_2^2),   M_2^2 = f_2^2 M_P^2 / 2,
so with M_P = 1 and M == M_2 the Weyl coefficient is 1/(4 M^2).
Check 6 below gives the beta'' Hessian -6 e^{3 alpha}/M^2 against the GR kinetic term
+3 e^{3 alpha} beta'^2, i.e. the Pais-Uhlenbeck form x'^2/2 - x''^2/(2M^2) with frequencies
0 and M: the ghost pole is at M and non-tachyonic in this normalisation.

Metric:  ds^2 = -N^2 dt^2 + a^2 s1^2 + b^2 s2^2 + c^2 s3^2,
         d s^i = (1/2) eps_{ijk} s^j ^ s^k   (so a=b=c=1 is the round S^3 of radius 2,
         R_3 = 3/2; the 16 pi^2 volume is dropped).
Misner:  a = e^{alpha + beta_+ + sqrt3 beta_-}, b = e^{alpha + beta_+ - sqrt3 beta_-},
         c = e^{alpha - 2 beta_+}.

Method:  Cartan structure equations in the orthonormal frame e^0 = N dt, e^i = a_i s^i.
         Everything is a function of t only, so d acts as (1/N) d/dt on coefficients
         plus the structure functions. Riemann -> R, Ric^2, Riem^2 -> C^2, Gauss-Bonnet.

Checks performed (all must print OK):
  1. torsion-free connection solves de + w^e = 0
  2. Ricci symmetric
  3. static round case: R = 3/(2a^2), C^2 = 0; isotropic a=b=c: C^2 = 0 for any a(t)
  4. Gauss-Bonnet density N a b c G is a total derivative on the ansatz (E-L vanish)
  5. C^2 piece contains no alpha-ddot (conformal invariance) and no N-ddot
  6. Hessian in beta-ddot has rank 2
  7. the constraint (delta S / delta N) is conserved by the N=1 equations
  8. GR (vacuum Bianchi IX) solutions solve the full equations: the Einstein manifold
     is invariant  (Bach tensor vanishes on Einstein spaces)

Output: bianchi9_system.pkl with the N=1 equations of motion as sympy expressions,
plus a numerical RHS (cse + lambdify) used by bianchi9_test.py.

Phase space (N = 1):  y = (alpha, alpha', b+, b-, b+', b-', b+'', b-'', b+''', b-''')
10-d, one constraint -> 9-d.  Einstein manifold: b'' and b''' at their GR values,
6-d, constraint -> 5-d.  Transverse: 4-d = the two homogeneous traceless-diagonal
polarisations of the Weyl ghost, one oscillator each.

Run:  %run bianchi9_derive.py          (about a minute)
"""
import sys, time, pickle, itertools
import sympy as sp

t0 = time.time()
def stamp(msg):
    print(f"[{time.time()-t0:6.1f}s] {msg}", flush=True)

t = sp.symbols('t')
M = sp.symbols('M', positive=True)
N, a, b, c = [sp.Function(s)(t) for s in ('N', 'a', 'b', 'c')]
A = [None, a, b, c]
eta = [-1, 1, 1, 1]
idx = range(4)

def d0(f):                      # frame derivative of a function of t
    return sp.diff(f, t) / N

# ---------------------------------------------------------------- structure functions
# de^a = (1/2) F^a_{bc} e^b ^ e^c
F = [[[sp.S(0)]*4 for _ in idx] for _ in idx]
for i in (1, 2, 3):
    F[i][0][i] = sp.diff(A[i], t) / (N * A[i])
    F[i][i][0] = -F[i][0][i]
    for j in (1, 2, 3):
        for k in (1, 2, 3):
            e = sp.LeviCivita(i, j, k)
            if e != 0:
                F[i][j][k] = e * A[i] / (A[j] * A[k])
Fl = [[[eta[p] * F[p][q][r] for r in idx] for q in idx] for p in idx]   # F_{abc}

# ---------------------------------------------------------------- connection
# w_{abc}: w_{ab} = w_{abc} e^c,  w_{abc} = (F_{abc} + F_{bca} - F_{cab})/2
w = [[[sp.simplify((Fl[p][q][r] + Fl[q][r][p] - Fl[r][p][q]) / 2) for r in idx]
      for q in idx] for p in idx]
# check 1: torsion-free  F_{abc} = w_{abc} - w_{acb},  antisymmetry w_{abc} = -w_{bac}
ok = all(sp.simplify(Fl[p][q][r] - (w[p][q][r] - w[p][r][q])) == 0
         and sp.simplify(w[p][q][r] + w[q][p][r]) == 0
         for p in idx for q in idx for r in idx)
stamp(f"check 1 (torsion-free, antisymmetric connection): {'OK' if ok else 'FAIL'}")
wu = [[[eta[p] * w[p][q][r] for r in idx] for q in idx] for p in idx]     # w^a_{bc}

# ---------------------------------------------------------------- curvature
# R^a_b = d w^a_b + w^a_c ^ w^c_b = (1/2) R^a_{bde} e^d ^ e^e
Riem = [[[[sp.S(0)]*4 for _ in idx] for _ in idx] for _ in idx]
for p in idx:
    for q in idx:
        for dd in idx:
            for ee in idx:
                val = sp.S(0)
                # d(w^a_{bc}) ^ e^c  with d f = (f'/N) e^0
                if dd == 0:
                    val += d0(wu[p][q][ee])
                if ee == 0:
                    val -= d0(wu[p][q][dd])
                # w^a_{bc} de^c
                val += sum(wu[p][q][r] * F[r][dd][ee] for r in idx)
                # w^a_c ^ w^c_b
                val += sum(wu[p][r][dd] * wu[r][q][ee] - wu[p][r][ee] * wu[r][q][dd]
                           for r in idx)
                Riem[p][q][dd][ee] = sp.simplify(val)
stamp("Riemann done")

Ric = [[sp.simplify(sum(Riem[p][q][p][r] for p in idx)) for r in idx] for q in idx]
Rs = sp.simplify(sum(eta[q] * Ric[q][q] for q in idx))
ok = all(sp.simplify(Ric[p][q] - Ric[q][p]) == 0 for p in idx for q in idx)
stamp(f"check 2 (Ricci symmetric): {'OK' if ok else 'FAIL'}")

Riem_l = [[[[eta[p] * Riem[p][q][r][s] for s in idx] for r in idx] for q in idx] for p in idx]
Riem2 = sp.simplify(sum(eta[p]*eta[q]*eta[r]*eta[s] * Riem_l[p][q][r][s]**2
                        for p in idx for q in idx for r in idx for s in idx))
Ric2 = sp.simplify(sum(eta[p]*eta[q] * Ric[p][q]**2 for p in idx for q in idx))
C2 = sp.simplify(Riem2 - 2*Ric2 + Rs**2/3)          # Weyl^2 in 4-d
GB = sp.simplify(Riem2 - 4*Ric2 + Rs**2)             # Gauss-Bonnet
stamp("scalars done")

# check 3
A0 = sp.symbols('A0', positive=True)
static = {N: 1, a: A0, b: A0, c: A0}
R_static = sp.simplify(Rs.subs(static).doit())
C2_static = sp.simplify(C2.subs(static).doit())
f = sp.Function('f')(t)
C2_iso = sp.simplify(C2.subs({a: f, b: f, c: f}).doit())
ok = (sp.simplify(R_static - sp.Rational(3, 2)/A0**2) == 0 and C2_static == 0 and C2_iso == 0)
stamp(f"check 3 (round S^3: R = {R_static}, C^2 = {C2_static}; isotropic C^2 = {C2_iso}): "
      f"{'OK' if ok else 'FAIL'}")

# ---------------------------------------------------------------- Misner variables
al, bp, bm = [sp.Function(s)(t) for s in ('alpha', 'bp', 'bm')]
s3 = sp.sqrt(3)
misner = {a: sp.exp(al + bp + s3*bm), b: sp.exp(al + bp - s3*bm), c: sp.exp(al - 2*bp)}
sqrtg = N * a * b * c

L_R = (sqrtg * Rs / 2).subs(misner).doit()
L_C = (-sqrtg * C2 / (4 * M**2)).subs(misner).doit()
L_GB = (sqrtg * GB).subs(misner).doit()
stamp("Misner substitution done")

# ---------------------------------------------------------------- jet symbols
# q_k = k-th time derivative, as plain symbols, so we can differentiate w.r.t. them
funcs = [al, bp, bm, N]
names = ['al', 'bp', 'bm', 'N']
KMAX = 6
J = {nm: [sp.symbols(f'{nm}{k}') for k in range(KMAX + 1)] for nm in names}
def to_jet(expr):
    expr = expr.doit()
    rep = {}
    for fn, nm in zip(funcs, names):
        for k in range(KMAX, 0, -1):
            rep[sp.Derivative(fn, (t, k))] = J[nm][k]
    expr = expr.subs(rep)
    for fn, nm in zip(funcs, names):
        expr = expr.subs(fn, J[nm][0])
    return expr
def D(expr):                    # total time derivative on jet expressions
    return sum(sp.diff(expr, J[nm][k]) * J[nm][k+1]
               for nm in names for k in range(KMAX))

L_R, L_C, L_GB = [sp.expand(to_jet(e)) for e in (L_R, L_C, L_GB)]
stamp("jet form done")

def EL(L, nm, order=4):
    """Euler-Lagrange expression  sum_k (-D)^k dL/dq_k  for variable nm."""
    out = sp.S(0)
    for k in range(order + 1):
        term = sp.diff(L, J[nm][k])
        for _ in range(k):
            term = D(term)
        out += (-1)**k * term
    return sp.expand(out)

# check 4: Gauss-Bonnet is a total derivative on the ansatz
ok = all(sp.simplify(EL(L_GB, nm)) == 0 for nm in names)
stamp(f"check 4 (Gauss-Bonnet density is a total derivative): {'OK' if ok else 'FAIL'}")

# check 5: C^2 piece has no alpha'' and no N''
ok = sp.diff(L_C, J['al'][2]) == 0 and sp.diff(L_C, J['N'][2]) == 0
stamp(f"check 5 (C^2 term free of alpha'' and N''): {'OK' if ok else 'FAIL'}")

# ---------------------------------------------------------------- remove alpha'' etc. from the GR piece
# L_R is linear in the second derivatives with coefficients free of velocities; subtract
# d/dt(sum q' dL/dq'') to get the usual first-order Misner Lagrangian.
g = {nm: sp.diff(L_R, J[nm][2]) for nm in ('al', 'bp', 'bm')}
assert all(sp.diff(g[nm], J[m][2]) == 0 and sp.diff(g[nm], J[m][1]) == 0
           for nm in g for m in ('al', 'bp', 'bm')), "GR piece not linear in q''"
L_R1 = sp.expand(sp.simplify(L_R - D(sum(J[nm][1] * g[nm] for nm in g))))
assert all(sp.diff(L_R1, J[nm][2]) == 0 for nm in names)
stamp("GR piece reduced to first order:")
sp.pprint(sp.factor(L_R1))

# check 6: Hessian in (bp'', bm'') of the Weyl piece
Hs = sp.Matrix([[sp.diff(L_C, J[p][2], J[q][2]) for q in ('bp', 'bm')] for p in ('bp', 'bm')])
Hs = Hs.applyfunc(sp.simplify)
stamp(f"check 6 (Hessian rank {Hs.rank()}):  H =")
sp.pprint(Hs)

L = sp.expand(L_R1 + L_C)

# ---------------------------------------------------------------- equations, N = 1
gauge = {J['N'][0]: 1, J['N'][1]: 0, J['N'][2]: 0, J['N'][3]: 0}
E_al = EL(L, 'al').subs(gauge)
E_bp = EL(L, 'bp').subs(gauge)
E_bm = EL(L, 'bm').subs(gauge)
Con = EL(L, 'N').subs(gauge)                         # Hamiltonian constraint
stamp("Euler-Lagrange expressions done")

def orders(expr, nm):
    return [k for k in range(KMAX + 1) if sp.diff(expr, J[nm][k]) != 0]
for lab, ex in (('E_alpha', E_al), ('E_bp', E_bp), ('E_bm', E_bm), ('constraint', Con)):
    print(f"   {lab:11s} contains derivative orders  al:{orders(ex,'al')}  "
          f"bp:{orders(ex,'bp')}  bm:{orders(ex,'bm')}")

# alpha'' from E_alpha (linear).  Work with numerators and denominators explicitly:
# substituting rational expressions and expanding is what makes sympy crawl here.
c1 = E_al.coeff(J['al'][2])                     # = 6 (M^2 + bp1^2 + bm1^2) e^{3 al} / M^2
c0 = sp.expand(E_al - c1 * J['al'][2])
assert sp.diff(c0, J['al'][2]) == 0 and sp.diff(c1, J['al'][2]) == 0
al2 = -c0 / c1
print("   alpha'' coefficient:", sp.factor(c1))
# alpha''' = D(al2) = -(D c0 c1 - c0 D c1)/c1^2 with al2 re-substituted; cleared -> n3 / c1^3
n3 = sp.expand(-(D(c0) * c1 - c0 * D(c1)))
n3 = sp.expand(n3.subs(J['al'][2], -c0 / c1) * c1)
def clear(E):
    """E with al2 -> -c0/c1, al3 -> n3/c1^3, multiplied through by c1^3."""
    E = sp.expand(E)
    out = sp.S(0)
    for term in sp.Add.make_args(E):
        p2 = sp.degree(term, J['al'][2]); p3 = sp.degree(term, J['al'][3])
        assert p2 <= 1 and p3 <= 1 and p2 + p3 <= 1
        base = term / J['al'][2]**p2 / J['al'][3]**p3
        out += base * (-c0)**p2 * c1**(3 - p2 - 3*p3) * n3**p3
    return sp.expand(out)
Eb = [clear(E) for E in (E_bp, E_bm)]
Amat, rhs = sp.linear_eq_to_matrix(Eb, [J['bp'][4], J['bm'][4]])
assert all(sp.diff(Amat, J[nm][4]) == sp.zeros(2, 2) for nm in ('bp', 'bm'))
print("   beta-4 coefficient matrix / c1^3:", sp.factor(Amat[0, 0] / c1**3), "|", sp.factor(Amat[0, 1]))
b4 = list(Amat.LUsolve(rhs))
stamp("solved for alpha'' and beta-4")

# ---- numerical RHS: y = (al0, al1, bp0, bm0, bp1, bm1, bp2, bm2, bp3, bm3)
state = [J['al'][0], J['al'][1], J['bp'][0], J['bm'][0], J['bp'][1], J['bm'][1],
         J['bp'][2], J['bm'][2], J['bp'][3], J['bm'][3]]
import numpy as np
f_al2 = sp.lambdify(state + [M], al2, 'numpy', cse=True)
f_b4 = sp.lambdify(state + [M], b4, 'numpy', cse=True)
def rhs_num(y, Mv):
    a2 = f_al2(*y, Mv); q4 = f_b4(*y, Mv)
    return np.array([y[1], a2, y[4], y[5], y[6], y[7], y[8], y[9], q4[0], q4[1]])
f_Con = sp.lambdify(state + [M], Con.subs(J['al'][2], al2), 'numpy', cse=True)
e1 = np.eye(10)[1]

import random
random.seed(1)
def rand_point(Mv=1.0):
    y = np.array([random.uniform(-1, 1), 0, random.uniform(-.5, .5), random.uniform(-.5, .5)]
                 + [random.uniform(-1, 1) for _ in range(6)])
    for start in (0.5, -0.5, 1.5, -1.5, 3, -3):
        y[1] = start
        for _ in range(80):
            c = f_Con(*y, Mv)
            d = (f_Con(*(y + 1e-7*e1), Mv) - f_Con(*(y - 1e-7*e1), Mv)) / 2e-7
            if abs(d) < 1e-14: break
            y[1] -= c / d
            if abs(c) < 1e-13: break
        if abs(f_Con(*y, Mv)) < 1e-11:
            return y
    return None

# check 7: constraint conserved: d/dt Con along the flow, by central difference
ok = True; tested = 0
for _ in range(8):
    y = rand_point()
    if y is None: continue
    yd = rhs_num(y, 1.0)
    h = 1e-5 / (1 + np.abs(yd).max())
    dC = (f_Con(*(y + h*yd), 1.0) - f_Con(*(y - h*yd), 1.0)) / (2*h)
    scale = abs(f_Con(*(y + 0.3*yd), 1.0)) / 0.3 + 1e-30     # typical off-shell size
    ok &= abs(dC) < 1e-6 * scale; tested += 1
    print(f"   dCon/dt = {dC:.2e}  (off-shell scale {scale:.2e})")
stamp(f"check 7 (constraint conserved, {tested} random on-shell points): {'OK' if ok else 'FAIL'}")

# check 8: Einstein manifold invariant.  GR equations from L_R1 alone.
E_al_GR = EL(L_R1, 'al', 2).subs(gauge)
E_bp_GR = EL(L_R1, 'bp', 2).subs(gauge)
E_bm_GR = EL(L_R1, 'bm', 2).subs(gauge)
Con_GR = EL(L_R1, 'N', 2).subs(gauge)
gr2 = sp.solve([E_al_GR, E_bp_GR, E_bm_GR], [J['al'][2], J['bp'][2], J['bm'][2]], dict=True)
assert len(gr2) == 1
gr2 = {k: sp.expand(v) for k, v in gr2[0].items()}
gr3 = {J[nm][3]: sp.expand(D(gr2[J[nm][2]]).subs(gr2)) for nm in ('al', 'bp', 'bm')}
gr4 = {J[nm][4]: sp.expand(D(gr3[J[nm][3]]).subs(gr3).subs(gr2)) for nm in ('al', 'bp', 'bm')}
print("   GR:  alpha'' =", sp.factor(gr2[J['al'][2]]))
print("   GR:  bp''    =", sp.factor(gr2[J['bp'][2]]))
f_ConGR = sp.lambdify(state[:6] + [M], Con_GR, 'numpy', cse=True)
f_dConGR = sp.lambdify(state[:6] + [M], sp.diff(Con_GR, J['al'][1]), 'numpy', cse=True)
f_gr2 = sp.lambdify(state[:6] + [M], [gr2[J['bp'][2]], gr2[J['bm'][2]]], 'numpy', cse=True)
f_gr3 = sp.lambdify(state[:6] + [M], [gr3[J['bp'][3]], gr3[J['bm'][3]]], 'numpy', cse=True)
f_al2GR = sp.lambdify(state[:6] + [M], gr2[J['al'][2]], 'numpy', cse=True)
f_gr4 = sp.lambdify(state[:6] + [M], [gr4[J['bp'][4]], gr4[J['bm'][4]]], 'numpy', cse=True)
def gr_point():
    y = np.array([random.uniform(-1, 1), 0, random.uniform(-.5, .5), random.uniform(-.5, .5),
                  random.uniform(-1, 1), random.uniform(-1, 1)])
    for start in (0.5, -0.5, 1.5, -1.5):
        y[1] = start
        for _ in range(80):
            c = f_ConGR(*y, 1.0); d = f_dConGR(*y, 1.0)
            if abs(d) < 1e-14: break
            y[1] -= c / d
            if abs(c) < 1e-13: break
        if abs(f_ConGR(*y, 1.0)) < 1e-11:
            return np.concatenate([y, f_gr2(*y, 1.0), f_gr3(*y, 1.0)])
    return None
ok = True; tested = 0
for _ in range(8):
    y = gr_point()
    if y is None: continue
    yd = rhs_num(y, 1.0)
    d1 = abs(yd[1] - f_al2GR(*y[:6], 1.0))
    d2 = np.abs(yd[8:] - np.array(f_gr4(*y[:6], 1.0))).max() / (1 + np.abs(yd[8:]).max())
    d3 = abs(f_Con(*y, 1.0))
    ok &= d1 < 1e-9 and d2 < 1e-8 and d3 < 1e-9
    tested += 1
stamp(f"check 8 (Einstein manifold invariant: full flow = GR flow on GR data, {tested} pts): {'OK' if ok else 'FAIL'}")

# ---------------------------------------------------------------- GR values of b'', b''' (define the ghost-free manifold)
gr_b2 = [gr2[J['bp'][2]], gr2[J['bm'][2]]]
gr_b3 = [gr3[J[nm][3]] for nm in ('bp', 'bm')]

# ---------------------------------------------------------------- save
out = dict(state=[str(s) for s in state], M=str(M),
           L=L, L_R1=L_R1, L_C=L_C, Con=Con, al2=al2, b4=b4,
           Con_GR=Con_GR, gr_b2=gr_b2, gr_b3=gr_b3, al2_GR=gr2[J['al'][2]])
with open('bianchi9_system.pkl', 'wb') as fh:
    pickle.dump(out, fh)
stamp("wrote bianchi9_system.pkl")
