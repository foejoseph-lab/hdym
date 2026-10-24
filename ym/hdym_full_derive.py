#!/usr/bin/env python3
"""
hdym_full_derive.py -- the FULL nonlinear homogeneous SU(2) Lee-Wick system, all nine components
A^a_i (spatial i, colour a), temporal gauge A_0 = 0, units g = M = 1.  Symbolic derivation,
self-checks, and C code for the numerical stage (csrc/hdym_full_escape.c via islands/code/offdiag_escape.py).

The Lagrangian is built character for character as in hdym_rung1.py (which linearises it about the
diagonal background): L = L_YM + L_HD, L_HD = (c_j J_0^2 - sum_i J_i^2)/2, J^nu = eta^{mu mu} D_mu F^{mu nu},
A_0 kept until the Gauss law has been varied.  What is new here is that nothing is linearised:
the Euler-Lagrange equations for the nine A^a_i are solved for the fourth derivative,

    A''''_m = R_m(A, A', A'', A''';  c_j),      m = 3 i + a,

the A_0 equation gives the three Gauss constraints G^a(A, A', A'', A''') = 0 (conserved by the flow,
to be imposed on initial data), and the Ostrogradsky energy
    E = sum_m [ A'_m (dL/dA'_m - d/dt dL/dA''_m) + A''_m dL/dA''_m ] - L
is the conserved quantity used as a drift diagnostic.

State vector for the numerical stage: s = (A, A', A'', A''') each as 9 components, s[9 k + m].
Diagonal q_i = A^i_i = s[3 i + i]; off-diagonal u^a_i = s[3 i + a] (hdym_rung1's ulist order is
(i, a) = (0,1), (0,2), (1,0), (1,2), (2,0), (2,1) -> m = 1, 2, 3, 5, 6, 7).

Checks (CHECK ... PASS/FAIL), every one a structural fact that the earlier scripts also use:
  1  the coefficient of A'''' in the nine equations is -I (Lee-Wick kinetic sign, as hdym_rung1 D4 = -I)
  2  with every off-diagonal component zero the equations reduce to the Euler-Lagrange equations of the
     diagonal closed-form Lagrangian  1/2|q'|^2 - V - 1/2|q'' + grad V|^2  (hdym_rung1 "O(eps^0)")
  3  the Jacobian of R with respect to the 24 off-diagonal state entries, on a ghost-free diagonal
     background (q'' = -grad V, q''' = -Hess V q'), equals the K(q, p, c_j) of rung1_system.pkl at
     random points (numerical, 1e-9 relative)
  4  the Gauss constraints are conserved by the flow: dG/dt = 0 along R at random points (1e-9)
  5  the energy is conserved: dE/dt = 0 along R at random points (1e-9)
  6  G vanishes identically on the diagonal sector (the ansatz is Gauss-trivial)

Writes
  full_system.pkl        the sympy expressions (R, G, E, state symbols)
  csrc/hdym_full_rhs.h   generated C: hdym_full_rhs(s, ds, cj) (and hdym_full_rhs_acc for the fourth derivatives alone),
                         hdym_full_gauss(s, G, cj), hdym_full_energy_v(s, E, cj); used by csrc/hdym_full_escape.c
                         (islands/code/offdiag_escape.py, island_full.py)
Runtime: a few minutes (sympy).  Added 2026-10-19.
"""
import os
import pickle
import sys
import time

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()
OK = True


def stamp(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def check(name, ok):
    global OK
    OK &= bool(ok)
    print(f"CHECK {name}: {'PASS' if ok else 'FAIL'}", flush=True)


t = sp.Symbol('t', real=True)
eps = sp.LeviCivita
eta = sp.diag(1, -1, -1, -1)
cj = sp.Symbol('c_j')

# ---------------------------------------------------------------- the Lagrangian, as in hdym_rung1.py
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
LHD = (cj * sum(J[0][c] ** 2 for c in range(3)) - sum(J[nu][c] ** 2 for nu in range(1, 4) for c in range(3))) / 2
L = sp.expand(LYM + LHD)
stamp(f"   {len(L.args)} terms")

# ---------------------------------------------------------------- Gauss law (A_0 varied first), then temporal gauge
stamp("A_0 equations (Gauss law)")
gauss_full = [sp.euler_equations(L, [A[0][c]], t)[0].lhs for c in range(3)]
zero_A0 = {}
for c in range(3):
    for k in range(4, 0, -1):
        zero_A0[A[0][c].diff(t, k)] = 0
    zero_A0[A[0][c]] = 0
gauss_tg = [sp.expand(gc.subs(zero_A0).doit()) for gc in gauss_full]
Ltg = sp.expand(L.subs(zero_A0).doit())
stamp(f"   temporal-gauge Lagrangian: {len(Ltg.args)} terms")

# ---------------------------------------------------------------- Euler-Lagrange equations of the nine spatial components
stamp("Euler-Lagrange equations of the nine A^a_i")
comps = [A[i + 1][a] for i in range(3) for a in range(3)]            # m = 3 i + a
eqs = [sp.expand(sp.euler_equations(Ltg, [f], t)[0].lhs) for f in comps]
D4 = sp.Matrix(9, 9, lambda r, c: eqs[r].coeff(comps[c].diff(t, 4)))
check("1 A'''' coefficient matrix is -I", D4 == -sp.eye(9))
R_fun = [sp.expand(eqs[m] + comps[m].diff(t, 4)) for m in range(9)]  # A'''' = R

# ---------------------------------------------------------------- Ostrogradsky energy
stamp("Ostrogradsky energy")
Efun = -Ltg
for f in comps:
    dL_d1 = Ltg.diff(f.diff(t))
    dL_d2 = Ltg.diff(f.diff(t, 2))
    Efun += f.diff(t) * (dL_d1 - dL_d2.diff(t)) + f.diff(t, 2) * dL_d2
Efun = sp.expand(Efun)

# ---------------------------------------------------------------- to symbols: s[9k + m] = d^k A_m / dt^k
S = [[sp.Symbol(f's{k}_{m}') for m in range(9)] for k in range(5)]
state = [S[k][m] for k in range(4) for m in range(9)]          # 36


def to_syms(expr):
    e = expr
    for k in range(4, -1, -1):
        e = e.subs({(f.diff(t, k) if k else f): S[k][m] for m, f in enumerate(comps)})
    return sp.expand(e)


stamp("substituting symbols")
R = [to_syms(e) for e in R_fun]
G = [to_syms(e) for e in gauss_tg]
E = to_syms(Efun)
check("R, G, E contain only the 36 state symbols and c_j",
      all(ex.free_symbols <= set(state) | {cj} for ex in R + G + [E]))
stamp(f"   R: {[len(r.args) for r in R]} terms; G: {[len(g.args) for g in G]}; E: {len(E.args)}")

# ---------------------------------------------------------------- check 2: diagonal reduction
stamp("check 2: diagonal reduction")
X = sp.symbols('x y z')
Vs = (X[0] ** 2 * X[1] ** 2 + X[1] ** 2 * X[2] ** 2 + X[2] ** 2 * X[0] ** 2) / 2
qf = [sp.Function(n)(t) for n in 'xyz']
Vq = Vs.subs(dict(zip(X, qf)))
Lclosed = (sum(v.diff(t) ** 2 for v in qf) / 2 - Vq
           - sum((qf[i].diff(t, 2) + Vq.diff(qf[i])) ** 2 for i in range(3)) / 2)
eq_closed = [sp.expand(sp.euler_equations(Lclosed, [qf[i]], t)[0].lhs) for i in range(3)]
R_closed = [sp.expand(eq_closed[i] + qf[i].diff(t, 4)) for i in range(3)]
diag_sub = {S[k][m]: 0 for k in range(4) for m in range(9) if m % 4 != 0}       # off-diagonal m in {1,2,3,5,6,7}
Qs = [[sp.Symbol(f'q{k}_{i}') for i in range(3)] for k in range(4)]
diag_sub.update({S[k][4 * i]: Qs[k][i] for k in range(4) for i in range(3)})


def closed_to_syms(expr):
    e = expr
    for k in range(4, -1, -1):
        e = e.subs({(qf[i].diff(t, k) if k else qf[i]): Qs[k][i] for i in range(3)} if k < 4 else
                   {qf[i].diff(t, 4): sp.Symbol('unused') for i in range(3)})
    return sp.expand(e)


ok2 = True
for i in range(3):
    lhs = sp.expand(R[4 * i].subs(diag_sub))
    ok2 &= sp.expand(lhs - closed_to_syms(R_closed[i])) == 0
check("2 off-diagonals = 0 reduces to the diagonal closed-form equations", ok2)
check("6 Gauss law vanishes identically on the diagonal sector",
      all(sp.expand(g.subs(diag_sub)) == 0 for g in G))

# ---------------------------------------------------------------- numerical functions
stamp("lambdifying")
args = state + [cj]
R_f = sp.lambdify(args, R, 'numpy')
G_f = sp.lambdify(args, G, 'numpy')
E_f = sp.lambdify(args, E, 'numpy')
rng = np.random.default_rng(20261019)


def rhs36(s, c):
    r = np.array(R_f(*s, c), dtype=float)
    return np.concatenate([s[9:36], r])


# check 4, 5: conservation along the flow, by the chain rule with symbolic gradients
stamp("check 4/5: dG/dt and dE/dt along the flow at random points")
dG = [[sp.diff(g, v) for v in state] for g in G]
dE = [sp.diff(E, v) for v in state]
dG_f = sp.lambdify(args, dG, 'numpy')
dE_f = sp.lambdify(args, dE, 'numpy')
worstG = worstE = 0.0
for trial in range(6):
    c = 1.0 if trial < 4 else 0.0
    s = rng.normal(size=36) * np.array([1.0] * 9 + [1.0] * 9 + [2.0] * 9 + [4.0] * 9)
    ds = rhs36(s, c)
    gdot = np.array(dG_f(*s, c), dtype=float) @ ds
    edot = float(np.array(dE_f(*s, c), dtype=float) @ ds)
    scale_g = np.abs(np.array(dG_f(*s, c), dtype=float)).sum() * np.abs(ds).max()
    scale_e = np.abs(np.array(dE_f(*s, c), dtype=float)).sum() * np.abs(ds).max()
    worstG = max(worstG, np.abs(gdot).max() / scale_g)
    worstE = max(worstE, abs(edot) / scale_e)
check(f"4 Gauss constraints conserved by the flow (worst relative {worstG:.1e}, tol 1e-9)", worstG < 1e-9)
check(f"5 energy conserved by the flow (worst relative {worstE:.1e}, tol 1e-9)", worstE < 1e-9)

# check 3: linearisation about a ghost-free diagonal background = rung1 K
stamp("check 3: Jacobian on the diagonal background vs rung1_system.pkl")
pkl = os.path.join(HERE, 'rung1_system.pkl')
if os.path.exists(pkl):
    dd = pickle.load(open(pkl, 'rb'))
    Kf = sp.lambdify(list(dd['X']) + list(dd['P']) + [dd['cj']], dd['K'], 'numpy')
    od = [1, 2, 3, 5, 6, 7]                                            # m of ulist order
    od_state = [9 * k + m for k in range(4) for m in od]               # 24 entries, rung1 order
    JR = [[sp.diff(R[m], state[c]) for c in od_state] for m in od]     # 6 x 24
    JR_f = sp.lambdify(args, JR, 'numpy')
    gV = [Vs.diff(v) for v in X]
    HV = sp.hessian(Vs, X)
    gV_f = sp.lambdify(X, gV, 'numpy')
    HV_f = sp.lambdify(X, HV, 'numpy')
    worst = 0.0
    for trial in range(5):
        q = rng.normal(size=3) * 1.5
        p = rng.normal(size=3)
        qdd = -np.array(gV_f(*q), dtype=float)
        qddd = -np.array(HV_f(*q), dtype=float) @ p
        s = np.zeros(36)
        for i in range(3):
            s[4 * i] = q[i]
            s[9 + 4 * i] = p[i]
            s[18 + 4 * i] = qdd[i]
            s[27 + 4 * i] = qddd[i]
        for c in (1.0, 0.0):
            Jnum = np.array(JR_f(*s, c), dtype=float)
            K = np.array(Kf(*q, *p, c), dtype=float)[18:24, :]
            worst = max(worst, np.abs(Jnum - K).max() / max(np.abs(K).max(), 1.0))
    check(f"3 Jacobian on the ghost-free diagonal background = rung1 K (worst relative {worst:.1e}, tol 1e-9)", worst < 1e-9)
else:
    print("CHECK 3 skipped: rung1_system.pkl not found next to this script (run hdym_rung1.py)")
    OK = False

# ---------------------------------------------------------------- write the pickle and the C code
stamp("writing full_system.pkl")
with open(os.path.join(HERE, 'full_system.pkl'), 'wb') as fh:
    pickle.dump(dict(R=R, G=G, E=E, state=state, cj=cj), fh)

stamp("generating C (cse)")
from sympy.printing.c import C99CodePrinter          # noqa: E402

printer = C99CodePrinter()


def _nopow(code):
    """pow(x, n) with small integer n -> explicit products (pow is ~50x a multiply and the kernel is all small powers)."""
    import re
    def rep(m):
        x, n = m.group(1), int(m.group(2))
        return "(" + "*".join([x] * n) + ")"
    return re.sub(r"pow\((\w+), ([2-6])\)", rep, code)


def cfun(name, exprs, outname, nout):
    reps, red = sp.cse(exprs, symbols=sp.numbered_symbols('w'), order='none')
    lines = [f"DEVFN void {name}(const double *s, double *{outname}, double c_j)", "{"]
    for k in range(4):
        for m in range(9):
            lines.append(f"    const double s{k}_{m} = s[{9 * k + m}];")
    for sym, ex in reps:
        lines.append(f"    const double {sym} = {_nopow(printer.doprint(ex))};")
    for j, ex in enumerate(red):
        lines.append(f"    {outname}[{j}] = {_nopow(printer.doprint(ex))};")
    lines.append("}")
    return "\n".join(lines), len(reps)


c_rhs, n1 = cfun('hdym_full_rhs_acc', R, 'r', 9)
c_gauss, n2 = cfun('hdym_full_gauss', G, 'G', 3)
c_energy, n3 = cfun('hdym_full_energy_v', [E], 'E', 1)
header = ("/* hdym_full_rhs.h -- GENERATED by ym/hdym_full_derive.py on " + time.strftime("%Y-%m-%d") + "; do not edit.\n"
          " *\n"
          " * The full nonlinear homogeneous SU(2) Lee-Wick system, nine components A^a_i, temporal gauge,\n"
          " * units g = M = 1, as a first-order system in s = (A, A1, A2, A3) = (A and its first three time\n"
          " * derivatives), 36 doubles, s[9k + m], m = 3i + a:\n"
          " *   hdym_full_rhs_acc(s, r, c_j)   r[m] = the fourth derivative of A_m (the other 27 components of ds/dt are s[9..35])\n"
          " *   hdym_full_gauss(s, G, c_j)     the three Gauss constraints (conserved; zero on physical data)\n"
          " *   hdym_full_energy_v(s, E, c_j)  E[0] = the Ostrogradsky energy (conserved)\n"
          " * c_j multiplies the J_0^2 piece of the Lee-Wick term (1 = the theory).  Temporaries after cse:\n"
          f" * rhs {n1}, gauss {n2}, energy {n3}.  DEVFN as in hdym_model.h.\n"
          " */\n"
          "#ifndef HDYM_FULL_RHS_H\n#define HDYM_FULL_RHS_H\n#ifndef DEVFN\n#define DEVFN static inline\n#endif\n\n"
          + c_rhs + "\n\n" + c_gauss + "\n\n" + c_energy + "\n\n"
          "DEVFN void hdym_full_rhs(const double *s, double *ds, double c_j)\n{\n"
          "    for (int j = 0; j < 27; ++j) ds[j] = s[j + 9];\n"
          "    hdym_full_rhs_acc(s, ds + 27, c_j);\n}\n#endif\n")
os.makedirs(os.path.join(HERE, 'csrc'), exist_ok=True)
with open(os.path.join(HERE, 'csrc', 'hdym_full_rhs.h'), 'w') as fh:
    fh.write(header)
stamp(f"wrote csrc/hdym_full_rhs.h ({n1}/{n2}/{n3} cse temporaries)")

# the generated C must agree with the sympy functions: compile a tiny test if gcc is available
try:
    import ctypes
    import subprocess
    import tempfile
    sys.path.insert(0, HERE)
    import cbuild
    dtmp = tempfile.mkdtemp()
    src = os.path.join(dtmp, 't.c')
    lib = os.path.join(dtmp, 't' + ('.dll' if sys.platform == 'win32' else '.so'))
    open(src, 'w').write('#include <math.h>\n#include "%s"\nvoid rhs(const double*s,double*ds,double c){hdym_full_rhs(s,ds,c);}\n'
                         'void gauss(const double*s,double*G,double c){hdym_full_gauss(s,G,c);}\n'
                         'double energy(const double*s,double c){double E;hdym_full_energy_v(s,&E,c);return E;}\n'
                         % os.path.join(HERE, 'csrc', 'hdym_full_rhs.h').replace('\\', '/'))
    cbuild.check_call([os.environ.get('CC', 'gcc'), '-O2', '-shared', '-fPIC', '-o', lib, src, '-lm'])
    Lb = ctypes.CDLL(lib)
    dp = ctypes.POINTER(ctypes.c_double)
    Lb.rhs.argtypes = [dp, dp, ctypes.c_double]
    Lb.gauss.argtypes = [dp, dp, ctypes.c_double]
    Lb.energy.argtypes = [dp, ctypes.c_double]
    Lb.energy.restype = ctypes.c_double
    worst = 0.0
    for trial in range(4):
        s = rng.normal(size=36)
        ds = np.zeros(36)
        Gc = np.zeros(3)
        Lb.rhs(s.ctypes.data_as(dp), ds.ctypes.data_as(dp), 1.0)
        Lb.gauss(s.ctypes.data_as(dp), Gc.ctypes.data_as(dp), 1.0)
        Ec = Lb.energy(s.ctypes.data_as(dp), 1.0)
        worst = max(worst, np.abs(ds - rhs36(s, 1.0)).max() / np.abs(ds).max(),
                    np.abs(Gc - np.array(G_f(*s, 1.0), dtype=float)).max() / max(np.abs(Gc).max(), 1.0),
                    abs(Ec - float(E_f(*s, 1.0))) / max(abs(Ec), 1.0))
    check(f"7 generated C agrees with sympy (worst relative {worst:.1e}, tol 1e-12)", worst < 1e-12)
except Exception as ex:                                           # noqa: BLE001
    print(f"CHECK 7 skipped (no compiler here): {ex}")

print(f"\n{'ALL CHECKS PASS' if OK else 'SOME CHECK FAILED'}")
print(f"DONE {time.time() - T0:.0f}")
