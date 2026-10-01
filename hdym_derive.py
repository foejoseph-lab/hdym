#!/usr/bin/env python3
"""
hdym_derive.py -- Step 1: symbolic derivation and consistency checks for the
spatially homogeneous sector of SU(2) Lee-Wick (higher-derivative) Yang-Mills.

Theory (metric eta = diag(+1,-1,-1,-1)):

    L = -1/4 F^a_{mu nu} F^{a mu nu} + 1/(2 M^2) (D^mu F_{mu nu})^a (D_rho F^{rho nu})^a

The sign of the M^{-2} term is the one that gives a massive (non-tachyonic)
spin-1 ghost at p^2 = M^2 in this metric signature; CHECK 0 below verifies it
from the linearised equations instead of trusting the comment.

What this script does (each CHECK prints PASS/FAIL):

  CHECK 0  Sign convention: linearised homogeneous equation has ghost frequency M
           (oscillatory), not a tachyon.
  CHECK 1  Full homogeneous reduction A^a_mu(t) (12 functions incl. A_0) ->
           Euler-Lagrange equations (4th order in t) + Gauss-law constraint.
  CHECK 2  Diagonal ansatz A^a_i = delta^a_i q_i(t), A_0 = 0 is a CONSISTENT
           TRUNCATION: all off-diagonal EOMs and the Gauss law vanish identically.
  CHECK 3  Reduced Lagrangian equals the closed form
               L_red = 1/2|qdot|^2 - V(q) - 1/(2M^2) |qddot + grad V|^2,
               V = g^2/2 (x^2 y^2 + y^2 z^2 + z^2 x^2)
           and its EOMs coincide with the restricted full EOMs
           (symmetric criticality holds here, verified rather than assumed).
  CHECK 4  Auxiliary-field (order-reduced) Lagrangian
               L' = 1/2|qdot|^2 + wdot.qdot - V - w.grad V + M^2/2 |w|^2
           reproduces the same 4th-order dynamics.  With p = q + w:
               L' = 1/2|pdot|^2 - 1/2|wdot|^2 - U(p,w)
           so p is the normal sector and w the ghost sector.
  CHECK 5  Scaling: t -> tau/M, q -> (M/g) u gives L = (M^4/g^2) L~(u) with
           g = M = 1.  The classical problem has NO free couplings; only the
           initial data (energy) matters.

Output: prints results; writes `reduced_system.txt` with the Hamiltonian used by
the numerics (in units g = M = 1).

Runtime is measured and printed per check.  Expected: CHECK 1-2 dominate.
"""
import time
import sympy as sp

T0 = time.time()


def stamp(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def report(name, ok):
    print(f"    {name}: {'PASS' if ok else 'FAIL'}", flush=True)
    return ok


t = sp.Symbol('t', real=True)
g, M = sp.symbols('g M', positive=True)
eta = sp.diag(1, -1, -1, -1)
eps = sp.LeviCivita

# ---------------------------------------------------------------- CHECK 0
stamp("CHECK 0: sign convention via abelian homogeneous limit")
a = sp.Function('a')(t)
# abelian, homogeneous, A_0 = 0, A_1 = a(t):  F_{01} = a', (d^mu F_{mu 1}) = a''
# -1/4 F^2 = +1/2 a'^2 ;  (D^mu F_{mu nu})(D_rho F^{rho nu}) = -(a'')^2
L0 = sp.Rational(1, 2) * a.diff(t)**2 + 1 / (2 * M**2) * (-(a.diff(t, 2))**2)
eom0 = sp.euler_equations(L0, [a], t)[0].lhs
# try a'' ~ exp(i w t)
w_ = sp.Symbol('omega')
trial = eom0.subs(a, sp.exp(sp.I * w_ * t)).doit()
roots = sp.solve(sp.simplify(trial / sp.exp(sp.I * w_ * t)), w_)
stamp(f"   linear frequencies: {roots}")
ok0 = report("ghost frequency is +-M (oscillatory)", set(roots) >= {M, -M})

# ---------------------------------------------------------------- CHECK 1
stamp("CHECK 1: full homogeneous reduction (12 functions)")
A = [[sp.Function(f'A{mu}{c}')(t) for c in range(3)] for mu in range(4)]


def d(mu, expr):
    return expr.diff(t) if mu == 0 else sp.Integer(0)


# F_{mu nu}^a
F = [[[None] * 3 for _ in range(4)] for _ in range(4)]
for mu in range(4):
    for nu in range(4):
        for c in range(3):
            F[mu][nu][c] = (d(mu, A[nu][c]) - d(nu, A[mu][c])
                            + g * sum(eps(c, b, e) * A[mu][b] * A[nu][e]
                                      for b in range(3) for e in range(3)))


def Dcov(mu, X):
    """adjoint covariant derivative of a colour vector X[c]"""
    return [d(mu, X[c]) + g * sum(eps(c, b, e) * A[mu][b] * X[e]
                                  for b in range(3) for e in range(3))
            for c in range(3)]


# J_nu^a = D^mu F_{mu nu}^a  (index raised with eta)
J = []
for nu in range(4):
    acc = [sp.Integer(0)] * 3
    for mu in range(4):
        col = [F[mu][nu][c] for c in range(3)]
        Dm = Dcov(mu, col)
        acc = [acc[c] + eta[mu, mu] * Dm[c] for c in range(3)]
    J.append(acc)

LYM = sp.Integer(0)
for mu in range(4):
    for nu in range(4):
        for c in range(3):
            LYM += -sp.Rational(1, 4) * eta[mu, mu] * eta[nu, nu] * F[mu][nu][c]**2
LHD = sum(eta[nu, nu] * J[nu][c]**2 for nu in range(4) for c in range(3)) / (2 * M**2)
L = sp.expand(LYM + LHD)
stamp(f"   Lagrangian built: {len(L.args)} terms")

funcs = [A[mu][c] for mu in range(4) for c in range(3)]
eqs = sp.euler_equations(L, funcs, t)
stamp("   Euler-Lagrange equations derived")
EOM = {(mu, c): eqs[3 * mu + c].lhs for mu in range(4) for c in range(3)}
gauss = [EOM[(0, c)] for c in range(3)]
ok1 = report("12 EOMs obtained", len(eqs) == 12)

# ---------------------------------------------------------------- CHECK 2
stamp("CHECK 2: diagonal ansatz is a consistent truncation")
x, y, z = [sp.Function(n)(t) for n in 'xyz']
q = [x, y, z]
ans = {}
for c in range(3):
    ans[A[0][c]] = sp.Integer(0)
    for i in range(3):
        ans[A[i + 1][c]] = q[i] if i == c else sp.Integer(0)


def on_ansatz(expr):
    return sp.expand(expr.subs(ans).doit())


ok2 = True
for c in range(3):
    ok2 &= report(f"Gauss law component a={c + 1} vanishes", on_ansatz(gauss[c]) == 0)
for i in range(3):
    for c in range(3):
        if i != c:
            ok2 &= report(f"off-diagonal EOM (i={i + 1},a={c + 1}) vanishes",
                          on_ansatz(EOM[(i + 1, c)]) == 0)
diag_eoms = [on_ansatz(EOM[(i + 1, i)]) for i in range(3)]
stamp("   restricted diagonal EOMs computed")

# ---------------------------------------------------------------- CHECK 3
stamp("CHECK 3: closed-form reduced Lagrangian and symmetric criticality")
V = g**2 / 2 * (x**2 * y**2 + y**2 * z**2 + z**2 * x**2)
gradV = [V.diff(v) for v in q]
Lclosed = (sp.Rational(1, 2) * sum(v.diff(t)**2 for v in q) - V
           - sum((q[i].diff(t, 2) + gradV[i])**2 for i in range(3)) / (2 * M**2))
Lred = on_ansatz(L)
ok3 = report("L|ansatz == closed form", sp.simplify(Lred - sp.expand(Lclosed)) == 0)
red_eqs = sp.euler_equations(Lclosed, q, t)
for i in range(3):
    # the full EOM for A_i^i is dS/dA_i^i; reduced EOM is dS/dq_i: same sign convention
    diff_ = sp.simplify(sp.expand(red_eqs[i].lhs) - diag_eoms[i])
    ok3 &= report(f"reduced EOM {q[i].func} == restricted full EOM", diff_ == 0)

# ---------------------------------------------------------------- CHECK 4
stamp("CHECK 4: auxiliary-field order reduction")
wv = [sp.Function(f'w{n}')(t) for n in 'xyz']
Lp = (sp.Rational(1, 2) * sum(v.diff(t)**2 for v in q)
      + sum(wv[i].diff(t) * q[i].diff(t) for i in range(3))
      - V - sum(wv[i] * gradV[i] for i in range(3))
      + M**2 / 2 * sum(wv[i]**2 for i in range(3)))
aux_eqs = sp.euler_equations(Lp, q + wv, t)
# w-equations are algebraic in w: solve, substitute into q-equations
wsol = sp.solve([e.lhs for e in aux_eqs[3:]], wv, dict=True)[0]
ok4 = True
for i in range(3):
    lhs = sp.expand(aux_eqs[i].lhs.subs(wsol).doit())
    # compare with reduced 4th-order EOM up to overall factor
    ref = sp.expand(red_eqs[i].lhs)
    ratio = sp.simplify(lhs / ref)
    ok4 &= report(f"aux q-equation {i} proportional to reduced EOM (ratio={ratio})",
                  ratio.free_symbols <= {g, M} and not ratio.has(t))

# ---------------------------------------------------------------- CHECK 5
stamp("CHECK 5: scaling to g = M = 1")
tau = sp.Symbol('tau', real=True)
u = [sp.Function(f'u{n}')(tau) for n in 'xyz']
# q_i(t) = (M/g) u_i(M t);  d/dt = M d/dtau
sub = {}
for i in range(3):
    qi = q[i]
    sub[qi.diff(t, 2)] = (M / g) * M**2 * u[i].diff(tau, 2)
    sub[qi.diff(t)] = (M / g) * M * u[i].diff(tau)
    sub[qi] = (M / g) * u[i]
Ls = sp.expand(Lclosed.subs(sub))
Vu = (u[0]**2 * u[1]**2 + u[1]**2 * u[2]**2 + u[2]**2 * u[0]**2) / 2
gVu = [Vu.diff(v) for v in u]
Lu = (sp.Rational(1, 2) * sum(v.diff(tau)**2 for v in u) - Vu
      - sum((u[i].diff(tau, 2) + gVu[i])**2 for i in range(3)) / 2)
ok5 = report("L = (M^4/g^2) * L~(u; g=M=1)",
             sp.simplify(Ls - M**4 / g**2 * sp.expand(Lu)) == 0)

# ---------------------------------------------------------------- summary
allok = all([ok0, ok1, ok2, ok3, ok4, ok5])
stamp(f"ALL CHECKS {'PASSED' if allok else 'NOT ALL PASSED'}")

txt = r"""
Reduced system used by the numerics (units g = M = 1, time in units of 1/M,
fields in units of M/g):

  q = (x, y, z) diagonal colour-space amplitudes, V(q) = 1/2 (x^2 y^2 + y^2 z^2 + z^2 x^2)
  normal coordinate  p = q + w,   ghost coordinate  w   (so q = p - w)

  H(p, w, P, W) = 1/2 |P|^2 - 1/2 |W|^2 + U(p, w)
  U(p, w)       = V(q) + w . grad V(q) - 1/2 |w|^2 ,    q = p - w

  pdot = P,  wdot = -W,  Pdot = -dU/dp,  Wdot = -dU/dw
  dU/dp = grad V(q) + Hess V(q) w
  dU/dw = -Hess V(q) w - w

  Free limit (V=0): p free particle, w ghost oscillator at frequency 1 (= M).
  On shell, w = qddot + grad V(q): the ghost amplitude is the violation of the
  ordinary Yang-Mills equation of motion.
"""
with open('reduced_system.txt', 'w') as fh:
    fh.write(txt)
print(txt)
