"""
bianchi9_kasner.py -- the Weyl ghost on a Kasner epoch (Bianchi I limit of the system in
bianchi9_system.pkl): linear perturbation theory with the onset as a function of H/M.

1. L_BI = the e^{3 alpha}-homogeneous part of L (spatial-curvature terms, which scale as
   e^{alpha} and e^{-alpha}, dropped).  Check: with the Weyl term off it is the Bianchi I
   Misner Lagrangian 3 e^{3a}(-a'^2 + b'^2).
2. Background: vacuum Kasner, alpha = ln(t)/3, beta = (ln t /3)(cos th, sin th), proper time,
   H = 1/(3t).  Check: solves the full quadratic-gravity equations (Einstein manifold).
3. Perturbations (da, db+, db-) -> quadratic Lagrangian -> linear equations.  Eliminate da
   with the linearised constraint and alpha equation.  Result: two coupled 4th-order linear
   ODEs for (db+, db-) with coefficients polynomial in 1/t and M.
4. t -> 0 (H >> M): keep the leading powers; Euler-type system; indicial exponents s in
   db ~ t^s as functions of th.  Compared with the GR modes (s = 0 and the Kasner-rotation
   mode) this says whether the Einstein branch is approached or left as t -> 0.
5. t -> inf (H << M): WKB; frequency M, amplitude ~ t^{-1/2} (= e^{-3 alpha/2}) -- the adiabatic
   ghost.
6. Numerics: integrate the linear system in x = M t from large x (adiabatic ghost state of
   unit action) to small x; record the ghost action and the projection on the GR modes.  Where
   the action stops being conserved is the onset; the growth law at small x is the indicial
   exponent.

Output: printed exponents, kasner_onset.npz, kasner_onset.png.
Run:   %run bianchi9_kasner.py
"""
import pickle, time
import numpy as np
import sympy as sp

t0 = time.time()
S = pickle.load(open('bianchi9_system.pkl', 'rb'))
M = sp.Symbol(S['M'], positive=True)
KMAX = 6
J = {nm: [sp.Symbol(f'{nm}{k}') for k in range(KMAX + 1)] for nm in ('al', 'bp', 'bm', 'N')}
def D(expr):
    return sum(sp.diff(expr, J[nm][k]) * J[nm][k+1] for nm in J for k in range(KMAX))
gauge = {J['N'][0]: 1, J['N'][1]: 0, J['N'][2]: 0, J['N'][3]: 0}
L = sp.expand(S['L'].subs(gauge))

# ---------------------------------------------------------------- 1. Bianchi I part
c = sp.Symbol('c')
Lc = sp.expand(L.subs(J['al'][0], J['al'][0] + c))
L_BI = sp.S(0)
for term in sp.Add.make_args(Lc):
    # terms scaling as exp(3c) exactly
    if sp.simplify(term / sp.exp(3 * c)).has(c):
        continue
    L_BI += term / sp.exp(3 * c)
L_BI = sp.expand(L_BI)
al0, al1, al2 = J['al'][:3]; bp0, bp1, bp2, bp3 = J['bp'][:4]; bm0, bm1, bm2, bm3 = J['bm'][:4]
print("   L_BI / e^{3 alpha} =")
sp.pprint(sp.collect(sp.expand(L_BI * sp.exp(-3 * al0)), M))
L_GR_BI = sp.limit(L_BI, M, sp.oo)
print("   GR part:", sp.factor(L_GR_BI))
print(f"[{time.time()-t0:5.1f}s] Bianchi I Lagrangian")

# ---------------------------------------------------------------- 2. Kasner background
t, th = sp.symbols('t theta', positive=True)
bg = {J['al'][0]: sp.log(t) / 3, J['bp'][0]: sp.cos(th) * sp.log(t) / 3, J['bm'][0]: sp.sin(th) * sp.log(t) / 3}
def jet_of(funcs):
    """dict jet symbol -> expression, for functions of t"""
    out = {}
    for nm, f in funcs.items():
        for k in range(KMAX + 1):
            out[J[nm][k]] = sp.diff(f, t, k)
    return out
bgjet = jet_of({'al': bg[J['al'][0]], 'bp': bg[J['bp'][0]], 'bm': bg[J['bm'][0]]})
def EL(Lx, nm, order=4):
    out = sp.S(0)
    for k in range(order + 1):
        term = sp.diff(Lx, J[nm][k])
        for _ in range(k): term = D(term)
        out += (-1)**k * term
    return sp.expand(out)
E = {nm: EL(L_BI, nm) for nm in ('al', 'bp', 'bm')}
Con = EL(L.subs(gauge), 'N') if False else None
# constraint for L_BI: dL/dN at N=1 -- rebuild from the full L (N-dependence) restricted to e^{3 alpha} part
LN = sp.expand(S['L'].subs(J['al'][0], J['al'][0] + c))
L_BI_N = sp.S(0)
for term in sp.Add.make_args(LN):
    if sp.simplify(term / sp.exp(3 * c)).has(c): continue
    L_BI_N += term / sp.exp(3 * c)
E['N'] = EL(sp.expand(L_BI_N), 'N').subs(gauge)
for nm in ('al', 'bp', 'bm', 'N'):
    val = sp.simplify(E[nm].subs(bgjet))
    print(f"   Kasner solves E_{nm}: {val}")
print(f"[{time.time()-t0:5.1f}s] background checked")

# ---------------------------------------------------------------- 3. linearise
eps = sp.Symbol('epsilon')
da, dbp, dbm = [sp.Function(s)(t) for s in ('da', 'dbp', 'dbm')]
pert = jet_of({'al': bg[J['al'][0]] + eps * da, 'bp': bg[J['bp'][0]] + eps * dbp, 'bm': bg[J['bm'][0]] + eps * dbm})
lin = {}
for nm in ('al', 'bp', 'bm', 'N'):
    ex = E[nm].subs(pert)
    lin[nm] = sp.expand(sp.diff(ex, eps).subs(eps, 0))
    lin[nm] = sp.simplify(lin[nm] * t**4 * sp.exp(-3 * bg[J['al'][0]]) )     # clear the t^-n and e^{3 alpha}
print(f"[{time.time()-t0:5.1f}s] linearised; e.g. linearised constraint:")
sp.pprint(sp.collect(sp.expand(lin['N']), [da, dbp, dbm]))

with open('kasner_linear.pkl', 'wb') as fh:
    pickle.dump(dict(L_BI=L_BI, lin={k: str(v) for k, v in lin.items()}, lin_expr=lin), fh)

# ---------------------------------------------------------------- 4. indicial exponents, t -> 0 (symbolic)
# The M -> 0 (pure Weyl) leading operator is degenerate: its determinant is s^3 (s-2)^3 Q(s,th)
# with Q == 0 identically (conformal invariance: the pure-C^2 operator does not see the Kasner
# background in the ghost directions at leading order). The exponents come from the next
# order in (M t)^2 -> degenerate perturbation problem; done numerically below.
s = sp.Symbol('s')
A, B, C = sp.symbols('A B C')
rows = []
for nm in ('al', 'bp', 'bm'):
    ex = lin[nm].subs({da: A * t**s, dbp: B * t**s, dbm: C * t**s}).doit()
    rows.append(sp.expand(sp.simplify(ex / t**s)))
Mat0 = sp.Matrix([[sp.diff(sp.expand(r * M**2).subs(M, 0), v) for v in (A, B, C)] for r in rows]).applyfunc(lambda e: sp.simplify(e.subs(t, 1)))
print("   t -> 0 (pure Weyl) indicial determinant =", sp.factor(sp.trigsimp(Mat0.det())))
MatGR = sp.Matrix([[sp.diff(sp.limit(r, M, sp.oo), v) for v in (A, B, C)] for r in rows]).applyfunc(lambda e: sp.simplify(e.subs(t, 1)))
print("   GR-limit indicial determinant =", sp.factor(MatGR.det()), " roots", sp.solve(MatGR.det(), s))
print(f"[{time.time()-t0:5.1f}s] indicial (symbolic) done")

# ---------------------------------------------------------------- 5. numerical linear system in x = M t
# unknowns: da (2nd order), dbp, dbm (4th order). Solve the three linearised equations for the
# highest derivatives, build the 10-d first-order system in x with M = 1.
x = sp.Symbol('x', positive=True)
u = {}
for nm_, f in (('a', da), ('p', dbp), ('m', dbm)):
    for k in range(5):
        u[(nm_, k)] = sp.Symbol(f'd{nm_}{k}')
def to_u(expr):
    rep = {}
    for nm_, f in (('a', da), ('p', dbp), ('m', dbm)):
        for k in range(4, 0, -1):
            rep[sp.Derivative(f, (t, k))] = u[(nm_, k)]
    expr = expr.subs(rep)
    for nm_, f in (('a', da), ('p', dbp), ('m', dbm)):
        expr = expr.subs(f, u[(nm_, 0)])
    return expr.subs({M: 1, t: x})
eqs = [to_u(lin[nm]) for nm in ('al', 'bp', 'bm')]
for nm, e in zip(('al', 'bp', 'bm'), eqs):
    print(f"   lin {nm} depends on:", sorted({str(v) for v in e.free_symbols if str(v).startswith('d')}))
da2 = sp.solve(eqs[0], u[('a', 2)])
assert len(da2) == 1; da2 = sp.simplify(da2[0])
def Dx(expr):
    out = sp.diff(expr, x)
    for nm_ in ('a', 'p', 'm'):
        for k in range(4):
            out += sp.diff(expr, u[(nm_, k)]) * u[(nm_, k + 1)]
    return out
da3 = Dx(da2).subs(u[('a', 2)], da2)
eb = [sp.expand(e.subs(u[('a', 3)], da3).subs(u[('a', 2)], da2)) for e in eqs[1:]]
Am, rv = sp.linear_eq_to_matrix(eb, [u[('p', 4)], u[('m', 4)]])
b4 = Am.LUsolve(rv)
b4 = [sp.simplify(sp.together(v).subs(u[('a', 2)], da2)) for v in b4]
state = [u[('a', 0)], u[('a', 1)], u[('p', 0)], u[('m', 0)], u[('p', 1)], u[('m', 1)], u[('p', 2)], u[('m', 2)], u[('p', 3)], u[('m', 3)]]
rhs = sp.Matrix([u[('a', 1)], da2, u[('p', 1)], u[('m', 1)], u[('p', 2)], u[('m', 2)], u[('p', 3)], u[('m', 3)], b4[0], b4[1]])
Jlin = rhs.jacobian(sp.Matrix(state))              # linear system: y' = Jlin(x) y
f_J = sp.lambdify([x, th], Jlin, 'numpy', cse=True)
print(f"[{time.time()-t0:5.1f}s] linear system in x built")

from scipy.integrate import solve_ivp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def run(theta, x_hi=60.0, x_lo=0.02, n=400):
    J0 = np.array(f_J(x_hi, theta), float)
    lam, V = np.linalg.eig(J0)
    ghost = np.where(np.abs(lam.imag) > 0.5)[0]
    grm = np.where(np.abs(lam.imag) <= 0.5)[0]
    F = lambda xx, Y: (np.array(f_J(xx, theta), float) @ Y.reshape(10, 10)).ravel()
    xs = np.geomspace(x_hi, x_lo, n)
    sol = solve_ivp(F, (x_hi, x_lo), np.eye(10).ravel(), t_eval=xs, method='DOP853', rtol=1e-10, atol=1e-13)
    Phi = sol.y.T.reshape(-1, 10, 10)
    Eg = np.column_stack([V[:, ghost].real, V[:, ghost].imag]); Eg = np.linalg.qr(Eg)[0]
    # scaled norm so that an adiabatic ghost reads as constant: with M = 1 the ghost has
    # components (b, b1, b2, b3) ~ t^(-1/2) (1,1,1,1); multiply by sqrt(x).
    def gnorm(Y, xx):
        return np.linalg.norm(Y[2:, :], axis=0) * np.sqrt(xx)
    gg = np.array([gnorm(Ph @ Eg, xx) for Ph, xx in zip(Phi, xs)])
    return xs, gg, lam[ghost], lam[grm]
fig, ax = plt.subplots(figsize=(7, 4.5))
out = {}
for theta in (np.pi, np.pi/2, 0.0, np.pi/3):
    xs, gg, lg, lgr = run(theta)
    HM = 1 / (3 * xs)
    gmax = gg.max(axis=1); gmin = gg.min(axis=1)
    ax.loglog(HM, gmax / gmax[0], label=f'theta={theta:.2f}')
    out[theta] = (HM, gmax, gmin)
    k = slice(-40, -1)
    slope = np.polyfit(np.log(xs[k]), np.log(gmax[k]), 1)[0]
    i2 = np.argmax(gmax / gmax[0] > 2); i10 = np.argmax(gmax / gmax[0] > 10)
    print(f"   theta={theta:.2f}: ghost freq at x_hi {lg.imag.max():.4f}; scaled ghost amplitude grows by "
          f"{gmax[-1]/gmax[0]:.2e} from H/M={HM[0]:.4f} to {HM[-1]:.2f}; small-x power law g ~ x^{slope:.2f}; "
          f"x2 at H/M={HM[i2]:.3f}, x10 at H/M={HM[i10]:.3f}")
ax.set_xlabel('H / M along the Kasner epoch (increasing toward the singularity)')
ax.set_ylabel('scaled ghost amplitude / initial (adiabatic = 1)'); ax.legend(); ax.axhline(1, color='k', lw=.5)
plt.tight_layout(); plt.savefig('kasner_onset.png', dpi=130)
np.savez('kasner_onset.npz', **{f'theta_{k:.3f}': np.array(v) for k, v in out.items()})
print(f"[{time.time()-t0:5.1f}s] wrote kasner_onset.png")
