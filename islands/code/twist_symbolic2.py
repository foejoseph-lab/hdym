"""
twist_symbolic2.py -- directional form of the exact twist formula, exported for the checker.

Return map P on the section {y=0, py>0, H=E} with coordinates xi = (x - x_s, px - px_s).
Lift to the energy surface:  z0 - z* = xi1 e_x + xi2 e_2 + c2 xi2^2 e_py + c3 xi2^3 e_py
  e_x = (1,0,0,0),  e_py = (0,0,0,1),  e_2 = (0,0,1,c1),
  c1 = -px_s/py_s,  c2 = -E/py_s^3 (= -(py_s^2+px_s^2)/(2 py_s^3)),  c3 = -px_s E/py_s^5 ... (from the sqrt series)
Flow map to third order in directional form (all 4-vectors, evaluated at the orbit's period T):
  Mx  = Dphi e_x,  M2 = Dphi e_2,  Mpy = Dphi e_py
  Bxx = D2phi(e_x,e_x), Bx2 = D2phi(e_x,e_2), B22 = D2phi(e_2,e_2), Bxpy = D2phi(e_x,e_py), B2py = D2phi(e_2,e_py)
  Cxxx, Cxx2, Cx22, C222 = D3phi on the corresponding triples
Then  Phi(z0) - z* = xi1 Mx + xi2 M2 + (c2 xi2^2 + c3 xi2^3) Mpy
                   + 1/2 [xi1^2 Bxx + 2 xi1 xi2 Bx2 + xi2^2 B22 + 2 c2 xi1 xi2^2 Bxpy + 2 c2 xi2^3 B2py]
                   + 1/6 [xi1^3 Cxxx + 3 xi1^2 xi2 Cxx2 + 3 xi1 xi2^2 Cx22 + xi2^3 C222]   + O(4),
and the return-time / section reduction as in twist_symbolic.py.

Output: twist_formula.json -- for each of P1, P2 the Taylor coefficients (monomials xi1^a xi2^b,
a+b in 1..3) as expression trees over the symbols above plus x_s, px_s, py_s, E.
Tree format: ["+", t1, t2, ...] | ["*", f1, f2, ...] | ["^", base, int] | ["q", "num/den"] | ["s", name]
Also prints the float evaluation of tau for the DR orbit as a cross-check (expect -2.140e5).
"""
import json, itertools, time
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

xi1, xi2 = sp.symbols('xi1 xi2')
XS, PXS, PYS, EE = sp.symbols('x_s px_s py_s E', positive=True)
VEC = ['Mx', 'M2', 'Mpy', 'Bxx', 'Bx2', 'B22', 'Bxpy', 'B2py', 'Cxxx', 'Cxx2', 'Cx22', 'C222']
V = {n: [sp.Symbol(f'{n}_{i}') for i in range(4)] for n in VEC}
ORDER = 3

class TS:
    def __init__(self, d=None): self.d = {m: c for m, c in (d or {}).items() if c != 0}
    @staticmethod
    def const(c): return TS({(0, 0): sp.sympify(c)})
    def __add__(self, o):
        o = o if isinstance(o, TS) else TS.const(o); d = dict(self.d)
        for m, c in o.d.items(): d[m] = sp.expand(d.get(m, 0) + c)
        return TS(d)
    __radd__ = __add__
    def __neg__(self): return TS({m: -c for m, c in self.d.items()})
    def __sub__(self, o): return self + (-(o if isinstance(o, TS) else TS.const(o)))
    def __mul__(self, o):
        if not isinstance(o, TS): return TS({m: sp.expand(c * o) for m, c in self.d.items()})
        d = {}
        for (a1, b1), c1 in self.d.items():
            for (a2, b2), c2 in o.d.items():
                if a1 + a2 + b1 + b2 <= ORDER:
                    m = (a1 + a2, b1 + b2); d[m] = d.get(m, 0) + c1 * c2
        return TS({m: sp.expand(c) for m, c in d.items()})
    __rmul__ = __mul__
    def __pow__(self, n):
        r = TS.const(1)
        for _ in range(n): r = r * self
        return r
    def c0(self): return self.d.get((0, 0), sp.Integer(0))
    def inv(self):
        c0 = self.c0(); r = (self - c0) * (1 / c0)
        acc = TS.const(0); term = TS.const(1)
        for n in range(ORDER + 1): acc = acc + term * ((-1) ** n); term = term * r
        return acc * (1 / c0)

x_, y_, px_, py_ = sp.symbols('x y px py')
def poly_at(expr, Zs):
    P = sp.Poly(sp.expand(expr), x_, y_, px_, py_); out = TS.const(0)
    for (a, b, c, d), coef in P.terms():
        t = TS.const(coef)
        for base, e in ((Zs[0], a), (Zs[1], b), (Zs[2], c), (Zs[3], d)):
            if e: t = t * (base ** e)
        out = out + t
    return out

X1 = TS({(1, 0): sp.Integer(1)}); X2 = TS({(0, 1): sp.Integer(1)})
# py(px_s + xi2) series -> c1, c2, c3 in terms of px_s, py_s, E
pyser = sp.series(sp.sqrt(2*EE - (PXS + xi2)**2), xi2, 0, ORDER + 1).removeO()
pyser = sp.expand(pyser.subs(sp.sqrt(2*EE - PXS**2), PYS))
pc = {b: sp.simplify(c) for (b,), c in sp.Poly(pyser, xi2).terms()}
c1, c2, c3 = pc[1], pc[2], pc[3]
print("c1 =", c1, "\nc2 =", c2, "\nc3 =", c3)

def vec_ts(name, coeff_ts):
    """4-vector of TS: coeff_ts * V[name]"""
    return [coeff_ts * V[name][i] for i in range(4)]
def vadd(*vs):
    out = [TS.const(0)] * 4
    for v in vs: out = [a + b for a, b in zip(out, v)]
    return out

half, sixth = sp.Rational(1, 2), sp.Rational(1, 6)
Phi = vadd(vec_ts('Mx', X1), vec_ts('M2', X2), vec_ts('Mpy', X2**2 * c2 + X2**3 * c3),
           vec_ts('Bxx', X1**2 * half), vec_ts('Bx2', X1 * X2), vec_ts('B22', X2**2 * half),
           vec_ts('Bxpy', X1 * X2**2 * c2), vec_ts('B2py', X2**3 * c2),
           vec_ts('Cxxx', X1**3 * sixth), vec_ts('Cxx2', X1**2 * X2 * half),
           vec_ts('Cx22', X1 * X2**2 * half), vec_ts('C222', X2**3 * sixth))
zstar = [XS, sp.Integer(0), PXS, PYS]
Z = [Phi[i] + zstar[i] for i in range(4)]
fvec = [px_, py_, -x_*y_**2, -x_**2*y_]
def Dflow(e): return sp.expand(sum(sp.diff(e, s) * fs for s, fs in zip((x_, y_, px_, py_), fvec)))
f1 = [sp.expand(c) for c in fvec]; f2 = [Dflow(c) for c in f1]; f3 = [Dflow(c) for c in f2]
F1 = [poly_at(c, Z) for c in f1]; F2 = [poly_at(c, Z) for c in f2]; F3 = [poly_at(c, Z) for c in f3]
invF1y = F1[1].inv()
dts = TS.const(0)
for _ in range(ORDER):
    dts = -(Z[1] + F2[1] * (dts ** 2) * half + F3[1] * (dts ** 3) * sixth) * invF1y
def flow_comp(i): return Z[i] + F1[i] * dts + F2[i] * (dts ** 2) * half + F3[i] * (dts ** 3) * sixth
P1 = flow_comp(0) - XS; P2 = flow_comp(2) - PXS
A1 = {m: c for m, c in P1.d.items() if m != (0, 0)}
A2 = {m: c for m, c in P2.d.items() if m != (0, 0)}
print("directional return map built:", len(A1), len(A2), "monomials")

# ---------------------------------------------------------- export as expression trees
def tree(e):
    e = sp.sympify(e)
    if e.is_Rational: return ["q", f"{e.p}/{e.q}"]
    if e.is_Symbol: return ["s", e.name]
    if e.is_Add: return ["+"] + [tree(a) for a in e.args]
    if e.is_Mul: return ["*"] + [tree(a) for a in e.args]
    if e.is_Pow and e.exp.is_Integer: return ["^", tree(e.base), int(e.exp)]
    raise ValueError(f"unsupported node: {e} ({type(e)})")
out = {"symbols": [f'{n}_{i}' for n in VEC for i in range(4)] + ['x_s', 'px_s', 'py_s', 'E'],
       "P1": {f"{a},{b}": tree(c) for (a, b), c in A1.items()},
       "P2": {f"{a},{b}": tree(c) for (a, b), c in A2.items()},
       "note": "P_i = sum coeff[a,b] xi1^a xi2^b; constant terms vanish at the fixed point and are omitted."}
json.dump(out, open("twist_formula.json", "w"))
print("wrote twist_formula.json")

# ---------------------------------------------------------- float cross-check (directional variational ODEs)
E = 0.5; xs, pxs = 3.14640122769837, 0.00179319341881; T = 31.053330268585066; pys = np.sqrt(2*E - pxs**2)
c1n, c2n, c3n = [float(c.subs({XS: xs, PXS: pxs, PYS: pys, EE: E})) for c in (c1, c2, c3)]
ex = np.array([1., 0, 0, 0]); epy = np.array([0, 0, 0, 1.]); e2 = np.array([0, 0, 1., c1n])
pairs = [('Bxx', 0, 0), ('Bx2', 0, 1), ('B22', 1, 1), ('Bxpy', 0, 2), ('B2py', 1, 2)]
trip = [('Cxxx', 0, 0, 0), ('Cxx2', 0, 0, 1), ('Cx22', 0, 1, 1), ('C222', 1, 1, 1)]
def A_of(z):
    x, y, px, py = z
    A = np.zeros((4, 4)); A[0, 2] = A[1, 3] = 1; A[2, 0] = -y*y; A[2, 1] = -2*x*y; A[3, 0] = -2*x*y; A[3, 1] = -x*x
    return A
def A2v(z, u, v):
    x, y = z[0], z[1]
    return np.array([0, 0, -2*y*(u[0]*v[1] + u[1]*v[0]) - 2*x*u[1]*v[1], -2*y*u[0]*v[0] - 2*x*(u[0]*v[1] + u[1]*v[0])])
def A3v(u, v, w):
    return np.array([0, 0, -2*(u[0]*v[1]*w[1] + u[1]*v[0]*w[1] + u[1]*v[1]*w[0]),
                     -2*(u[0]*v[0]*w[1] + u[0]*v[1]*w[0] + u[1]*v[0]*w[0])])
def rhs(t, s):
    z = s[:4]; m = s[4:16].reshape(3, 4); b = s[16:36].reshape(5, 4); c = s[36:52].reshape(4, 4)
    A = A_of(z)
    dz = [z[2], z[3], -z[0]*z[1]**2, -z[0]**2*z[1]]
    dm = (A @ m.T).T
    db = np.array([A @ b[k] + A2v(z, m[i], m[j]) for k, (_, i, j) in enumerate(pairs)])
    bidx = {(0, 0): 0, (0, 1): 1, (1, 1): 2}
    dc = []
    for k, (_, i, j, l) in enumerate(trip):
        bij, bil, bjl = b[bidx[tuple(sorted((i, j)))]], b[bidx[tuple(sorted((i, l)))]], b[bidx[tuple(sorted((j, l)))]]
        dc.append(A @ c[k] + A2v(z, m[i], bjl) + A2v(z, m[j], bil) + A2v(z, m[l], bij) + A3v(m[i], m[j], m[l]))
    return np.concatenate([dz, dm.ravel(), db.ravel(), np.array(dc).ravel()])
s0 = np.concatenate([[xs, 0, pxs, pys], np.array([ex, e2, epy]).ravel(), np.zeros(20), np.zeros(16)])
sol = solve_ivp(rhs, (0, T), s0, method='DOP853', rtol=1e-12, atol=1e-14)
sT = sol.y[:, -1]; m = sT[4:16].reshape(3, 4); b = sT[16:36].reshape(5, 4); c = sT[36:52].reshape(4, 4)
vals = {'x_s': xs, 'px_s': pxs, 'py_s': pys, 'E': E}
for i in range(4):
    vals[f'Mx_{i}'], vals[f'M2_{i}'], vals[f'Mpy_{i}'] = m[0, i], m[1, i], m[2, i]
    for k, (n, *_r) in enumerate(pairs): vals[f'{n}_{i}'] = b[k, i]
    for k, (n, *_r) in enumerate(trip): vals[f'{n}_{i}'] = c[k, i]
def ev(t):
    if t[0] == 'q': p, q = t[1].split('/'); return int(p) / int(q)
    if t[0] == 's': return vals[t[1]]
    if t[0] == '+': return sum(ev(a) for a in t[1:])
    if t[0] == '*':
        r = 1.0
        for a in t[1:]: r *= ev(a)
        return r
    if t[0] == '^': return ev(t[1]) ** t[2]
a1 = {tuple(map(int, k.split(','))): ev(t) for k, t in out["P1"].items()}
a2 = {tuple(map(int, k.split(','))): ev(t) for k, t in out["P2"].items()}
L = np.array([[a1[(1, 0)], a1[(0, 1)]], [a2[(1, 0)], a2[(0, 1)]]])
print("DP =", L.tolist(), " det", np.linalg.det(L), " HI", np.trace(L)/2)
# complexify exactly as twist_symbolic.py / as the checker will
tr = np.trace(L)/2; sn = np.sqrt(1 - tr*tr); a, bb = L[0, 0], L[0, 1]
u = np.array([bb, tr - a]); w = np.array([0.0, sn]); area = u[0]*w[1] - u[1]*w[0]
if area < 0: w = -w; area = -area
S = np.column_stack([u, w]) / np.sqrt(area); Sinv = np.linalg.inv(S)
zr, zi = sp.symbols('zr zi', real=True)
sub = {xi1: S[0, 0]*zr + S[0, 1]*zi, xi2: S[1, 0]*zr + S[1, 1]*zi}
P1n = sum(v * xi1**p * xi2**q for (p, q), v in a1.items()).subs(sub)
P2n = sum(v * xi1**p * xi2**q for (p, q), v in a2.items()).subs(sub)
zeta, zetab = sp.symbols('zeta zetabar')
zp = sp.expand(((Sinv[0, 0]*P1n + Sinv[0, 1]*P2n) + sp.I*(Sinv[1, 0]*P1n + Sinv[1, 1]*P2n))
               .subs({zr: (zeta + zetab)/2, zi: (zeta - zetab)/(2*sp.I)}))
G = {mm: complex(cc) for mm, cc in sp.Poly(zp, zeta, zetab).terms()}
lam = G[(1, 0)]
from math import factorial
def g(j, k): return G.get((j, k), 0) * factorial(j) * factorial(k)
lb = np.conj(lam)
c1v = (g(2, 0)*g(1, 1)*(lb - 3 + 2*lam) / (2*(lam**2 - lam)*(lb - 1)) + abs(g(1, 1))**2/(1 - lb)
       + abs(g(0, 2))**2/(2*(lam**2 - lb)) + g(2, 1)/2)
print(f"|lambda| = {abs(lam):.12f}   Re(conj(lam) c1) = {(lb*c1v).real:.3e}")
print(f"TWIST (directional formula, float) tau = {(lb*c1v).imag:.6e}    [twist_symbolic.py gave -2.140067e+05]")
