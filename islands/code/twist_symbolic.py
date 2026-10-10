"""
twist_symbolic.py -- exact (computer-algebra) formula for the Birkhoff twist of an elliptic
periodic orbit of  H = (px^2 + py^2)/2 + x^2 y^2 / 2  on the section {y = 0, py > 0, H = E},
and its numerical evaluation for the Dahlqvist-Russberg orbit at E = 1/2.

Symbolic layer (sympy):
  * the time-T flow map  Phi(z* + v) = z* + M v + B(v,v)/2 + C(v,v,v)/6 + O(4)   with
    M (4x4), B (4x4x4, symmetric), C (4x4x4x4, symmetric) as SYMBOLS;
  * the lift of section coordinates (xi1, xi2) = (x - x*, px - px*) to the energy surface;
  * the return time  T + dt(xi)  from  y(phi_{T+dt}(z0)) = 0, solved as a series to order 3;
  * the return map  P(xi) = (x, px)(phi_{T+dt}(z0))  to order 3;
  * the complex normal form and the first Birkhoff coefficient (Kuznetsov's formula), giving
        tau = Im( conj(lambda) * c1 )          in the canonical complex coordinate zeta.
  Everything up to here is exact algebra in the symbols; the result is a rational function.

Numerical layer (floats, NOT rigorous): integrate the 1st/2nd/3rd-order variational equations
along the orbit to get M, B, C, substitute, and print tau.  Cross-check against
twist_estimate.py (rotation number vs amplitude), which gave tau ~ -2.5e5.
"""
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
import itertools, time

# ============================================================ symbolic layer
xi1, xi2, dt = sp.symbols('xi1 xi2 dt')
XS, PXS, PYS, EE = sp.symbols('x_s px_s py_s E', positive=True)      # fixed point, y_s = 0
Msym = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'M{i}{j}'))
# symmetric B and C via sorted-index symbols
def Bsym(i, j, k):
    j, k = sorted((j, k)); return sp.Symbol(f'B{i}_{j}{k}')
def Csym(i, j, k, l):
    j, k, l = sorted((j, k, l)); return sp.Symbol(f'C{i}_{j}{k}{l}')

ORDER = 3
# ---- truncated bivariate series with sympy coefficients: {(a,b): coeff}, a+b <= ORDER
class TS:
    def __init__(self, d=None): self.d = {m: c for m, c in (d or {}).items() if c != 0}
    @staticmethod
    def const(c): return TS({(0, 0): sp.sympify(c)})
    @staticmethod
    def var(i): return TS({(1, 0) if i == 0 else (0, 1): sp.Integer(1)})
    def __add__(self, o):
        o = o if isinstance(o, TS) else TS.const(o)
        d = dict(self.d)
        for m, c in o.d.items(): d[m] = sp.expand(d.get(m, 0) + c)
        return TS(d)
    __radd__ = __add__
    def __neg__(self): return TS({m: -c for m, c in self.d.items()})
    def __sub__(self, o): return self + (-(o if isinstance(o, TS) else TS.const(o)))
    def __rsub__(self, o): return TS.const(o) - self
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
        for n in range(ORDER + 1):
            acc = acc + term * ((-1) ** n); term = term * r
        return acc * (1 / c0)
def poly_at(expr, Zs):
    """evaluate a sympy polynomial in (x,y,px,py) at series Zs (list of TS)."""
    P = sp.Poly(sp.expand(expr), x_, y_, px_, py_)
    out = TS.const(0)
    for (a, b, c, d), coef in P.terms():
        t = TS.const(coef)
        for base, e in ((Zs[0], a), (Zs[1], b), (Zs[2], c), (Zs[3], d)):
            if e: t = t * (base ** e)
        out = out + t
    return out

X1, X2 = TS.var(0), TS.var(1)
# lift: py(px_s + xi2) = sqrt(2E - px^2) as a series in xi2 to ORDER, written with py_s
pyser = sp.series(sp.sqrt(2*EE - (PXS + xi2)**2), xi2, 0, ORDER + 1).removeO()
pyser = sp.expand(pyser.subs(sp.sqrt(2*EE - PXS**2), PYS))
py_ts = TS({(0, b): c for (b,), c in sp.Poly(pyser, xi2).terms()})
v = [X1, TS.const(0), X2, py_ts - PYS]

Phi = []
for i in range(4):
    acc = TS.const(0)
    for j in range(4):
        if v[j].d: acc = acc + v[j] * Msym[i, j]
    for j in range(4):
        for k in range(4):
            if v[j].d and v[k].d: acc = acc + (v[j] * v[k]) * (sp.Rational(1, 2) * Bsym(i, j, k))
    for j in range(4):
        for k in range(4):
            for l in range(4):
                if v[j].d and v[k].d and v[l].d:
                    acc = acc + (v[j] * v[k] * v[l]) * (sp.Rational(1, 6) * Csym(i, j, k, l))
    Phi.append(acc)
    print(f"  Phi[{i}] built: {len(acc.d)} monomials", flush=True)
zstar = [XS, sp.Integer(0), PXS, PYS]
Z = [Phi[i] + zstar[i] for i in range(4)]

x_, y_, px_, py_ = sp.symbols('x y px py')
fvec = [px_, py_, -x_*y_**2, -x_**2*y_]
def Dflow(expr):
    return sp.expand(sum(sp.diff(expr, s) * fs for s, fs in zip((x_, y_, px_, py_), fvec)))
f1 = [sp.expand(c) for c in fvec]; f2 = [Dflow(c) for c in f1]; f3 = [Dflow(c) for c in f2]
F1 = [poly_at(c, Z) for c in f1]; F2 = [poly_at(c, Z) for c in f2]; F3 = [poly_at(c, Z) for c in f3]
print("  F1,F2,F3 at Z built", flush=True)

invF1y = F1[1].inv()
dts = TS.const(0)
for _ in range(ORDER):
    dts = -(Z[1] + F2[1] * (dts ** 2) * sp.Rational(1, 2) + F3[1] * (dts ** 3) * sp.Rational(1, 6)) * invF1y
print("  return time series built", flush=True)
def flow_comp(i):
    return Z[i] + F1[i] * dts + F2[i] * (dts ** 2) * sp.Rational(1, 2) + F3[i] * (dts ** 3) * sp.Rational(1, 6)
P1 = flow_comp(0) - XS
P2 = flow_comp(2) - PXS
print("symbolic return map built", flush=True)
def coeffs(ts): return dict(ts.d)
A1, A2 = coeffs(P1), coeffs(P2)
# fixed-point condition: constant terms must vanish (they do identically only if Phi(z*) = z*,
# which is encoded by M, B, C being derivatives at the fixed point; we drop them explicitly)
for A in (A1, A2):
    A.pop((0, 0), None)

# ----------------------------------------------------------- complex normal form
# linear part L = [[A1(1,0), A1(0,1)], [A2(1,0), A2(0,1)]].  Numerical from here on for the
# frame (eigenvectors are algebraic but ugly); the coefficient FORMULA stays exact in g_jk.
lam, g20, g11, g02, g21 = sp.symbols('lambda g20 g11 g02 g21')
lamb = sp.conjugate(lam)
c1_formula = (g20*g11*(lamb - 3 + 2*lam) / (2*(lam**2 - lam)*(lamb - 1))
              + g11*sp.conjugate(g11) / (1 - lamb)
              + g02*sp.conjugate(g02) / (2*(lam**2 - lamb))
              + g21 / 2)
print("c1 =", c1_formula)
print("tau = Im(conj(lambda) * c1)   [Kuznetsov (4.21), area-preserving case: Re(conj(lambda) c1) = 0]")

# ============================================================ numerical layer
E = 0.5
xs, pxs = 3.14640122769837, 0.00179319341881
T = 31.053330268585066
pys = np.sqrt(2*E - pxs**2)

def rhs(t, s):
    z = s[:4]; M = s[4:20].reshape(4, 4); B = s[20:84].reshape(4, 4, 4); C = s[84:340].reshape(4, 4, 4, 4)
    x, y, px, py = z
    A = np.zeros((4, 4)); A[0, 2] = A[1, 3] = 1
    A[2, 0] = -y*y; A[2, 1] = -2*x*y; A[3, 0] = -2*x*y; A[3, 1] = -x*x
    A2 = np.zeros((4, 4, 4))
    A2[2, 0, 1] = A2[2, 1, 0] = -2*y; A2[2, 1, 1] = -2*x
    A2[3, 0, 0] = -2*y;               A2[3, 0, 1] = A2[3, 1, 0] = -2*x
    A3 = np.zeros((4, 4, 4, 4))
    for perm in itertools.permutations((0, 1, 1)): A3[(2,) + perm] = -2
    for perm in itertools.permutations((0, 0, 1)): A3[(3,) + perm] = -2
    dz = [px, py, -x*y*y, -x*x*y]
    dM = A @ M
    dB = np.einsum('ia,ajk->ijk', A, B) + np.einsum('iab,aj,bk->ijk', A2, M, M)
    MB = np.einsum('iab,aj,bkl->ijkl', A2, M, B)
    dC = (np.einsum('ia,ajkl->ijkl', A, C) + MB + MB.transpose(0, 2, 1, 3) + MB.transpose(0, 3, 2, 1)
          + np.einsum('iabc,aj,bk,cl->ijkl', A3, M, M, M))
    return np.concatenate([dz, dM.ravel(), dB.ravel(), dC.ravel()])

s0 = np.concatenate([[xs, 0.0, pxs, pys], np.eye(4).ravel(), np.zeros(64), np.zeros(256)])
t0 = time.time()
sol = solve_ivp(rhs, (0, T), s0, method='DOP853', rtol=1e-12, atol=1e-14)
sT = sol.y[:, -1]
zT = sT[:4]; Mn = sT[4:20].reshape(4, 4); Bn = sT[20:84].reshape(4, 4, 4); Cn = sT[84:340].reshape(4, 4, 4, 4)
print(f"\nvariational integration done ({time.time()-t0:.0f}s); |Phi(z*) - z*| = {np.linalg.norm(zT - s0[:4]):.2e}")
print("tr M =", np.trace(Mn), " (HI =", (np.trace(Mn) - 2)/2, ")")

subs = {XS: xs, PXS: pxs, PYS: pys, EE: E}
for i in range(4):
    for j in range(4):
        subs[Msym[i, j]] = Mn[i, j]
        for k in range(j, 4):
            subs[Bsym(i, j, k)] = Bn[i, j, k]
            for l in range(k, 4):
                subs[Csym(i, j, k, l)] = Cn[i, j, k, l]
a1 = {m: complex(sp.N(c.subs(subs))) .real for m, c in A1.items()}
a2 = {m: complex(sp.N(c.subs(subs))) .real for m, c in A2.items()}
L = np.array([[a1[(1, 0)], a1[(0, 1)]], [a2[(1, 0)], a2[(0, 1)]]])
print("DP (from symbolic reduction) =\n", L, "\ndet =", np.linalg.det(L), " HI =", np.trace(L)/2)

# canonical complex coordinate: same construction as twist_estimate.py
theta = np.arccos(np.trace(L)/2)
lamv, V = np.linalg.eig(L)
k = np.argmin(abs(lamv - np.exp(1j*theta)))
u, w = V[:, k].real, V[:, k].imag
area = u[0]*w[1] - u[1]*w[0]
if area < 0: w = -w; area = -area
S = np.column_stack([u, w]) / np.sqrt(area); Sinv = np.linalg.inv(S)
# xi-vector = S (xi, eta),  zeta = xi + i eta  =>  (xi1, xi2) = S @ [Re zeta, Im zeta]
# Express P in zeta:  zeta' = xi' + i eta' with (xi', eta') = Sinv @ P(S @ (Re z, Im z)).
# Expand numerically via sympy in (zr, zi), then convert to (zeta, zetabar).
zr, zi = sp.symbols('zr zi', real=True)
sub_xi = {xi1: S[0, 0]*zr + S[0, 1]*zi, xi2: S[1, 0]*zr + S[1, 1]*zi}
P1n = sum(c * xi1**a * xi2**b for (a, b), c in a1.items()).subs(sub_xi)
P2n = sum(c * xi1**a * xi2**b for (a, b), c in a2.items()).subs(sub_xi)
xi_p = Sinv[0, 0]*P1n + Sinv[0, 1]*P2n
eta_p = Sinv[1, 0]*P1n + Sinv[1, 1]*P2n
zeta, zetab = sp.symbols('zeta zetabar')
zeta_p = sp.expand((xi_p + sp.I*eta_p).subs({zr: (zeta + zetab)/2, zi: (zeta - zetab)/(2*sp.I)}))
Pz = sp.Poly(zeta_p, zeta, zetab)
G = {m: complex(c) for m, c in Pz.terms()}
lam_num = G[(1, 0)]
print("\nlinear coefficient of zeta:", lam_num, " |lambda| =", abs(lam_num), " (zetabar coeff, should be ~0:", abs(G.get((0, 1), 0)), ")")
def g(j, k):   # Kuznetsov normalization: coefficient = g_jk / (j! k!)
    from math import factorial
    return G.get((j, k), 0) * factorial(j) * factorial(k)
c1 = complex(c1_formula.subs({lam: lam_num, g20: g(2, 0), g11: g(1, 1), g02: g(0, 2), g21: g(2, 1)}).evalf())
tau = (np.conj(lam_num) * c1).imag
print(f"\nc1 = {c1}")
print(f"Re(conj(lambda) c1) = {(np.conj(lam_num)*c1).real:.3e}   (must be ~0 for an area-preserving map)")
print(f"TWIST  tau = Im(conj(lambda) c1) = {tau:.6e}      (twist_estimate.py slope: ~ -2.5e5)")
