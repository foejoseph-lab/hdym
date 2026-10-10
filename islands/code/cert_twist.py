"""
cert_twist.py -- TRUSTED layer, continued: directional second/third variations along the
certified orbit and the Birkhoff twist of the return map, in exact rational interval
arithmetic.  Imports the interval primitives from cert_check.py.

Lemmas (function names):
  VAR2/VAR3  Taylor recurrences for the directional variational equations
               m' = A m,   b' = A b + A2(m_a, m_b),   c' = A c + A2(m_i,b_jl) + A2(m_j,b_il) + A2(m_l,b_ij) + A3(m_i,m_j,m_l)
             with A = Df, A2 = D2f, A3 = D3f along the orbit (all polynomial; A3 is constant).
             Same Lagrange-remainder lemma as LAG, applied to the augmented polynomial ODE.
  FRAME2     accumulation of b, c across steps in the Lohner frame:  b = A_k beta_k,
               beta_{k+1} = T_k beta_k + A_{k+1}^{-1} Delta_k,   Delta_k = b-increment of the step with b(0) = 0
             (linearity of the inhomogeneous linear ODE in its initial datum).
  TWIST      the exact return-map formula (twist_formula.json, generated once by sympy in
             twist_symbolic2.py) evaluated on interval data; complexification in the canonical
             frame; Kuznetsov's normal-form coefficient; tau = Im(conj(lambda) c1).
             Non-resonance lambda^q != 1 (q = 1,2,3) is enforced by the divisions.
"""
import json, sys, os
from fractions import Fraction as Fr
# The interval primitives are taken from the RUNNING cert_check module, which cert_check.main
# passes in explicitly via bind().  (A plain `import cert_check` could create a second copy of the
# module whose rounding grid is unset; bind() makes the choice explicit under any runner.)
cc = None
def bind(module):
    global cc, I, add, sub, mul, neg, scal, div_int, hull, inside, width, ipow, isqrt
    global vadd, vsub, mvec, mmul, madd, mscal, eye, _conv, CheckFailed
    cc = module
    if cc.GRID is None:
        raise RuntimeError("cert_twist.bind: the rounding grid of the given cert_check module is not set")
    I, add, sub, mul, neg, scal, div_int, hull, inside, width, ipow, isqrt = \
        cc.I, cc.add, cc.sub, cc.mul, cc.neg, cc.scal, cc.div_int, cc.hull, cc.inside, cc.width, cc.ipow, cc.isqrt
    vadd, vsub, mvec, mmul, madd, mscal, eye, _conv, CheckFailed = \
        cc.vadd, cc.vsub, cc.mvec, cc.mmul, cc.madd, cc.mscal, cc.eye, cc._conv, cc.CheckFailed

# ------------------------------------------------------------- directional series (VAR2 / VAR3)
PAIRS = [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2)]          # (x,x) (x,2) (2,2) (x,py) (2,py)
TRIPS = [(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)]      # on {x, 2}
PAIR_IDX = {pr: k for k, pr in enumerate(PAIRS)}

def A_series(ser, p):
    """Taylor coefficients A_k of A(z(t)) = Df(z(t)) (m = 0)."""
    x, y, px, py, xx, yy = ser
    Ak = []
    for k in range(p + 1):
        xy = _conv(x, y, k)
        d = I(1 if k == 0 else 0)
        a20 = neg(yy[k]) if k < len(yy) else I(0)
        a31 = neg(xx[k]) if k < len(xx) else I(0)
        Ak.append([[I(0), I(0), d, I(0)], [I(0), I(0), I(0), d],
                   [a20, scal(Fr(-2), xy), I(0), I(0)], [scal(Fr(-2), xy), a31, I(0), I(0)]])
    return Ak

def _conv3(a, b, c, n):
    """sum_{i+j+k=n} a_i b_j c_k."""
    s = I(0)
    for i in range(n + 1):
        for j in range(n + 1 - i):
            s = add(s, mul(a[i], mul(b[j], c[n - i - j])))
    return s

def _A2_coeff(x, y, u, v, n):
    """n-th Taylor coefficient of A2(z(t))(u(t), v(t)) = D2f(z)[u, v]:
       comp2 = -2 y (u0 v1 + u1 v0) - 2 x u1 v1 ;  comp3 = -2 y u0 v0 - 2 x (u0 v1 + u1 v0)."""
    t01 = add(_conv3(y, u[0], v[1], n), _conv3(y, u[1], v[0], n))
    c2 = scal(Fr(-2), add(t01, _conv3(x, u[1], v[1], n)))
    s01 = add(_conv3(x, u[0], v[1], n), _conv3(x, u[1], v[0], n))
    c3 = scal(Fr(-2), add(_conv3(y, u[0], v[0], n), s01))
    return [I(0), I(0), c2, c3]

def _A3_coeff(u, v, w, n):
    """A3 = D3f is constant:  comp2 = -2(u0 v1 w1 + u1 v0 w1 + u1 v1 w0),  comp3 = -2(u0 v0 w1 + u0 v1 w0 + u1 v0 w0)."""
    c2 = scal(Fr(-2), add(add(_conv3(u[0], v[1], w[1], n), _conv3(u[1], v[0], w[1], n)), _conv3(u[1], v[1], w[0], n)))
    c3 = scal(Fr(-2), add(add(_conv3(u[0], v[0], w[1], n), _conv3(u[0], v[1], w[0], n)), _conv3(u[1], v[0], w[0], n)))
    return [I(0), I(0), c2, c3]

def _Av_coeff(Ak, vser, n):
    """n-th coefficient of A(t) v(t)."""
    out = [I(0)] * 4
    for k in range(n + 1):
        out = vadd(out, mvec(Ak[k], [vser[i][n - k] for i in range(4)]))
    return out

def dir_series(ser, p, m0, b_full0):
    """Taylor coefficient series (component-major: v[i][k]) up to order p of
         m_a (3 directions),  b_zero_ab (b(0)=0),  b_full_ab (b(0)=b_full0),  c_abc (c(0)=0, forced by m and b_full).
    Pairwise product series are cached so each coefficient costs O(n).  Lemma VAR2/VAR3."""
    x, y = ser[0], ser[1]
    Ak = A_series(ser, p)
    def init(v0): return [[v0[i]] for i in range(4)]
    m = [init(v) for v in m0]
    bz = [init([I(0)] * 4) for _ in PAIRS]
    bf = [init(v) for v in b_full0]
    c = [init([I(0)] * 4) for _ in TRIPS]
    # cached product series for (m_a, m_b): P00, P01s = u0 v1 + u1 v0, P11
    mm = {pr: ([], [], []) for pr in PAIRS}
    # cached product series for (m_i, b_full_k) needed by c
    mb_keys = set()
    for (i1, i2, i3) in TRIPS:
        mb_keys |= {(i1, PAIR_IDX[(i2, i3)]), (i2, PAIR_IDX[(i1, i3)]), (i3, PAIR_IDX[(i1, i2)])}
    mb = {key: ([], [], []) for key in mb_keys}
    def extend_products(n):
        for (a, bq), (P00, P01, P11) in mm.items():
            u, v = m[a], m[bq]
            P00.append(_conv(u[0], v[0], n)); P11.append(_conv(u[1], v[1], n))
            P01.append(add(_conv(u[0], v[1], n), _conv(u[1], v[0], n)))
        for (a, k), (Q00, Q01, Q11) in mb.items():
            u, v = m[a], bf[k]
            Q00.append(_conv(u[0], v[0], n)); Q11.append(_conv(u[1], v[1], n))
            Q01.append(add(_conv(u[0], v[1], n), _conv(u[1], v[0], n)))
    def A2_from(P00, P01, P11, n):
        c2 = scal(Fr(-2), add(_conv(y, P01, n), _conv(x, P11, n)))
        c3 = scal(Fr(-2), add(_conv(y, P00, n), _conv(x, P01, n)))
        return [I(0), I(0), c2, c3]
    def A3_from(u, P00, P01, P11, n):
        c2 = scal(Fr(-2), add(_conv(u[0], P11, n), _conv(u[1], P01, n)))
        c3 = scal(Fr(-2), add(_conv(u[0], P01, n), _conv(u[1], P00, n)))
        return [I(0), I(0), c2, c3]
    for n in range(p):
        extend_products(n)
        for a in range(3):
            nxt = div_int_vec(_Av_coeff(Ak, m[a], n), n + 1)
            for i in range(4): m[a][i].append(nxt[i])
        for k, pr in enumerate(PAIRS):
            force = A2_from(*mm[pr], n)
            for store in (bz, bf):
                nxt = div_int_vec(vadd(_Av_coeff(Ak, store[k], n), force), n + 1)
                for i in range(4): store[k][i].append(nxt[i])
        for k, (i1, i2, i3) in enumerate(TRIPS):
            rhs = _Av_coeff(Ak, c[k], n)
            rhs = vadd(rhs, A2_from(*mb[(i1, PAIR_IDX[(i2, i3)])], n))
            rhs = vadd(rhs, A2_from(*mb[(i2, PAIR_IDX[(i1, i3)])], n))
            rhs = vadd(rhs, A2_from(*mb[(i3, PAIR_IDX[(i1, i2)])], n))
            rhs = vadd(rhs, A3_from(m[i1], *mm[(i2, i3)], n))
            nxt = div_int_vec(rhs, n + 1)
            for i in range(4): c[k][i].append(nxt[i])
    return m, bz, bf, c

def div_int_vec(v, n): return [div_int(a, n) for a in v]

def _horner_vec(vser, h, p):
    return [cc.horner(vser[i], h, p) for i in range(4)]

def lemma_dir_step(B, Z, hI, p, m_k, b_k):
    """One step of the directional systems.  Returns (Delta_b[5], Delta_c[4]) where
       Delta_b = b(h) with b(0) = 0  (forcing from m),       and
       Delta_c = c(h) with c(0) = 0  (forcing from m and the FULL b(t), b(0) = b_k).
    Enclosures: Taylor to order p on the box series + Lagrange remainder on the Z series (LAG)."""
    serB = cc.taylor_coeffs(B, p)
    serZ = cc.taylor_coeffs(Z, p + 1)
    hp1 = ipow(hI, p + 1)
    _, bzB, _, cB = dir_series(serB, p, m_k, b_k)          # Taylor part on the box
    _, bzZ, _, cZ = dir_series(serZ, p + 1, m_k, b_k)      # Lagrange remainder on Z
    Db = [vadd(_horner_vec(bzB[k], hI, p), [mul(hp1, bzZ[k][i][p + 1]) for i in range(4)]) for k in range(5)]
    Dc = [vadd(_horner_vec(cB[k], hI, p), [mul(hp1, cZ[k][i][p + 1]) for i in range(4)]) for k in range(4)]
    return Db, Dc

class DirState:
    """FRAME2: b = A beta, c = A gamma in the current Lohner frame."""
    def __init__(self, m0):
        self.m0 = m0                                   # 3 initial interval vectors (frame A_0 = I)
        self.beta = [[I(0)] * 4 for _ in range(5)]
        self.gamma = [[I(0)] * 4 for _ in range(4)]
    def vectors(self, A, R):
        m_k = [mvec(mmul(A, R), v) for v in self.m0]
        b_k = [mvec(A, be) for be in self.beta]
        c_k = [mvec(A, ga) for ga in self.gamma]
        return m_k, b_k, c_k
    def advance(self, T, Ainv_next, Db, Dc):
        self.beta = [vadd(mvec(T, be), mvec(Ainv_next, d)) for be, d in zip(self.beta, Db)]
        self.gamma = [vadd(mvec(T, ga), mvec(Ainv_next, d)) for ga, d in zip(self.gamma, Dc)]

def lemma_dir_time_interval(Z, delta, Jt, m0, b, c):
    """Enclose m, b, c for t in [T_lo, T_hi]: m(t) in Jt m0 (Jt from lemma_rough_time_interval);
       b(t), c(t) by Picard with the forcing evaluated on the enclosures."""
    D = (Fr(0), delta)
    x, y = Z[0], Z[1]
    AZ = [[I(0), I(0), I(1), I(0)], [I(0), I(0), I(0), I(1)],
          [neg(mul(y, y)), scal(Fr(-2), mul(x, y)), I(0), I(0)],
          [scal(Fr(-2), mul(x, y)), neg(mul(x, x)), I(0), I(0)]]
    mt = [mvec(Jt, v) for v in m0]
    def A2v(u, v):
        t01 = add(mul(u[0], v[1]), mul(u[1], v[0]))
        return [I(0), I(0), scal(Fr(-2), add(mul(y, t01), mul(x, mul(u[1], v[1])))),
                scal(Fr(-2), add(mul(y, mul(u[0], v[0])), mul(x, t01)))]
    def A3v(u, v, w):
        return [I(0), I(0),
                scal(Fr(-2), add(add(mul(u[0], mul(v[1], w[1])), mul(u[1], mul(v[0], w[1]))), mul(u[1], mul(v[1], w[0])))),
                scal(Fr(-2), add(add(mul(u[0], mul(v[0], w[1])), mul(u[0], mul(v[1], w[0]))), mul(u[1], mul(v[0], w[0]))))]
    def picard(v0, force):
        vt = v0
        for _ in range(60):
            vn = vadd(v0, [mul(D, a) for a in vadd(mvec(AZ, vt), force)])
            if all(inside(vn[i], vt[i]) for i in range(4)): return vn
            vt = [vt[i] if inside(vn[i], vt[i]) else add(hull(vt[i], vn[i]), (-(width(vn[i]) / 10 + cc.GRID), width(vn[i]) / 10 + cc.GRID)) for i in range(4)]
        raise CheckFailed("Picard for directional variations over the final interval did not close")
    bt = [picard(b[k], A2v(mt[a], mt[bq])) for k, (a, bq) in enumerate(PAIRS)]
    ct = []
    for k, (i1, i2, i3) in enumerate(TRIPS):
        force = vadd(vadd(A2v(mt[i1], bt[PAIR_IDX[(i2, i3)]]), A2v(mt[i2], bt[PAIR_IDX[(i1, i3)]])),
                     vadd(A2v(mt[i3], bt[PAIR_IDX[(i1, i2)]]), A3v(mt[i1], mt[i2], mt[i3])))
        ct.append(picard(c[k], force))
    return mt, bt, ct

# ------------------------------------------------------------- expression trees (TWIST formula)
def eval_tree(t, env):
    if t[0] == 'q':
        n, d = t[1].split('/'); return I(Fr(int(n), int(d)))
    if t[0] == 's': return env[t[1]]
    if t[0] == '+':
        s = I(0)
        for a in t[1:]: s = add(s, eval_tree(a, env))
        return s
    if t[0] == '*':
        s = I(1)
        for a in t[1:]: s = mul(s, eval_tree(a, env))
        return s
    if t[0] == '^':
        base = eval_tree(t[1], env); e = t[2]
        if e >= 0: return ipow(base, e)
        return ipow(cc._div(I(1), base), -e)
    raise CheckFailed(f"bad tree node {t[0]}")

# ------------------------------------------------------------- complex intervals and bivariate polynomials
def cadd(a, b): return (add(a[0], b[0]), add(a[1], b[1]))
def csub(a, b): return (sub(a[0], b[0]), sub(a[1], b[1]))
def cmul(a, b): return (sub(mul(a[0], b[0]), mul(a[1], b[1])), add(mul(a[0], b[1]), mul(a[1], b[0])))
def cconj(a):   return (a[0], neg(a[1]))
def cscal(s, a): return (scal(s, a[0]), scal(s, a[1]))
def cabs2(a):   return add(mul(a[0], a[0]), mul(a[1], a[1]))
def cdiv(a, b):
    d = cabs2(b)
    if not d[0] > 0: raise CheckFailed("complex division: denominator interval contains 0 (resonance?)")
    inv = cc.rnd(1 / d[1], 1 / d[0])
    n = cmul(a, cconj(b))
    return (mul(n[0], inv), mul(n[1], inv))
def _C0(): return (I(0), I(0))
def _C1(): return (I(1), I(0))

class CPoly:
    """bivariate polynomial in (s, t) with complex-interval coefficients, degree <= 3."""
    def __init__(self, d=None): self.d = dict(d or {})
    @staticmethod
    def const(c): return CPoly({(0, 0): c})
    @staticmethod
    def var(i): return CPoly({(1, 0) if i == 0 else (0, 1): _C1()})
    def __add__(self, o):
        d = dict(self.d)
        for m, c in o.d.items(): d[m] = cadd(d.get(m, _C0()), c)
        return CPoly(d)
    def __mul__(self, o):
        d = {}
        for (a1, b1), c1 in self.d.items():
            for (a2, b2), c2 in o.d.items():
                if a1 + a2 + b1 + b2 <= 3:
                    m = (a1 + a2, b1 + b2); d[m] = cadd(d.get(m, _C0()), cmul(c1, c2))
        return CPoly(d)
    def scale(self, c): return CPoly({m: cmul(c, v) for m, v in self.d.items()})
    def __pow__(self, n):
        r = CPoly.const(_C1())
        for _ in range(n): r = r * self
        return r

def lemma_twist(formula, env, DP_check=None):
    """TWIST: from the return-map coefficients (interval) to tau.  Returns (tau, lambda, costheta)."""
    def coeffs(key):
        return {tuple(int(v) for v in k.split(',')): eval_tree(t, env) for k, t in formula[key].items()}
    a1, a2 = coeffs("P1"), coeffs("P2")
    L = [[a1[(1, 0)], a1[(0, 1)]], [a2[(1, 0)], a2[(0, 1)]]]
    ct = scal(Fr(1, 2), add(L[0][0], L[1][1]))                 # cos theta
    one_minus = sub(I(1), mul(ct, ct))
    if not one_minus[0] > 0: raise CheckFailed("return map not elliptic (cos^2 theta >= 1)")
    sn = isqrt(one_minus)                                       # sin theta > 0
    # non-resonance q <= 4:  cos theta not in {1, -1, -1/2, 0}
    for bad in (Fr(1), Fr(-1), Fr(-1, 2), Fr(0)):
        if ct[0] <= bad <= ct[1]: raise CheckFailed(f"resonance: cos theta interval contains {bad}")
    a, b = L[0][0], L[0][1]
    u = [b, sub(ct, a)]; w = [I(0), sn]
    area = mul(b, sn)
    if area[1] < 0: w = [I(0), neg(sn)]; area = neg(area)
    if not area[0] > 0: raise CheckFailed("canonical frame: area interval contains 0")
    sa = isqrt(area); inv_sa = cc._div(I(1), sa)
    S = [[mul(u[0], inv_sa), mul(w[0], inv_sa)], [mul(u[1], inv_sa), mul(w[1], inv_sa)]]
    Sinv = [[mul(w[1], inv_sa), neg(mul(w[0], inv_sa))], [neg(mul(u[1], inv_sa)), mul(u[0], inv_sa)]]
    # zeta = s, zetabar = t ;  zr = (s + t)/2,  zi = -i (s - t)/2  ;  xi = S (zr, zi)
    s, t = CPoly.var(0), CPoly.var(1)
    zr = (s + t).scale((I(Fr(1, 2)), I(0)))
    zi = (s + t.scale((I(-1), I(0)))).scale((I(0), I(Fr(-1, 2))))
    def real_times(x, P): return P.scale((x, I(0)))
    xi1 = real_times(S[0][0], zr) + real_times(S[0][1], zi)
    xi2 = real_times(S[1][0], zr) + real_times(S[1][1], zi)
    pw1 = [CPoly.const(_C1()), xi1, xi1 * xi1, xi1 * xi1 * xi1]
    pw2 = [CPoly.const(_C1()), xi2, xi2 * xi2, xi2 * xi2 * xi2]
    def poly_of(coef):
        P = CPoly()
        for (p1, p2), cv in coef.items():
            P = P + (pw1[p1] * pw2[p2]).scale((cv, I(0)))
        return P
    P1, P2 = poly_of(a1), poly_of(a2)
    zeta_p = real_times(Sinv[0][0], P1) + real_times(Sinv[0][1], P2) \
             + (real_times(Sinv[1][0], P1) + real_times(Sinv[1][1], P2)).scale((I(0), I(1)))
    G = zeta_p.d
    lam = G[(1, 0)]
    from math import factorial
    def g(j, k): c = G.get((j, k), _C0()); return cscal(Fr(factorial(j) * factorial(k)), c)
    lamb = cconj(lam)
    lam2 = cmul(lam, lam)
    t1 = cdiv(cmul(cmul(g(2, 0), g(1, 1)), cadd(csub(lamb, (I(3), I(0))), cscal(Fr(2), lam))),
              cscal(Fr(2), cmul(csub(lam2, lam), csub(lamb, _C1()))))
    t2 = cdiv((cabs2(g(1, 1)), I(0)), csub(_C1(), lamb))
    t3 = cdiv((cabs2(g(0, 2)), I(0)), cscal(Fr(2), csub(lam2, lamb)))
    t4 = cscal(Fr(1, 2), g(2, 1))
    c1 = cadd(cadd(t1, t2), cadd(t3, t4))
    lc = cmul(lamb, c1)
    return lc[1], lam, ct, lc[0], G.get((0, 1), _C0())

def load_formula(path):
    return json.load(open(path))
