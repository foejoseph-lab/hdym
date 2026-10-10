#!/usr/bin/env python3
"""
cert_lemmas.py -- the trusted lemma layer of the periodic-orbit certificate checker
(split out of cert_check.py on 2026-10-10; cert_check.py is now the driver).

Only Python integers and fractions.Fraction; no floats enter any computation.  Every
function named lemma_* (and the recurrences TAY, VAR it relies on) is the shadow of one
lemma a proof assistant would have to state.  The rounding grid lives in this module and
nowhere else, so cert_check.py and cert_twist.py share it by a plain import.

Lemma list (function: statement):
  rnd / add / mul ...       outward rounding to the grid 2^-K preserves inclusion
  taylor_coeffs (TAY)       interval Taylor coefficients of z' = f(z) enclose the true ones
  variational_coeffs (VAR)  same for X' = A(z) X, X(0) = I  (linear in X(0): coefficients
                            for X(0) = X0 are those for X(0) = I times X0)
  lemma_rough (PIC)         B + [0,h] f(Z) subset Z  =>  solutions from B stay in Z on [0,h]
  lemma_var_rough (PICV)    J0 + [0,d] A(Z) Xt subset Xt  =>  X(t) in Xt on [0,d]
  lemma_taylor_point (LAG)  phi_h(c) enclosure: Taylor part + Lagrange remainder on Z
  lemma_jacobian (MVT)      D phi_h over B: Taylor part on B + h^{p+1} X_{p+1}(Z) Xt
  lemma_inverse (NEU)       Neumann-series enclosure of A^{-1}
  lemma_step                one Lohner step: phi_h(c + A r) subset c' + A' r'
  run_chain (FRAME)         accumulated Jacobian in the moving frame
  lemma_rough_time_interval enclosure of z(t), X(t) over the final interval [T_lo, T_hi]
  lemma_krawczyk (KRA)      unique zero of F in X
  lemma_trace_stability (T), lemma_symplectic_4x4 (SYMP)
"""
import sys, time
from fractions import Fraction as Fr

# ============================================================== dyadic interval arithmetic
# An interval is a tuple (lo, hi) of Fractions with lo <= hi.  After every operation the
# result is rounded OUTWARD to the grid 2^-K.  Rounding outward is sound: the true result
# is contained in the rounded interval.  [Lemma: outward rounding preserves inclusion.]
K = 200
GRID = None

def set_grid(k):
    global K, GRID
    K = k; GRID = Fr(1, 2 ** k)

def floor_grid(x):  return Fr((x / GRID).__floor__()) * GRID
def ceil_grid(x):   return Fr(-((-x / GRID).__floor__())) * GRID
def rnd(lo, hi):
    if not lo <= hi:                       # an explicit check, not assert: python -O must not drop it
        raise CheckFailed(f"interval with lo > hi: ({lo}, {hi})")
    return (floor_grid(lo), ceil_grid(hi))

def I(x):            x = Fr(x); return (x, x)
def add(a, b):       return rnd(a[0] + b[0], a[1] + b[1])
def sub(a, b):       return rnd(a[0] - b[1], a[1] - b[0])
def neg(a):          return (-a[1], -a[0])
def mul(a, b):
    ps = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return rnd(min(ps), max(ps))
def scal(s, a):      # s exact Fraction
    return rnd(min(s * a[0], s * a[1]), max(s * a[0], s * a[1]))
def div_int(a, n):   return rnd(a[0] / n, a[1] / n)
def hull(a, b):      return (min(a[0], b[0]), max(a[1], b[1]))
def inside(a, b):    return b[0] <= a[0] and a[1] <= b[1]          # a subset of b
def strictly_inside(a, b): return b[0] < a[0] and a[1] < b[1]
def width(a):        return a[1] - a[0]
def absmax(a):       return max(abs(a[0]), abs(a[1]))
def ipow(a, n):
    r = I(1)
    for _ in range(n): r = mul(r, a)
    return r
def isqrt(a):
    """Enclosure of sqrt on a >= 0 interval, via exact rational bounds."""
    if not a[0] >= 0:
        raise CheckFailed("isqrt: interval not >= 0")
    def lo_sqrt(x):   # largest grid point whose square <= x
        from math import isqrt as _isq
        n = x / GRID / GRID   # x = n * GRID^2 ; want floor(sqrt(n)) * GRID
        return Fr(_isq(int(n))) * GRID
    def hi_sqrt(x):
        from math import isqrt as _isq
        n = x / GRID / GRID
        m = _isq(int(n))
        if Fr(m) * m < n: m += 1
        return Fr(m) * GRID
    return (lo_sqrt(a[0]), hi_sqrt(a[1]))

# vectors / matrices of intervals
def vadd(u, v):      return [add(a, b) for a, b in zip(u, v)]
def vsub(u, v):      return [sub(a, b) for a, b in zip(u, v)]
def mvec(M, v):
    out = []
    for row in M:
        s = I(0)
        for m, x in zip(row, v): s = add(s, mul(m, x))
        out.append(s)
    return out
def mmul(A, B):
    n, m, k = len(A), len(B[0]), len(B)
    return [[_dot([A[i][l] for l in range(k)], [B[l][j] for l in range(k)]) for j in range(m)] for i in range(n)]
def _dot(u, v):
    s = I(0)
    for a, b in zip(u, v): s = add(s, mul(a, b))
    return s
def madd(A, B):      return [[add(a, b) for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]
def mscal(s, A):     return [[scal(s, a) for a in r] for r in A]
def mscal_iv(s, A):  return [[mul(s, a) for a in r] for r in A]
def eye(n):          return [[I(1 if i == j else 0) for j in range(n)] for i in range(n)]
def mat_inf_norm_hi(A):   # upper bound on ||A||_inf
    return max(sum(absmax(a) for a in row) for row in A)

# ============================================================== the vector field
# z = (x, y, px, py);   f(z) = (px, py, -x y^2, -x^2 y)
def f(z):
    x, y, px, py = z
    return [px, py, neg(mul(x, mul(y, y))), neg(mul(mul(x, x), y))]

# Lemma TAY (Taylor coefficients of a polynomial ODE):
#   if z(t) solves z' = f(z), z(0) = z0, then its Taylor coefficients z_k satisfy the
#   recurrences below (Cauchy products), and evaluating them on an interval box
#   containing z0 encloses the true coefficients.
def taylor_coeffs(z0, p):
    x = [z0[0]]; y = [z0[1]]; px = [z0[2]]; py = [z0[3]]
    yy = []; xx = []; xyy = []; xxy = []
    for k in range(p):
        yy.append(_conv(y, y, k)); xx.append(_conv(x, x, k))
        xyy.append(_conv(x, yy, k)); xxy.append(_conv(xx, y, k))
        x.append(div_int(px[k], k + 1)); y.append(div_int(py[k], k + 1))
        px.append(div_int(neg(xyy[k]), k + 1)); py.append(div_int(neg(xxy[k]), k + 1))
    return x, y, px, py, xx, yy

def _conv(a, b, k):
    s = I(0)
    for i in range(k + 1): s = add(s, mul(a[i], b[k - i]))
    return s

# Lemma VAR (variational equation): X' = A(z(t)) X, X(0) = I, with
#   A = [[0,0,1,0],[0,0,0,1],[-y^2 - m, -2xy, 0, 0],[-2xy, -x^2 - m, 0, 0]],
#   m = 0 for the base flow, m = 1 for the ghost (w'' = -(1 + Hess V) w).
#   Its Taylor coefficients follow the Cauchy-product recurrence X_{k+1} = (A X)_k / (k+1).
def variational_coeffs(series, p, m):
    x, y, px, py, xx, yy = series
    mI = I(m)
    Ak = []
    for k in range(p + 1):
        xy = _conv(x, y, k)
        d = I(1 if k == 0 else 0)
        a20 = neg(yy[k]) if k < len(yy) else I(0)
        a31 = neg(xx[k]) if k < len(xx) else I(0)
        if k == 0:
            a20 = sub(a20, mI); a31 = sub(a31, mI)
        Ak.append([[I(0), I(0), d, I(0)], [I(0), I(0), I(0), d],
                   [a20, scal(Fr(-2), xy), I(0), I(0)], [scal(Fr(-2), xy), a31, I(0), I(0)]])
    X = [eye(4)]
    for k in range(p):
        S = [[I(0)] * 4 for _ in range(4)]
        for i in range(k + 1): S = madd(S, mmul(Ak[i], X[k - i]))
        X.append(mscal(Fr(1, k + 1), S))
    return X

def horner(coeffs, h, p):
    if isinstance(coeffs[0], list) and isinstance(coeffs[0][0], list):   # matrix series
        acc = [[I(0)] * 4 for _ in range(4)]
        for k in range(p, -1, -1): acc = madd(mscal_iv(h, acc), coeffs[k])
        return acc
    acc = I(0)
    for k in range(p, -1, -1): acc = add(mul(acc, h), coeffs[k])
    return acc

# ============================================================== lemma-level checks
class CheckFailed(Exception): pass

def lemma_rough(B, h, Z):
    """Lemma PIC: if  B + [0,h] f(Z) subset Z  then every solution starting in B stays in Z
    on [0,h].  (Picard operator maps the ball into itself.)  We verify the hypothesis."""
    H = (Fr(0), h)
    fZ = f(Z)
    for i in range(4):
        img = add(B[i], mul(H, fZ[i]))
        if not inside(img, Z[i]):
            raise CheckFailed(f"rough enclosure: component {i} image {img} not inside Z {Z[i]}")

def lemma_taylor_point(c, h, p, Z):
    """Lemma LAG: phi_h(c) in sum_{k<=p} c_k h^k + h^{p+1} f^{[p+1]}(Z)  (Lagrange remainder,
    the (p+1)-th coefficient evaluated on an enclosure of the trajectory over [0,h])."""
    ser_c = taylor_coeffs(c, p)
    ser_Z = taylor_coeffs(Z, p + 1)
    hp1 = ipow(h, p + 1)
    return [add(horner(ser_c[i], h, p), mul(hp1, ser_Z[i][p + 1])) for i in range(4)]

def lemma_jacobian(B, h, p, Z, m):
    """Lemma MVT: for z in the convex box B,  phi_h(z) - phi_h(c) in J (z - c)  with J an
    enclosure of D phi_h over B  (B must contain c; lemma_step ensures it), where
        J = sum_{k<=p} X_k(B) h^k  +  h^{p+1} X_{p+1}(Z) Xt.
    D phi_t = X(t), X' = A(z(t)) X, X(0) = I.  Lagrange remainder of X at t = h, entrywise:
    X^{(p+1)}(xi)/(p+1)! is the (p+1)-th Taylor coefficient of X about xi, which by linearity
    (VAR) is X_{p+1}[z(xi), X(0)=I] * X(xi), with z(xi) in Z (PIC) and X(xi) in Xt (PICV).
    Status (2026-10-10): the factor Xt was missing (the remainder used X(xi) = I); found by
    the independent review of 10 Oct."""
    XB = variational_coeffs(taylor_coeffs(B, p), p, m)
    XZ = variational_coeffs(taylor_coeffs(Z, p + 1), p + 1, m)
    Xt = lemma_var_rough(Z, h[1], eye(4), m)
    hp1 = ipow(h, p + 1)
    return madd(horner(XB, h, p), mscal_iv(hp1, mmul(XZ[p + 1], Xt)))

def lemma_inverse(A, B):
    """Lemma NEU: if ||I - B A||_inf <= e < 1 then A^{-1} exists and
       A^{-1} in B + [-d, d]^{n x n},  d = ||B||_inf * e / (1 - e).   (Neumann series.)"""
    n = len(A)
    R = madd(eye(n), mscal(Fr(-1), mmul(B, A)))
    e = mat_inf_norm_hi(R)
    if not e < 1:
        raise CheckFailed(f"inverse: ||I - B A|| = {float(e):.3e} >= 1")
    d = mat_inf_norm_hi(B) * e / (1 - e)
    return [[add(B[i][j], (-d, d)) for j in range(n)] for i in range(n)]

def set_box(c, A, r):
    """Enclosing box of  c + A r ."""
    return vadd(c, mvec(A, r))

def lemma_step(step, nxt, h, p, m_list):
    """One Lohner step.  Verifies:  phi_h( c + A r )  subset  c' + A' r'.
    Returns the Jacobian enclosure(s) over the set (one per m in m_list)."""
    c, A, r, Z = step["c"], step["A"], step["r"], step["Z"]
    # B must contain the set c + A r AND the centre c: LAG encloses phi_h(c) with Z, and MVT
    # integrates D phi_h along the segment from c to z.  r comes from the certificate and need
    # not contain 0, so take the hull with c.  Status (2026-10-10): added after the review of
    # 10 Oct (before, B = set_box(c, A, r) only).
    B = [hull(b, ci) for b, ci in zip(set_box(c, A, r), c)]
    lemma_rough(B, h, Z)                       # h: exact Fraction
    hI = I(h)                                  # as a point interval for series evaluation
    phi_c = lemma_taylor_point(c, hI, p, Z)
    Js = [lemma_jacobian(B, hI, p, Z, m) for m in m_list]
    J = Js[0]                                            # base-flow Jacobian drives the set
    # image of the set:  phi_c + J A r  must lie in  c' + A' r',  i.e.  A'^{-1}(phi_c - c') + (A'^{-1} J A) r  subset r'.
    # The grouping matters: (A'^{-1} J A) r is evaluated as a product of matrices first, then
    # applied to r.  Interval arithmetic is not distributive, and this grouping is the tight one
    # (it is the point of the Lohner frame).  Both groupings are sound enclosures.
    Ainv = lemma_inverse(nxt["A"], nxt["Binv"])
    # frame transfer matrices, one per variational system: m=0 uses the set frame A -> A';
    # m>0 uses its own frames G_m -> G_m' (certificate fields G/Ginv), else falls back to A.
    Ts = [mmul(mmul(Ainv, Js[0]), A)]
    for mi in range(1, len(Js)):
        if "G" in step and "G" in nxt:
            Ginv = lemma_inverse(nxt["G"][mi - 1], nxt["Ginv"][mi - 1])
            Ts.append(mmul(mmul(Ginv, Js[mi]), step["G"][mi - 1]))
        else:
            Ts.append(mmul(mmul(Ainv, Js[mi]), A))
    w = vadd(mvec(Ainv, vsub(phi_c, nxt["c"])), mvec(Ts[0], r))
    for i in range(4):
        if not inside(w[i], nxt["r"][i]):
            raise CheckFailed(f"set update: component {i}: w=({float(w[i][0])!r}, {float(w[i][1])!r}) not inside r'=({float(nxt['r'][i][0])!r}, {float(nxt['r'][i][1])!r})")
    lemma_step.extras = dict(B=B, Z=Z, Ainv=Ainv, hI=hI)       # for the directional variations (cert_twist)
    return Ts

def lemma_var_rough(Z, delta, J, m):
    """Lemma PICV: if z(t) in Z for t in [0, delta] and  J + [0,delta] A(Z) Jt  subset Jt,
    then the solution of X' = A(z(t)) X, X(0) in J, satisfies X(t) in Jt, hence in
    J + [0,delta] A(Z) Jt, for all t in [0, delta].  Returns the latter (tighter) enclosure."""
    D = (Fr(0), delta)
    x, y = Z[0], Z[1]
    mI = I(m)
    AZ = [[I(0), I(0), I(1), I(0)], [I(0), I(0), I(0), I(1)],
          [sub(neg(mul(y, y)), mI), scal(Fr(-2), mul(x, y)), I(0), I(0)],
          [scal(Fr(-2), mul(x, y)), sub(neg(mul(x, x)), mI), I(0), I(0)]]
    # Picard for J: iterate Jt <- J + D * AZ * Jt, inflating until self-mapping
    Jt = J
    for _ in range(60):
        Jn = madd(J, mscal_iv(D, mmul(AZ, Jt)))
        if all(inside(Jn[i][j], Jt[i][j]) for i in range(4) for j in range(4)):
            return Jn
        # inflate only the entries that failed
        Jt = [[Jt[i][j] if inside(Jn[i][j], Jt[i][j])
               else add(hull(Jt[i][j], Jn[i][j]), (-(width(Jn[i][j]) / 10 + GRID), width(Jn[i][j]) / 10 + GRID))
               for j in range(4)] for i in range(4)]
    raise CheckFailed("Picard enclosure of the variational matrix did not close")

def lemma_rough_time_interval(B, delta, Z, J, m):
    """Enclose z(t) and J(t) for t in [T_lo, T_hi] = [0, delta] after the chain:
       z(t) in Z (Lemma PIC);  J(t) in J + [0,delta] A(Z) Jt  with Jt a Picard enclosure (PICV)."""
    lemma_rough(B, delta, Z)
    return lemma_var_rough(Z, delta, J, m), f(Z)

def lemma_krawczyk(xhat, X, Fc, DF, Y):
    """Lemma KRA: K(X) = xhat - Y F(xhat) + (I - Y DF(X))(X - xhat).
       If K(X) is strictly inside X then F has exactly one zero in X."""
    n = 3
    YF = mvec(Y, Fc)
    M = madd(eye(n), mscal(Fr(-1), mmul(Y, DF)))
    d = [sub(X[i], xhat[i]) for i in range(n)]
    Md = mvec(M, d)
    Kv = [add(sub(xhat[i], YF[i]), Md[i]) for i in range(n)]
    ok = all(strictly_inside(Kv[i], X[i]) for i in range(n))
    return ok, Kv

def lemma_trace_stability(trM):
    """Lemma T: a 4x4 real symplectic monodromy of a periodic orbit of an autonomous
    Hamiltonian flow has eigenvalue 1 with multiplicity >= 2; the other two have product 1
    and sum  tr M - 2.  If  tr M - 2 in (-2, 2)  they are a conjugate pair on the unit
    circle (elliptic); if |tr M - 2| > 2 they are real, one of modulus > 1 (hyperbolic)."""
    s = sub(trM, I(2))
    if s[0] > -2 and s[1] < 2:  return "ELLIPTIC (stable)", s
    if s[1] < -2 or s[0] > 2:   return "HYPERBOLIC (unstable)", s
    return "UNDECIDED (enclosure straddles |tr-2| = 2)", s

def lemma_symplectic_4x4(M):
    """Lemma SYMP: for a real symplectic 4x4 matrix with a = tr M, b = (a^2 - tr M^2)/2, the
    eigenvalues come in reciprocal pairs and s = lam + 1/lam satisfies s^2 - a s + (b - 2) = 0.
    Both roots real in [-2, 2]  <=> all eigenvalues on the unit circle (stable).
    A real root with |s| > 2, or complex roots (disc < 0)  =>  an eigenvalue off the circle."""
    a = _trace(M); t2 = _trace(mmul(M, M))
    b = scal(Fr(1, 2), sub(mul(a, a), t2))
    disc = sub(mul(a, a), scal(Fr(4), sub(b, I(2))))
    if disc[1] < 0:
        return "UNSTABLE (complex quadruple)", a, disc
    if disc[0] < 0:
        return "UNDECIDED (discriminant straddles 0)", a, disc
    sq = isqrt(disc)
    sp = scal(Fr(1, 2), add(a, sq)); sm = scal(Fr(1, 2), sub(a, sq))
    if sp[0] > 2 or sm[1] < -2:
        return "UNSTABLE (real multiplier off the unit circle)", a, disc
    if sp[1] <= 2 and sm[0] >= -2:
        return "STABLE (all multipliers on the unit circle)", a, disc
    return "UNDECIDED (a root enclosure straddles |s| = 2)", a, disc


# ============================================================== FRAME: one chain
def run_chain(ch, p, m_list, label, verbose, initial_box, dirstate=None):
    h, steps = ch["h"], ch["steps"]
    # steps[0] must contain the initial data
    B0 = set_box(steps[0]["c"], steps[0]["A"], steps[0]["r"])
    for i in range(4):
        if not inside(initial_box[i], B0[i]):
            raise CheckFailed(f"{label}: initial data component {i} not inside steps[0]")
    # Lemma FRAME: accumulated Jacobian  M_k = A_k R_k  with  R_{k+1} = (A_{k+1}^{-1} J_k A_k) R_k,
    # A_0 = I.  Keeps the product in the moving frame so it does not wrap.  (Sound because each
    # factor is an enclosure and matrix products of enclosures enclose the product.)
    Rs = [eye(4) for _ in m_list]
    frames0 = [steps[0]["A"]] + (steps[0]["G"] if "G" in steps[0] else [])
    for Fm in frames0[:len(m_list)]:
        for i in range(4):
            for j in range(4):
                if Fm[i][j] != I(1 if i == j else 0):
                    raise CheckFailed(f"{label}: all initial Jacobian frames must be the identity")
    t0 = time.time()
    n = len(steps) - 1
    for k in range(n):
        Ts = lemma_step(steps[k], steps[k + 1], h, p, m_list)
        if dirstate is not None:                      # VAR2/VAR3 + FRAME2 (cert_twist)
            import cert_twist as ctw
            ex = lemma_step.extras
            m_k, b_k, _ = dirstate.vectors(steps[k]["A"], Rs[0])
            Db, Dc = ctw.lemma_dir_step(ex["B"], ex["Z"], ex["hI"], p, m_k, b_k)
            dirstate.advance(Ts[0], ex["Ainv"], Db, Dc)
        Rs = [mmul(T, R) for T, R in zip(Ts, Rs)]
        if verbose and ((k + 1) % max(1, n // 10) == 0):
            w = max(width(b) for b in set_box(steps[k + 1]["c"], steps[k + 1]["A"], steps[k + 1]["r"]))
            print(f"  [{label}] step {k+1}/{n} ok  box width {float(w):.2e}  ({time.time()-t0:.0f}s)")
            sys.stdout.flush()
    final = steps[-1]
    framesN = [final["A"]] + (final["G"] if "G" in final else [final["A"]] * (len(m_list) - 1))
    Jacc = [mmul(Fm, R) for Fm, R in zip(framesN, Rs)]
    if dirstate is not None:
        run_chain.dir_final = dirstate.vectors(final["A"], Rs[0])
    return set_box(final["c"], final["A"], final["r"]), Jacc, h * n

def _div(a, b):
    if not (b[0] > 0 or b[1] < 0):
        raise CheckFailed("division: denominator interval contains 0")
    inv = rnd(min(1 / b[0], 1 / b[1]), max(1 / b[0], 1 / b[1]))
    return mul(a, inv)
def _trace(M):
    s = I(0)
    for i in range(4): s = add(s, M[i][i])
    return s
