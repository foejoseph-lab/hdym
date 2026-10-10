"""
cap_dr_orbit.py  (v2: certificate generator)

UNTRUSTED layer.  Finds a periodic orbit of Yang-Mills mechanics V = x^2 y^2 / 2 near the
Dahlqvist-Russberg initial data, runs a validated Lohner/Taylor integration to choose
frames, and writes a certificate (cert.json) for cert_check.py, which re-verifies every
enclosure in exact rational arithmetic.  Nothing computed here needs to be correct; a bad
certificate simply fails to check.

Design choices forced by the checker:
  * every frame centre c, frame matrix A and its approximate inverse Binv are doubles
    (exact dyadics), so the checker reproduces them exactly;
  * interval endpoints r, Z are rounded OUTWARD to doubles before emission;
  * step h is a double; chain P runs N steps so T_center = N h exactly (emitted as a
    rational); chain B runs the same N steps to T_lo = T_center and the Krawczyk box in T
    is one-sided, [T_center, T_center + 2 radT]  (Krawczyk only needs x^ in X).

Run:  python3 cap_dr_orbit.py [p=24] [N=800] [rad=1e-9] [out=cert.json]
Then: python3 cert_check.py cert.json --ghost --verbose
"""
import sys, time, json
import numpy as np
from fractions import Fraction
from mpmath import iv, mp
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

iv.prec = 113
E = iv.mpf(1) / 2

# ----------------------------------------------------------------------------- interval helpers
from mpmath.libmp import to_rational
def lo_frac(x):      return Fraction(*to_rational(x._mpi_[0]))      # exact lower endpoint
def hi_frac(x):      return Fraction(*to_rational(x._mpi_[1]))      # exact upper endpoint
def lo_f(x):         return float(lo_frac(x))
def hi_f(x):         return float(hi_frac(x))
def I(x):            return iv.mpf(x)
def hull(a, b):
    from mpmath.libmp import mpf_lt
    lo = a._mpi_[0] if mpf_lt(a._mpi_[0], b._mpi_[0]) else b._mpi_[0]
    hi = b._mpi_[1] if mpf_lt(a._mpi_[1], b._mpi_[1]) else a._mpi_[1]
    return iv.make_mpf((lo, hi))
def width(a):        return hi_f(a) - lo_f(a)
def mid_double(a):   return float((lo_frac(a) + hi_frac(a)) / 2)
def absmax_f(a):     return max(abs(lo_f(a)), abs(hi_f(a)))
def sym(r):          return iv.mpf([-r, r])
def contains_int(a, b):   return lo_frac(b) < lo_frac(a) and hi_frac(a) < hi_frac(b)
def mat_mul(A, B):
    n, m, k = len(A), len(B[0]), len(B)
    return [[sum((A[i][l] * B[l][j] for l in range(k)), I(0)) for j in range(m)] for i in range(n)]
def mat_vec(A, v):   return [sum((A[i][j] * v[j] for j in range(len(v))), I(0)) for i in range(len(A))]
def mat_add(A, B):   return [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]
def mat_scale(s, A): return [[s * a for a in r] for r in A]
def eye(n):          return [[I(1 if i == j else 0) for j in range(n)] for i in range(n)]
def to_iv_mat(M):    return [[I(float(M[i, j])) for j in range(M.shape[1])] for i in range(M.shape[0])]
def np_mid(A):       return np.array([[float(a.mid) for a in r] for r in A])

def rigorous_inverse(Qf):
    Q = to_iv_mat(Qf); QT = to_iv_mat(Qf.T)
    Eres = mat_add(eye(4), mat_scale(I(-1), mat_mul(QT, Q)))
    n = max(sum(absmax_f(e) for e in row) for row in Eres)
    nQT = max(sum(absmax_f(e) for e in row) for row in QT)
    assert n < 1
    b = (nQT * n / (1 - n)) * (1 + 1e-12)
    return [[QT[i][j] + sym(b) for j in range(4)] for i in range(4)]

# outward rounding of an iv interval to doubles (exact comparisons via Fraction)
def out_lo(x):
    lo = lo_frac(x); v = float(lo)
    return v if Fraction(v) <= lo else float(np.nextafter(v, -np.inf))
def out_hi(x):
    hi = hi_frac(x); v = float(hi)
    return v if Fraction(v) >= hi else float(np.nextafter(v, np.inf))
def hx(v):      return float(v).hex()
def hiv(x):     return [hx(out_lo(x)), hx(out_hi(x))]

# ----------------------------------------------------------------------------- vector field, Taylor
def f(z):
    x, y, px, py = z
    return [px, py, -x * y * y, -x * x * y]

def taylor_coeffs(z0, p):
    x = [z0[0]]; y = [z0[1]]; px = [z0[2]]; py = [z0[3]]
    yy = []; xx = []; xyy = []; xxy = []
    for k in range(p):
        yy.append(sum((y[i] * y[k - i] for i in range(k + 1)), I(0)))
        xx.append(sum((x[i] * x[k - i] for i in range(k + 1)), I(0)))
        xyy.append(sum((x[i] * yy[k - i] for i in range(k + 1)), I(0)))
        xxy.append(sum((xx[i] * y[k - i] for i in range(k + 1)), I(0)))
        x.append(px[k] / (k + 1)); y.append(py[k] / (k + 1))
        px.append(-xyy[k] / (k + 1)); py.append(-xxy[k] / (k + 1))
    return x, y, px, py, xx, yy

def variational_coeffs(series, p, m=0):
    x, y, px, py, xx, yy = series
    z0 = I(0); mI = I(m)
    Ak = []
    for k in range(p + 1):
        xy = sum((x[i] * y[k - i] for i in range(k + 1)), I(0))
        a20 = -yy[k] if k < len(yy) else z0
        a31 = -xx[k] if k < len(xx) else z0
        if k == 0: a20 = a20 - mI; a31 = a31 - mI
        d = I(1 if k == 0 else 0)
        Ak.append([[z0, z0, d, z0], [z0, z0, z0, d], [a20, -2 * xy, z0, z0], [-2 * xy, a31, z0, z0]])
    X = [eye(4)]
    for k in range(p):
        S = [[I(0)] * 4 for _ in range(4)]
        for i in range(k + 1): S = mat_add(S, mat_mul(Ak[i], X[k - i]))
        X.append(mat_scale(I(1) / (k + 1), S))
    return X

def horner(coeffs, h, p):
    if isinstance(coeffs[0], list):
        acc = [[I(0)] * 4 for _ in range(4)]
        for k in range(p, -1, -1): acc = mat_add(mat_scale(h, acc), coeffs[k])
        return acc
    acc = I(0)
    for k in range(p, -1, -1): acc = acc * h + coeffs[k]
    return acc

# ----------------------------------------------------------------------------- Lohner
class LohnerSet:
    def __init__(self, c, A, r, Rs=None, Fr=None):
        self.c, self.A, self.r = c, A, r             # c: list of iv (double-valued), A: float 4x4, r: iv
        self.Rs = Rs                                  # list of interval 4x4: Jacobian_m = Fr[m] @ Rs[m]
        self.Fr = Fr                                  # list of float 4x4 frames, Fr[0] is A
    def box(self):
        Ai = to_iv_mat(self.A)
        return [self.c[i] + sum((Ai[i][j] * self.r[j] for j in range(4)), I(0)) for i in range(4)]

def rough_enclosure(B, h, maxit=60):
    """Z with B + [0,h] f(Z) strictly inside Z, by Picard iteration with epsilon-inflation
    of the new iterate (never of the old one, which is what makes it diverge)."""
    H = I([0, h])
    Z = [b + sym(1e-12 * max(1.0, absmax_f(b))) for b in B]
    for _ in range(maxit):
        fZ = f(Z)
        Znew = [B[i] + H * fZ[i] for i in range(4)]
        if all(contains_int(Znew[i], Z[i]) for i in range(4)):
            # emit-safe margin so the checker's outward-rounded Z still satisfies the test
            return [z + sym(1e-13 * max(1.0, absmax_f(z)) + 1e-13 * width(z)) for z in Z]
        # inflate ONLY the components that failed; the others stay fixed so the iteration closes
        Z = [Z[i] if contains_int(Znew[i], Z[i])
             else hull(Z[i], Znew[i]) + sym(0.1 * width(Znew[i]) + 1e-12 * max(1.0, absmax_f(Znew[i])))
             for i in range(4)]
    raise RuntimeError("rough enclosure failed; increase N")

def lohner_step(S, h, p, m_list):
    hI = I(h)
    B = S.box()
    Z = rough_enclosure(B, h)
    ser_c = taylor_coeffs(S.c, p)
    ser_Z = taylor_coeffs(Z, p + 1)
    phi_c = [horner(ser_c[i], hI, p) + ser_Z[i][p + 1] * hI ** (p + 1) for i in range(4)]
    serB = taylor_coeffs(B, p)
    Js = []
    for m in m_list:
        XB = variational_coeffs(serB, p, m)
        XZ = variational_coeffs(ser_Z, p + 1, m)
        Js.append(mat_add(horner(XB, hI, p), mat_scale(hI ** (p + 1), XZ[p + 1])))
    J = Js[0]
    c_new = [I(mid_double(v)) for v in phi_c]              # DOUBLE centres
    JA = mat_mul(J, to_iv_mat(S.A))
    Qf, _ = np.linalg.qr(np_mid(JA))
    Qinv = rigorous_inverse(Qf)
    QJA = mat_mul(Qinv, JA)                                # = A'^{-1} J A  (nearly triangular)
    # accumulated Jacobians, each in ITS OWN moving frame:  M_m = G_m R_m,  R_m' = (G_m'^{-1} J_m G_m) R_m
    Rs_new, Fr_new = [], []
    for mi, (Jm, R) in enumerate(zip(Js, S.Rs)):
        if mi == 0:
            Rs_new.append(mat_mul(QJA, R)); Fr_new.append(Qf)
        else:
            JG = mat_mul(Jm, to_iv_mat(S.Fr[mi]))
            Gf, _ = np.linalg.qr(np_mid(JG))
            Ginv = rigorous_inverse(Gf)
            Rs_new.append(mat_mul(mat_mul(Ginv, JG), R)); Fr_new.append(Gf)
    r_new = mat_vec(QJA, S.r)
    shift = [phi_c[i] - c_new[i] for i in range(4)]
    r_new = [r_new[i] + sum((Qinv[i][j] * shift[j] for j in range(4)), I(0)) for i in range(4)]
    # slack so the checker's own (slightly different) rounding still passes
    # slack must NOT compound: 1e-5 relative over 1000 steps is a factor 1.01 in total
    r_new = [r + sym(1e-4 * width(r) + 1e-27 * max(1.0, absmax_f(r))) for r in r_new]
    return LohnerSet(c_new, Qf, r_new, Rs_new, Fr_new), Js, Z

def step_record(S, Z=None):
    rec = {"c": [hx(float(c.mid)) for c in S.c],
           "A": [[hx(S.A[i, j]) for j in range(4)] for i in range(4)],
           "Binv": [[hx(S.A.T[i, j]) for j in range(4)] for i in range(4)],
           "r": [hiv(r) for r in S.r]}
    if S.Fr is not None and len(S.Fr) > 1:            # extra frames (ghost etc.), with approx inverses
        rec["G"] = [[[hx(G[i, j]) for j in range(4)] for i in range(4)] for G in S.Fr[1:]]
        rec["Ginv"] = [[[hx(G.T[i, j]) for j in range(4)] for i in range(4)] for G in S.Fr[1:]]
    if Z is not None: rec["Z"] = [hiv(z) for z in Z]
    return rec

def run_chain(S, h, N, p, m_list, label):
    steps = []
    S.Rs = [eye(4) for _ in m_list]                  # all frames start at I, so R_0 = I
    S.Fr = [np.eye(4) for _ in m_list]
    t0 = time.time()
    for k in range(N):
        S_new, Js, Z = lohner_step(S, h, p, m_list)
        steps.append(step_record(S, Z))
        S = S_new
        if (k + 1) % max(1, N // 20) == 0:
            w = max(width(b) for b in S.box())
            print(f"  [{label}] step {k+1}/{N}  max box width={w:.2e}  ({time.time()-t0:.0f}s)"); sys.stdout.flush()
    steps.append(step_record(S))
    Jacc = [mat_mul(to_iv_mat(G), R) for G, R in zip(S.Fr, S.Rs)]
    return S, Jacc, steps

def initial_set(xb, pxb):
    """Set containing (x0, 0, px0, sqrt(2E - px0^2)) for x0 in xb, px0 in pxb; double centre."""
    py0 = iv.sqrt(2 * E - pxb * pxb)
    z0 = [xb, I(0), pxb, py0]
    c = [I(mid_double(v)) for v in z0]
    r = [(z0[i] - c[i]) + sym(1e-25) for i in range(4)]       # slack for checker's isqrt grid
    return LohnerSet(c, np.eye(4), r)

def newton_refine(x0, px0, T, E0=0.5):
    def F(v):
        x, px, T = v
        py = np.sqrt(2 * E0 - px ** 2)
        s = solve_ivp(lambda t, s: f(s), (0, T), [x, 0, px, py], rtol=1e-13, atol=1e-15)
        xT, yT, pxT, pyT = s.y[:, -1]
        return [xT - x, pxT - px, yT]
    return fsolve(F, [x0, px0, T], xtol=1e-14)

def main():
    p    = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    N    = int(sys.argv[2]) if len(sys.argv) > 2 else 800
    rad  = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-9
    out  = sys.argv[4] if len(sys.argv) > 4 else "cert.json"
    radT = 10 * rad

    x0, px0, T0 = 3.14640122769753, 0.0017931934191, 31.053330268627366
    print("Newton refine ...")
    xr, pxr, Tr = newton_refine(x0, px0, T0)
    h = float(Tr / N)                                  # a double
    Tc = Fraction(h) * N                               # exact dyadic rational = N h
    xr, pxr = float(xr), float(pxr)
    print(f"  x0={xr!r} px0={pxr!r} T={float(Tc)!r} h={h!r}")

    # chain P: center point
    print("chain P ...")
    SP = initial_set(I(xr), I(pxr))
    SPf, _, stepsP = run_chain(SP, h, N, p, [0], "P")
    zT = SPf.box()
    Fc = [zT[0] - I(xr), zT[2] - I(pxr), zT[1]]
    print("  F(center) widths:", [f"{width(v):.1e}" for v in Fc])

    # chain B: box
    print("chain B ...")
    # the Krawczyk box is defined with DOUBLE endpoints so the certificate carries it exactly
    Xb = [iv.mpf([float(xr - rad), float(xr + rad)]), iv.mpf([float(pxr - rad), float(pxr + rad)])]
    SB = initial_set(Xb[0], Xb[1])
    hB = float((Tc - Fraction(radT)) / N)             # a double; N hB = T_lo (exact rational) < Tc
    Tlo = Fraction(hB) * N
    assert Tlo < Tc
    Thi = 2 * Tc - Tlo                                 # symmetric box around the centre T_c
    deltaT = float(Thi - Tlo)
    SBf, Jaccs, stepsB = run_chain(SB, hB, N, p, [0, 1], "B")
    Bbox = SBf.box()
    Zf = rough_enclosure(Bbox, deltaT)
    # preview Krawczyk (non-rigorous midpoint linear algebra for Y)
    D = I([0, deltaT])
    x, y = Zf[0], Zf[1]
    AZ = [[I(0), I(0), I(1), I(0)], [I(0), I(0), I(0), I(1)], [-y*y, -2*x*y, I(0), I(0)], [-2*x*y, -x*x, I(0), I(0)]]
    J = Jaccs[0]
    Jt = [[a + sym(1e-12 * max(1.0, absmax_f(a))) for a in row] for row in J]
    for _ in range(60):
        Jn = mat_add(J, mat_scale(D, mat_mul(AZ, Jt)))
        if all(contains_int(Jn[i][j], Jt[i][j]) for i in range(4) for j in range(4)):
            Jt = Jn; break
        Jt = [[Jt[i][j] if contains_int(Jn[i][j], Jt[i][j])
               else hull(Jt[i][j], Jn[i][j]) + sym(0.1 * width(Jn[i][j]) + 1e-12 * max(1.0, absmax_f(Jn[i][j])))
               for j in range(4)] for i in range(4)]
    py0 = iv.sqrt(2 * E - Xb[1] * Xb[1])
    dz0 = [[I(1), I(0)], [I(0), I(0)], [I(0), I(1)], [I(0), -Xb[1] / py0]]
    JZ = mat_mul(Jt, dz0); fz = f(Zf)
    DF = [[JZ[r][0], JZ[r][1], fz[r]] for r in [0, 2, 1]]
    DF[0][0] -= 1; DF[1][1] -= 1
    Yf = np.linalg.inv(np_mid(DF))
    trM = sum((Jt[i][i] for i in range(4)), I(0))
    trG = sum((Jaccs[1][i][i] for i in range(4)), I(0))
    print(f"  preview: tr M in [{lo_f(trM):.6f}, {hi_f(trM):.6f}]  HI ~ {(mid_double(trM)-2)/2:.5f};  ghost tr (at T_lo) in [{lo_f(trG):.4f}, {hi_f(trG):.4f}]")

    cert = {"E": hx(0.5), "p": p, "K": 200,
            "box": {"x0": hiv(Xb[0]), "px0": hiv(Xb[1]), "T": [f"{Tlo.numerator}/{Tlo.denominator}", f"{Thi.numerator}/{Thi.denominator}"]},
            "center": {"x0": hx(xr), "px0": hx(pxr), "T": f"{Tc.numerator}/{Tc.denominator}"},
            "Y": [[hx(Yf[i, j]) for j in range(3)] for i in range(3)],
            "chainP": {"h": hx(h), "steps": stepsP},
            "chainB": {"h": hx(hB), "steps": stepsB, "Zfinal": [hiv(z) for z in Zf]}}
    json.dump(cert, open(out, "w"))
    print(f"wrote {out}  ({len(stepsP)} + {len(stepsB)} frames).  Now: python3 cert_check.py {out} --ghost --verbose")

if __name__ == "__main__":
    main()
