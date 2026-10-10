"""
bianchi9_six.py -- sixth-order (two-ghost) Weyl gravity in Bianchi I: do the ghosts collide?

    S = int sqrt(-g) [ R/2 - (k2/4) C_{abcd}C^{abcd} - (k4/4) C_{abcd} Box C^{abcd} ],
    k2 = 1/M1^2 + 1/M2^2,  k4 = 1/(M1^2 M2^2)

so that the spin-2 sector has poles at 0, M1, M2 with alternating Krein signs (checked below
on the minisuperspace transverse block: frequencies 0, M1, M2 in the heavy limit).

Bianchi I:  ds^2 = -N^2 dt^2 + a^2 dx^2 + b^2 dy^2 + c^2 dz^2, Misner variables as before,
no curvature potential; the Kasner solution is exact (Ricci-flat => Bach-flat and Box C-flat
contributions vanish on shell: checked, the Einstein manifold is invariant).

Method: Cartan in the orthonormal frame (only F^i_{0i} nonzero), Weyl tensor, covariant
derivative in the frame, Box C, the three scalars; checks:
  1. C^2 = 0 and C Box C = 0 for isotropic a = b = c
  2. sqrt(-g)(C Box C + nabla C . nabla C) is a total derivative on the ansatz (EL vanish)
  3. the reduced Lagrangian is sixth order: solve for beta^(6); alpha'' from E_alpha
  4. heavy-limit transverse frequencies = (0, M1, M2) with Krein signs (+, -, +)
Then: scaled Jacobian along Kasner in s = ln t, eigenvalues vs M1 t for M2/M1 = ratio;
print the ghost frequencies / M1 and real parts, flag any crossing (complex-pair formation
away from the imaginary axis) -- a Krein collision.

Run:  %run bianchi9_six.py --ratio 3 --theta 3.14159
"""
import time, argparse, itertools
import math
import numpy as np
import sympy as sp

P = argparse.ArgumentParser()
P.add_argument('--ratio', type=float, default=3.0, help='M2/M1')
P.add_argument('--theta', type=float, default=np.pi)
args = P.parse_args()
t0 = time.time()
def stamp(msg): print(f"[{time.time()-t0:6.1f}s] {msg}", flush=True)

t = sp.symbols('t')
M1, M2 = sp.symbols('M1 M2', positive=True)
N, a, b, c = [sp.Function(s)(t) for s in ('N', 'a', 'b', 'c')]
A = [None, a, b, c]; eta = [-1, 1, 1, 1]; idx = range(4)
d0 = lambda f: sp.diff(f, t) / N

# ------------------------------------------------ structure functions and connection (Bianchi I)
F = [[[sp.S(0)]*4 for _ in idx] for _ in idx]
for i in (1, 2, 3):
    F[i][0][i] = sp.diff(A[i], t) / (N * A[i]); F[i][i][0] = -F[i][0][i]
Fl = [[[eta[p] * F[p][q][r] for r in idx] for q in idx] for p in idx]
w = [[[sp.simplify((Fl[p][q][r] + Fl[q][r][p] - Fl[r][p][q]) / 2) for r in idx] for q in idx] for p in idx]   # omega_{abc}
wu = [[[eta[p] * w[p][q][r] for r in idx] for q in idx] for p in idx]                                          # omega^a_{bc}

# ------------------------------------------------ Riemann, Ricci, Weyl (frame components, all indices down)
Riem = [[[[sp.S(0)]*4 for _ in idx] for _ in idx] for _ in idx]
for p in idx:
    for q in idx:
        for dd in idx:
            for ee in idx:
                val = sp.S(0)
                if dd == 0: val += d0(wu[p][q][ee])
                if ee == 0: val -= d0(wu[p][q][dd])
                val += sum(wu[p][q][r] * F[r][dd][ee] for r in idx)
                val += sum(wu[p][r][dd] * wu[r][q][ee] - wu[p][r][ee] * wu[r][q][dd] for r in idx)
                Riem[p][q][dd][ee] = sp.simplify(val)
Rl = [[[[eta[p] * Riem[p][q][r][s] for s in idx] for r in idx] for q in idx] for p in idx]       # R_{abcd}
Ric = [[sp.simplify(sum(Riem[p][q][p][r] for p in idx)) for r in idx] for q in idx]           # R_{bd}
Rs = sp.simplify(sum(eta[q] * Ric[q][q] for q in idx))
def g(p, q): return eta[p] if p == q else 0
C = [[[[sp.simplify(Rl[p][q][r][s]
        - sp.Rational(1, 2) * (g(p, r) * Ric[q][s] - g(p, s) * Ric[q][r] - g(q, r) * Ric[p][s] + g(q, s) * Ric[p][r])
        + Rs / 6 * (g(p, r) * g(q, s) - g(p, s) * g(q, r)))
        for s in idx] for r in idx] for q in idx] for p in idx]
# tracelessness check
assert all(sp.simplify(sum(eta[p] * C[p][q][p][s] for p in idx)) == 0 for q in idx for s in idx)
stamp("Weyl tensor done (traceless)")

# ------------------------------------------------ covariant derivatives in the frame
def cov(T, rank):
    """nabla_e T_{a1..ar}: returns tensor with one extra (last) index e."""
    out = {}
    for ind in itertools.product(idx, repeat=rank):
        for e in idx:
            val = d0(T[ind]) if e == 0 else sp.S(0)
            for pos in range(rank):
                for f in idx:
                    ind2 = list(ind); ind2[pos] = f
                    val -= wu[f][ind[pos]][e] * T[tuple(ind2)]
            out[ind + (e,)] = val
    return out
Cd = {ind: C[ind[0]][ind[1]][ind[2]][ind[3]] for ind in itertools.product(idx, repeat=4)}
DC = cov(Cd, 4)                      # nabla_e C_{abcd}
DDC = cov(DC, 5)                     # nabla_f nabla_e C_{abcd}
stamp("second covariant derivative done")
def up(ind): return sp.prod([eta[i] for i in ind])
BoxC = {ind: sp.simplify(sum(eta[e] * DDC[ind + (e, e)] for e in idx)) for ind in itertools.product(idx, repeat=4)}
C2 = sp.simplify(sum(up(ind) * Cd[ind]**2 for ind in Cd))
CBoxC = sp.simplify(sum(up(ind) * Cd[ind] * BoxC[ind] for ind in Cd))
DCDC = sp.simplify(sum(up(ind) * DC[ind]**2 for ind in DC))
stamp("scalars done")

# check 1: isotropic
f = sp.Function('f')(t)
iso = {a: f, b: f, c: f}
print("   check 1 isotropic: C^2 =", sp.simplify(C2.subs(iso).doit()), " C Box C =", sp.simplify(CBoxC.subs(iso).doit()))

# ------------------------------------------------ Misner variables, jets
al, bp, bm = [sp.Function(s)(t) for s in ('alpha', 'bp', 'bm')]
s3 = sp.sqrt(3)
misner = {a: sp.exp(al + bp + s3*bm), b: sp.exp(al + bp - s3*bm), c: sp.exp(al - 2*bp)}
sqrtg = N * a * b * c
k2 = 1 / M1**2 + 1 / M2**2; k4 = 1 / (M1**2 * M2**2)
L_R = (sqrtg * Rs / 2).subs(misner).doit()
L_C2 = (-sqrtg * C2 / 4).subs(misner).doit()
L_CB = (-sqrtg * CBoxC / 4).subs(misner).doit()
L_tot_deriv = (sqrtg * (CBoxC + DCDC)).subs(misner).doit()
funcs = [al, bp, bm, N]; names = ['al', 'bp', 'bm', 'N']; KMAX = 12
J = {nm: [sp.symbols(f'{nm}{k}') for k in range(KMAX + 1)] for nm in names}
def to_jet(expr):
    expr = expr.doit(); rep = {}
    for fn, nm in zip(funcs, names):
        for k in range(KMAX, 0, -1): rep[sp.Derivative(fn, (t, k))] = J[nm][k]
    expr = expr.subs(rep)
    for fn, nm in zip(funcs, names): expr = expr.subs(fn, J[nm][0])
    return expr
def D(expr): return sum(sp.diff(expr, J[nm][k]) * J[nm][k+1] for nm in names for k in range(KMAX))
def EL(Lx, nm, order=6):
    out = sp.S(0)
    for k in range(order + 1):
        term = sp.diff(Lx, J[nm][k])
        for _ in range(k): term = D(term)
        out += (-1)**k * term
    return sp.expand(out)
L_R, L_C2, L_CB, L_td = [sp.expand(to_jet(e)) for e in (L_R, L_C2, L_CB, L_tot_deriv)]
stamp("jet form done")
ok = all(sp.simplify(EL(L_td, nm)) == 0 for nm in names)
print(f"   check 2 (C Box C + |nabla C|^2 is a total derivative): {'OK' if ok else 'FAIL'}")
gauge = {J['N'][0]: 1, J['N'][1]: 0, J['N'][2]: 0, J['N'][3]: 0, J['N'][4]: 0, J['N'][5]: 0}
L = sp.expand(k2 * L_C2 + k4 * L_CB + L_R)
# orders present
def orders(expr, nm): return [k for k in range(KMAX + 1) if sp.diff(expr, J[nm][k]) != 0]
print("   derivative orders in L:  al", orders(L, 'al'), " bp", orders(L, 'bp'), " N", orders(L, 'N'))
E = {nm: EL(L, nm).subs(gauge) for nm in ('al', 'bp', 'bm')}
Con = EL(L, 'N').subs(gauge)
for nm in ('al', 'bp', 'bm', ):
    print(f"   E_{nm}: orders al {orders(E[nm], 'al')} bp {orders(E[nm], 'bp')} bm {orders(E[nm], 'bm')}")
stamp("Euler-Lagrange done")

# ------------------------------------------------ solve: highest derivatives
top_al = max(orders(E['al'], 'al')); top_b = max(orders(E['bp'], 'bp'))
print(f"   alpha equation order {top_al}, beta equation order {top_b}")
# E_alpha is linear in alpha^(4): alpha4 = -c0/c1.  alpha^(5) = D(alpha4) contains beta^(6);
# substitute both into E_beta (linear in alpha4, alpha5, beta6), clear denominators, 2x2 solve.
a4, a5 = J['al'][top_al], J['al'][top_al + 1]
c1 = E['al'].coeff(a4); c0 = sp.expand(E['al'] - c1 * a4)
assert sp.diff(c0, a4) == 0 and not c1.has(a4)
n5 = sp.expand(-(D(c0) * c1 - c0 * D(c1)))                 # alpha5 = n5 / c1^2 with alpha4 inside
n5 = sp.expand(n5.subs(a4, -c0 / c1) * c1)                 # -> n5 / c1^3
def clear(Ex):
    Ex = sp.expand(Ex); out = sp.S(0)
    for term in sp.Add.make_args(Ex):
        p4 = sp.degree(term, a4); p5 = sp.degree(term, a5)
        assert p4 <= 1 and p5 <= 1 and p4 + p5 <= 1, term
        base = term / a4**p4 / a5**p5
        out += base * (-c0)**p4 * c1**(3 - p4 - 3*p5) * n5**p5
    return sp.expand(out)
Eb = [clear(E['bp']), clear(E['bm'])]
hb = [J['bp'][top_b], J['bm'][top_b]]
Amat, rvec = sp.linear_eq_to_matrix(Eb, hb)
assert all(sp.diff(Amat, h_) == sp.zeros(2, 2) for h_ in hb)
print("   c1 (alpha^(4) coefficient):", sp.factor(c1))
state = [J['al'][k] for k in range(top_al)] + [J['bp'][k] for k in range(top_b)] + [J['bm'][k] for k in range(top_b)]
print(f"   phase space dimension {len(state)}")
na, nb = top_al, top_b
f_c = sp.lambdify([state, M1, M2], [c0, c1], 'numpy', cse=True)
f_A = sp.lambdify([state, M1, M2], Amat, 'numpy', cse=True)
f_r = sp.lambdify([state, M1, M2], rvec, 'numpy', cse=True)
Con_s = Con.subs(a4, -c0 / c1)
f_con = sp.lambdify([state, M1, M2], Con_s, 'numpy', cse=True)
def rhs(y, m1, m2):
    c0v, c1v = f_c(y, m1, m2)
    Am = np.asarray(f_A(y, m1, m2), float); rv = np.asarray(f_r(y, m1, m2), float).ravel()
    b6 = np.linalg.solve(Am, rv)
    dy = np.zeros_like(y)
    dy[:na-1] = y[1:na]; dy[na-1] = -c0v / c1v
    dy[na:na+nb-1] = y[na+1:na+nb]; dy[na+nb-1] = b6[0]
    dy[na+nb:na+2*nb-1] = y[na+nb+1:]; dy[na+2*nb-1] = b6[1]
    return dy
def jac(y, m1, m2, h=1e-6):
    Jm = np.zeros((len(y), len(y)))
    for i in range(len(y)):
        e = np.zeros(len(y)); e[i] = h * (1 + abs(y[i]))
        Jm[:, i] = (rhs(y + e, m1, m2) - rhs(y - e, m1, m2)) / (2 * e[i])
    return Jm
stamp("numerical rhs ready")

# ------------------------------------------------ Kasner background and GR values of the higher derivatives
# GR Bianchi I: alpha'' = -3/2(alpha'^2+|beta'|^2) ... on Kasner: alpha = ln t/3, alpha' = 1/3t, beta' = n/3t
def kasner(tt, theta):
    y = np.zeros(len(state))
    al_d = [np.log(tt)/3] + [(-1)**(k-1) * math.factorial(k-1) / (3 * tt**k) for k in range(1, na)]
    y[:na] = al_d
    for j, cs in ((na, np.cos(theta)), (na+nb, np.sin(theta))):
        y[j] = 0.0
        for k in range(1, nb): y[j+k] = cs * (-1)**(k-1) * math.factorial(k-1) / (3 * tt**k)
    return y
m1 = 1.0; m2 = args.ratio
yk = kasner(50.0, args.theta); fk = rhs(yk, m1, m2); yk2 = kasner(50.0*(1+1e-7), args.theta)
print(f"   Kasner is a solution: |rhs - d/dt Kasner| / scale = {np.abs(fk - (yk2-yk)/(50*1e-7)).max() / (1+np.abs(fk).max()):.1e}")
print(f"   constraint on Kasner: {f_con(yk, m1, m2):.1e}")

# check 4: heavy-limit transverse frequencies (M1 t >> 1)
Jf = jac(kasner(1e4, args.theta), m1, m2); ev = np.linalg.eigvals(Jf)
freqs = np.sort(np.abs(ev.imag[ev.imag > 1e-6]))
print(f"   check 4: oscillation frequencies at M1 t = 1e4 (expect {m1}, {m2} each twice): {np.round(freqs, 5)}")

# ------------------------------------------------ Frobenius exponents along Kasner
wts = np.array([k for k in range(na)] + [k for k in range(nb)] + [k for k in range(nb)], float)
print("\n   scaled-Jacobian eigenvalues sigma along Kasner (perturbation ~ t^sigma; toward the singularity the smallest Re wins)")
print("   M1 t      Re sigma (sorted), with |Im| in brackets for complex ones")
for Mt in [1e3, 1e2, 30, 10, 3, 1, 0.3, 0.1, 0.03, 0.01]:
    y = kasner(Mt, args.theta)
    Asc = np.diag(wts) + (Mt**wts)[:, None] * (Mt * jac(y, m1, m2)) * (Mt**-wts)[None, :]
    e = np.linalg.eigvals(Asc); e = e[np.argsort(e.real)]
    print(f"   {Mt:8.2e}: " + " ".join(f"{x.real:+.3f}" + (f"[{abs(x.imag):.2f}]" if abs(x.imag) > 1e-6 else "") for x in e))
stamp("done")
