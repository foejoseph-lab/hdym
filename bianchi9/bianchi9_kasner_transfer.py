"""
bianchi9_kasner.py -- Bianchi I limit: linear stability of the Kasner solution of
R/2 - C^2/(4 M^2) gravity against the Weyl ghost, as a function of M t and Kasner angle.

Bianchi I Lagrangian: keep from L (bianchi9_system.pkl, N = 1) only the e^{3 alpha} terms
(the curvature potential of Bianchi IX carries e^{alpha} and e^{-alpha}).  Then
    L_I = 3 e^{3a}(-a'^2 + b'^2) - (3 e^{3a}/M^2) b''^2 + (e^{3a}/M^2) x (cubic/quartic in a', b', b'')
and the GR Kasner solution  e^{3 alpha} = t,  alpha' = 1/(3t),  beta' = n/(3t)  (|n| = 1, n the
Kasner angle; alpha'^2 = |beta'|^2 is the Bianchi I constraint) is an exact solution of the
full equations (Einstein manifold).  With s = ln(M t) the linearised equations about Kasner
have coefficients depending only on s, so the ghost's fate depends on nothing but M t and n.

Method: proper-time rhs from L_I by the same pieces/solve as bianchi9_codegen, Jacobian by
finite differences; integrate the 10x10 tangent matrix along Kasner with Radau (implicit) in
log t, from M t = 1e2 down to 1e-2; at each output, symplectically normalised frozen ghost
frame (first adiabatic correction) and the ghost transfer matrix from the start; report the
ghost action gain log s_max versus M t, and the frozen spectrum.

Run:  %run bianchi9_kasner.py --theta 3.14159
"""
import time, argparse, pickle
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

p = argparse.ArgumentParser()
p.add_argument('--theta', type=float, default=np.pi, help='Kasner angle of beta-prime in the beta plane')
p.add_argument('--Mt-start', type=float, default=100.0)
p.add_argument('--Mt-end', type=float, default=0.01)
p.add_argument('--rtol', type=float, default=1e-9)
p.add_argument('--out', default=None)
args = p.parse_args()
t0 = time.time()
M = 1.0                                              # units: M = 1, time measured in 1/M

# ---------------------------------------------------------------- Bianchi I Lagrangian and rhs (sympy, once)
S = pickle.load(open('bianchi9_system.pkl', 'rb'))
state = sp.symbols(S['state']); Msym = sp.Symbol(S['M'], positive=True)
al0, al1, bp0, bm0, bp1, bm1, bp2, bm2, bp3, bm3 = state
KMAX = 6
J = {nm: [sp.Symbol(f'{nm}{k}') for k in range(KMAX + 1)] for nm in ('al', 'bp', 'bm', 'N')}
def D(expr): return sum(sp.diff(expr, J[nm][k]) * J[nm][k+1] for nm in J for k in range(KMAX))
def EL(Lx, nm, order=4):
    out = sp.S(0)
    for k in range(order + 1):
        term = sp.diff(Lx, J[nm][k])
        for _ in range(k): term = D(term)
        out += (-1)**k * term
    return sp.expand(out)
gauge = {J['N'][0]: 1, J['N'][1]: 0, J['N'][2]: 0, J['N'][3]: 0}
L = sp.expand(S['L'])
# keep terms whose alpha-dependence is exactly e^{3 alpha} (and N-dependence arbitrary): Bianchi I
LIN = sp.S(0)
for term in sp.Add.make_args(L):
    # substitute alpha -> alpha + c and read the exponent
    c = sp.Symbol('c')
    ratio = sp.simplify(term.subs(J['al'][0], J['al'][0] + c) / term)
    if sp.simplify(ratio - sp.exp(3 * c)) == 0:
        LIN += term
LI = LIN.subs(gauge)
assert sp.diff(LI, bp0) == 0 and sp.diff(LI, bm0) == 0, "Bianchi I Lagrangian should not depend on beta"
print(f"[{time.time()-t0:5.1f}s] Bianchi I Lagrangian: {len(sp.Add.make_args(LI))} terms")
E_al = EL(LI, 'al'); E_bp = EL(LI, 'bp'); E_bm = EL(LI, 'bm')
Con = EL(LIN, 'N').subs(gauge)                       # Hamiltonian constraint (N kept in LIN)
c1 = E_al.coeff(J['al'][2]); c0 = sp.expand(E_al - c1 * J['al'][2])
n3 = sp.expand(-(D(c0) * c1 - c0 * D(c1))); n3 = sp.expand(n3.subs(J['al'][2], -c0 / c1) * c1)
def clear(E):
    E = sp.expand(E); out = sp.S(0)
    for term in sp.Add.make_args(E):
        p2 = sp.degree(term, J['al'][2]); p3 = sp.degree(term, J['al'][3])
        base = term / J['al'][2]**p2 / J['al'][3]**p3
        out += base * (-c0)**p2 * c1**(3 - p2 - 3*p3) * n3**p3
    return sp.expand(out)
Eb = [clear(E) for E in (E_bp, E_bm)]
Amat, rvec = sp.linear_eq_to_matrix(Eb, [J['bp'][4], J['bm'][4]])
pieces = sp.lambdify([state, Msym], [c0, c1, Amat[0,0], Amat[0,1], Amat[1,0], Amat[1,1], rvec[0], rvec[1]], 'numpy', cse=True)
f_con = sp.lambdify([state, Msym], Con.subs(J['al'][2], -c0 / c1), 'numpy', cse=True)
# GR Bianchi I values of beta'', beta''' for the Einstein manifold
LGR = sp.S(0)
for term in sp.Add.make_args(sp.expand(S['L_R1'])):
    c = sp.Symbol('c'); ratio = sp.simplify(term.subs(J['al'][0], J['al'][0] + c) / term)
    if sp.simplify(ratio - sp.exp(3 * c)) == 0: LGR += term
LGR = LGR.subs(gauge)
gr2 = sp.solve([EL(LGR, 'al', 2), EL(LGR, 'bp', 2), EL(LGR, 'bm', 2)], [J['al'][2], J['bp'][2], J['bm'][2]], dict=True)[0]
gr3 = {J[nm][3]: sp.expand(D(gr2[J[nm][2]]).subs(gr2)) for nm in ('al', 'bp', 'bm')}
f_gr = sp.lambdify([state, Msym], [gr2[J['bp'][2]], gr2[J['bm'][2]], gr3[J['bp'][3]], gr3[J['bm'][3]]], 'numpy', cse=True)
print("   GR Bianchi I:  alpha'' =", sp.factor(gr2[J['al'][2]]), "   bp'' =", sp.factor(gr2[J['bp'][2]]))
print(f"[{time.time()-t0:5.1f}s] rhs built")

def rhs(y):
    c0v, c1v, a11, a12, a21, a22, r1, r2 = pieces(y, M)
    det = a11*a22 - a12*a21
    return np.array([y[1], -c0v/c1v, y[4], y[5], y[6], y[7], y[8], y[9], (a22*r1 - a12*r2)/det, (a11*r2 - a21*r1)/det])
def jac(y, h=1e-6):
    Jm = np.zeros((10, 10))
    Jm[0, 1] = Jm[2, 4] = Jm[3, 5] = Jm[4, 6] = Jm[5, 7] = Jm[6, 8] = Jm[7, 9] = 1.0
    for i in range(10):
        e = np.zeros(10); e[i] = h * (1 + abs(y[i]))
        fp = rhs(y + e); fm = rhs(y - e)
        Jm[1, i] = (fp[1] - fm[1]) / (2*e[i]); Jm[8:, i] = (fp[8:] - fm[8:]) / (2*e[i])
    return Jm

# ---------------------------------------------------------------- Kasner background and checks
def kasner(t, theta):
    y = np.zeros(10)
    y[0] = np.log(t) / 3; y[1] = 1 / (3*t)                        # e^{3 alpha} = t; the singularity is at t -> 0
    y[4], y[5] = np.cos(theta) / (3*t), np.sin(theta) / (3*t)     # |beta'| = |alpha'|
    y[6:] = f_gr(y, M)
    return y
ts = args.Mt_start
y = kasner(ts, args.theta)
f = rhs(y); yk = kasner(ts * (1 + 1e-6), args.theta)
print(f"   Kasner check: constraint {f_con(y, M):.1e};  rhs vs d/dt of Kasner family: {np.abs(f - (yk - y)/(ts*1e-6)).max() / (1 + np.abs(f).max()):.1e}")
Jf = jac(y); lam = np.sort_complex(np.linalg.eigvals(Jf))
print("   frozen eigenvalues at M t = %g:" % ts, np.array2string(lam, precision=4))

# ---------------------------------------------------------------- symplectic form in the Bianchi I limit (same construction as bianchi9_symp)
p_al = sp.expand(sp.diff(LI, J['al'][1]))
P2p = sp.expand(sp.diff(LI, bp2)); P2m = sp.expand(sp.diff(LI, bm2))
al2s = sp.Symbol('al2')
P1p = sp.expand(sp.diff(LI, bp1) - D(P2p)).subs(al2s, -c0 / c1)
P1m = sp.expand(sp.diff(LI, bm1) - D(P2m)).subs(al2s, -c0 / c1)
z = sp.Matrix([al0, bp0, bm0, bp1, bm1, p_al, P1p, P1m, P2p, P2m])
f_dz = sp.lambdify([state, Msym], z.jacobian(sp.Matrix(state)), 'numpy', cse=True)
OMC = np.block([[np.zeros((5, 5)), np.eye(5)], [-np.eye(5), np.zeros((5, 5))]])
def omega(y):
    d = np.asarray(f_dz(y, M)); return d.T @ OMC @ d
print(f"[{time.time()-t0:5.1f}s] symplectic form built")

# ---------------------------------------------------------------- tangent flow in s = ln t
Sc = np.array([1, 1, 1, 1, 1, 1, M, M, M**2, M**2]); Sinv = 1 / Sc
def F(s, Y):
    t = np.exp(s); yb = Y[:10]; Ph = Y[10:].reshape(10, 10)
    fb = rhs(yb); A = t * jac(yb)                        # d/ds = t d/dt
    Ah = (Sinv[:, None] * A) * Sc[None, :]
    return np.concatenate([t * fb, (Ah @ Ph).ravel()])
def unscale(Ph): return (Sc[:, None] * Ph) * Sinv[None, :]

def ghost_frame(y, h=1e-4):
    Jm = jac(y); lam, V = np.linalg.eig(Jm); W = np.linalg.inv(V)
    pos = np.where(lam.imag > 0)[0]; pos = pos[np.argsort(-np.abs(lam[pos].imag))][:2]
    Om = omega(y); E = V[:, pos]
    # first adiabatic correction (as in bianchi9_action --wkb)
    def proj(yy):
        l2, V2 = np.linalg.eig(jac(yy)); W2 = np.linalg.inv(V2)
        q = np.where(l2.imag > 0)[0]; q = q[np.argsort(-np.abs(l2[q].imag))][:2]
        return V2[:, q] @ W2[q, :]
    fb = rhs(y); Edot = ((proj(y + h*fb) - proj(y - h*fb)) / (2*h)) @ E
    a = W @ Edot; lp = lam[pos].mean(); delta = np.zeros_like(E)
    for k in range(10):
        if k in pos: continue
        delta += np.outer(V[:, k], a[k, :]) / (lam[k] - lp)
    E = E + delta
    G = 1j * E.conj().T @ Om @ E; w, U = np.linalg.eigh(G)
    return E @ U / np.sqrt(np.abs(w)), np.sign(w), lam[pos], lam, Om
def amplitudes(E, sgn, Om, V): return sgn[:, None] * (1j * E.conj().T @ Om @ V)

E0, sg0, lg0, lam0, Om0 = ghost_frame(y)
B0 = np.column_stack([E0[:, 0].real, E0[:, 0].imag, E0[:, 1].real, E0[:, 1].imag])
R0 = np.vstack([amplitudes(E0, sg0, Om0, B0).real, amplitudes(E0, sg0, Om0, B0).imag])
print(f"   ghost frozen frequency / M at M t = {ts:g}: {lg0.imag.mean():.5f}, Re = {lg0.real.mean():.4f} (3H/2 = {0.5/ts:.4f})")

s_grid = np.linspace(np.log(ts), np.log(args.Mt_end), 160)
Y = np.concatenate([y, np.eye(10).ravel()])
rec = dict(Mt=[], smax=[], smin=[], smean=[], lam=[], con=[], delta=[])
for k in range(1, len(s_grid)):
    sol = solve_ivp(F, (s_grid[k-1], s_grid[k]), Y, method='Radau', rtol=args.rtol, atol=1e-12)
    Y = sol.y[:, -1]; yb = Y[:10]; Phi = unscale(Y[10:].reshape(10, 10)); t = np.exp(s_grid[k])
    E, sg, lg, lam, Om = ghost_frame(yb)
    C = amplitudes(E, sg, Om, Phi @ B0); R = np.vstack([C.real, C.imag]); T = R @ np.linalg.inv(R0)
    sv = np.linalg.svd(T, compute_uv=False)
    gr_ = np.asarray(f_gr(yb, M))
    rec['Mt'].append(t); rec['smax'].append(np.log(sv[0])); rec['smin'].append(np.log(sv[-1])); rec['smean'].append(np.log(sv).mean())
    rec['lam'].append(lam); rec['con'].append(f_con(yb, M)); rec['delta'].append(np.abs(yb[6:] - gr_).max() / (1 + np.abs(gr_).max()))
    if k % 16 == 0 or t < 0.3:
        print(f"   M t = {t:8.4f}  H/M = {1/(3*t):8.4f}  log s: max {rec['smax'][-1]:+.4f} mean {rec['smean'][-1]:+.4f} min {rec['smin'][-1]:+.4f}   "
              f"ghost Re lam - 3H/2 = {lg.real.mean() - 0.5/t:+.3e}  Im/M = {lg.imag.mean():.4f}   con {rec['con'][-1]:.0e} delta {rec['delta'][-1]:.0e} ({time.time()-t0:.0f}s)", flush=True)
    if not np.isfinite(Y).all(): break
rec = {k_: np.array(v) for k_, v in rec.items()}
out = args.out or f"kasner_theta{args.theta:.3f}.npz"
np.savez(out, theta=args.theta, **rec)
print(f"[{time.time()-t0:5.1f}s] wrote {out}")
