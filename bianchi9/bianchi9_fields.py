"""
bianchi9_fields.py -- healthy test fields on the (unchanged) Bianchi IX background:
a massive scalar (boson) and a massive Dirac field (fermion), homogeneous (k = 0), linear,
no back-reaction.  Same GR trajectory as the ghost runs (beta0 = -0.2, head-on into the
beta_+ -> -inf wall).  Observable: adiabatic particle number n(t) produced by the bounce,
versus m/H_b, to compare with the Weyl ghost's n = sinh^2 r.

Scalar:  phi'' + 3 alpha' phi' + m^2 phi = 0  (proper time); chi = e^{3alpha/2} phi obeys
         chi'' + W^2 chi = 0,  W^2 = m^2 - (9/4) alpha'^2 - (3/2) alpha''.
         n = (|chi'|^2 + W^2 |chi|^2)/(2W) - 1/2  for Wronskian-normalised chi.
         The scalar couples to the volume only: it does not see the anisotropy at k = 0.

Dirac:   (gamma^a e_a^mu nabla_mu + m) psi = 0 (mostly-plus Clifford algebra) with the orthonormal frame e^0 = N dt,
         e^i = a_i sigma^i and the spin connection Gamma_mu = (1/4) omega_{mu ab} gamma^a gamma^b.
         For homogeneous psi(t) this is  d psi/dt = X(t) psi  (4x4).  With chi = e^{3alpha/2} psi,
         i d chi/dt = H_D(t) chi, H_D Hermitian (checked numerically): H_D = m beta_D + A(t) Sigma,
         where A = (1/4)(a/bc + b/ca + c/ab) is the axial coupling (checked against the eigenvalues) from the Bianchi IX structure
         constants -- the fermion sees the anisotropy through the SQUARE ROOTS of the curvature
         walls (a/bc = e^{-alpha + 2 beta_+ + 2 sqrt3 beta_-}, etc.).  Eigenvalues +-E(t),
         E = sqrt(m^2 + A^2).  n = norm^2 of the projection of chi onto the instantaneous
         negative-energy subspace, starting from a positive-energy eigenstate.

Run:  %run bianchi9_fields.py --ms 20 28 40 56 80
"""
import time, argparse, pickle
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import bianchi9_rhs as m_

p = argparse.ArgumentParser()
p.add_argument('--ms', type=float, nargs='+', default=[20, 28, 40, 56, 80])
p.add_argument('--beta0', type=float, default=-0.2)
p.add_argument('--theta', type=float, default=np.pi)
p.add_argument('--alpha-end', type=float, default=-0.7)
p.add_argument('--rtol', type=float, default=1e-10)
p.add_argument('--out', default='fields.npz')
args = p.parse_args()
t0 = time.time()

# ---------------------------------------------------------------- background (GR, M-independent)
Mbg = 300.0
y = np.zeros(10); y[2] = args.beta0; y[4], y[5] = np.cos(args.theta), np.sin(args.theta); y[1] = -1.0
for _ in range(100):
    c = m_.conGR(y, Mbg); d = m_.dconGR(y, Mbg); y[1] -= c / d
    if abs(c) < 1e-15: break
y[6:] = m_.gr(y, Mbg)
Fbg = lambda tau, yy: np.exp(3 * yy[0]) * m_.rhs(yy, Mbg)
stop = lambda t, yy: yy[0] - args.alpha_end; stop.terminal = True
bg = solve_ivp(Fbg, (0, 50), y, method='DOP853', rtol=1e-11, atol=1e-13, events=stop, dense_output=True, max_step=0.02)
tau_end = bg.t[-1]
Yb = lambda tau: bg.sol(tau)
tb = np.linspace(0, tau_end, 3000); Yg = bg.sol(tb).T; ib = np.argmin(Yg[:, 2]); tau_b = tb[ib]; Hb = abs(Yg[ib, 1])
print(f"[{time.time()-t0:5.1f}s] background: bounce at tau = {tau_b:.3f}, alpha_b = {Yg[ib,0]:.3f}, H_b = {Hb:.3f}, tau_end = {tau_end:.3f}")

# ---------------------------------------------------------------- Dirac reduction (sympy, frame connection as in bianchi9_derive)
t = sp.symbols('t'); N_, a_, b_, c_ = [sp.Function(s)(t) for s in ('N', 'a', 'b', 'c')]
A = [None, a_, b_, c_]; eta = [-1, 1, 1, 1]; idx = range(4)
F = [[[sp.S(0)]*4 for _ in idx] for _ in idx]
for i in (1, 2, 3):
    F[i][0][i] = sp.diff(A[i], t) / (N_ * A[i]); F[i][i][0] = -F[i][0][i]
    for j in (1, 2, 3):
        for k in (1, 2, 3):
            e = sp.LeviCivita(i, j, k)
            if e != 0: F[i][j][k] = e * A[i] / (A[j] * A[k])
Fl = [[[eta[p_] * F[p_][q][r] for r in idx] for q in idx] for p_ in idx]
w = [[[sp.simplify((Fl[p_][q][r] + Fl[q][r][p_] - Fl[r][p_][q]) / 2) for r in idx] for q in idx] for p_ in idx]   # omega_{abc}
# gamma matrices with {g^a, g^b} = 2 eta^{ab}, eta = diag(-1,1,1,1):  g^a = i * (standard Dirac gamma^a)
I2 = np.eye(2); Z2 = np.zeros((2, 2))
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
g0s = np.block([[I2, Z2], [Z2, -I2]]).astype(complex)
gis = [np.block([[Z2, s], [-s, Z2]]) for s in (sx, sy, sz)]
gam = [1j * g0s] + [1j * g for g in gis]
for a1 in range(4):
    for b1 in range(4):
        assert np.allclose(gam[a1] @ gam[b1] + gam[b1] @ gam[a1], 2 * eta[a1] * (a1 == b1) * np.eye(4))
# numerical connection along the background: omega_{abc}(a, b, c, a', b', c'), N = 1
syms = sp.symbols('a b c ad bd cd')
rep = {sp.Derivative(a_, t): syms[3], sp.Derivative(b_, t): syms[4], sp.Derivative(c_, t): syms[5]}
w_num = sp.lambdify(syms, [[[w[p_][q][r].subs(N_, 1).subs(rep).subs({a_: syms[0], b_: syms[1], c_: syms[2]}) for r in idx] for q in idx] for p_ in idx], 'numpy')
s3 = np.sqrt(3)
def abc(yy):
    al, bp, bm = yy[0], yy[2], yy[3]; alp, bpp, bmp = yy[1], yy[4], yy[5]
    a = np.exp(al + bp + s3*bm); b = np.exp(al + bp - s3*bm); c = np.exp(al - 2*bp)
    return a, b, c, a*(alp + bpp + s3*bmp), b*(alp + bpp - s3*bmp), c*(alp - 2*bpp)
def dirac_X(yy, m):
    """psi' = X psi (proper time, N = 1):  g^0 psi' + (1/4) omega_{abc} g^c g^a g^b psi + m psi = 0"""
    om = np.array(w_num(*abc(yy)), float)
    S = np.zeros((4, 4), complex)
    for a1 in idx:
        for b1 in idx:
            for c1 in idx:
                if om[a1, b1, c1] != 0:
                    S += om[a1, b1, c1] * (gam[c1] @ gam[a1] @ gam[b1])
    # mostly-plus Clifford algebra: (g^a nabla_a + m) psi = 0  ->  g^0 psi' = -(m + S/4) psi
    return np.linalg.solve(gam[0], -(m * np.eye(4) + 0.25 * S))
def dirac_H(yy, m):
    """H_D for chi = e^{3 alpha/2} psi:  i chi' = H chi,  H = i X + (3/2) i alpha' * 1"""
    return 1j * dirac_X(yy, m) + 1.5j * yy[1] * np.eye(4)
H0 = dirac_H(Yb(0.0), 1.0)
herm = np.abs(H0 - H0.conj().T).max() / np.abs(H0).max()
ev0 = np.linalg.eigvalsh((H0 + H0.conj().T) / 2)
a0, b0, c0 = abc(Yb(0.0))[:3]
Aax = (a0/(b0*c0) + b0/(c0*a0) + c0/(a0*b0)) / 8
print(f"   Dirac reduction: hermiticity of H_D {herm:.1e};  eigenvalues at tau=0 (m=1): {np.round(ev0, 5)};  "
      f"sqrt(m^2 + A^2) with A = (1/4) sum a/bc = {np.sqrt(1 + 4*Aax**2):.5f}")
print(f"[{time.time()-t0:5.1f}s] Dirac operator built")

# ---------------------------------------------------------------- integrations
def scalar_run(m):
    # chi'' + W^2 chi = 0 in proper time; integrate in tau: d/dtau = N d/dt, N = e^{3 alpha}
    def W2(yy):
        al2 = m_.rhs(yy, Mbg)[1]
        return m**2 - 2.25 * yy[1]**2 - 1.5 * al2
    y0 = Yb(0.0); W0 = np.sqrt(W2(y0))
    chi0 = np.array([1/np.sqrt(2*W0), -1j*np.sqrt(W0/2)])      # positive-frequency WKB state
    def Fs(tau, Z):
        yy = Yb(tau); Nn = np.exp(3*yy[0]); chi, dchi = Z[0] + 1j*Z[1], Z[2] + 1j*Z[3]
        ddchi = -W2(yy) * chi
        return np.array([ (Nn*dchi).real, (Nn*dchi).imag, (Nn*ddchi).real, (Nn*ddchi).imag ])
    Z0 = np.array([chi0[0].real, chi0[0].imag, chi0[1].real, chi0[1].imag])
    taus = np.linspace(0, tau_end, 400)
    sol = solve_ivp(Fs, (0, tau_end), Z0, method='DOP853', rtol=args.rtol, atol=1e-14, t_eval=taus, max_step=0.01)
    n = []
    for k, tau in enumerate(sol.t):
        yy = Yb(tau); W = np.sqrt(W2(yy)); chi = sol.y[0, k] + 1j*sol.y[1, k]; dchi = sol.y[2, k] + 1j*sol.y[3, k]
        n.append((abs(dchi)**2 + W**2 * abs(chi)**2) / (2*W) - 0.5)
    return sol.t, np.array(n)

def dirac_run(m):
    y0 = Yb(0.0); H = dirac_H(y0, m); ev, U = np.linalg.eigh((H + H.conj().T)/2)
    chi0 = U[:, np.argmax(ev)]                                  # a positive-energy eigenstate
    def Fd(tau, Z):
        yy = Yb(tau); Nn = np.exp(3*yy[0]); chi = Z[:4] + 1j*Z[4:]
        d = -1j * Nn * (dirac_H(yy, m) @ chi)
        return np.concatenate([d.real, d.imag])
    Z0 = np.concatenate([chi0.real, chi0.imag])
    taus = np.linspace(0, tau_end, 400)
    sol = solve_ivp(Fd, (0, tau_end), Z0, method='DOP853', rtol=args.rtol, atol=1e-14, t_eval=taus, max_step=0.01)
    n = []
    for k, tau in enumerate(sol.t):
        yy = Yb(tau); H = dirac_H(yy, m); ev, U = np.linalg.eigh((H + H.conj().T)/2)
        chi = sol.y[:4, k] + 1j*sol.y[4:, k]
        neg = U[:, ev < 0]
        n.append(np.linalg.norm(neg.conj().T @ chi)**2 / np.linalg.norm(chi)**2)
    return sol.t, np.array(n)

alpha_of = lambda taus: np.array([Yb(tau)[0] for tau in taus])
res = {}
print(f"\n   m      m/H_b     scalar n(alpha=-0.5)  n(end)      Dirac n(alpha=-0.5)  n(end)     [ghost r^2 at -0.5, frozen-frame upper bounds: 0.049 -> <1e-5, 0.035 -> <3e-7]")
for m in args.ms:
    ts_, ns = scalar_run(m); td_, nd = dirac_run(m)
    al_s = alpha_of(ts_); al_d = alpha_of(td_)
    js = np.argmin(np.abs(al_s + 0.5)); jd = np.argmin(np.abs(al_d + 0.5))
    res[m] = dict(tau=ts_, alpha=al_s, n_scalar=ns, n_dirac=nd)
    print(f"   {m:5g}   {Hb/m:.4f}     {ns[js]:.3e}             {ns[-1]:.3e}   {nd[jd]:.3e}            {nd[-1]:.3e}      ({time.time()-t0:.0f}s)", flush=True)
np.savez(args.out, ms=np.array(args.ms), Hb=Hb, tau_b=tau_b,
         **{f"tau_{m:g}": res[m]['tau'] for m in args.ms}, **{f"alpha_{m:g}": res[m]['alpha'] for m in args.ms},
         **{f"ns_{m:g}": res[m]['n_scalar'] for m in args.ms}, **{f"nd_{m:g}": res[m]['n_dirac'] for m in args.ms})
print(f"   wrote {args.out}")
