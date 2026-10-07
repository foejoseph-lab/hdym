"""
bianchi9_bogoliubov.py -- first-order (Bogoliubov / adiabatic perturbation theory) squeeze of
the Weyl ghost across one mixmaster bounce, as a quadrature along the exact GR background.

Along the Einstein-manifold trajectory y(tau) (Misner time, dt = e^{3 alpha} d tau) the
linearised flow has, at each instant, a ghost pair of eigenmodes E+ (eigenvalues
lambda = gamma + i omega, two polarisations) and the conjugate pair E- = conj(E+).  With the
symplectic normalisation i E^H Omega E = +-1, first-order adiabatic perturbation theory gives
the amplitude transferred from the + pair into the - pair:

    B = int d tau  C(tau) exp(-i phi(tau)),     C = i E-^H Omega dE+/dtau,
    phi(tau) = int 2 omega_tau d tau,  omega_tau = e^{3 alpha} Im lambda  (Misner-time frequency)

(real parts of lambda cancel between + and -).  dE+/dtau is taken through the projector onto
the pair, so only the part of the motion that leaves the pair subspace counts (within-pair
rotations are the adiabatic phase).  The squeeze is r = |B| (largest singular value of the
2x2 B); this is the perturbative version of the transfer-matrix r of bianchi9_action.py and is
valid when r << 1 -- exactly the regime of interest.  There is no frame artefact of the
frozen-mode type here because the formula is first order in the coupling by construction;
its own error is O(r^2) and O(H/M) relative.

Also returned: the integrand |C| and the local phase rate, so the saddle-point exponent can be
read off, and the position of the bounce.

Run:  %run bianchi9_bogoliubov.py --Ms 28 40 56 80 160 320 640
"""
import time, argparse
import numpy as np
from scipy.integrate import solve_ivp
import bianchi9_rhs as m
import bianchi9_symp_gen as g

p = argparse.ArgumentParser()
p.add_argument('--Ms', type=float, nargs='+', default=[28, 40, 56, 80, 160, 320])
p.add_argument('--beta0', type=float, default=-0.2)
p.add_argument('--theta', type=float, default=np.pi)
p.add_argument('--alpha0', type=float, default=0.0)
p.add_argument('--alpha-end', type=float, default=-0.7)
p.add_argument('--tau-start', type=float, default=0.0)
p.add_argument('--n', type=int, default=4000, help='quadrature points in tau')
p.add_argument('--h', type=float, default=1e-4, help='finite-difference step for dE/dtau (proper time)')
p.add_argument('--out', default='bogoliubov.npz')
args = p.parse_args()
t0 = time.time()

# ---------------------------------------------------------------- background (M-independent on the Einstein manifold)
Mbg = 300.0
def gr_initial(alpha0, theta, beta0, Mv):
    y = np.zeros(10); y[0] = alpha0; y[2] = beta0
    y[4], y[5] = np.cos(theta), np.sin(theta); y[1] = -1.0
    for _ in range(100):
        c = m.conGR(y, Mv); d = m.dconGR(y, Mv); y[1] -= c / d
        if abs(c) < 1e-15: break
    y[6:] = m.gr(y, Mv)
    return y
y0 = gr_initial(args.alpha0, args.theta, args.beta0, Mbg)
Fbg = lambda tau, y: np.exp(3 * y[0]) * m.rhs(y, Mbg)
stop = lambda t, y: y[0] - args.alpha_end; stop.terminal = True
sol = solve_ivp(Fbg, (0, 50), y0, method='DOP853', rtol=1e-11, atol=1e-13, events=stop, dense_output=True, max_step=0.02)
tau_end = sol.t[-1]
taus = np.linspace(args.tau_start, tau_end, args.n)
Y = sol.sol(taus).T                                   # n x 10, GR trajectory
ib = np.argmin(Y[:, 2])
print(f"[{time.time()-t0:5.1f}s] background: tau_end = {tau_end:.3f}, bounce at tau = {taus[ib]:.3f}, alpha_b = {Y[ib,0]:.3f}, "
      f"H_b = {abs(Y[ib,1]):.3f}", flush=True)

# ---------------------------------------------------------------- frames
def pair(y, Mv):
    """eigen-decomposition; returns lam, V, W=V^-1, index of the + pair, and projector P"""
    J = m.jac(y, Mv); lam, V = np.linalg.eig(J)
    pos = np.where(lam.imag > 0)[0]
    pos = pos[np.argsort(-np.abs(lam[pos].imag))][:2]
    W = np.linalg.inv(V)
    return lam, V, W, pos, V[:, pos] @ W[pos, :]

def normalise(E, Om):
    G = 1j * E.conj().T @ Om @ E
    w, U = np.linalg.eigh(G)
    return E @ U / np.sqrt(np.abs(w)), np.sign(w)

results = {}
for Mv in args.Ms:
    # along the grid: ghost frame, coupling, phase rate
    C = np.zeros((args.n, 2, 2), complex); om = np.zeros(args.n); gam = np.zeros(args.n)
    Eprev = None
    for k, y in enumerate(Y):
        # on the Einstein manifold the state is M-independent, but y[6:] (beta'', beta''') were
        # computed with Mbg; they are GR values, identical for any M (invariance). Fine.
        lam, V, W, pos, P = pair(y, Mv)
        Om = g.omega(y, Mv)
        Ep = V[:, pos]
        Ep, sg = normalise(Ep, Om)
        # parallel transport of the within-pair frame: the symplectic normalisation leaves a
        # 2x2 unitary free at every point; fix it so the overlap with the previous frame
        # (in the Hermitian form i E^H Om E) is Hermitian positive (polar decomposition).
        if Eprev is not None:
            O = 1j * Eprev.conj().T @ Om @ Ep
            Uo, _, Vho = np.linalg.svd(O)
            Ep = Ep @ (Uo @ Vho).conj().T
        Eprev = Ep
        Em = Ep.conj()                                      # conjugate pair, eigenvalues conj(lam)
        f = m.rhs(y, Mv)
        Pp = pair(y + args.h * f, Mv)[4]; Pm = pair(y - args.h * f, Mv)[4]
        Edot_t = ((Pp - Pm) / (2 * args.h)) @ Ep            # d/dt (proper time) of the pair, applied to E+
        N = np.exp(3 * y[0])
        Edot = N * Edot_t                                   # d/dtau
        C[k] = 1j * Em.conj().T @ Om @ Edot                 # coupling + -> -, per unit tau
        om[k] = N * lam[pos].imag.mean()                    # Misner-time frequency
        gam[k] = lam[pos].real.mean()
    phi = np.concatenate([[0], np.cumsum(0.5 * (2 * om[1:] + 2 * om[:-1]) * np.diff(taus))])
    integrand = C * np.exp(-1j * phi)[:, None, None]
    B = np.trapezoid(integrand, taus, axis=0)
    s = np.linalg.svd(B, compute_uv=False)
    # where does the integrand live?  cumulative |B| and the modulus of the coupling
    cum = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(taus)[:, None, None], axis=0)
    cumr = np.array([np.linalg.svd(c, compute_uv=False)[0] for c in cum])
    Cn = np.array([np.linalg.norm(c) for c in C])
    kmax = np.argmax(Cn)
    results[Mv] = dict(B=B, r=s[0], r2=s[1], cumr=cumr, Cn=Cn, om=om, C=C)
    print(f"   M={Mv:6g}  H_b/M={abs(Y[ib,1])/Mv:.4f}   r_Bog = {s[0]:.3e} (second {s[1]:.3e})   "
          f"|C| peak at tau={taus[kmax]:.3f} (bounce {taus[ib]:.3f}), |C|max={Cn[kmax]:.3e}, "
          f"phase rate at bounce 2w={2*om[ib]:.1f}   ({time.time()-t0:.0f}s)", flush=True)

np.savez(args.out, taus=taus, Y=Y, Ms=np.array(args.Ms), ib=ib,
         r=np.array([results[M]['r'] for M in args.Ms]),
         cumr=np.array([results[M]['cumr'] for M in args.Ms]),
         Cn=np.array([results[M]['Cn'] for M in args.Ms]),
         om=np.array([results[M]['om'] for M in args.Ms]),
         C=np.array([results[M]['C'] for M in args.Ms]))
print(f"   wrote {args.out}")
