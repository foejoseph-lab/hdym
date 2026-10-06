"""
bianchi9_action.py -- exact ghost action and ghost transfer matrix across a mixmaster bounce

Along a trajectory on the Einstein manifold (GR data, full quadratic-gravity flow), integrate
the 10x10 tangent matrix Phi (Misner time, dt = e^{3 alpha} d tau).  At each output point:

  * frozen Jacobian J(y) of the proper-time flow, eigen-decomposition; the two eigenvectors
    with Im(lambda) > 0 of largest |Im| are the Weyl-ghost polarisations (frequency ~ M);
  * symplectic normalisation with the Ostrogradsky form Omega(y) (bianchi9_symp_gen):
    G = i E^H Omega E (Hermitian 2x2), E -> E G^{-1/2}; then for a real tangent vector v the
    complex amplitudes c = i E^H Omega v and the ghost action is I = |c|^2 (exact adiabatic
    invariant of each polarisation; no phase wobble);
  * ghost transfer matrix: T(tau) = real 4x4 map from (Re c, Im c) at tau=0 to the same at tau,
    built from Phi.  Its singular values s_i are the action-gain factors (sqrt) for the best
    and worst ghost phase; log(s_max) and mean log s are the per-bounce observables.
  * frozen spectrum: all 10 eigenvalues, saved, for the Kasner-epoch onset analysis.

Checks: Omega(Phi v1, Phi v2) = Omega(v1, v2) along the run (symplecticity of the tangent flow),
and the invariance/constraint numbers as in bianchi9_test.py.

Run:   %run bianchi9_action.py --M 20 --beta0 -0.2 --alpha-end -0.9
"""
import time, argparse
import numpy as np
from scipy.integrate import solve_ivp
import bianchi9_rhs as m
import bianchi9_symp_gen as g

p = argparse.ArgumentParser()
p.add_argument('--M', type=float, default=20.0)
p.add_argument('--alpha0', type=float, default=0.0)
p.add_argument('--beta0', type=float, default=-0.2)
p.add_argument('--theta', type=float, default=np.pi)
p.add_argument('--alpha-end', type=float, default=-0.9)
p.add_argument('--HM-max', type=float, default=1e9)
p.add_argument('--rtol', type=float, default=1e-10)
p.add_argument('--dtau', type=float, default=0.02)
p.add_argument('--out', default=None)
args = p.parse_args()
Mv = args.M
t0 = time.time()

def gr_initial(alpha0, theta, beta0):
    y = np.zeros(10); y[0] = alpha0; y[2] = beta0
    y[4], y[5] = np.cos(theta), np.sin(theta); y[1] = -1.0
    for _ in range(100):
        c = m.conGR(y, Mv); d = m.dconGR(y, Mv); y[1] -= c / d
        if abs(c) < 1e-15: break
    y[6:] = m.gr(y, Mv)
    return y

# The tangent matrix is integrated in scaled coordinates, Phi = Sc PhiHat Sc^-1 with
# Sc = diag(1,..,1, M, M, M^2, M^2): a ghost mode has delta beta'' ~ M delta beta and
# delta beta''' ~ M^2 delta beta, so the scaled entries are O(1) and a uniform atol is
# meaningful. Unscaled, the entries span 1e-6..1e4 and the integrator crawls.
Sc = np.array([1, 1, 1, 1, 1, 1, Mv, Mv, Mv**2, Mv**2]); Sinv = 1 / Sc
def F(tau, Y):
    y = Y[:10]; Ph = Y[10:].reshape(10, 10)
    N = np.exp(3 * y[0]); f = m.rhs(y, Mv); J = m.jac(y, Mv)
    A = N * J + np.outer(f, 3 * N * np.eye(10)[0])        # d(N f)/dy
    Ah = (Sinv[:, None] * A) * Sc[None, :]                 # Sc^-1 A Sc
    return np.concatenate([N * f, (Ah @ Ph).ravel()])
def unscale(Ph):
    return (Sc[:, None] * Ph) * Sinv[None, :]

def ghost_frame(y):
    """symplectically normalised ghost eigenvectors E (10x2 complex), frequencies, and all eigenvalues"""
    J = m.jac(y, Mv); lam, V = np.linalg.eig(J)
    Om = g.omega(y, Mv)
    pos = np.where(lam.imag > 0)[0]
    pos = pos[np.argsort(-np.abs(lam[pos].imag))][:2]
    E = V[:, pos]
    G = 1j * E.conj().T @ Om @ E                           # Hermitian
    w, U = np.linalg.eigh(G)
    sgn = np.sign(w)
    E = E @ U / np.sqrt(np.abs(w))                         # now i E^H Om E = diag(sgn)
    return E, sgn, lam[pos], lam, Om

def amplitudes(E, sgn, Om, V):
    """complex ghost amplitudes of the columns of V (real 10xk):  c = sgn * i E^H Om V"""
    return (sgn[:, None]) * (1j * E.conj().T @ Om @ V)

y0 = gr_initial(args.alpha0, args.theta, args.beta0)
E0, sgn0, lam0, _, Om0 = ghost_frame(y0)
print(f"   M = {Mv}, initial H/M = {abs(y0[1])/Mv:.4f}; ghost frozen frequencies / M = "
      f"{lam0.imag[0]/Mv:.5f}, {lam0.imag[1]/Mv:.5f}; Krein signs {sgn0}; Re lambda = {lam0.real}")
# real basis of the ghost plane: columns Re e1, Im e1, Re e2, Im e2
B0 = np.column_stack([E0[:, 0].real, E0[:, 0].imag, E0[:, 1].real, E0[:, 1].imag])
C0 = amplitudes(E0, sgn0, Om0, B0)                          # 2x4 complex
R0 = np.vstack([C0.real, C0.imag])                          # 4x4 real: basis -> (Re c, Im c)
def dH(y, h=1e-6):
    return np.array([(m.con(y + h*np.eye(10)[i], Mv) - m.con(y - h*np.eye(10)[i], Mv)) / (2*h) for i in range(10)])
# two random tangent vectors IN THE CONSTRAINT SURFACE for the symplecticity check: the
# Misner-time flow is Hamiltonian only on H = 0, so dH.v must vanish for Omega to be preserved
rng = np.random.default_rng(1)
g0 = dH(y0)
def project(v, g): return v - g * (g @ v) / (g @ g)
v1, v2 = project(rng.normal(size=10), g0), project(rng.normal(size=10), g0)
om12_0 = v1 @ Om0 @ v2
print(f"   ghost modes in the constraint surface? |dH.e|/|dH||e| = "
      f"{np.abs(g0 @ E0).max() / np.linalg.norm(g0):.1e}")

Y = np.concatenate([y0, np.eye(10).ravel()]); tau = 0.0
rec = dict(tau=[], alpha=[], HM=[], smax=[], smin=[], slogmean=[], lam=[], symp=[], delta=[], con=[], bp=[], bm=[])
stop = lambda t, Y: Y[0] - args.alpha_end; stop.terminal = True
while True:
    sol = solve_ivp(F, (tau, tau + args.dtau), Y, method='DOP853', rtol=args.rtol, atol=1e-12,
                    events=stop, max_step=args.dtau)
    Y = sol.y[:, -1]; tau = sol.t[-1]
    y = Y[:10]; Phi = unscale(Y[10:].reshape(10, 10))
    E, sgn, lamg, lam, Om = ghost_frame(y)
    C = amplitudes(E, sgn, Om, Phi @ B0)                    # 2x4: image of the initial ghost basis
    R = np.vstack([C.real, C.imag])
    T = R @ np.linalg.inv(R0)                               # 4x4 real transfer matrix on (Re c, Im c)
    s = np.linalg.svd(T, compute_uv=False)
    gr_ = m.gr(y, Mv)
    rec['tau'].append(tau); rec['alpha'].append(y[0]); rec['HM'].append(abs(y[1]) / Mv)
    rec['smax'].append(np.log(s[0])); rec['smin'].append(np.log(s[-1])); rec['slogmean'].append(np.log(s).mean())
    rec['lam'].append(lam); rec['symp'].append((Phi @ v1) @ Om @ (Phi @ v2) / om12_0 - 1)
    rec['delta'].append(np.abs(y[6:] - gr_).max() / (1 + np.abs(gr_).max())); rec['con'].append(m.con(y, Mv))
    rec['bp'].append(y[2]); rec['bm'].append(y[3])
    if len(rec['tau']) % 25 == 0:
        print(f"   tau={tau:5.2f} alpha={y[0]:7.3f} H/M={rec['HM'][-1]:.4f}  log s: max {rec['smax'][-1]:+.3f} "
              f"mean {rec['slogmean'][-1]:+.3f} min {rec['smin'][-1]:+.3f}  Re lam_ghost/M={lamg.real.max()/Mv:+.1e} "
              f"symp err {abs(rec['symp'][-1]):.1e} ({time.time()-t0:.0f}s)", flush=True)
    if sol.status == 1 or rec['HM'][-1] > args.HM_max or not np.isfinite(y).all():
        break
rec = {k: np.array(v) for k, v in rec.items()}
ib = np.argmin(rec['bp'])
print(f"[{time.time()-t0:5.1f}s] done: alpha {rec['alpha'][0]:.3f} -> {rec['alpha'][-1]:.3f}, H/M {rec['HM'][0]:.4f} -> {rec['HM'][-1]:.4f}, "
      f"bounce at alpha = {rec['alpha'][ib]:.3f} (H_b/M = {rec['HM'][ib]:.4f})")
print(f"   invariance {rec['delta'].max():.1e}   constraint {np.abs(rec['con']).max():.1e}   symplecticity {np.abs(rec['symp']).max():.1e}")
print(f"   ghost transfer across the run:  log s_max = {rec['smax'][-1]:+.4f}   mean log s = {rec['slogmean'][-1]:+.4f}   log s_min = {rec['smin'][-1]:+.4f}")
out = args.out or f"action_M{Mv:g}.npz"
np.savez(out, M=Mv, **rec)
print(f"   wrote {out}")
