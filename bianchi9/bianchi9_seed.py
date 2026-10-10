"""
bianchi9_seed.py -- direct repeller test in the strongly coupled regime (H_b ~ M).

Start on GR data (same trajectory as all sweeps: beta0 = -0.2, head-on into the beta_+ wall),
add a tiny seed off the Einstein manifold (delta beta'' = eps * M * H0, i.e. a ghost of relative
amplitude eps), integrate the FULL nonlinear equations with Radau (implicit; DOP853 crawls for
H/M > 0.2), and record the deviation from the manifold

    dev(t) = |(beta'', beta''') - gr(y)|  in units  (M H, M^2 H)   (ghost-amplitude units)

relative to the seed, together with the on-manifold reference (eps = 0: dev stays 0 to
roundoff; this is the invariance check in the same run).  Also recorded: the deviation of the
LOW variables (alpha, beta, their velocities) from the unseeded trajectory, which is the
physical back-reaction of the ghost on the geometry.

Run:  %run bianchi9_seed.py --Ms 1 2 3 5 --eps 1e-6 --alpha-end -1.5
"""
import time, argparse
import numpy as np
from scipy.integrate import solve_ivp
import bianchi9_rhs as m

p = argparse.ArgumentParser()
p.add_argument('--Ms', type=float, nargs='+', default=[1, 2, 3, 5])
p.add_argument('--eps', type=float, default=1e-6)
p.add_argument('--beta0', type=float, default=-0.2)
p.add_argument('--theta', type=float, default=np.pi)
p.add_argument('--alpha-end', type=float, default=-1.5)
p.add_argument('--HM-max', type=float, default=8.0)
p.add_argument('--rtol', type=float, default=1e-9)
p.add_argument('--out', default='seed.npz')
p.add_argument('--seed-pol', choices=['p','m'], default='p', help='seed polarisation: beta+ or beta-')
p.add_argument('--quiet', action='store_true')
args = p.parse_args()
t0 = time.time()

def gr_initial(Mv):
    y = np.zeros(10); y[2] = args.beta0
    y[4], y[5] = np.cos(args.theta), np.sin(args.theta); y[1] = -1.0
    for _ in range(100):
        c = m.conGR(y, Mv); d = m.dconGR(y, Mv); y[1] -= c / d
        if abs(c) < 1e-15: break
    y[6:] = m.gr(y, Mv)
    return y

def integrate(y0, Mv):
    F = lambda tau, y: np.exp(3 * y[0]) * m.rhs(y, Mv)
    def Jf(tau, y):
        N = np.exp(3 * y[0]); f = m.rhs(y, Mv); J = m.jac(y, Mv)
        return N * J + np.outer(f, 3 * N * np.eye(10)[0])
    ev1 = lambda t, y: y[0] - args.alpha_end; ev1.terminal = True
    ev2 = lambda t, y: abs(y[1]) / Mv - args.HM_max; ev2.terminal = True
    sol = solve_ivp(F, (0, 200), y0, method='Radau', jac=Jf, rtol=args.rtol, atol=1e-12,
                    events=[ev1, ev2], dense_output=True, max_step=0.02)
    return sol

results = {}
for Mv in args.Ms:
    y0 = gr_initial(Mv); H0 = abs(y0[1])
    ref = integrate(y0, Mv)
    ys = y0.copy(); ys[6 if args.seed_pol == 'p' else 7] += args.eps * Mv * H0   # seed in ghost units
    sd = integrate(ys, Mv)
    tend = min(ref.t[-1], sd.t[-1])
    taus = np.linspace(0, tend, 600)
    Yr = ref.sol(taus).T; Ys = sd.sol(taus).T
    H = np.abs(Yr[:, 1]); HM = H / Mv; alpha = Yr[:, 0]
    unit = np.array([Mv * H, Mv * H, Mv**2 * H, Mv**2 * H]).T
    dev_ref = np.linalg.norm((Yr[:, 6:] - np.array([m.gr(y, Mv) for y in Yr])) / unit, axis=1)
    dev_sd = np.linalg.norm((Ys[:, 6:] - np.array([m.gr(y, Mv) for y in Ys])) / unit, axis=1)
    lowunit = np.column_stack([np.ones_like(H), H, np.ones_like(H), np.ones_like(H), H, H])
    low = np.linalg.norm((Ys[:, :6] - Yr[:, :6]) / lowunit, axis=1)
    ib = np.argmin(Yr[:, 2])
    results[Mv] = dict(taus=taus, alpha=alpha, HM=HM, dev_ref=dev_ref, dev_sd=dev_sd, low=low, ib=ib)
    g = dev_sd / args.eps
    def at(a):
        j = np.argmin(np.abs(alpha - a)); return g[j], HM[j]
    print(f"M={Mv:4g}  H_b/M={HM[ib]:.3f}  bounce alpha={alpha[ib]:.3f}  run to alpha={alpha[-1]:.2f} (H/M={HM[-1]:.2f})  "
          f"ref-dev max {dev_ref.max():.1e}  | seed growth: pre-bounce(-0.2) {at(-0.2)[0]:.2e}  "
          f"post(-0.5) {at(-0.5)[0]:.2e}  (-0.8) {at(-0.8)[0]:.2e}  (-1.1) {at(-1.1)[0]:.2e}  end {g[-1]:.2e}  "
          f"| geometry back-reaction end {low[-1]:.1e}   ({time.time()-t0:.0f}s)", flush=True)
np.savez(args.out, eps=args.eps, **{f'M{Mv:g}_{k}': v for Mv, r in results.items() for k, v in r.items()})
print(f"   wrote {args.out}")
