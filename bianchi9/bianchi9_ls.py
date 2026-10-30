"""
bianchi9_ls.py -- Lehners & Stelle (arXiv:1909.01169) Fig. 1 initial data through the bianchi9 equations.

LS evolve EH + C^2 (+ R^2, which vacuum solutions do not see) on Bianchi IX from a(0) = 100, adot(0) = -1/3,
beta_+(0) = 1, beta_- = 0, higher-derivative data "taken from the two-derivative solution", and report that the
quadratic terms do not drastically alter BKL behaviour but speed up the crunch.  They give no value for the Weyl
coupling sigma, so their H/M is unknown.  This script runs their data at chosen H0/M.

Convention map (derive.log check 3: round S^3 with a = b = c = A0 has R = 3/(2 A0^2), so our one-forms obey
d sigma^1 = sigma^2 ^ sigma^3 and e^alpha is HALF the sphere radius; LS write (l_m/2)^2 sigma_m^2 with the same
one-forms, so their a is the radius):
    e^alpha = a_LS / 2,   beta_+ = beta_+^LS / 2,   beta_- = beta_-^LS / 2,   alpha' = adot/a = H,
    sigma_LS = 2 M^2,  kappa^2 = 2,  spin-2 mass^2 = sigma/kappa^2 = M^2.
Fig. 1 data in our variables: alpha0 = ln 50, alpha' = -1/300, beta_+ = 0.5, beta_- = beta_-' = 0, and beta_+'
from the GR constraint (two signs; LS do not say which, both are run).  Higher data on the Einstein branch,
y[6:] = gr(y), which is exactly "taken from the two-derivative solution".

Two facts the run tests:
  1. A vacuum GR solution is an EXACT solution of EH + C^2 + R^2 (Bach vanishes on Einstein spaces; the R^2
     variation vanishes at R = 0), and the Einstein branch is invariant in this code to 1e-11.  So an on-branch
     run must coincide with the GR run to roundoff at every H/M.  LS's solid curves differing from their dashed
     ones can only come from initial data off the branch or from their numerics; this run shows what on-branch
     data do.
  2. With a seed off the branch (as in bianchi9_seed.py), where in the collapse H/M crosses 0.2 and what the
     ghost does from there, at the H0/M values that bracket what LS could have used.

Outputs (written here, move to results\\): ls_rerun.npz with every trajectory (signed), and the log.
Run (Spyder):
    import sys
    %cd C:\\GitHub\\hdym\\bianchi9
    !"{sys.executable}" -B -u bianchi9_ls.py --H0M 0.01 0.1 0.3 1.0 --eps 1e-6 > results\\ls_rerun.log
Predicted wall clock: seconds per (H0/M, sign) pair with Radau; well under a minute in all.
"""
import time, argparse
import numpy as np
from scipy.integrate import solve_ivp
import bianchi9_rhs as m

p = argparse.ArgumentParser()
p.add_argument('--H0M', type=float, nargs='+', default=[0.01, 0.1, 0.3, 1.0], help='H0/M values to run')
p.add_argument('--a0', type=float, default=100.0, help='LS a(0) (sphere radius)')
p.add_argument('--adot0', type=float, default=-1.0 / 3.0, help='LS adot(0)')
p.add_argument('--bpLS', type=float, default=1.0, help='LS beta_+(0)')
p.add_argument('--eps', type=float, default=1e-6, help='seed off the Einstein branch, ghost units (0: none)')
p.add_argument('--HM-max', type=float, default=8.0)
p.add_argument('--alpha-end', type=float, default=-6.0)
p.add_argument('--rtol', type=float, default=1e-10)
p.add_argument('--out', default='ls_rerun.npz')
args = p.parse_args()
t0 = time.time()
name = 'ls_rerun'
print(f"# {name}  {time.strftime('%Y-%m-%d %H:%M:%S')}  LS data a0={args.a0} adot0={args.adot0} beta+LS={args.bpLS}", flush=True)

alpha0 = np.log(args.a0 / 2.0)
H0 = args.adot0 / args.a0                    # alpha' = adot/a, signed (negative: collapsing)
bp0 = args.bpLS / 2.0


def ls_initial(Mv, sign):
    """GR data from the LS numbers; beta_+' from the GR constraint, root of the chosen sign."""
    y = np.zeros(10)
    y[0], y[1], y[2] = alpha0, H0, bp0
    y[4] = sign * 0.01                       # Newton start for beta_+'
    for _ in range(200):
        c = m.conGR(y, Mv)
        # d conGR / d beta_+' = -6 beta_+' e^{3 alpha}
        d = -6.0 * y[4] * np.exp(3.0 * y[0])
        y[4] -= c / d
        if abs(c) < 1e-14:
            break
    y[6:] = m.gr(y, Mv)
    return y


def integrate(y0, Mv):
    F = lambda tau, y: np.exp(3 * y[0]) * m.rhs(y, Mv)
    def Jf(tau, y):
        N = np.exp(3 * y[0]); f = m.rhs(y, Mv); J = m.jac(y, Mv)
        return N * J + np.outer(f, 3 * N * np.eye(10)[0])
    ev1 = lambda t, y: y[0] - args.alpha_end; ev1.terminal = True
    ev2 = lambda t, y: abs(y[1]) / Mv - args.HM_max; ev2.terminal = True
    return solve_ivp(F, (0, 1e3), y0, method='Radau', jac=Jf, rtol=args.rtol, atol=1e-13,
                     events=[ev1, ev2], dense_output=True, max_step=0.02)


results = {}
for H0M in args.H0M:
    Mv = abs(H0) / H0M
    for sign in (+1, -1):
        y0 = ls_initial(Mv, sign)
        c0 = m.conGR(y0, Mv)
        ref = integrate(y0, Mv)                                 # on-branch: LS's "higher data from GR"
        taus = np.linspace(0, ref.t[-1], 1200)
        Yr = ref.sol(taus).T
        dev_ref = np.linalg.norm((Yr[:, 6:] - np.array([m.gr(y, Mv) for y in Yr])) /
                                 np.column_stack([Mv * np.abs(Yr[:, 1])] * 2 + [Mv**2 * np.abs(Yr[:, 1])] * 2), axis=1)
        HM = np.abs(Yr[:, 1]) / Mv
        cross = np.argmax(HM > 0.2) if np.any(HM > 0.2) else -1
        # Mixmaster bookkeeping on the GR branch: bounces = sign changes of beta_+' (head-on wall) and of beta_-'
        scale = np.abs(Yr[:, 4:6]).max()
        def nflips(v):  # sign changes of a velocity, ignoring roundoff about an identically zero one
            v = v[np.abs(v) > 1e-9 * scale]
            return int(np.sum(np.diff(np.sign(v)) != 0)) if v.size else 0
        nb_p, nb_m = nflips(Yr[:, 4]), nflips(Yr[:, 5])
        line = (f"H0/M={H0M:5g}  M={Mv:.4e}  sign={sign:+d}  beta+'(0)={y0[4]:+.5f}  con0={c0:.1e}  "
                f"ref: alpha {Yr[0,0]:.3f}->{Yr[-1,0]:.3f}  H/M {HM[0]:.3f}->{HM[-1]:.3f}  "
                f"H/M crosses 0.2 at alpha={Yr[cross,0]:.3f}" if cross >= 0 else
                f"H0/M={H0M:5g}  M={Mv:.4e}  sign={sign:+d}  beta+'(0)={y0[4]:+.5f}  con0={c0:.1e}  "
                f"ref: alpha {Yr[0,0]:.3f}->{Yr[-1,0]:.3f}  H/M {HM[0]:.3f}->{HM[-1]:.3f}  never crosses 0.2")
        line += f"  wall bounces (beta+', beta-' sign changes) {nb_p}, {nb_m}  on-branch dev max {dev_ref.max():.1e}"
        rec = dict(taus=taus, Y=Yr, HM=HM, dev_ref=dev_ref, M=Mv, sign=sign, y0=y0)
        if args.eps > 0:
            ys = y0.copy(); ys[6] += args.eps * Mv * abs(H0)
            sd = integrate(ys, Mv)
            tend = min(ref.t[-1], sd.t[-1]); taus2 = np.linspace(0, tend, 1200)
            Yr2 = ref.sol(taus2).T; Ys = sd.sol(taus2).T
            unit = np.column_stack([Mv * np.abs(Yr2[:, 1])] * 2 + [Mv**2 * np.abs(Yr2[:, 1])] * 2)
            dev_sd = np.linalg.norm((Ys[:, 6:] - np.array([m.gr(y, Mv) for y in Ys])) / unit, axis=1) / args.eps
            lowunit = np.column_stack([np.ones(len(taus2)), np.abs(Yr2[:, 1]), np.ones(len(taus2)),
                                       np.ones(len(taus2)), np.abs(Yr2[:, 1]), np.abs(Yr2[:, 1])])
            low = np.linalg.norm((Ys[:, :6] - Yr2[:, :6]) / lowunit, axis=1)
            HM2 = np.abs(Yr2[:, 1]) / Mv
            def at(hm):
                j = np.argmin(np.abs(HM2 - hm)); return dev_sd[j]
            line += (f"  | seeded (eps={args.eps:g}): growth at H/M=0.1 {at(0.1):.2e}, 0.2 {at(0.2):.2e}, "
                     f"0.5 {at(0.5):.2e}, 1 {at(1.0):.2e}, end {dev_sd[-1]:.2e}; geometry back-reaction end {low[-1]:.1e}")
            rec.update(taus_seed=taus2, Y_seed=Ys, dev_seed=dev_sd, low=low)
        results[(H0M, sign)] = rec
        print(line, flush=True)

np.savez(args.out, **{f"H0M{k[0]}_s{'p' if k[1] > 0 else 'm'}": np.array(v, dtype=object) for k, v in results.items()},
         alpha0=alpha0, H0=H0, bp0=bp0, eps=args.eps)
print(f"# END {name} elapsed={time.time() - t0:.1f}s status=ok", flush=True)
