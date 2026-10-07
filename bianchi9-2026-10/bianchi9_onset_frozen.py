"""
bianchi9_onset_frozen.py -- frozen spectrum of the ghost along the collapse vs H/M.

Along the GR background (same trajectory as the sweeps, extended to alpha_end) and for a list
of M, compute at each point the eigenvalues of the proper-time Jacobian J(y). The ghost pair
has Re lambda ~ 3H/2 kinematically (amplitude ~ a^{-3/2} in a collapse). Report
   excess = Re lambda_ghost - 1.5 H      (per unit proper time, in units of H)
and the ghost frequency / M, as functions of H/M. Also the largest real eigenvalue of the
whole spectrum (mixmaster + ghost) for reference.
Run: %run bianchi9_onset_frozen.py --Ms 2 5 10 20 --alpha-end -2.0
"""
import argparse, numpy as np
from scipy.integrate import solve_ivp
import bianchi9_rhs as m
p = argparse.ArgumentParser()
p.add_argument('--Ms', type=float, nargs='+', default=[2, 5, 10, 20])
p.add_argument('--alpha-end', type=float, default=-2.0)
p.add_argument('--beta0', type=float, default=-0.2)
p.add_argument('--n', type=int, default=400)
p.add_argument('--out', default='onset_frozen.npz')
a = p.parse_args()
Mbg = 300.0
y = np.zeros(10); y[2] = a.beta0; y[4] = -1.0; y[1] = -1.0
for _ in range(100):
    c = m.conGR(y, Mbg); d = m.dconGR(y, Mbg); y[1] -= c / d
    if abs(c) < 1e-15: break
y[6:] = m.gr(y, Mbg)
F = lambda t, y: np.exp(3*y[0]) * m.rhs(y, Mbg)
stop = lambda t, y: y[0] - a.alpha_end; stop.terminal = True
sol = solve_ivp(F, (0, 200), y, method='DOP853', rtol=1e-11, atol=1e-13, events=stop, dense_output=True, max_step=0.02)
taus = np.linspace(0, sol.t[-1], a.n); Y = sol.sol(taus).T
H = np.abs(Y[:, 1]); alpha = Y[:, 0]
beta = np.hypot(Y[:, 2], Y[:, 3]); turns = np.where(np.diff(np.sign(np.diff(Y[:, 2]))) != 0)[0]
print(f"background to alpha={alpha[-1]:.2f}: H {H[0]:.3f} -> {H[-1]:.1f}; beta+ turning points at alpha = {np.round(alpha[turns+1], 3)}")
out = {}
for M in a.Ms:
    rows = []
    for yk in Y:
        lam = np.linalg.eigvals(m.jac(yk, M))
        pos = np.where(lam.imag > 0)[0]
        if len(pos) == 0:
            rows.append((np.nan, np.nan, lam.real.max())); continue
        g = pos[np.argsort(-np.abs(lam[pos].imag))][:2]
        rows.append((lam[g].real.max(), np.abs(lam[g].imag).mean(), lam.real.max()))
    rows = np.array(rows)
    out[M] = rows
    HM = H / M
    print(f"\nM={M:g}:  H/M     omega/M   (Re lam - 1.5H)/H   max Re lam / H   alpha")
    for q in np.linspace(0, len(taus)-1, 14).astype(int):
        print(f"       {HM[q]:7.4f}  {rows[q,1]/M:8.4f}   {(rows[q,0]-1.5*H[q])/H[q]:+9.4f}      {rows[q,2]/H[q]:8.3f}      {alpha[q]:6.3f}")
np.savez(a.out, taus=taus, Y=Y, Ms=np.array(a.Ms), rows=np.array([out[M] for M in a.Ms]))
