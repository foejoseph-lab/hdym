"""Measure nu0: rate density of channel entries per unit transverse-action area at E=1.
Pure YM, sub-stepped.  A visit to channel i starts when |p_i| crosses A_in upward (with
|p_j|,|p_k| < A_in) and ends when it drops back below.  Records D_max and the transverse
action J = E_perp/|p_i| at entry.  Predictions (uniform action measure):
   rate(any channel, D_max > D) = 3 nu0 / (2 D^2)      (from J = 1/D_max)
   rate density in J = J_y + J_z  :  3 nu0 J  per unit time  (any channel)
"""
import sys, numpy as np
sys.path.insert(0, '.')
from hdym_core import ym_initial, ym_step_substepped, V
N, E, dt, T, tburn, A_in = 512, 1.0, 0.01, 3000.0, 100.0, 3.0
rng = np.random.default_rng(5)
p, P = ym_initial(N, E, rng)
nsteps = int(T/dt); nburn = int(tburn/dt)
inch = -np.ones(N, int); dmax = np.zeros(N); Jent = np.zeros(N)
D_list, J_list = [], []
for k in range(nsteps):
    p, P = ym_step_substepped(p, P, dt)
    if k < nburn: continue
    ap = np.abs(p); i = np.argmax(ap, 1); d = ap[np.arange(N), i]
    # entries
    ent = (inch < 0) & (d > A_in)
    if ent.any():
        idx = np.nonzero(ent)[0]; ii = i[idx]
        x = p[idx, ii]
        m = np.ones((idx.size, 3), bool); m[np.arange(idx.size), ii] = False
        pperp = p[idx][m].reshape(-1, 2); Pperp = P[idx][m].reshape(-1, 2)
        Eperp = 0.5*np.sum(Pperp**2, 1) + 0.5*x*x*np.sum(pperp**2, 1)
        Jent[idx] = Eperp/np.abs(x); inch[idx] = ii; dmax[idx] = d[idx]
    inv = inch >= 0
    dmax[inv] = np.maximum(dmax[inv], d[inv])
    ex = inv & (d < A_in)
    if ex.any():
        D_list.extend(dmax[ex]); J_list.extend(Jent[ex]); inch[ex] = -1
Ttot = (nsteps-nburn)*dt*N
D = np.array(D_list); J = np.array(J_list)
print(f"{len(D)} visits in {Ttot:.0f} trajectory-time units; visit rate {len(D)/Ttot:.4f}/unit time")
print("rate(any channel, Dmax > D) * D^2  -> should be 3 nu0/2 (constant):")
for Dc in (4, 6, 8, 12, 16, 24, 32, 48, 64):
    r = np.sum(D > Dc)/Ttot
    print(f"  D={Dc:3d}  rate={r:.2e}  rate*D^2={r*Dc*Dc:.3f}  nu0={2*r*Dc*Dc/3:.3f}   (n={np.sum(D>Dc)})")
print("from J at entry: rate(J < j)/j^2 -> should be 3 nu0/2:")
for jc in (0.01, 0.02, 0.05, 0.1, 0.2):
    r = np.sum(J < jc)/Ttot
    print(f"  j={jc:.2f}  rate={r:.2e}  rate/j^2={r/jc/jc:.3f}  nu0={2*r/jc/jc/3:.3f}  (n={np.sum(J<jc)})")
good = J > 0
print(f"check Dmax*J (should be ~1 - E_perp ~ 1): median {np.median(D[good]*J[good]):.3f}")
np.savez("nu0_E1.npz", D=D, J=J, Ttot=Ttot)
