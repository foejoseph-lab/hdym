"""Pure Yang-Mills mechanics at E=1: channel-visit statistics (re-derivation of Fig. 2 / Sec. 5.3).
A visit to the x-channel starts when |x| crosses D_in upward while |y|,|z|<|x|, ends when |x| falls below D_in.
Records D_max and the transverse action J = J_y+J_z at entry, J_i = (1/2 v_i^2 + 1/2 x^2 q_i^2)/|x|.
"""
import numpy as np, time
rng = np.random.default_rng(2)
N, dt, T, D_in = 512, 0.01, 2900.0, 3.0
def gradV(q):
    x, y, z = q.T
    return np.stack([x*(y*y+z*z), y*(x*x+z*z), z*(x*x+y*y)], 1)
def V(q):
    x, y, z = q.T; return 0.5*(x*x*y*y + y*y*z*z + z*z*x*x)
c = 2**(1/3); w1 = 1/(2-c); w0 = -c/(2-c)
D = [w1/2, (w0+w1)/2, (w0+w1)/2, w1/2]; K = [w1, w0, w1]
def step(q, p, h):
    q += D[0]*h*p; p -= K[0]*h*gradV(q); q += D[1]*h*p; p -= K[1]*h*gradV(q)
    q += D[2]*h*p; p -= K[2]*h*gradV(q); q += D[3]*h*p
# random phase points scaled to E=1
q = rng.normal(size=(N, 3)); p = rng.normal(size=(N, 3))
lo, hi = np.full(N, 1e-3), np.full(N, 10.0)
for _ in range(80):
    m = 0.5*(lo+hi); e = 0.5*((m[:,None]*p)**2).sum(1) + V(m[:,None]*q)
    lo = np.where(e < 1, m, lo); hi = np.where(e < 1, hi, m)
a = 0.5*(lo+hi); q *= a[:, None]; p *= a[:, None]
# burn-in
for _ in range(int(300/dt)): step(q, p, dt)
in_visit = np.zeros((N, 3), bool); dmax = np.zeros((N, 3)); Jent = np.zeros((N, 3))
Dlist, Jlist = [], []
t = 0.0; t0 = time.time(); nsteps = int(T/dt)
for k in range(nsteps):
    depth = np.abs(q).max()
    nsub = int(max(1, np.ceil(depth*dt/0.25)))
    for _ in range(nsub): step(q, p, dt/nsub)
    absq = np.abs(q)
    for i in range(3):
        j, l = (i+1) % 3, (i+2) % 3
        deep = (absq[:, i] > D_in) & (absq[:, j] < absq[:, i]) & (absq[:, l] < absq[:, i])
        enter = deep & ~in_visit[:, i]
        if enter.any():
            X = absq[enter, i]
            J = (0.5*p[enter, j]**2 + 0.5*q[enter, i]**2*q[enter, j]**2)/X + \
                (0.5*p[enter, l]**2 + 0.5*q[enter, i]**2*q[enter, l]**2)/X
            Jent[enter, i] = J; dmax[enter, i] = X; in_visit[enter, i] = True
        cont = deep & in_visit[:, i]
        dmax[cont, i] = np.maximum(dmax[cont, i], absq[cont, i])
        leave = in_visit[:, i] & ~deep
        if leave.any():
            Dlist.extend(dmax[leave, i]); Jlist.extend(Jent[leave, i]); in_visit[leave, i] = False
    if k % (nsteps//10) == 0:
        print(f"t={k*dt:6.0f}  visits={len(Dlist)}  E err={abs(0.5*(p*p).sum(1)+V(q)-1).max():.1e}  {time.time()-t0:.0f}s", flush=True)
Dm = np.array(Dlist); J = np.array(Jlist); Ttot = N*T
np.savez('nu0_E1_rederived.npz', D=Dm, J=J, Ttot=Ttot)
print(f"\nvisits: {len(Dm)}  (paper: 35585 for D_in unspecified)")
print("adiabatic relation D_max*J: median %.3f, 10-90%%: %.3f-%.3f" % (np.median(Dm*J), *np.percentile(Dm*J, [10, 90])))
print("rate of visits with D_max > D, and implied nu0 = 2 D^2 rate / 3:")
for Dc in [4, 6, 8, 12, 16, 24, 32, 48]:
    r = np.sum(Dm > Dc)/Ttot
    print(f"  D>{Dc:3d}: rate {r:.3e}   nu0 {2*Dc*Dc*r/3:.4f}   (pred 0.1417)")
