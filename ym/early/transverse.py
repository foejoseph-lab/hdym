"""Linearised ghost on pure YM trajectories: transverse Lyapunov exponent (vector and diagonal ghost),
tachyonic fraction, and (at E=1) the spectrum of the Hessian diagonal entry K_xx=y^2+z^2.
usage: python3 transverse.py E T N
"""
import numpy as np, sys, time
E, T, N = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(7); dt = 0.02
def gradV(q):
    x, y, z = q.T; return np.stack([x*(y*y+z*z), y*(x*x+z*z), z*(x*x+y*y)], 1)
def V(q):
    x, y, z = q.T; return 0.5*(x*x*y*y + y*y*z*z + z*z*x*x)
def hess(q):
    x, y, z = q.T
    H = np.empty((len(q), 3, 3))
    H[:,0,0] = y*y+z*z; H[:,1,1] = x*x+z*z; H[:,2,2] = x*x+y*y
    H[:,0,1] = H[:,1,0] = 2*x*y; H[:,1,2] = H[:,2,1] = 2*y*z; H[:,0,2] = H[:,2,0] = 2*x*z
    return H
c = 2**(1/3); w1 = 1/(2-c); w0 = -c/(2-c)
Dc = [w1/2, (w0+w1)/2, (w0+w1)/2, w1/2]; Kc = [w1, w0, w1]
def step(q, p, wv, wd, h):
    # q,p: YM; wv=(w,v) vector ghost; wd=(w,v) diagonal ghost. Ghost force evaluated at current q.
    for s in range(4):
        q += Dc[s]*h*p; wv[0] += Dc[s]*h*wv[1]; wd[0] += Dc[s]*h*wd[1]
        if s < 3:
            H = hess(q); p -= Kc[s]*h*gradV(q)
            wv[1] -= Kc[s]*h*(wv[0] + np.einsum('nij,nj->ni', H, wv[0]))
            d = np.stack([H[:,0,0], H[:,1,1], H[:,2,2]], 1)
            wd[1] -= Kc[s]*h*(wd[0] + d*wd[0])
q = rng.normal(size=(N, 3)); p = rng.normal(size=(N, 3))
lo, hi = np.full(N, 1e-3), np.full(N, 10.0)
for _ in range(80):
    m = 0.5*(lo+hi); e = 0.5*((m[:,None]*p)**2).sum(1) + V(m[:,None]*q)
    lo = np.where(e < E, m, lo); hi = np.where(e < E, hi, m)
a = 0.5*(lo+hi); q *= a[:, None]; p *= a[:, None]
wv = [rng.normal(size=(N,3)), rng.normal(size=(N,3))]; wd = [rng.normal(size=(N,3)), rng.normal(size=(N,3))]
burn = 300*E**-0.25
for _ in range(int(burn/dt)):
    depth = np.abs(q).max(); n = int(max(1, np.ceil(depth*dt/0.25)))
    for _ in range(n): step(q, p, wv, wd, dt/n)
def renorm(w):
    nrm = np.sqrt((w[0]**2).sum(1) + (w[1]**2).sum(1)); w[0] /= nrm[:,None]; w[1] /= nrm[:,None]; return np.log(nrm)
renorm(wv); renorm(wd)
nsteps = int(T/dt); lv = np.zeros(N); ld = np.zeros(N); lv_half = np.zeros(N); ld_half = np.zeros(N)
tach = 0; tach_n = 0; rec = [] if E == 1.0 else None; t0 = time.time()
for k in range(1, nsteps+1):
    depth = np.abs(q).max(); n = int(max(1, np.ceil(depth*dt/0.25)))
    for _ in range(n): step(q, p, wv, wd, dt/n)
    if k % 50 == 0:
        gv = renorm(wv); gd = renorm(wd); lv += gv; ld += gd
        if k > nsteps//2: lv_half += gv; ld_half += gd
    if k % 10 == 0:
        ev = np.linalg.eigvalsh(hess(q))[:,0]; tach += np.sum(ev <= -1.0); tach_n += N
    if rec is not None: rec.append(q[:, 1]**2 + q[:, 2]**2)
    if k % (nsteps//5) == 0:
        print(f"E={E} t={k*dt:8.0f} lam_v={lv.mean()/(k*dt):.3e} lam_d={ld.mean()/(k*dt):.3e} {time.time()-t0:.0f}s", flush=True)
Eerr = abs(0.5*(p*p).sum(1)+V(q)-E).max()/E
res = dict(E=E, T=T, N=N, lam_v=lv.mean()/T, lam_v_se=lv.std()/T/np.sqrt(N), lam_d=ld.mean()/T, lam_d_se=ld.std()/T/np.sqrt(N),
           lam_v_half=lv_half.mean()/(T/2), lam_d_half=ld_half.mean()/(T/2), tach_frac=tach/tach_n, Eerr=Eerr)
print(res)
np.savez(f"transverse_E{E}.npz", **res)
if rec is not None:
    K = np.array(rec)  # (steps, N)
    K -= K.mean(0)
    F = np.fft.rfft(K, axis=0); S = 2*(np.abs(F)**2).mean(1)*dt/len(K)  # two-sided spectrum estimate, per unit time
    u = 2*np.pi*np.fft.rfftfreq(len(K), dt)
    np.savez("spectrum_E1.npz", u=u, S=S)
    for lo_, hi_ in [(6,12),(12,20),(20,30),(6,30)]:
        m = (u>=lo_)&(u<=hi_); sl = np.polyfit(np.log(u[m]), np.log(S[m]), 1)[0]
        print(f"spectral tail exponent over u={lo_}-{hi_}: {sl:.2f}")
    for uu in [2.94, 3.56, 5.32, 5.72]:
        print(f"  S({uu}) = {np.interp(uu, u, S):.3e}")
