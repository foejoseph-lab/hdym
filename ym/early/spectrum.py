"""Spectrum of the Hessian diagonal entry K_xx=y^2+z^2 along pure YM trajectories at E=1:
raw, Hann-windowed, and conditioned on channel occupancy. Compares to the resonant prediction 32 pi I nu0 / u^6."""
import numpy as np, time
rng = np.random.default_rng(11); N, dt, T = 256, 0.02, 800.0
def gradV(q):
    x, y, z = q.T; return np.stack([x*(y*y+z*z), y*(x*x+z*z), z*(x*x+y*y)], 1)
def V(q):
    x, y, z = q.T; return 0.5*(x*x*y*y + y*y*z*z + z*z*x*x)
c = 2**(1/3); w1 = 1/(2-c); w0 = -c/(2-c); Dc = [w1/2, (w0+w1)/2, (w0+w1)/2, w1/2]; Kc = [w1, w0, w1]
def step(q, p, h):
    q += Dc[0]*h*p; p -= Kc[0]*h*gradV(q); q += Dc[1]*h*p; p -= Kc[1]*h*gradV(q); q += Dc[2]*h*p; p -= Kc[2]*h*gradV(q); q += Dc[3]*h*p
q = rng.normal(size=(N, 3)); p = rng.normal(size=(N, 3))
lo, hi = np.full(N, 1e-3), np.full(N, 10.0)
for _ in range(80):
    m = 0.5*(lo+hi); e = 0.5*((m[:,None]*p)**2).sum(1) + V(m[:,None]*q); lo = np.where(e < 1, m, lo); hi = np.where(e < 1, hi, m)
a = 0.5*(lo+hi); q *= a[:, None]; p *= a[:, None]
for _ in range(int(300/dt)):
    n = int(max(1, np.ceil(np.abs(q).max()*dt/0.25)))
    for _ in range(n): step(q, p, dt/n)
n_s = int(T/dt); Q = np.empty((n_s, N, 3), np.float32)
for k in range(n_s):
    n = int(max(1, np.ceil(np.abs(q).max()*dt/0.25)))
    for _ in range(n): step(q, p, dt/n)
    Q[k] = q
x, y, z = Q[...,0].astype(float), Q[...,1].astype(float), Q[...,2].astype(float)
K = y*y + z*z
u = 2*np.pi*np.fft.rfftfreq(n_s, dt)
def spec(sig, win=None):
    s = sig - sig.mean(0)
    if win is not None: s = s*win[:, None]; norm = (win**2).sum()*dt
    else: norm = n_s*dt
    F = np.fft.rfft(s, axis=0); return (np.abs(F)**2).mean(1)*dt*dt/norm
hann = np.hanning(n_s)
ax = np.abs(x); ay = np.abs(y); az = np.abs(z)
in_x = (ax > 3) & (ay < ax) & (az < ax); in_yz = ((ay > 3) & (ax < ay) & (az < ay)) | ((az > 3) & (ax < az) & (ay < az))
def cond(sig, mask):   # keep signal only on mask, mean-subtract within mask, zero elsewhere (spectrum per unit total time)
    s = np.where(mask, sig, np.nan); mu = np.nanmean(s, 0); s = np.where(mask, sig - mu, 0.0); return s
I = 64/(105*np.sqrt(2)); nu0 = 0.1417; pred = 32*np.pi*I*nu0/u[1:]**6
res = {'raw': spec(K), 'hann': spec(K, hann), 'not_yz': spec(cond(K, ~in_yz), hann), 'x_only': spec(cond(K, in_x), hann)}
print(f"fraction of time in x-channel {in_x.mean():.3f}, in y/z channel {in_yz.mean():.3f}")
for name, S in res.items():
    line = f"{name:8s}:"
    for lo_, hi_ in [(6,12),(12,20),(20,30)]:
        m = (u>=lo_)&(u<=hi_); line += f"  slope[{lo_}-{hi_}]={np.polyfit(np.log(u[m]), np.log(S[m]), 1)[0]:+.2f}"
    for uu in [2.94, 3.56, 5.32, 5.72, 8, 12]:
        line += f"  S({uu})={np.interp(uu, u, S):.2e}"
    print(line)
print("resonant prediction 32 pi I nu0/u^6:", "  ".join(f"S({uu})={32*np.pi*I*nu0/uu**6:.2e}" for uu in [2.94,3.56,5.32,5.72,8,12]))
np.savez('spectrum_E1_full.npz', u=u, **res, in_x=in_x.mean(), in_yz=in_yz.mean())
