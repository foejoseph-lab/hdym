"""hdym_spectrum_xy.py -- two-sided spectrum of the off-diagonal Hessian entry xy along 3-d YM
trajectories at E = 1 (Welch/Hann), against the derived tail u^2 s_xy(u) = (16 sqrt2 pi/15) nu0 = 0.672
(hdym_analytic.py Sec. 6: chirps sqrt(2 J_y |x|) cos(phi_y) at frequency |x| in the x-channel, and
the same in the y-channel).  Also the entry z^2 for reference (u^-6 tail, A_res = 6.14 resonant share).
Usage: python hdym_spectrum_xy.py --N 128 --T 3000
"""
import argparse, numpy as np
YC = (0.6756035959798288, -0.17560359597982886, -0.17560359597982886, 0.6756035959798288)
YD = (1.3512071919596576, -1.7024143839193153, 1.3512071919596576, 0.0)
def V(q):
    x, y, z = q; return 0.5 * (x*x*y*y + y*y*z*z + z*z*x*x)
def grad(q):
    x, y, z = q; return np.array([x*(y*y+z*z), y*(z*z+x*x), z*(x*x+y*y)])
def step(q, p, h):
    for c, d in zip(YC, YD):
        q = q + c*h*p
        if d != 0.0: p = p - d*h*grad(q)
    return q, p
ap = argparse.ArgumentParser(); ap.add_argument('--N', type=int, default=128); ap.add_argument('--T', type=float, default=3000.0)
ap.add_argument('--tburn', type=float, default=300.0); ap.add_argument('--dt', type=float, default=0.005); ap.add_argument('--every', type=int, default=4)
ap.add_argument('--seg', type=float, default=163.84); ap.add_argument('--seed', type=int, default=3)
a = ap.parse_args()
rng = np.random.default_rng(a.seed)
q = rng.standard_normal((3, a.N)); p = rng.standard_normal((3, a.N))
lam = (1.0/(0.5*np.sum(p*p, 0) + V(q)))**0.25; q *= lam; p *= lam**2
for _ in range(int(a.tburn/a.dt)):
    # sub-step in channels
    om = np.max(np.abs(q)); n = max(1, int(np.ceil(om*a.dt/0.25)))
    for _ in range(n): q, p = step(q, p, a.dt/n)
nrec = int(a.T/a.dt/a.every); Kxy = np.empty((a.N, nrec)); Kzz = np.empty((a.N, nrec))
for k in range(nrec):
    for _ in range(a.every):
        om = np.max(np.abs(q)); n = max(1, int(np.ceil(om*a.dt/0.25)))
        for _ in range(n): q, p = step(q, p, a.dt/n)
    Kxy[:, k] = q[0]*q[1]; Kzz[:, k] = q[0]**2 + q[1]**2
dts = a.dt*a.every; M = int(round(a.seg/dts)); win = np.hanning(M); Teff = dts*np.sum(win**2); nseg = nrec//M
u = 2*np.pi*np.fft.rfftfreq(M, dts)
def spec(K):
    S = np.zeros(len(u))
    for i in range(a.N):
        for j in range(nseg):
            s = K[i, j*M:(j+1)*M]; s = s - s.mean(); X = dts*np.fft.rfft(win*s); S += np.abs(X)**2/Teff
    return S/(a.N*nseg)
Sxy = spec(Kxy); Szz = spec(Kzz)
nu0 = 8*np.pi**2/557.1; Axy = 16*np.sqrt(2)*np.pi/15*nu0; Ares = 32*np.pi*(64/(105*np.sqrt(2)))*nu0
print(f"E1 err {np.max(np.abs(0.5*np.sum(p*p,0)+V(q)-1)):.1e};  predicted u^2 s_xy = {Axy:.3f};  A_res(zz) = {Ares:.2f}, total 7/2 A_res = {3.5*Ares:.1f}")
print("   u-band    u^2 s_xy   /pred     u^6 s_zz")
for lo, hi in ((2,3),(3,4),(4,5),(5,6),(6,8),(8,10),(10,12),(12,16),(16,20),(20,30),(30,45)):
    m = (u >= lo) & (u < hi)
    print(f"  {lo:3d}-{hi:3d}   {np.mean(u[m]**2*Sxy[m]):8.3f}   {np.mean(u[m]**2*Sxy[m])/Axy:6.3f}    {np.mean(u[m]**6*Szz[m]):8.2f}")
np.savez("spec_xy_E1.npz", u=u, Sxy=Sxy, Szz=Szz)
