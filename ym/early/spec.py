# Power spectrum of Hess V(p(t)) along pure-YM trajectories at E=1 (scaled units).
import numpy as np, time, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))   # ym/hdym_core.py
from hdym_core import ym_initial, gradV
from scipy.signal.windows import blackmanharris
rng=np.random.default_rng(3)
N=512; dt=0.005; sample=10; nseg=8192   # sampling dt_s=0.05 -> Nyquist 62.8 rad
E=1.0
p,P=ym_initial(N,E,rng)
C=(0.6756035959798288,-0.17560359597982886,-0.17560359597982886,0.6756035959798288)
D=(1.3512071919596576,-1.7024143839193153,1.3512071919596576,0.0)
def step(p,P,h):
    for c,d in zip(C,D):
        p=p+c*h*P
        if d: P=P-d*h*gradV(p)
    return p,P
def sstep(p,P):
    om=np.max(np.abs(p),1); n=np.maximum(1,np.ceil(om*dt/0.1)).astype(int)
    if np.all(n==1): return step(p,P,dt)
    p,P=p.copy(),P.copy()
    for k in np.unique(n):
        i=n==k; a,b=p[i],P[i]
        for _ in range(k): a,b=step(a,b,dt/k)
        p[i],P[i]=a,b
    return p,P
for _ in range(int(200/dt)): p,P=sstep(p,P)     # burn-in
E0=0.5*(P**2).sum(1)+0.5*((p[:,0]*p[:,1])**2+(p[:,1]*p[:,2])**2+(p[:,2]*p[:,0])**2)
nseg_total=4
win=blackmanharris(nseg); U=(win**2).sum()
Sdiag=np.zeros(nseg//2+1); Soff=np.zeros(nseg//2+1)
t0=time.time()
for s in range(nseg_total):
    buf=np.empty((nseg,N,6))
    for j in range(nseg):
        for _ in range(sample): p,P=sstep(p,P)
        x,y,z=p.T
        buf[j]=np.stack([y*y+z*z,z*z+x*x,x*x+y*y,2*x*y,2*x*z,2*y*z],1)
    buf-=buf.mean(0)
    F=np.fft.rfft(buf*win[:,None,None],axis=0)
    S=(np.abs(F)**2)/U
    Sdiag+=S[:,:,:3].mean((1,2)); Soff+=S[:,:,3:].mean((1,2))
    print("seg",s,time.time()-t0,flush=True)
E1=0.5*(P**2).sum(1)+0.5*((p[:,0]*p[:,1])**2+(p[:,1]*p[:,2])**2+(p[:,2]*p[:,0])**2)
print("max rel E err",np.max(abs(E1-E0)/E0))
om=2*np.pi*np.fft.rfftfreq(nseg,dt*sample)
np.savez("spec.npz",om=om,Sd=Sdiag/nseg_total,So=Soff/nseg_total)
