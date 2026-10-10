"""Nonlinear escape at E=0.3: escape times at eps=1e-4 (compare eps=1e-2 median ~95-106) and the
mechanism near blow-up: transverse displacement rho off the channel axis and 1+lambda_min(Hess V) there."""
import numpy as np, sys, time, os
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lyap_full.py')).read().split('# --- initial conditions')[0]; src=src.replace('int(sys.argv[1]) if len(sys.argv) > 1 else 1','3').replace('float(sys.argv[2]) if len(sys.argv) > 2 else 0.01','float(sys.argv[1])'); exec(src)
N, dt, T = 400, 0.004, 400.0
q0 = rng.normal(size=(N, 3)); P0 = rng.normal(size=(N, 3)); gw = rng.normal(size=(N, 3)); gW = rng.normal(size=(N, 3))
s = np.sqrt(2*eps*E / ((gw*gw).sum(1) + (gW*gW).sum(1)))[:, None]; w0_ = gw*s; W0_ = gW*s
def H_of_scale(a):
    return energy([a[:, None]*q0 + w0_, w0_.copy(), a[:, None]*P0, W0_.copy()])
lo, hi = np.full(N, 1e-3), np.full(N, 10.0)
for _ in range(80):
    mid = 0.5*(lo+hi); m = H_of_scale(mid) < E; lo = np.where(m, mid, lo); hi = np.where(m, hi, mid)
a = 0.5*(lo+hi); S = [a[:, None]*q0 + w0_, w0_.copy(), a[:, None]*P0, W0_.copy()]
def lam_min_hess(q):
    x, y, z = q.T
    H = np.empty((len(q), 3, 3))
    H[:,0,0] = y*y+z*z; H[:,1,1] = x*x+z*z; H[:,2,2] = x*x+y*y
    H[:,0,1] = H[:,1,0] = 2*x*y; H[:,1,2] = H[:,2,1] = 2*y*z; H[:,0,2] = H[:,2,0] = 2*x*z
    return np.linalg.eigvalsh(H)[:, 0]
active = np.ones(N, bool); t_esc = np.full(N, np.nan); t = 0.0
snaps = {10: {}, 100: {}, 1000: {}}   # at first crossing of max|component| = level: (rho, 1+lam_min, |w|, |q|)
t0 = time.time()
while t < T and active.any():
    depth = np.abs(S[0][active]).max(); nsub = int(max(1, np.ceil(depth*dt/0.25)))
    for _ in range(nsub): step(S, dt/nsub)
    t += dt
    comp = np.max([np.abs(S[i]).max(1) for i in range(4)], axis=0)
    q = S[0] - S[1]
    for lev in snaps:
        hit = active & (comp > lev)
        for i in np.where(hit)[0]:
            if i not in snaps[lev]:
                ax = np.argmax(np.abs(q[i])); rho = np.sqrt((q[i]**2).sum() - q[i][ax]**2)
                snaps[lev][i] = (rho, 1+lam_min_hess(q[i:i+1])[0], np.linalg.norm(S[1][i]), np.linalg.norm(q[i]))
    newly = active & (comp > 1e4); t_esc[newly] = t; active &= ~newly
print(f"eps={eps}: escaped {np.isfinite(t_esc).sum()}/{N} by T={T}; median t_esc {np.nanmedian(t_esc):.1f}, "
      f"10-90% {np.nanpercentile(t_esc,10):.0f}-{np.nanpercentile(t_esc,90):.0f}  ({time.time()-t0:.0f}s)")
for lev in snaps:
    v = np.array(list(snaps[lev].values()))
    if len(v):
        print(f"at max|component|={lev}: n={len(v)}  rho median {np.median(v[:,0]):.2f}  1+lam_min median {np.median(v[:,1]):.1f} "
              f"(10-90%: {np.percentile(v[:,1],10):.1f} to {np.percentile(v[:,1],90):.1f})  |w| median {np.median(v[:,2]):.2g}  |q| median {np.median(v[:,3]):.2g}")
np.savez(f'escape_eps{eps}.npz', t_esc=t_esc)
