"""Maximal Lyapunov exponent of the full 12-d Lee-Wick YM system at E=0.3, eps=0.01.
Benettin two-trajectory method, Yoshida-4 symplectic integrator for
H = 1/2|P|^2 - 1/2|W|^2 + V(q) + w.gradV(q) - 1/2|w|^2,  q = p - w.
Also records escape times (max|component|>1e4) for comparison with the paper (median ~95-106).
"""
import numpy as np, sys, time
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
E, eps = 0.3, float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
N, dt, T, ESC = 400, 0.004, 160.0, 1e4
d0, renorm_every = 1e-8, 1.0

def V(q):
    x, y, z = q.T
    return 0.5 * (x*x*y*y + y*y*z*z + z*z*x*x)
def gradV(q):
    x, y, z = q.T
    return np.stack([x*(y*y+z*z), y*(x*x+z*z), z*(x*x+y*y)], 1)
def hessV_w(q, w):
    x, y, z = q.T
    wx, wy, wz = w.T
    return np.stack([(y*y+z*z)*wx + 2*x*y*wy + 2*x*z*wz,
                     2*x*y*wx + (x*x+z*z)*wy + 2*y*z*wz,
                     2*x*z*wx + 2*y*z*wy + (x*x+y*y)*wz], 1)
def kick(S, h):
    p, w, P, W = S
    q = p - w
    g = gradV(q); Hw = hessV_w(q, w)
    P -= h * (g + Hw)
    W -= h * (-(Hw + w))
def drift(S, h):
    p, w, P, W = S
    p += h * P
    w -= h * W
c = 2**(1/3); w1 = 1/(2-c); w0 = -c/(2-c)
D = [w1/2, (w0+w1)/2, (w0+w1)/2, w1/2]; K = [w1, w0, w1]
def step(S, h):
    drift(S, D[0]*h); kick(S, K[0]*h); drift(S, D[1]*h); kick(S, K[1]*h)
    drift(S, D[2]*h); kick(S, K[2]*h); drift(S, D[3]*h)
def energy(S):
    p, w, P, W = S; q = p - w
    return 0.5*(P*P).sum(1) - 0.5*(W*W).sum(1) + V(q) + (w*gradV(q)).sum(1) - 0.5*(w*w).sum(1)

# --- initial conditions: random YM point, ghost with 1/2|w|^2+1/2|W|^2 = eps*E, total H = E by scaling YM part
q0 = rng.normal(size=(N, 3)); P0 = rng.normal(size=(N, 3))
gw = rng.normal(size=(N, 3)); gW = rng.normal(size=(N, 3))
if eps > 0:
    s = np.sqrt(2*eps*E / ((gw*gw).sum(1) + (gW*gW).sum(1)))[:, None]; w0_ = gw*s; W0_ = gW*s
else:
    w0_ = np.zeros((N, 3)); W0_ = np.zeros((N, 3))
def H_of_scale(a):
    S = [a[:, None]*q0 + w0_, w0_.copy(), a[:, None]*P0, W0_.copy()]  # p = q + w
    return energy(S)
lo, hi = np.full(N, 1e-3), np.full(N, 10.0)
for _ in range(80):
    mid = 0.5*(lo+hi); m = H_of_scale(mid) < E; lo = np.where(m, mid, lo); hi = np.where(m, hi, mid)
a = 0.5*(lo+hi)
S = [a[:, None]*q0 + w0_, w0_.copy(), a[:, None]*P0, W0_.copy()]
assert np.allclose(energy(S), E, atol=1e-6)
# burn-in of the YM part is skipped: the paper's escape times are measured from a random phase point too.

# perturbed copy
pert = rng.normal(size=(N, 12)); pert /= np.linalg.norm(pert, axis=1)[:, None]
S2 = [S[i] + d0*pert[:, 3*i:3*i+3] for i in range(4)]

active = np.ones(N, bool); t_esc = np.full(N, np.nan)
logsum = np.zeros(N); t_last = np.zeros(N)   # accumulated log growth and time over which it was accumulated
lam_hist = [[] for _ in range(N)]             # per-renorm-interval local exponents
t = 0.0; t_ren = 0.0; t0 = time.time()
while t < T and active.any():
    # sub-stepping by channel depth over the active set
    depth = max(np.abs(S[0][active]).max(), np.abs(S[0][active] - S[1][active]).max())
    nsub = int(max(1, np.ceil(depth*dt/0.25)))
    for _ in range(nsub):
        step(S, dt/nsub); step(S2, dt/nsub)
    t += dt
    # escape test
    comp = np.max([np.abs(S[i]).max(1) for i in range(4)], axis=0)
    newly = active & (comp > ESC)
    t_esc[newly] = t; active &= ~newly
    if t - t_ren >= renorm_every:
        diff = np.concatenate([S2[i]-S[i] for i in range(4)], 1)
        dist = np.linalg.norm(diff, axis=1)
        g = np.log(dist/d0)
        for i in np.where(active)[0]:
            lam_hist[i].append(g[i]/(t - t_ren))
        logsum[active] += g[active]; t_last[active] = t
        for i in range(4):
            S2[i] = S[i] + d0*diff[:, 3*i:3*i+3]/dist[:, None]
        t_ren = t
        if int(t) % 20 == 0:
            print(f"t={t:6.1f} active={active.sum():4d} {time.time()-t0:6.0f}s", flush=True)

# exponent: per trajectory, mean of local exponents excluding the last 10 time units before escape
lams = []
for i in range(N):
    h = np.array(lam_hist[i])
    if np.isfinite(t_esc[i]):
        h = h[:max(1, int(t_esc[i]) - 10)]
    else:
        h = h[5:]
    if len(h) >= 20: lams.append(h.mean())
lams = np.array(lams)
print(f"\nE={E} eps={eps} N={N}")
print(f"escaped: {np.isfinite(t_esc).sum()}/{N}; median t_esc = {np.nanmedian(t_esc):.1f}, "
      f"10-90%: {np.nanpercentile(t_esc,10):.0f}-{np.nanpercentile(t_esc,90):.0f}")
print(f"lambda_max (per-trajectory time average, excl. last 10 units): mean {lams.mean():.3f}, "
      f"median {np.median(lams):.3f}, SE {lams.std()/np.sqrt(len(lams)):.3f}, n={len(lams)}")
print(f"ensemble-average (sum log / sum t): {logsum.sum()/t_last.sum():.3f}")
np.savez(f"lyap_E{E}_eps{eps}.npz", t_esc=t_esc, lams=lams, logsum=logsum, t_last=t_last)
