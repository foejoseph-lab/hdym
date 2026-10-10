import numpy as np, sys
from scipy.integrate import solve_ivp
x0, vx0, y0, vy0 = 3.14640122769753, 0.0017931934191, 0.0, 0.99999839222738
E0 = 0.5

def run(E, dx, Ttot=3000.0):
    lam = (E/E0)**0.25
    def rhs(t, s):
        x, y, vx, vy = s[:4]
        H = np.array([[y*y, 2*x*y],[2*x*y, x*x]])
        # base tangent (for base Lyapunov) and ghost vector (wx,wy)
        dq = s[4:6]; dp = s[6:8]
        w = s[8:10]; wd = s[10:12]
        return [vx, vy, -x*y*y, -y*x*x,
                dp[0], dp[1], *(-(H@dq)),
                wd[0], wd[1], *(-(w + H@w))]
    # perturb x0 and fix energy by rescaling vy
    xs = lam*(x0+dx); ys = 0.0; vxs = lam**2*vx0
    vys = np.sqrt(2*E - vxs**2 - xs**2*ys**2)
    s = np.array([xs, ys, vxs, vys, 1,0,0,0, 1,0,0,0], float)
    s[4:8] /= np.linalg.norm(s[4:8]); s[8:12] /= np.linalg.norm(s[8:12])
    lb = lg = 0.0; t = 0.0; dt = 10.0
    while t < Ttot:
        r = solve_ivp(rhs, (t, t+dt), s, rtol=1e-10, atol=1e-12)
        s = r.y[:,-1]; t += dt
        nb = np.linalg.norm(s[4:8]); ng = np.linalg.norm(s[8:12])
        lb += np.log(nb); lg += np.log(ng); s[4:8] /= nb; s[8:12] /= ng
    return lb/Ttot, lg/Ttot

for E in [0.2, 0.3, 1.0]:
    for dx in [0.0, 2e-4, 5e-4, 1e-3, 3e-3]:
        lb, lg = run(E, dx)
        print(f"E={E:5.2f} dx={dx:7.1e}  base lambda={lb:9.2e}  ghost lambda={lg:9.3e}  {'(regular)' if lb<2e-3 else '(chaotic)'}")
    print()
