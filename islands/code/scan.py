import numpy as np, sys
from scipy.integrate import solve_ivp
x0, vx0, y0, vy0 = 3.14640122769753, 0.0017931934191, 0.0, 0.99999839222738
E0 = 0.5; T = 31.053330268627366

def rate(E):
    lam = (E/E0)**0.25; TE = T/lam
    def rhs(t, s):
        x, y, vx, vy = s[:4]
        H = np.array([[y*y, 2*x*y],[2*x*y, x*x]])
        A = np.zeros((4,4)); A[0,2]=A[1,3]=1; A[2:,:2] = -(np.eye(2)+H)
        P = s[4:20].reshape(4,4)
        return [vx, vy, -x*y*y, -y*x*x] + list((A@P).ravel())
    s0 = [lam*x0, lam*y0, lam**2*vx0, lam**2*vy0] + list(np.eye(4).ravel())
    r = solve_ivp(rhs, (0,TE), s0, rtol=1e-11, atol=1e-13)
    P = r.y[4:20,-1].reshape(4,4)
    ev = np.linalg.eigvals(P)
    # Floquet multipliers for 4x4 symplectic: two pairs; report traces/2 of each pair via eigs
    return np.log(max(abs(ev)))/TE, ev

a, b, n = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
for E in np.linspace(a, b, n):
    g, ev = rate(E)
    flag = "UNSTABLE" if g > 1e-6 else ""
    print(f"{E:9.5f}  {g:11.4e}  {flag}")
