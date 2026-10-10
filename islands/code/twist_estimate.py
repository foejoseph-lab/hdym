"""
twist_estimate.py -- NON-RIGOROUS estimate of the Birkhoff twist of the DR orbit at E = 1/2.

Return map P on the section {y = 0, py > 0, H = E}, coordinates (x, px): area-preserving.
Fixed point z* (certified).  Linear part DP has eigenvalues e^{+-i theta}, cos theta = HI.
Choose real symplectic S (det S = 1) with S^{-1} DP S = rotation by theta; in zeta = xi + i eta
the map is  zeta -> e^{i theta} zeta + O(2),  and the Birkhoff normal form is
      zeta -> exp(i (theta + tau |zeta|^2 + ...)) zeta .
Here |zeta|^2 / 2 = action.  Estimate tau from the slope of the mean angular increment per
iterate against the mean |zeta|^2 along quasi-periodic orbits inside the island.
"""
import numpy as np
from scipy.integrate import solve_ivp

E = 0.5
xs, pxs = 3.14640122769837, 0.00179319341881          # certified fixed point
T = 31.053330268585066

def f(t, s):
    x, y, px, py = s
    return [px, py, -x*y*y, -x*x*y]

def lift(x, px):
    return [x, 0.0, px, np.sqrt(2*E - px*px)]

def section_orbit(x, px, n_cross, rtol=1e-11):
    """n_cross successive crossings of y=0 with py>0 (excluding the initial point)."""
    ev = lambda t, s: s[1]; ev.direction = 1
    s0 = lift(x, px); out = []
    t_end = (n_cross + 1) * T / 9 * 1.2       # 9 upward crossings per period
    sol = solve_ivp(f, (0, t_end), s0, rtol=rtol, atol=1e-14, events=ev, dense_output=False)
    pts = sol.y_events[0]
    pts = pts[sol.t_events[0] > 1e-9]
    return pts[:n_cross, [0, 2]]              # (x, px) at crossings

def return_map(x, px):
    """P = 9th crossing (one period of the DR orbit)."""
    return section_orbit(x, px, 9)[-1]

# ---- linear part by central differences
d = 1e-6
DP = np.zeros((2, 2))
for j, dv in enumerate([(d, 0), (0, d)]):
    zp = return_map(xs + dv[0], pxs + dv[1]); zm = return_map(xs - dv[0], pxs - dv[1])
    DP[:, j] = (zp - zm) / (2*d)
print("DP =\n", DP, "\ndet DP =", np.linalg.det(DP), "  HI = tr/2 =", np.trace(DP)/2)
theta = np.arccos(np.trace(DP)/2)
# symplectic normalizing frame: eigenvector v = u + i w of DP for e^{i theta}; omega(u,w) = u x w
lam, V = np.linalg.eig(DP)
k = np.argmin(abs(lam - np.exp(1j*theta)))
v = V[:, k]; u, w = v.real, v.imag
area = u[0]*w[1] - u[1]*w[0]
if area < 0: w = -w; area = -area; theta = -theta           # orient so omega(u,w) > 0
S = np.column_stack([u, w]) / np.sqrt(area)                   # det S = 1
Sinv = np.linalg.inv(S)
print("theta =", theta, " rotation number theta/2pi =", theta/(2*np.pi), " check S^-1 DP S =\n", Sinv @ DP @ S)

# ---- rotation number vs amplitude
print("\n  d (x-offset)   mean |zeta|^2    mean dphi      dphi - theta    regular?")
rows = []
for d0 in [2e-5, 5e-5, 1e-4, 1.5e-4, 2.5e-4]:
    n_iter = 120
    pts = section_orbit(xs + d0, pxs, 9*n_iter)[8::9]          # every 9th crossing = one P-iterate
    zeta = (Sinv @ (pts - [xs, pxs]).T)
    z = zeta[0] + 1j*zeta[1]
    r2 = np.mean(abs(z)**2)
    dphi = np.angle(z[1:] / z[:-1])
    regular = np.std(abs(z)) / np.mean(abs(z)) < 0.3
    rows.append((d0, r2, dphi.mean()))
    print(f"  {d0:9.1e}   {r2:12.4e}   {dphi.mean():12.8f}   {dphi.mean()-theta:+.3e}   {regular}  (|z| spread {np.std(abs(z))/np.mean(abs(z)):.2f})")
rows = np.array(rows)
good = rows[:, 1] < 1e-7
if good.sum() >= 2:
    tau = np.polyfit(rows[good, 1], rows[good, 2] - theta, 1)[0]
    print(f"\ntwist estimate  tau = d(dphi)/d|zeta|^2  ~ {tau:.4g}   (action units: |zeta|^2 = 2 I)")
