import numpy as np, sympy as sp, pickle, time
exec(open('bianchi9_kasner.py').read().split("from scipy.integrate import solve_ivp\nimport matplotlib")[0])
from scipy.integrate import solve_ivp
theta = np.pi
# frozen spectrum of Jlin(x): eigenvalues scaled by x (so Euler-type exponents appear as x*lambda)
print("frozen spectrum: x, H/M, eigenvalues*x (sorted by |Im|)")
for xx in (30, 10, 3, 1.5, 1.0, 0.7, 0.5, 0.3, 0.1, 0.03):
    lam = np.linalg.eigvals(np.array(f_J(xx, theta), float)) * xx
    lam = lam[np.argsort(-np.abs(lam.imag))]
    print(f"  x={xx:5.2f} H/M={1/(3*xx):6.3f}  ", "  ".join(f"{l.real:+.2f}{l.imag:+.2f}i" for l in lam[:4]), " | GR-ish:", "  ".join(f"{l.real:+.2f}{l.imag:+.2f}i" for l in lam[4:]))
# clean ghost preparation: integrate from x_hi = 200 and read growth of each component
x_hi, x_lo = 200.0, 0.02
J0 = np.array(f_J(x_hi, theta), float); lam, V = np.linalg.eig(J0)
ghost = np.where(np.abs(lam.imag) > 0.5)[0]
v0 = V[:, ghost[0]].real; v0 /= np.linalg.norm(v0)
F = lambda xx, Y: np.array(f_J(xx, theta), float) @ Y
xs = np.geomspace(x_hi, x_lo, 600)
sol = solve_ivp(F, (x_hi, x_lo), v0, t_eval=xs, method='DOP853', rtol=1e-11, atol=1e-14)
Y = sol.y
names = ['da','da1','bp','bm','bp1','bm1','bp2','bm2','bp3','bm3']
print("\nsmall-x power laws per component (fit over x in [0.02, 0.1]):")
k = xs < 0.1
for i, nm in enumerate(names):
    yy = np.abs(Y[i, k]); 
    if yy.max() < 1e-12: continue
    sl = np.polyfit(np.log(xs[k]), np.log(yy + 1e-300), 1)[0]
    print(f"  {nm:4s} ~ x^{sl:+.2f}")
# envelope of beta+ perturbation vs H/M, normalised to adiabatic t^{-1/2}
env = np.abs(Y[2]) * np.sqrt(xs)
import numpy as np
for hm in (0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 1, 3):
    j = np.argmin(np.abs(1/(3*xs) - hm)); w = slice(max(j-15,0), j+15)
    print(f"  H/M={hm:4.2f}: |delta beta| sqrt(x) envelope (max over window) = {env[w].max():.3e}")

print("\ngeneric and GR-mode initial data, growth exponents of |delta beta+| at small x:")
rng = np.random.default_rng(0)
def expo(v0, x_hi=3.0):
    xs = np.geomspace(x_hi, 0.01, 400)
    sol = solve_ivp(F, (x_hi, 0.01), v0, t_eval=xs, method='DOP853', rtol=1e-11, atol=1e-14)
    k = xs < 0.1
    sl = [np.polyfit(np.log(xs[k]), np.log(np.abs(sol.y[i, k]) + 1e-300), 1)[0] for i in (2, 3, 0)]
    return sl, np.abs(sol.y[2, -1]) / np.abs(sol.y[2, 0])
for trial in range(3):
    v = rng.normal(size=10); print(f"  random #{trial}: (bp, bm, da) ~ x^{expo(v)[0]}")
# GR time-shift mode: delta y = f_GR(y) dt -> components (alpha', alpha'', beta', beta'',...) of the background
# background: alpha = ln t/3 -> alpha' = 1/(3t), alpha'' = -1/(3t^2); beta+ = cos th ln t/3 etc.
tt = 3.0; cth, sth = np.cos(theta), np.sin(theta)
v_shift = np.array([1/(3*tt), -1/(3*tt**2), cth/(3*tt), sth/(3*tt), -cth/(3*tt**2), -sth/(3*tt**2), 2*cth/(3*tt**3), 2*sth/(3*tt**3), -6*cth/(3*tt**4), -6*sth/(3*tt**4)])
print(f"  GR time-shift mode: (bp, bm, da) ~ x^{expo(v_shift)[0]}")
# a 'pure GR' perturbation: change of Kasner angle (delta theta): beta+ -> -sin th ln t/3, beta- -> cos th ln t /3
lt = np.log(tt)
v_ang = np.array([0, 0, -sth*lt/3, cth*lt/3, -sth/(3*tt), cth/(3*tt), sth/(3*tt**2), -cth/(3*tt**2), -2*sth/(3*tt**3), 2*cth/(3*tt**3)])
print(f"  GR Kasner-angle mode: (bp, bm, da) ~ x^{expo(v_ang)[0]}")
