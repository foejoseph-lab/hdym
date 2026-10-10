import numpy as np
from scipy.integrate import solve_ivp

# Dahlqvist-Russberg stable period-9 orbit in V = x^2 y^2 / 2 at E = 1/2
# (Contopoulos & Harsoula 2023, arXiv:2302.12071)
x0, vx0, y0, vy0 = 3.14640122769753, 0.0017931934191, 0.0, 0.99999839222738
E0 = 0.5

def rhs_base(t, s):
    x, y, vx, vy = s
    return [vx, vy, -x*y*y, -y*x*x]

# --- find period: 9 upward crossings of y=0 ---
def ev(t, s): return s[1]
ev.direction = 1
sol = solve_ivp(rhs_base, (0, 200), [x0, y0, vx0, vy0], rtol=1e-12, atol=1e-14,
                events=ev, dense_output=True)
tc = sol.t_events[0]
# first event at t=0 may or may not register; take crossings strictly > 1e-6
tc = tc[tc > 1e-6]
T = tc[8]
s_T = sol.sol(T)
print("period T =", T)
print("return error |s(T)-s(0)| =", np.linalg.norm(s_T - [x0, y0, vx0, vy0]))

# --- monodromy of base flow in-plane and z-direction, scale free, at E0 ---
def rhs_full(t, s):
    x, y, vx, vy = s[:4]
    out = [vx, vy, -x*y*y, -y*x*x]
    # in-plane variation (4x4 fundamental matrix, columns stacked)
    H = np.array([[y*y, 2*x*y], [2*x*y, x*x]])
    A = np.zeros((4, 4)); A[0, 2] = A[1, 3] = 1; A[2:, :2] = -H
    Phi = s[4:20].reshape(4, 4)
    out += list((A @ Phi).ravel())
    # z-direction: zdd = -(x^2+y^2) z, 2x2 fundamental
    Az = np.array([[0, 1], [-(x*x + y*y), 0]])
    Pz = s[20:24].reshape(2, 2)
    out += list((Az @ Pz).ravel())
    return out

s_init = [x0, y0, vx0, vy0] + list(np.eye(4).ravel()) + list(np.eye(2).ravel())
r = solve_ivp(rhs_full, (0, T), s_init, rtol=1e-12, atol=1e-14)
M = r.y[4:20, -1].reshape(4, 4)
Mz = r.y[20:24, -1].reshape(2, 2)
eig = np.linalg.eigvals(M)
print("in-plane monodromy eigenvalues:", np.round(eig, 6))
# Henon index: non-trivial pair mu + 1/mu = 2 HI
nt = [e for e in eig if abs(e - 1) > 1e-3]
if nt:
    print("Henon index HI =", (nt[0] + 1/nt[0]).real / 2)
print("z-direction trace Tr Mz =", np.trace(Mz), " (|Tr|<2 stable)")
print("z-direction growth rate per unit time =", np.log(max(abs(np.linalg.eigvals(Mz)))) / T)

# --- ghost Floquet along the orbit at energy E ---
# scaling: q -> lam q, p -> lam^2 p, t -> t/lam, lam = (E/E0)^(1/4)
def ghost_floquet(E):
    lam = (E / E0) ** 0.25
    TE = T / lam
    def rhs(t, s):
        x, y, vx, vy = s[:4]
        out = [vx, vy, -x*y*y, -y*x*x]
        # in-plane ghost pair (wx, wy): wdd = -(1 + H) w
        H = np.array([[y*y, 2*x*y], [2*x*y, x*x]])
        A = np.zeros((4, 4)); A[0, 2] = A[1, 3] = 1; A[2:, :2] = -(np.eye(2) + H)
        P = s[4:20].reshape(4, 4)
        out += list((A @ P).ravel())
        # z ghost: wdd = -(1 + x^2 + y^2) w
        Az = np.array([[0, 1], [-(1 + x*x + y*y), 0]])
        Pz = s[20:24].reshape(2, 2)
        out += list((Az @ Pz).ravel())
        return out
    s_init = [lam*x0, lam*y0, lam**2*vx0, lam**2*vy0] + list(np.eye(4).ravel()) + list(np.eye(2).ravel())
    r = solve_ivp(rhs, (0, TE), s_init, rtol=1e-11, atol=1e-13)
    P = r.y[4:20, -1].reshape(4, 4)
    Pz = r.y[20:24, -1].reshape(2, 2)
    g_pair = np.log(max(abs(np.linalg.eigvals(P)))) / TE
    g_z = np.log(max(abs(np.linalg.eigvals(Pz)))) / TE
    return g_pair, g_z, TE

print("\n   E        T(E)     ghost pair rate   z-ghost rate   (rate < ~1e-6 means stable)")
Es = np.concatenate([np.logspace(-2, 0, 41), np.linspace(1.2, 4.0, 15)])
rows = []
for E in Es:
    gp, gz, TE = ghost_floquet(E)
    rows.append((E, TE, gp, gz))
    print(f"{E:8.4f}  {TE:8.3f}   {gp:12.3e}   {gz:12.3e}")
np.save("dr_ghost_floquet.npy", np.array(rows))
