"""dr_offdiag_floquet.py -- the off-diagonal ghost on the Dahlqvist-Russberg orbit (float).

Everything else in islands/ is about the colour-diagonal ghost.  The full homogeneous theory
also has six off-diagonal components u^a_i (a != i); ym/hdym_rung1.py linearises the full
Lee-Wick system about the diagonal background and exports the 24-dim system xdot = K(q, p) x,
which splits exactly into three 8-dim blocks (u12,u21), (u13,u31), (u23,u32).  Each block holds
a gauge zero mode, a conserved Gauss charge, one Yang-Mills-type degree of freedom (leaving the
diagonal ansatz) and two ghost degrees of freedom.

Here: the monodromy of each block over one period of the DR orbit (planar, z = pz = 0), for a
range of energies (the orbit is the E = 1/2 one scaled: q -> lam q, p -> lam^2 p, t -> t/lam,
lam = (2E)^{1/4}; the ghost mass breaks the scaling, the background does not).  Reported per block:
the largest |multiplier| and the Floquet exponent ln|mu|_max / T.  Separately, the pure
Yang-Mills off-diagonal block (uddot = -m(q) u, 4-dim) on the same orbit, so the YM-type
multipliers can be told apart from the ghost ones.

Checks printed: closure of the background orbit; det M = 1; multipliers come in (mu, 1/mu)
pairs; the multiplier 1 (zero mode, Gauss charge) present in every block; with c_j = 0 the block
reduces to the scalar Hill equations wddot = -(1 + z^2 -+ xy) w (paper, off-diagonal section),
whose monodromies are computed independently and must match.

Float only (DOP853, rtol 1e-12); nothing here is certified.
Usage (from islands/):  python code/dr_offdiag_floquet.py > results/dr_offdiag_floquet.txt
Added 2026-10-12.
"""
import os
import pickle
import sys
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
PKL = os.path.join(HERE, '..', '..', 'ym', 'rung1_system.pkl')

# DR orbit at E = 1/2 (islands/README.md, Conventions; certified box in results/cert16_check_output.txt)
X0, PX0, T_HALF = 3.14640122769837, 0.00179319341881, 31.053330268585066
BLOCKS = {'12': [0, 2], '13': [1, 4], '23': [3, 5]}      # u order (u12,u13,u21,u23,u31,u32)
RT = dict(method='DOP853', rtol=1e-12, atol=1e-13)


def block_idx(name):
    return [k * 6 + jj for k in range(4) for jj in BLOCKS[name]]


def load_K():
    d = pickle.load(open(PKL, 'rb'))
    X, P, cj = d['X'], d['P'], d['cj']
    return sp.lambdify(list(X) + list(P) + [cj], d['K'], 'numpy')


def gradV(q):
    x, y, z = q
    return np.array([x * (y * y + z * z), y * (z * z + x * x), z * (x * x + y * y)])


def ym_mass(q, name):
    x, y, z = q
    if name == '12':
        return np.array([[z * z, -x * y], [-x * y, z * z]])
    if name == '13':
        return np.array([[y * y, -x * z], [-x * z, y * y]])
    return np.array([[x * x, -y * z], [-y * z, x * x]])


def orbit(E):
    lam = (2.0 * E) ** 0.25
    q0 = np.array([lam * X0, 0.0, 0.0])
    p0 = lam * lam * np.array([PX0, np.sqrt(1.0 - PX0 ** 2), 0.0])
    return q0, p0, T_HALF / lam


def monodromies(Kf, E, cj):
    q0, p0, T = orbit(E)
    nb = len(BLOCKS)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, cj), dtype=float)
        out = [p, -gradV(q)]
        off = 6
        for b in BLOCKS:
            idx = block_idx(b)
            M = s[off:off + 64].reshape(8, 8)
            out.append((K[np.ix_(idx, idx)] @ M).ravel())
            off += 64
        for b in BLOCKS:
            M = s[off:off + 16].reshape(4, 4)
            A = np.zeros((4, 4))
            A[0, 2] = A[1, 3] = 1.0
            A[2:, :2] = -ym_mass(q, b)
            out.append((A @ M).ravel())
            off += 16
        return np.concatenate(out)

    s0 = np.concatenate([q0, p0] + [np.eye(8).ravel()] * nb + [np.eye(4).ravel()] * nb)
    sol = solve_ivp(rhs, (0.0, T), s0, **RT)
    sT = sol.y[:, -1]
    closure = np.max(np.abs(sT[:6] - s0[:6]))
    Ms, Ys = {}, {}
    off = 6
    for b in BLOCKS:
        Ms[b] = sT[off:off + 64].reshape(8, 8)
        off += 64
    for b in BLOCKS:
        Ys[b] = sT[off:off + 16].reshape(4, 4)
        off += 16
    return T, closure, Ms, Ys


def hill_monodromy(E, sign):
    """c_j = 0 check: wddot = -(1 + z^2 - sign*x*y) w on the planar orbit, scalar."""
    q0, p0, T = orbit(E)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        x, y, z = q
        w, wd, w2, wd2 = s[6:]
        om2 = 1.0 + z * z - sign * x * y
        return np.concatenate([p, -gradV(q), [wd, -om2 * w, wd2, -om2 * w2]])
    s0 = np.concatenate([q0, p0, [1.0, 0.0, 0.0, 1.0]])
    sT = solve_ivp(rhs, (0.0, T), s0, **RT).y[:, -1]
    return np.array([[sT[6], sT[8]], [sT[7], sT[9]]])


def summary(M):
    mu = np.linalg.eigvals(M)
    a = np.sort(np.abs(mu))[::-1]
    return mu, a


def ghost_multipliers(M):
    """The 8 multipliers minus the 4 nearest 1 (gauge zero mode, Gauss charge and their Jordan
    partners: a defective eigenvalue 1, which floating point splits by ~1e-5).  What is left
    belongs to the two ghost degrees of freedom; sorted by modulus."""
    mu = np.linalg.eigvals(M)
    keep = np.argsort(np.abs(mu - 1.0))[4:]
    g = mu[keep]
    return g[np.argsort(-np.abs(g))]


def main():
    t0 = time.time()
    print(f"dr_offdiag_floquet.py  python {sys.version.split()[0]}  numpy {np.__version__}  "
          f"sympy {sp.__version__}  ({sys.executable})")
    Kf = load_K()

    # ---- c_j = 0 check against the scalar Hill equations, at E = 1/2 and 0.3
    print("\nCheck (c_j = 0): block multipliers vs the scalar Hill equations w'' = -(1 + z^2 -+ xy) w")
    for E in (0.5, 0.3):
        T, clo, Ms, _ = monodromies(Kf, E, 0.0)
        hill = sorted([abs(m) for s in (+1, -1) for m in np.linalg.eigvals(hill_monodromy(E, s))], reverse=True)
        _, a12 = summary(Ms['12'])
        print(f"  E={E}: block 12 |mu| = {np.round(a12, 6)}")
        print(f"          Hill  |mu| = {np.round(hill, 6)}   (must appear among the block's)")

    # ---- the theory, c_j = 1
    print("\nc_j = 1 (the theory).  Per block, the ghost multipliers (the 4 not near 1): max |mu|, and the")
    print("Floquet rate ln|mu|/T when |mu| > 1.  'YM' = largest |mu| of the pure-YM off-diagonal block.")
    print(f"{'E':>6} {'T':>8} {'closure':>8} | {'blk12':>8} {'blk13':>8} {'blk23':>8} | {'YM max':>7} | "
          f"{'det-1':>8} {'pair':>8}  verdict")

    def row(E):
        T, clo, Ms, Ys = monodromies(Kf, E, 1.0)
        mx, rate, dets, pairs = [], [], [], []
        for b in BLOCKS:
            g = ghost_multipliers(Ms[b])
            mx.append(abs(g[0]))
            rate.append(np.log(abs(g[0])) / T)
            dets.append(abs(np.linalg.det(Ms[b]) - 1.0))
            pairs.append(abs(np.prod(np.abs(g)) - 1.0))       # (mu, 1/mu) pairing: product of moduli 1
        ym = max(summary(Ys[b])[1][0] for b in BLOCKS)
        unst = [m > 1.0 + 1e-6 for m in mx]
        verdict = ("UNSTABLE " + ", ".join(f"blk{b} rate {r:.4f}" for b, r, u in zip(BLOCKS, rate, unst) if u)
                   if any(unst) else "stable")
        print(f"{E:6.4f} {T:8.4f} {clo:8.1e} | " + " ".join(f"{m:8.4f}" for m in mx)
              + f" | {ym:7.5f} | {max(dets):8.1e} {max(pairs):8.1e}  {verdict}", flush=True)
        return any(unst)

    Es = [1.0, 0.75, 0.6, 0.5, 0.45, 0.4, 0.35, 0.315, 0.3, 0.28, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05]
    for E in Es:
        row(E)

    print("\nFine scan, E = 0.05 ... 1 (400 points, log-spaced): unstable intervals of the off-diagonal ghost")
    grid = np.logspace(np.log10(0.05), 0.0, 400)
    flags = []
    for E in grid:
        T, clo, Ms, Ys = monodromies(Kf, E, 1.0)
        flags.append(max(abs(ghost_multipliers(Ms[b])[0]) for b in BLOCKS) > 1.0 + 1e-6)
    flags = np.array(flags)
    k = 0
    while k < len(grid):
        if flags[k]:
            j = k
            while j + 1 < len(grid) and flags[j + 1]:
                j += 1
            print(f"  unstable  E in [{grid[k]:.4f}, {grid[j]:.4f}]  ({j - k + 1} points)")
            k = j + 1
        else:
            k += 1
    print(f"  unstable fraction of the grid: {flags.mean():.2f}")
    i3 = np.argmin(np.abs(grid - 0.3))
    print(f"  grid point nearest E = 0.3: {grid[i3]:.4f} -> {'UNSTABLE' if flags[i3] else 'stable'}")

    print("\nTolerance check at E = 0.3 (block 23 ghost multipliers):")
    global RT
    for rt in (1e-9, 1e-11, 1e-12):
        RT = dict(method='DOP853', rtol=rt, atol=rt * 0.1)
        T, clo, Ms, Ys = monodromies(Kf, 0.3, 1.0)
        g = ghost_multipliers(Ms['23'])
        print(f"  rtol {rt:.0e}: " + "  ".join(f"{m.real:+.6f}{m.imag:+.6f}i" for m in g))

    print("\nFull multipliers, c_j = 1:")
    for E in (0.5, 0.3):
        T, clo, Ms, Ys = monodromies(Kf, E, 1.0)
        for b in BLOCKS:
            mu = np.linalg.eigvals(Ms[b])
            mu = mu[np.argsort(-np.abs(mu))]
            print(f"  E={E} block {b}: " + "  ".join(f"{m.real:+.5f}{m.imag:+.5f}i" for m in mu))
    print(f"\n({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
