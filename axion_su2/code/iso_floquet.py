"""iso_floquet.py -- every linear sector on the isotropic Yang-Mills orbit x = y = z (float).

Background.  The homogeneous isotropic SU(2) field A^a_i = delta^a_i s(t) is the background of
axion-SU(2) inflation; after inflation, once the axion decouples, it oscillates under the quartic
self-coupling (Fujita, Murai, Nakayama, Yin, arXiv:2312.06889: Qddot + 2 g^2 Q^3 = 0 in conformal
time; flat space here, no expansion, no axion).  In this repository's units (g = M = 1,
V = (x^2 y^2 + y^2 z^2 + z^2 x^2)/2) it is the orbit x = y = z = s(t) of the 12-d Lee-Wick system:
    s'' = -2 s^3,   E = (3/2)(s'^2 + s^4),   s = A cn(sqrt2 A t | k^2 = 1/2),   A = (2E/3)^{1/4},
    period T = 4 K(1/2) / (sqrt2 A).

Sectors, each by its monodromy over one period:
  D-YM   diagonal Yang-Mills perturbations  dq'' = -Hess V dq.  Hess V = 2 s^2 (all-ones matrix):
         the isotropic direction is the orbit's own tangent; the two anisotropic directions feel no
         force (marginal, linear drift).  Reported for completeness.
  D-gh   diagonal ghost  w'' = -(1 + Hess V) w.  Isotropic mode: w'' = -(1 + 6 s^2) w, which in
         tau = sqrt2 A t is the Lame equation w_tt + (h - 6 k^2 sn^2) w = 0 with n = 2, k^2 = 1/2,
         h = 3 + 1/(2A^2).  Band edges h = 1.268, 1.5, 3, 4.5, 4.732, so for h >= 3 the only gap is
         (4.5, 4.732), i.e. E in (1/8, 1/6) exactly.  Anisotropic modes: w'' = -w, free.
  O-YM   pure Yang-Mills off-diagonal block (u_ab, u_ba), 4-dim, mass m(q) as in hdym_rung1_lyap.py
  O-gh   the off-diagonal 8-dim blocks of ym/rung1_system.pkl (full Lee-Wick, c_j = 1, the theory;
         c_j = 0 as a diagnostic).  By the isotropy all three blocks must give the same multipliers
         (checked).  "Ghost multipliers" = the 8 minus the 4 nearest 1 (gauge zero mode, Gauss
         charge, Jordan partners), as in islands/code/dr_offdiag_floquet.py.

Checks printed: orbit closure; the period against 4K(1/2)/(sqrt2 A); D-gh unstable interval
against (1/8, 1/6); the three O-gh blocks against each other; det M = 1.

Float only (DOP853, rtol 1e-11); nothing here is certified.
Usage (from the repository root):
    python axion_su2/code/iso_floquet.py > axion_su2/results/iso_floquet.txt
Added 2026-10-12.
"""
import os
import pickle
import sys
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.special import ellipk

HERE = os.path.dirname(os.path.abspath(__file__))
PKL = os.path.join(HERE, '..', '..', 'ym', 'rung1_system.pkl')
BLOCKS = {'12': [0, 2], '13': [1, 4], '23': [3, 5]}      # u order (u12,u13,u21,u23,u31,u32)
RT = dict(method='DOP853', rtol=1e-11, atol=1e-13)


def block_idx(name):
    return [k * 6 + jj for k in range(4) for jj in BLOCKS[name]]


def load_K():
    d = pickle.load(open(PKL, 'rb'))
    return sp.lambdify(list(d['X']) + list(d['P']) + [d['cj']], d['K'], 'numpy')


def gradV(q):
    x, y, z = q
    return np.array([x * (y * y + z * z), y * (z * z + x * x), z * (x * x + y * y)])


def hessV(q):
    x, y, z = q
    return np.array([[y * y + z * z, 2 * x * y, 2 * x * z],
                     [2 * x * y, z * z + x * x, 2 * y * z],
                     [2 * x * z, 2 * y * z, x * x + y * y]])


def ym_mass(q, name):
    x, y, z = q
    if name == '12':
        return np.array([[z * z, -x * y], [-x * y, z * z]])
    if name == '13':
        return np.array([[y * y, -x * z], [-x * z, y * y]])
    return np.array([[x * x, -y * z], [-y * z, x * x]])


def orbit(E):
    A = (2.0 * E / 3.0) ** 0.25
    T = 4.0 * ellipk(0.5) / (np.sqrt(2.0) * A)
    return np.array([A, A, A]), np.zeros(3), T


def second_order(M2):
    """4x4 (or 6x6) first-order matrix for w'' = -M2 w."""
    n = M2.shape[0]
    A = np.zeros((2 * n, 2 * n))
    A[:n, n:] = np.eye(n)
    A[n:, :n] = -M2
    return A


def monodromies(Kf, E, cj):
    q0, p0, T = orbit(E)
    sizes = [('DYM', 6), ('Dgh', 6)] + [('OYM' + b, 4) for b in BLOCKS] + [('Ogh' + b, 8) for b in BLOCKS]

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        H = hessV(q)
        K = np.array(Kf(*q, *p, cj), dtype=float)
        out = [p, -gradV(q)]
        off = 6
        for name, n in sizes:
            M = s[off:off + n * n].reshape(n, n)
            if name == 'DYM':
                A = second_order(H)
            elif name == 'Dgh':
                A = second_order(np.eye(3) + H)
            elif name.startswith('OYM'):
                A = second_order(ym_mass(q, name[3:]))
            else:
                idx = block_idx(name[3:])
                A = K[np.ix_(idx, idx)]
            out.append((A @ M).ravel())
            off += n * n
        return np.concatenate(out)

    s0 = np.concatenate([q0, p0] + [np.eye(n).ravel() for _, n in sizes])
    sT = solve_ivp(rhs, (0.0, T), s0, **RT).y[:, -1]
    closure = np.max(np.abs(sT[:6] - s0[:6]))
    Ms, off = {}, 6
    for name, n in sizes:
        Ms[name] = sT[off:off + n * n].reshape(n, n)
        off += n * n
    return T, closure, Ms


def ghost_multipliers(M):
    mu = np.linalg.eigvals(M)
    g = mu[np.argsort(np.abs(mu - 1.0))[4:]]
    return g[np.argsort(-np.abs(g))]


def maxmod(M):
    return np.max(np.abs(np.linalg.eigvals(M)))


def intervals(grid, flags):
    out, k = [], 0
    while k < len(grid):
        if flags[k]:
            j = k
            while j + 1 < len(grid) and flags[j + 1]:
                j += 1
            out.append((grid[k], grid[j], j - k + 1))
            k = j + 1
        else:
            k += 1
    return out


def main():
    t0 = time.time()
    print(f"iso_floquet.py  python {sys.version.split()[0]}  numpy {np.__version__}  sympy {sp.__version__}"
          f"  ({sys.executable})")
    Kf = load_K()

    print("\nSpot values (c_j = 1).  Columns: max |mu| per sector; O-gh = largest ghost multiplier of each block.")
    print(f"{'E':>7} {'T':>8} {'closure':>8} | {'D-YM':>8} {'D-gh':>8} | {'O-YM':>8} | "
          f"{'O-gh12':>8} {'O-gh13':>8} {'O-gh23':>8} | {'det-1':>7}")
    for E in (0.01, 0.05, 0.1, 0.14, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0):
        T, clo, Ms = monodromies(Kf, E, 1.0)
        og = [abs(ghost_multipliers(Ms['Ogh' + b])[0]) for b in BLOCKS]
        oy = max(maxmod(Ms['OYM' + b]) for b in BLOCKS)
        dets = max(abs(np.linalg.det(Ms[k]) - 1) for k in Ms)
        print(f"{E:7.3f} {T:8.4f} {clo:8.1e} | {maxmod(Ms['DYM']):8.4f} {maxmod(Ms['Dgh']):8.4f} | {oy:8.4f} | "
              + " ".join(f"{v:8.4f}" for v in og) + f" | {dets:7.1e}", flush=True)

    for cj in (1.0, 0.0):
        print(f"\nScan E = 0.005 ... 10 (500 points, log-spaced), c_j = {cj:g}: unstable intervals per sector")
        grid = np.logspace(np.log10(0.005), 1.0, 500)
        F = {k: [] for k in ('Dgh', 'OYM', 'Ogh')}
        rate = {k: [] for k in F}
        spread = 0.0
        for E in grid:
            T, clo, Ms = monodromies(Kf, E, cj)
            dg = maxmod(Ms['Dgh'])
            oy = max(maxmod(Ms['OYM' + b]) for b in BLOCKS)
            og3 = [abs(ghost_multipliers(Ms['Ogh' + b])[0]) for b in BLOCKS]
            spread = max(spread, max(og3) - min(og3))
            for k, v in (('Dgh', dg), ('OYM', oy), ('Ogh', max(og3))):
                F[k].append(v > 1.0 + (1e-3 if k == 'OYM' else 1e-5))   # O-YM is parabolic (Jordan): float splits it by ~1e-5
                rate[k].append(np.log(v) / T)
        for k, label in (('Dgh', 'diagonal ghost'), ('OYM', 'off-diagonal Yang-Mills'), ('Ogh', 'off-diagonal ghost')):
            iv = intervals(grid, np.array(F[k]))
            print(f"  {label}: {len(iv)} interval(s); unstable fraction of grid {np.mean(F[k]):.2f}")
            for a, b, n in iv:
                sel = (grid >= a) & (grid <= b)
                print(f"     E in [{a:.4f}, {b:.4f}]  ({n} pts)  max rate {np.max(np.array(rate[k])[sel]):.4g}")
        print(f"  max spread of the largest ghost multiplier across the three blocks: {spread:.1e}  (isotropy check)")
        if cj == 1.0:
            print("  diagonal ghost, closed form (Lame n = 2, k^2 = 1/2): unstable for E in (1/8, 1/6) = (0.1250, 0.1667)")
    print(f"\n({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
