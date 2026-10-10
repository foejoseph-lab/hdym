"""hdym2d_nu0.py -- channel statistics of the regulated two-field model at E = 1.

    V = 1/2 x^2 y^2 + 1/2 mu1^2 (x^2 + y^2)

Pure (ghost-free) dynamics, Yoshida-4, N trajectories.  A visit to the x-channel starts
when |x| crosses A_in upward with |y| < A_in and ends when it drops back below; same for y.
Records D_max and the transverse action J = E_perp / sqrt(x^2 + mu^2) at entry.

Predictions (microcanonical flux, uniform action measure, Sec. 2 of NOTE_2d.md):
    nu0 = 4 pi / Z2(mu1),   Z2 = 2 pi * area{ V < 1 }
    rate(any channel, D_max > D) = 2 nu0 * J_max(D),  J_max(D) = (1 - mu^2 D^2 / 2) / D
    rate(J_entry < j)            = 2 nu0 * j
    D_max * J_entry + mu^2 D_max^2 / 2 = 1

Usage:  python hdym2d_nu0.py --mu1 0.1 --N 256 --T 3000
"""
import argparse
import numpy as np
from scipy.integrate import quad

YC = (0.6756035959798288, -0.17560359597982886, -0.17560359597982886, 0.6756035959798288)
YD = (1.3512071919596576, -1.7024143839193153, 1.3512071919596576, 0.0)


def Z2(mu):
    f = lambda x: np.sqrt((2 - mu * mu * x * x) / (x * x + mu * mu))
    v, _ = quad(f, 0, np.sqrt(2) / mu, limit=500, points=[1e-3, 1e-2, 0.1, 1, 10])
    return 2 * np.pi * 4 * v


def V(q, mu2):
    x, y = q
    return 0.5 * x * x * y * y + 0.5 * mu2 * (x * x + y * y)


def grad(q, mu2):
    x, y = q
    return np.array([x * y * y + mu2 * x, y * x * x + mu2 * y])


def step(q, p, h, mu2):
    for c, d in zip(YC, YD):
        q = q + c * h * p
        if d != 0.0:
            p = p - d * h * grad(q, mu2)
    return q, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mu1', type=float, default=0.1)
    ap.add_argument('--N', type=int, default=256)
    ap.add_argument('--T', type=float, default=3000.0)
    ap.add_argument('--tburn', type=float, default=100.0)
    ap.add_argument('--dt', type=float, default=0.01)
    ap.add_argument('--A-in', type=float, default=3.0)
    ap.add_argument('--seed', type=int, default=5)
    a = ap.parse_args()
    mu2 = a.mu1 ** 2
    rng = np.random.default_rng(a.seed)
    q = rng.standard_normal((2, a.N))
    p = rng.standard_normal((2, a.N))
    E0 = 0.5 * np.sum(p * p, 0) + V(q, mu2)
    # rescale to E = 1 by bisection on the overall amplitude (not scale-free with mu fixed)
    s = np.ones(a.N)
    for _ in range(60):
        E = 0.5 * s ** 4 * np.sum(p * p, 0) + V(s * q, mu2)
        s *= (1.0 / E) ** 0.25 * 0.5 + 0.5
    q *= s
    p *= s ** 2
    E = 0.5 * np.sum(p * p, 0) + V(q, mu2)
    print(f"mu1={a.mu1}  N={a.N}  E in [{E.min():.6f}, {E.max():.6f}]")
    nsteps, nburn = int(a.T / a.dt), int(a.tburn / a.dt)
    inch = -np.ones(a.N, int)
    dmax = np.zeros(a.N)
    Jent = np.zeros(a.N)
    D_list, J_list = [], []
    for k in range(nsteps):
        q, p = step(q, p, a.dt, mu2)
        if k < nburn:
            continue
        aq = np.abs(q)
        i = np.argmax(aq, 0)
        d = aq[i, np.arange(a.N)]
        ent = (inch < 0) & (d > a.A_in)
        if ent.any():
            idx = np.nonzero(ent)[0]
            ii = i[idx]
            x = q[ii, idx]
            yperp = q[1 - ii, idx]
            pperp = p[1 - ii, idx]
            Om = np.sqrt(x * x + mu2)
            Eperp = 0.5 * pperp ** 2 + 0.5 * Om ** 2 * yperp ** 2
            Jent[idx] = Eperp / Om
            inch[idx] = ii
            dmax[idx] = d[idx]
        inv = inch >= 0
        dmax[inv] = np.maximum(dmax[inv], d[inv])
        ex = inv & (d < a.A_in)
        if ex.any():
            D_list.extend(dmax[ex])
            J_list.extend(Jent[ex])
            inch[ex] = -1
    Ttot = (nsteps - nburn) * a.dt * a.N
    D = np.array(D_list)
    J = np.array(J_list)
    E1 = 0.5 * np.sum(p * p, 0) + V(q, mu2)
    print(f"max |E1-E0|/E0 = {np.max(np.abs(E1 - E)):.1e}")
    Z = Z2(a.mu1)
    nu0 = 4 * np.pi / Z
    print(f"Z2 = {Z:.2f}, predicted nu0 = 4pi/Z2 = {nu0:.4f}")
    print(f"{len(D)} visits in {Ttot:.0f} trajectory-time units; visit rate {len(D) / Ttot:.4f}")
    print("rate(any channel, Dmax > D) / (2 Jmax(D))  -> nu0:")
    Dcut = np.sqrt(2) / a.mu1
    for Dc in (4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40):
        if Dc >= 0.9 * Dcut:
            break
        r = np.sum(D > Dc) / Ttot
        Jmax = (1 - 0.5 * mu2 * Dc * Dc) / Dc
        print(f"  D={Dc:3d}  rate={r:.3e}  nu0_meas={r / (2 * Jmax):.4f}  ratio={r / (2 * Jmax) / nu0:.3f}  (n={np.sum(D > Dc)})")
    print("rate(J_entry < j) / (2 j) -> nu0:")
    for jc in (0.02, 0.05, 0.1, 0.2, 0.3):
        r = np.sum(J < jc) / Ttot
        print(f"  j={jc:.2f}  nu0_meas={r / (2 * jc):.4f}  ratio={r / (2 * jc) / nu0:.3f}  (n={np.sum(J < jc)})")
    good = (J > 0) & (D > 4)
    adiab = D[good] * J[good] + 0.5 * mu2 * D[good] ** 2
    print(f"adiabatic check D*J + mu^2 D^2/2 (should be 1): median {np.median(adiab):.3f}, "
          f"IQR {np.percentile(adiab, 25):.3f}-{np.percentile(adiab, 75):.3f}")
    np.savez(f"nu0_2d_mu{a.mu1}.npz", D=D, J=J, Ttot=Ttot, mu1=a.mu1, Z2=Z)


if __name__ == '__main__':
    main()
