#!/usr/bin/env python3
"""
hdym_rung1_lyap.py -- numerical stage of Rung 1.

Loads rung1_system.pkl (from hdym_rung1.py): the 24-dim linear system xdot = K(q, p) x for the
six off-diagonal components u^a_i (a != i) and their first three time derivatives, on a
ghost-free diagonal Yang-Mills background q(t).  The system splits exactly into three 8-dim
blocks (u12,u21), (u13,u31), (u23,u32), each carrying one gauge zero mode, one conserved
Gauss charge, one YM-type degree of freedom and TWO ghost degrees of freedom (the second is
the longitudinal polarisation of the massive Lee-Wick vector, which has no diagonal analogue).

1. Frozen-background scan (p = 0 turning points on the shell): eigenvalues of K.  Real
   eigenvalues = instantaneous exponential directions.  Reports the threshold energy at which
   non-YM (ghost-type) real eigenvalues first appear.
2. Lyapunov spectrum of each block along chaotic YM trajectories (QR method), at several E.
   Expect per block {+lam_YM', +g1, +g2, 0, 0, -g2, -g1, -lam_YM'} with lam_YM' the rate at
   which the background leaves the diagonal ansatz, and g1, g2 the off-diagonal ghost rates.

Usage: python hdym_rung1_lyap.py --E 1.0 0.3 --N 16 --T 3000
"""
import argparse
import pickle
import time

import numpy as np
import sympy as sp

YC = (0.6756035959798288, -0.17560359597982886, -0.17560359597982886, 0.6756035959798288)
YD = (1.3512071919596576, -1.7024143839193153, 1.3512071919596576, 0.0)
BLOCKS = {'12': [0, 2], '13': [1, 4], '23': [3, 5]}      # indices into the u order (u12,u13,u21,u23,u31,u32)


def load_K():
    d = pickle.load(open('rung1_system.pkl', 'rb'))
    X, P = d['X'], d['P']
    Kf = sp.lambdify(list(X) + list(P), d['K'], 'numpy')
    Gf = sp.lambdify(list(X) + list(P), d['G'], 'numpy')
    return Kf, Gf


def V(q):
    x, y, z = q
    return 0.5 * (x * x * y * y + y * y * z * z + z * z * x * x)


def gradV(q):
    x, y, z = q
    return np.array([x * (y * y + z * z), y * (z * z + x * x), z * (x * x + y * y)])


def ym_step(q, p, h):
    for c, dcoef in zip(YC, YD):
        q = q + c * h * p
        if dcoef != 0.0:
            p = p - dcoef * h * gradV(q)
    return q, p


def block_idx(name):
    j = BLOCKS[name]
    return [k * 6 + jj for k in range(4) for jj in j]


def frozen_scan(Kf, E, nsamp, rng):
    """p = 0 points on the shell V = E: real parts of eig(K) per block, max over samples."""
    out = {b: 0.0 for b in BLOCKS}
    out_arg = {b: None for b in BLOCKS}
    for _ in range(nsamp):
        q = rng.standard_normal(3)
        q *= (E / V(q)) ** 0.25
        K = np.array(Kf(*q, 0.0, 0.0, 0.0), dtype=float)
        for b in BLOCKS:
            idx = block_idx(b)
            ev = np.linalg.eigvals(K[np.ix_(idx, idx)])
            r = ev.real.max()
            if r > out[b]:
                out[b] = r
                out_arg[b] = (q.copy(), np.sort_complex(ev))
    return out, out_arg


def lyap_block(Kf, Gf, E, T, dt, rng, name, nqr=10, tburn=100.0):
    idx = block_idx(name)
    q = rng.standard_normal(3)
    p = rng.standard_normal(3)
    s = (E / (0.5 * p @ p + V(q))) ** 0.25
    q *= s
    p *= s * s
    for _ in range(int(tburn / dt)):
        q, p = ym_step(q, p, dt)
    n = len(idx)
    Q = np.linalg.qr(rng.standard_normal((n, n)))[0]
    acc = np.zeros(n)
    nsteps = int(T / dt)
    Gini = None
    for k in range(1, nsteps + 1):
        K0 = np.array(Kf(*q, *p), dtype=float)[np.ix_(idx, idx)]
        qh, ph = ym_step(q, p, dt / 2)
        Kh = np.array(Kf(*qh, *ph), dtype=float)[np.ix_(idx, idx)]
        q, p = ym_step(qh, ph, dt / 2)
        K1 = np.array(Kf(*q, *p), dtype=float)[np.ix_(idx, idx)]
        k1 = K0 @ Q
        k2 = Kh @ (Q + 0.5 * dt * k1)
        k3 = Kh @ (Q + 0.5 * dt * k2)
        k4 = K1 @ (Q + dt * k3)
        Q = Q + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if k % nqr == 0:
            Q, R = np.linalg.qr(Q)
            acc += np.log(np.abs(np.diag(R)))
    lam = acc / (nsteps * dt)
    return np.sort(lam)[::-1], 0.5 * p @ p + V(q)


def ym_mass(q, name):
    x, y, z = q
    if name == '12':
        return np.array([[z * z, -x * y], [-x * y, z * z]])
    if name == '13':
        return np.array([[y * y, -x * z], [-x * z, y * y]])
    return np.array([[x * x, -y * z], [-y * z, x * x]])


def lyap_ym_block(E, T, dt, rng, name, nqr=10, tburn=100.0):
    """Pure Yang-Mills: uddot = -m(t) u for the off-diagonal pair; 4-dim; the top exponent is the
    rate at which the background leaves the diagonal ansatz (includes the gauge zero mode and
    the conserved Gauss charge, which give two zero exponents)."""
    q = rng.standard_normal(3)
    p = rng.standard_normal(3)
    s = (E / (0.5 * p @ p + V(q))) ** 0.25
    q *= s
    p *= s * s
    for _ in range(int(tburn / dt)):
        q, p = ym_step(q, p, dt)
    Q = np.linalg.qr(rng.standard_normal((4, 4)))[0]
    acc = np.zeros(4)
    nsteps = int(T / dt)

    def Kmat(q):
        K = np.zeros((4, 4))
        K[0, 2] = K[1, 3] = 1.0
        K[2:, :2] = -ym_mass(q, name)
        return K
    for k in range(1, nsteps + 1):
        K0 = Kmat(q)
        qh, ph = ym_step(q, p, dt / 2)
        Kh = Kmat(qh)
        q, p = ym_step(qh, ph, dt / 2)
        K1 = Kmat(q)
        k1 = K0 @ Q
        k2 = Kh @ (Q + 0.5 * dt * k1)
        k3 = Kh @ (Q + 0.5 * dt * k2)
        k4 = K1 @ (Q + dt * k3)
        Q = Q + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if k % nqr == 0:
            Q, R = np.linalg.qr(Q)
            acc += np.log(np.abs(np.diag(R)))
    return np.sort(acc / (nsteps * dt))[::-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--E', type=float, nargs='+', default=[1.0, 0.3])
    ap.add_argument('--N', type=int, default=8)
    ap.add_argument('--T', type=float, default=2000.0, help='run length in YM units (scaled by E^-1/4)')
    ap.add_argument('--dt', type=float, default=0.005)
    ap.add_argument('--scan', type=int, default=4000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--no-scan', action='store_true')
    ap.add_argument('--ym-only', action='store_true')
    a = ap.parse_args()
    Kf, Gf = load_K()
    rng = np.random.default_rng(a.seed)
    if not a.no_scan:
      print("=== frozen (p = 0) scan: largest real eigenvalue of K per block, max over the shell ===")
      print("    diagonal-sector reference: ghost tachyon iff E >= 1/2; YM-type instability of the ansatz: sqrt(|xy| - z^2) etc.")
      for E in (0.1, 0.3, 0.45, 0.49, 0.51, 0.6, 1.0):
          out, arg = frozen_scan(Kf, E, a.scan, rng)
          q, ev = arg['12']
          re = np.sort(ev.real)[::-1]
          print(f"  E={E:5.2f}: max Re per block {out['12']:.4f} {out['13']:.4f} {out['23']:.4f} | block 12 at q={np.round(q, 3)}: "
                f"Re(eig) = {np.round(re, 4)}")
    print("=== Lyapunov spectra along chaotic trajectories (per block, sorted) ===")
    for E in a.E:
        t0 = time.time()
        ym = {b: [] for b in BLOCKS}
        for i in range(a.N):
            for b in BLOCKS:
                ym[b].append(lyap_ym_block(E, a.T * E ** -0.25, a.dt / max(1.0, E ** 0.25), rng, b))
        for b in BLOCKS:
            S = np.array(ym[b])
            print(f"  E={E:5.3g} pure-YM block {b}: " + "  ".join(f"{v:+.4f}({e:.4f})" for v, e in zip(S.mean(0), S.std(0, ddof=1) / np.sqrt(a.N))))
        np.save(f'rung1_ym_E{E:g}.npy', np.concatenate([np.array(ym[b]) for b in BLOCKS], 0))
        if a.ym_only:
            continue
        spec = {b: [] for b in BLOCKS}
        for i in range(a.N):
            for b in BLOCKS:
                lam, E1 = lyap_block(Kf, Gf, E, a.T * E ** -0.25, a.dt / max(1.0, E ** 0.25), rng, b)
                spec[b].append(lam)
        for b in BLOCKS:
            S = np.array(spec[b])
            m, se = S.mean(0), S.std(0, ddof=1) / np.sqrt(a.N)
            print(f"  E={E:5.3g} block {b}: " + "  ".join(f"{v:+.4f}({e:.4f})" for v, e in zip(m, se)))
        allS = np.concatenate([np.array(spec[b]) for b in BLOCKS], 0)
        print(f"  E={E:5.3g} all blocks: lam_YM(E) ref = {0.38 * E ** 0.25:.4f}; top exponent {allS[:, 0].mean():.4f}, "
              f"second {allS[:, 1].mean():.4f}, third {allS[:, 2].mean():.4f}   ({time.time() - t0:.0f}s)")
        np.save(f'rung1_lyap_E{E:g}.npy', allS)


if __name__ == '__main__':
    main()
