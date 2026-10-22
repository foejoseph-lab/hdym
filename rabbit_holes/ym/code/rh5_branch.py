"""rh5_branch.py -- which branch?  A Floquet exponent from M is defined mod Omega = 2 pi / T; the reviewer's test to fix
the physical branch of each ghost pair: integrate the Floquet solution x(t) = Phi_t v over one period, Fourier-analyse
it, and read the dominant frequency.  The solution e^{i omega t} p(t) has lines at omega + k Omega weighted by p's
harmonics; the strongest line is the physical frequency.  Spectrum = sum over the 8 block components of |FFT|^2
(Euclidean weighting -- not canonical, but the dominant line is what is asked).  Both pairs (the eigenvectors of M with
Im mu > 0, sorted by |arg|), at nE energies.
Usage: python rabbit_holes/ym/code/rh5_branch.py [--nE 40] > rabbit_holes/ym/results/rh5_branch.txt
Added 2026-10-21 (late).
"""
import argparse
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402

RES = os.path.join(HERE, '..', 'results')
IDX = base.block_idx('12')


def spectrum(Kf, E, nsamp=4096):
    q0, p0, T = base.orbit(E)
    # monodromy and eigenvectors
    def rhsM(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, 1.0), dtype=float)[np.ix_(IDX, IDX)]
        return np.concatenate([p, -base.gradV(q), (K @ s[6:].reshape(8, 8)).ravel()])
    sT = solve_ivp(rhsM, (0.0, T), np.concatenate([q0, p0, np.eye(8).ravel()]), **base.RT).y[:, -1]
    M = sT[6:].reshape(8, 8)
    mu, V = np.linalg.eig(M)
    keep = np.argsort(np.abs(mu - 1.0))[4:]
    g, Vg = mu[keep], V[:, keep]
    sel = [k for k in range(4) if g[k].imag > 0] if (g.imag > 1e-9).sum() == 2 else list(np.argsort(-np.abs(g))[:2])
    sel.sort(key=lambda k: abs(np.angle(g[k])))
    # propagate the selected eigenvectors densely over one period
    nv = len(sel)
    Vs = Vg[:, sel]

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, 1.0), dtype=float)[np.ix_(IDX, IDX)]
        X = s[6:].reshape(2 * nv, 8)
        return np.concatenate([p, -base.gradV(q), (X @ K.T).ravel()])
    s0 = np.concatenate([q0, p0, np.concatenate([Vs.real.T, Vs.imag.T]).ravel()])
    te = np.linspace(0, T, nsamp, endpoint=False)
    sol = solve_ivp(rhs, (0.0, T), s0, t_eval=te, **base.RT)
    X = sol.y[6:].reshape(2 * nv, 8, nsamp)
    out = []
    freqs = np.fft.fftfreq(nsamp, d=T / nsamp) * 2 * np.pi
    for k in range(nv):
        x = X[k] + 1j * X[nv + k]                                  # (8, nsamp) complex solution
        # field-level components only (u, u'): weight by the Euclidean norm of (u12, u21, u12', u21')
        P = np.abs(np.fft.fft(x[:4], axis=1)) ** 2
        P = P.sum(axis=0)
        order = np.argsort(-P)
        f1, f2, f3 = freqs[order[0]], freqs[order[1]], freqs[order[2]]
        w1, w2 = P[order[0]] / P.sum(), P[order[1]] / P.sum()
        out.append((g[sel[k]], f1, w1, f2, w2, f3))
    return T, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nE', type=int, default=40)
    ap.add_argument('--emin', type=float, default=0.05)
    ap.add_argument('--emax', type=float, default=1.0)
    a = ap.parse_args()
    Kf = base.load_K()
    grid = np.logspace(np.log10(a.emin), np.log10(a.emax), a.nE)
    print("rh5_branch.py: dominant Fourier line of each ghost pair's Floquet solution over one period (field components u, u')")
    print(" E       Omega   | pair 1 (smaller |arg mu|): mu-angle/pi  f_dom  (weight)  f_2nd (weight)  f_3rd | pair 2: same")
    rows = []
    for E in grid:
        T, out = spectrum(Kf, E)
        Om = 2 * np.pi / T
        s = f"{E:.4f}  {Om:.4f}  |"
        for (m, f1, w1, f2, w2, f3) in out:
            s += f" {np.angle(m) / np.pi:+.3f}  {f1:+.3f} ({w1:.2f})  {f2:+.3f} ({w2:.2f})  {f3:+.3f} |"
        print(s)
        rows.append((E, T, [(o[1], o[2], o[3], o[4]) for o in out]))
    np.savez(os.path.join(RES, 'rh5_branch.npz'), E=grid, T=np.array([r[1] for r in rows]),
             dom=np.array([[o[0] for o in r[2]] for r in rows]), wdom=np.array([[o[1] for o in r[2]] for r in rows]))


if __name__ == '__main__':
    main()
