"""rh6_polarisation.py -- rabbit hole 6: which of block 12's two ghost-sector pairs is "the ghost's", and in which
polarisation does each live, measured from the eigenvectors B005 already computes.

The repository defines no polarisation basis for the off-diagonal components (hdym_rung1_lyap's docstring names
"the longitudinal polarisation of the massive Lee-Wick vector" without writing it down), so three measurable
contents are taken for each pair, from the eigenvector v of N = S~ Phi_T/2 at the orbit's reference point (t = 0,
the turning point y = 0, the point at which the sector basis is the reflection's eigenbasis) and, since the
contents move along the orbit, also at t = T/4 (propagated):

  component content  f12 = |u12|^2 / (|u12|^2 + |u21|^2) using the field-level entries (u and u') of v;
  sector content     f+  = |u+|^2 / (|u+|^2 + |u-|^2), u+- = (u12 +- u21)/sqrt 2 -- the reflection's eigenbasis;
  local frequency    w^2 = |u''| / |u| (field entries): a heavy mode (the massive Lee-Wick partner, mass 1 plus
                     the background's contribution) has |u''| ~ w^2 |u| with w^2 of order one or more, a light
                     mode has w^2 ~ |xy|; this is the Pais-Uhlenbeck light/heavy split read off the eigenvector.

For each pair B005's Krein sign of the member at positive half-period angle is also printed.  The sign of a mode
is the sign of its conserved quadratic form, i.e. of energy times frequency; a pair whose POSITIVE-frequency member
carries NEGATIVE sign is a ghost mode, so "which pair is the ghost's" needs the true (unwrapped) frequency sign,
which the mod-2pi angle does not give (see notes/RH2_ladder.md on the rotation number).  What can be said without it
is printed: the contents, and whether the sign of the member at |angle| < pi/2 (where the branch is unambiguous)
is + or -.

Usage (from the repo root): python rabbit_holes/ym/code/rh6_polarisation.py [--nE 100] > rabbit_holes/ym/results/rh6_polarisation.txt
Cost: ~2.5 s per energy.  Added 2026-10-21 (rabbit-hole session).
"""
import argparse
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'islands', 'code'))
import dr_offdiag_floquet as base            # noqa: E402
import krein_b005 as kb                      # noqa: E402

RES = os.path.join(HERE, '..', 'results')
IDX = base.block_idx('12')
STILDE = np.diag(np.array(kb.STILDE['12'] * 4, dtype=float))
# block state order: (u12, u21, u12', u21', u12'', u21'', u12''', u21''')


def contents(v):
    u = v[0:2]
    ud = v[2:4]
    udd = v[4:6]
    w2 = np.linalg.norm(udd) / max(np.linalg.norm(u), 1e-300)
    f12 = (abs(u[0]) ** 2 + abs(ud[0]) ** 2) / (np.linalg.norm(u) ** 2 + np.linalg.norm(ud) ** 2)
    up, um = (u[0] + u[1]) / np.sqrt(2), (u[0] - u[1]) / np.sqrt(2)
    upd, umd = (ud[0] + ud[1]) / np.sqrt(2), (ud[0] - ud[1]) / np.sqrt(2)
    fp = (abs(up) ** 2 + abs(upd) ** 2) / (np.linalg.norm(u) ** 2 + np.linalg.norm(ud) ** 2)
    return f12, fp, w2


def propagate(Kf, E, V, tfrac):
    q0, p0, T = base.orbit(E)
    nv = V.shape[1]

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, 1.0), dtype=float)[np.ix_(IDX, IDX)]
        X = s[6:].reshape(2 * nv, 8)
        return np.concatenate([p, -base.gradV(q), (X @ K.T).ravel()])
    s0 = np.concatenate([q0, p0, np.concatenate([V.real.T, V.imag.T]).ravel()])
    sT = solve_ivp(rhs, (0.0, tfrac * T), s0, **base.RT).y[:, -1]
    X = sT[6:].reshape(2 * nv, 8)
    return (X[:nv] + 1j * X[nv:]).T, sT[:3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nE', type=int, default=100)
    ap.add_argument('--emin', type=float, default=0.05)
    ap.add_argument('--emax', type=float, default=1.0)
    a = ap.parse_args()
    t0 = time.time()
    print(f"rh6_polarisation.py  python {sys.version.split()[0]}  numpy {np.__version__}; nE={a.nE} log-spaced in [{a.emin}, {a.emax}], RT={base.RT}")
    print("columns: E | per pair (sorted by |angle| of N, pair 1 = smaller): angle/pi of the member with positive angle, its Krein sign,"
          " f12, f+ and w^2 at t=0 ; f12, f+, w^2 at t=T/4 | background x,y at T/4")
    Kf = base.load_K()
    Of, _ = kb.build_omega()
    grid = np.logspace(np.log10(a.emin), np.log10(a.emax), a.nE)
    rows = []
    for E in grid:
        Ph, M, (q0, p0), _, T = kb.half_and_full(Kf, E, '12')
        N = STILDE @ Ph
        Om = np.array(Of(*q0, *p0), dtype=float)[np.ix_(IDX, IDX)]
        lam, V = np.linalg.eig(N)
        keep = kb.ghost_four(lam)
        lam4, V4 = lam[keep], V[:, keep]
        # pairs: conjugate partners; take the member with positive angle (or, off circle, the larger |lambda|)
        pairs = []
        used = set()
        for i in range(4):
            if i in used:
                continue
            j = min((abs(lam4[k] - np.conj(lam4[i])), k) for k in range(4) if k != i and k not in used)[1]
            used |= {i, j}
            pos = i if np.angle(lam4[i]) >= np.angle(lam4[j]) else j
            pairs.append(pos)
        pairs.sort(key=lambda k: abs(np.angle(lam4[k])))
        V2 = V4[:, pairs]
        Vq, qq = propagate(Kf, E, V2, 0.25)
        row = [E]
        txt = []
        for n, k in enumerate(pairs):
            v = V4[:, k]
            on = abs(abs(lam4[k]) - 1) < 1e-6
            sgn = np.sign(kb.krein_sign(Om, v)[0]) if on else 0
            f12, fp, w2 = contents(v)
            f12q, fpq, w2q = contents(Vq[:, n])
            row += [np.angle(lam4[k]) / np.pi, abs(lam4[k]), sgn, f12, fp, w2, f12q, fpq, w2q]
            txt.append(f"ang {np.angle(lam4[k]) / np.pi:+.3f} {'|l|=%.3f' % abs(lam4[k]) if not on else 'sign %+d' % sgn:10s} f12 {f12:.2f} f+ {fp:.2f} w2 {w2:6.2f} ; T/4: f12 {f12q:.2f} f+ {fpq:.2f} w2 {w2q:6.2f}")
        rows.append(row)
        print(f"E={E:.5f} | " + " || ".join(txt) + f" | x,y {qq[0]:+.2f},{qq[1]:+.2f}")
    R = np.array(rows)
    np.savez(os.path.join(RES, 'rh6_polarisation.npz'), rows=R)

    # summary statistics on the circle
    print("\n== summary over the scan (both pairs on the circle) ==")
    on = (np.abs(R[:, 2] - 1) < 1e-6) & (np.abs(R[:, 11] - 1) < 1e-6)
    for n, off in ((1, 1), (2, 10)):
        f12, fp, w2 = R[on, off + 3], R[on, off + 4], R[on, off + 5]
        f12q, fpq, w2q = R[on, off + 6], R[on, off + 7], R[on, off + 8]
        print(f"  pair {n} ({"the one with the smaller |angle|" if n == 1 else "the larger"}): at t=0  f12 median {np.median(f12):.2f} [{np.percentile(f12, 10):.2f},{np.percentile(f12, 90):.2f}],"
              f" f+ median {np.median(fp):.2f} [{np.percentile(fp, 10):.2f},{np.percentile(fp, 90):.2f}], w^2 median {np.median(w2):.2f} [{np.percentile(w2, 10):.2f},{np.percentile(w2, 90):.2f}];"
              f" at T/4 f12 {np.median(f12q):.2f}, f+ {np.median(fpq):.2f}, w^2 {np.median(w2q):.2f}")
    s1, s2 = R[on, 3], R[on, 12]
    a1, a2 = np.abs(R[on, 1]), np.abs(R[on, 10])
    print(f"  Krein sign of the positive-angle member: pair 1 is + at {100 * (s1 > 0).mean():.0f} % of on-circle energies, pair 2 at {100 * (s2 > 0).mean():.0f} %")
    m = a1 < 0.5
    m2 = a2 < 0.5
    pct = lambda x: f"{100 * (x > 0).mean():.0f} %" if len(x) else "n/a"
    print(f"  restricted to |angle| < pi/2 (branch unambiguous): pair 1 + at {pct(s1[m])} ({m.sum()} energies); pair 2 + at {pct(s2[m2])} ({m2.sum()} energies)")
    print(f"DONE {time.time() - t0:.0f}")


if __name__ == '__main__':
    main()
