#!/usr/bin/env python3
"""
hdym_island.py -- escape times and base-flow regularity in the neighbourhood of the
Dahlqvist-Russberg (DR) stable periodic orbit, in the full 12-dimensional Lee-Wick
Yang-Mills system (three fields + three ghosts), compared with the chaotic sea.

Conventions: exactly those of hdym_escape_map.py / hdym_precision.py.
  C1  state (12, N): p, w, P, W;  q = p - w;  pdot = P, wdot = -W;  kick_nl of hdym_escape_map
  C2  E = |P|^2/2 + V(q), V = (x^2y^2 + y^2z^2 + z^2x^2)/2;  ghost |w|^2 + |W|^2 = 2 eps E;
      p = q + w, P = P_q - W
  C3  DR orbit: q* = (3.14640122769837, 0, 0), P* = (0.00179319341881, sqrt(1 - px*^2), 0)
      at E = 1/2 (certified, cert16.json); scaled to E by q -> lam q, P -> lam^2 P,
      lam = (E/0.5)^(1/4).  Symmetry copies: cyclic permutations of (x,y,z), rotations by
      90 degrees in the plane, and time reversal.
  C4  escape: any |component| > R (default 1e4), dt = 0.005, omega_dt_max = 0.25.

For each radius r (log-spaced) we sample N points in a ball of radius r around the DR phase
point (q and P both perturbed, P then rescaled to E), seed the ghost as in C2 with a random
orientation, and integrate a PAIR (partner displaced by delta in q, shared sub-step schedule,
hdym_precision.CPairs) to get
    t_esc(a), t_esc(b)       escape times (C4)
    sep(t_k)                 phase-space separation at checkpoint times  -> FTLE proxy
A control sample of random sea points at the same E (hdym_escape_map.base_point with
different seeds) is run with identical settings.

Predictions (from the 2-dof analysis):
  * eps small enough that sqrt(eps E) << island width (~1e-4):
      - E in a Floquet tongue (e.g. 1.0): inside the island t_esc is finite, SHORT and smooth
        (plateau), FTLE ~ 0;  the halo (r ~ 2e-4 .. 1e-3) has the broadest t_esc spread.
      - E between tongues (e.g. 0.3): inside the island the ghost is bounded -> t_esc = inf
        (hole of positive measure), FTLE ~ 0.
  * eps = 0.01 (paper): the ghost kicks the base by ~0.07 >> island width; the island may be
    invisible.  Whether any trace survives is the question.

Usage
  python hdym_island.py --E 1.0 --eps 1e-10 --N 200 --backend c
  python hdym_island.py --E 0.3 --eps 1e-10 --N 200 --tmax 3000
  python hdym_island.py --E 1.0 --eps 0.01  --N 200
  python hdym_island.py --E 1.0 --eps 1e-10 --copy 7      # a different symmetry copy
"""
import argparse, json, os, sys, time
import numpy as np
# the production modules live at the repository root; this file lives in islands/code/
_here = os.path.dirname(os.path.abspath(__file__))
for _p in (_here, os.path.join(_here, '..'), os.path.join(_here, '..', '..')):
    if _p not in sys.path: sys.path.insert(0, _p)
import hdym_escape_map as em
import hdym_precision as hp

# ---------------------------------------------------------------- C3: the DR orbit and its copies
XS, PXS = 3.14640122769837, 0.00179319341881
PYS = np.sqrt(1.0 - PXS**2)

def dr_point(E, copy=0):
    """One of the 24 symmetry copies of the DR phase point at energy E.
    copy = 8*perm + 2*rot + rev :  perm in 0..2 (cyclic (x,y,z)), rot in 0..3 (90-degree
    rotations in the orbit plane), rev in 0..1 (time reversal P -> -P)."""
    lam = (E / 0.5) ** 0.25
    q = np.array([XS, 0.0, 0.0]) * lam
    P = np.array([PXS, PYS, 0.0]) * lam**2
    perm, rot, rev = copy // 8, (copy % 8) // 2, copy % 2
    for _ in range(rot):                        # (x, y) -> (-y, x) in the plane
        q = np.array([-q[1], q[0], q[2]]); P = np.array([-P[1], P[0], P[2]])
    if rev: P = -P
    q = np.roll(q, perm); P = np.roll(P, perm)  # cyclic permutation of the fields
    return q, P

def rescale_to_E(q, P, E):
    """scale P so that |P|^2/2 + V(q) = E (as section_ym does); returns P, ok."""
    Vq = em.V1(q.T)
    K = 0.5 * np.sum(P * P, axis=1)
    ok = (E - Vq) > 0
    s = np.sqrt(np.clip(E - Vq, 0, None) / K)
    return P * s[:, None], ok

def seed_ghost(N, E, eps, rng):
    """C2: random ghost orientation with |w|^2 + |W|^2 = 2 eps E."""
    nw = rng.standard_normal((N, 3)); nW = rng.standard_normal((N, 3))
    nrm = np.sqrt(np.sum(nw * nw, 1) + np.sum(nW * nW, 1))
    amp = np.sqrt(2.0 * eps * E) / nrm
    return nw * amp[:, None], nW * amp[:, None]

def assemble(q, P, w, W):
    """C1/C2 layout (12, N)."""
    p = q + w
    Pt = P - W
    return np.concatenate([p.T, w.T, Pt.T, W.T], 0)

PLANAR = [False]      # --planar: restrict to the exact invariant subsystem z = p_z = w_z = W_z = 0

def ball_sample(qc, Pc, r, N, E, eps, rng, delta):
    """N points in a ball of radius r around (qc, Pc), partner displaced by delta in q."""
    u = rng.standard_normal((N, 6)); u /= np.linalg.norm(u, axis=1)[:, None]
    rad = r * rng.random(N) ** (1 / 6)            # uniform in the 6-ball
    dq, dP = u[:, :3] * rad[:, None], u[:, 3:] * rad[:, None]
    q = qc + dq; P = Pc + dP
    w, W = seed_ghost(N, E, eps, rng)
    n = rng.standard_normal((N, 3)); n /= np.linalg.norm(n, axis=1)[:, None]
    if PLANAR[0]:
        for arr in (q, P, w, W, n): arr[:, 2] = 0.0
        n /= np.linalg.norm(n, axis=1)[:, None]
    P, ok = rescale_to_E(q, P, E)
    Sa = assemble(q, P, w, W)
    qb = q + delta * n
    Pb, okb = rescale_to_E(qb, P, E)
    Sb = assemble(qb, Pb, w, W)
    return Sa, Sb, ok & okb

def sea_sample(N, E, eps, rng, delta, seed0=1000):
    qs, Ps = [], []
    for s in range(N):
        q, P, _, _ = em.base_point(E, eps, seed0 + s)
        qs.append(q); Ps.append(P)
    q, P = np.array(qs), np.array(Ps)
    w, W = seed_ghost(N, E, eps, rng)
    Sa = assemble(q, P, w, W)
    n = rng.standard_normal((N, 3)); n /= np.linalg.norm(n, axis=1)[:, None]
    qb = q + delta * n
    Pb, okb = rescale_to_E(qb, P, E)
    Sb = assemble(qb, Pb, w, W)
    return Sa, Sb, okb

def summarize(label, out, ok, chk_t, delta):
    ta, tb = out[0], out[1]
    esc = ok & (ta > 0)
    unesc = ok & (ta < 0)
    seps = [out[2 + j][ok & (out[2 + j] > 0)] for j in range(len(chk_t))]
    ftle = []
    for j, t in enumerate(chk_t):
        s = seps[j]
        ftle.append(np.median(np.log(s / delta)) / t if s.size else np.nan)
    te = ta[esc]
    med = np.median(te) if te.size else np.nan
    p10 = np.percentile(te, 10) if te.size else np.nan
    p90 = np.percentile(te, 90) if te.size else np.nan
    spread = (np.log10(p90 / p10) if te.size and p10 > 0 else np.nan)
    print(f"  {label:>14s}  n={ok.sum():4d}  escaped {esc.sum()/max(1,ok.sum()):.3f}  "
          f"t_esc med {med:8.1f} [p10 {p10:7.1f}, p90 {p90:8.1f}] log10 spread {spread:5.2f}   "
          f"FTLE@{chk_t[0]:g} {ftle[0]:8.4f}  @{chk_t[-1]:g} {ftle[-1]:8.4f}", flush=True)
    return dict(label=label, n=int(ok.sum()), frac_esc=float(esc.sum() / max(1, ok.sum())),
                frac_unesc=float(unesc.sum() / max(1, ok.sum())), t_med=float(med),
                t_p10=float(p10), t_p90=float(p90), ftle=[float(x) for x in ftle],
                t_esc_a=ta[ok].tolist())

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--E', type=float, default=1.0)
    ap.add_argument('--eps', type=float, default=1e-10)
    ap.add_argument('--N', type=int, default=200, help='points per radius shell')
    ap.add_argument('--radii', type=float, nargs='+',
                    default=[1e-6, 1e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 1e-2])
    ap.add_argument('--copy', type=int, default=0, help='symmetry copy of the DR orbit (0..23)')
    ap.add_argument('--delta', type=float, default=1e-10, help='pair separation (FTLE proxy)')
    ap.add_argument('--t-chk', type=float, nargs='+', default=[31.05, 155.3, 621.1],
                    help='checkpoint times for the separation (1, 5, 20 DR periods at E=1/2)')
    ap.add_argument('--dt', type=float, default=0.005)
    ap.add_argument('--tmax', type=float, default=None)
    ap.add_argument('--R', type=float, default=1e4)
    ap.add_argument('--omega-dt-max', type=float, default=0.25)
    ap.add_argument('--backend', choices=['c', 'cupy'], default='c')
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--out', default='results_island')
    ap.add_argument('--tag', default=None)
    ap.add_argument('--planar', action='store_true',
                    help='restrict everything to the exact invariant 2-dof subsystem z = p_z = w_z = W_z = 0 '
                         '(ghost in the plane too); the sea control stays 3-dof')
    a = ap.parse_args()
    PLANAR[0] = a.planar
    if a.tmax is None:
        a.tmax = min(5000.0, 40.0 / (0.027 * a.E ** 2.24))        # hdym_escape_map default
    lam = (a.E / 0.5) ** 0.25
    chk_t = [t / lam for t in a.t_chk]                              # DR period scales as E^-1/4
    chk = [int(round(t / a.dt)) for t in chk_t]
    nmax = int(np.ceil(a.tmax / a.dt))
    tag = a.tag or f"island_E{a.E:g}_eps{a.eps:g}_c{a.copy}_N{a.N}" + ("_planar" if a.planar else "")
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    be = hp.CupyPairs() if a.backend == 'cupy' else hp.CPairs()
    qc, Pc = dr_point(a.E, a.copy)
    print(f"{tag}: DR copy {a.copy}: q*={qc}, P*={Pc}, E check={0.5*Pc@Pc + em.V1(qc):.6f}; "
          f"ghost amplitude sqrt(2 eps E)={np.sqrt(2*a.eps*a.E):.2e}; tmax={a.tmax:.0f}; "
          f"FTLE checkpoints t={[round(t,1) for t in chk_t]}", flush=True)
    print(f"  island half-width (2-dof, E=1/2) ~ 2e-4 in x; scaled ~ {2e-4*lam:.1e}")
    results = []
    t0 = time.time()
    for r in a.radii:
        Sa, Sb, ok = ball_sample(qc, Pc, r, a.N, a.E, a.eps, rng, a.delta)
        out = be.run(np.concatenate([Sa, Sb], 0), a.dt, nmax, a.R, a.omega_dt_max, 0, 1, chk)
        results.append(dict(r=r, **summarize(f"r={r:.0e}", out, ok, chk_t, a.delta)))
        print(f"      ({time.time()-t0:.0f}s)", flush=True)
    Sa, Sb, ok = sea_sample(a.N, a.E, a.eps, rng, a.delta)
    out = be.run(np.concatenate([Sa, Sb], 0), a.dt, nmax, a.R, a.omega_dt_max, 0, 1, chk)
    results.append(dict(r=np.inf, **summarize("sea control", out, ok, chk_t, a.delta)))
    with open(os.path.join(a.out, tag + '.json'), 'w') as fh:
        json.dump(dict(args=vars(a), chk_t=chk_t, results=results), fh)
    print("wrote", os.path.join(a.out, tag + '.json'))
    print("\nReading guide: island => FTLE ~ 0 and (tongue) short smooth t_esc / (off-tongue) "
          "unescaped;  halo => FTLE between island and sea, widest t_esc spread;  sea => FTLE ~ "
          "lambda_YM ~ 0.38 E^1/4 and the paper's t_esc distribution.")

if __name__ == '__main__':
    main()
