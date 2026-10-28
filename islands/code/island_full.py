"""island_full.py -- rabbit hole 9: is the E = 0.3 island a refuge in the FULL nine-component Lee-Wick system?  (No.)

hdym_island.py (12-d, colour-diagonal ghost only) finds the island around the DR point to be a non-escaping hole at
E = 0.3 for eps <= 1e-8 (rabbit hole 10).  The off-diagonal ghost is Floquet-unstable on that orbit (block 23, real
multiplier 1.552 at E = 0.3, rate 0.0125; islands/README off-diagonal item), so the island is "not a refuge in the
full theory" -- until this script, a linear statement.  This script makes it nonlinear: the same ball sampling
(hdym_island.ball_sample / sea_sample, same shells, same rescaling to E), each point converted to the 36-component
state of ym/hdym_full_derive.py and integrated with ym/csrc/hdym_full_escape.c (offdiag_escape.Kernel), with the
six off-diagonal components seeded at a tiny Gauss-projected amplitude --delta-od (default 1e-8, three decades below
the diagonal ghost's sqrt(2 eps E) = 7.7e-5 at eps = 1e-8) or, with --delta-od 0, left to the integrator's truncation.

Conversion (hdym_core conventions: q = p - w, w = q'' + grad V(q), w' = -W, q' = P + W = P_q):
    q0 = q,  q1 = P_q,  q2 = -grad V(q) + w,  q3 = -Hess V(q) P_q - W,
placed at the diagonal entries m = 0, 4, 8 of each derivative block s[9 k + m].

Per shell: escaped fraction by tmax, median t_esc, the sector largest at escape, the first time |w| reaches 1 and the
first time each off-diagonal block norm reaches 1 (tcross), with the sea control.  Resumable per shell (--resume).
Output <results>/<tag>.npz + the log.  Checks: the kernel's own CHECK lines; the 36-state's energy and Gauss residual
per point (printed per shell: max |G| at t = 0 must be O(delta_od^2), and the Ostrogradsky energy must equal
E + E_ghost of the 12-d state to 1e-10 -- the conversion check).

Usage (Spyder, from the repo root; ym/hdym_full_derive.py run once):
    !"{sys.executable}" -B -u islands\\code\\island_full.py --E 0.3 --eps 1e-8 --tag island_full_E0.3_eps1e-8 --resume >> islands\\results\\island_full_E0.3_eps1e-8.txt 2>&1
    !"{sys.executable}" -B -u islands\\code\\island_full.py --E 0.3 --eps 1e-8 --delta-od 0 --tag island_full_E0.3_eps1e-8_od0 --resume >> islands\\results\\island_full_E0.3_eps1e-8_od0.txt 2>&1
Cost: 8 shells x 200 + 200 sea = 1800 trajectories, tmax 5000, ~1 s each single-threaded: ~10-20 min on 24 threads
per run (estimate; the committed runs took 2999 s and 632 s, DONE lines of the two logs).  Added 2026-10-22 (rabbit-hole session); runs of 23 Oct (results/island_full_E0.3_eps1e-8.*, _od0.*): with
the off-diagonal seed every core shell r = 1e-6 .. 1e-4 escapes, 200 of 200, median t_esc 1340-1360 (ln(1e8)/0.0125
~ 1500), the diagonal ghost the escaping sector in 91-96 % of cells and block 23 in the rest; with --delta-od 0 the
island is intact, 0 of 200 in the four core shells, halo and sea as in hole 10.  Not a refuge (islands/README.md,
Results; rabbit_holes/ym/notes/RH9_island_full.md).
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
for _p in (HERE, os.path.join(ROOT, 'ym')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import hdym_island as isl                     # noqa: E402  (dr_point, ball_sample, sea_sample, rescale_to_E)
import hdym_core as hc                        # noqa: E402
import offdiag_escape as OE                   # noqa: E402  (Kernel, Model, gauss_projector, od_index)

SECTOR_NAMES = OE.SECTOR_NAMES


def hessV(q):
    x, y, z = q
    return np.array([[y * y + z * z, 2 * x * y, 2 * x * z], [2 * x * y, z * z + x * x, 2 * y * z], [2 * x * z, 2 * y * z, x * x + y * y]])


def to_full(Sa):
    """(12, N) hdym_core state (p, w, P, W) -> (N, 36) full state with the off-diagonal sector zero."""
    p, w, Pt, W = Sa[0:3].T, Sa[3:6].T, Sa[6:9].T, Sa[9:12].T
    q = p - w
    Pq = Pt + W
    N = q.shape[0]
    S = np.zeros((N, 36))
    for n in range(N):
        g = hc.gradV(q[n:n + 1])[0]
        H = hessV(q[n])
        q2 = -g + w[n]
        q3 = -(H @ Pq[n]) - W[n]
        for i in range(3):
            S[n, 4 * i] = q[n, i]
            S[n, 9 + 4 * i] = Pq[n, i]
            S[n, 18 + 4 * i] = q2[i]
            S[n, 27 + 4 * i] = q3[i]
    return S


def seed_offdiag(S, model, delta, rng):
    """add a Gauss-projected off-diagonal seed of norm delta to every row (projector built at each point's own q, p)."""
    if delta <= 0:
        return S
    S = S.copy()
    idx = OE.od_index()
    for n in range(S.shape[0]):
        P = OE.gauss_projector(model, S[n], 'all')
        v = P @ rng.normal(size=24)
        v *= delta / np.linalg.norm(v)
        S[n, idx] += v
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--E', type=float, default=0.3)
    ap.add_argument('--eps', type=float, default=1e-8)
    ap.add_argument('--delta-od', type=float, default=1e-8)
    ap.add_argument('--N', type=int, default=200)
    ap.add_argument('--radii', type=float, nargs='+', default=[1e-6, 1e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 1e-2])
    ap.add_argument('--copy', type=int, default=0)
    ap.add_argument('--tmax', type=float, default=5000.0)
    ap.add_argument('--R', type=float, default=1e4)
    ap.add_argument('--rtol', type=float, default=1e-11)
    ap.add_argument('--atol', type=float, default=1e-13)
    ap.add_argument('--hmax', type=float, default=0.25)
    ap.add_argument('--nchk', type=int, default=100)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--tag', default='island_full')
    ap.add_argument('--results', default=os.path.join(HERE, '..', 'results'))
    ap.add_argument('--resume', action='store_true')
    a = ap.parse_args()
    out = os.path.join(os.path.normpath(a.results), a.tag + '.npz')
    t0 = time.time()
    print(f"island_full.py  python {sys.version.split()[0]}  numpy {np.__version__}")
    print("run parameters: " + " ".join(f"{k}={v}" for k, v in vars(a).items() if k != 'results'))
    kern = OE.Kernel()
    model = OE.Model()
    rng = np.random.default_rng(a.seed)
    qc, Pc = isl.dr_point(a.E, a.copy)
    lam = (a.E / 0.5) ** 0.25
    print(f"DR copy {a.copy} at E={a.E}: q*={qc}, P*={Pc}; diagonal ghost amplitude sqrt(2 eps E) = {np.sqrt(2 * a.eps * a.E):.2e}; off-diagonal seed {a.delta_od:.1e}; "
          f"island half-width (2-dof, E=1/2) ~ 2e-4 in x, scaled ~ {2e-4 * lam:.1e}; tmax={a.tmax:.0f}")
    labels = [f"r={r:g}" for r in a.radii] + ['sea control']
    nrow = len(labels)
    tesc = np.full((nrow, a.N), np.nan)
    sector = np.full((nrow, a.N), -1, dtype=np.int32)
    status = np.full((nrow, a.N), -1, dtype=np.int32)
    tcross = np.full((nrow, a.N, 5), -1.0)
    ok_all = np.zeros((nrow, a.N), dtype=bool)
    start = 0
    if a.resume and os.path.exists(out):
        d = np.load(out, allow_pickle=True)
        same = all(str(d[k]) == str(getattr(a, k)) for k in ('E', 'eps', 'delta_od', 'N', 'tmax', 'seed', 'copy') if k in d.files)
        if same:
            start = int(d['rows_done'])
            for k, arr in (('tesc', tesc), ('sector', sector), ('status', status), ('tcross', tcross), ('ok', ok_all)):
                arr[:start] = d[k][:start]
            print(f"--resume: continuing from row {start} of {nrow}", flush=True)
    for i, lab in enumerate(labels):
        # the sampling consumes the rng in the same order whether or not the row is skipped, so --resume is exact
        if i < len(a.radii):
            Sa, Sb, ok = isl.ball_sample(qc, Pc, a.radii[i], a.N, a.E, a.eps, rng, 1e-10)
        else:
            Sa, Sb, ok = isl.sea_sample(a.N, a.E, a.eps, rng, 1e-10)
        S = to_full(Sa)
        S = seed_offdiag(S, model, a.delta_od, rng)
        if i < start:
            continue
        # conversion checks on this shell
        E12 = hc.hamiltonian((Sa[0:3].T, Sa[3:6].T, Sa[6:9].T, Sa[9:12].T))
        Efull = np.array([model.energy(S[n]) for n in range(min(a.N, 8))])
        gres = max(np.linalg.norm(model.gauss(S[n])) for n in range(min(a.N, 8)))
        dE = np.abs(Efull - E12[:len(Efull)]).max() / a.E
        te, se, ch, dg, _, st = kern.escape(S, 1.0, a.tmax, a.R, a.rtol, a.atol, a.hmax, a.nchk)
        tesc[i], sector[i], status[i], tcross[i], ok_all[i] = te, se, st, kern.last_tcross, ok
        esc = ok & (st >= 0)
        secs = ",".join(f"{SECTOR_NAMES[g]}:{int((se[esc] == g).sum())}" for g in range(4) if (se[esc] == g).any())
        t1 = kern.last_tcross[:, 0]
        tb = np.where(kern.last_tcross[:, 1:4] > 0, kern.last_tcross[:, 1:4], np.inf).min(axis=1)
        print(f"  {lab:>12s}  n={int(ok.sum()):4d}  escaped {esc.sum() / max(1, ok.sum()):.3f}  t_esc med {np.median(te[esc]) if esc.any() else float('nan'):8.1f} "
              f"[p10 {np.percentile(te[esc], 10) if esc.any() else float('nan'):8.1f}, p90 {np.percentile(te[esc], 90) if esc.any() else float('nan'):8.1f}]  sectors [{secs}]  "
              f"|w|=1 reached {int((t1 > 0).sum()):3d} (med t {np.median(t1[t1 > 0]) if (t1 > 0).any() else float('nan'):7.0f})  a block norm=1 reached {int(np.isfinite(tb).sum()):3d} (med t {np.median(tb[np.isfinite(tb)]) if np.isfinite(tb).any() else float('nan'):7.0f})  "
              f"CHECK conversion: |E_full - E_12|/E {dE:.1e}, |G| {gres:.1e}  ({time.time() - t0:.0f}s)", flush=True)
        np.savez(out, E=a.E, eps=a.eps, delta_od=a.delta_od, N=a.N, tmax=a.tmax, seed=a.seed, copy=a.copy, radii=np.array(a.radii), labels=np.array(labels),
                 tesc=tesc, sector=sector, status=status, tcross=tcross, ok=ok_all, rows_done=i + 1, rows=nrow)
    fails = int((tesc == -2).sum())
    print(f"\nCHECK no integrator failures: {fails} {'PASS' if fails == 0 else 'FAIL'}")
    print("Reading: hdym_island's hole (0 escaped in the core shells at eps = 1e-8 for the diagonal ghost) is a refuge in the full theory iff the core shells")
    print("still do not escape here; the off-diagonal block norms' first crossing of 1 and the escaping sector say which mechanism empties it if they do.")
    print(f"DONE {time.time() - t0:.0f}")


if __name__ == '__main__':
    main()
