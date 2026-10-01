#!/usr/bin/env python3
"""
hdym_uncertainty.py -- uncertainty-exponent estimate of the escape-time
landscape on a 2D section (Grebogi-McDonald-Ott-Yorke).

For M random points x in the section and partners x + delta n (n a random unit
vector in section coordinates), f(delta) = fraction of pairs whose escape
times differ by more than a threshold.  If the level sets of t_esc have a
fractal boundary of dimension d on the section, f(delta) ~ delta^alpha with
alpha = 2 - d for small delta.  alpha = 1: smooth (d = 1).  alpha < 1: fractal.

The threshold matters: t_esc has a smooth component of slope ~ 1/lambda_perp,
so tiny delta always gives tiny differences on a smooth landscape.  The tool
uses several thresholds at once (default: 2, 5, 10 time units and 5 % of the
median t_esc) and reports alpha for each; a fractal boundary gives alpha < 1
for every threshold, a smooth one gives alpha -> 1 as delta -> 0 for every
threshold.

Section definitions and the integrator come from hdym_escape_map.py (keep it
next to this file, and use the same --seed, --dt, --R as the map).

Examples
  python hdym_uncertainty.py --E 0.3 --eps 0.01 --section ym --M 4000 --backend c
  python hdym_uncertainty.py --E 0.3 --eps 0.01 --section ym --M 20000 --backend cupy \
         --delta-min 1e-6 --delta-max 1e-1 --n-delta 11
  # optional figure: map + zoom + f(delta)
  python hdym_uncertainty.py ... --figure results_escape/ym_E0.3_eps0.01_n1024_s1.npz \
                                 results_escape/zoom10x.npz
"""
import argparse
import os
import time

import numpy as np

import hdym_escape_map as em


def sample_pairs(section, E, eps, seed, plane, center, width, M, delta, rng):
    """Return (S_a, S_b, ok) for M random pairs at separation delta."""
    q0, P0, w0, W0 = em.base_point(E, eps, seed)
    e1, e2 = em.unit_pair(plane)
    cu, cv = center
    ua = cu + (rng.random(M) - 0.5) * width
    va = cv + (rng.random(M) - 0.5) * width
    ang = rng.random(M) * 2 * np.pi
    ub, vb = ua + delta * np.cos(ang), va + delta * np.sin(ang)

    def build(u, v):
        if section == 'ghost':
            a = np.sqrt(eps * E)
            w = a * (np.cos(u)[:, None] * e1 + np.sin(u)[:, None] * e2)
            W = a * (np.cos(v)[:, None] * e1 + np.sin(v)[:, None] * e2)
            q = np.broadcast_to(q0, w.shape)
            Pq = np.broadcast_to(P0, w.shape)
            ok = np.ones(M, bool)
        else:
            q = q0 + u[:, None] * e1 + v[:, None] * e2
            Vq = em.V1(q.T)
            K0 = 0.5 * P0 @ P0
            ok = (E - Vq) > 0
            Pq = np.sqrt(np.clip(E - Vq, 0, None) / K0)[:, None] * P0
            w = np.broadcast_to(w0, q.shape)
            W = np.broadcast_to(W0, q.shape)
        p = q + w
        P = Pq - W
        return np.concatenate([p.T, w.T, P.T, W.T], 0), ok

    Sa, oka = build(ua, va)
    Sb, okb = build(ub, vb)
    return Sa, Sb, oka & okb


def fit_alpha(deltas, f, fmin=0.02, fmax=0.9):
    m = (f > fmin) & (f < fmax) & np.isfinite(f)
    if m.sum() < 3:
        return np.nan, np.nan, m
    x, y = np.log(deltas[m]), np.log(f[m])
    A = np.stack([x, np.ones_like(x)], 1)
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    dof = m.sum() - 2
    s2 = (res[0] / dof) if (len(res) and dof > 0) else np.nan
    cov = s2 * np.linalg.inv(A.T @ A) if np.isfinite(s2) else np.full((2, 2), np.nan)
    return coef[0], np.sqrt(cov[0, 0]), m


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--E', type=float, required=True)
    ap.add_argument('--eps', type=float, required=True)
    ap.add_argument('--section', choices=['ghost', 'ym'], default='ym')
    ap.add_argument('--plane', default='xy')
    ap.add_argument('--center', type=float, nargs=2, default=None)
    ap.add_argument('--width', type=float, default=None)
    ap.add_argument('--M', type=int, default=4000, help='pairs per delta')
    ap.add_argument('--delta-min', type=float, default=1e-5)
    ap.add_argument('--delta-max', type=float, default=1e-1)
    ap.add_argument('--n-delta', type=int, default=9)
    ap.add_argument('--thresholds', type=float, nargs='+', default=[2.0, 5.0, 10.0],
                    help='|t_a - t_b| thresholds in time units; 5%% of median t_esc is added')
    ap.add_argument('--dt', type=float, default=0.005)
    ap.add_argument('--tmax', type=float, default=None)
    ap.add_argument('--R', type=float, default=1e4)
    ap.add_argument('--omega-dt-max', type=float, default=0.25)
    ap.add_argument('--backend', choices=['cupy', 'c', 'numpy'], default='c')
    ap.add_argument('--seed', type=int, default=1, help='base-point seed (match the map)')
    ap.add_argument('--pair-seed', type=int, default=12345)
    ap.add_argument('--out', default='results_escape')
    ap.add_argument('--tag', default=None)
    ap.add_argument('--figure', nargs='*', default=None,
                    help='npz files of a map and (optionally) its zoom: makes the combined figure')
    a = ap.parse_args()

    if a.center is None:
        a.center = [np.pi, np.pi] if a.section == 'ghost' else [0.0, 0.0]
    if a.width is None:
        a.width = em.default_width(a.section, a.E)
    if a.tmax is None:
        a.tmax = min(5000.0, 40.0 / (0.027 * a.E ** 2.24))
    # deltas are in section units; for the ghost section (angles) scale relative to 2pi
    deltas = np.logspace(np.log10(a.delta_min), np.log10(a.delta_max), a.n_delta) * a.width
    tag = a.tag or f"unc_{a.section}_E{a.E:g}_eps{a.eps:g}_M{a.M}_s{a.seed}"
    os.makedirs(a.out, exist_ok=True)

    be = em.make_backend(a.backend)
    rng = np.random.default_rng(a.pair_seed)
    rows = []
    all_diff = []
    t0 = time.time()
    for d in deltas:
        Sa, Sb, ok = sample_pairs(a.section, a.E, a.eps, a.seed, a.plane, a.center, a.width,
                                  a.M, d, rng)
        S = np.concatenate([Sa, Sb], 1)
        t = em.integrate_escape(be, S, a.dt, a.tmax, a.R, a.omega_dt_max)
        ta, tb = t[:a.M], t[a.M:]
        good = ok & np.isfinite(ta) & np.isfinite(tb)
        diff = np.abs(ta - tb)[good]
        all_diff.append(diff)
        rows.append((d, good.sum(), np.median(np.concatenate([ta[good], tb[good]])), diff))
        print(f"  delta={d:.3e}  pairs={good.sum()}  median|dt|={np.median(diff):.3g}  "
              f"p90|dt|={np.percentile(diff, 90):.3g}  ({time.time() - t0:.0f}s)", flush=True)

    med_t = np.median([r[2] for r in rows])
    thresholds = list(a.thresholds) + [0.05 * med_t]
    print(f"\nmedian t_esc = {med_t:.1f};  thresholds = {['%.3g' % x for x in thresholds]}")
    F = np.zeros((len(thresholds), len(deltas)))
    for j, th in enumerate(thresholds):
        F[j] = [np.mean(r[3] > th) if len(r[3]) else np.nan for r in rows]
    print("\n  delta        " + "  ".join(f"f(>{th:.3g})" for th in thresholds))
    for i, d in enumerate(deltas):
        print(f"  {d:.3e}    " + "  ".join(f"{F[j, i]:8.4f}" for j in range(len(thresholds))))
    print("\nuncertainty exponent alpha (fit over 0.02 < f < 0.9), boundary dimension d = 2 - alpha:")
    alphas = []
    for j, th in enumerate(thresholds):
        al, se, m = fit_alpha(deltas, F[j])
        alphas.append((al, se))
        note = 'smooth' if al > 0.9 else ('FRACTAL' if al < 0.8 else 'marginal')
        print(f"  threshold {th:7.3g}:  alpha = {al:.3f} +- {se:.3f}   d = {2 - al:.3f}   "
              f"[{m.sum()} points]  {note}" if np.isfinite(al) else
              f"  threshold {th:7.3g}:  not enough points in the scaling window")

    np.savez(os.path.join(a.out, tag + '.npz'), deltas=deltas, F=F, thresholds=thresholds,
             alphas=np.array(alphas), med_t=med_t, E=a.E, eps=a.eps, section=a.section,
             seed=a.seed, M=a.M, width=a.width, center=np.array(a.center))
    print("wrote", os.path.join(a.out, tag + '.npz'))

    # -------- figure
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        maps = [np.load(f, allow_pickle=True) for f in (a.figure or [])]
        ncol = 1 + len(maps)
        fig, axs = plt.subplots(1, ncol, figsize=(5.5 * ncol, 5))
        axs = np.atleast_1d(axs)
        for k, mp in enumerate(maps):
            ax = axs[k]
            te, x1, x2 = mp['t_esc'], mp['ax1'], mp['ax2']
            im = ax.imshow(np.log10(te), origin='lower', cmap='magma', aspect='auto',
                           extent=[x1[0], x1[-1], x2[0], x2[-1]], interpolation='nearest')
            ax.set_title(f"log10 t_esc, {str(mp['section'])} section, {int(mp['n'])}^2, "
                         f"width {float(mp['width']):.3g}", fontsize=9)
            plt.colorbar(im, ax=ax)
            if k + 1 < len(maps):        # draw the next map's frame on this one
                nx = maps[k + 1]
                c, wd = nx['center'], float(nx['width'])
                ax.add_patch(plt.Rectangle((c[0] - wd / 2, c[1] - wd / 2), wd, wd,
                                           fill=False, ec='cyan', lw=1.2))
        ax = axs[-1]
        for j, th in enumerate(thresholds):
            al, se = alphas[j]
            ax.loglog(deltas, F[j], 'o-', ms=4,
                      label=f"|dt|>{th:.3g}: alpha={al:.2f}+-{se:.2f}" if np.isfinite(al)
                      else f"|dt|>{th:.3g}")
        ax.set_xlabel('delta (section units)'); ax.set_ylabel('f(delta)')
        ax.set_title(f"uncertainty exponent  E={a.E:g} eps={a.eps:g} {a.section}", fontsize=9)
        ax.legend(fontsize=7)
        plt.tight_layout()
        png = os.path.join(a.out, tag + '.png')
        plt.savefig(png, dpi=150)
        print("wrote", png)
    except Exception as e:                        # noqa
        print("plot skipped:", e)


if __name__ == '__main__':
    main()
