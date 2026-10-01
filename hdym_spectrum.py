#!/usr/bin/env python3
"""
hdym_spectrum.py (v2.1, parallel) -- lambda_perp(E) from the power spectrum of Hess V along
pure Yang-Mills trajectories (no ghost integrated; no long runs).

Why.  Linearised, the ghost is a unit-frequency oscillator with a fluctuating
mass matrix,  wdd = -(1 + K(t)) w,  K = Hess V(p(t)).  For weak driving the
parametric-resonance result for a scalar oscillator
      xdd + (1 + k(t)) x = 0     =>     lambda = S_k(2) / 8,
with S_k the two-sided spectral density  S_k(om) = int <k(t)k(0)> e^{i om t} dt.
The YM scaling symmetry makes one spectrum serve every energy: with p(t) at
energy E written as E^{1/4} P(E^{1/4} t) for a trajectory P at E = 1,
      K_E(t) = E^{1/2} K_1(E^{1/4} t),
so with s(u) a spectral density computed once at E = 1
      lambda(E) ~ E^{3/4} s(2 E^{-1/4}) / 8 .
The ghost frequency in scaled units is u = 2 E^{-1/4}: lowering E moves the
resonance out along the tail of s(u).

Two mechanisms, two spectra.  In the instantaneous eigenframe of K the ghost
is three oscillators with frequencies sqrt(1 + lam_i(t)), coupled by the
frame rotation Omega_ij = e_i . de_j/dt.
  * Single-mode parametric driving: spectrum of the tracked eigenvalues
    lam_i(t) at u = 2 E^{-1/4}.  This is the low-E mechanism and gives the
    E^2.2 law (the tail of s(u) is ~ u^-6).  The basis-dependent diagonal
    entries K_xx = y^2 + z^2 approximate it because in a channel the diagonal
    entry IS the lowest eigenvalue.
  * Mode mixing: spectrum of Omega_ij, resonant when it has power at the
    difference/sum of two mode frequencies -- i.e. in the interior where the
    three modes are degenerate.  This is the high-E mechanism: measured
    vec/diag ratio ~1.5 at E <= 0.03, rising through 0.05-0.2 to ~20.
Which basis is physical depends on E.  The eigenframe description holds when
the ghost follows it adiabatically, i.e. when the frame rotates slowly
compared with the mode SPLITTING (~ E^{1/2}); interior rotation rates are
~ E^{1/4}, so at low E the ghost polarisation is effectively frozen in a fixed
basis while the interior Hessian rotates under it, and only in channels does
the frame lock (to the axis).  So the fixed-basis diagonal entries are not
merely a stand-in at low E; the tracked-mode spectrum is the alternative
description, and --calibrate against the --ghost diag scan decides which fits.
Only the lowest mode is tracked (see track_lowest); sorted eigenvalues have
kinks at crossings that manufacture a u^-2 tail and are kept only for
comparison.

Estimator.  --taper dpss (default): Thomson multitaper, but NOT the textbook
NW = 4, K = 2NW-1: the tail at u ~ 20 sits 6-8 decades below the peak, and the
higher-order Slepians at NW = 4 leak 1e-4 .. 6e-2 of the peak (their
concentration eigenvalues fall short of 1 by that much), which buries the
tail under a spurious u^-2 (step-discontinuity) spectrum.  Default NW = 16
with K = 8 = NW/2 tapers: every taper concentrated to machine precision, edge
values < 1e-15, resolution bandwidth NW/T_seg ~ 0.04 in u -- irrelevant for a
smooth spectrum.  Gain over a single Blackman-Harris taper: the 8 tapers use
the whole segment instead of ~half and do not down-weight bursts near segment
edges (the tail IS bursts: channel passages), roughly 3-4x in effective data.
--taper bh is the cross-check; they must agree over u = 2-20.  Neither fixes
the real limit past u ~ 20, which is the number of deep channel passages in
the ensemble (exceedance ~ D^-2): only N x T does.

Outputs (in --out):
  spectrum.npz    om; S_diag, S_off (entries); S_eig (3 columns: tracked
                  lowest mode, mean of the other two, sorted lowest); S_rot
                  (rotation rate |de0/dt| of the tracked mode);
                  S_diag_shallow / S_diag_deep (depth split); depth
  prediction.csv  E, u, lam_diag, lam_eig, lam_eigmin, (calibrated column)
  spectrum.png

Usage
  python hdym_spectrum.py --N 2048 --T 20000 --t-burn 500 --calibrate results_diag/transverse_gpu.csv
"""
import argparse
import csv
import os
import time

import numpy as np

from hdym_core import ym_initial, ym_step_substepped as _step


def hess(p):
    x, y, z = p.T
    H = np.empty((p.shape[0], 3, 3))
    H[:, 0, 0] = y * y + z * z
    H[:, 1, 1] = z * z + x * x
    H[:, 2, 2] = x * x + y * y
    H[:, 0, 1] = H[:, 1, 0] = 2 * x * y
    H[:, 0, 2] = H[:, 2, 0] = 2 * x * z
    H[:, 1, 2] = H[:, 2, 1] = 2 * y * z
    return H


def track_lowest(evals, evecs, prev_e0):
    """Follow ONE mode continuously: the eigenvector with the largest overlap
    with the previously tracked one (sign-fixed).  Returns (lam0, e0, lam_pair)
    with lam_pair = (tr H - lam0)/2 the mean of the other two.

    Why only one.  In a channel the two transverse modes are degenerate to
    O(2xy) and their eigenvectors rotate at frequency |x|, undersampled at any
    reasonable output rate, so labelling them is meaningless and produces
    jumps (a u^-2 spectral tail).  The lowest mode is well separated there
    (~ y^2+z^2 vs x^2) and slow in the interior, so it tracks cleanly, and the
    rotation rate |de0/dt| is invariant under rotations inside the pair."""
    N = evals.shape[0]
    if prev_e0 is None:
        j = np.zeros(N, int)
    else:
        j = np.abs(np.einsum('ni,nik->nk', prev_e0, evecs)).argmax(1)
    idx = np.arange(N)
    e0 = evecs[idx, :, j]
    if prev_e0 is not None:
        sgn = np.sign(np.einsum('ni,ni->n', prev_e0, e0))
        sgn[sgn == 0] = 1.0
        e0 = e0 * sgn[:, None]
    lam0 = evals[idx, j]
    lam_pair = 0.5 * (evals.sum(1) - lam0)
    return lam0, e0, lam_pair


def spectra(buf, tapers, dts):
    """Averaged multitaper spectral density along axis 0 of buf (nseg, ...).
    Returns (nfreq, ...).  S(om) = int <k(t)k(0)> e^{i om t} dt convention."""
    buf = buf - buf.mean(0)
    out = 0.0
    for w in tapers:
        F = np.fft.rfft(buf * w.reshape((-1,) + (1,) * (buf.ndim - 1)), axis=0)
        out = out + np.abs(F) ** 2 * (dts / np.sum(w * w))
    return out / len(tapers)


def _worker(args):
    """Integrate N trajectories for T after burn-in, return spectra summed over
    segments (not averaged), the per-segment depths, and the segment count."""
    (N, T, dt, sample, t_burn, seed, omega_dt_max, nseg_pow, taper, nw, ktapers, wid) = args
    import os
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    from scipy.signal.windows import blackmanharris, dpss
    rng = np.random.default_rng(seed)
    p, P = ym_initial(N, 1.0, rng)
    for _ in range(int(t_burn / dt)):
        p, P = _step(p, P, dt, omega_dt_max)
    nseg = 2 ** nseg_pow
    dts = dt * sample
    nsegs = max(1, int(T / (nseg * dts)))
    tapers = dpss(nseg, NW=nw, Kmax=ktapers) if taper == 'dpss' else blackmanharris(nseg)[None, :]
    acc = {}
    depth_all = []
    prev = None
    for s in range(nsegs):
        ent = np.empty((nseg, N, 6))
        eig = np.empty((nseg, N, 3))          # lam0 (tracked), lam_pair, lam_min (sorted)
        e0s = np.empty((nseg, N, 3))
        mx = np.zeros(N)
        for j in range(nseg):
            for _ in range(sample):
                p, P = _step(p, P, dt, omega_dt_max)
            H = hess(p)
            ent[j] = np.stack([H[:, 0, 0], H[:, 1, 1], H[:, 2, 2],
                               H[:, 0, 1], H[:, 0, 2], H[:, 1, 2]], 1)
            ev, vc = np.linalg.eigh(H)
            lam0, e0, lpair = track_lowest(ev, vc, prev)
            eig[j, :, 0], eig[j, :, 1], eig[j, :, 2] = lam0, lpair, ev[:, 0]
            e0s[j] = e0
            prev = e0
            mx = np.maximum(mx, np.abs(p).max(1))
        de = (e0s[2:] - e0s[:-2]) / (2 * dts)
        rot = np.sqrt(np.sum(de * de, -1))
        rot = np.concatenate([rot[:1], rot, rot[-1:]], 0)[:, :, None]
        S = spectra(ent, tapers, dts)
        Se = spectra(eig, tapers, dts)
        Sr = spectra(rot, tapers, dts)[:, :, 0]
        med = np.median(mx)
        sh = mx <= med
        upd = dict(diag=S[:, :, :3].mean((1, 2)), off=S[:, :, 3:].mean((1, 2)),
                   eig=Se.mean(1), rot=Sr.mean(1),
                   diag_shallow=S[:, sh, :3].mean((1, 2)), diag_deep=S[:, ~sh, :3].mean((1, 2)))
        for k, v in upd.items():
            acc[k] = acc.get(k, 0) + v
        depth_all.append(mx)
        if wid == 0:
            print(f"  worker 0: segment {s + 1}/{nsegs}", flush=True)
    return acc, np.concatenate(depth_all), nsegs, dts, nseg


def run(N, T, dt, sample, t_burn, seed, omega_dt_max, nseg_pow, taper, nw=16, ktapers=8,
        workers=1, log=print):
    """Parallel over trajectory chunks (multiprocessing, spawn-safe).  The
    per-trajectory numpy work is single-threaded, so this is a near-linear
    speed-up up to the physical core count."""
    import time as _time
    t0 = _time.time()
    workers = max(1, min(workers, N))
    chunks = [N // workers + (1 if i < N % workers else 0) for i in range(workers)]
    jobs = [(n, T, dt, sample, t_burn, seed + 1000 * i, omega_dt_max, nseg_pow, taper, nw, ktapers, i)
            for i, n in enumerate(chunks) if n > 0]
    if len(jobs) == 1:
        results = [_worker(jobs[0])]
    else:
        import multiprocessing as mp
        ctx = mp.get_context('spawn')
        with ctx.Pool(len(jobs)) as pool:
            results = pool.map(_worker, jobs)
    acc = {}
    depth_all = []
    wsum = 0
    for (a, dep, nsegs, dts, nseg), n in zip(results, [j[0] for j in jobs]):
        for k, v in a.items():
            acc[k] = acc.get(k, 0) + v * n          # segment-sum, weighted by trajectories
        depth_all.append(dep)
        wsum += n * nsegs
    for k in acc:
        acc[k] = acc[k] / wsum
    om = 2 * np.pi * np.fft.rfftfreq(nseg, dts)
    log(f"{len(jobs)} workers, {N} trajectories, {nsegs} segments each ({_time.time() - t0:.0f}s)")
    return om, acc, np.concatenate(depth_all)


def predict(om, s, energies):
    """lambda(E) = E^{3/4} s(2 E^{-1/4}) / 8, s interpolated in log-log."""
    u = 2.0 * np.asarray(energies) ** -0.25
    ok = (om > 0) & (s > 0)
    ls = np.interp(np.log(u), np.log(om[ok]), np.log(s[ok]))
    return u, np.asarray(energies) ** 0.75 * np.exp(ls) / 8.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--N', type=int, default=512)
    ap.add_argument('--T', type=float, default=3000.0, help='time per trajectory after burn-in')
    ap.add_argument('--dt', type=float, default=0.005)
    ap.add_argument('--sample', type=int, default=10, help='record every `sample` steps')
    ap.add_argument('--t-burn', type=float, default=500.0, help='YM units (= absolute at E = 1)')
    ap.add_argument('--nseg-pow', type=int, default=13, help='FFT segment length = 2**this samples')
    ap.add_argument('--omega-dt-max', type=float, default=0.1)
    ap.add_argument('--taper', choices=['dpss', 'bh'], default='dpss')
    ap.add_argument('--nw', type=float, default=16.0, help='DPSS time-bandwidth product')
    ap.add_argument('--ktapers', type=int, default=8, help='number of Slepian tapers (<= NW/2 here!)')
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--energies', type=float, nargs='*', default=list(np.logspace(-4, 0, 33)))
    ap.add_argument('--calibrate', default=None,
                    help='CSV with columns E, lam_mean -- use the --ghost diag scan, since the '
                         'formula is for the single-mode mechanism; one constant fitted on E <= --cal-emax')
    ap.add_argument('--cal-emax', type=float, default=0.05)
    ap.add_argument('--entry', choices=['diag', 'eig', 'eigmin'], default='diag',
                    help='spectrum for the prediction: diag = mean diagonal entry (fixed basis); '
                         'eig = continuously tracked lowest mode; eigmin = sorted lowest eigenvalue '
                         '(has kinks at crossings).  See docstring on which basis is physical.')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) // 2),
                    help='processes (default: half the logical CPUs = physical cores)')
    ap.add_argument('--out', default='results_spectrum')
    a = ap.parse_args()
    print(f"hdym_spectrum v2.1  taper={a.taper} NW={a.nw} K={a.ktapers} workers={a.workers}", flush=True)
    os.makedirs(a.out, exist_ok=True)

    om, acc, depth = run(a.N, a.T, a.dt, a.sample, a.t_burn, a.seed, a.omega_dt_max,
                         a.nseg_pow, a.taper, a.nw, a.ktapers, a.workers)
    np.savez(os.path.join(a.out, 'spectrum.npz'), om=om, S_diag=acc['diag'], S_off=acc['off'],
             S_eig=acc['eig'], S_rot=acc['rot'], S_diag_shallow=acc['diag_shallow'],
             S_diag_deep=acc['diag_deep'], depth=depth, taper=a.taper)
    print("channel depth per segment: median %.1f, p99 %.1f, max %.1f"
          % tuple(np.percentile(depth, [50, 99, 100])))

    curves = dict(diag=acc['diag'], eig=acc['eig'][:, 0], eigmin=acc['eig'][:, 2])
    s = curves[a.entry]
    ok = (om > 1) & (s > 0)
    lo, ls = np.log(om[ok]), np.log(s[ok])
    print(f"tail of s(u) [{a.entry}]: local log-log slope")
    for u0, u1 in [(2, 3), (3, 5), (5, 8), (8, 12), (12, 20), (20, 30), (30, 45)]:
        m = (lo >= np.log(u0)) & (lo <= np.log(u1))
        if m.sum() > 3:
            print(f"   u = {u0:>3}-{u1:<3}: {np.polyfit(lo[m], ls[m], 1)[0]:6.2f}")

    E = np.array(a.energies)
    preds = {k: predict(om, v, E)[1] for k, v in curves.items()}
    u = 2.0 * E ** -0.25
    cal = np.nan
    data = None
    if a.calibrate:
        with open(a.calibrate) as fh:
            rows = [r for r in csv.DictReader(fh)]
        data = np.array([(float(r['E']), float(r['lam_mean'])) for r in rows])
        m = data[:, 0] <= a.cal_emax
        if m.sum():
            _, lp = predict(om, s, data[m, 0])
            cal = np.exp(np.mean(np.log(data[m, 1] / lp)))
            print(f"calibration factor (data / prediction[{a.entry}], geometric mean over "
                  f"{m.sum()} points with E <= {a.cal_emax}): {cal:.3f}")
    with open(os.path.join(a.out, 'prediction.csv'), 'w', newline='') as fh:
        wr = csv.writer(fh)
        wr.writerow(['E', 'u', 'lam_diag', 'lam_eig', 'lam_eigmin', f'lam_{a.entry}_cal'])
        for i, e in enumerate(E):
            wr.writerow([e, u[i], preds['diag'][i], preds['eig'][i], preds['eigmin'][i],
                         preds[a.entry][i] * cal if np.isfinite(cal) else ''])
    print(f"{'E':>10} {'u':>6} {'diag':>10} {'eig':>10} {'eigmin':>10}" + ("   calibrated" if np.isfinite(cal) else ""))
    for i, e in enumerate(E):
        line = f"{e:10.3e} {u[i]:6.2f} {preds['diag'][i]:10.3e} {preds['eig'][i]:10.3e} {preds['eigmin'][i]:10.3e}"
        if np.isfinite(cal):
            line += f"   {preds[a.entry][i] * cal:.3e}"
        print(line)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
        ok = om > 0.3
        ax[0].loglog(om[ok], acc['diag'][ok], label='diag entries')
        ax[0].loglog(om[ok], acc['off'][ok], label='off-diag entries', alpha=.6)
        ax[0].loglog(om[ok], acc['eig'][ok, 0], '-', lw=.8, label='tracked lowest mode')
        ax[0].loglog(om[ok], acc['eig'][ok, 1], '-', lw=.8, label='mean of other two')
        ax[0].loglog(om[ok], acc['eig'][ok, 2], '-', lw=.6, label='sorted lowest (kinks)')
        ax[0].loglog(om[ok], acc['diag_deep'][ok], '--', label='diag, deep half', alpha=.6)
        ax[0].loglog(om[ok], acc['diag_shallow'][ok], ':', label='diag, shallow half', alpha=.6)
        ax[0].set_xlabel('u'); ax[0].set_ylabel('s(u)'); ax[0].legend(fontsize=7)
        ax[0].set_title(f'Hess V spectra, E = 1 ({a.taper})')
        ax[1].loglog(om[ok], acc['rot'][ok], lw=.8)
        ax[1].set_xlabel('u'); ax[1].set_ylabel('s(u)')
        ax[1].set_title('rotation rate |de0/dt| of the lowest mode (mixing)')
        for k in curves:
            ax[2].loglog(E, preds[k], label=f'from {k}', lw=.9)
        if np.isfinite(cal):
            ax[2].loglog(E, preds[a.entry] * cal, 'k--', label=f'{a.entry} x {cal:.2f}')
        if data is not None:
            ax[2].loglog(data[:, 0], data[:, 1], 'o', label='measured (calib. file)', ms=4)
        ax[2].set_xlabel('E'); ax[2].set_ylabel('lambda'); ax[2].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(a.out, 'spectrum.png'), dpi=130)
    except Exception as e:  # noqa
        print("plot skipped:", e)
    print("wrote", a.out)


if __name__ == '__main__':
    main()
