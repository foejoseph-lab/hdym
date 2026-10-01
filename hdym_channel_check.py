#!/usr/bin/env python3
"""
hdym_channel_check.py -- four numerical checks of the channel-passage
derivation of the u^-6 tail (pure Yang-Mills, E = 1, no ghost).

  Step 1  In a channel along axis i the transverse action J = E_perp/|x_i|
          is adiabatically invariant and x_i(t) is a parabola (linear
          effective potential J|x|), turning at depth D = E_par/J ~ 1/J.
  Step 2  K_ii = (sum of the other two squares) is a chirp: amplitude ~ J/|x|,
          instantaneous frequency 2|x(t)|.
  Step 5  The transverse action at channel entry has density P(J) ~ J near 0
          (two transverse degrees of freedom), equivalently
          P(depth > D) ~ D^-2.  This is the step that turns u^-5 into u^-6.

Panels in channel_check.png
  (a) the deepest visit found: x_i(t) with a parabola fit, J(t) on a twin axis
  (b) spectrogram of K_ii over that visit with 2|x_i(t)| overlaid
  (c) peak depth vs 1/J_entry for all visits (expect slope 1)
  (d) CDF of J_entry (expect slope 2 near 0) and exceedance of depth
      (expect slope -2)

Usage:  python hdym_channel_check.py --N 256 --T 2000 --A-in 3
"""
import argparse
import time

import numpy as np

from hdym_core import ym_initial, ym_step_substepped as step


def transverse_action(p, P, i):
    """J = E_perp / |x_i| for the channel along axis i (vectorised over rows)."""
    x = p[:, i]
    j = [k for k in range(3) if k != i]
    Eperp = 0.5 * (P[:, j[0]] ** 2 + P[:, j[1]] ** 2) + 0.5 * x * x * (p[:, j[0]] ** 2 + p[:, j[1]] ** 2)
    return Eperp / np.abs(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=256)
    ap.add_argument('--T', type=float, default=2000.0)
    ap.add_argument('--dt', type=float, default=0.005)
    ap.add_argument('--t-burn', type=float, default=200.0)
    ap.add_argument('--A-in', type=float, default=3.0, help='channel entry level in scaled units')
    ap.add_argument('--seed', type=int, default=5)
    ap.add_argument('--out', default='channel_check.png')
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    N, dt = a.N, a.dt
    p, P = ym_initial(N, 1.0, rng)
    for _ in range(int(a.t_burn / dt)):
        p, P = step(p, P, dt, 0.1)

    # ---- visit bookkeeping over the ensemble
    inside = np.zeros(N, bool)
    axis = np.zeros(N, int)
    J_entry = np.zeros(N)
    peak = np.zeros(N)
    J_at_peak = np.zeros(N)
    t_entry = np.zeros(N)
    visits = []                      # (J_entry, peak, J_at_peak, duration, traj)
    # full record of one trajectory (the one that ends up with the deepest visit
    # is not known in advance, so record everything for a subset and pick later)
    nrec = min(N, 64)
    nsteps = int(a.T / dt)
    rec = np.empty((nsteps, nrec, 6), np.float32)
    t0 = time.time()
    for k in range(nsteps):
        p, P = step(p, P, dt, 0.1)
        rec[k, :, :3] = p[:nrec]
        rec[k, :, 3:] = P[:nrec]
        amax = np.abs(p).max(1)
        imax = np.abs(p).argmax(1)
        # entries
        e = (~inside) & (amax > a.A_in)
        if e.any():
            idx = np.nonzero(e)[0]
            inside[idx] = True
            axis[idx] = imax[idx]
            J_entry[idx] = transverse_action(p[idx], P[idx], 0) * 0  # placeholder shape
            for i in np.unique(imax[idx]):
                sel = idx[imax[idx] == i]
                J_entry[sel] = transverse_action(p[sel], P[sel], i)
            peak[idx] = amax[idx]
            J_at_peak[idx] = J_entry[idx]
            t_entry[idx] = k * dt
        # update peaks
        ins = np.nonzero(inside)[0]
        if ins.size:
            cur = np.abs(p[ins, axis[ins]])
            better = cur > peak[ins]
            if better.any():
                b = ins[better]
                peak[b] = cur[better]
                for i in np.unique(axis[b]):
                    sel = b[axis[b] == i]
                    J_at_peak[sel] = transverse_action(p[sel], P[sel], i)
            # exits
            x = np.abs(p[ins, axis[ins]])
            out = x < a.A_in
            for n in ins[out]:
                visits.append((J_entry[n], peak[n], J_at_peak[n], k * dt - t_entry[n], n, t_entry[n]))
                inside[n] = False
    print(f"{len(visits)} visits in {N} trajectories x T={a.T} ({time.time() - t0:.0f}s)")
    V = np.array([v[:4] for v in visits])
    Je, Dpk, Jpk, dur = V.T

    # ---- Step 1: adiabatic invariance of J
    ratio = Jpk / Je
    print(f"Step 1  J_peak/J_entry: median {np.median(ratio):.3f}, "
          f"IQR {np.percentile(ratio, 25):.3f}-{np.percentile(ratio, 75):.3f} (expect ~1)")
    # ---- Step 1: depth ~ 1/J
    m = Dpk > 2 * a.A_in
    sl, ic = np.polyfit(np.log(1 / Je[m]), np.log(Dpk[m]), 1)
    print(f"Step 1  log(depth) vs log(1/J_entry) for depth > {2 * a.A_in}: slope {sl:.2f} "
          f"(expect 1), depth ~ {np.exp(ic):.2f}/J")
    # ---- Step 5: P(J) ~ J  <=>  CDF ~ J^2 ; exceedance of depth ~ D^-2
    Js = np.sort(Je)
    cdf = np.arange(1, len(Js) + 1) / len(Js)
    lo = (Js > np.percentile(Js, 1)) & (Js < np.percentile(Js, 20))
    s_cdf = np.polyfit(np.log(Js[lo]), np.log(cdf[lo]), 1)[0]
    Ds = np.sort(Dpk)[::-1]
    exc = np.arange(1, len(Ds) + 1) / len(Ds)
    hi = (Ds > np.percentile(Dpk, 80)) & (Ds < np.percentile(Dpk, 99.5))
    s_exc = np.polyfit(np.log(Ds[hi]), np.log(exc[hi]), 1)[0]
    print(f"Step 5  CDF of J_entry, low-J slope {s_cdf:.2f} (expect 2 => P(J) ~ J, tail u^-6; "
          f"1 => P(J) ~ const, tail u^-5)")
    print(f"Step 5  exceedance of peak depth, slope {s_exc:.2f} (expect -2)")

    # ---- Steps 1-2 on the deepest recorded visit
    deep = [v for v in visits if v[4] < nrec]
    v = max(deep, key=lambda v: v[1])
    n, te, du = v[4], v[5], v[3]
    k0, k1 = int(te / dt), int((te + du) / dt)
    seg = rec[k0:k1, n].astype(float)
    i = int(np.abs(seg[:, :3]).mean(0).argmax())
    j = [k for k in range(3) if k != i]
    t = np.arange(k1 - k0) * dt
    x = seg[:, i]
    J = transverse_action(seg[:, :3], seg[:, 3:], i)
    Kii = seg[:, j[0]] ** 2 + seg[:, j[1]] ** 2
    c = np.polyfit(t, x, 2)
    print(f"Steps 1-2 on deepest visit: traj {n}, axis {i}, depth {v[1]:.1f}, duration {du:.1f}, "
          f"J_entry {v[0]:.4f}; parabola fit x = {c[0]:.3f} t^2 + ..., "
          f"-> J_fit = {2 * abs(c[0]):.4f} (expect ~J_entry); rms residual {np.std(x - np.polyval(c, t)):.3f}")

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.signal import spectrogram
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
    a0 = ax[0, 0]
    a0.plot(t, x, lw=.8, label='x_i(t)')
    a0.plot(t, np.polyval(c, t), 'k--', lw=.8, label='parabola fit')
    a0.set_xlabel('t'); a0.set_ylabel('x_i'); a0.legend(loc='upper left', fontsize=8)
    a0b = a0.twinx(); a0b.plot(t, J, 'r', lw=.6, alpha=.7); a0b.set_ylabel('J (red)', color='r')
    a0b.set_ylim(0, 3 * np.median(J))
    a0.set_title('(a) Step 1: parabola and adiabatic invariant')
    nper = 256
    f, tt, S = spectrogram(Kii - Kii.mean(), fs=1 / dt, nperseg=nper, noverlap=nper * 3 // 4, scaling='density')
    a1 = ax[0, 1]
    a1.pcolormesh(tt, 2 * np.pi * f, np.log10(S + 1e-30), shading='auto', vmin=np.log10(S.max()) - 8)
    a1.plot(t, 2 * np.abs(x), 'w--', lw=.8, label='2|x_i(t)|')
    a1.set_ylim(0, 2.5 * np.abs(x).max()); a1.set_xlabel('t'); a1.set_ylabel('frequency (rad)')
    a1.legend(fontsize=8); a1.set_title('(b) Step 2: K_ii is a chirp at 2|x|')
    a2 = ax[1, 0]
    a2.loglog(1 / Je, Dpk, '.', ms=2, alpha=.4)
    xx = np.logspace(np.log10((1 / Je).min()), np.log10((1 / Je).max()), 2)
    a2.loglog(xx, np.exp(ic) * xx ** sl, 'k--', label=f'slope {sl:.2f}')
    a2.set_xlabel('1 / J_entry'); a2.set_ylabel('peak depth'); a2.legend(fontsize=8)
    a2.set_title('(c) Step 1: depth ~ 1/J')
    a3 = ax[1, 1]
    a3.loglog(Js, cdf, lw=.8, label=f'CDF of J_entry (low-J slope {s_cdf:.2f})')
    a3.loglog(Ds, exc, lw=.8, label=f'exceedance of depth (slope {s_exc:.2f})')
    a3.set_xlabel('J_entry  |  depth'); a3.set_ylabel('probability'); a3.legend(fontsize=8)
    a3.set_title('(d) Step 5: P(J) ~ J  <=>  P(depth > D) ~ D^-2')
    fig.tight_layout(); fig.savefig(a.out, dpi=130)
    print("wrote", a.out)


if __name__ == '__main__':
    main()
