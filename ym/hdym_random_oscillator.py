"""
hdym_random_oscillator.py -- the paper's spectral formula (eq:spectral, paper/main.tex: the randomly modulated
oscillator of Arnold, Papanicolaou & Wihstutz 1986 and Jarzynski 1993) re-evaluated with a fresh measured spectrum,
on the diagonal ghost and on the c_j = 0 off-diagonal blocks (notes/NOTE_jeans_sea.md, 29 Oct 2026).

For w'' = -(Omega^2 + h(t)) w with a weak stationary drive h, lambda = S_h(2 Omega) / (8 Omega^2), S_h the two-sided
spectral density.  In scaled units (q = E^{1/4} q^, tau = E^{1/4} t, E^ = 1) Omega = E^{-1/4} and
lambda_t = E^{1/4} lambda_tau, i.e. eq:spectral, lambda = E^{3/4} s^(2 E^{-1/4}) / 8.  Two sectors are exactly scalar:
  - the diagonal ghost of `hdym_gpu.py --ghost diag` (h = H_xx = y^2 + z^2); rates in data/transverse_diag_percomp
    {,_hi}.csv.  The paper predicts the ratio measured / (eq:spectral with the full spectrum) to be the resonant share
    2/7 (the y- and z-channel power is detuned) and measures 0.276-0.290 at E = 0.03 ... 0.008;
  - the c_j = 0 blocks, w_pm'' = -(1 + z^2 -+ xy) w_pm, with no detuning (paper, "The mechanism can be isolated"), so
    S_h = S_{z^2} + S_{xy} (the cross term is odd under y -> -y and averages out); rates lambda_inc of the three blocks
    in data/jeans_sea.txt (hdym_jeans_sea.py, same ensemble recipe); the closed-form law 0.0210 E^{5/4}, whose excess
    the paper attributes to the knee of the xy spectrum at u <= 3.6.
A power law S ~ omega^(-p) gives lambda_t ~ E^{(p+3)/4}: p = 6 gives 9/4, p = 2 gives 5/4.
The "unlifted" rows are a coarse, chunk-level version of the paper's resonant share (sea or x-led chunks only); a
20-unit chunk can hold several channel visits, so it is not the paper's accounting.

The spectrum: ghost-free YM at E^ = 1, hdym's GPU initial recipe (q, p standard normal, scaled; burn-in 300), DOP853
rtol 1e-9, chunks of 20.48 sampled at 1024 points (Kaiser beta 14 window, mean removed), each chunk tagged by its
maximum depth max|q^| and by which coordinate went deepest.  Spectra are per unit total time (sums over a subset of
chunks divided by the number of all chunks), for: all chunks; sea (max depth < 3); and, for H_xx, "unlifted" chunks
(sea, or the x-channel), i.e. without the y- and z-channels, where H_xx is large and lifts the ghost component's own
frequency (the detuned piece; paper, after eq:seven).

    python ym/hdym_random_oscillator.py --workers 22 > ym/data/random_oscillator.txt
The committed log and npz are from --workers 2 (48 trajectories, 594 s, on Claude's side); the log was then regenerated
with --from-npz after a column fix in part 3 (the spectra are those of the integration).  Added 2026-10-29.
"""
import argparse, csv, os, re, sys, time
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))


def rhs(t, Y):
    Y = Y.reshape(-1, 6); x, y, z = Y[:, 0], Y[:, 1], Y[:, 2]
    return np.concatenate([Y[:, 3:6], -np.stack([x * (y * y + z * z), y * (x * x + z * z), z * (x * x + y * y)], 1)], 1).ravel()


ENT = ['H_xx = y^2+z^2', 'H_xy = 2xy', 'z^2', 'xy']


def work(args):
    wk, n, tau, chunk, seed, rtol = args
    rng = np.random.default_rng(seed + 1000 * wk)
    q = rng.standard_normal((n, 3)); p = rng.standard_normal((n, 3))
    lam = (1 / (0.5 * (p * p).sum(1) + 0.5 * ((q[:, 0] * q[:, 1]) ** 2 + (q[:, 1] * q[:, 2]) ** 2 + (q[:, 0] * q[:, 2]) ** 2))) ** 0.25
    Y = np.concatenate([q * lam[:, None], p * (lam ** 2)[:, None]], 1)
    E0 = 0.5 * (Y[:, 3:] ** 2).sum(1) + 0.5 * ((Y[:, 0] * Y[:, 1]) ** 2 + (Y[:, 1] * Y[:, 2]) ** 2 + (Y[:, 0] * Y[:, 2]) ** 2)
    Y = solve_ivp(rhs, (0, 300), Y.ravel(), method='DOP853', rtol=rtol, atol=rtol * 1e-2).y[:, -1]
    nseg = 1024; dt = chunk / nseg; te = np.arange(nseg + 1) * dt; te[-1] = chunk
    win = np.kaiser(nseg, 14.0); wn = (win ** 2).sum(); nf = nseg // 2 + 1
    S = np.zeros((4, len(ENT), nf)); cnt = np.zeros(4); t = 0.0     # rows: all, sea, unlifted (sea or x-led), lifted
    while t < tau - 1e-9:
        s = solve_ivp(rhs, (t, t + chunk), Y, method='DOP853', rtol=rtol, atol=rtol * 1e-2, t_eval=t + te)
        Y = s.y[:, -1]; Q = s.y[:, :nseg].reshape(n, 6, nseg); x, y, z = Q[:, 0], Q[:, 1], Q[:, 2]
        amp = np.abs(Q[:, :3]).max(axis=2); dep = amp.max(1); lead = amp.argmax(1)
        sea = dep < 3; unl = sea | (lead == 0); lif = ~unl
        for k, h in enumerate([y * y + z * z, 2 * x * y, z * z, x * y]):
            hh = h - h.mean(1, keepdims=True)
            P = np.abs(np.fft.rfft(hh * win, axis=1)) ** 2 * dt / wn
            for r, m in enumerate((np.ones(n, bool), sea, unl, lif)):
                S[r, k] += P[m].sum(0)
        for r, m in enumerate((np.ones(n, bool), sea, unl, lif)):
            cnt[r] += m.sum()
        t += chunk
    Yf = Y.reshape(n, 6)
    Ef = 0.5 * (Yf[:, 3:] ** 2).sum(1) + 0.5 * ((Yf[:, 0] * Yf[:, 1]) ** 2 + (Yf[:, 1] * Yf[:, 2]) ** 2 + (Yf[:, 0] * Yf[:, 2]) ** 2)
    return S, cnt, float(np.abs(Ef - E0).max())


def read_diag():
    out = {}
    for f in ('transverse_diag_percomp_hi.csv', 'transverse_diag_percomp.csv'):
        with open(os.path.join(HERE, 'data', f)) as fh:
            for row in csv.DictReader(fh):
                out[float(row['E'])] = (float(row['lam_mean']), float(row['lam_se']))
    return dict(sorted(out.items(), reverse=True))


def read_blocks():
    out = {}
    with open(os.path.join(HERE, 'data', 'jeans_sea.txt'), encoding='utf-8') as fh:
        for ln in fh:
            m = re.match(r'\s*([\d.]+)\s+(b12|b13|b23)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)', ln)
            if m:
                E = float(m.group(1)); out.setdefault(E, []).append((float(m.group(3)), float(m.group(4)), float(m.group(5))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=22)
    ap.add_argument('--per', type=int, default=24)
    ap.add_argument('--tau', type=float, default=20000.0)
    ap.add_argument('--chunk', type=float, default=20.48)
    ap.add_argument('--seed', type=int, default=27)
    ap.add_argument('--rtol', type=float, default=1e-9)
    ap.add_argument('--out', default=os.path.join(HERE, 'data', 'random_oscillator.npz'))
    ap.add_argument('--from-npz', action='store_true', help='reuse the spectra in --out; no integration')
    a = ap.parse_args()
    T0 = time.time()
    print(f'# START hdym_random_oscillator workers={a.workers} per={a.per} tau={a.tau:g} chunk={a.chunk} seed={a.seed} rtol={a.rtol:g}'
          f' from_npz={a.from_npz}  python {sys.version.split()[0]} numpy {np.__version__}', flush=True)
    dt = a.chunk / 1024; om = 2 * np.pi * np.fft.rfftfreq(1024, dt)
    if a.from_npz:
        Z = np.load(a.out); S, cnt, dE = Z['S'], Z['cnt'], float(Z['dE'])
        print(f'# spectra from {os.path.relpath(a.out)} (written by an earlier run of this script with the parameters above)')
    else:
        from multiprocessing import Pool
        with Pool(min(a.workers, os.cpu_count() or 1)) as pool:
            R = pool.map(work, [(k, a.per, a.tau, a.chunk, a.seed, a.rtol) for k in range(a.workers)])
        S = sum(r[0] for r in R); cnt = sum(r[1] for r in R); dE = max(r[2] for r in R)
        np.savez(a.out, S=S, cnt=cnt, omega=om, dE=dE, entries=np.array(ENT), rows=np.array(['all', 'sea', 'unlifted', 'lifted']))
    tot = cnt[0]; Sp = S / tot                                       # per unit total time
    print(f'# {a.workers * a.per} trajectories, tau = {a.tau:g} each; time fractions: sea {cnt[1] / tot:.3f}, unlifted {cnt[2] / tot:.3f}, '
          f'lifted {cnt[3] / tot:.3f}; max |dE|/E {dE:.1e}', flush=True)
    print('# 1. spectra S(omega) per unit total time, all chunks (Kaiser 14)')
    print(f"{'omega':>6} " + ' '.join(f'{e:>15}' for e in ENT) + f" {'H_xx unlifted':>15} {'H_xx sea':>12}")
    for w in (2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 25, 30):
        i = np.argmin(np.abs(om - w))
        print(f'{om[i]:6.2f} ' + ' '.join(f'{Sp[0, k, i]:15.3e}' for k in range(len(ENT))) + f' {Sp[2, 0, i]:15.3e} {Sp[1, 0, i]:12.3e}')
    print('# power-law fits of ln S on ln omega (slope -p; the law exponent is (p + 3)/4)')
    for lab, r, k in (('H_xx all', 0, 0), ('H_xx unlifted', 2, 0), ('H_xy all', 0, 1), ('z^2 all', 0, 2), ('xy all', 0, 3)):
        for lo, hi in ((4, 8), (6, 14), (8, 20)):
            sel = (om >= lo) & (om <= hi)
            pp = np.polyfit(np.log(om[sel]), np.log(Sp[r, k, sel]), 1)
            rms = np.sqrt(np.mean((np.polyval(pp, np.log(om[sel])) - np.log(Sp[r, k, sel])) ** 2))
            print(f'  {lab:14s} omega {lo:2d}-{hi:2d}: p = {-pp[0]:.2f}, rms {rms:.3f}, law exponent {(-pp[0] + 3) / 4:.3f}')
    L = lambda Sx, w: float(np.exp(np.interp(w, om[1:], np.log(Sx[1:]))))
    print('# 2. diagonal ghost (scalar h = H_xx): lambda_pred = E^{1/4} S(2 Omega)/(8 Omega^2) against hdym_gpu --ghost diag')
    print(f"{'E':>6} {'2Omega':>7} {'measured':>11} {'pred all':>10} {'ratio':>6} {'pred unl.':>10} {'ratio':>6} {'0.0120E^9/4':>12}")
    for E, (m, se) in read_diag().items():
        O = E ** -0.25; pa = E ** 0.25 * L(Sp[0, 0], 2 * O) / (8 * O * O); pu = E ** 0.25 * L(Sp[2, 0], 2 * O) / (8 * O * O)
        print(f'{E:6.3f} {2 * O:7.3f} {m:11.3e} {pa:10.3e} {m / pa:6.3f} {pu:10.3e} {m / pu:6.3f} {0.0120 * E ** 2.25:12.3e}')
    print('# 3. c_j = 0 blocks (scalar h = z^2 -+ xy, S_h = S_z2 + S_xy): against data/jeans_sea.txt lambda_inc (mean of 3 blocks)')
    print(f"{'E':>6} {'2Omega':>7} {'measured':>11} {'+-':>8} {'pred':>10} {'ratio':>6} {'0.0210E^5/4':>12}")
    Sh = Sp[0, 2] + Sp[0, 3]
    for E, rows in sorted(read_blocks().items(), reverse=True):
        m = np.mean([r[0] for r in rows]); se = np.mean([r[1] for r in rows])          # the blocks share trajectories: no 1/sqrt(3)
        O = E ** -0.25; pr = E ** 0.25 * L(Sh, 2 * O) / (8 * O * O)
        print(f'{E:6.3f} {2 * O:7.3f} {m:11.3e} {se:8.1e} {pr:10.3e} {m / pr:6.3f} {0.0210 * E ** 1.25:12.3e}')
    print(f'# END hdym_random_oscillator status=OK wall={time.time() - T0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
