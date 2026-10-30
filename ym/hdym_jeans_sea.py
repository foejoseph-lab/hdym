"""
hdym_jeans_sea.py -- where the naive ghost grows (sea or channels) and whether it tracks the sea-only spectrum of the
drive: the "if resumed" items of the parked Boltzmann-Jeans idea (27 Oct 2026, sandbox notes, not committed).  Moved
into ym/ on 29 Oct 2026 with its run (data/jeans_sea.txt, 22 workers, 5125 s, on the author's machine); the result and
the random-oscillator reading are in notes/NOTE_random_oscillator.md (hdym_random_oscillator.py).

Model (as in the parked note): SU(2) colour-diagonal sector, V = (x^2 y^2 + y^2 z^2 + z^2 x^2)/2, naive ghost
w'' = -(1 + Hess V) w (the vector ghost, 3x3) and, separately, the three off-diagonal 2x2 blocks
u'' = -(1 + [[q_k^2, q_i q_j], [q_i q_j, q_k^2]]) u (hdym's c_j = 0 Hill blocks; the sign of the off-diagonal entry
is a change of basis and does not matter here).  Every sector has its OWN tangent (the su2_blocks lesson: one
tangent through several sectors measures their maximum).
Scaling: q = E^(1/4) q^(tau), tau = E^(1/4) t, so ONE pure-YM ensemble at E^ = 1 serves every E, the ghost seeing
w_tautau = -(Omega^2 + H^(tau)) w with Omega = E^(-1/4); physical rate lambda_t = E^(1/4) lambda_tau.

Items measured, on the same trajectories:
 1. Spectrum of the drive restricted to sea time.  Each chunk of --chunk tau (1024 samples, d tau = chunk/1024) is
    one Hann-windowed segment; its power spectrum S(omega) (two-sided density, mean over the matrix entries of the
    sector's drive, mean removed) is added to the class of its maximum depth, depth = max_i |q^_i| (the channel's
    transverse frequency at E^ = 1): classes [0,2), [2,3), [3,5), [5,inf).  "Sea" = max depth < 3 over the chunk.
    Reported: S_sea and S_all at omega = 2 .. 14, and fits of ln S_sea on omega (exponential) and on ln omega (power)
    over --fit-band.
 2. Sea-time-weighted rate.  The log-growth of every tangent is accumulated per chunk and booked to the chunk's depth
    class: lambda_sea = (growth in sea chunks)/(sea time), lambda_chan likewise, f_sea = sea time fraction.
 3. Tracking: R = lambda / (S(2 Omega) / (8 Omega^2)) per sector and E, with S = S_sea and S = S_all.  The weak-
    modulation formula predicts R = const (prefactor not derived here for the matrix case); "lambda tracks the sea
    spectrum" means R_sea is flat in E where R_all is not.  If it is flat, S_sea predicts lambda below the energies
    reachable by direct integration ([est] lines at E = 0.03, 0.01, 0.003).
Estimators per sector and E: lambda_T = L(T)/T and lambda_inc (second-half increment), SE over trajectories, and
lambda*T so a reader can see which E are trustworthy (lambda*T >= 20 rule, ym_run's convergence rule).

Initial data: hdym's GPU recipe (ym/hdym_rung1_gpu.py initial_state: q, p standard normal, scaled to the energy),
burn-in --tau-burn = 300 (the GPU's t_burn = 300 E^(-1/4) in scaled time), so the E = 0.1 block value can be checked
against hdym's 1.543e-3 (GPU one block, offset-free slope, t = 20,000 i.e. tau = 11,247) and E = 0.3 vector ghost
against 0.016506 (hdym ym/data/transverse_E0.3.npz).  Integrator: DOP853, rtol --rtol, vectorised over each
worker's trajectories (the step follows the deepest channel among them; that is the cost).  Sampling d tau = 0.02
(Nyquist 157): deep channels (max|q^| up to a few hundred) alias into the "all" spectra only; sea chunks have
max|q^| < 3 and are clean.  Windows: Kaiser beta = 14 for the numbers and fits, Hann alongside as a leakage check.

    python -B -u ym/hdym_jeans_sea.py --workers 22 --per 24 --tau 20000 --out jeans_out > ym/data/jeans_sea.txt 2>&1
    python -B -u ym/hdym_jeans_sea.py ... --resume >> ym/data/jeans_sea.txt 2>&1     (per-worker checkpoint)
    python -B -u ym/hdym_jeans_sea.py --combine-only --out jeans_out > table.txt 2>&1
The per-worker npz (jeans_out/) are working files and are not committed; the log is.  The run of data/jeans_sea.txt was
made from C:\\temp\\jeans as `python -B -u jeans_sea.py --workers 22 --per 24 --tau 20000 --out jeans_out`.
Written 29 Oct 2026.
"""
import argparse, glob, os, sys, time
import numpy as np
from scipy.integrate import solve_ivp

NCLS = 4
CLS_EDGES = (2.0, 3.0, 5.0)                      # depth classes [0,2) [2,3) [3,5) [5,inf)
SECT = ['vec', 'b12', 'b13', 'b23']
BLOCKS = [((0, 1), 2), ((0, 2), 1), ((1, 2), 0)]


def depth_class(d):
    return int(np.searchsorted(CLS_EDGES, d, side='right'))


def make_rhs(n, Om2):
    K = Om2.size

    def rhs(t, Y):
        Y = Y.reshape(n, -1)
        q, p = Y[:, 0:3], Y[:, 3:6]
        x, y, z = q.T
        g = np.stack([x * (y * y + z * z), y * (x * x + z * z), z * (x * x + y * y)], 1)
        Hd = np.empty((n, 3, 3))
        Hd[:, 0, 0] = y * y + z * z; Hd[:, 1, 1] = x * x + z * z; Hd[:, 2, 2] = x * x + y * y
        Hd[:, 0, 1] = Hd[:, 1, 0] = 2 * x * y; Hd[:, 0, 2] = Hd[:, 2, 0] = 2 * x * z; Hd[:, 1, 2] = Hd[:, 2, 1] = 2 * y * z
        G = Y[:, 6:].reshape(n, K, 18)
        out = np.empty_like(G)
        w, wd = G[:, :, 0:3], G[:, :, 3:6]
        out[:, :, 0:3] = wd
        out[:, :, 3:6] = -(Om2[None, :, None] * w + np.einsum('nij,nkj->nki', Hd, w))
        for b, ((i, j), k) in enumerate(BLOCKS):
            o = 6 + 4 * b
            u, ud = G[:, :, o:o + 2], G[:, :, o + 2:o + 4]
            d = (q[:, k] ** 2)[:, None]; c = (q[:, i] * q[:, j])[:, None]
            out[:, :, o:o + 2] = ud
            out[:, :, o + 2] = -(Om2[None, :] * u[:, :, 0] + d * u[:, :, 0] + c * u[:, :, 1])
            out[:, :, o + 3] = -(Om2[None, :] * u[:, :, 1] + c * u[:, :, 0] + d * u[:, :, 1])
        return np.concatenate([p, -g, out.reshape(n, -1)], 1).ravel()
    return rhs


SLICES = [slice(0, 6)] + [slice(6 + 4 * b, 10 + 4 * b) for b in range(3)]


def worker(args):
    wk, a = args
    a = argparse.Namespace(**a)
    Es = np.array(a.E); Om2 = Es ** -0.5; K = Es.size; n = a.per
    ck = os.path.join(a.out, f'w{wk:02d}.npz')
    rng = np.random.default_rng(a.seed + 1000 * wk)
    nseg = 1024; dtau = a.chunk / nseg
    wins = [np.kaiser(nseg, 14.0), np.hanning(nseg)]          # Kaiser beta 14 (main: ~110 dB sidelobes), Hann (check)
    wnorms = [(w_ ** 2).sum() for w_ in wins]
    nf = nseg // 2 + 1
    rhs = make_rhs(n, Om2)
    if a.resume and os.path.exists(ck):
        c = np.load(ck)
        Y = c['Y']; tau = float(c['tau']); L = c['L']; book = c['book']; btime = c['btime']
        Ssum = c['Ssum']; Scnt = c['Scnt']; Lh = list(c['Lh']); th = list(c['th']); wall0 = float(c['wall']); E0 = c['E0']
        dmax_all = float(c['dmax'])
    else:
        q = rng.standard_normal((n, 3)); p = rng.standard_normal((n, 3))
        lam = (1.0 / (0.5 * (p * p).sum(1) + 0.5 * ((q[:, 0] * q[:, 1]) ** 2 + (q[:, 1] * q[:, 2]) ** 2 + (q[:, 0] * q[:, 2]) ** 2))) ** 0.25
        q *= lam[:, None]; p *= (lam ** 2)[:, None]
        Y0 = np.concatenate([q, p], 1)
        Y0 = solve_ivp(lambda t, y: make_rhs(n, np.zeros(0))(t, y), (0, a.tau_burn), Y0.ravel(), method='DOP853',
                       rtol=a.rtol, atol=a.rtol * 1e-2).y[:, -1].reshape(n, 6)
        T = rng.standard_normal((n, K, 18))
        for s in SLICES:
            T[:, :, s] /= np.linalg.norm(T[:, :, s], axis=2, keepdims=True)
        Y = np.concatenate([Y0, T.reshape(n, -1)], 1)
        tau = 0.0; L = np.zeros((n, K, 4)); book = np.zeros((NCLS, K, 4)); btime = np.zeros(NCLS)
        Ssum = np.zeros((NCLS, 4, nf)); Scnt = np.zeros(NCLS); Lh = [L.copy()]; th = [0.0]; wall0 = 0.0
        E0 = 0.5 * (Y[:, 3:6] ** 2).sum(1) + 0.5 * ((Y[:, 0] * Y[:, 1]) ** 2 + (Y[:, 1] * Y[:, 2]) ** 2 + (Y[:, 0] * Y[:, 2]) ** 2)
        dmax_all = 0.0
    t0 = time.time(); last_save = time.time()
    teval = np.arange(nseg + 1) * dtau          # nseg samples for the spectrum, plus the chunk's end
    teval[-1] = a.chunk
    while tau < a.tau - 1e-9:
        sol = solve_ivp(rhs, (tau, tau + a.chunk), Y.ravel(), method='DOP853', rtol=a.rtol, atol=a.rtol * 1e-2,
                        t_eval=tau + teval)
        if sol.status != 0:
            raise RuntimeError(f'worker {wk}: solve_ivp failed at tau = {tau}: {sol.message}')
        Yend = sol.y[:, -1].reshape(n, -1).copy()
        S = sol.y[:, :nseg].reshape(n, -1, nseg)
        qs = S[:, 0:3, :]                                   # (n, 3, nseg)
        dep = np.abs(qs).max(axis=(1, 2)); dmax_all = max(dmax_all, float(dep.max()))
        cls = np.array([depth_class(d) for d in dep])
        x, y, z = qs[:, 0], qs[:, 1], qs[:, 2]
        hv = [y * y + z * z, x * x + z * z, x * x + y * y, 2 * x * y, 2 * x * z, 2 * y * z]   # 3 diag, 3 off (x2 weight)
        wv = np.array([1, 1, 1, 2, 2, 2]) / 9.0
        hb = [z * z, x * y, y * y, x * z, x * x, y * z]                                       # per block: diag, off
        wb = np.array([2, 2, 2, 2, 2, 2]) / 12.0
        for kind, (hs, ws) in enumerate(((hv, wv), (hb, wb))):
            for iw, (win, wnorm) in enumerate(zip(wins, wnorms)):
                P = np.zeros((n, nf))
                for h_, w_ in zip(hs, ws):
                    hh = h_ - h_.mean(axis=1, keepdims=True)
                    P += w_ * np.abs(np.fft.rfft(hh * win[None, :], axis=1)) ** 2 * dtau / wnorm
                for c_ in range(NCLS):
                    m = cls == c_
                    if m.any():
                        Ssum[c_, kind + 2 * iw] += P[m].sum(0)
        for c_ in range(NCLS):
            Scnt[c_] += np.sum(cls == c_)
        G = Yend[:, 6:].reshape(n, K, 18)
        inc = np.empty((n, K, 4))
        for si, s in enumerate(SLICES):
            nr = np.linalg.norm(G[:, :, s], axis=2)
            inc[:, :, si] = np.log(nr); G[:, :, s] /= nr[:, :, None]
        Yend[:, 6:] = G.reshape(n, -1)
        L += inc
        for c_ in range(NCLS):
            m = cls == c_
            if m.any():
                book[c_] += inc[m].sum(0); btime[c_] += m.sum() * a.chunk
        tau += a.chunk; Y = Yend
        Lh.append(L.copy()); th.append(tau)
        if time.time() - last_save > a.save_every or tau >= a.tau - 1e-9:
            tmp = ck + '.tmp.npz'
            np.savez(tmp, Y=Y, tau=tau, L=L, book=book, btime=btime, Ssum=Ssum, Scnt=Scnt, Lh=np.array(Lh), th=np.array(th),
                     wall=wall0 + time.time() - t0, E0=E0, dmax=dmax_all, E=Es, per=n, chunk=a.chunk, nseg=nseg)
            os.replace(tmp, ck); last_save = time.time()
    Ef = 0.5 * (Y[:, 3:6] ** 2).sum(1) + 0.5 * ((Y[:, 0] * Y[:, 1]) ** 2 + (Y[:, 1] * Y[:, 2]) ** 2 + (Y[:, 0] * Y[:, 2]) ** 2)
    return wk, wall0 + time.time() - t0, float(np.abs(Ef - E0).max()), dmax_all


def combine(a):
    Z = [np.load(p) for p in sorted(glob.glob(os.path.join(a.out, 'w*.npz')))]
    Z = [z for z in Z if not str(z.fid.name if hasattr(z, 'fid') else '').endswith('.tmp.npz')]
    if not Z:
        print('# no worker files'); return 'FAILED(no data)'
    Es = Z[0]['E']; K = Es.size; Om = Es ** -0.25
    nr = min(len(z['th']) for z in Z); th = Z[0]['th'][:nr]
    LH = np.concatenate([z['Lh'][:nr] for z in Z], axis=1)          # (nr, N, K, 4)
    N = LH.shape[1]; T = th[-1]
    book = sum(z['book'] for z in Z); btime = sum(z['btime'] for z in Z)
    Ssum = sum(z['Ssum'] for z in Z); Scnt = sum(z['Scnt'] for z in Z)
    chunk = float(Z[0]['chunk']); nseg = int(Z[0]['nseg']); dtau = chunk / nseg
    om = 2 * np.pi * np.fft.rfftfreq(nseg, dtau)
    print(f'# combined {len(Z)} workers, {N} trajectories, common tau <= {T:g} (scaled time at E^ = 1); '
          f'deepest channel max|q^| = {max(float(z["dmax"]) for z in Z):.0f}', flush=True)
    sea = (Scnt[0] + Scnt[1]); tot = Scnt.sum()
    print(f'# time fractions by max depth per chunk: [0,2) {Scnt[0] / tot:.3f}  [2,3) {Scnt[1] / tot:.3f}  '
          f'[3,5) {Scnt[2] / tot:.3f}  [5,inf) {Scnt[3] / tot:.3f}; sea (< 3) {sea / tot:.3f}', flush=True)
    S_sea = (Ssum[0] + Ssum[1]) / max(sea, 1); S_all = Ssum.sum(0) / tot
    print('# 1. drive spectrum S(omega) (two-sided density, mean over the sector\'s drive entries)', flush=True)
    print('#    Kaiser(14) window; the Hann column is the check: where they part, the Hann floor (window leakage) is reached', flush=True)
    print(f"{'omega':>6} {'vec sea':>11} {'vec all':>11} {'blk sea':>11} {'blk all':>11} {'vec sea Hann':>13}", flush=True)
    for w_ in (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20):
        i = np.argmin(np.abs(om - w_))
        print(f'{om[i]:6.2f} {S_sea[0, i]:11.3e} {S_all[0, i]:11.3e} {S_sea[1, i]:11.3e} {S_all[1, i]:11.3e} {S_sea[2, i]:13.3e}', flush=True)
    lo, hi = a.fit_band
    sel = (om >= lo) & (om <= hi)
    for kind, name in ((0, 'vec'), (1, 'blk')):
        for lab, Sx in (('sea', S_sea[kind]), ('all', S_all[kind])):
            ok = sel & (Sx > 0)
            pe = np.polyfit(om[ok], np.log(Sx[ok]), 1); re = np.sqrt(np.mean((np.polyval(pe, om[ok]) - np.log(Sx[ok])) ** 2))
            pp = np.polyfit(np.log(om[ok]), np.log(Sx[ok]), 1); rp = np.sqrt(np.mean((np.polyval(pp, np.log(om[ok])) - np.log(Sx[ok])) ** 2))
            print(f'  {name} {lab}, omega in [{lo:g}, {hi:g}]: exponential S ~ exp(-omega/{-1 / pe[0]:.3f}), rms {re:.3f};  '
                  f'power S ~ omega^{pp[0]:.2f}, rms {rp:.3f}', flush=True)
    h = np.nonzero(th >= T / 2)[0][0]
    print('# 2-3. rates per sector and E (physical lambda_t = E^(1/4) lambda_tau), sea booking, tracking ratio R', flush=True)
    print(f"{'E':>6} {'sect':>4} {'lam_T':>11} {'lam_inc':>11} {'+-':>9} {'lam*T':>7} {'f_sea':>6} {'lam_sea':>11} {'lam_chan':>11} "
          f"{'lam/f_sea':>11} {'R_sea':>9} {'R_all':>9}", flush=True)
    rows = []
    for k, E in enumerate(Es):
        for si, sn in enumerate(SECT):
            Lk = LH[:, :, k, si]
            lT = Lk[-1].mean() / T
            inc = (Lk[-1] - Lk[h]) / (T - th[h]); li, se = inc.mean(), inc.std() / np.sqrt(N)
            lsea = (book[0, k, si] + book[1, k, si]) / max(btime[0] + btime[1], 1)
            lch = (book[2, k, si] + book[3, k, si]) / max(btime[2] + btime[3], 1)
            fsea = (btime[0] + btime[1]) / btime.sum()
            kind = 0 if si == 0 else 1
            i2 = np.argmin(np.abs(om - 2 * Om[k]))
            pred_sea = S_sea[kind, i2] / (8 * Om[k] ** 2); pred_all = S_all[kind, i2] / (8 * Om[k] ** 2)
            s4 = E ** 0.25
            print(f'{E:6.3f} {sn:>4} {s4 * lT:11.4e} {s4 * li:11.4e} {s4 * se:9.1e} {li * T:7.1f} {fsea:6.3f} '
                  f'{s4 * lsea:11.4e} {s4 * lch:11.4e} {s4 * li / fsea:11.4e} {li / pred_sea:9.3g} {li / pred_all:9.3g}', flush=True)
            rows.append((E, si, li, li * T, pred_sea, pred_all))
    print('# [est] if R_sea is flat over the trustworthy rows (lam*T >= 20), S_sea extrapolates lambda:', flush=True)
    for si, sn in enumerate(SECT):
        good = [r for r in rows if r[1] == si and r[3] >= 20]
        if len(good) < 2:
            print(f'  {sn}: fewer than two trustworthy rows', flush=True); continue
        Rs = np.array([r[2] / r[4] for r in good]); Ra = np.array([r[2] / r[5] for r in good])
        kind = 0 if si == 0 else 1
        line = f'  {sn}: R_sea over E = {", ".join(f"{r[0]:g}" for r in good)}: {np.ptp(np.log(Rs)):.2f} spread in ln R ' \
               f'(R_all: {np.ptp(np.log(Ra)):.2f});  R_sea median {np.median(Rs):.3g} -> '
        pr = []
        for E in (0.03, 0.01, 0.003):
            O = E ** -0.25; i2 = np.argmin(np.abs(om - 2 * O))
            pr.append(f'E = {E:g}: {E ** 0.25 * np.median(Rs) * S_sea[kind, i2] / (8 * O * O):.2e}')
        print(line + '; '.join(pr) + ' (physical rate)', flush=True)
    return 'OK'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--E', type=float, nargs='+', default=[1.0, 0.5, 0.3, 0.2, 0.1, 0.07, 0.05])
    ap.add_argument('--workers', type=int, default=22)
    ap.add_argument('--per', type=int, default=24, help='trajectories per worker')
    ap.add_argument('--tau', type=float, default=20000.0, help='scaled time (E^ = 1) after burn-in')
    ap.add_argument('--tau-burn', type=float, default=300.0)
    ap.add_argument('--chunk', type=float, default=20.48, help='spectral segment and booking unit (1024 samples)')
    ap.add_argument('--rtol', type=float, default=1e-9)
    ap.add_argument('--seed', type=int, default=27)
    ap.add_argument('--fit-band', type=float, nargs=2, default=[4.0, 14.0])
    ap.add_argument('--save-every', type=float, default=300.0, help='seconds between worker checkpoints')
    ap.add_argument('--out', default='jeans_out')
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--combine-only', action='store_true')
    a = ap.parse_args()
    T0 = time.time()
    print(f'# START jeans_sea E={a.E} workers={a.workers} per={a.per} tau={a.tau:g} tau_burn={a.tau_burn:g} chunk={a.chunk} '
          f'rtol={a.rtol:g} seed={a.seed} resume={a.resume} combine_only={a.combine_only}  python {sys.version.split()[0]} '
          f'numpy {np.__version__}', flush=True)
    os.makedirs(a.out, exist_ok=True)
    status = 'OK'
    if not a.combine_only:
        os.environ.setdefault('OMP_NUM_THREADS', '1')
        from multiprocessing import Pool
        jobs = [(k, vars(a)) for k in range(a.workers)]
        with Pool(a.workers) as pool:
            for wk, wall, dE, dmax in pool.imap_unordered(worker, jobs):
                print(f'  worker {wk:2d} done: wall {wall:.0f}s, max |dE|/E {dE:.1e}, deepest max|q^| {dmax:.0f}  '
                      f'({time.time() - T0:.0f}s)', flush=True)
    status = combine(a)
    print(f'# END jeans_sea status={status} wall={time.time() - T0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
