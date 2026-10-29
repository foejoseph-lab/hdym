"""cantori_residence.py -- is the halo of the DR island sticky in the cantorus sense?  Pure Yang-Mills, 2 dof.

Question (29 Oct 2026; rabbit_holes/ym/notes/RH9_island_full.md, last section). The islands README calls the zone
around the DR island a "sticky halo", and the usual reading is cantori: remnants of broken KAM tori that slow the flux
across the island boundary.  The signature is in the base flow alone (no ghost): chaotic orbits started just outside
the last KAM curve leave the island's neighbourhood with an ALGEBRAIC tail in their exit times, P(t) ~ t^(-beta) over
decades, instead of the exponential tail of a hyperbolic sea.  This script measures that tail.  If it is absent, the
halo's intermediate FTLE and wide escape spread (hdym_island at E = 1, r = 1e-2) need another explanation, and the
script measures the obvious one too: plain proximity, i.e. how many periods a point started at distance r takes to
leave a neighbourhood of an elliptic orbit at all.

Model: the planar (z = 0, invariant) colour-diagonal system, H = (vx^2 + vy^2)/2 + x^2 y^2/2 at E = 1/2 (scale-free:
every E is this one), the DR orbit of islands/README (published data of the certified orbit, period T = 31.0533).
Section: x = 0 with vx < 0, coordinates (y, vy), vx fixed by the energy.  The orbit crosses it twice per period, at
c = (0.7813288, 0.4813174) and -c (the symmetry (y, vy) -> (-y, -vy)); the section map P sends c -> -c -> c.
Distance from the island: the normalised radius rho = sqrt(dz^T S dz) about the nearer of c, -c, with S the invariant
quadratic form of the linearised period map P^2 at c, scaled to det S = 1 (computed below, printed, checked).  The island is very elongated
in (y, vy) (aspect ~27) and the rotation number of P^2 is ~0.455, close to 5/11, so rho, not Euclidean distance, is
the radius; the island edge in rho is measured by `--section edge`.

Integration: Yoshida-4, fixed step --h (default 0.01; omega h <= 0.035 inside the tube), vectorised over trajectories,
section crossings by cubic Hermite interpolation.  A trajectory leaves when its section point has rho > the largest
--rho-out, or when max(|x|, |y|) > --wander (3.5; the DR orbit reaches 3.147), which counts as leaving at every
threshold.  Exit times are recorded only at crossings (every ~15.5 time units).  Pure YM: no ghost anywhere.

What counts as a result, fixed before the production run (29 Oct 2026, after a 4,000-start smoke run to t = 20,000):
  positive (cantori): exits past the bulk continue over at least two decades of t (from ~10^3 to ~10^5) with a roughly
    constant local log-log slope, far above what an exponential at the bulk's rate leaves; the long-lived orbits sit in
    a thin layer just outside the last KAM curve (`--section trace`).
  negative: after the bulk (exponential, a few hundred time units) exits stop, apart from a number compatible with
    the bulk's exponential, and the plateau is flat; then the halo's FTLE and spread are proximity, which `--section
    shells` measures directly.
  The smoke run already shows the bulk (local slopes rising 1.5 -> 5 by t ~ 700, i.e. exponential) and a residue of a
  few starts in 10^3 that live 10^3-1.5 x 10^4 inside rho ~ 3-7 x 10^-4; the production run decides which it is.

Sections
  --section edge      normalised-radius scan from the centre, 22 radii x 8 angles, to --tmax-edge: the island edge.
  --section survive   N starts uniform in area in the annulus rho in [--ra, --rb] about c (the sticky layer), first
                      exit time past each --rho-out, to --tmax; chunks of --chunk starts over --workers processes,
                      one npz per chunk (resumable with --resume).  Starts that never exit by tmax are the regular
                      plateau (KAM curves and satellite island chains); the tail of the rest is the test.
  --section shells    the same for Euclidean discs (y, vy) of radius --shells about c: the proximity control, the
                      2-dof analogue of hdym_island's shells.
  --section trace     re-integrates the --trace-k longest-lived exiting starts of a finished survive run and saves
                      their section points (for a picture of where the long-lived orbits sit).
CHECKs: the fixed point of P^2 (scipy DOP853, rtol 1e-13) and det M = 1; the same fixed point under the Yoshida
stepper at --h; energy drift on every chunk; an h/2 rerun of the first chunk with the same starts (--check-h), whose
exit-time counts, to min(tmax, --tmax-check), must agree with the h run in every bin within 3 binomial sigma + 2
(chaotic orbits differ one by one between step sizes; their ensemble statistics must not).

Usage (Spyder, from the repo root):
    !"{sys.executable}" -B -u islands\\code\\cantori_residence.py --section edge > islands\\results\\cantori_edge.txt
    !"{sys.executable}" -B -u islands\\code\\cantori_residence.py --section survive --N 200000 --tmax 300000 --workers 22 --resume >> islands\\results\\cantori_survive.txt
    !"{sys.executable}" -B -u islands\\code\\cantori_residence.py --section shells --N 20000 --tmax 30000 --workers 22 --resume >> islands\\results\\cantori_shells.txt
    !"{sys.executable}" -B -u islands\\code\\cantori_residence.py --section trace --trace-k 12 > islands\\results\\cantori_trace.txt
Cost (estimates from smoke runs here: 2,000 starts to t = 20,000 in 77 s per chunk, the plateau of ~14 % regular
starts dominating): survive at the defaults ~1.5-2 h on 22 workers plus ~5 min for the h/2 check; shells, edge and
trace a few minutes each.  Added 2026-10-29.
"""
import argparse, glob, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
X0, VX0, VY0 = 3.14640122769753, 0.0017931934191, 0.99999839222738     # interior.py / scan.py (published DR data)
T_DR = 31.053330268585064
CBRT = 2 ** (1 / 3)
W1 = 1 / (2 - CBRT); W0 = -CBRT / (2 - CBRT)
YC = (W1 / 2, (W0 + W1) / 2, (W0 + W1) / 2, W1 / 2)
YD = (W1, W0, W1, 0.0)


# ------------------------------------------------------------------ geometry of the section
def _rhs(t, s):
    return [s[2], s[3], -s[0] * s[1] ** 2, -s[1] * s[0] ** 2]


def _ev(t, s):
    return s[0]
_ev.direction = -1


def p2_scipy(z):
    """Second return of (y, vy) on x = 0, vx < 0 (one DR period), with DOP853 at rtol 1e-13."""
    from scipy.integrate import solve_ivp
    y, vy = z
    r = solve_ivp(_rhs, (0, 1.4 * T_DR), [0.0, y, -np.sqrt(1 - vy * vy), vy], rtol=1e-13, atol=1e-15, events=_ev)
    ev = [e for te, e in zip(r.t_events[0], r.y_events[0]) if te > 1.0]
    return np.array([ev[1][1], ev[1][3]])


def geometry():
    from scipy.integrate import solve_ivp
    r = solve_ivp(_rhs, (0, T_DR), [X0, 0.0, VX0, VY0], rtol=1e-13, atol=1e-15, events=_ev)
    z = np.array([r.y_events[0][0][1], r.y_events[0][0][3]])
    for _ in range(4):
        F = p2_scipy(z) - z
        J = np.zeros((2, 2)); hh = 1e-7
        for k in range(2):
            d = np.zeros(2); d[k] = hh
            J[:, k] = (p2_scipy(z + d) - p2_scipy(z - d)) / (2 * hh)
        z = z - np.linalg.solve(J - np.eye(2), F)
    res = np.abs(p2_scipy(z) - z).max()
    a, b, c, d = J.ravel()
    S = np.array([[c, (d - a) / 2], [(d - a) / 2, -b]])
    if S[0, 0] < 0:
        S = -S
    S /= np.sqrt(np.linalg.det(S))                      # det S = 1: rho is an area-normalised radius
    return z, J, S, res


# ------------------------------------------------------------------ the stepper
def step(x, y, vx, vy, h):
    for c, d in zip(YC, YD):
        x = x + c * h * vx; y = y + c * h * vy
        if d:
            vx = vx - d * h * x * y * y; vy = vy - d * h * y * x * x
    return x, y, vx, vy


def _hermite_root(x0, x1, v0, v1, h):
    s = x0 / (x0 - x1)
    for _ in range(4):
        f = (2*s**3 - 3*s**2 + 1)*x0 + (s**3 - 2*s**2 + s)*h*v0 + (-2*s**3 + 3*s**2)*x1 + (s**3 - s**2)*h*v1
        fp = (6*s**2 - 6*s)*x0 + (3*s**2 - 4*s + 1)*h*v0 + (-6*s**2 + 6*s)*x1 + (3*s**2 - 2*s)*h*v1
        s = np.clip(s - f / fp, 0.0, 1.0)
    return s


def _herm(a0, a1, b0, b1, h, s):
    return (2*s**3 - 3*s**2 + 1)*a0 + (s**3 - 2*s**2 + s)*h*b0 + (-2*s**3 + 3*s**2)*a1 + (s**3 - s**2)*h*b1


class Geo:
    def __init__(self, z, S):
        self.z = np.asarray(z); self.S = np.asarray(S)
        self.w, self.U = np.linalg.eigh(self.S)

    def rho(self, ys, vys):
        out = np.full(ys.shape, np.inf)
        for c in (self.z, -self.z):
            d0, d1 = ys - c[0], vys - c[1]
            out = np.minimum(out, np.sqrt(self.S[0, 0]*d0*d0 + 2*self.S[0, 1]*d0*d1 + self.S[1, 1]*d1*d1))
        return out

    def annulus(self, ra, rb, n, rng):
        r = np.sqrt(rng.uniform(ra * ra, rb * rb, n)); th = rng.uniform(0, 2 * np.pi, n)
        u = np.stack([np.cos(th) / np.sqrt(self.w[0]), np.sin(th) / np.sqrt(self.w[1])], 1) * r[:, None]
        dz = u @ self.U.T
        return self.z[0] + dz[:, 0], self.z[1] + dz[:, 1]

    def disc(self, r, n, rng):
        rr = r * np.sqrt(rng.uniform(0, 1, n)); th = rng.uniform(0, 2 * np.pi, n)
        return self.z[0] + rr * np.cos(th), self.z[1] + rr * np.sin(th)


def evolve(geo, y0, vy0, thr, tmax, h, wander, keep_points=False):
    """First time each trajectory's section point has rho > thr[k]; inf if never by tmax.  Returns (first, rho_max,
    wandered, dE_max, points) -- points: list of (y, vy) arrays per trajectory, only with keep_points."""
    ok = np.abs(vy0) < 1
    n = y0.size
    x = np.zeros(n); y = y0.copy(); vy = vy0.copy(); vx = -np.sqrt(np.where(ok, 1 - vy0 ** 2, 1.0))
    thr = np.asarray(thr, float)
    first = np.full((n, thr.size), np.inf); first[~ok] = 0.0
    rmax = np.zeros(n); wand = np.zeros(n, bool)
    E0 = 0.5 * (vx ** 2 + vy ** 2); dE = np.zeros(n)
    pts = [[] for _ in range(n)] if keep_points else None
    idx = np.nonzero(ok)[0]; x, y, vx, vy, E0s = x[idx], y[idx], vx[idx], vy[idx], E0[idx]
    t = 0.0; nst = 0
    while t < tmax - 1e-9 and idx.size:
        xn, yn, vxn, vyn = step(x, y, vx, vy, h); t += h; nst += 1
        cr = (x > 0) & (xn <= 0)
        if cr.any():
            s = _hermite_root(x[cr], xn[cr], vx[cr], vxn[cr], h)
            ys = _herm(y[cr], yn[cr], vy[cr], vyn[cr], h, s)
            vys = _herm(vy[cr], vyn[cr], -y[cr] * x[cr] ** 2, -yn[cr] * xn[cr] ** 2, h, s)
            r = geo.rho(ys, vys); ci = idx[cr]
            rmax[ci] = np.maximum(rmax[ci], r)
            sub = first[ci]; sub[(r[:, None] > thr[None, :]) & np.isinf(sub)] = t; first[ci] = sub
            if keep_points:
                for k, a_, b_ in zip(ci, ys, vys):
                    pts[k].append((t, a_, b_))
        wd = np.maximum(np.abs(xn), np.abs(yn)) > wander
        if wd.any():
            wi = idx[wd]; wand[wi] = True
            sub = first[wi]; sub[np.isinf(sub)] = t; first[wi] = sub
        alive = np.isinf(first[idx, -1])
        if not alive.all():
            Eg = 0.5 * (vxn[~alive] ** 2 + vyn[~alive] ** 2) + 0.5 * xn[~alive] ** 2 * yn[~alive] ** 2
            dE[idx[~alive]] = np.abs(Eg - E0s[~alive])
        idx = idx[alive]; x, y, vx, vy, E0s = xn[alive], yn[alive], vxn[alive], vyn[alive], E0s[alive]
    if idx.size:
        Eg = 0.5 * (vx ** 2 + vy ** 2) + 0.5 * x ** 2 * y ** 2
        dE[idx] = np.abs(Eg - E0s)
    return first, rmax, wand, dE, pts


# ------------------------------------------------------------------ chunks
def _chunk_job(args):
    k, kind, a_dict, z, S = args
    a = argparse.Namespace(**a_dict)
    geo = Geo(z, S)
    out = os.path.join(a.work, f'{kind}_c{k:04d}.npz')
    if a.resume and os.path.exists(out):
        return k, out, 0.0, True
    t0 = time.time()
    rng = np.random.default_rng(a.seed + 1000 * k)
    if kind == 'survive':
        y0, vy0 = geo.annulus(a.ra, a.rb, a.chunk, rng); tag = np.zeros(a.chunk)
    else:
        per = a.chunk // len(a.shells)
        ys, vs, tg = [], [], []
        for j, r in enumerate(a.shells):
            yy, vv = geo.disc(r, per, rng); ys.append(yy); vs.append(vv); tg.append(np.full(per, j))
        y0, vy0, tag = np.concatenate(ys), np.concatenate(vs), np.concatenate(tg)
    first, rmax, wand, dE, _ = evolve(geo, y0, vy0, a.rho_out, a.tmax, a.h, a.wander)
    tmp = out + '.tmp.npz'
    np.savez(tmp, y0=y0, vy0=vy0, tag=tag, first=first, rmax=rmax, wand=wand, dE=dE, h=a.h, tmax=a.tmax,
             rho_out=np.array(a.rho_out), seed=a.seed + 1000 * k)
    os.replace(tmp, out)
    return k, out, time.time() - t0, False


TBINS = np.array([30, 50, 100, 200, 300, 500, 1e3, 2e3, 3e3, 5e3, 1e4, 2e4, 3e4, 5e4, 1e5, 2e5, 3e5, 5e5, 1e6])


def table(first, tmax, rho_out, label):
    n = first.shape[0]
    ts = TBINS[TBINS <= tmax]
    print(f'# {label}: {n} starts; counts with exit time > t (never = no exit by tmax = {tmax:g})', flush=True)
    print('  rho_out   ' + ' '.join(f'{t:>8g}' for t in ts) + '    never', flush=True)
    for j, r in enumerate(rho_out):
        f = first[:, j]
        print(f'  {r:8.1e} ' + ' '.join(f'{int(np.sum(f > t)):8d}' for t in ts) + f' {int(np.isinf(f).sum()):8d}', flush=True)


def tail_fits(first, tmax, rho_out):
    """For each rho_out: among starts that exit by tmax, survival S(t) = P(exit > t | exits); fits over the window
    where S has between 300 and 10 counts, of log S on log t (power) and on t (exponential), with rms residuals,
    and the local log-log slope per factor of 2-3 in t.  Also the plateau (never exits)."""
    print('# tail fits (exiting starts only; window: 300 >= count >= 10)', flush=True)
    for j, r in enumerate(rho_out):
        f = first[:, j]; f = np.sort(f[np.isfinite(f)]); n = f.size
        if n < 50:
            print(f'  rho_out {r:.1e}: {n} exiting starts, too few', flush=True); continue
        cnt = n - np.arange(n)                         # number with exit >= f[i]
        sel = (cnt <= 300) & (cnt >= 10)
        if sel.sum() < 5:
            print(f'  rho_out {r:.1e}: window too short', flush=True); continue
        tt, ss = f[sel], cnt[sel] / n
        pp = np.polyfit(np.log(tt), np.log(ss), 1); rp = np.sqrt(np.mean((np.polyval(pp, np.log(tt)) - np.log(ss)) ** 2))
        pe = np.polyfit(tt, np.log(ss), 1); re = np.sqrt(np.mean((np.polyval(pe, tt) - np.log(ss)) ** 2))
        print(f'  rho_out {r:.1e}: n_exit {n}, window t = {tt[0]:.0f} .. {tt[-1]:.0f} ({tt[-1] / tt[0]:.1f}x); '
              f'power: beta = {-pp[0]:.3f}, rms {rp:.3f};  exponential: rate = {-pe[0]:.3e}, rms {re:.3f}', flush=True)
        edges = np.unique(np.concatenate([[f[0]], np.geomspace(max(f[0], 30), f[-1], 12)]))
        Sv = np.array([np.mean(f > e) for e in edges])
        loc = [f'{edges[i]:.0f}:{-np.log(Sv[i + 1] / Sv[i]) / np.log(edges[i + 1] / edges[i]):.2f}'
               for i in range(len(edges) - 1) if Sv[i + 1] > 0 and Sv[i] > 0 and n * Sv[i + 1] >= 5]
        print(f'           local log-log slopes (t: -dlnS/dlnt): ' + ' '.join(loc), flush=True)
        dec = [(lo, hi, int(np.sum((f >= lo) & (f < hi)))) for lo, hi in ((1e3, 1e4), (1e4, 1e5), (1e5, 1e6)) if lo < tmax]
        rate = -pe[0]; ex = np.sum(f >= 1e3)
        print(f'           exits per decade: ' + ', '.join(f'[{lo:g}, {min(hi, tmax):g}): {c}' for lo, hi, c in dec)
              + f'; an exponential at the fitted rate would leave {n * np.exp(np.polyval(pe, 1e3)):.2g} past t = 1000 '
              f'(observed {ex}); plateau (no exit by tmax) {int(np.sum(np.isinf(first[:, j])))}', flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--section', choices=['edge', 'survive', 'shells', 'trace'], required=True)
    ap.add_argument('--N', type=int, default=200000)
    ap.add_argument('--chunk', type=int, default=2000)
    ap.add_argument('--tmax', type=float, default=300000.0)
    ap.add_argument('--tmax-edge', type=float, default=30000.0)
    ap.add_argument('--tmax-check', type=float, default=30000.0, help='the h/2 check runs chunk 0 to min(tmax, this)')
    ap.add_argument('--h', type=float, default=0.01)
    ap.add_argument('--ra', type=float, default=2e-4, help='annulus inner rho (survive)')
    ap.add_argument('--rb', type=float, default=8e-4, help='annulus outer rho (survive)')
    ap.add_argument('--rho-out', type=float, nargs='+', default=[1.5e-3, 3e-3, 1e-2, 3e-2, 1e-1])
    ap.add_argument('--shells', type=float, nargs='+', default=[2e-4, 5e-4, 1e-3, 3e-3, 1e-2],
                    help='Euclidean (y, vy) disc radii about c (shells)')
    ap.add_argument('--wander', type=float, default=3.5)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--workers', type=int, default=22)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--check-h', action='store_true', default=True)
    ap.add_argument('--no-check-h', dest='check_h', action='store_false')
    ap.add_argument('--trace-k', type=int, default=12)
    ap.add_argument('--results', default=os.path.join(HERE, '..', 'results'))
    a = ap.parse_args()
    a.rho_out = sorted(a.rho_out)
    T0 = time.time()
    print(f'# START cantori_residence section={a.section} N={a.N} chunk={a.chunk} tmax={a.tmax:g} h={a.h} '
          f'ra={a.ra:g} rb={a.rb:g} rho_out={a.rho_out} shells={a.shells} wander={a.wander} seed={a.seed} '
          f'workers={a.workers} resume={a.resume}  python {sys.version.split()[0]} numpy {np.__version__}', flush=True)
    status = 'OK'
    z, M, S, res = geometry()
    tr = np.trace(M) / 2
    print(f'geometry: section x = 0, vx < 0; fixed point of P^2 c = ({z[0]:.12f}, {z[1]:.12f}); |P^2(c) - c| = {res:.1e}; '
          f'M = [[{M[0,0]:.6f}, {M[0,1]:.6f}], [{M[1,0]:.6f}, {M[1,1]:.6f}]], det M = {np.linalg.det(M):.9f}, '
          f'tr M / 2 = {tr:.6f}, rotation number {np.arccos(tr) / (2 * np.pi):.6f}; S (det 1) = '
          f'[[{S[0,0]:.6f}, {S[0,1]:.6f}], [{S[1,0]:.6f}, {S[1,1]:.6f}]], aspect {np.sqrt(np.linalg.eigvalsh(S)[1] / np.linalg.eigvalsh(S)[0]):.1f}',
          flush=True)
    ok1 = res < 1e-10 and abs(np.linalg.det(M) - 1) < 1e-5
    print(f'CHECK fixed point (DOP853) and area preservation: {"PASS" if ok1 else "FAIL"}', flush=True)
    status = status if ok1 else 'CHECK-FAIL'
    geo = Geo(z, S)
    f0, rm0, _, dE0, pts = evolve(geo, np.array([z[0]]), np.array([z[1]]), [1.0], 20 * T_DR, a.h, a.wander, keep_points=True)
    P0 = np.array(pts[0]); dev = float(geo.rho(P0[:, 1], P0[:, 2]).max())
    ok2 = dev < 1e-5
    print(f'CHECK the fixed point under the Yoshida stepper at h = {a.h}: max rho of its section points over 20 periods '
          f'{dev:.1e} (tol 1e-5, a twentieth of the island edge; the O(h^4) shift of the orbit), |dE| {dE0[0]:.1e}: '
          f'{"PASS" if ok2 else "FAIL"}', flush=True)
    status = status if ok2 else 'CHECK-FAIL'
    os.makedirs(a.results, exist_ok=True)

    if a.section == 'edge':
        rs = np.logspace(-6, -2.5, 22); th = np.linspace(0, 2 * np.pi, 8, endpoint=False)
        R, TH = np.meshgrid(rs, th); R, TH = R.ravel(), TH.ravel()
        u = np.stack([np.cos(TH) / np.sqrt(geo.w[0]), np.sin(TH) / np.sqrt(geo.w[1])], 1) * R[:, None]
        dz = u @ geo.U.T
        first, rmax, wand, dE, _ = evolve(geo, z[0] + dz[:, 0], z[1] + dz[:, 1], [0.05], a.tmax_edge, a.h, a.wander)
        print(f'# edge scan to t = {a.tmax_edge:g}: per rho0, 8 angles; left = section rho > 0.05 or wander', flush=True)
        for r in rs:
            m = R == r
            print(f'  rho0 {r:.2e}: left {int(np.isfinite(first[m, 0]).sum())}/8, median exit {np.median(first[m, 0]):9.1f}, '
                  f'median rho_max/rho0 {np.median(rmax[m] / r):9.2f}, max |dE| {dE[m].max():.1e}', flush=True)
        np.savez(os.path.join(a.results, 'cantori_edge.npz'), rho0=R, theta=TH, first=first[:, 0], rmax=rmax, dE=dE,
                 z=z, M=M, S=S, tmax=a.tmax_edge, h=a.h)

    elif a.section in ('survive', 'shells'):
        kind = a.section
        a.work = os.path.join(a.results, f'cantori_{kind}')
        os.makedirs(a.work, exist_ok=True)
        nchunk = int(np.ceil(a.N / a.chunk))
        ad = vars(a).copy()
        jobs = [(k, kind, ad, z, S) for k in range(nchunk)]
        from multiprocessing import Pool
        with Pool(a.workers) as pool:
            for k, path, wall, skipped in pool.imap_unordered(_chunk_job, jobs):
                print(f'  chunk {k:4d} {"(resumed)" if skipped else f"{wall:7.0f}s"}  ({time.time() - T0:.0f}s)', flush=True)
        Z = [np.load(p) for p in sorted(glob.glob(os.path.join(a.work, f'{kind}_c*.npz')))]
        Z = [q for q in Z if float(q['tmax']) == a.tmax and float(q['h']) == a.h]
        first = np.concatenate([q['first'] for q in Z]); tag = np.concatenate([q['tag'] for q in Z])
        wand = np.concatenate([q['wand'] for q in Z]); dE = np.concatenate([q['dE'] for q in Z])
        print(f'# {len(Z)} chunks, {first.shape[0]} starts; wandered (max(|x|,|y|) > {a.wander}) {int(wand.sum())}; '
              f'max |dE|/E {dE.max() / 0.5:.1e}', flush=True)
        if kind == 'survive':
            table(first, a.tmax, a.rho_out, f'annulus rho in [{a.ra:g}, {a.rb:g}]')
            tail_fits(first, a.tmax, a.rho_out)
        else:
            for j, r in enumerate(a.shells):
                m = tag == j
                table(first[m], a.tmax, a.rho_out, f'Euclidean disc r = {r:g} about c')
        if a.check_h and Z:
            q = Z[0]; tch = min(a.tmax, a.tmax_check)
            f2, _, _, _, _ = evolve(geo, q['y0'], q['vy0'], a.rho_out, tch, a.h / 2, a.wander)
            f1 = np.where(q['first'] > tch, np.inf, q['first']); n0 = f1.shape[0]
            ts = TBINS[TBINS < tch]; worst, wtxt, okh = 0.0, '', True
            for j in range(len(a.rho_out)):
                for t in list(ts) + [np.inf]:
                    c1 = np.sum(f1[:, j] > t) if np.isfinite(t) else np.isinf(f1[:, j]).sum()
                    c2 = np.sum(f2[:, j] > t) if np.isfinite(t) else np.isinf(f2[:, j]).sum()
                    tol = 3 * np.sqrt(c1 * (1 - c1 / n0) + c2 * (1 - c2 / n0)) + 2
                    if abs(c1 - c2) / tol > worst:
                        worst, wtxt = abs(c1 - c2) / tol, f'rho_out {a.rho_out[j]:.1e}, t > {t:g}: {c1} vs {c2}'
                    okh &= abs(c1 - c2) <= tol
            print(f'CHECK h/2 rerun of chunk 0 ({n0} starts, to t = {tch:g}): exit counts per bin agree within 3 sigma '
                  f'binomial + 2 (chaotic orbits differ one by one, ensembles must not), worst {worst:.2f} of tol '
                  f'({wtxt}): {"PASS" if okh else "FAIL"}', flush=True)
            status = status if okh else 'CHECK-FAIL'
        np.savez(os.path.join(a.results, f'cantori_{kind}.npz'), first=first, tag=tag, wand=wand, dE=dE,
                 y0=np.concatenate([q['y0'] for q in Z]), vy0=np.concatenate([q['vy0'] for q in Z]),
                 rho_out=np.array(a.rho_out), shells=np.array(a.shells), ra=a.ra, rb=a.rb, tmax=a.tmax, h=a.h, z=z, M=M, S=S)

    elif a.section == 'trace':
        src = np.load(os.path.join(a.results, 'cantori_survive.npz'))
        f = src['first'][:, -1]; fin = np.isfinite(f)
        order = np.argsort(-np.where(fin, f, -1))[:a.trace_k]
        tmax = float(f[order[0]]) + 100
        _, _, _, _, pts = evolve(geo, src['y0'][order], src['vy0'][order], src['rho_out'], tmax, float(src['h']), a.wander,
                                 keep_points=True)
        out = {}
        for i, k in enumerate(order):
            P = np.array(pts[i]); rr = geo.rho(P[:, 1], P[:, 2])
            out[f'traj{i}'] = P
            q = np.percentile(rr, [10, 50, 90])
            print(f'  start {k}: exit {f[k]:.0f}, {len(P)} crossings, rho percentiles 10/50/90 = '
                  f'{q[0]:.2e}/{q[1]:.2e}/{q[2]:.2e}', flush=True)
        np.savez(os.path.join(a.results, 'cantori_trace.npz'), index=order, z=z, S=S, **out)
    print(f'# END cantori_residence section={a.section} status={status} wall={time.time() - T0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
