"""offdiag_escape.py -- nonlinear escape times of the FULL nine-component Lee-Wick system seeded with
off-diagonal perturbations of the Dahlqvist-Russberg orbit: the nonlinear counterpart of the typed
Floquet scan (dr_offdiag_floquet_typed.py, B001).

Model: ym/hdym_full_derive.py (all nine A^a_i, temporal gauge, g = M = 1; its CHECK 3 shows the
linearisation about the diagonal background is exactly rung1_system.pkl's K, so this is the same system
the Floquet scans analysed, now without the linearisation).  Integrator: ym/csrc/hdym_full_escape.c
(Dormand-Prince 5(4), adaptive, OpenMP).  Nothing here is symplectic or certified; the Ostrogradsky
energy and the Gauss constraints are recorded per trajectory as the drift diagnostics.

Initial data.  The diagonal sector is placed ON the DR orbit at energy E (same point and scaling as
dr_offdiag_floquet.orbit) with q2 = -grad V + w, q3 = -Hess V q1 + w1 (qk = k-th time derivative), where
the diagonal ghost (w, w1) is a random vector of amplitude delta_diag (default = delta; --delta-diag 0
switches it off, but then the integrator's truncation error seeds it at ~rtol anyway: the ghost-free
manifold is Floquet-unstable on the orbit for E >~ 0.315 and in tongues below, and an unseeded run at
E = 1 escapes through the diagonal ghost at t ~ 120 from a 1e-10 truncation seed).  So the experiment is
a race between the two sectors seeded at equal amplitude; the escaping sector and the growth curves say
who won.  The diagonal ghost's own Floquet rate on the orbit is computed per E (6-dim Hill monodromy) and
stored as rate_diag.  The off-diagonal state (u and its first three derivatives for the six components,
24 numbers) is a random vector of norm delta, projected onto the linearised Gauss constraint (G of
rung1_system.pkl at the DR point); the full nonlinear Gauss residual is then O(delta^2) and is printed.
--block restricts the seed to one block.

Sections
  --section E       nE energies (log-spaced in [Emin, Emax]) x nseed random seeds each.  Observable:
                    t_esc(E, seed), the escaping sector, and the sector norms at nchk checkpoints (the
                    nonlinear growth curves).  Compare with the linear prediction
                    t_lin = ln(R / delta) / rate from the typed log (drawn by --plot if the log is present).
  --section phase   fixed E and --block; an ensemble (nalpha phases x nseed random admixtures and diagonal seeds)
                    over the phase alpha of the seed's projection on a Krein quadruplet's unstable subspace,
                    |c| fixed (rabbit hole 1; run_phase's docstring).  Added 2026-10-21.
  --section orient  fixed E; a 2D grid over the sphere angles (theta, phi) of the seed direction in the
                    span of three fixed random Gauss-projected directions e1, e2, e3 (seeded by --seed):
                    seed = delta (cos theta e1 + sin theta cos phi e2 + sin theta sin phi e3).  The
                    escape-time map in the style of hdym_escape_map's --section ghost, for the off-diagonal
                    sector.

Linear transients can be large: near a Krein band edge the two block-12 modes are nearly parallel and their beat
amplifies the seed by up to ~10^3 without instability (E = 0.05: seed 1e-4 -> 0.23 at t ~ 2150 and back, energy
conserved).  The growth curves (chk) show it; the "driver" of a runaway is therefore read at norm 1, not at the seed.
A blow-up too sharp to follow to R (step underflow or overflow with the state already beyond |A| = 10) counts as
an escape at the last resolved time and is flagged status = 1 ("singular"; the Lee-Wick runaway can reach a
finite-time singularity, cf. hdym_core.confirm_runaway); status = 2 is a real integrator failure and fails the run.

Since 2026-10-21 every section also stores tcross (N, 5): the first times |w| >= 1, |u_12|, |u_13|, |u_23| >= 1 and
|w| >= 0.1 (-1 = never; rabbit hole 4), and the orient section keeps its chk array (float32; rabbit hole 12).
Output: <results>/<tag>.npz, written after every row (E value, or theta row) so a killed run keeps what
it has (a partial file carries rows_done < rows); `--resume` with the same parameters continues it from there
exactly (the seed basis, the fixed diagonal seed and the per-cell seeds are all regenerated from --seed).  Use
`>>` on the log when resuming, or the restart overwrites it.  Log: CHECK lines first (build, integrator
cross-check against scipy DOP853 on one trajectory, Gauss residual of the seeds), then one line per row,
then DONE <seconds>.  --plot <tag> draws the figures from an existing npz (no run).

Usage (Spyder, from the repo root, kernel's own interpreter; the derivation must have been run once):
    import sys; %cd C:\\GitHub\\hdym
    !"{sys.executable}" -B -u ym\\hdym_full_derive.py > islands\\results\\hdym_full_derive.txt 2>&1
    !"{sys.executable}" -B -u islands\\code\\offdiag_escape.py --section E --nE 100 --nseed 32 --tag esc_E > islands\\results\\esc_E.txt 2>&1
    !"{sys.executable}" -B -u islands\\code\\offdiag_escape.py --section orient --E 0.156 --n 128 --tag esc_orient_0156 > islands\\results\\esc_orient_0156.txt 2>&1
    !"{sys.executable}" -B -u islands\\code\\offdiag_escape.py --section orient --E 0.30 --n 128 --tag esc_orient_030 > islands\\results\\esc_orient_030.txt 2>&1
  (orient only pays where the off-diagonal sector leads the race, i.e. where the diagonal ghost is Floquet-stable on the
  orbit: E = 0.156 is block 12's Krein band, E = 0.30 block 23's parametric tongue.  At E = 0.412 and 0.642 the diagonal
  ghost's rate is 2-4x the off-diagonal one and the map is structureless: 20 Oct run at 0.412, median 183, 96 % diagonal.)
    !"{sys.executable}" -B -u islands\\code\\offdiag_escape.py --plot esc_E

Run parameters and why: delta = delta_diag = 1e-4.  Convergence test (E = 0.549, both seeds, same initial
state, escape time): delta 1e-4 -> 151.64 / 151.58 / 151.59 at rtol 1e-10 / 1e-11 / 1e-12 against 151.61 from
scipy DOP853 at 1e-12 (converged, 4e-4); delta 1e-5 -> 198.0 / 200.2 / 200.2 vs 200.1 (needs 1e-11);
delta 1e-6 -> 231.3 / 238.8 / 240.9 vs 238.3 (not converged at any practical tolerance: the ghost-free manifold is
unstable and the integrator's truncation error competes with the seed).  Hence delta = 1e-4 with rtol 1e-11,
atol 1e-13 (the escape-time CHECK must pass to 1e-3); ln(1/delta) = 9.2, ln(R/delta) = 18.4 e-foldings.  The
paper's regime (ghost fraction eps = 0.01, |w| ~ 0.1) is larger still.  R = 1e4 (hdym_escape_map); tmax = 3000
(hdym_island at E = 0.3); hmax 0.25; nchk 60.  Cost: ~0.1-1 s per trajectory single-threaded.
Added 2026-10-19.
"""
import argparse
import ctypes
import os
import pickle
import subprocess
import sys
import tempfile
import time

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'ym'))
import dr_offdiag_floquet as base            # noqa: E402  (orbit, PKL, BLOCKS)
import cbuild                                # noqa: E402

OD_M = [1, 2, 3, 5, 6, 7]                    # m = 3 i + a of the off-diagonal components, rung1 ulist order
BLOCK_M = {'12': [1, 3], '13': [2, 6], '23': [5, 7]}
SECTOR_NAMES = ['diagonal', 'block 12', 'block 13', 'block 23']


# ---------------------------------------------------------------- the kernel
class Kernel:
    def __init__(self):
        src = os.path.join(ROOT, 'ym', 'csrc', 'hdym_full_escape.c')
        hdr = os.path.join(ROOT, 'ym', 'csrc', 'hdym_full_rhs.h')
        if not os.path.exists(hdr):
            raise SystemExit("ym/csrc/hdym_full_rhs.h missing: run ym/hdym_full_derive.py first")
        d = tempfile.mkdtemp()
        lib = os.path.join(d, 'hdym_full_escape' + ('.dll' if sys.platform == 'win32' else '.so'))
        cc = os.environ.get('CC', 'gcc')
        inc = '-I' + os.path.join(ROOT, 'ym', 'csrc')
        try:
            cbuild.check_call([cc, '-O3', '-fopenmp', '-shared', '-fPIC', inc, '-o', lib, src, '-lm'])
            omp = 'with OpenMP'
        except subprocess.CalledProcessError:
            cbuild.check_call([cc, '-O3', '-shared', '-fPIC', inc, '-o', lib, src, '-lm'])
            omp = 'single-threaded (no OpenMP)'
        self.lib = ctypes.CDLL(lib)
        dp, ip = ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int)
        self.lib.escape.argtypes = [dp, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, dp, ip, ip, dp, dp, dp, ip]
        self.lib.integrate_to.argtypes = [dp, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        self.lib.omp_threads.restype = ctypes.c_int
        self.lib.fp_mode_report.argtypes = [ctypes.POINTER(ctypes.c_uint)]
        rep = (ctypes.c_uint * 3)()
        self.lib.fp_mode_report(rep)
        print(f"CHECK kernel built with {cc} {omp}, {self.lib.omp_threads()} threads: PASS", flush=True)
        mx = (rep[0] & 0xFFFF, rep[1] & 0xFFFF)
        x87 = (rep[0] >> 16, rep[1] >> 16)
        print(f"INFO  FP state the threads had before the kernel pinned it: MXCSR control bits min 0x{mx[0] & 0xFFC0:04X} max 0x{mx[1] & 0xFFC0:04X}"
              f" (status flags ignored); x87 control word min 0x{x87[0]:04X} max 0x{x87[1]:04X} (Windows only; 0 elsewhere)."
              f"  Now MXCSR 0x{rep[2]:04X}, x87 53-bit, in every thread; the step controller no longer calls the C library.", flush=True)

    def escape(self, S0, cj, tmax, R, rtol, atol, hmax, nchk):
        S = np.array(S0, dtype=np.float64, order='C', copy=True)   # (N, 36); the kernel overwrites it with the final state
        N = S.shape[0]
        tesc = np.zeros(N)
        sector = np.zeros(N, dtype=np.int32)
        status = np.zeros(N, dtype=np.int32)
        thread_of = np.zeros(N, dtype=np.int32)
        chk = np.zeros((N, nchk, 6))
        diag = np.zeros((N, 4))
        tcross = np.zeros((N, 5))                               # first crossings: |w|>=1, |u_12|>=1, |u_13|>=1, |u_23|>=1, |w|>=0.1 (-1 = never)
        dp, ip = ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int)
        self.lib.escape(S.ctypes.data_as(dp), N, cj, tmax, R, rtol, atol, hmax, nchk,
                        tesc.ctypes.data_as(dp), sector.ctypes.data_as(ip), status.ctypes.data_as(ip), chk.ctypes.data_as(dp), diag.ctypes.data_as(dp),
                        tcross.ctypes.data_as(dp), thread_of.ctypes.data_as(ip))
        self.last_thread_of = thread_of
        self.last_tcross = tcross                               # (N, 5); the sections save it as 'tcross'
        return tesc, sector, chk, diag, S, status

    def integrate_to(self, s0, cj, T, rtol, atol, hmax):
        s = np.ascontiguousarray(s0, dtype=np.float64).copy()
        self.lib.integrate_to(s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), cj, T, rtol, atol, hmax)
        return s


# ---------------------------------------------------------------- the model's numpy side (checks and seeds)
class Model:
    def __init__(self):
        d = pickle.load(open(os.path.join(ROOT, 'ym', 'full_system.pkl'), 'rb'))
        args = list(d['state']) + [d['cj']]
        self.R = sp.lambdify(args, d['R'], 'numpy')
        self.G = sp.lambdify(args, d['G'], 'numpy')
        self.E = sp.lambdify(args, d['E'], 'numpy')
        r1 = pickle.load(open(base.PKL, 'rb'))
        self.Glin = sp.lambdify(list(r1['X']) + list(r1['P']) + [r1['cj']], r1['G'], 'numpy')   # 3 x 24, rung1 state order

    def rhs(self, t, s, cj=1.0):
        return np.concatenate([s[9:36], np.array(self.R(*s, cj), dtype=float)])

    def gauss(self, s, cj=1.0):
        return np.array(self.G(*s, cj), dtype=float)

    def energy(self, s, cj=1.0):
        return float(self.E(*s, cj))


def hessV(q):
    x, y, z = q
    return np.array([[y * y + z * z, 2 * x * y, 2 * x * z], [2 * x * y, z * z + x * x, 2 * y * z], [2 * x * z, 2 * y * z, x * x + y * y]])


def diag_ghost_rate(E):
    """Floquet rate of the colour-diagonal ghost on the DR orbit at E: the Hill system w2 = -(1 + Hess V(q(t))) w,
    6-dim monodromy over one period (as dr_orbit_ghost.py's coarse scan; float)."""
    from scipy.integrate import solve_ivp
    q0, p0, T = base.orbit(E)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        M = s[6:].reshape(6, 6)
        A = np.zeros((6, 6))
        A[:3, 3:] = np.eye(3)
        A[3:, :3] = -(np.eye(3) + hessV(q))
        return np.concatenate([p, -base.gradV(q), (A @ M).ravel()])
    sT = solve_ivp(rhs, (0.0, T), np.concatenate([q0, p0, np.eye(6).ravel()]), **base.RT).y[:, -1]
    mu = np.linalg.eigvals(sT[6:].reshape(6, 6))
    return float(np.log(np.abs(mu).max()) / T)


def seed_diag_ghost(s, delta_diag, rng):
    """add a diagonal ghost of amplitude delta_diag: w = q2 + grad V(q), so q2 += w and q3 += w1, random
    orientation, |w|^2 + |w1|^2 = delta_diag^2 (the (w, W) seeding of hdym_core.full_initial, scaled)."""
    if delta_diag <= 0:
        return s
    v = rng.normal(size=6)
    v *= delta_diag / np.linalg.norm(v)
    s = s.copy()
    for i in range(3):
        s[18 + 4 * i] += v[i]
        s[27 + 4 * i] += v[3 + i]
    return s


def diagonal_state(E):
    """s (36,) on the DR orbit at energy E, ghost-free, no off-diagonal part."""
    q0, p0, T = base.orbit(E)
    g = base.gradV(q0)
    H = hessV(q0)
    s = np.zeros(36)
    for i in range(3):
        s[4 * i] = q0[i]
        s[9 + 4 * i] = p0[i]
        s[18 + 4 * i] = -g[i]
        s[27 + 4 * i] = -(H @ p0)[i]
    return s, T


def od_index():
    """indices into s (36) of the 24 off-diagonal state entries, in rung1 state order (k major, ulist minor)."""
    return np.array([9 * k + m for k in range(4) for m in OD_M])


def gauss_projector(model, s_diag, block):
    """P (24 x 24) projecting an off-diagonal state vector onto ker G_lin, restricted to one block if asked."""
    q = s_diag[[0, 4, 8]]
    p = s_diag[[9, 13, 17]]
    G = np.array(model.Glin(*q, *p, 1.0), dtype=float)      # 3 x 24
    mask = np.ones(24, dtype=bool)
    if block != 'all':
        keep = [j for j, m in enumerate(OD_M) if m in BLOCK_M[block]]
        mask[:] = False
        for k in range(4):
            for j in keep:
                mask[6 * k + j] = True
    Gm = G[:, mask]
    Gm = Gm[np.abs(Gm).sum(axis=1) > 0]                       # rows that touch this block
    P = np.eye(24)
    if len(Gm):
        Pm = np.eye(mask.sum()) - Gm.T @ np.linalg.solve(Gm @ Gm.T, Gm)
        P = np.zeros((24, 24))
        P[np.ix_(mask, mask)] = Pm
    else:
        P = np.diag(mask.astype(float))
    return P


def seeded_state(s_diag, vec24):
    s = s_diag.copy()
    s[od_index()] += vec24
    return s


def pre_escape_drift(ch, dg, lin=1.0):
    """|E(t_k) - E(0)| / |E(0)| at the last checkpoint still in the linear regime (every sector norm and the ghost
    amplitude below `lin`): the integrator's accuracy where it matters for the escape time.  0 if no such checkpoint."""
    E0 = dg[:, 0]
    out = np.zeros(len(E0))
    for n in range(len(E0)):
        okk = ~np.isnan(ch[n, :, 5]) & (np.nanmax(ch[n, :, :5], axis=1) < lin)
        if okk.any():
            k = np.where(okk)[0][-1]
            out[n] = abs(ch[n, k, 5] - E0[n]) / abs(E0[n])
    return out


def driver_sector(ch, tesc=None):
    """who STARTED the runaway: among {block 12, block 13, block 23, diagonal ghost} (chk columns 1..4, all of which
    start at the seed size), the largest at the last checkpoint before escape; returned as 1, 2, 3, or 4
    (= diagonal ghost); -1 for a trajectory that did not escape.  Column 0 (the diagonal positions) is the orbit
    itself, not a deviation, and is not used; `sector` from the kernel says who was largest at |A| = R, a different
    question.  Caveat: linear transients can be large (near a Krein band edge the two block-12 modes are nearly
    parallel and their beat amplifies a 1e-4 seed to 0.2 while Floquet-stable, E = 0.05 in the 19 Oct smoke run),
    so "largest before escape" is the honest reading, not "first past a threshold"."""
    N = ch.shape[0]
    out = np.full(N, -1, dtype=int)
    for n in range(N):
        if tesc is not None and not tesc[n] > 0:
            continue
        valid = ~np.isnan(ch[n, :, 1])
        if valid.any():
            k = np.where(valid)[0][-1]
            out[n] = 1 + int(np.argmax(ch[n, k, 1:5]))
    return out


# ---------------------------------------------------------------- checks
def cross_check(kern, model, E, delta, rng, rtol, atol, hmax):
    from scipy.integrate import solve_ivp
    s_diag, T = diagonal_state(E)
    P = gauss_projector(model, s_diag, 'all')
    v = P @ rng.normal(size=24)
    v *= delta / np.linalg.norm(v)
    s0 = seeded_state(s_diag, v)
    Tc = 2 * T
    sc = kern.integrate_to(s0, 1.0, Tc, rtol, atol, hmax)
    sol = solve_ivp(model.rhs, (0.0, Tc), s0, method='DOP853', rtol=1e-12, atol=1e-14)
    sp_ = sol.y[:, -1]
    rel = np.abs(sc - sp_).max() / np.abs(sp_).max()
    od = od_index()
    rel_od = np.abs(sc[od] - sp_[od]).max() / max(np.abs(sp_[od]).max(), 1e-300)
    print(f"CHECK integrator vs scipy DOP853 at E={E}, t = 2T = {Tc:.2f}: relative difference {rel:.1e} (whole state), "
          f"{rel_od:.1e} (off-diagonal part, |u| = {np.abs(sp_[od]).max():.2e}) (tol 1e-6) {'PASS' if max(rel, rel_od) < 1e-6 else 'FAIL'}", flush=True)
    E0, E1 = model.energy(s0), model.energy(sc)
    G0, G1 = np.linalg.norm(model.gauss(s0)), np.linalg.norm(model.gauss(sc))
    print(f"CHECK energy drift over 2T: {abs(E1 - E0) / abs(E0):.1e} (tol 1e-8) {'PASS' if abs(E1 - E0) / abs(E0) < 1e-8 else 'FAIL'}", flush=True)
    print(f"CHECK Gauss residual of the seed: {G0:.1e} (O(delta^2) = {delta ** 2:.0e}; tol 100 delta^2) "
          f"{'PASS' if G0 < 100 * delta ** 2 else 'FAIL'}; after 2T: {G1:.1e}", flush=True)
    # escape time against scipy, one trajectory with both seeds at an energy where the diagonal ghost is fast
    s_d, _ = diagonal_state(0.549)
    Pd = gauss_projector(model, s_d, 'all')
    v = Pd @ rng.normal(size=24)
    v *= delta / np.linalg.norm(v)
    s1 = seed_diag_ghost(seeded_state(s_d, v), delta, rng)
    te, se, ch, dg, _, st = kern.escape(s1[None, :], 1.0, 600.0, 1e4, rtol, atol, hmax, 24)

    def ev(t, s):
        return 1e4 - np.abs(s[:18]).max()
    ev.terminal = True
    sol = solve_ivp(model.rhs, (0.0, 600.0), s1, method='DOP853', rtol=1e-12, atol=1e-14, events=ev)
    t_sp = sol.t_events[0][0] if len(sol.t_events[0]) else np.nan
    dt_rel = abs(te[0] - t_sp) / t_sp
    print(f"CHECK escape time vs scipy DOP853 (rtol 1e-12) at E=0.549, both seeds {delta:.0e}: kernel {te[0]:.3f}, scipy {t_sp:.3f}, "
          f"relative {dt_rel:.1e} (tol 1e-3) {'PASS' if dt_rel < 1e-3 else 'FAIL'}", flush=True)
    # determinism of the parallel kernel: identical initial states must give identical results.  Not a formality:
    # the 20 Oct orient run on Windows (MinGW gcc, OpenMP, 24 threads) gave TWO escape times across the theta = 0
    # row, whose 128 cells are bit-identical by construction (123 cells at 2965.95, 5 at 1818.84 = one thread's
    # share): numpy/MKL had set flush-to-zero on the main thread, the workers had the default mode, and a denormal
    # in the O(delta^2) z-sector was flushed on one thread only; the chaotic loiter made the last bit into 1100
    # time units.  The kernel now pins the FP control word in every thread (hdym_full_escape.c, fp_mode_set), the
    # INFO line above reports what the threads had before, and this check runs 48 copies so every thread gets some.
    # A short run can still pass with a latent difference (no denormal arises in 300 time units at E = 0.549), so the
    # INFO line is the real evidence; if this FAILS, do not trust any map until it passes.
    Scp = np.repeat(s1[None, :], 48, axis=0)      # more copies than threads, so every thread gets some
    tec, sec, chc, dgc, Sc, stc = kern.escape(Scp, 1.0, 300.0, 1e4, rtol, atol, hmax, 4)
    same = np.all(tec == tec[0]) and np.all(Sc == Sc[0]) and np.all(np.nan_to_num(chc) == np.nan_to_num(chc[0]))
    print(f"CHECK parallel kernel is deterministic (48 identical copies, E=0.549 to t=300): "
          f"t_esc spread {np.ptp(tec):.3e}, final-state spread {np.abs(Sc - Sc[0]).max():.3e} {'PASS' if same else 'FAIL'}", flush=True)
    if not same:
        thr = kern.last_thread_of
        vals, inv = np.unique(np.round(tec, 6), return_inverse=True)
        for k, v in enumerate(vals):
            print(f"      t_esc = {v:.6f}: {np.sum(inv == k)} copies, on threads {sorted(set(thr[inv == k].tolist()))}", flush=True)
    # the diagonal orbit alone must close: ghost-free, no seed
    s_end = kern.integrate_to(s_diag, 1.0, T, rtol, atol, hmax)
    clo = np.abs(s_end[[0, 4, 8, 9, 13, 17]] - s_diag[[0, 4, 8, 9, 13, 17]]).max()
    print(f"CHECK unseeded DR orbit closes after one period in the full system: {clo:.1e} (tol 1e-7) {'PASS' if clo < 1e-7 else 'FAIL'}", flush=True)
    return max(rel, rel_od) < 1e-6 and G0 < 100 * delta ** 2 and clo < 1e-7 and dt_rel < 1e-3 and same


def resume_from(out, a, keys):
    """--resume: load an existing npz with the same run parameters and return it (rows_done tells where to continue);
    None if absent or if any parameter differs.  The seed basis and the fixed diagonal seed are regenerated from
    --seed, so continuing is exact; only the random per-cell seeds of section E are regenerated too (same rng, so the
    rows already done are skipped by drawing and discarding the same numbers)."""
    if not a.resume or not os.path.exists(out):
        return None
    d = np.load(out, allow_pickle=True)
    for k in keys:
        if k not in d.files:
            continue                                   # files from before --resume existed lack a few keys
        if str(d[k]) != str(getattr(a, k)):
            print(f"--resume: parameter {k} differs ({d[k]} in file, {getattr(a, k)} now); starting afresh", flush=True)
            return None
    print(f"--resume: continuing {out} from row {int(d['rows_done'])} of {int(d['rows'])}", flush=True)
    return d


# ---------------------------------------------------------------- sections
def run_E(a, kern, model, rng, out):
    Es = np.logspace(np.log10(a.Emin), np.log10(a.Emax), a.nE)
    N = a.nseed
    tesc = np.full((a.nE, N), np.nan)
    sector = np.full((a.nE, N), -1, dtype=np.int32)
    status = np.full((a.nE, N), -1, dtype=np.int32)
    chk = np.full((a.nE, N, a.nchk, 6), np.nan)
    diag = np.full((a.nE, N, 4), np.nan)
    tcross = np.full((a.nE, N, 5), -1.0)
    Tper = np.zeros(a.nE)
    rate_diag = np.full(a.nE, np.nan)
    start = 0
    old = resume_from(out, a, ('nE', 'nseed', 'Emin', 'Emax', 'delta', 'delta_diag', 'block', 'tmax', 'rtol', 'seed', 'nchk'))
    if old is not None:
        start = int(old['rows_done'])
        for k, arr in (('tesc', tesc), ('sector', sector), ('chk', chk), ('diag', diag), ('T', Tper), ('rate_diag', rate_diag), ('status', status), ('tcross', tcross)):
            if k in old.files:
                arr[:start] = old[k][:start]
    t0 = time.time()
    for i, E in enumerate(Es):
        s_diag, T = diagonal_state(E)
        Tper[i] = T
        P = gauss_projector(model, s_diag, a.block)
        S0 = np.empty((N, 36))
        for n in range(N):
            v = P @ rng.normal(size=24)
            v *= a.delta / np.linalg.norm(v)
            S0[n] = seed_diag_ghost(seeded_state(s_diag, v), a.delta_diag, rng)
        if i < start:
            continue                                       # rows already in the file; the rng has been advanced identically
        te, se, ch, dg, _, st = kern.escape(S0, a.cj, a.tmax, a.R, a.rtol, a.atol, a.hmax, a.nchk)
        tesc[i], sector[i], chk[i], diag[i], status[i] = te, se, ch, dg, st
        tcross[i] = kern.last_tcross
        rate_diag[i] = diag_ghost_rate(E)
        esc = te > 0
        drift = pre_escape_drift(ch, dg)
        secs = ",".join(f"{SECTOR_NAMES[g]}:{(se[esc] == g).sum()}" for g in range(4) if (se[esc] == g).any())
        fail = (te == -2).sum()
        sing = (st == 1).sum()
        print(f"E={E:.5f} T={T:7.3f}  escaped {esc.sum():3d}/{N}{' (%d singular)' % sing if sing else ''}  t_esc median {np.median(te[esc]) if esc.any() else float('nan'):8.1f} "
              f"min {te[esc].min() if esc.any() else float('nan'):8.1f}  sectors [{secs}]  diag-ghost rate {rate_diag[i]:.4f}  max |dE/E| {drift.max():.1e}  "
              f"max |G|_end {dg[:, 3].max():.1e}{'  INTEGRATOR-FAIL x%d' % fail if fail else ''}  ({time.time() - t0:.0f}s)", flush=True)
        np.savez(out, section='E', E=Es, T=Tper, tesc=tesc, sector=sector, status=status, chk=chk, diag=diag, tcross=tcross, rate_diag=rate_diag, rows_done=i + 1, rows=a.nE,
                 nE=a.nE, nseed=a.nseed, Emin=a.Emin, Emax=a.Emax,
                 delta=a.delta, delta_diag=a.delta_diag, R=a.R, tmax=a.tmax, rtol=a.rtol, atol=a.atol, block=a.block, cj=a.cj, nchk=a.nchk, seed=a.seed)
    return tesc, sector


def run_orient(a, kern, model, rng, out):
    s_diag, T = diagonal_state(a.E)
    P = gauss_projector(model, s_diag, a.block)
    basis = []
    for _ in range(3):
        v = P @ rng.normal(size=24)
        for b in basis:
            v -= (v @ b) * b
        basis.append(v / np.linalg.norm(v))
    e1, e2, e3 = basis
    rng_diag = np.random.default_rng(a.seed + 1)
    if a.diag_orient == 'fixed':
        # ONE diagonal-ghost seed for the whole map, so that structure in the off-diagonal seed direction can show;
        # re-randomising it per cell decorrelates the map wherever the diagonal ghost wins the race (seen at E = 0.412,
        # 20 Oct: salt and pepper, 96 % of cells escaping through the diagonal ghost)
        w_fixed = rng_diag.normal(size=6)
        w_fixed *= a.delta_diag / np.linalg.norm(w_fixed)
    th = np.linspace(0, np.pi, a.n)
    ph = np.linspace(0, 2 * np.pi, a.n, endpoint=False)
    tesc = np.full((a.n, a.n), np.nan)
    sector = np.full((a.n, a.n), -1, dtype=np.int32)
    status = np.full((a.n, a.n), -1, dtype=np.int32)
    diag = np.full((a.n, a.n, 4), np.nan)
    chk = np.full((a.n, a.n, a.nchk, 6), np.nan, dtype=np.float32)   # kept since 2026-10-21 (rabbit hole 12); float32 to halve the file
    tcross = np.full((a.n, a.n, 5), -1.0)
    start = 0
    old = resume_from(out, a, ('n', 'E', 'delta', 'delta_diag', 'block', 'tmax', 'rtol', 'seed', 'diag_orient'))
    if old is not None:
        start = int(old['rows_done'])
        for k, arr in (('tesc', tesc), ('sector', sector), ('status', status), ('diag', diag), ('chk', chk), ('tcross', tcross)):
            if k in old.files:
                arr[:start] = old[k][:start]
    t0 = time.time()
    for i, t_ in enumerate(th):
        if i < start:
            continue
        S0 = np.empty((a.n, 36))
        for j, p_ in enumerate(ph):
            v = a.delta * (np.cos(t_) * e1 + np.sin(t_) * np.cos(p_) * e2 + np.sin(t_) * np.sin(p_) * e3)
            if a.diag_orient == 'fixed' and a.delta_diag > 0:
                s_ = seeded_state(s_diag, v)
                for ii in range(3):
                    s_[18 + 4 * ii] += w_fixed[ii]
                    s_[27 + 4 * ii] += w_fixed[3 + ii]
                S0[j] = s_
            else:
                S0[j] = seed_diag_ghost(seeded_state(s_diag, v), a.delta_diag, rng_diag)
        te, se, ch, dg, _, st = kern.escape(S0, a.cj, a.tmax, a.R, a.rtol, a.atol, a.hmax, a.nchk)
        tesc[i], sector[i], diag[i], status[i] = te, se, dg, st
        chk[i] = ch
        tcross[i] = kern.last_tcross
        esc = te > 0
        print(f"theta={t_:.4f} ({i + 1}/{a.n}): escaped {esc.sum()}/{a.n} ({(st == 1).sum()} singular, {(st == 2).sum()} failed), t_esc in [{te[esc].min() if esc.any() else float('nan'):.1f}, "
              f"{te[esc].max() if esc.any() else float('nan'):.1f}], max |dE/E| {pre_escape_drift(ch, dg).max():.1e}  ({time.time() - t0:.0f}s)", flush=True)
        np.savez(out, section='orient', E=a.E, T=T, theta=th, phi=ph, basis=np.array(basis), tesc=tesc, sector=sector, status=status, diag=diag, chk=chk, tcross=tcross, n=a.n,
                 rows_done=i + 1, rows=a.n, delta=a.delta, delta_diag=a.delta_diag, diag_orient=a.diag_orient, R=a.R, tmax=a.tmax, rtol=a.rtol, atol=a.atol, block=a.block, cj=a.cj, seed=a.seed)
    return tesc, sector


def run_phase(a, kern, model, rng, out):
    """--section phase (added 2026-10-21, rabbit hole 1): an ensemble over the PHASE of the seed's projection on a
    Krein quadruplet's unstable subspace at fixed E.  With v the right eigenvector of the block's period map for the
    unstable multiplier mu with Im mu > 0 and w its left eigenvector (w . v = 1, w . conj v = 0), the direction
    d(alpha) = Re(e^{i alpha} v) has projection c = w . d = e^{i alpha} / 2: fixed modulus, phase alpha.  Each cell
    is s = delta * unit(d(alpha) + eta r) with r a random Gauss-projected unit vector (eta = --eta, default 0.3)
    and an independent random diagonal-ghost seed, so the ensemble at each alpha samples the loiter.  The actual
    c = w . s of every cell is stored (its phase differs from alpha by O(eta)).  For a parametric (real) pair the
    construction degenerates to the sign of alpha and the run is refused.  Requires --block."""
    if a.block == 'all':
        raise SystemExit("--section phase needs --block 12 (or 13, 23)")
    s_diag, T = diagonal_state(a.E)
    P = gauss_projector(model, s_diag, a.block)
    Kf = base.load_K()
    _, clo, Ms, _ = base.monodromies(Kf, a.E, a.cj)
    M = Ms[a.block]
    mu, V = np.linalg.eig(M)
    W = np.linalg.inv(V)
    ghost = np.argsort(np.abs(mu - 1.0))[4:]                   # drop the zero-mode sector (a defective 1, split by ~1e-5 in floating point), as B001 does
    unst = np.array([k for k in ghost if abs(mu[k]) > 1 + 1e-6])
    if len(unst) != 2 or np.abs(mu[unst].imag).max() < 1e-6:
        raise SystemExit(f"--section phase: block {a.block} at E={a.E} has ghost multipliers {np.round(mu[ghost], 4)}: not a Krein quadruplet")
    k = unst[np.argmax(mu[unst].imag)]
    v8, w8 = V[:, k], W[k]
    idx = base.block_idx(a.block)
    v24 = np.zeros(24, dtype=complex)
    v24[idx] = v8
    w24 = np.zeros(24, dtype=complex)
    w24[idx] = w8
    gres = np.abs(np.array(model.Glin(*s_diag[[0, 4, 8]], *s_diag[[9, 13, 17]], 1.0), dtype=float) @ v24).max() / np.abs(v24).max()
    print(f"phase section: block {a.block}, E={a.E}, T={T:.3f}, Krein multiplier {np.round(mu[k], 5)} (rate {np.log(abs(mu[k])) / T:.4f}); "
          f"eigenvector Gauss residual {gres:.1e}; orbit closure {clo:.1e}; eta={a.eta}", flush=True)
    alphas = np.linspace(0, 2 * np.pi, a.nalpha, endpoint=False)
    N = a.nseed
    tesc = np.full((a.nalpha, N), np.nan)
    sector = np.full((a.nalpha, N), -1, dtype=np.int32)
    status = np.full((a.nalpha, N), -1, dtype=np.int32)
    diag = np.full((a.nalpha, N, 4), np.nan)
    chk = np.full((a.nalpha, N, a.nchk, 6), np.nan, dtype=np.float32)
    tcross = np.full((a.nalpha, N, 5), -1.0)
    cproj = np.full((a.nalpha, N), np.nan, dtype=complex)
    start = 0
    old = resume_from(out, a, ('nalpha', 'nseed', 'E', 'eta', 'delta', 'delta_diag', 'block', 'tmax', 'rtol', 'seed'))
    if old is not None:
        start = int(old['rows_done'])
        for kk, arr in (('tesc', tesc), ('sector', sector), ('status', status), ('diag', diag), ('chk', chk), ('tcross', tcross), ('cproj', cproj)):
            arr[:start] = old[kk][:start]
    t0 = time.time()
    for i, al in enumerate(alphas):
        d24 = P @ np.real(np.exp(1j * al) * v24)
        S0 = np.empty((N, 36))
        for n in range(N):
            r = P @ rng.normal(size=24)
            r /= np.linalg.norm(r)
            v = d24 / np.linalg.norm(d24) + a.eta * r
            v *= a.delta / np.linalg.norm(v)
            cproj[i, n] = w24 @ v / a.delta
            S0[n] = seed_diag_ghost(seeded_state(s_diag, v), a.delta_diag, rng)
        if i < start:
            continue
        te, se, ch, dg, _, st = kern.escape(S0, a.cj, a.tmax, a.R, a.rtol, a.atol, a.hmax, a.nchk)
        tesc[i], sector[i], diag[i], status[i], chk[i], tcross[i] = te, se, dg, st, ch, kern.last_tcross
        esc = te > 0
        print(f"alpha={al:.4f} ({i + 1}/{a.nalpha}): escaped {esc.sum()}/{N} ({(st == 1).sum()} singular, {(st == 2).sum()} failed), t_esc median {np.median(te[esc]) if esc.any() else float('nan'):.1f} "
              f"[{te[esc].min() if esc.any() else float('nan'):.1f}, {te[esc].max() if esc.any() else float('nan'):.1f}], |c| {np.abs(cproj[i]).mean():.3f}+-{np.abs(cproj[i]).std():.3f}, "
              f"arg c - alpha sd {np.std(np.angle(cproj[i] * np.exp(-1j * al))):.3f}, max |dE/E| {pre_escape_drift(ch, dg).max():.1e}  ({time.time() - t0:.0f}s)", flush=True)
        np.savez(out, section='phase', E=a.E, T=T, alpha=alphas, mu=mu[k], v24=v24, w24=w24, cproj=cproj, tesc=tesc, sector=sector, status=status, diag=diag, chk=chk, tcross=tcross,
                 rows_done=i + 1, rows=a.nalpha, nalpha=a.nalpha, nseed=a.nseed, eta=a.eta, delta=a.delta, delta_diag=a.delta_diag, R=a.R, tmax=a.tmax, rtol=a.rtol, atol=a.atol,
                 block=a.block, cj=a.cj, nchk=a.nchk, seed=a.seed)
    return tesc, sector


# ---------------------------------------------------------------- plots (from an npz; no run)
def linear_rates(results_dir):
    """rate per (E, block) from dr_offdiag_floquet_typed.txt, if present: {block: (E array, rate array, type array)}."""
    import re
    f = os.path.join(results_dir, 'dr_offdiag_floquet_typed.txt')
    if not os.path.exists(f):
        return None
    out = {b: [] for b in base.BLOCKS}
    E = T = None
    fine = False
    for l in open(f):
        if l.startswith('== Fine'):
            fine = True
        if not fine:
            continue
        m = re.match(r"E=([0-9.]+) T=\s*([0-9.]+)", l)
        if m:
            E, T = float(m.group(1)), float(m.group(2))
            continue
        m = re.match(r"\s+blk(\d\d)\s+([0-9.]+)\s+(\S+)", l)
        if m and E is not None:
            out[m.group(1)].append((E, np.log(max(float(m.group(2)), 1.0)) / T, m.group(3)))
    return {b: np.array(v, dtype=object) for b, v in out.items()}


def plot(tag, results_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = np.load(os.path.join(results_dir, tag + '.npz'), allow_pickle=True)
    sec = str(d['section'])
    cols = {0: '#52514e', 1: '#2a78d6', 2: '#eb6834', 3: '#1baf7a'}
    if sec == 'E':
        Es, tesc, sector = d['E'], d['tesc'], d['sector']
        nd = int(d['rows_done'])
        fig, ax = plt.subplots(figsize=(10, 5.5), facecolor='#fcfcfb')
        ax.set_facecolor('#fcfcfb')
        chk = d['chk']
        cols5 = dict(cols)
        cols5[4] = '#9085e9'
        names5 = SECTOR_NAMES + ['diagonal ghost']
        for i in range(nd):
            te = tesc[i]
            drv = driver_sector(chk[i], te)
            for g in range(5):
                m = (drv == g) & (te > 0)
                if m.any():
                    ax.scatter(np.full(m.sum(), Es[i]), te[m], s=9, color=cols5[g], alpha=0.7, lw=0)
            m = te == -1
            if m.any():
                ax.scatter([Es[i]], [float(d['tmax']) * 1.15], s=14, marker='^', color='#c3c2b7', lw=0)
        lr = linear_rates(results_dir)
        if lr:
            ln = np.log(float(d['R']) / float(d['delta']))
            for b, c in (('12', 1), ('13', 2), ('23', 3)):
                arr = lr[b]
                E_ = np.array([r[0] for r in arr], dtype=float)
                rate = np.array([r[1] for r in arr], dtype=float)
                ok = rate > 0
                tmax_ = float(d['tmax'])
                for ln_, ms_, al_ in ((ln, 2.5, 0.5), (np.log(1.0 / float(d['delta'])), 1.2, 0.35)):
                    tp = ln_ / rate[ok]
                    below = tp <= tmax_
                    ax.plot(E_[ok][below], tp[below], '.', ms=ms_, color=cols[c], alpha=al_)
                    # predictions beyond tmax: the run could not confirm them; drawn hollow and faint
                    ax.plot(E_[ok][~below], tp[~below], 'o', ms=ms_ + 1, mfc='none', mec=cols[c], mew=0.5, alpha=0.35)
            ax.plot([], [], '.', color='#0b0b0b', label='linear: ln(R/δ)/rate and ln(1/δ)/rate per block (typed log); hollow = beyond tmax, untestable here')
        ax.axhline(float(d['tmax']), color='#52514e', lw=0.8, ls='--')
        ax.text(float(d['E'][0]) * 1.02, float(d['tmax']) * 1.08, 'tmax', fontsize=8, color='#52514e')
        if 'rate_diag' in d.files:
            rd = d['rate_diag'][:nd]
            dd_ = max(float(d['delta_diag']), 1e-10)
            # NaN where the diagonal ghost is Floquet-stable, so the curve BREAKS there instead of bridging the gap
            # (an earlier version joined the neighbours and drew a spurious prediction across E = 0.27-0.31)
            rsafe = np.where(rd > 1e-6, rd, np.nan)
            ax.plot(Es[:nd], np.log(float(d['R']) / dd_) / rsafe, '-', color=cols5[4], lw=1,
                    label='linear, diagonal ghost: ln(R/δ_diag)/rate_diag (and ln(1/δ_diag), dotted); absent where Floquet-stable')
            ax.plot(Es[:nd], np.log(1.0 / dd_) / rsafe, ':', color=cols5[4], lw=1)
        for g in (1, 2, 3, 4):
            ax.scatter([], [], color=cols5[g], label=f'runaway led by {names5[g]} (largest before escape)')
        ax.scatter([], [], marker='^', color='#c3c2b7', label=f'not escaped by tmax = {float(d["tmax"]):.0f} (drawn at 1.15 tmax)')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('E')
        ax.set_ylabel('escape time')
        ax.set_title(f'Nonlinear escape times, off-diagonal seed δ = {float(d["delta"]):.0e} on the DR orbit, block {str(d["block"])} '
                     f'({nd}/{int(d["rows"])} rows)', fontsize=10, loc='left')
        ax.legend(frameon=False, fontsize=8, loc='lower left')
        ax.grid(True, color='#eeede9')
        fig.savefig(os.path.join(results_dir, tag + '.png'), dpi=150, bbox_inches='tight')
        print(f"wrote {tag}.png")
    else:
        te, se = d['tesc'], d['sector']
        fig, axs = plt.subplots(1, 2, figsize=(12, 5), facecolor='#fcfcfb')
        im = axs[0].imshow(np.where(te > 0, np.log10(np.where(te > 0, te, 1)), np.nan), origin='lower', aspect='auto',
                           extent=[0, 2, 0, 1], cmap='viridis')
        plt.colorbar(im, ax=axs[0], label='log10 t_esc')
        axs[0].set_xlabel('φ/π')
        axs[0].set_ylabel('θ/π')
        axs[0].set_title(f'escape time, E = {float(d["E"])}, δ = {float(d["delta"]):.0e}', fontsize=10)
        cm = matplotlib.colors.ListedColormap(['#c3c2b7', cols[0], cols[1], cols[2], cols[3]])
        axs[1].imshow(se + 1, origin='lower', aspect='auto', extent=[0, 2, 0, 1], cmap=cm, vmin=0, vmax=4)
        axs[1].set_xlabel('φ/π')
        axs[1].set_title('escaping sector (light grey: not escaped; dark grey: diagonal; blue/orange/green: blocks 12/13/23)', fontsize=9)
        fig.savefig(os.path.join(results_dir, tag + '.png'), dpi=150, bbox_inches='tight')
        print(f"wrote {tag}.png")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--section', choices=['E', 'orient', 'phase'], default='E')
    ap.add_argument('--plot', default=None, help='tag of an existing npz to draw; no run')
    ap.add_argument('--tag', default='offdiag_escape')
    ap.add_argument('--results', default=os.path.join(HERE, '..', 'results'))
    ap.add_argument('--Emin', type=float, default=0.05)
    ap.add_argument('--Emax', type=float, default=1.0)
    ap.add_argument('--nE', type=int, default=100)
    ap.add_argument('--nseed', type=int, default=32)
    ap.add_argument('--E', type=float, default=0.412)
    ap.add_argument('--n', type=int, default=96)
    ap.add_argument('--block', choices=['all', '12', '13', '23'], default='all')
    ap.add_argument('--delta', type=float, default=1e-4)
    ap.add_argument('--delta-diag', type=float, default=1e-4, help='diagonal ghost seed amplitude; 0 = none (then truncation error seeds it at ~rtol)')
    ap.add_argument('--R', type=float, default=1e4)
    ap.add_argument('--tmax', type=float, default=3000.0)
    ap.add_argument('--rtol', type=float, default=1e-11)
    ap.add_argument('--atol', type=float, default=1e-13)
    ap.add_argument('--hmax', type=float, default=0.25)
    ap.add_argument('--nchk', type=int, default=60)
    ap.add_argument('--cj', type=float, default=1.0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--resume', action='store_true', help='continue an existing <tag>.npz with identical parameters from its rows_done')
    ap.add_argument('--diag-orient', choices=['fixed', 'random'], default='fixed', help='orient section: one diagonal-ghost seed for the whole map (default) or one per cell')
    ap.add_argument('--nalpha', type=int, default=12, help='phase section: number of phases over [0, 2 pi)')
    ap.add_argument('--eta', type=float, default=0.3, help='phase section: random admixture to the eigen-direction')
    a = ap.parse_args()
    results = os.path.normpath(a.results)
    os.makedirs(results, exist_ok=True)
    if a.plot:
        plot(a.plot, results)
        return
    t0 = time.time()
    print(f"offdiag_escape.py  python {sys.version.split()[0]}  numpy {np.__version__}  sympy {sp.__version__}  ({sys.executable})")
    print("run parameters: " + " ".join(f"{k}={v}" for k, v in vars(a).items() if k not in ('plot', 'results')))
    print("orbit: DR, X0=%r PX0=%r T_HALF=%r (E=1/2, scaled); diagonal ghost not seeded; seed Gauss-projected (linear G)" % (base.X0, base.PX0, base.T_HALF))
    kern = Kernel()
    model = Model()
    rng = np.random.default_rng(a.seed)
    ok = cross_check(kern, model, a.E if a.section in ('orient', 'phase') else 0.412, a.delta, rng, a.rtol, a.atol, a.hmax)
    out = os.path.join(results, a.tag + '.npz')
    if a.section == 'E':
        tesc, sector = run_E(a, kern, model, rng, out)
    elif a.section == 'phase':
        tesc, sector = run_phase(a, kern, model, rng, out)
    else:
        tesc, sector = run_orient(a, kern, model, rng, out)
    fails = int((tesc == -2).sum())
    print(f"\nCHECK no integrator failures: {fails} {'PASS' if fails == 0 else 'FAIL'}")
    print(f"{'ALL CHECKS PASS' if ok and fails == 0 else 'SOME CHECK FAILED'}")
    print(f"DONE {time.time() - t0:.0f}")


if __name__ == "__main__":
    main()
