"""rh7_selection_rule.py -- rabbit hole 7: the selection rule of block 12 (folds only at M = +1, Krein
collisions only at M = -1) derived from the half-period symmetry, with the numerical checks of each step.

The derivation (notes/RH7_selection_rule.md) has three steps; this script checks each on the DR orbit.

  1. PARITY.  In the (u12, u21) basis the half-period symmetry S~ K(t + T/2) S~ = K(t), S~ = diag(+1, -1) on each
     derivative order, says the self-coupling blocks K_aa are T/2-periodic (even harmonics of Omega = 2 pi / T only)
     and the cross blocks K_ab flip sign (odd harmonics only).  In the sector basis u+- = (u12 +- u21)/sqrt 2 the
     same symmetry is the SWAP: K++(t + T/2) = K--(t), K+-(t + T/2) = K-+(t).  CHECK: both, pointwise along the orbit,
     and the harmonic content by FFT.
  2. SECTOR DEGENERACY.  Zero the cross blocks K+-, K-+ (the "decoupled symmetric model"): the two sectors are then
     half-period translates of each other, their period-T monodromies Phi_- Phi_+ and Phi_+ Phi_- are conjugate, so
     every multiplier of M_dec is double and N_dec = swap * Phi_T/2 has exactly +-sqrt(mu).  CHECK at every E.
     In this model a single-sector ghost pair at e^{+-i theta} gives N-eigenvalues at +-theta/2 and pi -+ theta/2,
     so the only coincidences of the two copies are at N = +-1 (theta = 0, the sector's +1 fold) and at N = +-i
     (theta = pi, the sector's period-doubling fold); at N = +-i the two coincident members carry OPPOSITE Krein
     sign (one descends from e^{i theta}, the other from e^{-i theta}).  Hence, in the decoupled model, parametric
     bands = the sector's +1 tongues, Krein bands = the sector's -1 tongues; the true cross-coupling (c_j = 1)
     splits the double real pair at -1 into a genuine complex quadruplet.  SCAN: the sector's own tongues over E,
     typed +1 / -1, against B001's block-12 intervals (results/dr_offdiag_floquet_typed.txt).
  3. LEMMA (exact, no model): M = N^2 with N real symplectic.  A real negative multiplier of M off the circle has
     square roots +-i r, non-real with |lambda| != 1, hence a full quadruplet {i r, -i r, i/r, -i/r} of N, whose
     squares are {-r^2, -r^2, -1/r^2, -1/r^2}: it is at least double.  So one pair can never leave through -1
     (B005's "a + 2 >= 0"), and a -1 instability is necessarily two pairs (a Krein quadruplet or, in the decoupled
     model, the degenerate limit of one).  CHECK on B005's log: every off-circle event is either a +-1 fold or a
     four-member collision, never a lone pair at -1.

Cost: step 2's scan is nE half-period integrations, ~1.7 s each on one core.  Default nE = 150: ~5 min.
Usage (from the repo root):
    python rabbit_holes/ym/code/rh7_selection_rule.py [--nE 150] > rabbit_holes/ym/results/rh7_selection_rule.txt
Writes results/rh7_sector_scan.npz (the decoupled-sector multipliers vs E) and results/rh7_sector_tongues.png.
Float, DOP853 rtol 1e-12 (dr_offdiag_floquet.RT).  Added 2026-10-21 (rabbit-hole session).
"""
import argparse
import os
import re
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
B001 = os.path.join(ROOT, 'islands', 'results', 'dr_offdiag_floquet_typed.txt')
B005 = os.path.join(ROOT, 'islands', 'results', 'krein_b005.txt')

# block-12 state order from base.block_idx('12'): (u12, u21, u12', u21', u12'', u21'', u12''', u21''')
IDX = base.block_idx('12')
STILDE = np.diag(np.array([1.0, -1.0] * 4))
# sector basis u+- = (u12 +- u21)/sqrt2 at every derivative order; R maps (u12,u21,...) -> (u+,u-,...)
R = np.kron(np.eye(4), np.array([[1, 1], [1, -1]]) / np.sqrt(2))
SWAP = np.kron(np.eye(4), np.array([[0, 1], [1, 0]]))
PLUS = [0, 2, 4, 6]       # sector-basis rows of u+ and its derivatives
MINUS = [1, 3, 5, 7]


def K12(Kf, q, p, cj=1.0):
    return np.array(Kf(*q, *p, cj), dtype=float)[np.ix_(IDX, IDX)]


def K_sector(Kf, q, p, cj=1.0):
    """K in the sector basis (u+, u-, u+', u-', ...)."""
    return R @ K12(Kf, q, p, cj) @ R.T


def K_decoupled(Kf, q, p):
    K = K_sector(Kf, q, p)
    Kd = K.copy()
    Kd[np.ix_(PLUS, MINUS)] = 0.0
    Kd[np.ix_(MINUS, PLUS)] = 0.0
    return Kd


def orbit_samples(E, n):
    q0, p0, T = base.orbit(E)
    sol = solve_ivp(lambda t, s: np.concatenate([s[3:], -base.gradV(s[:3])]), (0.0, T),
                    np.concatenate([q0, p0]), t_eval=np.linspace(0, T, n, endpoint=False), **base.RT)
    return sol.t, sol.y[:3].T, sol.y[3:].T, T


# ---------------------------------------------------------------- step 1: parity
def step1_parity(Kf, E, n=256):
    t, q, p, T = orbit_samples(E, n)
    Ks = np.array([K12(Kf, q[i], p[i]) for i in range(n)])          # (n, 8, 8)
    Ksec = np.array([K_sector(Kf, q[i], p[i]) for i in range(n)])
    h = n // 2
    scale = np.abs(Ks).max()
    # pointwise: S~ K(t+T/2) S~ = K(t)
    r_par = max(np.abs(STILDE @ Ks[(i + h) % n] @ STILDE - Ks[i]).max() for i in range(n)) / scale
    # sector swap: K++(t+T/2) = K--(t), K+-(t+T/2) = K-+(t)
    r_swap = max(np.abs(SWAP @ Ksec[(i + h) % n] @ SWAP - Ksec[i]).max() for i in range(n)) / scale
    # harmonic content in the (u12,u21) basis: self entries (a,a) and cross entries (a,b), rows 6,7 (the u'''' rows)
    F = np.fft.rfft(Ks, axis=0) / n                                  # (n/2+1, 8, 8), harmonic index k of Omega
    self_mask = np.zeros((8, 8), bool)
    cross_mask = np.zeros((8, 8), bool)
    for i in range(8):
        for j in range(8):
            (self_mask if (i % 2) == (j % 2) else cross_mask)[i, j] = True
    amp_self = np.abs(F[:, self_mask]).max(axis=1)
    amp_cross = np.abs(F[:, cross_mask]).max(axis=1)
    odd = np.arange(F.shape[0]) % 2 == 1
    even = ~odd
    return dict(r_par=r_par, r_swap=r_swap,
                self_odd=amp_self[odd].max() / amp_self.max(), self_even=amp_self[even].max() / amp_self.max(),
                cross_even=amp_cross[even].max() / amp_cross.max(), cross_odd=amp_cross[odd].max() / amp_cross.max(),
                amp_self=amp_self, amp_cross=amp_cross, T=T)


# ---------------------------------------------------------------- step 2: decoupled symmetric model
def decoupled_maps(Kf, E):
    """Phi_T/2 and Phi_T of the decoupled sector model (16+6 ODE), plus the true block maps for comparison."""
    q0, p0, T = base.orbit(E)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        Kd = K_decoupled(Kf, q, p)
        Kt = K_sector(Kf, q, p)
        return np.concatenate([p, -base.gradV(q), (Kd @ s[6:70].reshape(8, 8)).ravel(), (Kt @ s[70:].reshape(8, 8)).ravel()])
    s0 = np.concatenate([q0, p0, np.eye(8).ravel(), np.eye(8).ravel()])
    sol = solve_ivp(rhs, (0.0, T), s0, t_eval=[T / 2, T], **base.RT)
    sh, sT = sol.y[:, 0], sol.y[:, 1]
    Phd, Md = sh[6:70].reshape(8, 8), sT[6:70].reshape(8, 8)
    Pht, Mt = sh[70:].reshape(8, 8), sT[70:].reshape(8, 8)
    return Phd, Md, Pht, Mt, T


def sector_ghost_pair(Md):
    """The decoupled model's M is block-diagonal in the sector basis; each 4x4 sector has a zero-mode pair near 1
    and one ghost pair.  Return the + sector's ghost pair (the - sector's is the same set, step 2) and the
    conjugacy residual between the sectors."""
    Mp = Md[np.ix_(PLUS, PLUS)]
    Mm = Md[np.ix_(MINUS, MINUS)]
    ep, em = np.linalg.eigvals(Mp), np.linalg.eigvals(Mm)
    conj = np.abs(np.sort_complex(ep) - np.sort_complex(em)).max()
    keep = np.argsort(np.abs(ep - 1.0))[2:]                           # drop the 2 nearest 1 (zero-mode pair)
    g = ep[keep]
    return g[np.argsort(-np.abs(g))], conj, np.abs(Md[np.ix_(PLUS, MINUS)]).max() / np.abs(Md).max()


def classify(g, tol=1e-6):
    """type of the sector's ghost pair: 'circle', 'fold+1' (real positive, off circle), 'fold-1' (real negative)."""
    m = g[0]
    if abs(abs(m) - 1) < tol:
        return 'circle', float(np.angle(m) / np.pi)
    if abs(m.imag) < 1e-6 * abs(m):
        return ('fold+1' if m.real > 0 else 'fold-1'), float(np.log(abs(m)))
    return 'cplx', float(np.log(abs(m)))


def b001_intervals():
    """block-12 unstable intervals and types from the typed log."""
    txt = open(B001).read().split('== Unstable intervals')[1]
    out = []
    cur = None
    for line in txt.splitlines():
        m = re.match(r'\s*interval \d+: E in \[([\d.]+), ([\d.]+)\]', line)
        if m:
            cur = (float(m.group(1)), float(m.group(2)))
            continue
        m = re.match(r'\s*blk12: max \|mu\| = ([\d.]+) at E=([\d.]+); type sequence (.*?);', line)
        if m and cur:
            typ = 'CPLX' if 'CPLX' in m.group(3) else 'REAL'
            out.append((cur[0], cur[1], typ, float(m.group(1))))
    return out


def b005_events():
    """every LEAVE/RE-ENTER line of B005: kind (FOLD/COLLISION) and the angle for folds."""
    out = []
    for line in open(B005):
        m = re.match(r'\s*E=([\d.]+) (LEAVE|RE-ENTER): (.*?)(FOLD|COLLISION)', line)
        if m:
            ang = re.search(r'meet at arg ([-+\d.]+)pi', line)
            out.append((float(m.group(1)), m.group(2), m.group(4), float(ang.group(1)) if ang else np.nan))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nE', type=int, default=150)
    ap.add_argument('--emin', type=float, default=0.05)
    ap.add_argument('--emax', type=float, default=1.0)
    a = ap.parse_args()
    t0 = time.time()
    print(f"rh7_selection_rule.py  python {sys.version.split()[0]}  numpy {np.__version__}")
    print(f"run parameters: nE={a.nE} log-spaced in [{a.emin}, {a.emax}], RT={base.RT}, c_j = 1, block 12")
    Kf = base.load_K()

    print("\n== step 1: parity of K under the half-period reflection ==")
    for E in (0.1, 0.3, 0.6):
        r = step1_parity(Kf, E)
        ok = r['r_par'] < 1e-9 and r['r_swap'] < 1e-9
        print(f"E={E:.2f}: S~ K(t+T/2) S~ = K(t) to {r['r_par']:.1e}; sector swap K++(t+T/2) = K--(t), K+-(t+T/2) = K-+(t) to {r['r_swap']:.1e}  (tol 1e-9) {'PASS' if ok else 'FAIL'}")
        print(f"        harmonics of Omega = 2pi/T in the (u12,u21) basis, relative to the largest: self blocks odd {r['self_odd']:.1e} / even {r['self_even']:.1e};"
              f" cross blocks even {r['cross_even']:.1e} / odd {r['cross_odd']:.1e}")
        ks = np.argsort(-r['amp_self'])[:4]
        kc = np.argsort(-r['amp_cross'])[:4]
        print(f"        leading harmonics: self k = {list(ks)} (amp {np.round(r['amp_self'][ks], 3)}), cross k = {list(kc)} (amp {np.round(r['amp_cross'][kc], 3)})")

    print("\n== step 2: the decoupled symmetric model (sector cross blocks zeroed), scanned in E ==")
    print("columns: E | T | sector ghost pair of M_dec: type, angle/pi or rate | sector conjugacy | N_dec = swap Phi_T/2: +-sqrt(mu) residual | true block-12 ghost pair type")
    grid = np.logspace(np.log10(a.emin), np.log10(a.emax), a.nE)
    rows = []
    worst = dict(conj=0.0, sqrt=0.0, N2=0.0, offblock=0.0)
    for E in grid:
        Phd, Md, Pht, Mt, T = decoupled_maps(Kf, E)
        g, conj, offb = sector_ghost_pair(Md)
        typ, val = classify(g)
        Nd = SWAP @ Phd
        lam = np.linalg.eigvals(Nd)
        # +-sqrt(mu): the spectrum of N_dec must be symmetric under lambda -> -lambda
        sq = max(min(abs(l + lam)) for l in lam) / max(abs(lam))
        n2 = np.abs(Nd @ Nd - Md).max() / np.abs(Md).max()
        worst['conj'] = max(worst['conj'], conj)
        worst['sqrt'] = max(worst['sqrt'], sq)
        worst['N2'] = max(worst['N2'], n2)
        worst['offblock'] = max(worst['offblock'], offb)
        gt = base.ghost_multipliers(Mt)
        tt = 'circle' if abs(abs(gt[0]) - 1) < 1e-6 else ('CPLX' if abs(gt[0].imag) > 1e-6 * abs(gt[0]) else ('REAL+' if gt[0].real > 0 else 'REAL-'))
        rows.append((E, T, typ, val, abs(g[0]), conj, sq, tt, abs(gt[0]), g[0], gt[0]))
        print(f"E={E:.5f} T={T:8.4f} | sector pair {typ:7s} {val:+8.4f} |mu|={abs(g[0]):.4f} | conj {conj:.1e} | +-sqrt {sq:.1e} N^2 {n2:.1e} | true {tt:6s} |mu|={abs(gt[0]):.4f}")
    print(f"CHECK the two sectors of the decoupled model have the same period-T spectrum, worst {worst['conj']:.1e} (tol 1e-7) {'PASS' if worst['conj'] < 1e-7 else 'FAIL'}")
    print(f"CHECK M_dec is block-diagonal in the sector basis, worst off-block {worst['offblock']:.1e} (tol 1e-10) {'PASS' if worst['offblock'] < 1e-10 else 'FAIL'}")
    print(f"CHECK spec(N_dec) = +-sqrt(spec(M_dec)) (symmetric under lambda -> -lambda), worst {worst['sqrt']:.1e} (tol 1e-6) {'PASS' if worst['sqrt'] < 1e-6 else 'FAIL'}")
    print(f"CHECK N_dec^2 = M_dec, worst {worst['N2']:.1e} (tol 1e-8) {'PASS' if worst['N2'] < 1e-8 else 'FAIL'}")

    # tongues of the decoupled sector vs B001's block-12 intervals
    print("\n== the decoupled sector's tongues against B001's block-12 intervals ==")
    typs = np.array([r[2] for r in rows])
    Es = np.array([r[0] for r in rows])
    tongues = []
    i = 0
    while i < len(rows):
        if typs[i] in ('fold+1', 'fold-1'):
            j = i
            while j + 1 < len(rows) and typs[j + 1] == typs[i]:
                j += 1
            rate = max(r[3] for r in rows[i:j + 1])
            tongues.append((Es[i], Es[j], typs[i], rate))
            i = j + 1
        else:
            i += 1
    print("  decoupled sector (prediction):")
    for lo, hi, ty, rate in tongues:
        print(f"    E in [{lo:.4f}, {hi:.4f}]  {ty}  max rate per period {rate:.3f}   -> predicts a {'PARAMETRIC (REAL+)' if ty == 'fold+1' else 'KREIN (CPLX)'} band of the true block")
    print("  true block 12 (B001, measured):")
    for lo, hi, ty, mu in b001_intervals():
        print(f"    E in [{lo:.4f}, {hi:.4f}]  {ty}  max |mu| {mu:.4f}")
    print("  (the decoupled model has no complex quadruplets by construction: every instability is a double real pair; "
          "the true cross-coupling turns the -1 ones into quadruplets and splits the +1 ones between the two pairs)")

    print("\n== step 3: the exact lemma against B005's events ==")
    lone_minus = 0
    for E, what, kind, ang in b005_events():
        tag = ''
        if kind == 'FOLD':
            near = min(abs(ang), abs(abs(ang) - 1.0))
            tag = f"fold of N at arg {ang:+.3f}pi -> M at +1 (distance of |arg| from {{0, 1}}: {near:.3f}pi)"
            if 0.25 < abs(ang) < 0.75:
                lone_minus += 1
        else:
            tag = "four members leave together: two-pair collision (M at -1 by B005's angles +-0.44..0.59 pi of N)"
        print(f"  E={E:.5f} {what:8s} {tag}")
    print(f"CHECK no single-pair fold of N near +-i (which would be a lone pair of M leaving at -1, forbidden by the lemma): {lone_minus} found {'PASS' if lone_minus == 0 else 'FAIL'}")

    np.savez(os.path.join(RES, 'rh7_sector_scan.npz'), E=Es, T=np.array([r[1] for r in rows]),
             sector_mu=np.array([r[9] for r in rows]), true_mu=np.array([r[10] for r in rows]),
             sector_type=typs, true_type=np.array([r[7] for r in rows]))

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
        mu_s = np.array([r[9] for r in rows])
        mu_t = np.array([r[10] for r in rows])
        for lo, hi, ty, mu in b001_intervals():
            for axx in ax:
                axx.axvspan(lo, hi, color=('tab:red' if ty == 'CPLX' else 'tab:blue'), alpha=0.15, lw=0)
        ax[0].plot(Es, np.abs(mu_s), 'k.-', ms=3, lw=0.8, label='decoupled sector |mu|')
        ax[0].plot(Es, np.abs(mu_t), 'g.', ms=3, label='true block 12 |mu| (leading ghost)')
        ax[0].set_ylabel('|mu| over one period')
        ax[0].legend(loc='upper left', fontsize=8)
        ang = np.angle(mu_s) / np.pi
        col = np.where(typs == 'fold-1', 'tab:red', np.where(typs == 'fold+1', 'tab:blue', 'k'))
        ax[1].scatter(Es, ang, c=col, s=8)
        ax[1].set_ylabel('arg(mu_sector)/pi')
        ax[1].set_xlabel('E')
        ax[1].set_xscale('log')
        ax[0].set_title('rabbit hole 7: decoupled sector tongues (points: red = fold at -1, blue = fold at +1) vs B001 block-12 bands (shading: red Krein, blue parametric)', fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, 'rh7_sector_tongues.png'), dpi=130)
        print("figure: results/rh7_sector_tongues.png")
    except Exception as e:                                           # noqa: BLE001
        print(f"(no figure: {e})")
    print(f"DONE {time.time() - t0:.0f}")


if __name__ == '__main__':
    main()
