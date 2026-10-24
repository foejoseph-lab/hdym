"""krein_b005.py -- B005: the Krein signature of the two ghost-sector pairs of block 12 on the DR orbit, measured.

Everything before this measured |mu| and arg mu and INFERRED the Krein signs from whether pairs pass through each
other or collide (WORKORDER v3, response to the 19 Oct review).  Here the sign is computed.

1. The symplectic form.  L2 (ym/rung1_system.pkl) is the quadratic Lagrangian of the six off-diagonal components,
   L2(u, u', u''; q, p), with u''^2 coefficient -1/2 and no u'' cross terms.  Ostrogradsky: Q1 = u, Q2 = u',
   P2 = dL2/du'', P1 = dL2/du' - d/dt P2, with d/dt acting on u (u -> u', ..., u'' -> u''') and on the background
   (qdot = p, pdot = -grad V, ghost-free).  All four are linear in the state x = (u, u', u'', u''') with
   background-dependent coefficients, (Q1, Q2, P1, P2) = J(q, p) x, and the physical symplectic form in the state
   coordinates is Omega(q, p) = J^T J_can J, J_can the canonical form in (Q, P) ordering.  Restricted to a block it
   is an 8 x 8 antisymmetric matrix.  CHECK: M^T Omega(0) M = Omega(0) for the period map of every block (the form
   is periodic with the background) -- this is the test that Omega is the right form.
2. The half-period map.  q(t + T/2) = S q(t), S = diag(1, -1, 1) on the planar DR orbit, and the block-12 state
   transforms with S~ = diag(+1, -1) on (u12, u21) at every derivative order; N = S~ Phi_{T/2}, CHECK N^2 = M and
   N^T Omega N = Omega.  N's angles are half of M's: M's -1 (period doubling) is N's +-i, which is not a fold, so a
   single pair of M cannot go real at -1 (the tangency of a = mu + 1/mu at -2 seen in B001, never below).
3. Krein sign.  For an eigenvalue of N on the unit circle with eigenvector v, sigma = sign(i v^* Omega v) (real for
   antisymmetric real Omega; the two members of a conjugate pair have opposite sign).  The gauge zero mode and the
   rotation zero mode (the defective 1 of M, +-1 of N) are left out: the four remaining eigenvalues are the two
   ghost-sector pairs.  Each pair is described by psi+ in (-pi, pi]: the angle of its POSITIVE-sign member.
   Two pairs coincide as sets when psi1+ = +-psi2+ (mod 2 pi): SAFE if psi1+ = psi2+ (same-sign members meet,
   they pass), DANGEROUS if psi1+ = -psi2+ (opposite-sign members meet; generically they leave the circle).
4. Tracking.  The four ghost eigenvalues of N are followed continuously in E (slope-aware matching in the complex
   plane; inside a band where they are off the circle the sign is undefined and the tracking is by eigenvalue only).
   For one member of each pair, the sign along the track is printed: a change of sign of a tracked member between
   two on-circle intervals is the "reassignment" that a tongue can produce.
5. CHECK lines over the 400 points: N^2 = M; M^T Omega M = Omega; N^T Omega N = Omega; the trace formula
   A = tr M - 4, B = (A^2 - tr M^2 + 4)/2 against the drop-4 invariants; a + 2 >= 0 for every real a of M.

Predictions tested (WORKORDER v3): (i) every -1 collision of M is a meeting of N at +-i with psi1+ = -psi2+;
(ii) the +1 tongues change the sign relation (the pass-through at E ~ 0.76 is psi1+ = psi2+ while the -1-fold
meetings before 0.49 were psi1+ = -psi2+); (iii) tongues only at +1 and collisions only at -1 at every E.
Found (results/krein_b005.txt, 21 Oct; islands/README.md, "Scope and checks", B005): of eleven M-coincidences of
the two pairs, six are true meetings in N and five are fake degeneracies (N has lambda and -lambda: E = 0.054,
0.056, 0.094, 0.199, 0.283, 0.798); the six are the edges of the three Krein bands, at each the pairs carry opposite
sign (psi1+ = -psi2+ to 0.03-0.09 pi) and meet at N ~ +-i; every parametric event is one pair's own fold at
N = +-1.  So (i) holds; (ii)'s "pass-through at E ~ 0.76" of the 19 Oct brief is the fake degeneracy at 0.798,
where nothing meets.
Not done here: the segment-transfer-matrix localisation (claim 4 of the brief) -- it needs a positive norm adapted
to the block and is left open.

Usage (from the repo root): python islands/code/krein_b005.py [--npts 400] [--block 12] > islands/results/krein_b005.txt
Float, DOP853 rtol 1e-12 (dr_offdiag_floquet.RT).  Added 2026-10-21.
"""
import argparse
import os
import pickle
import sys
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dr_offdiag_floquet as base                                   # noqa: E402

UL = ['u12', 'u13', 'u21', 'u23', 'u31', 'u32']                     # rung1 ulist order
STILDE = {'12': [1, -1], '13': [1, 1], '23': [1, 1]}                # sign of (u_a, u_b) under y -> -y, each derivative order


# ---------------------------------------------------------------- the symplectic form from L2
def build_omega():
    d = pickle.load(open(base.PKL, 'rb'))
    L2, X, P, cj = d['L2'], d['X'], d['P'], d['cj']
    L2 = L2.subs(cj, 1)
    U = [[sp.Symbol(f'u{k}_{j}') for j in range(6)] for k in range(4)]
    state = [U[k][j] for k in range(4) for j in range(6)]              # 24, k-major: matches rung1 'state'
    Vs = (X[0] ** 2 * X[1] ** 2 + X[1] ** 2 * X[2] ** 2 + X[2] ** 2 * X[0] ** 2) / 2
    bg_rate = {X[i]: P[i] for i in range(3)}
    bg_rate.update({P[i]: -Vs.diff(X[i]) for i in range(3)})

    def ddt(expr):
        out = 0
        for k in range(3):
            for j in range(6):
                out += expr.diff(U[k][j]) * U[k + 1][j]
        for s, r in bg_rate.items():
            out += expr.diff(s) * r
        return sp.expand(out)

    Q1 = [U[0][j] for j in range(6)]
    Q2 = [U[1][j] for j in range(6)]
    P2 = [sp.expand(L2.diff(U[2][j])) for j in range(6)]
    P1 = [sp.expand(L2.diff(U[1][j]) - ddt(P2[j])) for j in range(6)]
    canon = Q1 + Q2 + P1 + P2                                           # 24 canonical coordinates
    Jm = sp.Matrix(24, 24, lambda r, c: canon[r].coeff(state[c]))
    # every canonical coordinate must be linear in the state
    lin_ok = all(sp.expand(canon[r] - sum(Jm[r, c] * state[c] for c in range(24))) == 0 for r in range(24))
    Jcan = sp.zeros(24, 24)
    for i in range(12):
        Jcan[i, 12 + i] = 1
        Jcan[12 + i, i] = -1
    Omega = (Jm.T * Jcan * Jm).applyfunc(sp.expand)
    f = sp.lambdify(list(X) + list(P), Omega, 'numpy')
    return f, lin_ok


# ---------------------------------------------------------------- maps
def half_and_full(Kf, E, b):
    """Phi_{T/2}, Phi_T of the block, the background at 0 and T/2."""
    q0, p0, T = base.orbit(E)
    idx = base.block_idx(b)

    def rhs(t, s):
        q, p = s[:3], s[3:6]
        K = np.array(Kf(*q, *p, 1.0), dtype=float)[np.ix_(idx, idx)]
        return np.concatenate([p, -base.gradV(q), (K @ s[6:].reshape(8, 8)).ravel()])
    sol = solve_ivp(rhs, (0.0, T), np.concatenate([q0, p0, np.eye(8).ravel()]), t_eval=[T / 2, T], **base.RT)
    sh, sT = sol.y[:, 0], sol.y[:, 1]
    return sh[6:].reshape(8, 8), sT[6:].reshape(8, 8), (q0, p0), (sh[:3], sh[3:6]), T


def krein_sign(Om, v):
    val = 1j * (np.conj(v) @ Om @ v)
    return float(val.real), float(abs(val.imag))


def ghost_four(lam):
    """the four eigenvalues of N farthest from {+1, -1} (the gauge and rotation sectors sit at +-1)"""
    dist = np.minimum(np.abs(lam - 1.0), np.abs(lam + 1.0))
    return np.argsort(dist)[4:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--npts', type=int, default=400)
    ap.add_argument('--emin', type=float, default=0.05)
    ap.add_argument('--emax', type=float, default=1.0)
    ap.add_argument('--block', default='12')
    a = ap.parse_args()
    t0 = time.time()
    b = a.block
    print(f"krein_b005.py  python {sys.version.split()[0]}  numpy {np.__version__}  sympy {sp.__version__}  ({sys.executable})")
    print(f"run parameters: block {b}, E in [{a.emin}, {a.emax}] x {a.npts} log-spaced, RT={base.RT}, c_j = 1, S~ = {STILDE[b]} on (u_a, u_b)")
    Kf = base.load_K()
    Of, lin_ok = build_omega()
    print(f"CHECK Ostrogradsky coordinates (Q1, Q2, P1, P2) are linear in the state (u, u', u'', u'''): {'PASS' if lin_ok else 'FAIL'}")
    idx = base.block_idx(b)
    St = np.diag(np.array(STILDE[b] * 4, dtype=float))

    grid = np.logspace(np.log10(a.emin), np.log10(a.emax), a.npts)
    rows = []
    worst = dict(N2=0.0, MOM=0.0, NON=0.0, trace=0.0, amin=np.inf, herm=0.0)
    # first E: check Omega on all three blocks with the full period map
    for bb in base.BLOCKS:
        Ph, M, (q0, p0), _, T = half_and_full(Kf, 0.3, bb)
        Om = np.array(Of(*q0, *p0), dtype=float)[np.ix_(base.block_idx(bb), base.block_idx(bb))]
        r = np.abs(M.T @ Om @ M - Om).max() / np.abs(Om).max()
        print(f"CHECK block {bb}, E=0.3: M^T Omega M = Omega to {r:.1e} (tol 1e-7) {'PASS' if r < 1e-7 else 'FAIL'}; "
              f"Omega antisymmetric to {np.abs(Om + Om.T).max():.1e}; rank {np.linalg.matrix_rank(Om)}")

    prev = None
    print("\n== per energy: the four ghost-sector eigenvalues of N = S~ Phi_T/2 (tracked), with Krein sign on the circle ==")
    print("columns: E | for each tracked member k=0..3: |lambda| arg/pi sign | pair structure: psi1+ psi2+ (angle/pi of the positive member of each pair) relation")
    for i, E in enumerate(grid):
        Ph, M, (q0, p0), (qh, ph), T = half_and_full(Kf, E, b)
        N = St @ Ph
        Om = np.array(Of(*q0, *p0), dtype=float)[np.ix_(idx, idx)]
        worst['N2'] = max(worst['N2'], np.abs(N @ N - M).max() / np.abs(M).max())
        worst['MOM'] = max(worst['MOM'], np.abs(M.T @ Om @ M - Om).max() / np.abs(Om).max())
        worst['NON'] = max(worst['NON'], np.abs(N.T @ Om @ N - Om).max() / np.abs(Om).max())
        lam, V = np.linalg.eig(N)
        keep = ghost_four(lam)
        lam4, V4 = lam[keep], V[:, keep]
        # trace formula vs drop-4 on M
        mu = np.linalg.eigvals(M)
        g = base.ghost_multipliers(M)
        A_d = np.sum(g).real
        B_d = ((np.sum(g)) ** 2 - np.sum(g ** 2)).real / 2
        A_t = np.trace(M) - 4
        B_t = (A_t ** 2 - np.trace(M @ M) + 4) / 2
        worst['trace'] = max(worst['trace'], abs(A_d - A_t), abs(B_d - B_t))
        # a + 2 for real a = mu + 1/mu of M's ghost pairs
        for m in g:
            am = m + 1 / m
            if abs(am.imag) < 1e-9 and abs(abs(m) - 1) < 1e-6:
                worst['amin'] = min(worst['amin'], am.real + 2)
        # tracking: match the 4 eigenvalues to the previous step (slope-aware)
        if prev is None:
            order = np.argsort(np.angle(lam4))
        else:
            pred = prev['lam'] if prev.get('lam2') is None else 2 * prev['lam'] - prev['lam2']
            order = np.full(4, -1)
            used = set()
            for k in np.argsort([-abs(p) for p in pred]):
                cands = [(abs(lam4[j] - pred[k]), j) for j in range(4) if j not in used]
                j = min(cands)[1]
                order[k] = j
                used.add(j)
        lam4, V4 = lam4[order], V4[:, order]
        signs, herms = [], []
        for k in range(4):
            if abs(abs(lam4[k]) - 1) < 1e-6:
                s_, h_ = krein_sign(Om, V4[:, k])
                signs.append(np.sign(s_))
                herms.append(h_ / max(abs(s_), 1e-300))
            else:
                signs.append(0.0)
        if herms:
            worst['herm'] = max(worst['herm'], max(herms))
        # pairs: members with conjugate eigenvalues
        pairs = []
        taken = set()
        for k in range(4):
            if k in taken:
                continue
            j = min((abs(lam4[m] - np.conj(lam4[k])), m) for m in range(4) if m != k and m not in taken)[1]
            pairs.append((k, j))
            taken |= {k, j}
        psi = []
        for (k, j) in pairs:
            if signs[k] > 0:
                psi.append(np.angle(lam4[k]))
            elif signs[j] > 0:
                psi.append(np.angle(lam4[j]))
            else:
                psi.append(np.nan)
        rel = '-'
        if not any(np.isnan(psi)) and len(psi) == 2:
            dsum = abs((psi[0] + psi[1] + np.pi) % (2 * np.pi) - np.pi)
            ddif = abs((psi[0] - psi[1] + np.pi) % (2 * np.pi) - np.pi)
            rel = f"d(psi1+psi2)={dsum / np.pi:.3f}pi d(psi1-psi2)={ddif / np.pi:.3f}pi"
        off = [k for k in range(4) if abs(abs(lam4[k]) - 1) >= 1e-6]
        print(f"E={E:.5f} | " + " | ".join(f"{abs(lam4[k]):.4f} {np.angle(lam4[k]) / np.pi:+.4f} {int(signs[k]):+d}" for k in range(4))
              + f" | psi+ = ({', '.join('nan' if np.isnan(p) else f'{p / np.pi:+.4f}' for p in psi)})pi {rel}"
              + (f"  OFF-CIRCLE members {off} (|lambda| {', '.join(f'{abs(lam4[k]):.4f}' for k in off)})" if off else ""), flush=True)
        rows.append(dict(E=E, lam=lam4, signs=signs, psi=psi, off=bool(off)))
        prev = dict(lam=lam4, lam2=None if prev is None else prev['lam'])

    print("\n== checks over the scan ==")
    print(f"CHECK N^2 = M, worst relative residual {worst['N2']:.1e} (tol 1e-8) {'PASS' if worst['N2'] < 1e-8 else 'FAIL'}")
    print(f"CHECK M^T Omega M = Omega, worst {worst['MOM']:.1e} (tol 1e-7) {'PASS' if worst['MOM'] < 1e-7 else 'FAIL'}")
    print(f"CHECK N^T Omega N = Omega (S~ is canonical), worst {worst['NON']:.1e} (tol 1e-7) {'PASS' if worst['NON'] < 1e-7 else 'FAIL'}")
    print(f"CHECK trace formula vs drop-4 invariants, worst |diff| {worst['trace']:.1e} (tol 1e-8) {'PASS' if worst['trace'] < 1e-8 else 'FAIL'}")
    print(f"CHECK a + 2 >= 0 for every real a = mu + 1/mu of M's ghost pairs: min a + 2 = {worst['amin']:.2e} {'PASS' if worst['amin'] > -1e-9 else 'FAIL'}")
    print(f"INFO  Krein form i v* Omega v: worst |imag|/|real| = {worst['herm']:.1e} (should be rounding)")

    # transitions: where a pair leaves / re-enters the circle, what were the signs and relation just before
    print("\n== transitions (a member leaves or re-enters the unit circle), with the sign relation just before ==")
    for i in range(1, len(rows)):
        if rows[i]['off'] != rows[i - 1]['off']:
            r0 = rows[i - 1] if rows[i]['off'] else rows[i]
            kind = 'LEAVE' if rows[i]['off'] else 'RE-ENTER'
            lam = r0['lam']
            on = [k for k in range(4) if abs(abs(lam[k]) - 1) < 1e-6]
            angs = {k: np.angle(lam[k]) / np.pi for k in on}
            psi = r0['psi']
            relation = '-'
            if len(psi) == 2 and not any(np.isnan(psi)):
                dsum = abs((psi[0] + psi[1] + np.pi) % (2 * np.pi) - np.pi)
                ddif = abs((psi[0] - psi[1] + np.pi) % (2 * np.pi) - np.pi)
                relation = 'DANGEROUS (psi1+ = -psi2+)' if dsum < ddif else 'SAFE (psi1+ = psi2+)'
                relation += f"  [d(sum) {dsum / np.pi:.3f}pi, d(diff) {ddif / np.pi:.3f}pi]"
            # which members are off the circle inside the band, and do they form one conjugate pair (a fold of that
            # pair: parametric, M's tongue at +1 = N's eigenvalue at +-1) or come from two different pairs (a collision)
            r1 = rows[i] if rows[i]['off'] else rows[i - 1]
            offk = [k for k in range(4) if abs(abs(r1['lam'][k]) - 1) >= 1e-6]
            lam_on = r0['lam']
            same_pair = (len(offk) == 2 and abs(lam_on[offk[0]] - np.conj(lam_on[offk[1]])) < 1e-3 * 2)
            if len(offk) == 4:
                what = 'ALL FOUR off: two-pair COLLISION (Krein quadruplet)'
            elif same_pair:
                what = f"one pair's own members meet at arg {angs.get(offk[0], float('nan')):+.3f}pi: FOLD (parametric)"
            else:
                what = f"members {offk} from different pairs: COLLISION"
            print(f"  E={rows[i]['E']:.5f} {kind}: {what}; on-circle just outside the band: "
                  + ", ".join(f"k{k} arg {angs[k]:+.3f}pi sign {int(r0['signs'][k]):+d}" for k in on)
                  + f" ; two-pair relation {relation}")
    # meetings of the two pairs: a coincidence of M's multipliers (N-angles equal mod pi) is a TRUE meeting only if
    # N's eigenvalues coincide (equal mod 2 pi); if N has lambda and -lambda the degeneracy of M is fake and nothing
    # meets.  True meetings are DANGEROUS (psi1+ = -psi2+) or SAFE (psi1+ = psi2+).
    print("\n== meetings of the two pairs (local minima of the angular distance, < 0.12 pi) ==")

    def dist(x, per):
        return abs((x + per / 2) % per - per / 2)
    cand = []
    for i, r in enumerate(rows):
        if r['off'] or any(np.isnan(r['psi'])) or len(r['psi']) != 2:
            cand.append((np.inf, np.inf, np.inf, np.inf))
            continue
        p1, p2 = r['psi']
        cand.append((dist(p1 + p2, 2 * np.pi), dist(p1 - p2, 2 * np.pi), dist(p1 + p2, np.pi), dist(p1 - p2, np.pi)))
    cand = np.array(cand)
    dM = np.minimum(cand[:, 2], cand[:, 3])
    for i in range(1, len(rows) - 1):
        if dM[i] <= dM[i - 1] and dM[i] <= dM[i + 1] and dM[i] < 0.12 * np.pi:
            dN_sum, dN_dif = cand[i, 0], cand[i, 1]
            if min(dN_sum, dN_dif) < 0.12 * np.pi:
                kind = 'TRUE meeting in N: ' + ('DANGEROUS (psi1+ = -psi2+)' if dN_sum < dN_dif else 'SAFE (psi1+ = psi2+)')
            else:
                kind = 'FAKE degeneracy of M (N has lambda and -lambda): no meeting'
            print(f"  E={rows[i]['E']:.4f}: M-coincidence to {dM[i] / np.pi:.3f}pi; N distances sum {dN_sum / np.pi:.3f}pi diff {dN_dif / np.pi:.3f}pi -> {kind}")
    # band edges are meetings too (the pairs leave the circle at once); list the dangerous ones with their angles
    print("  (the collision band edges listed under 'transitions' are the dangerous meetings that did collide)")

    # sign reassignment along each tracked member
    print("\n== sign of each tracked member over the on-circle intervals ==")
    for k in range(4):
        segs = []
        cur = None
        for r in rows:
            s_ = int(r['signs'][k])
            if s_ == 0:
                if cur:
                    segs.append(cur)
                    cur = None
                continue
            if cur is None or cur[2] != s_:
                if cur:
                    segs.append(cur)
                cur = [r['E'], r['E'], s_]
            else:
                cur[1] = r['E']
        if cur:
            segs.append(cur)
        print(f"  member k{k}: " + "  ".join(f"[{s[0]:.4f},{s[1]:.4f}]:{s[2]:+d}" for s in segs))
    print(f"\nDONE {time.time() - t0:.0f}")


if __name__ == "__main__":
    main()
