"""rh2_ladder.py -- rabbit hole 2: the resonance ladder of the block-12 tongue edges as a formula check.

Input: islands/results/krein_b005.txt (the four tracked ghost-sector eigenvalues of N = S~ Phi_T/2 at 400 energies)
and islands/results/dr_offdiag_floquet_typed.txt (the measured unstable intervals).  No integration.

The check.  Unwrap the tracked N-angles psi_k(E) continuously in E (inside a band the off-circle members are
interpolated across the band from the two edges).  Theta_k(E) = 2 psi_k, unwrapped, is the accumulated Floquet
phase of pair k over one period, an integer multiple of pi at every tongue edge: parametric folds at Theta = 2 pi m
(M = +1, N = +-1), Krein collisions where Theta_1 + Theta_2 = 2 pi m' with both near an odd multiple of pi (M = -1,
N = +-i; rabbit hole 7).  The naive ladder "omega_ghost T / 2 pi = n / 2" with the Lee-Wick mass omega = 1 is the
straight line Theta = T(E) = T_half (2E)^{-1/4}; the measured Theta_k(E)/T(E) is the effective ghost frequency
omega_k(E) of each pair, the quantity the ladder needs.

Output: results/rh2_ladder.txt (this print), results/rh2_ladder.png, results/rh2_ladder.npz.
Usage (from the repo root): python rabbit_holes/ym/code/rh2_ladder.py > rabbit_holes/ym/results/rh2_ladder.txt
Added 2026-10-21 (rabbit-hole session).
"""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
RES = os.path.join(HERE, '..', 'results')
B005 = os.path.join(ROOT, 'islands', 'results', 'krein_b005.txt')
B001 = os.path.join(ROOT, 'islands', 'results', 'dr_offdiag_floquet_typed.txt')
T_HALF = 31.053330268585066


def parse_b005():
    E, mod, ang, sgn = [], [], [], []
    pat = re.compile(r'E=([\d.]+) \| (.*?) \| psi\+')
    for line in open(B005):
        m = pat.match(line.strip())
        if not m:
            continue
        E.append(float(m.group(1)))
        cols = [c.split() for c in m.group(2).split(' | ')]
        mod.append([float(c[0]) for c in cols])
        ang.append([float(c[1]) for c in cols])
        sgn.append([int(c[2]) for c in cols])
    return np.array(E), np.array(mod), np.array(ang), np.array(sgn)


def intervals_b001(block='12'):
    txt = open(B001).read().split('== Unstable intervals')[1]
    out, cur = [], None
    for line in txt.splitlines():
        m = re.match(r'\s*interval \d+: E in \[([\d.]+), ([\d.]+)\]', line)
        if m:
            cur = (float(m.group(1)), float(m.group(2)))
            continue
        m = re.match(r'\s*blk(\d\d): max \|mu\| = ([\d.]+) at E=([\d.]+); type sequence (.*?);', line)
        if m and cur:
            out.append((m.group(1), cur[0], cur[1], 'CPLX' if 'CPLX' in m.group(4) else 'REAL', float(m.group(2))))
    return out


def unwrap_member(E, mod, ang, k):
    """continuous unwrapped angle (units of pi) of tracked member k; off-circle stretches interpolated linearly
    between the last on-circle value before and the first after (the fold/collision point itself is at an
    integer or half-integer, so the interpolation passes through it)."""
    on = np.abs(mod[:, k] - 1) < 1e-4
    a = ang[:, k].copy()
    # unwrap only across on-circle neighbours
    out = np.full(len(E), np.nan)
    last = None
    for i in range(len(E)):
        if not on[i]:
            continue
        if last is None:
            out[i] = a[i]
        else:
            d = a[i] - a[last]
            d = (d + 1) % 2 - 1
            out[i] = out[last] + d
        last = i
    # fill the bands by interpolation in log E
    idx = np.where(~np.isnan(out))[0]
    out = np.interp(np.log(E), np.log(E[idx]), out[idx])
    return out, on


def main():
    E, mod, ang, sgn = parse_b005()
    T = T_HALF * (2 * E) ** -0.25
    print(f"rh2_ladder.py: {len(E)} energies from krein_b005.txt; T(E) = T_half (2E)^(-1/4), T_half = {T_HALF}")
    th = {}
    for k in range(4):
        u, on = unwrap_member(E, mod, ang, k)
        th[k] = u
    # the two pairs: (k0, k3) and (k1, k2) are conjugate members (B005's pairing at E = 0.05: +-0.8426, +-0.2861)
    # unwrapped positive-member angle of each pair, in units of pi, 2 psi = Theta/pi
    psiA = (th[3] - th[0]) / 2        # pair A half-difference = |angle| continued
    psiB = (th[2] - th[1]) / 2
    ThA, ThB = 2 * psiA, 2 * psiB     # accumulated phase per period / pi
    # pin the branch: choose the integer offsets so that Theta/T is positive and of order the Lee-Wick frequency
    wA, wB = np.pi * ThA / T, np.pi * ThB / T
    print(f"raw unwrapped Theta/pi at E = {E[0]:.3f}: A {ThA[0]:+.3f}, B {ThB[0]:+.3f}; at E = 1: A {ThA[-1]:+.3f}, B {ThB[-1]:+.3f}")
    print(f"total winding over the scan: A {ThA[-1] - ThA[0]:+.3f} pi, B {ThB[-1] - ThB[0]:+.3f} pi; T goes {T[0]:.2f} -> {T[-1]:.2f}")
    # the ladder line omega = 1
    n_line = T / np.pi
    print("\n== naive ladder: omega_ghost = 1, tongues at T(E) = n pi ==")
    print("  n   E_n (T = n pi)   nearest measured block-12 band (B001)        type expected (even n: parametric, odd n: Krein)")
    ivs = intervals_b001()
    for n in range(8, 20):
        En = (T_HALF / (n * np.pi)) ** 4 / 2
        best = min(ivs, key=lambda iv: abs(np.log(En / np.sqrt(iv[1] * iv[2]))))
        inside = best[1] <= En <= best[2]
        print(f"  {n:2d}   {En:.4f}          blk{best[0]} [{best[1]:.4f}, {best[2]:.4f}] {best[3]:4s} {'(E_n inside)' if inside else '            '}   {'parametric' if n % 2 == 0 else 'Krein'}")
    np.savez(os.path.join(RES, 'rh2_ladder.npz'), E=E, T=T, ThetaA=ThA, ThetaB=ThB, thetas=np.array([th[k] for k in range(4)]), mod=mod, ang=ang, sgn=sgn)

    # per-band: the unwrapped phase at the band edges
    print("\n== measured block-12 bands against the unwrapped phases (Theta/pi of each pair at the edges) ==")
    for blk, lo, hi, ty, mu in ivs:
        if blk != '12':
            continue
        i0, i1 = np.searchsorted(E, lo), np.searchsorted(E, hi)
        i0 = max(i0 - 1, 0)
        i1 = min(i1 + 1, len(E) - 1)
        print(f"  [{lo:.4f}, {hi:.4f}] {ty}: Theta_A/pi {ThA[i0]:+.3f} -> {ThA[i1]:+.3f}, Theta_B/pi {ThB[i0]:+.3f} -> {ThB[i1]:+.3f}, (A+B)/pi {ThA[i0] + ThB[i0]:+.3f} -> {ThA[i1] + ThB[i1]:+.3f}, T/pi {T[i0]:.2f} -> {T[i1]:.2f}")

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
        for blk, lo, hi, ty, mu in ivs:
            if blk == '12':
                for axx in ax:
                    axx.axvspan(lo, hi, color=('tab:red' if ty == 'CPLX' else 'tab:blue'), alpha=0.15, lw=0)
        for k, c in zip(range(4), ['C0', 'C1', 'C2', 'C3']):
            ax[0].plot(E, th[k], c, lw=1, label=f'member k{k}')
        ax[0].plot(E, T / np.pi / 2, 'k--', lw=1, label='T/(2 pi): omega = 1 line (N-angle)')
        ax[0].plot(E, -T / np.pi / 2, 'k--', lw=1)
        ax[0].set_ylabel('unwrapped N-angle / pi')
        ax[0].legend(fontsize=7, ncol=3)
        ax[1].plot(E, wA, 'C0', lw=1.2, label='omega_A = Theta_A / T')
        ax[1].plot(E, wB, 'C1', lw=1.2, label='omega_B = Theta_B / T')
        ax[1].axhline(1.0, color='k', ls='--', lw=1, label='Lee-Wick mass 1')
        ax[1].set_ylabel('effective frequency per period')
        ax[1].set_xlabel('E')
        ax[1].set_xscale('log')
        ax[1].legend(fontsize=7)
        ax[0].set_title('rabbit hole 2: unwrapped half-period Floquet angles of block 12 (B005) vs its bands (red Krein, blue parametric)', fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(RES, 'rh2_ladder.png'), dpi=130)
        print("figure: results/rh2_ladder.png")
    except Exception as e:                                            # noqa: BLE001
        print(f"(no figure: {e})")


if __name__ == '__main__':
    main()
