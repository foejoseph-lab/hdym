# Rabbit hole 2 — the resonance ladder of block 12's bands (a measurement; not a frequency)

21 Oct 2026. Scripts `code/rh2_ladder.py` (parses B005 and B001, no integration), `code/rh2_predict.py` (72
half-period maps, 41 s), `code/rh2_rotation.py` (a rotation-number attempt, see below). Logs `results/rh2_ladder.txt`,
`results/rh2_predict.txt`; figure `results/rh2_ladder.png`; data `results/rh2_ladder.npz`. Scope: block 12 on the
planar DR orbit, c_j = 1. Float.

## The ladder

B001's block-12 bands alternate parametric / Krein going down in E. Their centres (geometric mean of the edges),
written as the orbit period T(E) = T_half (2E)^{−1/4}, are equally spaced:

| rung n | band | type | E_c | T_c | T_c − 4.057 (n + 0.210) |
|---|---|---|---|---|---|
| 14 | [0.0427, 0.0435] (out of sample, below) | K | — | — | predicted at E = 0.0421 |
| 13 | [0.0547, 0.0551] | P (pair k0/k3 at N = −1) | 0.0549 | 53.95 | +0.35 (+0.09 rung) |
| 13 | [0.0577, 0.0594] | P (pair k1/k2 at N = +1) | 0.0585 | 53.09 | −0.51 (−0.13 rung) |
| 12 | [0.0761, 0.0773] | K | 0.0767 | 49.62 | +0.08 |
| 11 | [0.1020, 0.1133] | P | 0.1075 | 45.60 | +0.12 |
| 10 | [0.1519, 0.1613] | K | 0.1565 | 41.52 | +0.09 |
| 9 | [0.2228, 0.2648] | P | 0.2429 | 37.20 | −0.17 |
| 8 | [0.2811, 0.4937] | K | 0.3725 | 33.42 | +0.12 |
| 7 | [0.5694, 0.7241] | P | 0.6421 | 29.17 | −0.08 |

Least squares on n = 7–13: **T_c(n) = 4.057 (n + 0.210)**, rms residual 0.11 time units (≤ 0.04 rung for n = 7–12). (The fit is in
no committed script; refitting the eight rows of this table gives 4.057 / 0.210 with rms 0.24 over all eight and 0.115 over n = 7–12.)
Odd n parametric, even n Krein (the naive ω = 1 version, the first table of `rh2_ladder.txt`, lands within a rung but with
the parity wrong below E ≈ 0.1). Rung 13 appears twice, once per pair, the two copies 0.21 rung apart — the pairs'
phase splitting Δ(E) at that energy (`RH27_averaging.md`, item 5).

**Out of sample** (`rh2_predict.py`, predictions written before scanning): rung 14 (Krein) predicted at
E = 0.0421 — a Krein quadruplet found at E = 0.0427–0.0435, rung 13.9; rung 15 (parametric) predicted at 0.0321 — a
parametric fold at E = 0.0350, rung 14.67, one copy of the split pair; nothing else down to rung 16. Rung 6 (Krein
by alternation) predicted at E = 1.154 — found instead a **parametric** band from E = 1.014 to 1.271 (rungs 6.20 →
5.85) and a Krein band from E = 1.286 on (rung 5.83 → 5.6, still open at 1.5). So the ladder is a clean formula
from rung 14 to rung 7 (E ≈ 0.04–0.72) and the alternation of types ends above E ≈ 1, where the two pairs — 0.2
rung apart at the bottom — have separated by a full rung (Δ(E) = 0.95 π at E = 1: one pair at M = +1, the other at
M = −1, B005) and the "wrong" pair folds at the odd rung. Band widths grow meanwhile (Krein widths 0.05, 0.15,
1.16 rung for n = 12, 10, 8) while the centres stay on the ladder until the splitting reaches a rung.

## What the ladder is and is not

The content is the slope and the count: consecutive bands alternate M = +1 and M = −1, so the per-period phase of
the ghost pair advances by π per rung, and the rungs n = 7–14 are the eight bands of E ∈ [0.04, 0.72], one each.
The number π/4.057 = 0.774 that the alternation gives as "ω_eff ΔT = π" is **not a frequency**:
the ladder uses ωT mod 2π only, and ω → ω + kΩ shifts ωT by 2πk, leaving every rung position, the count and the
parity unchanged and moving only the label n by 2k. The physical frequency of the pairs was measured otherwise
— spectra of the Floquet solutions' periodic parts (`code/rh5_branch.py`, `results/rh5_branch.txt`): a line at
|f| ≈ 0.85–1.2 plus a quasi-static line, the same two lines in both pairs; and the instantaneous frequency along
the orbit (`code/rh5_instfreq.py`, `results/rh5_instfreq.png`), a swept oscillator. The "constant effective
frequency 0.77, 23 % below the Lee–Wick mass" of the first reading (21 Oct) is withdrawn; do not quote it
(`FOR_THE_REWRITE_2026-10-25.md`, item 6 and section D).

## What did not work: a canonical rotation number

The ladder uses band centres only. Between bands B005 gives each pair's N-angle mod 2π, and unwrapping it in E
across a band is ambiguous (inside the wide Krein band the pairs could wind by up to 2π unseen; `rh2_ladder.py`'s
unwrapped totals, −1.7 π and −3.5 π over the scan against the 9.3 π of the ω = 1 line, are not trustworthy for that
reason). `rh2_rotation.py` tried the phase along the orbit instead — the unwrapped argument of the complex u12
component of the Floquet solution over a half period, whose exponential equals the eigenvalue of N exactly
(residual 2 × 10⁻¹¹, no committed log) — but u12 passes within 10⁻³ of zero along the orbit (the ghost mode lives in u21 for half the
period), so its winding count is not canonical: ±13.2 π and ±1.7 π for the two pairs at E = 0.05, ±0.52 π and
±0.006 π at E = 1 (no committed log). Signature-oriented unwrapping (`code/rh5_unwrap.py`, `results/rh5_unwrap.txt`) gives a slope
that is one branch of the line near 1 (0.76–0.83 in three windows below E ≈ 0.22; the branch is fixed by the
spectrum, not by the unwrapping). A Conley–Zehnder / Maslov index of the symplectic path restricted to the pair's
invariant 2-plane would be the canonical object; not attempted.

## For the paper

One paragraph, one figure if wanted (`FOR_THE_REWRITE_2026-10-25.md`, item 6): the band centres of block 12 lie
on T_c = 4.057 (n + 0.210), n = 7–14, parametric at odd n and Krein at even n, rms 0.11 time units; the
out-of-sample Krein band at E = 0.043; the alternation ends above E ≈ 1 (the band at E = 1.014–1.271, rung 6, is parametric where alternation would put a Krein band). Numbers from
`results/rh2_predict.txt` and the table above.
