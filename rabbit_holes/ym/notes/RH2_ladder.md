# Rabbit hole 2 — the resonance ladder of block 12's tongue edges, as a formula check

Session of 21 Oct 2026 (Claude, cloud). Scripts `code/rh2_ladder.py` (parses B005 and B001, no integration),
`code/rh2_predict.py` (72 half-period maps, 41 s), `code/rh2_rotation.py` (tried, see below). Logs
`results/rh2_ladder.txt`, `results/rh2_predict.txt`; figures `results/rh2_ladder.png`. Scope: block 12 on the
planar DR orbit, c_j = 1. Float.

## The ladder, measured

B001's seven block-12 bands alternate parametric / Krein going down in E. Their centres (geometric mean of the
edges), written as the orbit period T(E) = T_half (2E)^{−1/4}, are equally spaced:

| rung n | band | type | E_c | T_c | T_c − 4.057 (n + 0.210) |
|---|---|---|---|---|---|
| 13 | [0.0547, 0.0551] | P (pair k0/k3 at N = −1) | 0.0549 | 53.95 | +0.35 (+0.09 rung) |
| 13 | [0.0577, 0.0594] | P (pair k1/k2 at N = +1) | 0.0585 | 53.09 | −0.51 (−0.13 rung) |
| 12 | [0.0761, 0.0773] | K | 0.0767 | 49.62 | +0.08 |
| 11 | [0.1020, 0.1133] | P | 0.1075 | 45.60 | +0.12 |
| 10 | [0.1519, 0.1613] | K | 0.1565 | 41.52 | +0.09 |
| 9 | [0.2228, 0.2648] | P | 0.2429 | 37.20 | −0.17 |
| 8 | [0.2811, 0.4937] | K | 0.3725 | 33.42 | +0.12 |
| 7 | [0.5694, 0.7241] | P | 0.6421 | 29.17 | −0.08 |

Least squares: **T_c(n) = 4.057 (n + 0.210)**, rms residual 0.11 time units (≤ 0.04 rung for n = 7–12). Odd n
parametric, even n Krein. Rung 13 appears twice — once per pair, the two copies 0.21 rung apart (the pairs'
phase splitting, hole 7).

Reading. Consecutive bands alternate M = +1 and M = −1, so the ghost pair's phase per period advances by π
from one rung to the next: ω_eff ΔT = π with ΔT = 4.06 gives **ω_eff = 0.774**, constant across E = 0.055–0.64
(T from 54 to 29): over this range the ghost on the DR orbit accumulates phase like a free mode of frequency
0.77, i.e. 23 % below the Lee–Wick mass, independent of E. This is the formula the hole asked for:
"ω_ghost(E) T(E)/2π = n/2" holds with ω_ghost = 0.774 and a phase offset of 0.21 rung. The naive version with
ω = 1 (the first table of `rh2_ladder.txt`) lands within a rung but with the parity wrong below E ≈ 0.1.

## Out-of-sample test (`rh2_predict.py`, made before looking)

Predictions: rung 14 (Krein) at E = 0.0421, rung 15 (parametric) at E = 0.0321, rung 6 (Krein) at E = 1.154.

- Low E, scanned [0.025, 0.05]: **a Krein quadruplet at E = 0.0427–0.0435, rung 13.9** — predicted 0.0421 at
  rung 14.0: confirmed to 0.1 rung. A parametric fold at E = 0.0350, rung 14.67 (predicted 15.0; one copy of
  the split pair, cf. rung 13's two copies at 13.2 and 12.9). Nothing else down to rung 16.
- High E, scanned [1.0, 1.5]: a **parametric** band from E = 1.014 to 1.271 (rungs 6.20 → 5.85) and a Krein band
  from E = 1.286 on (rung 5.83 → 5.6, still open at 1.5). The ladder's rung 6 is a parametric band, not the
  Krein one the alternation predicts; the bands are now wider than a quarter rung and adjacent.

So the ladder is a clean formula from rung 14 to rung 7 (E ≈ 0.04–0.72) and **breaks above E ≈ 1**, where the
two pairs — 0.2 rung apart at the bottom of the ladder — have separated by a full rung (at E = 1 one pair sits
at M = +1 and the other at M = −1, B005) and the alternation of types with n no longer holds. That is the
quantitative version of the list's "tongues overlap above E ≈ 0.31": the bands widen (Krein band widths 0.05,
0.15, 1.16 rung for n = 12, 10, 8) while their centres stay on the ladder until the splitting reaches a rung.

## What did not work: a canonical rotation number

The ladder uses band centres only. Between bands B005 gives each pair's N-angle mod 2π, and unwrapping it
in E across a band is ambiguous (inside the wide Krein band the pairs could wind by up to 2π unseen;
`rh2_ladder.py`'s unwrapped totals, −1.7 π and −3.5 π over the scan against the 9.3 π of the ω = 1 line, are not
trustworthy for that reason). `rh2_rotation.py` tried the phase along the orbit instead — the unwrapped
argument of the complex u12 component of the Floquet solution over a half period, whose exponential equals
the eigenvalue of N exactly (residual 2 × 10⁻¹¹) — but u12 passes within 10⁻³ of zero along the orbit (the
ghost mode lives in u21 for half the period), so its winding count is not canonical; at E = 0.05 it gives
±13.2 π and ±1.7 π for the two pairs, at E = 1 ±0.52 π and ±0.006 π, which is the mod-2π angle with an
arbitrary integer. A Conley–Zehnder / Maslov index of the symplectic path restricted to the pair's invariant
2-plane would be the canonical object; not attempted (it is a day, and the band-centre ladder already gives
the frequency).

## For the paper (if wanted)

One figure and two sentences: the band centres of block 12 lie on T_c = 4.06 (n + 0.21), n = 7–14, alternating
parametric (odd n) and Krein (even n) — the ghost accumulates phase at the constant effective frequency 0.77
over a factor 16 in E — with the out-of-sample Krein band at E = 0.043 as the check; above E ≈ 1 the two pairs
are a rung apart and the alternation ends. Numbers from `results/rh2_predict.txt` and the table above.

## Open

- Why 0.774? A WKB phase ∫ω(t) dt of the frozen ghost frequency over the orbit would test whether the 23 %
  redshift is the time average of the background's effect; the frozen block has frequencies 1.85 and 2.83
  and a real pair ±1.32 at the turning point (E = 0.3), so "the" ghost frequency along the orbit needs the
  dressed-sector picture of hole 7 first.
- The pair splitting as a function of E (0.2 rung at the bottom, 1 rung at E = 1): a curve, readable from the
  two copies of each parametric rung where both are visible (rungs 13 and 15 so far).
