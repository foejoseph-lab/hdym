# Rabbit hole 10 — the ε-crossover of the E = 0.3 island (22 Oct 2026; four `hdym_island.py` runs)

Logs `islands/results/island_E0.3_eps1e-{9,8,7,6}.log` (ball sampling around the DR point, 200 points per shell,
8 shells r = 10⁻⁶ … 10⁻², tmax 5000, Yoshida-4 pair kernel; ~15 s each). The diagonal ghost's amplitude is
|w| = √(2εE); the island half-width at E = 0.3 is 1.8 × 10⁻⁴ (scaled from the certified E = ½ value).

| ε | |w| | |w| / half-width | core shells (r ≤ 10⁻⁴, 800 points): escaped | core FTLE | r = 2 × 10⁻⁴ shell |
|---|---|---|---|---|---|
| 10⁻⁹ | 2.4 × 10⁻⁵ | 0.14 | 0 | −7 × 10⁻⁴ … 2 × 10⁻⁴ | 0.5 % |
| 10⁻⁸ | 7.7 × 10⁻⁵ | 0.43 | 0 | ≈ 0 | 5 % |
| 10⁻⁷ | 2.4 × 10⁻⁴ | 1.4 | 21–25 % | 4–10 × 10⁻⁴ | 27 % |
| 10⁻⁶ | 7.7 × 10⁻⁴ | 4.3 | 90–92 % | 0.01–0.03 | 90 % |
| (sea control) | — | — | 100 % | 0.29–0.31 | — |

- Intact hole up to ε = 10⁻⁸; begins to leak at ε = 10⁻⁷, exactly where |w| crosses the island half-width (ratio
  1.4) — with the core FTLE still at zero: the trajectories are still on tori and a quarter of them are carried out
  by the ghost's back-reaction, not by chaos; erased at ε = 10⁻⁶ (ratio 4.3), where the core escapes like the halo
  and the FTLE is halo-like (0.01–0.03).
- The halo shells r = 5 × 10⁻⁴, 10⁻³ (47–60 % and 78–81 % escape, FTLE 0.002–0.01) are ε-independent below 10⁻⁷:
  they are the base flow's, as they should be.
- Crossover criterion: |w| ≈ island half-width, i.e. ε_c ≈ (half-width)²/(2E) ≈ 5 × 10⁻⁸ at E = 0.3; erasure one
  decade above. This scopes the paper's "no interior at ε = 0.01" sentence: the island survives the diagonal ghost
  only for ε ≲ 10⁻⁷, six decades below the paper's regime (and it is not a refuge for the full theory at any ε,
  islands/README off-diagonal item).

Hole closed. (All diagonal-ghost. The nine-component version is hole 9, run 23 Oct, `islands/results/island_full_E0.3_eps1e-8.*`: with the off-diagonal sector seeded at 10⁻⁸ every core shell escapes, with none seeded none does; reading in `../README.md`, row 9.)
