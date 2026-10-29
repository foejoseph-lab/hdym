# Rabbit hole 8 — the loophole at the paper's seed size, measured (22 Oct 2026)

`islands/results/esc_E_delta1e-2_b12.{txt,npz}`: section E, 48 energies in [0.045, 0.2] × 64 seeds, **off-diagonal
seed δ = 10⁻² in block 12 only, diagonal ghost not seeded** (so that any escape at a Floquet-stable energy is the
off-diagonal transient's doing, not the diagonal ghost's own tongue), tmax 3000, 7963 s on 24 threads.

**Two CHECK lines fail and are reported as such**: the kernel-vs-scipy cross-checks at δ = 10⁻² miss their
tolerances (state at 2T: 8 × 10⁻³ relative against tol 10⁻⁶; escape time at E = 0.549: 1.8 × 10⁻² against 10⁻³).
The tolerances were set for 10⁻⁴ seeds; at 10⁻² the orbit's neighbourhood is sensitive enough that integrators at
rtol 10⁻¹¹ and 10⁻¹² part at the percent level within two periods. The escape *fractions* below are robust to that
(a 2 % shift of individual times does not move a 64-cell fraction); the individual times are at 1 % accuracy, not
the paper's 10⁻³. No integrator failures; energy drift ≤ 2 × 10⁻⁴; Gauss residual of the seed 7 × 10⁻⁵ = O(δ²), conserved.

## Result: at δ = 10⁻² the escape does not follow the linear bands at all

| E range | block 12 linear | diagonal ghost | escaped by 3000 at δ = 10⁻² (blk-12 seed) | at δ = 10⁻⁴ (both seeds, `esc_E`) |
|---|---|---|---|---|
| 0.045–0.075 | stable (two narrow bands inside) | stable | 0–1.6 % | 0–3 % |
| 0.080–0.100 | **stable** | stable (tongue at 0.089–0.092) | **6 → 22 %, rising smoothly** | 0 % off the tongue |
| 0.102–0.113 | Krein/param. band | stable | 23–39 % | 19–34 % |
| **0.117–0.150** | **stable** | **stable** (0.124 tongue, rate 0.0035) | **44 → 83 %, rising smoothly** | **0 %** (except the tongue) |
| 0.152–0.161 | Krein band | stable | 88–92 % | 86 % |
| 0.165–0.188 | stable | unstable tongue | 89–97 % | 100 % |
| **0.194–0.200** | **stable** | **stable** | **89–98 %** | **0 %** |

- The escape fraction at δ = 10⁻² is a smooth, monotone function of E — ≈ 0 below 0.08, 50 % at E ≈ 0.12, > 90 %
  above 0.15 — with no feature at any band edge. The Floquet-stable gaps 0.117–0.150 and 0.194–0.200, where at
  δ = 10⁻⁴ *nothing* escapes by 3000, escape at 44–98 %.
- The mechanism is the one the hole named: the diagonal ghost reaches |w| = 1 at t ≈ 60–150 at **every** energy,
  stable or not (median of t_cross(|w| = 1) over the 48 rows: 56–154) — the band-edge beat amplifies the seed by
  ~10³ (docstring: 10⁻⁴ → 0.23 at E = 0.05) and a 10⁻² seed is carried to order one within two periods. After that
  the trajectory is in the nonlinear regime and whether it escapes by 3000 is a smooth function of E, unrelated to
  the linear classification; below E ≈ 0.08 it sits at order one for the whole run (the τ = 100-orbit Poisson regime
  of the 12000 run, or longer).
- At δ = 10⁻⁴ the linear bands predict the escapes essentially perfectly (`RH4_loiter.md`); at δ = 10⁻²
  they predict nothing. The crossover seed is between, where the transient gain × δ reaches order one: ~10⁻³.

## Reading for the paper

"Floquet-stable means safe" holds on this orbit only for seeds below ~10⁻³ — two decades under the paper's
ε = 0.01. At the paper's own regime the escape from the orbit is governed by the nonlinear dynamics after a linear
transient has already delivered the ghost to order one, and the linear typing (Krein vs parametric, the bands, the
ladder) describes the δ → 0 limit, not the physical seed. This is consistent with the paper's main statement (no
threshold, no exponential suppression) and sharpens the scope of its off-diagonal section: the band structure is a
property of the linearised orbit, and the nonlinear escape at ε = 0.01 is smooth in E. The one new number: on the
DR orbit at δ = 10⁻², the escape-by-3000 fraction crosses 50 % at E ≈ 0.12 — a finite-time statement (the 12000 run
says the low-E side is a slow Poisson tail, not a stable regime).

Hole closed. (The δ-crossover scan, 10⁻⁴ → 10⁻² in half-decades at two energies, would pin the ~10⁻³; one hour. Done
28–29 Oct 2026: next section.)

## The δ-crossover scan (28–29 Oct 2026)

`islands/results/esc_dscan_<δ>_b12.{txt,npz}`, δ = 10⁻⁴, 3.1623 × 10⁻⁴, 10⁻³, 3.1623 × 10⁻³: section E at two
energies, the Floquet-stable rows 34 and 46 of `esc_E_delta1e-2_b12` (E = 0.13239 in the gap 0.117–0.150 and
E = 0.19375 in 0.194–0.200), 64 cells each, block-12 seed only, diagonal ghost not seeded, tmax 3000, 6 threads. The
same `--seed 1` draws the same 64 directions at every δ, scaled to δ, so the rows differ only in δ. The δ = 10⁻² row is
the 22 Oct run's.

| δ | E = 0.13239: escaped | \|w\| ≥ 1 | E = 0.19375: escaped | \|w\| ≥ 1 |
|---|---|---|---|---|
| 10⁻⁴ | 0/64 | 0/64 | 0/64 | 0/64 |
| 3.2 × 10⁻⁴ | 0/64 | 0/64 | 0/64 | 0/64 |
| 10⁻³ | 2/64 (3 %) | 21/64 | 14/64 (22 %) | 21/64 |
| 3.2 × 10⁻³ | 24/64 (38 %) | 57/64 | 41/64 (64 %) | 52/64 |
| 10⁻² | 46/64 (72 %) | 64/64 | 57/64 (89 %) | 64/64 |

- **The crossover is pinned: escape by t = 3000 starts between δ = 3.2 × 10⁻⁴ and 10⁻³ at both energies** [grid]. Half
  the cells escape at δ ≈ 5 × 10⁻³ (E = 0.132) and ≈ 2 × 10⁻³ (E = 0.194) [est, log-linear between grid points]. The
  "~10⁻³" above is the onset, and "Floquet-stable means safe" holds on these 128 cells for δ ≤ 3.2 × 10⁻⁴.
- \|w\| = 1 is necessary [grid]: no cell escapes without the diagonal ghost reaching 1 first (0 of 128 at each δ), and
  26 (δ = 10⁻³) and 44 (3.2 × 10⁻³) cells reach it and do not escape by 3000. Below the onset the transient still
  pumps: \|w\| reaches 0.1 in 56–66 % of cells at δ = 10⁻⁴ and 84–86 % at 3.2 × 10⁻⁴ (median t 1495–1753 and 531–713),
  a gain of 10³ or more without instability.
- The median time to \|w\| = 1 falls from 340–355 (δ = 10⁻³) to 118–169 (3.2 × 10⁻³) and 73–82 (10⁻²).
- Checks: the integrator, energy, Gauss, determinism and no-failure CHECKs pass at every δ. The single-trajectory
  escape-time CHECK at E = 0.549 passes at 10⁻⁴ (8.1 × 10⁻⁴) and 10⁻³ (8.2 × 10⁻⁵) and fails at 3.2 × 10⁻⁴
  (1.2 × 10⁻¹) and 3.2 × 10⁻³ (1.7 × 10⁻¹); the logs at those two δ end SOME CHECK FAILED. Not monotone in δ, so it is
  the sensitivity of one chaotic escape time, as at 10⁻² above; the fractions do not depend on it [est].

The physics reading of the onset is the author's.
