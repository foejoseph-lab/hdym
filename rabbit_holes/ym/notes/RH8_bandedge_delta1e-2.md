# Rabbit hole 8 — the loophole at the paper's seed size, measured (22 Oct 2026; Joseph's run)

`islands/results/esc_E_delta1e-2_b12.{txt,npz}`: section E, 48 energies in [0.045, 0.2] × 64 seeds, **off-diagonal
seed δ = 10⁻² in block 12 only, diagonal ghost not seeded** (so that any escape at a Floquet-stable energy is the
off-diagonal transient's doing, not the diagonal ghost's own tongue), tmax 3000, 7963 s.

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
- At δ = 10⁻⁴ the linear bands predict the escapes essentially perfectly (`RH_overnight_2026-10-22.md`); at δ = 10⁻²
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

Hole closed. (The δ-crossover scan, 10⁻⁴ → 10⁻² in half-decades at two energies, would pin the ~10⁻³; one hour.)
