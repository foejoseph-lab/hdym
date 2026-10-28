# The gauge-reduced 4-dim ghost block of block 12 — built and verified on the Floquet side (22 Oct 2026)

Script `code/rh_reduced_block.py` (numerical; the Ostrogradsky form, the Gauss row and their orbit derivatives are
exact sympy objects, the reduction and integration are float). Logs `results/rh_reduced_block_cj1.txt` (E = 0.1,
0.3, 0.6; the E = 0.3 row rerun with the corrected Ω in `_cj1_E03.txt`), `_cj025.txt`, `_cj05.txt`; frozen spectra in
`results/rh_reduced_block_E*_cj*.npz`. This is the zeroth order `RH27_averaging.md` ended by asking for.

## Construction

Block 12 (8-dim state x = (u12, u21, u′, u″, u‴)) carries one linearised Gauss row g(t)·x = 0 (first class: a
random Gauss-clean initial condition keeps g·x at 10⁻¹⁵ over a period). With Ω(t) the Ostrogradsky symplectic form
of the block (c_j kept symbolic — the c_j = 1 form gives a wrong gauge direction at c_j ≠ 1, found by a wrong first result at c_j ≠ 1), the
gauge direction is v_G = Ω⁻¹g.

**Constant gauge fixing is not admissible on the DR orbit.** The unit gauge direction v_G(t) sweeps the
whole 8-space over a period (singular values of its samples 0.74, 0.39, 0.36, 0.30, 0.18, 0.16, 0.09, 0.09 at E = 0.1; from an earlier run, the committed logs
print only a fixed string "0.74/0.39/.../0.09"),
and the best fixed functional f has min |f·v_G| ≈ 10⁻³ somewhere on the orbit, so the gauge-fixed flow x′ = Kx + λv_G
with λ = −f·Kx/f·v_G blows up there (a multiplier 4 × 10⁸ in the first attempt; no committed log). The moving gauge f(t) = v_G(t) is
used instead, with the connection term from the exact orbit derivatives: λ = −(v̇_G·x + v_G·Kx)/|v_G|², v̇_G =
Ω⁻¹(ġ − Ω̇v_G). Off-surface components g·x and v_G·x are then conserved (10⁻¹¹ over a period) and the surface
S(t) = {g·x = 0, v_G·x = 0} is invariant (leak 10⁻¹⁴).

**Second reduction.** The two leftover dimensions are the global colour-rotation zero mode z(t) = (−y, x, −p_y, p_x,
…) (an exact solution, Gauss-clean to 10⁻¹⁸) and its conserved charge c(t) = Ωz. Quotienting by z (made orthogonal to
v_G) and restricting to c·x = 0, with a second moving gauge term λ₂z, gives the 4-dim block.

## Checks (every E, every c_j in {0.25, 0.5, 1})

- CHECK A (6-dim): the four ghost multipliers of the full block are reduced multipliers to ≤ 8 × 10⁻¹¹; the other two
  are the unipotent pair at 1.
- CHECK A4 (4-dim): the reduced monodromy's spectrum *is* the four ghost multipliers, to ≤ 8 × 10⁻¹¹, including the
  c_j = 0.25, E = 0.6 case with |μ| = 34 and 16.5 (the double +1 fold).
- The reduced 4 × 4 Ω is nondegenerate and the reduced monodromy is symplectic for it to 3 × 10⁻¹³.

So the 4-dim reduced block is exact (not a model): it is the ghost dynamics with the constraint and symmetry sectors
removed, and anything computed from it is a statement about the theory.

## First outputs

**Frozen spectrum along the orbit** (eigenvalues of the restricted generator at fixed t; the gauge terms make it not
exactly Hamiltonian, real parts ≤ 0.1 at the printed times and up to 0.9 elsewhere, so read the imaginary parts as instantaneous frequencies; a
symplectic frame in which it is exactly Hamiltonian, and in which these numbers change, is in `RH_reduced_block_transfer.md`), c_j = 1:

| E | turning point t = 0 | bounce t = T/4 | t = T/2 |
|---|---|---|---|
| 0.1 | 0.71, 1.70 | 0.81, 1.44 | 0.71, 1.70 |
| 0.3 | 0.34, 1.29 | 0.62, 1.57 | 0.33, 1.29 |
| 0.6 | 0.13, 1.02 | 0.73, 1.64 | real pair ±0.17 and 1.02 |

Two modes: a light one that goes soft at the turning points (0.13–0.7, hyperbolic at E = 0.6, t = T/2) and a heavy
one at 1.0–1.7. This is the sweep of `rh5_instfreq.png` seen from the Hamiltonian side, with one disagreement with
that reading, not resolved here (resolved 23 Oct: `RH27_averaging.md` item 10): the frozen spectrum puts the light mode's softening at the **turning points** (x extremal,
y = 0), not at the bounce (x = 0), where both modes are stiffest, while the instantaneous-frequency estimator of
`rh5_instfreq.py` puts the dips at the bounce (`RH27_averaging.md`, item 10). The quasi-static spectral line is the light mode; the main line near 1 is the heavy mode dressed
down by the sweep.

**c_j family** (frozen, E = 0.6, t = T/2): c_j = 1 real pair ±0.17; c_j = 0.5 real ±0.06/0.23; c_j = 0.25 real
±0.12/0.29 — the light mode is hyperbolic at the turning points for every c_j, and the static instability of the
c_j = 0.25 family is **not** a mean-stiffness effect: the averaged generator at c_j = 0.25 has purely imaginary
spectrum (0.59, 0.66 at E = 0.6; 0.82, 0.85 at 0.3; 0.98, 1.07 at 0.1) while the exact multipliers are 34 and 16.5.
Of the two alternatives — negative mean stiffness, or transient-negative stiffness (a Kapitza-type effect the averaged
block cannot see) — only the second is compatible with the data, and it is not the whole story. The hyperbolic
fraction of the light mode and its strength both grow as c_j falls (`rh_reduced_block_hyp.txt`; the criterion is
|Im| < 10⁻³ and |Re| > 0.1, the gauge terms' real parts staying below 0.1 where the spectrum is elliptic):

| E | c_j | hyperbolic fraction | max Re (from the `_E*_cj*.npz` spectra; not printed in the hyp log) | ∫Re dt over the hyperbolic stretches | measured ln\|μ\| per period |
|---|---|---|---|---|---|
| 0.1 | 0.25 / 1 | 0 / 0 | 0.30 / 0.26 | 0 / 0 | 0 / 0 |
| 0.3 | 0.25 / 0.5 / 1 | 0.29 / 0.24 / 0.22 | 2.03 / 1.52 / 0.90 | 11.6 / 7.4 / 5.1 | 0 / 0 / 0 (stable at all three) |
| 0.6 | 0.25 / 0.5 / 1 | 0.48 / 0.47 / 0.44 | 4.04 / 3.74 / 2.84 | 27.2 / 24.5 / 17.3 | 3.5 / ≈ 0 / 0.76 (ln 2.132, `_cj1.txt`; stated as 0.9 until 28 Oct 2026) |

So the frozen hyperbolicity is an order of magnitude too large to be the growth: at E = 0.3 the generator is
hyperbolic for a fifth to a third of every period, with ∫Re dt = 5–12, and the block is stable at every c_j; at
E = 0.6, c_j = 0.5, the integral is 24 and the block is (just) stable. What the frozen stretches expand, the rest of
the period undoes — RH3's lesson (large per-segment stretching does not decide stability) at the level of the
generator. The c_j ordering of the static instability does correlate with the hyperbolic strength (max Re 0.90 →
2.03 at E = 0.3, 2.84 → 4.04 at 0.6), but the mechanism is in how much of the arm's expansion survives the bounce,
which is the transfer matrix (open item (ii) below), not the frozen spectrum.

**Averaged block (parallel-transported frame).** Not Floquet-accurate: at c_j = 1 the averaged
frequencies (1.16/1.61, 1.43/1.31, 1.12/1.50 at E = 0.1/0.3/0.6) give ωT mod 2π at 1.14/1.82, 0.05/0.73, 0.62/0.14 π
against the measured 0.10/0.26, 0.40/0.56, 0.00/0.77 π (no committed log). The modulation is too strong for first-order averaging even on
the exact reduced block; the frozen spectrum and the bounce transfer matrix are the objects to use.

## What this settles and what it opens

- Settled: the 4-dim block exists, is exact, and carries two modes — light (swept, soft at the turning points) and
  heavy (≈ 1–1.7). The "near 1" line of the spectra is the heavy mode; the quasi-static line the light one.
- Settled: the low-c_j static instability is not negative mean stiffness; the light mode is hyperbolic on the arms
  for 22–48 % of the period at E ≥ 0.3, more strongly at lower c_j, but the frozen expansion overshoots the measured
  growth by an order of magnitude and does not by itself decide stability (table above).
- Done (23 Oct): (i) the hyperbolic fraction and strength per (E, c_j), `results/rh_reduced_block_hyp.txt`.
- Closed 24 Oct (`RH_reduced_block_transfer.md`): (ii) is a frame convention (no invariant answer), (iii) is not
  first order, (iv) is built and shows the frozen spectrum above to be one convention among many — the hyperbolic
  fractions in the table are 0.29 / 0.51 / 0.73 at E = 0.3 / 0.412 / 0.6 in the symplectic frame. As posed on 22 Oct: (ii) the arm/turning-point transfer matrix: integrate the 4-dim
  reduced flow across one hyperbolic arm (t/T ∈ [−0.13, 0.13], say) and across one bounce, and read how much of the
  arm's expansion the bounce returns and how the sectors mix (`RH27_averaging.md`, part E); (iii) the Krein width at
  [0.151, 0.163] across c_j as the one first-order test; (iv) a proper Hamiltonian form H_red(t) = −Ω_red K_red in a
  symplectic frame (the current frame is orthonormal, which is why the frozen generator is not exactly Hamiltonian).
