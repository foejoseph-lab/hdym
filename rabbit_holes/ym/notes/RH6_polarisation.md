# Rabbit hole 6 — which polarisation carries which sign (a negative result, with the reason)

21 Oct 2026. Script `code/rh6_polarisation.py` (100 energies, 68 s), log
`results/rh6_polarisation.txt`, data `results/rh6_polarisation.npz`. Scope: block 12 on the planar DR orbit,
c_j = 1, E ∈ [0.05, 1]. Float.

## The question, sharpened

The question as posed (`../README.md`, hole 6) is which of the two ghost-sector pairs is "the ghost's" and whether that changes with E, expecting
a projection onto the transverse / longitudinal polarisations of the massive Lee–Wick vector. Two things stand
in the way, one practical and one structural.

Practical: the repository has no polarisation basis for the off-diagonal components — `hdym_rung1_lyap.py`'s
docstring names "the longitudinal polarisation" without defining it, and in the homogeneous sector (zero spatial
momentum) transverse/longitudinal has no momentum to refer to. So three measurable contents were used instead:
the component split f12 = |u12|²/(|u12|² + |u21|²), the reflection-sector split f+ for u± = (u12 ± u21)/√2, and a
local frequency w² = |u″|/|u| (the Pais–Uhlenbeck heavy/light split read off the eigenvector), each at the
turning point t = 0 and, propagated, at t = T/4.

Structural: "the ghost's pair" means the pair whose positive-frequency member carries negative Krein sign
(sign = energy × frequency). B005's sign is attached to the member at positive *half-period angle*, which is the
frequency sign only modulo 2π: the half-period phase is ≈ 0.77 × T/2 ≈ 10–20 rad, so which member is at positive
frequency is exactly the rotation-number question hole 2 could not settle (`RH2_ladder.md`). That is why B005's
"sign reassignments" across tongues exist: the label flips when an angle crosses ±π, not the physics.

## What the contents say

- **At the turning point both pairs are an exact 50/50 mixture of the two reflection sectors** (f+ = 0.50 ±
  0.02 at every on-circle energy, both pairs). This is the symmetry: at y = 0 the eigenvectors of N = S̃Φ_{T/2}
  carry equal weight in u+ and u−. No sector label is available at t = 0.
- A quarter period later the pairs often sit in opposite sectors (E = 0.05: f+ = 0.81 vs 0.19; 0.21: 0.69 vs
  0.07; 0.31: 0.57 vs 0.20) but not always (0.149: 0.14 vs 0.64; 0.44, inside the Krein band: 0.35 vs 0.34), and
  the component split f12 moves between 0.02 and 0.99 with E. The contents reshuffle along the orbit and with
  energy; there is no E-independent assignment of a pair to a component or a sector.
- The heavy/light split is not clean either: w² at t = 0 is 0.93 [0.27, 1.98] for the smaller-angle pair and
  1.13 [0.74, 2.90] for the larger, overlapping; at T/4 both are ≈ 0.65. Neither pair is "the light YM-like
  mode" and neither is "the massive ghost": both carry order-one |u″|/|u|, as the frozen K at the bounce
  (frequencies 1.85 and 2.83 and a real pair ±1.32) already said.
- The Krein sign of the positive-angle member is + for 59–61 % of on-circle energies for either pair — i.e.
  uninformative, as the branch argument predicts.

## Reading

This is what hole 7's decoupled picture says it should be: the two pairs are *the same mode in the two
reflection sectors*, dressed by an O(1) cross-coupling. There is no "ghost pair" and "other pair" — each pair
is a ghost pair of its sector, and what distinguishes them is a phase splitting (0.2 rung at low E, a rung at
E ≈ 1), not a polarisation. The question the paper can answer is therefore not "which mode is the ghost's" but
"the two sector ghosts, degenerate by the reflection symmetry up to the splitting δ(E)"; δ(E) is the quantity
to tabulate (from the two copies of each parametric rung, hole 2).

If the polarisation question is still wanted, it needs (i) a definition of the basis on the diagonal
background (the natural candidates are the gauge/rotation zero-mode directions (−q₂, q₁) and (q₁, −q₂) in
(u12, u21) and their Ω-complements), and (ii) the rotation number, to fix the frequency sign. Neither was
attempted here.

## Status

Negative result, recorded; nothing for the paper beyond the sentence above. Hole closed unless the rotation
number gets done.
