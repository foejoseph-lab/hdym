# Rabbit hole 3 — where on the orbit the growth happens, answered on the exact reduced block (22 Oct 2026)

Script `code/rh3_segments.py`, log `results/rh3_segments.txt`. A first attempt on 19 Oct with segment transfer matrices
in the Euclidean norm of the 8-dim block was inconclusive (local rates ±3 around 0.03). Two things were
missing, and both are now in hand: the constraint and symmetry sectors are gone (the 4-dim reduced block of
`RH_reduced_block.md`), and the norm is canonical — each quarter-period segment's transfer matrix is expressed
between Darboux frames of the reduced symplectic form at its two ends (symplectic Gram–Schmidt; Dᵀ Ω_red D = J to
10⁻¹⁶), so its singular values are the frame-independent expansion of a symplectic map. Segments: TP1 = [−T/8, T/8]
and TP2 = [3T/8, 5T/8] around the two turning points (x extremal, y = 0, where the light reduced mode is soft), B1 =
[T/8, 3T/8] and B2 = [5T/8, 7T/8] around the two bounces (x = 0).

Checks: every segment map symplectic to ≤ 2 × 10⁻¹²; the product of the four reproduces the block's ghost
multipliers to ≤ 4 × 10⁻¹⁰; the half-period symmetry shows as TP1 = TP2 and B1 = B2 to all printed digits.

## Result

Per-segment ln growth of the **actual unstable Floquet mode** (real and imaginary representatives; the max-singular-
value gains, which are worst-case bounds, are in the log):

| E | band | TP1 | B1 | TP2 | B2 | net per period (ln|μ|) |
|---|---|---|---|---|---|---|
| 0.412 | Krein, |μ| = 2.87 | −0.95 / −1.11 | +1.40 / +1.72 | −1.10 / −0.94 | +1.58 / +1.51 | +1.05 |
| 0.642 | parametric, |μ| = 2.44 | +0.17 | +0.28 | +0.17 | +0.28 | +0.89 |
| 0.156 | Krein, |μ| = 1.29 | −0.11 / +0.17 | −0.00 / +0.20 | +0.18 / −0.12 | +0.18 / +0.02 | +0.25 |
| 0.30 | stable | (worst-case gains 1.9 / 2.2 per segment cancel exactly) | | | | 0 |

- **In the strong Krein band the instability belongs to the bounce**, in the strong sense: the unstable mode grows
  by e^{1.4–1.7} at each bounce and *shrinks* by e^{−1} at each turning point; the net factor 2.87 per period is the
  residue of bounce growth minus turning-point decay. "The instability belongs to the bounce" is a result at E = 0.412, not a hypothesis.
- **In the parametric band it is not localised at all**: growth is spread 0.17 / 0.28 / 0.17 / 0.28 over the four
  segments. A parametric resonance is a whole-period interference, as it should be.
- **In the weak Krein band at 0.156** the per-segment numbers are ±0.1–0.2 with cancellations and the net 0.25 is an
  interference residue; "where" has no sharp answer there.
- At the stable energy each segment still stretches its worst-case vector by 1.9–2.2 (ln 0.65–0.81, bounce larger
  than turning point), and the product is exactly elliptic: large per-segment stretching does not decide
  stability, and the bounce stretches more at every E, unstable or not. Segment gains alone (what the 19 Oct attempt
  measured) cannot answer the question; the eigenvector's per-segment growth can.

Caveat: the canonical norm between two *different* Darboux frames is not frame-invariant for a single vector (the
singular values are; a vector's growth is not), so the per-segment eigenvector numbers carry a frame convention —
the same convention at every boundary, which the TP1 = TP2 and B1 = B2 symmetry confirms is applied consistently.

## Bearing on the paper

Assumption A2 (channel entry) and B002's "the DR orbit's instabilities belong to the channel bounce, not the
adiabatic interior" — supported for the Krein band at its peak, by the strongest possible version (the turning
points *undo* growth); not applicable to the parametric band, which is delocalised. Hole closed.
