# The reduced block's open items (ii)–(iv), worked 24 Oct 2026: the transfer matrix is not an invariant, the Krein width is not first order, a symplectic frame exists

Scripts `code/rh_reduced_block_transfer.py` (items (ii) and (iv)) and `code/rh_krein_width_cj.py` (item (iii)); logs
`results/rh_reduced_block_transfer_cj{1.0,0.5,0.25}_a0.13.txt`, `_cj1.0_a0.10.txt` (with `.npz`, named `_cj1_a0.13`, `_cj1_a0.1`, `_cj0.5_a0.13`, `_cj0.25_a0.13`), and
`results/rh_krein_width_cj.{txt,npz}`. Everything is computed on the exact 4-dim gauge-reduced block of
`RH_reduced_block.md` (moving gauge, colour-rotation quotient), whose checks pass here as there: every segment map
symplectic for Ω_red to ≤ 2 × 10⁻¹², the product of the four segments equal to the four ghost multipliers to
≤ 3 × 10⁻¹⁰, and the half-period map N = S̃Φ_{T/2} on the surface (S̃ maps S(t + T/2) onto S(t) to 10⁻¹⁰) squaring to
the monodromy to ≤ 3 × 10⁻¹⁰, at every E and c_j.

## (ii) The arm / bounce transfer matrix — the question has no frame-free answer

Segments: ARM = [−aT, aT] around a turning point (the light mode's hyperbolic window at E ≥ 0.3 is |t/T| ≤ 0.13,
`rh_reduced_block_hyp.txt`), BOUNCE = [aT, (½ − a)T], a = 0.13 (and 0.10 as a check). Each segment map is expressed
between Darboux frames of Ω_red at its ends (symplectic Gram–Schmidt, as `rh3_segments.py`), and its largest
singular value read as the segment's gain; the "return" is the fraction of the arm's ln-stretch that its
most-stretched output keeps after the bounce.

| E (c_j = 1) | arm ln gain, a = 0.13 / 0.10 | bounce ln gain | return | ln\|μ\|/2 (exact, per half period) |
|---|---|---|---|---|
| 0.1 | +0.40 / +0.71 | +0.19 / +0.44 | −0.39 / −0.55 | 0 |
| 0.156 (Krein) | +0.17 / +0.42 | +0.44 / +0.74 | −0.78 / +0.72 | +0.13 |
| 0.3 (stable) | +0.42 / +1.49 | +0.70 / +1.92 | +0.44 / −0.01 | 0 |
| 0.412 (Krein) | +0.53 / +1.99 | +1.07 / +2.60 | +0.09 / −0.10 | +0.53 |
| 0.6 | +1.39 / +2.73 | +1.65 / +2.35 | +0.04 / −0.34 | +0.38 (ln 2.132 / 2, multiplier in `_cj1.0_a0.13.txt`; stated as +0.19 until 28 Oct 2026) |
| 0.642 (parametric) | +1.35 / +2.51 | +1.59 / +1.63 | +0.12 / −0.14 | +0.45 |

Moving the cut from a = 0.13 to 0.10 changes the arm's gain by up to a factor e^1.46 (ln gain +0.53 → +1.99 at E = 0.412; stated as "e³" until
28 Oct 2026) and flips the sign of the return at five of six energies (all but E = 0.1; stated as four until 28 Oct 2026). That is not a sensitivity to where the hyperbolic window ends; it is the frame. **A Darboux
frame is defined up to Sp(4), and the singular values of a symplectic map are invariant only under its unitary
subgroup U(2).** The script's REFRAME check makes it explicit: with three random symplectic frames at the three
boundaries of a half period (at E = 0.3, c_j = 1: arm ln gain 0.42 → 0.98, 3.38, 2.55; return +0.44 → +0.65, −0.31,
+0.18) the segment gains and the return take any value, while the spectrum of N is unchanged. The same applies to
the per-segment growth of the unstable eigenvector of `RH3_segments.md` (it noted the frame convention for vectors
but called the singular values frame-independent; they are not). What is invariant on a time-periodic block is the
monodromy, hence N and its eigenvalues (Δ(E) among them); any split of a period into "what the arm does" and "what the
bounce returns" needs a frame convention, and the two conventions used so far (Gram–Schmidt from the unit vectors at
a = 0.125 in RH3, and here at 0.13 and 0.10) disagree at the level of the sign. RH3's localisation statement
("the growth belongs to the bounce" in the strong Krein band) is therefore a statement in RH3's frame, not a
property of the block; the frame-free content is the energy dependence of \|μ\| and Δ. Item (ii) closes as ill-posed
in the form asked; a well-posed version would fix the frame physically (e.g. the frozen normal modes of a chosen
symplectic frame, normalised iEᴴΩE = ±1, as the Bianchi IX transfer matrices of `bianchi9/README.md` do), and item
(iv) below shows that the frozen modes themselves depend on the frame.

**Δ(E) against the sector mixing** (the conjecture of `RH27_averaging.md`, part E, that Δ is the inter-sector mixing
of Φ_{T/2}). With P± = (1 ± S̃)/2 restricted to the surface and the mixing of a surface map A defined as
‖A − (W₊AW₊ + W₋AW₋)‖_F/‖A‖_F (Euclidean, the convention of hole 6's f₊), the mixing of Φ_{T/2} is 0.51–0.77 at every
(E, c_j) run, varying by ±0.1 across the whole table, while Δ runs from 0.09π to 0.77π (c_j = 1; the convention
of `rh5_splitting.py` reproduced: Δ = 0.161π at E = 0.3 against its 0.164π at 0.3031; off-circle rows are
marked and their Δ not meaningful). The surface itself is not sector-pure (the colour-rotation mode z = (−y, x, …)
and the gauge direction mix u12 and u21), so the sector decomposition on the reduced block is a convention too;
in this convention Δ and the mixing are unrelated, and the conjecture is withdrawn. The exact fact that survives is the
one RH27 started from: if Φ_{T/2} were sector-diagonal, N would have eigenvalues ±√(ab) and Δ = 0.

## (iii) The Krein band at [0.151, 0.163] across c_j — not first order

`rh_krein_width_cj.py`: block 12 monodromy on 71 energies in [0.12, 0.19] (ln-E step 0.0066) at c_j = 0.25, 0.5,
0.75, 1. The two-mode normal form of RH27 part B predicts, at first order in the ghost–background coupling, a band
of half-width 2\|g\| in detuning and peak rate \|g\| with g ∝ c_j, at a centre fixed by the unperturbed frequencies.

| c_j | Krein band in the window | width in ln E / in rungs (T_c/4.057) | peak rate (per unit time) | peak \|μ\| |
|---|---|---|---|---|
| 1 | [0.1520, 0.1612] | 0.059 / 0.15 | 0.0061 | 1.29 |
| 0.75 | none; +1 folds at [0.127, 0.144] (\|μ\| 1.77) and [0.163, 0.168] | — | — | — |
| 0.5 | [0.1316, 0.1377] | 0.046 / 0.12 | 0.0053 | 1.26 |
| 0.25 | [0.1490, 0.1634] | 0.092 / 0.24 | 0.0122 | 1.67 |

Neither the width nor the rate is monotone in c_j: the band is 1.6× wider and 2× stronger at c_j = 0.25 than at 1,
has moved to 0.132–0.138 at c_j = 0.5, and at c_j = 0.75 is replaced by a parametric fold of \|μ\| 1.77 where the two
pairs no longer meet. The first-order test fails; the low Krein band is not a weakly coupled resonance of two
fixed modes dressed by c_j — its two pairs' positions (Δ(E; c_j), watch-list item 9) move with c_j as fast as the
coupling does. This is the last of RH27's first-order predictions and it goes the way of the others.

## (iv) A symplectic frame along the orbit — exists, exactly Hamiltonian, and the frozen spectrum in it differs

Frame: in the polar-transported orthonormal basis B(t) of `run_4dim`, Ω_red(t) has two invariant 2-planes (the
eigen-planes of −Ω_red², strengths w₁ ≥ w₂; their ratio stays ≥ 0.095 on every orbit run, so the planes never
cross); e_k the unit vector of plane k rotated as little as possible from the previous sample, f_k = −Ω_red e_k/w_k²,
D = [e₁, e₂, f₁, f₂] gives DᵀΩ_red D = J to 10⁻¹⁵. The generator in this frame, K_z = D⁻¹(K_red D − Ḋ) (Ḋ by central
differences on 2000 samples), has J K_z symmetric to a median 5 × 10⁻⁵ – 2 × 10⁻⁴ (max 10⁻² where the frame turns
fastest), i.e. it is Hamiltonian, and its frozen spectrum is of the three Hamiltonian kinds only (imaginary pairs,
real pairs, quadruplets; the "non-Hamiltonian fuzz" — complex pairs with a real part — drops from 0.3–1.6 in the
orthonormal frame to ≤ 0.07 here, except at E = 0.642 where the frame's own rotation produces a quadruplet).
The remaining freedom is a rotation within each plane, fixed by minimal rotation; the connection term is not small
(‖Ḋ‖/‖K_red D‖ ≈ 0.2–0.5), which is the point: the frozen spectrum changes with it.

| E (c_j = 1) | hyperbolic fraction, orthonormal frame (`rh_reduced_block_hyp.txt`) | symplectic frame |
|---|---|---|
| 0.1 | 0.00 | 0.00 |
| 0.3 | 0.23 | 0.29 |
| 0.412 | 0.37 | 0.51 |
| 0.6 | 0.44 | 0.73 |
| 0.642 | 0.45 | 0.67 |

(c_j = 0.5: 0.24 → 0.41 at E = 0.3, 0.46 → 0.68 at 0.6; c_j = 0.25: 0.29 → 0.31, 0.48 → 0.52.) A time-dependent
symplectic change of frame D(t) sends K → D⁻¹KD − D⁻¹Ḋ, and the second term is an arbitrary Hamiltonian matrix, so no
frozen quantity of a time-dependent block is frame-free; the orthonormal frame of `RH_reduced_block.md` is one
convention (and not a Hamiltonian one), this one another. The qualitative picture holds in both — soft and
hyperbolic near the turning points at E ≥ 0.3, nothing at E = 0.1 — but the numbers 22 % / 44 % (and 29 % / 48 % at
c_j = 0.25) quoted in watch-list items 2 and 8 are convention-laden at the 30–60 % level and should be quoted with
the frame named, or as "a fifth to three quarters of the period depending on the frame".

## What this closes

- (ii) closed, negative: segment transfer matrices and "how much the bounce returns" are frame conventions on a
  periodic block; only N is invariant. RH3's per-segment numbers and its "frame-independent singular values" need
  that caveat. The Δ-equals-mixing conjecture of RH27 part E is withdrawn.
- (iii) closed, negative: the Krein band's width and rate are not ∝ c_j and the band does not persist across c_j.
- (iv) done: a symplectic frame with an exactly Hamiltonian frozen generator; its use is to show that frozen
  hyperbolic fractions are frame-dependent, not to replace the earlier numbers.
- For the rewrite (`FOR_THE_REWRITE_2026-10-25.md`): items 2 and 8 should not quote hyperbolic fractions without a
  frame; RH3's "growth belongs to the bounce" is frame-bound and should be dropped or stated as such; nothing
  first-order about the bands should be claimed.
