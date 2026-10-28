# Rabbit hole 7 — the selection rule of block 12, derived (what is exact and what is not)

21 Oct 2026. Script `code/rh7_selection_rule.py`,
log `results/rh7_selection_rule.txt` (150 energies, 184 s), figure `results/rh7_sector_tongues.png`.
Scope: the off-diagonal ghost block (u12, u21) on the planar DR orbit, c_j = 1, E ∈ [0.05, 1]. Float.

## What B005 measured

Single-mode folds only at M = +1 (the half-period map N = S̃ Φ_{T/2} at ±1), two-pair collisions only at
M = −1 (N at ≈ ±i), and a + 2 ≥ 0 for every real a = μ + 1/μ of a single pair. Measured at every E, not derived.

## Setting

x′ = K(t) x, K real and T-periodic, Hamiltonian for the Ostrogradsky form Ω(t) (B005 verified Mᵀ Ω M = Ω).
S̃ = diag(+1, −1) on (u12, u21) at every derivative order is a constant canonical involution (B005: Nᵀ Ω N = Ω)
with **S̃ K(t + T/2) S̃ = K(t)** (step 1 below: residual 2.6 × 10⁻¹¹). Hence N := S̃ Φ_{T/2} satisfies N² = M and is
Ω-symplectic (B005: both checked to 10⁻⁸).

## Lemma (exact): a single pair of M cannot leave the unit circle at −1

The spectrum of M is the set of squares of the spectrum of N, and the spectrum of a real symplectic N is closed
under λ → λ̄ and λ → 1/λ. Let μ be a multiplier of M with μ real, μ < 0, |μ| ≠ 1. Its square roots are ±i r with
r = √|μ| ≠ 1: non-real, off the circle. So one of them, say i r, is an eigenvalue of N and brings with it the
whole quadruplet {i r, −i r, i/r, −i/r}, whose squares are {−r², −r², −1/r², −1/r²}. Therefore **μ has
algebraic multiplicity at least two in M, and so does 1/μ**: a real negative multiplier of M off the circle is
never an isolated pair. Equivalently, for an isolated ghost pair a = μ + 1/μ can touch −2 but never cross it
(B005's "a + 2 ≥ 0", min 1.4 × 10⁻⁵), and a −1 instability of the block is necessarily a two-pair event: a Krein
quadruplet of N near ±i (the generic case) or its degenerate limit, a double real negative pair.

Corollary (same argument read the other way): a fold of a single pair of N at ±1 is a fold of M at +1 — the
period-doubling instability of a single pair, which a generic 2-dof Hill system shows, is **forbidden on any
background with a half-period reflection symmetry**. This is the derivation of two of the three measured rules.
It is the Swift–Wiesenfeld mechanism (PRL 52, 705, 1984: period doubling suppressed for orbits with a
half-period symmetry, M = N²) in its symplectic multiplicity form; see `RH7_literature.md` for what is theirs
and what is ours.

Check (step 3 of the script): all 16 LEAVE/RE-ENTER events of B005 are either single-pair folds of N at
|arg| < 0.05 π from {0, π} (ten events, all M = +1) or four-member collisions (six) (`rh7_selection_rule.txt`, step 3;
stated as eight and eight until 28 Oct 2026); no single pair ever
folds near ±i. Scanning B005's angles for single-pair passages through ±0.5 π (M = −1) outside the bands finds
exactly one, at E = 0.520 (pair k0/k3, with the other pair 0.27 π away: |ψ| 0.508 vs 0.781 π at E = 0.52038 of
`islands/results/krein_b005.txt`; stated as 0.17 π until 28 Oct 2026): the pair crosses M = −1 and nothing
happens, as the lemma requires. Every single-pair passage through 0 or ±π (M = +1) in the scan opens a tongue.

## Parity of the coupling (exact, and what it does and does not give)

Step 1 of the script: in the (u12, u21) basis the self-coupling blocks are T/2-periodic (only even harmonics
of Ω = 2π/T; the odd ones are 4 × 10⁻¹² relative) and the cross blocks flip sign over T/2 (only odd harmonics;
even ones 8 × 10⁻¹²). In the sector basis u± = (u12 ± u21)/√2 the same symmetry is the swap: K₊₊(t + T/2) =
K₋₋(t), K₊₋(t + T/2) = K₋₊(t) (residual 2 × 10⁻¹¹). This is the "self even, cross odd" parity statement conjectured in the 19 Oct review,
written out and verified.

The textbook reading — self-coupling at 2kΩ pumps a single mode at 2ω = 2kΩ, i.e. ωT = 2πk, M = +1; cross-
coupling at (2k + 1)Ω pumps the sum resonance ω₁ + ω₂ = (2k + 1)Ω — explains the *pairing of mechanisms with
harmonics* but does **not** by itself put the collisions at M = −1: a sum resonance fixes ω₁ + ω₂, not the angle
at which the two pairs meet. The collision lands at −1 only if the two pairs are degenerate in phase per period
(ω₁ T ≈ ω₂ T mod 2π), which is the content of hole 5's "near-degenerate pairs". The leading cross harmonics are
k = 13, 15 at E = 0.3–0.6, i.e. (ω₁ + ω₂)/Ω — the cross-coupling is resonant with the ghost pair, as it should be.

## The decoupled symmetric model (where the rule is exact, and why that is not enough)

Zeroing the sector cross blocks K₊₋, K₋₊ keeps the swap symmetry and makes the two sectors half-period
translates of each other: their period-T monodromies Φ₋Φ₊ and Φ₊Φ₋ are conjugate (same spectrum, checked to
2 × 10⁻⁵, the DOP853 tolerance on an 8 × 8 product; the log's CHECK line reports this as FAIL against its tol 10⁻⁷), every multiplier of M_dec is double, and N_dec has exactly
±√μ (checked to 7 × 10⁻⁷ over 150 energies). In this model a sector ghost pair at e^{±iθ} appears in N at
±θ/2 and π ∓ θ/2, with the member at π − θ/2 descending from e^{−iθ} (opposite Krein sign to the one at θ/2):
the two copies coincide only at N = ±1 (θ = 0, a same-pair fold, M = +1) and at N = ±i (θ = π, an opposite-sign
double eigenvalue, M = −1). So **in the decoupled symmetric model the rule is exact and its Krein bands are
the sector's period-doubling tongues**, dressed into complex quadruplets by any symmetric cross-coupling.

But the model is not quantitative. Its tongues (scan, step 2) alternate −1/+1 down the ladder like the real
bands, yet sit one rung off: decoupled fold at −1 in [0.245, 0.265] where the real block is parametric
(+1) in [0.223, 0.265]; decoupled fold at +1 in [0.447, 0.526] inside the real Krein band [0.281, 0.494];
decoupled fold at −1 from 0.83 to 1 where the real block is stable. The cross-coupling is O(1) (cross harmonic
amplitudes 2.8 against self 5.8 at E = 0.6) and shifts the per-period phase by about π. So the decoupled model
explains the *form* of the rule, not the positions, and the measured collision angles of N, ±(0.44–0.68) π
against the model's exactly ±0.5 π, are the size of the splitting the coupling produces.

## Status of the three rules

| rule (B005) | status after this hole |
|---|---|
| a + 2 ≥ 0; no single pair leaves at −1 | **derived, exact** (Lemma); checked on all 16 events and on the one safe passage at E = 0.520 |
| single-mode folds only at M = +1 | **derived, exact** (Corollary of the Lemma: a single pair can only fold at N = ±1) |
| collisions only at M = −1 | exact in the decoupled symmetric limit; in the real block it is "collisions at N = ±(0.44–0.68) π", i.e. the two pairs are degenerate in period-phase to 0.1–0.2 π at every collision (Δ(E), `RH27_averaging.md` item 5). Not a theorem; an explicit obstruction is the E > 1 scan of hole 2 (`results/rh2_predict.txt`), where the pairs have separated by a full rung and a parametric band sits at rung 6 where the alternation would put a Krein one. |

## What this changes in the paper (if anything)

One sentence in the off-diagonal section can now be a statement rather than an observation: *on a background
with a half-period reflection symmetry the ghost's period-doubling instability is forbidden for a single
pair; every multiplier at −1 is a two-pair Krein event, and every single-pair instability is a fold at +1
(Lemma, two lines from M = N²).* The "collisions at −1" half stays measured, with the N-angles quoted.

## Not done / next entrances

- Literature: done, `RH7_literature.md` (Swift–Wiesenfeld 1984; Lamb–Roberts 1998; Frauenfelder–Moreno 2023).
- The sector-basis picture suggests the natural object is not the decoupled model but the *dressed* sector:
  the cross-coupling's resonant harmonics k = 13, 15. A two-harmonic Mathieu-type reduction (self 2kΩ, cross
  (2k+1)Ω, two modes) would be the first model that could predict rung positions. Not attempted.
- Hole 5 (why block 12 and not 13, 23): blocks 13 and 23 have S̃ = +1 on both components (B005's STILDE), so
  for them N = Φ_{T/2} with no twist and the Lemma applies verbatim: no single-pair period doubling in any
  block (B001: "no pair ever leaves at −1" anywhere). What they lack is the two-pair event: hole 5's
  observation that their K is exactly block-diagonal on the planar orbit means their two pairs never couple,
  so no collision can open — only +1 folds, which is what B001 sees. Block 12 is the one whose two components
  mix under the orbit's reflection (S̃ = diag(+1, −1)), i.e. the one with a cross-coupling at all. That is the
  derivation of hole 5's "why" in one sentence; the non-planar question (do all blocks couple off the plane)
  is answered in `RH5_nonplanar.md`: on the x = y = z orbit none does.
