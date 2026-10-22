# Rabbit holes 2 + 7 — first-order averaging on the canonical block, the sympy part done (21 Oct 2026, late)

Script `code/rh27_averaging.py` (sympy for the canonical form and the normal form, numpy for the scan; 120
energies, 12 s), log `results/rh27_averaging.txt`, figure `results/rh27_averaging.png`. Written to be rerun;
nothing here is at size.

## A. Canonical form (sympy, exact)

From L2 of `ym/rung1_system.pkl` restricted to (u12, u21), Ostrogradsky gives z = (Q1, Q2, P1, P2) = Jx(q, p) x and
a quadratic Hamiltonian H = ½ zᵀ Hm(q, p) z, 30 nonzero entries, with the canonical J. Checked along the orbit at
E = 0.1, 0.3, 0.6: J Hm = (dJx/dt + Jx K) Jx⁻¹ to 2.6 × 10⁻¹⁶ — the canonical equations are exactly the K the
Floquet scans integrate. (This is the symmetric, time-dependent Hamiltonian the paper's off-diagonal section can
quote; krein_b005's Ω = Jxᵀ J Jx is its pull-back.)

## B. Two-mode resonance normal form (sympy, exact)

Modes with signed frequencies ν = σω (σ the energy sign) and a resonant coupling 2 Re(g e^{inΩt} ζ_a ζ_b): in the
frame rotating at nΩ/2 the Hamiltonian is constant and sympy's eigenvalues, with Δ = ν_a + ν_b − nΩ, reduce to
λ = ±(i/2)(Δ ∓ √(Δ² − 4|g|²)) — purely imaginary for |Δ| > 2|g|, real part ±√(|g|² − Δ²/4) inside: the
**sum resonance of signed frequencies is unstable in a band of half-width 2|g| in detuning**. The difference
coupling 2 Re(g e^{inΩt} ζ_a ζ̄_b) gives λ = ±(i/2)(Δ ∓ √(Δ² + 4|g|²)): always stable. Parametric: a = b,
2ν = nΩ. Two negative-energy modes therefore collide at ω_a + ω_b = |n|Ω (n < 0 in the convention e^{inΩt}), and
a negative-energy mode with a positive one at ω₊ − ω₋ = nΩ — which is the Krein-signature statement in formula
form.

## C. What the averaged block says (numpy)

The time-averaged H₀ commutes with S̃ (4 × 10⁻¹²), so it is block-diagonal in the *component* basis (u12, u21)
— the cross-coupling has odd harmonics only and averages to zero — and each component carries one oscillatory
mode and one hyperbolic pair (±0.15–0.38: the gauge/rotation sector, which the average does not keep on the
constraint surface; dropped). Both oscillatory modes are **negative-energy** (ghost-like), ω₀ = 1.05 → 1.33 and
ω₁ = 1.50 → 2.69 over E = 0.04 → 1. They are not degenerate: the half-period symmetry maps each component to
itself (u12 → u12, u21 → −u21), so the two components are *different* T/2-periodic systems; the "same mode in two
sectors" picture of RH6/RH7 was the ± basis, and the averaged system says the natural basis is the component one.

**Mode 0's parametric rungs reproduce the whole measured ladder.** 2ν₀ = nΩ, i.e. ω₀(E) T(E) = |n| π:

| n | predicted E | type if odd n is read as "Krein" | measured block-12 band (B001) |
|---|---|---|---|
| 18 | 0.0580 | parametric | [0.0577, 0.0594] P — inside |
| 17 | 0.0759 | Krein | [0.0761, 0.0773] K — on the edge |
| 16 | 0.1020 | parametric | [0.1020, 0.1133] P — inside |
| 15 | 0.1416 | Krein | [0.1519, 0.1613] K — 7 % low |
| 14 | 0.2052 | parametric | [0.2228, 0.2648] P — 8 % low |
| 13 | 0.3159 | Krein | [0.2811, 0.4937] K — inside |
| 12 | 0.5328 | parametric | [0.5694, 0.7241] P — 6 % low |
| 19 | 0.0452 | Krein | the rung-14 Krein band of `rh2_predict.txt` at 0.0427–0.0435 — 5 % high |

Eight for eight, positions within 8 % in E, and the alternation of types comes out as **even n = parametric,
odd n = period doubling of mode 0** — which the Lemma of RH7 forbids for a single pair and which the true system
realises as the two-pair Krein collision. So the ladder's constant "ω_eff = 0.774" was never a frequency: the
mode's frequency is ω₀ ≈ 1.1–1.3 (the Lee–Wick mass dressed by the background), rising with E, and T_c linear
in n is what ω₀(E) T(E) = nπ looks like when ω₀ rises more slowly than 1/T. The combination resonances of modes
0 and 1 (ω₀ + ω₁ = |n|Ω) do not line up with any band, and their first-order couplings are mostly ≈ 0.

What first order does **not** give: the widths (|g| ≈ 10⁻⁴–10⁻³ against measured widths of 0.05–1.2 rung —
the instability is not driven by a single harmonic of H at first order) and the mechanism of the odd rungs,
which first-order theory calls single-mode period doubling. The next sympy step is clear and bounded:
second-order degenerate perturbation theory on the 4-dim subspace {mode 0 of u12, mode 0's partner in u21}
at an odd rung, with the cross-coupling's odd harmonics (k = 13, 15 at E = 0.3–0.6, i.e. the resonant ones)
entering twice — that is where the Krein quadruplet, its width and the collision angle's drift from ±i
should come from. It is also where the exact statement for "collisions only at −1" would be proved or
bounded. One sitting with sympy, no runs.

## Status

Hole 2 is now a formula with a mechanism: ω₀(E) from the averaged canonical block, rungs at ω₀T = nπ, odd
rungs Krein. Hole 7's third rule has its model (mode-0 period doubling dressed by the odd-harmonic cross-
coupling), pending the second-order calculation. The paper sentence: "the band centres of block 12 are the
resonances ω₀(E) T(E) = nπ of the ghost mode of the time-averaged block, parametric at even n and, by the
half-period symmetry, Krein at odd n".

---

## Review (21 Oct, 20:50–21:10): an outside reader's objections, and what the data said

An AI reader with only the four-paragraph summary (no code, no paper) was asked to judge the "odd n = Krein"
claim. Its objections, and the checks they prompted (`results/rh27_cj_scan.txt`, `rh27_mu_scan.txt`,
`rh5_splitting.txt`, and the scratch nulls recorded here):

1. **"The parity is a theorem, not a prediction."** Accepted. Given an instability at every odd rung, the lemma
   makes it two-pair; the content is the existence of a band at every odd rung, and the one non-trivial
   observation is that no Krein collision happens at +1.
2. **"Within 8 % is not evidence."** Confirmed, worse than it guessed: the half-integer rungs ω₀T = (n + ½)π also
   fall within 8 % of a band 7 times in 9, because bands plus an 8 % halo cover 79 % of the log-E range; inside-
   band hits for the integer rungs are 4 of 9 against 44 % coverage. **Rung 11 lands at E = 1.049 inside the
   parametric band [1.014, 1.271]: the parity fails at the ninth rung.** What survives of the positional claim is
   the slope d(ω₀T)/dE and the count (nine rungs, nine bands): one number and an anchor.
3. **"Odd rungs are first order in the 0–1 coupling, not second order."** Accepted; the width comparison in part C
   was against the wrong matrix element. But the two-mode combination resonances with the averaged ω₁ do not
   locate the Krein bands either, and the averaged mode 0 tracks only one of the two measured pairs' phases
   (correlation 0.62, mean error 0.10 π over 400 energies) while the averaged mode 1 tracks nothing (−0.36 / 0.08).
   The time-averaged Hamiltonian is a correct zeroth order for one pair and not the other.
4. **Its test 4 (scale the component cross-coupling by μ).** Run, and invalid as posed: at μ < 1 the block grows by
   |μ| ~ 10²–10³ per period, but the unstable eigenvectors violate the linearised Gauss law at the 0.9 level — the
   deformation breaks the gauge structure (B002's F118 trap). The gauge-consistent dial is c_j itself (Gauss
   violation 10⁻¹² at every c_j): at **c_j = 0.25 the block is far more unstable than at c_j = 1** (|μ| up to 34, a
   Krein band over [0.41, 1.0]), at 0.5 still more than at 1 (|μ| 6.5 over [0.65, 1.0]), and the band positions move
   with c_j. So the c_j = 1 theory is the *least* unstable member of the family except c_j = 0 (pure YM, stable):
   the Krein bands at c_j = 1 are what remains after a strong cancellation, and nothing about them is
   perturbative in the ghost–background coupling. The reviewer's dichotomy ("pinned at the rung = kinematic,
   drifting = coupling-shifted") comes out on the drifting side.
5. **Its second note: "the pairs are mirror-related, φ_B ≈ π − φ_A in M".** Measured (`rh5_splitting.py`): the
   splitting Δ(E) = (φ_A + φ_B) mod 2π is **0.07–0.25 π over E = 0.05–0.35** — the pairs are near-degenerate as
   sets, the decoupled-symmetric-model structure of RH7 (N-angles summing to π), not mirror-related — grows
   through the wide Krein band to 0.5–0.8 π at E = 0.5–0.72, passes **exactly through 0 at E = 0.80** (B005's fake
   degeneracy: λ and −λ in N, no interaction), and reaches 0.95 π at E = 1, where the mirror relation does hold and
   the parametric bands alternate between the pairs as the reviewer predicted (rung 11). This curve is hole 5's
   δ(E), the quantity RH6 said to tabulate; it is also why the ladder works below E ≈ 0.35 and breaks above:
   the rungs locate the bands only while the pairs are degenerate to a fraction of a rung.
6. **Signatures at the collisions.** Already in B005: at every true meeting the colliding members carry opposite
   Krein sign (ψ₁⁺ = −ψ₂⁺), so the reviewer's "either B is positive-energy or S̃ reverses an orientation" is
   answered by measurement; the averaged system calls both modes negative-energy, consistent with the second
   branch.

**Standing after review.** Exact: the lemma, the canonical form, the normal form. Measured: the ladder in T
(slope and count), δ(E), the Krein signatures, the c_j family. Withdrawn: "mode 0's rungs predict the bands" as
more than slope-plus-count, and "odd n = Krein" as a prediction. The open problem is now stated correctly: the
right zeroth order is the gauge-reduced four-dimensional block (two ghost pairs, no constraint sector), whose
two frequencies are near-degenerate below E ≈ 0.35 with the splitting δ(E) above; averaging the full
eight-dimensional canonical form does not produce it. Sympy can build the reduced system from part A plus the
constraint; that is the next sitting, before any more resonance theory.

### Review, round 3 (21:15–21:35)

7. **Gauss residual of the averaged modes** (t = 0; E = 0.1 / 0.3 / 0.6): hyperbolic pairs 0.5–1.2 (constraint-
   contaminated, as the reviewer predicted); averaged mode 1 (1.7–2.4) **clean**, 0.01–0.02; averaged mode 0
   (1.1–1.25) the **dirty** oscillatory one, 0.04 → 0.23 with E. The mode that tracked a measured pair was the
   contaminated one. Orbit-averaged residual not yet computed.
8. **Signature-oriented unwrapping** (`rh5_unwrap.py`, the reviewer's recipe): below E ≈ 0.22 both pairs advance at
   the same rate — effective ω_A ≈ 0.78 / 0.76 / 0.82 and ω_B ≈ 0.83 / 0.79 / 0.81 in three windows — and split above
   0.27 where Δ(E) grows. Positive-member identity switches at fold re-entries (E = 0.056, 0.109, 0.267 for A;
   0.248, 0.730 for B), read as within-pair orientation flips.
9. **The branch** (`rh5_branch.py`, revised): the first version FFT'd x(t) over one period with a rectangular window —
   Joseph's spectral-leakage objection: e^{iωt}p(t) is not periodic over [0, T), and the phase jump at the wrap leaks
   a single line into 1/k sidebands, which is what the "sidebands to ±8Ω" of the first run were. Revised: divide out
   the Floquet phase (principal branch), FFT the periodic part (periodicity residual ≤ 2 × 10⁻²). Result: the spectra
   are genuinely **two-line** — a line at |f| ≈ 0.85–1.2 (0.93 at E = 0.05, 1.1 at 0.1–0.3, 0.84–0.87 at 0.36–0.55,
   1.15–1.2 at E ≥ 0.85; the dressed Lee–Wick mass) and a quasi-static line at |f| ≲ 0.1, each with 20–50 % of the
   power, plus a weaker line at ±(1–2)Ω; nothing near 1.7–2.4. Both pairs carry the same lines. So the physical
   frequency of the doublet is ≈ 1, E-dependent at the 20 % level, riding on a large quasi-static component; the
   slope-derived 0.8 of item 8 was a branch artefact, and the dirty averaged mode's 1.1–1.25 was close to the real
   line after all.
10. **c_j = 0.25, corrected** (`rh5_splitting_cj.py`, and the classifier of `rh27_mu_scan.py` fixed: four real
    eigenvalues off the circle are a double +1 fold, not a collision): Δ is 0.47–0.57 π where the pairs are on the
    circle, two narrow Krein bands survive ([0.151, 0.163], [0.304, 0.321] — the first at the same place as at
    c_j = 1), and the giant band [0.41, 1.0] is a **double +1 fold** whose exponent is smooth and single-humped across
    the band (ln|μ| 2.0 → 3.5 → 2.3, no rung structure): a static instability, ω_eff² < 0, as the reviewer
    suggested. So Δ controls whether collisions occur; the c_j stability ordering is a static +1 phenomenon the
    reduced block can address as ω_eff²(c_j) at fixed E.

**Target for the reduced block, as it now stands:** two ghost pairs with multipliers degenerate (Δ ≲ 0.25 π)
below E ≈ 0.22, each a line near 1 (0.85–1.2) plus a comparable quasi-static line; no mode at 1.7–2.4; ω_eff²(c_j) with a negative
minimum near c_j ≈ 0.3–0.4 and marginal at c_j = 1. Checks against ωT mod 2π, not ω.
