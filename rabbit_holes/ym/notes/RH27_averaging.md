# Rabbit holes 2 + 7 — the canonical block, first-order averaging, and what it took to find the right zeroth order

21 Oct 2026 (one evening, in four rounds); rewritten 23 Oct 2026 from the state at that date. Scripts
`code/rh27_averaging.py` (sympy canonical and normal forms; numpy scan of the averaged block, 120 energies, 12 s),
`code/rh27_mu_scan.py` (the μ and c_j scans), `code/rh5_splitting.py` and `rh5_splitting_cj.py` (Δ(E) at c_j = 1 and
0.25), `code/rh5_unwrap.py` (signature-oriented unwrapping), `code/rh5_branch.py` (spectra of the Floquet solutions),
`code/rh5_instfreq.py` (instantaneous frequency along the orbit). Logs and figures `results/rh27_*`, `results/rh5_*`.
Scope: block 12 on the planar DR orbit; c_j = 1 unless stated. Nothing here is at size.

**Where this ended.** The canonical form and the resonance normal form are exact and stand. The claim that
first-order averaging of the 8-dim block predicts the bands — "ω₀(E)T(E) = nπ, parametric at even n, Krein at
odd n" — was made and withdrawn the same evening (the null test below). What the evening established instead is
measured: the two pairs are near-degenerate at low E and antipodal at E = 1 (Δ(E)); their physical frequency is a
line near 1 with a quasi-static companion, and the oscillator is swept, soft once per half period; the c_j
family's stability ordering is a static +1 effect; and the right zeroth order is the gauge-reduced 4-dim block,
built the next day (`RH_reduced_block.md`). The paper watch-list carries the retired claims (`FOR_THE_REWRITE_2026-10-25.md`,
items 6, 8, 9, 10, section D).

## A. Canonical form (sympy, exact)

From L2 of `ym/rung1_system.pkl` restricted to (u12, u21), Ostrogradsky gives z = (Q1, Q2, P1, P2) = Jx(q, p) x and a
quadratic Hamiltonian H = ½ zᵀ Hm(q, p) z, 30 nonzero entries, with the canonical J. Checked along the orbit at
E = 0.1, 0.3, 0.6: J Hm = (dJx/dt + Jx K) Jx⁻¹ to 2.6 × 10⁻¹⁶ — the canonical equations are exactly the K the Floquet
scans integrate. krein_b005's Ω = Jxᵀ J Jx is its pull-back. This is the symmetric time-dependent Hamiltonian the
off-diagonal section can quote.

## B. Two-mode resonance normal form (sympy, exact)

Modes with signed frequencies ν = σω (σ the energy sign) and a resonant coupling 2 Re(g e^{inΩt} ζ_a ζ_b): in the
frame rotating at nΩ/2 the Hamiltonian is constant and its eigenvalues, with Δ = ν_a + ν_b − nΩ, are
λ = ±(i/2)(Δ ∓ √(Δ² − 4|g|²)) — purely imaginary for |Δ| > 2|g|, real part ±√(|g|² − Δ²/4) inside: the **sum
resonance of signed frequencies is unstable in a band of half-width 2|g| in detuning**. The difference coupling
2 Re(g e^{inΩt} ζ_a ζ̄_b) gives λ = ±(i/2)(Δ ∓ √(Δ² + 4|g|²)): always stable. Parametric: a = b, 2ν = nΩ. Two
negative-energy modes therefore collide at ω_a + ω_b = |n|Ω (n < 0 in the convention e^{inΩt}), and a
negative-energy mode with a positive one at ω₊ − ω₋ = nΩ — the Krein-signature statement in formula form.

## C. What the averaged 8-dim block said, and why it was withdrawn

The time-averaged H₀ commutes with S̃ (4 × 10⁻¹²), so it is block-diagonal in the component basis (u12, u21): the
cross-coupling has odd harmonics only and averages to zero. Each component carries one oscillatory mode and one
hyperbolic pair (±0.15–0.38, dropped). Both oscillatory modes are negative-energy, ω₀ = 1.05 → 1.33 and ω₁ = 1.50 →
2.69 over E = 0.04 → 1, not degenerate. Mode 0's rungs ω₀T = |n|π, n = 12–19, fall within 8 % in E of the eight
measured block-12 bands with even n at parametric and odd n at Krein bands (predicted E = 0.0580, 0.0759, 0.1020,
0.1416, 0.2052, 0.3159, 0.5328 for n = 18 … 12 and 0.0452 for n = 19; `results/rh27_averaging.txt`).

Why that is not a prediction (the checks, same evening):

1. **The parity is a theorem given existence.** The lemma of `RH7_selection_rule.md` makes any −1 instability
   two-pair; "odd n = Krein" carries no information beyond "there is a band at every odd rung". The one
   non-trivial observation in it is that no Krein collision occurs at +1.
2. **Within 8 % is not evidence.** The half-integer rungs ω₀T = (n + ½)π also fall within 8 % of a band 7 times in
   9, because bands plus an 8 % halo cover 79 % of the log-E range; inside-band hits for the integer rungs are
   4 of 9 against 44 % coverage. Rung 11 lands at E = 1.049 inside the parametric band [1.014, 1.271]: the parity
   fails at the ninth rung. What survives of the positional claim is the slope d(ω₀T)/dE and the count (nine
   rungs, nine bands): one number and an anchor.
3. **The odd rungs are first order in the 0–1 coupling, not "second order".** The width comparison of the first
   reading (|g| ≈ 10⁻⁴–10⁻³ against widths of 0.05–1.2 rung) was against the wrong matrix element; but the
   two-mode combination resonances with the averaged ω₁ do not locate the Krein bands either.
4. **The averaged modes are not the pairs.** The averaged mode 0 tracks one of the two measured pairs' phases
   (correlation 0.62, mean error 0.10 π over 400 energies); mode 1 tracks nothing (−0.36 / 0.08). Gauss residual of
   the averaged eigenvectors at t = 0, E = 0.1 / 0.3 / 0.6: hyperbolic pairs 0.5–1.2 (constraint-contaminated),
   mode 1 (1.7–2.4) clean at 0.01–0.02, mode 0 (1.1–1.25) dirty, 0.04 → 0.23 with E — the mode that tracked a
   pair was the contaminated one. Averaging in fixed coordinates mixes the constraint sector into the physical
   modes because the Gauss direction rotates along the orbit; the fix is to reduce first (`RH_reduced_block.md`).

So: withdrawn as more than slope-plus-count; the normal form of B is right about the mechanism of a sum
resonance but the averaged 8-dim block is the wrong object to feed it.

## D. Measured on the way

5. **Δ(E), the splitting of the two pairs** (`rh5_splitting.py`, from B005): Δ = (φ_A + φ_B) mod 2π is **0.07–0.25 π
   over E = 0.05–0.35** — the pairs are near-degenerate as sets, the decoupled-symmetric structure of RH7 (N-angles
   summing to π) — grows through the wide Krein band to 0.5–0.8 π at E = 0.5–0.72, passes **exactly through 0 at
   E = 0.80** (B005's fake degeneracy: λ and −λ in N, no interaction), and reaches 0.95 π at E = 1, where the pairs
   are antipodal (φ_B ≈ π − φ_A in M: when one is at +1 the other is at −1) and the parametric bands alternate
   between the pairs, which is why rung 11 is parametric. This is hole 5's δ(E) and the quantity RH6 said to
   tabulate; the ladder works below E ≈ 0.35 because the rungs locate the bands only while the pairs are
   degenerate to a fraction of a rung. Figure `results/rh5_splitting.png`.
6. **Signatures at the collisions.** From B005: at every true meeting the colliding members carry opposite Krein
   sign (ψ₁⁺ = −ψ₂⁺); the averaged system calls both modes negative-energy, which is consistent with the half-period
   symmetry reversing one pair's orientation rather than with a ghost–physical collision.
7. **Which dial is gauge-consistent.** Scaling the component cross-coupling by μ (`rh27_mu_scan.py`) is invalid: at
   μ < 1 the block grows by |μ| ~ 10²–10³ per period but the unstable eigenvectors violate the linearised Gauss law
   at the 0.9 level (B002's F118 trap). The gauge-consistent dial is c_j (Gauss violation 10⁻¹² at every value).
8. **The c_j family** (`results/rh27_cj_scan.txt`, `rh5_splitting_cj025.txt`): the physical c_j = 1 is the *least*
   unstable member except the ghost-free c_j = 0 (pure YM, stable). c_j = 0.25: |μ| up to 34, a giant band over
   E ∈ [0.41, 1.0] — a **double +1 fold** (four real eigenvalues off the circle; the first classifier called it Krein
   and was corrected) whose exponent is smooth and single-humped across the band (ln|μ| 2.0 → 3.5 → 2.3, no rung
   structure): a static instability, not a resonance. Δ at c_j = 0.25 is 0.47–0.57 π where the pairs are on the
   circle (the degeneracy is lifted); two narrow Krein bands survive, [0.151, 0.163] (at the same place as at c_j = 1)
   and [0.304, 0.321]. c_j = 0.5: |μ| 6.5 over [0.65, 1.0]. Band positions and even types move with c_j (the
   0.58–0.70 window is Krein at 0.75 and parametric at 1). So Δ controls whether collisions occur at all, and the
   stability ordering of the family is a +1 phenomenon; the reduced block identifies it as transient hyperbolicity
   of the light mode (`RH_reduced_block.md`), not negative mean stiffness.
9. **The physical frequency** (`rh5_branch.py`, revised): an FFT of x(t) over one period with a rectangular window
   leaks — e^{iωt}p(t) is not periodic over [0, T) and the phase jump at the wrap puts a single line into 1/k
   sidebands, which is what a first run's "sidebands to ±8Ω" were. Dividing out the Floquet phase (principal
   branch) and transforming the periodic part (periodicity residual ≤ 2 × 10⁻²): the spectra are **two-line** — a
   line at |f| ≈ 0.85–1.2 (0.93 at E = 0.05, 1.1 at 0.1–0.3, 0.84–0.87 at 0.36–0.55, 1.15–1.2 at E ≥ 0.85; the
   dressed Lee–Wick mass) and a quasi-static line at |f| ≲ 0.1, each with 20–50 % of the power, plus a weaker line
   at ±(1–2)Ω at the few-percent level, near the 10⁻² periodicity floor; nothing near 1.7–2.4. Both pairs carry the
   same lines. The slope-derived "0.8" of the signature-oriented unwrapping (`rh5_unwrap.py`: ω_A ≈ 0.78 / 0.76 /
   0.82 and ω_B ≈ 0.83 / 0.79 / 0.81 in three windows below E ≈ 0.22, splitting above 0.27) is one branch of this
   line, and the ladder's 0.774 another; neither is the frequency. Positive-member identity switches at fold
   re-entries (E = 0.056, 0.109, 0.267 for A; 0.248, 0.730 for B) are within-pair orientation flips.
10. **A swept oscillator, not a beat** (`rh5_instfreq.py`, energy-weighted phase rate of the field components along
    one period; `results/rh5_instfreq.png`): at E = 0.3 one pair sits near |ω| ≈ 1 for most of the period and sweeps
    to 0 twice per period while the other peaks where it dips — the two reflection sectors trading roles every half
    period; at E = 0.1 the same with one pair mostly soft (0.2–0.3) and spiking to ≈ 1 twice; at E = 0.6 one pair
    sweeps 0–1 four times (the other is real there, inside the parametric band). So the two-line spectrum is a
    swept oscillator; the zeroth order is adiabatic/WKB with a stiffness that drops once per half period, the phase
    advance per period is ∫ω_inst dt, and that is where the frequency of ⟨H⟩ goes wrong. **Where the soft point
    is**: this estimator puts the dips at t = T/4 and 3T/4, where x crosses zero and |y| peaks (the bounce); the
    frozen spectrum of the exact reduced block puts the light mode's softening at t = 0 and T/2, the turning points
    (x extremal, y = 0), with both reduced modes stiffest at the bounce, and the per-segment growth of RH3 (growth at
    the bounces, decay at the turning points) agrees with the reduced block. The two measure different things — the
    phase rate of a Floquet solution's field components against the eigenfrequencies of the frozen reduced
    generator — and the discrepancy is not resolved here. The paper text uses "turning points" (watch-list items
    1 and 8).

## E. What this set as the target, and what the reduced block did with it

Target for the reduced block as stated at the end of the evening: two ghost pairs with multipliers degenerate
(Δ ≲ 0.25 π) below E ≈ 0.22, each a line near 1 (0.85–1.2) plus a comparable quasi-static line, nothing at
1.7–2.4, ω_eff²(c_j) with a negative minimum near c_j ≈ 0.3–0.4 and marginal at c_j = 1; checks against ωT mod 2π,
not ω. `RH_reduced_block.md` (22 Oct) delivered: the 4-dim block, exact to 8 × 10⁻¹¹ on the Floquet side; a light
mode (the quasi-static line) and a heavy mode at 1.0–1.7 (the line near 1); the c_j = 0.25 instability as transient
hyperbolicity rather than negative mean stiffness; and the finding that first-order averaging is not Floquet-accurate
even on the exact reduced block (ωT mod 2π off by up to a half-turn). Two structural remarks from the evening carry
over as its open items: if Φ_{T/2} were sector-diagonal, N = S̃Φ_{T/2} would have eigenvalues ±√(ab) — antipodal in N,
degenerate in M, the low-E doublet with the orientation reversal that makes the ±i meetings Krein-unstable — so
**Δ(E) is the off-diagonal part of Φ_{T/2} in the sector basis, the inter-sector mixing within one half period**,
computable from the reduced block without averaging; and the natural model of N is quiet-phase WKB propagation ×
a transfer matrix across the soft point × S̃, with that 4 × 4 symplectic matrix carrying Δ, the Krein widths and, at
low c_j, the static instability.

Cross-reference, off this list by design: a ghost going soft at a bounce is what the spin-2 ghost does on Bianchi IX
(`bianchi9/README.md`). Two backgrounds, one behaviour; the transfer matrix is the object that would say whether that
is the mechanism or a coincidence.
