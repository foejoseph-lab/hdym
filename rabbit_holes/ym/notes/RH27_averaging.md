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
