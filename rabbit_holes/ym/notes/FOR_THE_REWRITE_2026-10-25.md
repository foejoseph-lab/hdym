# For the rewrite instance (Saturday 25 Oct 2026): what the rabbit-hole session changed, in watch-list form

You are helping Joseph rewrite `paper/main.tex` in his own words: grammar, flow, missing-content checks, placing
formulas and tables in TeX. You will see the submitted paper and his new text. You will not see the session that
produced the results below; this file is your only source for them. Every number here comes from a committed log
in `rabbit_holes/ym/results/` or `islands/results/` (named per item) — quote those, not memory. Nothing in the
submitted paper is contradicted by this work; the entries are refinements, additions and three labels to retire.

## A. Claims in the submitted paper to watch for, and what they should become

1. **"the second ghost being the longitudinal polarisation of the massive vector"** (Sec. off-diagonal sector, the
   mode count after eq. L2). Keep the count (one gauge zero mode, one Gauss charge, one YM-type and two ghost
   degrees of freedom per block). Retire the polarisation label: on the DR orbit the two ghost pairs are the same
   sector mode twice, split by a phase δ(E) (`rh5_splitting.txt`; `RH6_polarisation.md`), and the exact gauge-
   reduced block shows the two physical ghost modes to be a light one, soft at the turning points, and a heavy one
   at 1–1.7 (`rh_reduced_block_cj1.txt`). The channel-limit frequencies 1 and √(1 + x²) stand as the channel limit.
   If the new text says "longitudinal", flag it.

2. **"a transiently tachyonic mass is not an instability on a periodic orbit"** (Sec. landing points, DR orbit).
   True as stated for the diagonal ghost at the energies given. Do not let it generalise: in the c_j family of the
   off-diagonal block the +1 instability at c_j = 0.25–0.5 *is* transient hyperbolicity of the light reduced mode
   (frozen spectrum real for part of each period; the averaged block is purely oscillatory while the exact
   multipliers reach 34), `rh_reduced_block_cj025.txt`, `RH_reduced_block.md`. If the new text turns the sentence
   into a principle, flag it.

3. **"The only candidate for a positive-measure set on which the ghost stays bounded is the regular islands … whether
   any survive … is open"** and open item (iii). Still open as a positive-measure statement, but two measured facts
   bound it now and should appear if the passage is rewritten: the E = 0.3 island is not a refuge in the full
   homogeneous theory (islands/README, off-diagonal item, pre-session); and at linearly stable energies between the
   off-diagonal bands a transiently pumped ghost of order one sits for ≥ 150–230 orbits without escaping
   (`rh4_lowE_long.txt`, E = 0.0611–0.0869, 3–13 cells per energy, all censored at tmax = 12000).

4. **Open item (iv), "the rate law of the off-diagonal ghost sector, whose derivation is the natural sequel".** Still
   not derived. The apparatus now exists (A.6 below). Do not let the new text claim a derivation.

## B. Additions the new text may make, at the strength the data supports (and no stronger)

5. **The selection rule as a lemma** (off-diagonal section, two paragraphs). On a background with a half-period
   reflection symmetry S̃K(t + T/2)S̃ = K(t), M = N² with N real symplectic, so a real negative multiplier of M off the
   unit circle has square roots ±ir that force the full quadruplet {ir, −ir, i/r, −i/r} in N: it is never isolated.
   Hence a single ghost pair cannot leave the circle at −1; every −1 instability is a two-pair event (a Krein
   collision of N near ±i between members of opposite signature, measured in `krein_b005.txt`), and every single-pair
   instability is a fold at +1. This is the Swift–Wiesenfeld mechanism in symplectic form; cite them, and
   Frauenfelder–Moreno for the symplectic version and its failure beyond two degrees of freedom (BibTeX in C).
   **Do not** let the text say "collisions only at −1" as a theorem: measured, the collisions are at N-angles
   ±(0.44–0.68)π, and the rule is exact only in the decoupled symmetric limit (`RH7_selection_rule.md`).

6. **The ladder, as a measurement** (off-diagonal section, one paragraph, one figure if wanted; `RH2_ladder.md`,
   `rh2_predict.txt`). The block-12 band centres lie on T_c = 4.057(n + 0.210), n = 7–14, parametric at odd n and
   Krein at even n, rms 0.11 time units; an out-of-sample Krein band predicted at E = 0.0421 was found at
   0.0427–0.0435. Above E ≈ 1 the alternation ends (rung 11 is parametric). **Do not** call 0.774 (= π/4.057) a
   frequency: the ladder is branch-independent (ω → ω + kΩ shifts ωT by 2πk) and cannot fix the physical frequency;
   the Floquet solutions' spectra put the physical line at 0.85–1.2 with a quasi-static companion
   (`rh5_branch.txt`, `rh5_instfreq.png`). **Do not** present "ω₀T = nπ of the averaged block predicts the bands":
   withdrawn (half-rung null: bands plus 8 % halo cover 79 % of log E; `RH27_averaging.md`, Review).

7. **The loiter, as a measurement** (fractal / nonlinear section, one paragraph; `RH_overnight_2026-10-22.md`,
   `rh_overnight_2026-10-22.txt`, `rh4_lowE_long.txt`). After the diagonal ghost reaches |w| = 1 the escape is
   memoryless: where the diagonal ghost is itself Floquet-unstable the loiter has median 0.4–0.6 orbits (pooled
   0.58, scaling with the period: Spearman −0.05 in periods vs −0.46 in time units), stretched tail (p90/median
   4.3–5.7 vs 3.3 exponential); where an off-diagonal band pumps it, the same process with τ = 43–106 orbits
   (E = 0.0581: τ = 106, Kaplan–Meier on the exponential after a ~1000-unit dead time; 0.106–0.112: τ = 43–45, exponential
   to ±0.1). The linear phase is the clock: t(0.1 → 1) matches ln 10/rate where the diagonal ghost is unstable
   (68.7 vs 65.9 at E = 0.33, 10.8 vs 11.2 at 0.85) and is far faster than the ghost's own rate where an off-diagonal
   mode leads (41.6 vs 183 at E = 0.24). **Do not** say the loiter is E-independent in absolute time (the earlier
   checkpoint reading; wrong). **Do not** cite a Krein seed-phase effect: closed negative by the 3072-cell ensemble
   (z = 0.7–1.3 after controlling |c|; `rh_overnight_2026-10-22.txt`).

8. **The gauge-reduced ghost block** (appendix). One Gauss row per block, first class; the gauge direction v_G = Ω⁻¹g
   sweeps the whole 8-space over a period, so constant gauge fixing is inadmissible on the DR orbit and a moving
   gauge f(t) = v_G(t) with the exact connection term is used; the second reduction is by the colour-rotation zero
   mode z = (−y, x, …) and its charge Ωz. The 4-dim reduced monodromy reproduces the four ghost multipliers to
   ≤ 8 × 10⁻¹¹ at every E and c_j ∈ {0.25, 0.5, 1} (`rh_reduced_block_cj*.txt`). Its frozen spectrum: a light mode
   soft at the turning points (0.13–0.7, hyperbolic at E = 0.6 for c_j = 1 and at every E for c_j = 0.25) and a heavy
   mode at 1.0–1.7. First-order averaging of the exact block is not Floquet-accurate (ωT mod 2π off by up to a
   half-turn): state this, do not quote averaged frequencies as predictions. Ω must keep c_j symbolic.

9. **The c_j family** (appendix, with 8): the physical c_j = 1 block is the least unstable member of the family
   except the ghost-free c_j = 0 (c_j = 0.25: |μ| up to 34, a double +1 fold over [0.41, 1] with a smooth
   single-humped exponent, i.e. static; `rh27_cj_scan.txt`). Splitting Δ at c_j = 0.25 is 0.47–0.57π (no near-
   degenerate pairs), two narrow Krein bands survive ([0.151, 0.163] at the same place as at c_j = 1).

10. **δ(E)** (appendix, one figure `rh5_splitting.png`): Δ = (φ_A + φ_B) mod 2π is 0.07–0.25π for E = 0.05–0.35
    (near-degenerate pairs), grows through the wide Krein band to 0.5–0.8π at 0.5–0.72, passes exactly through 0
    at E = 0.80 (B005's fake degeneracy), reaches 0.95π at E = 1 (antipodal pairs). The ladder works where Δ is a
    fraction of a rung.

## C. References to add

```
@article{SwiftWiesenfeld1984, author = {Swift, J. W. and Wiesenfeld, K.}, title = {Suppression of period doubling in symmetric systems},
  journal = {Phys. Rev. Lett.}, volume = {52}, pages = {705--708}, year = {1984}, doi = {10.1103/PhysRevLett.52.705}}
@article{FrauenfelderMoreno2023, author = {Frauenfelder, U. and Moreno, A.}, title = {On doubly symmetric periodic orbits},
  journal = {Celest. Mech. Dyn. Astron.}, volume = {135}, pages = {20}, year = {2023}, doi = {10.1007/s10569-023-10135-6}}
@article{LambRoberts1998, author = {Lamb, J. S. W. and Roberts, J. A. G.}, title = {Time-reversal symmetry in dynamical systems: a survey},
  journal = {Physica D}, volume = {112}, pages = {1--39}, year = {1998}}
```

## D. Things an instance will be tempted to write, and must not

- "the ghost frequency is 0.774" — see 6.
- "odd rungs are Krein because …" as a prediction — the parity is a theorem given existence (5); the content is the
  count and slope.
- "the two ghost polarisations" — see 1.
- "the loiter is independent of energy" — see 7.
- "a Krein instability imprints the seed phase on the escape" — closed negative (7).
- "the averaged block predicts …" — see 8.
- "we derive the off-diagonal rate law" — not done (4).
- Any number not traceable to a named log above.
