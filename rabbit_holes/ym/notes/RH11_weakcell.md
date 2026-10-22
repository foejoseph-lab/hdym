# Rabbit hole 11 — the frozen-channel model's weak Krein cell: an artefact, and which one (22 Oct 2026)

Script `code/rh11_weakcell.py`, log `results/rh11_weakcell.txt`. B002's one in-model Krein cell (D = 3.8946,
J = 0.0315, ψ = 0, c_j = 1; |μ| = 1.002 at arg 0.065π, blocks 12 and 13) was downgraded to "probably an artefact"
without a decisive test. Two tests, seconds each.

**A. The cell is the split gauge Jordan pair, not a ghost collision.** B002's verdict used all eight multipliers
(correctly, since on the non-exact background the drop-4 rule has no exact meaning). Reproducing the cell and
applying the drop-4 labelling: the four ghost-sector multipliers are on the circle to 10⁻¹⁵ (−0.056 ± 0.999i,
0.979 ± 0.204i), and the thing at |μ| = 1.002, arg ±0.065π is the gauge pair, whose Jordan block the frozen
background splits by |μ − 1| = 0.20 ≈ √ε (ε = 0.021) — exactly the suspicion of the 19 Oct response. Scanning
J = J₀/8 … 4J₀ at fixed D: the ghost pairs stay on the circle at every J (|μ| − 1 at 10⁻¹⁵), the split grows
0.07 → 0.21 with √ε as it should, and no threshold in J appears (a real collision would switch on when the pairs
meet).

**B. The exact background agrees.** Integrating the real Yang–Mills equations from the same initial data for one
transverse period (depth drifts by 2.1 %, matching the model's ε = 2.1 %): the linearised Gauss law is then
conserved to 6 × 10⁻¹⁵ (the frozen background leaks 0.39), the only Gauss-clean eigen-pair of the transfer matrix is
on the circle (−0.053 ± 0.999i, |μ| = 1.0002, residual 0.012), and the off-circle eigenvalues (1.035 ± 0.21i) have
Gauss residual 1.15 — constraint-sector content of a non-periodic transfer matrix, not growth. The singular values
of the Gauss-clean restriction over the period agree between exact and frozen backgrounds to 2–7 % (3.39 vs 3.32,
1.66 vs 1.79, …), so the frozen model is adequate for the physical sector at this cell. (Block 23's 1:1 resonance
cluster, singular values ~190, is B002's F118 and INFO only.)

**Verdict.** Artefact, mechanism identified: the broken constraint's Jordan splitting, |μ − 1| ≈ √ε, mislabelled
as a weak Krein quadruplet by the all-8 rule. B002's "the c_j = 1 off-diagonal ghost is Floquet-stable in every
block at both phases within the frozen-channel model" holds without the exception. Hole closed.
