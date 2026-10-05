# Rung 1 — off-ansatz linearisation of homogeneous SU(2) Lee–Wick Yang–Mills

*2026-10-01/02. Files: `hdym_rung1.py` (sympy, writes `rung1_system.pkl`), `hdym_rung1_lyap.py`
(CPU Lyapunov spectra, QR), `hdym_rung1_gpu.py` (CuPy/C kernel for the off-diagonal ghost
exponents, large ensembles).*

## 1. Setup

Full homogeneous field A^a_i(t), i, a = 1..3, temporal gauge A_0 = 0, units g = M = 1.
Background: the colour-diagonal solution A = diag(x, y, z)(t) on its ghost-free manifold
(q̈ = −∇V). Perturbation: the six off-diagonal components u^a_i (a ≠ i). The Lagrangian is
expanded to second order in u including the (D^μF_{μ0})² = J_0² piece of the Lee–Wick
term, which vanishes identically on the ansatz and is therefore absent from the paper, but is
quadratic in u off it. The A_0 equation (Gauss law) is expanded to first order.

**Structural results (all sympy checks PASS):**

1. O(u⁰) reproduces the paper's reduced Lagrangian; O(u¹) vanishes (consistent truncation).
2. The six u's split exactly into three 2-component blocks (u₁₂,u₂₁), (u₁₃,u₃₁), (u₂₃,u₃₂),
   related by the permutation symmetry. For block 12:

       L₂ = ½|u̇|² − ½ uᵀm u − ½|ü + m u|² + ½ j²

   with m = [[z², −xy], [−xy, z²]] the (12,21) block of the Hessian of the 9-dof YM potential
   V₉ = ½Σ_{i<j}|A_i × A_j|², and j = x u̇₁₂ − y u̇₂₁ − ẋ u₁₂ + ẏ u₂₁ the linearised colour
   charge (A_i × Ȧ_i)₃. The first three terms are the diagonal-sector structure (normal
   kinetic term, YM potential, Lee–Wick term −½|EOM violation|²); the ½j² term is new.
3. Each block has an exact gauge zero mode u = θ × A (constant θ; verified to solve the linear
   equations on the ghost-free background) and one linearised Gauss constraint G·x = 0, which
   involves u up to u''' and is preserved by the linear flow for *every* state (verified).
4. Counting per block: 2 × 4th order = 8 phase dims, minus gauge, minus constraint = 6 = 3 dof:
   one YM-type and **two** ghost. Total 6 YM + 9 ghost = (9 − 3 gauge) + 3 colours × 3
   polarisations of the massive vector. The second ghost in each block is the longitudinal
   polarisation; it has no diagonal analogue. In a channel (|x| ≫ |y|, |z|) the two ghost
   modes of block 12 sit at frequencies 1 and √(1+x²): the j² term lifts the longitudinal one.
5. The 4th-derivative coefficient is −1 (Lee–Wick sign as in the diagonal sector). The ü
   coefficient matrix is −[[1+x²+2z², −3xy], [−3xy, 1+y²+2z²]].

## 2. Result (b): the diagonal ansatz is transversally *neutral* within homogeneous YM

The SVD A^a_i = (O_colour · diag(x,y,z) · O_space)^a_i shows that the six off-diagonal
directions at a diagonal point are exactly the 3 colour-rotation (gauge) and 3 spatial-rotation
directions. Gauss's law sets the colour-rotation momenta to zero; spatial angular momentum
J_ij = A_i·Ȧ_j − A_j·Ȧ_i is conserved. So the diagonal ansatz is the J = 0 sector of the full
homogeneous dynamics, and off-diagonal perturbations are symmetry directions plus J-changing
directions: polynomial, not exponential, growth. In the linearised YM block (ü = −m u, 4-dim)
there are two conserved linear functionals (j and δJ₁₂) and one zero mode, which forces all
four exponents to zero.

Measured (`hdym_rung1_lyap.py --ym-only`, 18 samples per E): all four exponents
±0.004 at T = 10³, ±0.002 at T = 3×10³ (E = 1), ±0.0013 at T = 5.3×10³ (E = 0.1) — i.e.
ln T / T, the finite-time signature of neutral directions. Zero to within the method.

**Consequences.** (i) The prior expectation "the background leaves the ansatz at
~0.38 E^{1/4}" is wrong; it leaves it at most algebraically. (ii) **Rung 2 (9+9 dof) is moot:**
at J = 0 it *is* the diagonal system; at J ≠ 0 it is the same three-dof potential plus a
centrifugal term. The homogeneous theory has no "generic background" beyond what the paper
already treats. (iii) The scope paragraph of the paper can say this in one sentence, and the
lattice (Rung 3) becomes the only place where the background changes qualitatively.

## 3. Result (a): off-diagonal ghost exponents

Lyapunov spectrum per block along chaotic diagonal backgrounds (QR, `hdym_rung1_lyap.py`,
N = 18, T = 3000 E^{-1/4}):

| E | exponents (block, sorted) | noise floor ln T/T |
|---|---|---|
| 1 | +0.0126(9), +0.0070(5), +0.003, +0.0015, −0.0015, −0.003, −0.0070, −0.0127 | 0.003 |
| 0.1 | +0.0018(2), +0.0013(1), +0.0005, +0.0001, …, −0.0018 | 0.0016 |

Four neutral exponents (gauge, Gauss, YM pair) and two resolved positive ghost exponents at
E = 1; at E = 0.1 the ghost exponents are at the noise floor of the QR method.

Dedicated kernel (`hdym_rung1_gpu.py`), GPU, N = 2×10⁴ (E = 1, 0.3), 10⁵ (0.1), 2×10⁵ (0.03),
`--target-rel 0.02`, 2026-10-02:

| E | T | g₁ | extrap. | g₂ | g₁/g₂ | g₁/λ_diag | local slope of g₁ |
|---|---|---|---|---|---|---|---|
| 1 | 12 000 | 1.180×10⁻² ± 0.1% | 1.178×10⁻² | 5.882×10⁻³ | 2.01 | 0.98 | |
| 0.3 | 12 000 | 3.371×10⁻³ ± 0.2% | 3.398×10⁻³ | 1.732×10⁻³ | 1.95 | 4.2 | 1.04 |
| 0.1 | 60 000 | 6.155×10⁻⁴ ± 0.1% | 6.271×10⁻⁴ | 3.943×10⁻⁴ | 1.56 | 9.1 | 1.55 |
| 0.03 | 120 000 | 2.625×10⁻⁴ ± 0.1% | 2.635×10⁻⁴ | 1.461×10⁻⁴ | 1.80 | **58** | **0.71** |
| 0.01 | 240 000 | 1.326×10⁻⁴ ± 0.2% | 1.316×10⁻⁴ | 7.753×10⁻⁵ | 1.71 | **349** | **0.62** |

(λ_diag here is 0.0120 E^{9/4}; at E = 1 the measured diagonal single-mode rate is 0.004 and the
physical three-vector rate 0.087.)

**Reading.** Neither E^{9/4} nor E^{5/4}. The local slope is non-monotonic — 1.04, 1.55, 0.71, 0.62 — so
two mechanisms cross near E ≈ 0.1: a steep one (channel-tail physics) above, a shallow one
below. The 0.1 → 0.03 slope of 0.71 is within error of **3/4**, the exponent of ε²τ for a
perturbation of strength ε² = E^{1/2} and correlation time τ = E^{-1/4} acting at *zero*
frequency, i.e. without any spectral tail. Zero-frequency response needs either a
velocity-dependent coupling or coupling between modes of opposite energy sign. The j² term
has both features: in the (p, w) variables it is the only term in L₂ that couples the normal
(positive-energy, slow) and ghost (negative-energy, frequency 1) sectors at linear order —
the diagonal sector has no such term, which is why there "the normal sector decouples". With
c_j = 0 the ghost obeys ẅ = −(1 + m)w, the diagonal-type Hill equation with the off-diagonal
Hessian block, for which the paper's machinery predicts the resonant single-mode rate
½·0.0120 E^{9/4} (only J_z² drives w₁₂ through z²) times a mixing factor ≲ 2 from −xy:
roughly 0.006–0.013 E^{9/4}, i.e. 4–9×10⁻⁴ at E = 0.3 and 3–8×10⁻⁵ at E = 0.1.

Predictions for g₁ at E = 0.01 from the 0.03 point were 1.15×10⁻⁴ (slope 3/4), 6.7×10⁻⁵ (5/4),
2.2×10⁻⁵ (9/4); measured 1.33×10⁻⁴, slope 0.62 — shallower than 3/4 and still falling
(1.04 → 1.55 → 0.71 → 0.62). g₂ has slopes 1.02, 1.35, 0.82, 0.58. A slope tending to 1/2 would
mean a rate first order in the O(E^{1/2}) coupling, i.e. a Krein-type (opposite-signature)
instability rather than a second-order noise-driven one; this is to be settled by the c_j = 0
run and the growth-location split, not by fitting four points.

**Growth location.** The kernel now also splits the log-growth of the leading vector by
max|q| ≶ 3E^{1/4} (channel mouth). At E = 1, N = 512: 16% of the growth in the 26% of time
spent in channels, i.e. the off-diagonal ghost grows mostly in the interior, unlike the
diagonal single-mode ghost, whose growth is entirely at depth ~E^{-1/4}. At low E this split
is the cleanest fingerprint of the mechanism: a channel-tail mechanism puts ~100% of the
growth in the channels; a zero-frequency interior mechanism puts it where the time is.

## 3b. c_j = 0 results and the derived E^{5/4} law (2026-10-02 evening)

| E | g₁ | g₂ | mean | predicted 0.0210 E^{5/4} | ratio | growth in channels |
|---|---|---|---|---|---|---|
| 0.3 | 5.189×10⁻³ | 4.945×10⁻³ | 5.07×10⁻³ | 4.66×10⁻³ | 1.09 | 37% of growth in 27% of time (resonant depth D = 2.7 is below the D = 3 cut) |
| 0.1 | 1.541×10⁻³ | 1.403×10⁻³ | 1.47×10⁻³ | 1.18×10⁻³ (1.32 with the measured knee) | 1.25 (1.11) | **75% of growth in 27% of time** |

With c_j = 0, m = z²·1 − xy σ_x has fixed eigenvectors (1, ±1): the pair decouples into
ẅ_± = −(1 + z² ∓ xy)w_±. The xy entry chirps at |x| (not 2|x|) with amplitude √(2J|x|); both
channels resonate (no longitudinal lift). Passage integral with one power of J:
u²ŝ_xy = 16√2πν₀/15 = 0.672 (measured 0.68–0.72 over u = 4–20, `hdym_spectrum_xy.py`), and
λ_± = (4√2π³/15) E^{5/4}/Z = 0.0210 E^{5/4} (`hdym_analytic.py` §6). Confirmed by flux, tail and
exponents. In the real theory (c_j = 1) the j² term lifts the longitudinal mode to √(1+x²) in
channels and detunes this resonance (g₁ = 6.2×10⁻⁴ vs 15.4×10⁻⁴ at E = 0.1), but leaves the
interior pair degenerate and coupled to the YM mode; the shallow low-E law is that interior
contribution. Its derivation is open.

## 3c. The theory at E = 0.1 with the split, and the Q statistic

E = 0.1, c_j = 1 repeat: g₁ = 6.155×10⁻⁴ (identical, same seed), **79% of the growth in the 27% of
time in channels** — same as c_j = 0. So the j² suppression (×2.5) acts within the channels; the
"interior mechanism" reading was wrong at this energy.

Q = σ²T/λ (variance of per-trajectory finite-time exponents × T / mean; the log-gain per growth
event; T-independent), from every `lam_i` file (`make_fig_Q.py`, `fig_Q.pdf`):

| ghost | Q |
|---|---|
| diagonal (3-d) | 1.0–2.5 over E = 0.015–1 (6.4 at E = 0.008) |
| physical 3-vector | 3.5 → 1.0 over E = 0.008 → 0.03 |
| two-field | 1.2–1.4 |
| off-diagonal, c_j = 0 | 3.6 (E = 0.3), 2.7 (0.1) |
| off-diagonal, the theory | 2.1, 2.1, 6.1, 11.9, 18.9 at E = 1, 0.3, 0.1, 0.03, 0.01 (≈ 1.9 E^{-1/2}) |

Parametric mechanisms: Q ~ 1–3, flat. The theory's off-diagonal ghost: growth in rare deep
visits with log-gain ~19 at E = 0.01, <0.2% negative trajectories — an instability lasting the
visit, not a kick per passage. Candidate: the frozen coefficients of the Gauss-charge term change
sign in channels (the frozen scan at E = 0.1 found a real eigenvalue 0.54 at x = 3.8, z = 0.12;
its time average over the fast z oscillation vanishes by the virial relation, so the effect is
not a static tachyon; the derivation is open). Paper text: block D of `additions.tex` is complete.

## 4. GPU runs

Done (c_j = 1): E = 1, 0.3, 0.1, 0.03; E = 0.01 running. Next, in order:

```
%run hdym_rung1_gpu.py --E 0.3 0.1 --cj 0 --N 100000 --t-block 5000 --target-rel 0.02 --target-abs 1e-6
%run hdym_rung1_gpu.py --E 0.1 --N 100000 --t-block 10000 --target-rel 0.02 --target-abs 1e-6
```

The first (~1.5 h) drops the (D^μF_{μ0})² term: if g₁ falls to the 10⁻⁴/10⁻⁵ level predicted
above, the Gauss-charge coupling between the normal and ghost sectors is the mechanism. The
second repeats the c_j = 1 point at E = 0.1 with the new kernel to get the channel/interior
split of the growth (the earlier runs predate it). Output rows append to
`results_rung1/rung1_offdiag.csv` with a `cj` column and the two split fractions.

## 5. For the paper

- Scope paragraph: replace the implied "the ansatz is a measure-zero slice of chaotic
  homogeneous dynamics" by: the diagonal ansatz is the J = 0, gauge-fixed sector of homogeneous
  SU(2) Yang–Mills mechanics (SVD); off-diagonal perturbations are neutral; the full homogeneous
  ghost sector has 9 dof, 3 of which are treated in the paper, and the other 6 (including the
  longitudinal polarisations, which couple through the (D^μF_{μ0})² term) have exponents
  [to be filled from §4] — then one sentence on which of the two scalings holds.
- If E^{5/4}: the abstract's λ_⊥ ∝ E^{9/4} must be qualified as the diagonal-sector law, and
  the physically relevant transverse exponent of the homogeneous ghost-free manifold is the
  off-diagonal one. That would be a bigger paper, not a smaller one.
