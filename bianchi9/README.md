# bianchi9/ — Weyl ghost in Bianchi IX minisuperspace of quadratic gravity

Working notes 5–6 Oct 2026, closure 10 Oct; rewritten 23 Oct 2026 from the state at that date. Status: closed 10 Oct for the one-ghost theory (single wall, one trajectory family, θ = π, 2.9, 2.6, 2.3, 2.0); two-ghost case unfinished; the README-vs-code mismatches of 12 Oct verified against the scripts on 23 Oct (all stand; the Dirac one resolved in substance) and held.

## Question and answer

The Lee–Wick Yang–Mills question (`paper/main.tex`, `ym/README.md`), asked of the Weyl ghost. In the diagonal Bianchi IX minisuperspace of

    S = ∫ d⁴x √−g [ R/2 − C_{abcd}C^{abcd}/(4M²) ]        (M̄_P = 1, no R² term)

the vacuum-GR mixmaster solutions form an invariant submanifold (Einstein branch); the theory carries a spin-2 ghost of mass M. Is the Einstein branch a repeller for the ghost in the chaotic collapse, at what rate per bounce vs curvature/M²?

**Answer (closure 10 Oct; four methods agree: transfer matrix, Bogoliubov quadrature, Kasner exponents, nonlinear seeds).**

1. H < 0.2M — adiabatic. A bounce excites the ghost by ∝ e^{−(π/2)M/H_b}, no more than a healthy scalar of the same mass (computed as a control): the coupling is even under the bounce (β̇², α̇) and vanishes at the turning point — opposite parity to the Yang–Mills coupling w·∇V(q), linear and largest at channel entry, which kicks at every energy.
2. H ≈ 0.2M — the ghost wakes in a Kasner epoch, not at a bounce: amplitude law from the adiabatic t^{−1/2} to log-periodic.
3. H > 0.2M — the Einstein branch is a repeller by a **power law**, δβ ∝ 1/t (deviation ∝ (H/M)² in ghost units; physical GR mode t^{−0.3}, relative t^{−0.7}), the same for all bounces; a bounce at H_b ≳ M only sets the prefactor, O(10). No exponential runaway anywhere.
4. H ≫ M — Querella (thesis, gr-qc/9902044), Cotsakis–Demaret–De Rop–Querella 1993, Barrow–Cotsakis 1989: mixmaster chaos is structurally unstable in fourth-order gravity, the generic approach a monotonic, isotropising power law. No billiard at H ≫ M; the departure found here is the transition to the quadratic theory's own asymptotics. The Weyl ghost does not ride the mixmaster; it ends it.

## Layout and reproduction

Scripts and generated inputs (`bianchi9_system.pkl`, `bianchi9_rhs.py`, `bianchi9_symp_gen.py`, `kasner_linear.pkl`) here; committed outputs in `results/`; scripts write to the working directory; run commands in each docstring; reproduces on two machines (no committed record). DOP853 crawls for H/M ≳ 0.2 — use Radau (3–15 s per seeded run vs an hour).

| Script (section) | Outputs in `results/` |
|---|---|
| `bianchi9_derive.py` (Derivation; eight internal checks pass, 35.6 s in `derive.log` (stated as 37 s until 28 Oct 2026); check 7 tests the random on-shell points found (4 of 8 in the log) against |dCon/dt| < 1e-6 × an off-shell scale, `bianchi9_derive.py`:272-283) | `derive.log`; writes `bianchi9_system.pkl` |
| `bianchi9_codegen.py` → `bianchi9_rhs.py` (85 kB; 0.4 ms per rhs call (no committed log), matches lambdify to 1e-14; do not lambdify the closed-form `b4` — ten times slower (no committed log)) | `codegen.log` |
| `bianchi9_test.py` (1) | `sweep_M{5,10,20,50,100,200}.log/.npz`, `sweep_first_look.png` |
| `bianchi9_symp.py` → `bianchi9_symp_gen.py` (Derivation) | — |
| `bianchi9_action.py` [`--wkb`] (2) | `action_M{10,20,40,80,160}.log/.npz`, `action_M{20,28,40,56,80,160}_wkb.log/.npz`, `control_M80_wkb.log/.npz`, `action_sweep.png`, `action_sweep_wkb.png` |
| `bianchi9_bogoliubov.py` (3) | `bog.log`, `bog.npz`, `bog_diag.png` |
| `bianchi9_fields.py` (4) | `fields.log`, `fields.npz`, `fields_light.npz`, `fields_n.png` |
| `bianchi9_kasner.py`, `bianchi9_kasner_check.py`, `bianchi9_kasner_transfer.py`, `bianchi9_onset_frozen.py` (5) | `kasner_onset.npz/.png`, `kasner_pi.log`, `onset_frozen.npz`; writes `kasner_linear.pkl` |
| `bianchi9_seed.py` [`--theta`, `--seed-pol p|m`] (6) | `seed.log`, `seed.npz`, `seed_growth.png`, `seed_th{3.1416,2.9,2.6,2.3,2.0}_{p,m}.npz` |
| `bianchi9_six.py` (7) | `six.log` |

## Derivation (`bianchi9_derive.py`, 5 Oct)

Cartan calculus on ds² = −N²dt² + a²σ₁² + b²σ₂² + c²σ₃², Misner variables (α, β₊, β₋); N = 1 equations solved for α″, β⁗; constraint; state y = (α, α′, β₊, β₋, β₊′, β₋′, β₊″, β₋″, β₊‴, β₋‴).

1. Weyl² vanishes for isotropic a = b = c (any a(t)); Gauss–Bonnet is a total derivative on the ansatz: C² is the whole non-GR content.
2. No α″, no N″ in the C² Lagrangian (conformal invariance); Hessian in (β₊″, β₋″) = −6e^{3α}/M² · 𝟙 (polarisations degenerate, isotropic). With the GR term +3e^{3α}|β′|² the transverse block is Pais–Uhlenbeck, ½ẋ² − ẍ²/(2M²), frequencies 0 and exactly M: non-tachyonic, pole at M with −C²/(4M²).
3. α equation second order, α″ coefficient 6(M² + β₊′² + β₋′²)e^{3α}/M² ≠ 0; β⁗ coefficient matrix (M² + |β′|²)𝟙 − β′β′ᵀ times −1296(M² + |β′|²)²e^{12α}/M⁸ (`derive.log`, "beta-4 coefficient matrix", read off its (0,1) entry; "up to positive factors" until 28 Oct 2026), determinant ∝ M²(M² + |β′|²) ≠ 0: regular on the whole phase space.
4. Phase space 10-d, one constraint → 9-d; Einstein manifold 6-d → 5-d; transverse 4-d, two oscillators.
5. Einstein manifold invariant (Bach tensor vanishes on Einstein spaces), checked on random constraint-satisfying points to 1e-9.

`bianchi9_symp.py`: Ostrogradsky map z(y) (α, β, β′ | p_α, P₁, P₂; α″ in P₁ on shell), pulled-back Ω(y); the linearised Misner-time flow preserves Ω to 1.4e-9 for tangent vectors with dH·v = 0, off the constraint surface Ω drifts by O(1).

## Results

One GR bounce off the β₊ wall, head-on (θ = π) unless stated, α_b = −0.38 (`sweep_M*.npz`; −0.36 in `action_*.log`, `bog.log`, `fields.log`, `seed.log`), H_b ≈ 2.0 in units with initial |β′| = 1; M = 5…200, i.e. H_b/M = 0.4…0.01. Einstein-manifold invariance 1e-11 (relative), constraint 1e-10, every run.

### 1. Bounce, first look (`bianchi9_test.py`, 5 Oct)

DOP853 in Misner time dt = e^{3α}dτ, transverse projection δv = v_hi − (∂gr/∂y_lo)v_lo, observable ½ log J, J = e^{3α}(δβ″² + δβ‴²/M²)/M³ (crude, phase wobble ±0.2); on-manifold tangent tracked, unused.

| H_b/M | 0.40 | 0.20 | 0.10 | 0.04 | 0.02 | 0.01 |
|---|---|---|---|---|---|---|
| transverse log-growth over Δα = −0.85, α −0.05 → −0.90 (`sweep_M*.log` line 8, "transverse": log √(ghost action), `bianchi9_test.py:102`) | +3.14 | +2.76 | +0.64 | +0.20 | −0.17 | +0.05 |
| "log action gain across the bounce", as stated until 28 Oct 2026 (source unknown; the logged quantity is log √action, so "log action" is off by a factor 2) | +3.0 | +2.6 | +0.4 | +0.1 | −0.1 | 0.0 |

Constant within the wobble for H_b/M ≤ 0.1; the growth at 0.2, 0.4 is not a step but sets in after the bounce as H/M → 1 along the following Kasner epoch (H ∝ e^{−3α}).

### 2. Bounce, transfer matrix (`bianchi9_action.py`, 5–6 Oct)

10×10 tangent matrix in scaled coordinates (δβ″ ~ Mδβ, δβ‴ ~ M²δβ; unscaled, atol makes the integrator crawl), frozen ghost modes symplectically normalised (iEᴴΩE = ±1), transfer matrix T (real 4×4), singular values e^{±r} → squeeze r (action gain cosh 2r; mean log s = 0). Frozen frequencies 0.97 M; frozen Re λ ≈ 3H/2 is kinematic (a^{−3/2}); |dH·e|/|dH| ∝ (H/M)² (2e-1 … 1e-3 for M = 10 … 160), a frozen-mode error.

| H_b/M | 0.012 | 0.025 | 0.049 | 0.099 | 0.198 |
|---|---|---|---|---|---|
| r at α = −0.5 (just past the bounce), frozen frame | 1.6e-4 | 1.0e-3 | 8.8e-3 | 0.13 | 1.2 |

The three smallest are upper bounds (frame adiabatic only to O((H/M)²)): r ≤ 1e-3 for H_b/M ≤ 0.025, r ≤ 1e-2 at 0.05, r = 0.13 at 0.1.

`--wkb` (6 Oct): first adiabatic correction ẽ = e + Σ_k v_k (w_k·ė)/(λ_k − λ_pair), ė from the projector onto the near-degenerate pair (central difference); |dH·e|/|dH| 8e-3 → 1e-3 at M = 40, 3e-3 → 1e-4 at M = 80. Artefact floor from a control run away from the wall (`control_M80_wkb.*`, θ = 0): r_floor ≈ 1.7e-5 (H/M / 0.0127)^1.8 at matched H/M.

| H_b/M | 0.099 | 0.071 | 0.049 | 0.035 | 0.025 |
|---|---|---|---|---|---|
| r, corrected, at α = −0.5 | ~0.3 | 0.040 | 2.9e-3 | 5.6e-4 | 1.4e-4 |
| H/M at α = −0.5 (`action_M*_wkb.npz`) | 0.119 | 0.085 | 0.060 | 0.043 | 0.030 |
| artefact floor, the fit above at that H/M (the row as stated until 28 Oct 2026) | 1e-3 | 5e-4 | 3e-4 | 1.5e-4 | 8e-5 |
| artefact floor, the control run itself at that H/M (log-interpolated in `control_M80_wkb.npz`, smax) | — (0.119 is beyond the control's 0.107) | 9.4e-3 | 1.35e-3 | 4.0e-4 | 1.2e-4 |
| r / control floor | — | 4.2 | 2.1 | 1.4 | 1.1 |

0.099 only marginally adiabatic (mean log s = +0.13 ≠ 0); below 0.035 the floor dominates. (28 Oct 2026: the floor row
was the power-law fit evaluated at each run's H/M at α = −0.5; the fit follows the control only for H/M ≲ 0.03 and
falls below it by ×1.5–18 above that, so against the control run itself r is within ×2.1 of the floor already at
H_b/M = 0.049 and ×4.2 at 0.071. This strengthens the withdrawal below; whether to restate the reading is the
author's call.) Withdrawn (6 Oct, late): the reading of 0.071 → 0.049 → 0.035 (×14 then ×5 per factor 1.4; local exponent 7 → 5; as e^{−cM/H_b}, c ≈ 0.44 → 0.22) as bounce kicks — the 0.049 and 0.035 values are second-order frame terms: the two first-order-removed estimators (`--wkb`; Bogoliubov with first endpoint term removed) disagree by ×2.5–5 there. A fixed spectral-tail power law (the YM mechanism) stays excluded.

### 3. Bounce, Bogoliubov quadrature and Taub estimate (`bianchi9_bogoliubov.py`, 6 Oct)

B = ∫ C e^{−iφ} dτ, C = iE₋ᴴ Ω Ė₊, in a parallel-transported, symplectically normalised frame (without parallel transport the within-pair unitary is arbitrary and the phase-sensitive result meaningless). |C(τ)| has a **dip** at the bounce, not a peak (its rate 2β̇₊β̈₊ vanishes at the turning point; same feature as the flattening of the M = 10, 20 curves at α ≈ −0.4 in `action_sweep.png`); cumulative |B| has **no step at the bounce**, only the oscillating endpoint term of a smooth integrand growing with H/M. Taub estimate: the bounce is β̈ ∝ sech²(4vt), Fourier transform ∝ ω/sinh(πω/8v); phase rate at the bounce 2ω = 19.5 (M = 28), 27.3, 37.9, 53.8 (M = 80) (`bog.log`; stated until 28 Oct 2026 as "phase per
bounce 2ω_τ/a ≈ 10 … 27", half these values, the factor ½ not recorded), so the kick is O(e^{−(π/2)M/H_b}) — 10⁻¹⁰ at H_b/M = 0.07, 4×10⁻⁴ at 0.2, negligible for H_b < 0.15M, below anything measured (c ≈ π/2, O(1) corrections from Δα ≈ ½). Measurable frame terms are reversible dressing; the planned three-decade r(H_b/M) sweep was cancelled.

### 4. Bounce, healthy probes as control (`bianchi9_fields.py`, 6 Oct)

Massive scalar (sees only the volume at k = 0) and massive Dirac field (sees the anisotropy through the spin connection, the square roots of the walls), homogeneous, no back-reaction. Dirac axial coupling A = ¼(a/bc + b/ca + c/ab): the committed `fields.log` has H_D eigenvalues ±1.33972 at τ = 0 (m = 1), which is √(1 + A²) with A = ¼Σ exactly (0.89156) and not with ⅛Σ (the log's own check line, written by an earlier printout, says "A = (1/8) sum a/bc = 1.09486" and so appears to fail; the script now prints √(1 + 4·Aax²) with `Aax` = Σ/8, the half-coupling it stores, and the Dirac integration itself uses the spin connection, not `Aax`). The ⅛ form that appeared in the parallel-branch section was that variable. Checks: Hermiticity of H_D for χ = e^{3α/2}ψ (4e-16; the wrong Clifford convention gives 2.0 and n ≡ ½), E = √(m² + A²). Observable: adiabatic particle number n at α = −0.5.

| m/H_b | 0.099 | 0.071 | 0.050 | 0.035 | 0.025 |
|---|---|---|---|---|---|
| scalar n | 2.0e-5 | 2.5e-6 | 2.4e-7 | 3.3e-8 | 3.6e-9 |
| Dirac n | 5.1e-5 | 8.9e-6 | 1.9e-6 | 1.6e-7 | 1.1e-7 |

No step at the bounce for either (scalar: smooth dressing with the same turning-point dip; fermion: a beat pattern passing through unchanged). The scalar n falls ×8.3, 10.1, 7.3, 9.3 per step of the table, the Dirac n ×5.8, 4.6, 12.5, 1.5 (`fields.log`, n(alpha=−0.5) columns, successive ratios, e.g. 2.029e-5/2.450e-6; stated as "both fall ~×8–10 per factor 1.4" until 28 Oct 2026); n ∝ (H/m)⁶ (scalar, second adiabatic order) and (H/m)⁴ (fermion, first order, A ~ H): reversible dressing. The ghost's n_ghost = sinh²r, r ∝ e^{−(π/2)M/H_b}, is at or below these: the quiet bounce is the background's smoothness, not the ghost's sign.

Sanity check (`fields_light.npz`), m ≪ H_b: boson production saturates (n ~ 0.7, no adiabatic vacuum for m < 3H/2 — the sqrt warning); fermion production peaks at m ≈ H_b, falling both sides (conformal protection at m → 0, adiabatic at m → ∞): Zel'dovich–Starobinsky 1971 reproduced. Cosmetic: Dirac large-m curves in `fields_n.png` under-sampled (400 points); scalar lower-left fuzz is the rtol = 1e-10 floor.

### 5. Kasner epoch, onset and exponents (`bianchi9_kasner.py`, `bianchi9_kasner_check.py`, 6 Oct)

Bianchi I limit of `L` (e^{3α}-homogeneous part; GR part checked = 3e^{3α}(−α′²+β′²)); Kasner α = ⅓ln t, β = ⅓ln t (cos θ, sin θ) verified to solve the full equations (Ricci-flat ⇒ Bach-flat); by scale invariance the linearised 10-d problem depends only on x = Mt (H/M = 1/3x) and θ; integrated from x = 200 to 0.01 (committed script: 60 → 0.02, Open).

* Symbolic t → 0: pure-Weyl indicial determinant s³(s−2)³·Q(s,θ) with Q ≡ 0 (conformal invariance); exponents come from the next order, numerically.
* (No committed log for this item and the next: `bianchi9_kasner.py` and `bianchi9_kasner_check.py` print them.) Frozen spectrum in Euler form xλ: ghost pair −0.5 ± ix for H ≪ M (adiabatic t^{−1/2}); Re crosses 0 at **H/M ≈ 0.2**; for H ≫ M the pair becomes 0.27 ± 0.56i (log-periodic, t^{0.27} cos(0.56 ln t)) and 0.45 ± 0.35i. The other frozen exponents (−4.8, −2.5) are not realised as growth (frozen eigenvalues mislead in a non-autonomous Euler system).
* Integrated exponents, H ≫ M (fit over x ∈ [0.02, 0.1], all θ): generic data δβ₊ ∝ t^{−1.0}, δα ∝ t^{−1.0}, δβ₋ ∝ t^{−0.2…−0.3}. GR time-shift mode t^{−0.99999999999} (exact check). GR Kasner-angle mode (the physical GR perturbation): δβ ∝ t^{−0.30}, δα ∝ t^{−0.81}.
* Scaled ghost amplitude |δβ|√x flat for 0.01 < H/M < 0.1 within the preparation transient (the `bianchi9_kasner_check.py` envelope, no committed log; the committed `kasner_onset.npz` gmax at θ = π rises ×5 from H/M = 0.01 to 0.1, 298 → 1494), knee at 0.2–0.3 for all Kasner angles run (recorded as "all four"; the seeded map below uses five).

Local Frobenius exponents σ (perturbation ∼ t^σ; toward t → 0 the smallest Re σ dominates), scaled Jacobian in s = ln t, θ = π:

| Mt | 10⁴ | 10 | 1 | 0.3 | 0.1 | 10⁻² |
|---|---|---|---|---|---|---|
| ghost pair(s) | 2 ± iMt (×2) | 2.006 ± 9.93i, 2 ± 10i | 2.20 ± 0.77i, 2 ± i | 2 ± 0.3i, 1.85 ± 0.16i | 2 ± 0.1i, 1.99, 1.12 | 2 ± 0.01i, 2.00, 1.00 |
| on-manifold | 1, 0,0,0,0, −1 | 0.99, 0⁴, −1 | 0.41, 0⁴, −1 | 0.20, 0⁴, −1 | 0.05, 0⁴, −1 | 0.00, 0⁴, −1 |

For Mt ≫ 1 the WKB correction (−½ d ln ω/d ln t) brings the ghost's Re σ from 2 to 3/2 (adiabatic invariant t M A² = const); −1 is the time-translation gauge mode; physical Kasner perturbations (angle, p) sit at σ = 0. At Mt ≲ 1 the ghost turns from oscillatory to power-law; the Radau run crawls below Mt ≈ 1, its frozen-frame observable meaningless there (mode crossing). Withdrawn (10 Oct): the reading "no transverse exponent below 1, so Kasner is relatively attracting and the Einstein branch not a classical repeller" — replaced by the integrated exponents above (relative growth is also direction-of-time dependent); the bounce half stands. Literature: Toporensky & Müller, arXiv:1603.02851 (R² + Ric², χ = α/3β): Kasner locally stable toward the singularity for χ > 4; this theory is χ → ∞.

Independent observable (`bianchi9_kasner_transfer.py`, `kasner_pi.log`, "log s: max" column): r ≤ 1e-3 at H/M = 0.054 and < 1e-4 for H/M ≤ 0.021 (stated as "r < 1e-4 for H/M ≤ 0.054" until 28 Oct 2026), r = 0.12 at H/M = 0.136, r = 1.03 at 0.34 — same onset ≈ 0.2; frozen frequency 0.9999 M at Mt = 100 with |Re λ| = 3H/2 to 4 digits (log: Re = −0.0050, 3H/2 = 0.0050). `bianchi9_onset_frozen.py` (`onset_frozen.npz`): Re λ_ghost, |Im λ_ghost|, max Re λ vs H/M along the real collapse, M = 2, 5, 10, 20; data only. Exponential growth could only come from bounces at H ≳ M (κ ~ H ~ M, kick O(1)) with chaos re-entering there — excluded by the literature (Answer 4).

### 6. Nonlinear seeded runs through bounces at H_b ~ M (`bianchi9_seed.py`, 6 Oct)

GR data + ghost seed δβ₊″ = 10⁻⁶ M H₀, full nonlinear equations, Radau, M = 1, 2, 3, 5 so the bounce sits at H_b/M = 2, 1, 0.66, 0.4; run to H/M = 8. Deviation from the Einstein manifold in ghost units (δβ″/MH, δβ‴/M²H) relative to the seed; unseeded reference stays on the manifold to 3e-7; back-reaction linear in the seed (no committed log: `seed.log` has one ε). All four curves fall on **one universal track, deviation ∝ (H/M)^{2.0}** (fits 2.03, 2.04, 2.06, 2.34; no committed log) — exactly δβ ∝ 1/t, the Kasner exponent reproduced through a real bounce. The bounce only lifts the trajectory onto the track: ×2.5 at H_b/M = 0.4, ×6 at 0.66, ×13 at 1, ×50 at 2; after it every trajectory forgets it. Ignition at H/M ≈ 0.2–0.3 on every curve (no committed log).

Prefactor map (`seed_th*_*.npz`), M = 2 and 5, θ = π, 2.9, 2.6, 2.3, 2.0 (head-on to glancing; glancing bounces occur later, at larger H_b/M). Kick at matched H_b/M ≈ 1–1.3 (growth from α = −0.2 to −0.5): ×13 (θ = π), ×12 (2.9), ×9 (2.6) — **O(10) at H_b ≈ M**. Polarisation: head-on, β₋ (odd under β₋ → −β₋) is barely excited (×9 vs ×800 for β₊ at M = 2; ×5 vs ×420 at M = 5) and drifts with t^{−0.3}, not t^{−1}. The (H/M)² track is that of the polarisation **aligned with the shear**; oblique bounces mix the two (θ = 2.0: the β₋ seed ends above the β₊ one), so in a real mixmaster sequence the aligned law governs after the first oblique bounce.

### 7. Two-ghost (sixth-order) theory (`bianchi9_six.py`, `six.log`, 6 Oct) — unfinished

S = ∫√−g [R/2 − (k₂/4)C² − (k₄/4)C□C], k₂ = 1/M₁² + 1/M₂², k₄ = 1/M₁²M₂² (spin-2 poles 0, M₁, M₂, Krein signs +, −, +), Bianchi I. Checked: Weyl traceless; C² = C□C = 0 isotropic; √−g(C□C + |∇C|²) a total derivative on the ansatz (EL vanish, validating the frame covariant derivative); Einstein manifold invariant; phase space 16-d; Noether identity α′E_α + β′·E_β = d(Con)/dt (in the derivation; not in the committed script, Open). The α-equation is **fourth** order with coefficient −6|β′|²e^{3α}/M₁²M₂², vanishing at zero shear (the one-ghost coefficient never did): a Kasner point with β′ = 0 is singular for the solve. The β⁽⁶⁾ coefficient matrix is exactly the transverse projector |β′|²𝟙 − β′β′ᵀ (longitudinal shear mode of lower order): a constrained higher-order system. The script runs to the singular-matrix point; everything before it is reusable. Recipe on the β₋ ≡ 0 branch: α⁗ from E_α, β₊⁽⁵⁾ from Con, β₋⁽⁶⁾ from E_{β₋}; 15-d closed system; then the exponent table vs M₁t for M₂/M₁ = 3 and whether the two ghosts Krein-collide along Kasner.

## Open

* Two-ghost theory: constrained reduction (recipe above), exponent table, Krein collision yes/no.
* Probe suite (CPU; doubles as the Python 3.13 acceptance suite): control; ghost + light; control + light; control + ghost; all three — runs 4–5 (control + ghost on one trajectory, m = M resonance) not done. "Light" is a probe, not a background wall: a homogeneous magnetic field on diagonal Bianchi IX did not pass d⋆F = 0 in a quick check; arXiv:1910.11970 adds 2b²e^{2α−4β₊} to the potential without discussing the constraint — to be understood first. Head-on, wall-centre geometry (worst case by the Taub argument); corner-era run and Kasner-angle eigenvalue loop separately.
* Direct tangent-matrix integration through Mt ~ 1 (Radau, scaled variables) with the transverse / on-manifold split against physical directions (dH·v = 0, modulo the time-translation mode); Kasner-angle sweep.
* Other wall types / Kasner angles for the prefactor (cheap with Radau); n_ghost vs n_scalar on one axis.
* (Optional) superadiabatic series to optimal order to exhibit the e^{−(π/2)M/H_b} kick explicitly.
* Port the rhs to C for ensembles over bounce parameters (u, wall type).
* Comparison of the Cotsakis–Demaret–De Rop–Querella / Barrow–Cotsakis minisuperspace Lagrangian with `L` in the pkl (literature read 10 Oct; comparison not recorded).
* Independent review of the conclusion by a reader who has not seen the work.

### Known README-vs-code mismatches (found 12 Oct, held; verified 23 Oct)

Each checked against the committed scripts on 23 Oct; all nine stand as written. Not fixed (held by the decision of 12 Oct). (28 Oct 2026: the Dirac check line's label is fixed in `bianchi9_fields.py`, which now prints "E = sqrt(m^2 + A^2)";
the committed `fields.log` keeps the old line.)

* six figures in results/ with no writer — `savefig` occurs only in `bianchi9_kasner.py` (`kasner_onset.png`); `action_sweep.png`, `action_sweep_wkb.png`, `bog_diag.png`, `fields_n.png`, `seed_growth.png`, `sweep_first_look.png` were made outside the committed scripts;
* `bianchi9_bogoliubov.py` default `--out bogoliubov.npz` vs the committed `bog.npz` (renamed on save), and no endpoint term removed — the script is the plain first-order quadrature; the "first endpoint term removed" estimator mentioned under Results 2 is not in the repository;
* `bianchi9_symp.py` computes the Ω-conservation residual JᵀΩ + ΩJ + dΩ/dt for the full 10-d Jacobian at a point, with no restriction to dH·v = 0, so the Derivation's "1.4e-9 for tangent vectors in the constraint surface" is not what the committed script tests;
* Dirac axial coefficient: ¼ is the coupling and Σ/8 the stored half-variable `Aax`; the committed `fields.log` check line is from an earlier printout and reads as a failure — resolved above (Results 4), log not regenerated;
* `bianchi9_kasner.py` integrates x = 60 → 0.02 (`run()` defaults, no command-line override) vs the 200 → 0.01 quoted under Results 5;
* `bianchi9_kasner_check.py` fixes θ = π; `bianchi9_kasner.py` itself loops over θ = π, π/2, 0, π/3 (the "all four" of Results 5), so "all θ" means those four and the check script one of them;
* the Frobenius σ table (Results 5) has no writer: the script that produced it was overwritten by the exponent version of `bianchi9_kasner.py`, and the restored `bianchi9_kasner_transfer.py` still calls itself `bianchi9_kasner.py` in its docstring and run line;
* `bianchi9_six.py` has no Noether check (the identity α′E_α + β′·E_β = d(Con)/dt was verified in the derivation session, not in the script) and solves for β⁽⁶⁾ numerically without printing its coefficient;
* `kasner_linear.pkl` is written by `bianchi9_kasner.py` (line 101) and read by nothing.

## History

* 5 Oct 2026 — derivation, codegen, first-look sweep; symplectic form; frozen-frame transfer matrix.
* 6 Oct — `--wkb` frame and "clean points" (withdrawn the same evening by the Bogoliubov quadrature); Frobenius table and "Kasner relatively attracting" (withdrawn 10 Oct); healthy probes; two-ghost derivation; Kasner-epoch onset with integrated exponents. A parallel branch merged at 21:58 added `bianchi9_kasner_transfer.py` (its `bianchi9_kasner.py` had been overwritten) and `bianchi9_onset_frozen.py`. 22:10 nonlinear seeds, literature resolved; 22:25 prefactor map.
* 10 Oct — closure (the four-item answer); the 6 Oct Kasner conclusion superseded.
* 12 Oct — README-against-code pass; mismatches held.
* 23 Oct — this rewrite.
* 24 Oct — code comments brought up to the README.
