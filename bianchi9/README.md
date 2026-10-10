# bianchi9-2026-10 — Weyl ghost in Bianchi IX minisuperspace of quadratic gravity

Started 2026-10-05. Same question as the Lee–Wick Yang–Mills paper (root of this
repository), asked of the physical Weyl ghost: in the diagonal Bianchi IX
minisuperspace of

    S = ∫ d⁴x √−g [ R/2 − C_{abcd}C^{abcd}/(4M²) ]        (M̄_P = 1, no R² term yet)

the vacuum-GR mixmaster solutions form an invariant submanifold. Is it a repeller
for the ghost, and at what rate per bounce as a function of curvature/M²?

## Summary (read this first)

**Question.** Quadratic gravity, R/2 − C²/(4M²), carries a massive spin-2 ghost of mass M.
Its Einstein solutions remain exact solutions. In the chaotic Bianchi IX (mixmaster)
collapse, is the Einstein branch a repeller for the ghost — and if so, when, and how fast?

**Answer (one trajectory family, four independent methods, all code and logs here).**
1. For H < 0.2M the ghost is adiabatic. A mixmaster bounce excites it by a factor
   ∝ e^{−(π/2)M/H_b} — no more than the same bounce excites an ordinary massive scalar of the
   same mass (computed as a control). The coupling in fact *vanishes* at the turning point
   by parity: the ghost sees the background through β̇² and α̇, both even under the bounce.
2. At H ≈ 0.2M the ghost wakes, not at a bounce but in a Kasner epoch: its amplitude law
   changes from the adiabatic t^{−1/2} to log-periodic.
3. For H > 0.2M the Einstein branch is a repeller by a **power law**, δβ ∝ 1/t
   (deviation ∝ (H/M)² in ghost units), identical for all bounces; a bounce at H_b ≳ M only
   sets the prefactor, by a factor O(10). No exponential runaway anywhere.
4. For H ≫ M the literature (Querella gr-qc/9902044; Cotsakis–Demaret–De Rop–Querella 1993;
   Barrow–Cotsakis 1989) shows mixmaster chaos is structurally unstable in fourth-order
   gravity: the quadratic theory has its own monotonic asymptotics. The Weyl ghost does not
   ride the mixmaster; it ends it.

Contrast with the Yang–Mills Lee–Wick case (root of this repository), where channel visits
kick the ghost with a power law at every energy: there the coupling is linear and sharpest at
channel entry; here it is even and quietest at the bounce. Parity and derivative order of the
coupling about the turning point decide which case a model falls in.

**Figures.** `seed_growth.png` (the universal (H/M)² track through bounces at H_b/M = 0.4–2),
`kasner_onset.png` (the knee at H/M ≈ 0.2), `bog_diag.png` (the coupling dip at the bounce),
`fields_n.png` (healthy scalar and Dirac particle numbers for comparison).

**Status.** Working notes, 5–6 Oct 2026, in chronological order below, including two
intermediate claims that were later withdrawn (marked where they occur). Single wall,
single trajectory family, five Kasner angles; the two-ghost (sixth-order) case is set up but
unfinished. Everything reproduces on two machines; run commands are in each script's docstring.

## Files

* `bianchi9_derive.py` — sympy derivation. Orthonormal-frame Cartan calculus on
  ds² = −N²dt² + a²σ₁² + b²σ₂² + c²σ₃², Misner variables (α, β₊, β₋), fourth-order
  Lagrangian, N = 1 equations of motion solved for α″ and β⁗, Hamiltonian
  constraint. Eight internal checks, all passing (log: `derive.log`, 37 s).
  Output `bianchi9_system.pkl` (sympy expressions; state
  y = (α, α′, β₊, β₋, β₊′, β₋′, β₊″, β₋″, β₊‴, β₋‴)).

## Structural facts established by the script (2026-10-05)

1. Weyl² vanishes identically for isotropic a = b = c (any a(t)) and the
   Gauss–Bonnet density is a total derivative on the ansatz — so the C² term is the
   whole non-GR content; R² would only add the healthy scalaron.
2. The C² Lagrangian contains no α″ and no N″ (conformal invariance); its Hessian in
   (β₊″, β₋″) is −6e^{3α}/M² · 𝟙: both ghost polarisations are degenerate and the
   ghost kinetic matrix is isotropic in the β-plane. Together with the GR kinetic term
   +3e^{3α}|β′|² the transverse block has the Pais–Uhlenbeck form ½ẋ² − ẍ²/(2M²), whose
   frequencies are 0 and exactly M: the ghost is massive and non-tachyonic with mass M, so
   the normalisation −C²/(4M²) in the action is the one that puts the Weyl pole at M.
3. The α equation stays second order; its α″ coefficient is
   6(M² + β₊′² + β₋′²)e^{3α}/M², never zero. The β⁗ coefficient matrix has the form
   (M² + |β′|²)𝟙 − β′β′ᵀ up to positive factors, determinant ∝ M²(M² + |β′|²) ≠ 0.
   The N = 1 fourth-order system is regular on the whole phase space.
4. Phase space 10-d, one constraint → 9-d. Einstein manifold (β″, β‴ at their GR
   values) 6-d → 5-d after the constraint. Transverse: 4-d, two oscillators.
5. The Einstein manifold is invariant: on GR data the full flow reproduces the GR
   flow (Bach tensor vanishes on Einstein spaces), checked numerically on random
   constraint-satisfying points to 1e-9.

* `bianchi9_codegen.py` → `bianchi9_rhs.py` (generated, 85 kB): flat numpy/math code for
  the proper-time rhs (polynomial pieces + a 2×2 solve at run time), finite-difference
  Jacobian, Einstein-manifold values and their Jacobian, constraints. 0.4 ms per rhs call in
  pure Python; checked against lambdify to 1e-14. The closed-form β⁗ from `LUsolve` is ten
  times slower — do not lambdify `b4` directly.
* `bianchi9_test.py` — DOP853 in Misner time dt = e^{3α}dτ, tangent vector alongside,
  transverse projection δv = v_hi − (∂gr/∂y_lo)v_lo, observable ½ log J with
  J = e^{3α}(δβ″² + δβ‴²/M²)/M³ (heavy-ghost adiabatic action, crude: phase wobble ±0.2).
  `sweep_M*.log/.npz`, `sweep_first_look.png`: one GR bounce (β₊ wall, α_b = −0.38,
  H_b ≈ 2.0 in units with the initial |β′| = 1) and M = 5…200, i.e. H_b/M = 0.4…0.01.

## First numerical result (2026-10-05 evening, one bounce, one trajectory)

| H_b/M | 0.40 | 0.20 | 0.10 | 0.04 | 0.02 | 0.01 |
|---|---|---|---|---|---|---|
| log action gain across the bounce | +3.0 | +2.6 | +0.4 | +0.1 | −0.1 | 0.0 |

Invariance of the Einstein manifold holds to 1e-11 (relative) and the constraint to 1e-10 in
every run. For H_b/M ≤ 0.1 the ghost action is constant through the bounce within the ±0.2
wobble of the observable: **the mixmaster bounce is adiabatic for a heavy Weyl ghost**. The
growth at H_b/M = 0.2, 0.4 is not a step at the bounce; it is continuous and sets in after
it, as H/M → 1 along the following Kasner epoch (H ∝ e^{−3α} while M is fixed). This is the
opposite of the Yang–Mills case, where channel visits give a power-law kick at every energy.
Expected physics: the walls are smooth (exponential), so the bounce spectrum at 2M is
exponentially small in M/H_b; the Weyl ghost switches on at a definite epoch H ~ M rather
than being excited at all scales. Both halves of that sentence need to be nailed down before
anything is written: the exponential law needs a cleaner action and ≥ 3 more decades of
dynamic range, and the onset needs the Kasner-epoch linear analysis (frozen-coefficient
spectrum of the transverse 4×4 block as a function of shear/M).

Caveats: single trajectory, single wall, crude action, the on-manifold tangent (mixmaster
Lyapunov) is tracked but not yet used for anything.

## Exact ghost action and transfer matrix (2026-10-05, late)

* `bianchi9_symp.py` → `bianchi9_symp_gen.py`: Ostrogradsky canonical map z(y) (α, β, β′ | p_α,
  P₁, P₂; α″ in P₁ replaced on shell) and the pulled-back symplectic form Ω(y). Check: the
  linearised Misner-time flow preserves Ω to 1.4e-9 **for tangent vectors in the constraint
  surface** (dH·v = 0); off it the τ-flow is not Hamiltonian and Ω drifts by O(1). Same lesson
  as Rung 1.
* `bianchi9_action.py`: 10×10 tangent matrix in scaled coordinates (δβ″ ~ Mδβ, δβ‴ ~ M²δβ;
  unscaled, atol makes the integrator crawl), frozen ghost modes symplectically normalised
  (iEᴴΩE = ±1), ghost transfer matrix T (real 4×4 on the complex amplitudes). Its singular
  values e^{±r} give the squeeze r across the bounce; a symplectic ghost block has mean log s = 0,
  so the physical gain is r (random-phase action gain cosh 2r), not the mean.
  `action_M*.log/.npz`, `action_sweep.png`.

| H_b/M | 0.012 | 0.025 | 0.049 | 0.099 | 0.198 |
|---|---|---|---|---|---|
| squeeze r at α = −0.5 (just past the bounce) | 1.6e-4 | 1.0e-3 | 8.8e-3 | 0.13 | 1.2 |

Successive-doubling ratios 6, 9, 15: rising toward large H/M, the signature of a power law
(≈ 3) at small x, not of e^{−cM/H}. **But the three smallest values are upper bounds, not
measurements**: a frozen-mode frame is adiabatic only to O((H/M)²), and a symmetric ±r band of
that order is what the frame error alone produces. Defensible tonight: r ≤ 1e-3 for
H_b/M ≤ 0.025, r ≤ 1e-2 at 0.05, r = 0.13 at 0.1. The functional form below 0.05 is open.
Ghost modes satisfy |dH·e|/|dH| ∝ (H/M)² (2e-1 … 1e-3 for M = 10 … 160), as a frozen-mode
error should. Frozen frequencies 0.97 M; frozen Re λ ≈ 3H/2 is the kinematic a^{−3/2} factor,
not an instability. DOP853 crawls once H/M ≳ 0.2 (the onset regime); use an implicit method
there.

## Corrected adiabatic frame (2026-10-06)

`bianchi9_action.py --wkb`: first adiabatic correction to the ghost frame, ẽ = e + Σ_k v_k
(w_k·ė)/(λ_k − λ_pair), with ė from the projector onto the near-degenerate pair (central
difference along the flow). Check: corrected modes lie in the constraint surface an order
better (|dH·e|/|dH| 8e-3 → 1e-3 at M = 40, 3e-3 → 1e-4 at M = 80). Artefact floor measured
with a control run pointed away from the wall (`control_M80_wkb.*`, θ = 0): r_floor ≈
1.7e-5 (H/M / 0.0127)^1.8 at matched H/M. `action_M*_wkb.*`, `action_sweep_wkb.png`.

| H_b/M | 0.099 | 0.071 | 0.049 | 0.035 | 0.025 |
|---|---|---|---|---|---|
| r, corrected, at α = −0.5 | ~0.3 | 0.040 | 2.9e-3 | 5.6e-4 | 1.4e-4 |
| artefact floor | 1e-3 | 5e-4 | 3e-4 | 1.5e-4 | 8e-5 |

Clean points 0.071 → 0.049 → 0.035: r falls ×14 then ×5 per factor 1.4 in H_b/M. As a power
law the local exponent drifts 7 → 5; as e^{−cM/H_b}, c ≈ 0.44 → 0.22. A fixed power law of
spectral-tail type (the YM mechanism) is excluded. Below 0.035 the floor dominates; the
0.099 point is only marginally adiabatic (mean log s = +0.13 ≠ 0). **Theory next**: the
single-wall bounce is the Taub solution, β̈ ∝ sech²(4vt), Fourier transform ∝ ω/sinh(πω/8v),
so r should be exponential in M/v with a computable c and a power-law prefactor that explains
the drift. The three numbers above are what that formula has to hit.

## Bogoliubov quadrature and the correction of the 06 Oct table (2026-10-06, late)

`bianchi9_bogoliubov.py`: first-order adiabatic perturbation theory along the exact GR
background, B = ∫ C e^{−iφ} dτ with C = iE₋ᴴ Ω Ė₊ in a parallel-transported, symplectically
normalised ghost frame (without the parallel transport the within-pair unitary is arbitrary
and the phase-sensitive result is garbage). `bog.npz`, `bog_diag.png`.

Result: |C(τ)| has a **dip** at the bounce, not a peak (the frame depends on α̇ and β̇², which
vary least at the turning point), and the cumulative |B| shows **no step at the bounce** —
only the oscillating endpoint term of a smooth integrand that grows with H/M. The Taub
estimate agrees: the phase accumulated during one bounce is 2ω_τ/a ≈ 10 (M = 28) … 27
(M = 80), so the bounce kick is O(e^{−(π/2)M/H_b}) — 10⁻¹⁰ at H_b/M = 0.07, 4×10⁻⁴ at 0.2 —
below anything measured. **Therefore the 06 Oct "clean points" at H_b/M = 0.049 and 0.035
were second-order adiabatic frame terms, not bounce kicks, and the exponent/c values quoted
from them are withdrawn.** Evidence: the two first-order-removed estimators (transfer matrix
with --wkb; Bogoliubov with the first endpoint term removed) disagree by ×2.5–5 at the same
points, the signature of second-order residuals.

What stands, and is stronger: mixmaster bounces do not excite a heavy Weyl ghost at any
measurable level (c ≈ π/2, with O(1) corrections from Δα ≈ ½ across the bounce). The
measurable frame terms are reversible dressing, not accumulated growth. The onset of the
instability is a Kasner-epoch phenomenon at H ~ M set by shear, not by walls; the chaos is
beside the point for the ghost. The planned three-decade r(H_b/M) sweep is cancelled.

## Kasner (Bianchi I) linear stability: the "onset" is not an instability (2026-10-06, late)

`bianchi9_kasner.py`: Bianchi I limit of `L` (the e^{3α} terms), Kasner is an exact solution
(Ricci-flat ⇒ Bach-flat), scale invariance makes the linearised problem depend only on Mt and
the Kasner angle. Local Frobenius exponents σ (perturbation ∼ t^σ; toward the singularity
t → 0 the smallest Re σ dominates), from the scaled Jacobian in s = ln t, head-on angle θ = π:

| Mt | 10⁴ | 10 | 1 | 0.3 | 0.1 | 10⁻² |
|---|---|---|---|---|---|---|
| ghost pair(s) | 2 ± iMt (×2) | 2.006 ± 9.93i, 2 ± 10i | 2.20 ± 0.77i, 2 ± i | 2 ± 0.3i, 1.85 ± 0.16i | 2 ± 0.1i, 1.99, 1.12 | 2 ± 0.01i, 2.00, 1.00 |
| on-manifold | 1, 0,0,0,0, −1 | 0.99, 0⁴, −1 | 0.41, 0⁴, −1 | 0.20, 0⁴, −1 | 0.05, 0⁴, −1 | 0.00, 0⁴, −1 |

For Mt ≫ 1 the WKB amplitude correction (−½ d ln ω/d ln t) brings the ghost's Re σ from 2 to
3/2, which is the adiabatic invariant t M A² = const; −1 is the time-translation gauge mode;
physical Kasner perturbations (angle, p) sit at σ = 0. **No transverse exponent drops below 1
at any Mt**, so the ghost decays relative to physical on-manifold perturbations toward the
singularity in every regime: Kasner is relatively attracting, not repelling, for the Weyl ghost.
At Mt ≲ 1 the ghost turns from oscillatory to power-law and adiabaticity is lost — this is
what the frozen-frame observables registered as "onset" — but it is not an instability.

Literature: Toporensky & Müller, arXiv:1603.02851 (R² + Ric² gravity, χ = α/3β): Kasner
locally stable toward the singularity for χ > 4; our theory (no R²) is χ → ∞, the stable side.
Consistent. Our additions: the ghost-specific symplectic/adiabatic structure, and the walls.

Caveats: linear; one Kasner angle; "relative growth" is direction-of-time dependent (in the
expanding direction the ghost falls off more slowly than GR perturbations); the direct
tangent-matrix integration through Mt ~ 1 with a physical (not frame-based) transverse split is
still to be done; the Radau run of `bianchi9_kasner.py` crawls below Mt ≈ 1 and its
frozen-frame observable is meaningless there (mode crossing).

**Standing conclusion of the two days (to be reviewed independently before anything is written):
in R/2 − C²/4M² gravity the Einstein branch is not a classical repeller for the Weyl ghost on
the approach to the singularity: bounces do not excite it (e^{−(π/2)M/H_b}) and Kasner epochs do
not amplify it relative to ordinary perturbations. The opposite of the Yang–Mills result.**

## Healthy probes on the same background: boson and fermion (2026-10-06, late)

`bianchi9_fields.py`: massive scalar (couples to the volume only at k = 0) and massive Dirac
field (couples to the anisotropy through the spin connection: axial term A = ¼(a/bc + b/ca +
c/ab), i.e. the square roots of the curvature walls) as linear test fields on the unchanged GR
trajectory. Dirac reduction checked by Hermiticity of H_D for χ = e^{3α/2}ψ (4e-16; the wrong
Clifford convention gives 2.0 and n ≡ ½) and by E = √(m² + A²). Observable: instantaneous
adiabatic particle number n. `fields.npz`, `fields_n.png`.

| m/H_b | 0.099 | 0.071 | 0.050 | 0.035 | 0.025 |
|---|---|---|---|---|---|
| scalar n (α = −0.5) | 2.0e-5 | 2.5e-6 | 2.4e-7 | 3.3e-8 | 3.6e-9 |
| Dirac n (α = −0.5) | 5.1e-5 | 8.9e-6 | 1.9e-6 | 1.6e-7 | 1.1e-7 |

No step at the bounce for either; the scalar's n is smooth dressing with the same turning-
point dip as the ghost's coupling, the fermion's is a beat pattern passing through the bounce
unchanged. Scalings n ∝ (H/m)⁶ (scalar, second adiabatic order) and (H/m)⁴ (fermion, first
order, since A ~ H) identify both as reversible dressing, not production. **The bounce is
quiet for a healthy boson, for a healthy fermion that does see the walls, and for the ghost:
the result is a property of the background's smoothness, not of the ghost's sign.**

Planned suite (Joseph's ordering; CPU; doubles as the Python 3.13 acceptance suite):
control; ghost + light; control + light; control + ghost; all three — with "light" now
understood as a *probe*, not a background wall (a homogeneous magnetic field on diagonal
Bianchi IX did not pass d⋆F = 0 in a quick check; arXiv:1910.11970 adds 2b²e^{2α−4β₊} to the
potential without discussing the constraint — to be understood before it is used). Head-on,
wall-centre geometry throughout (worst case by the Taub argument); one m = M resonance run;
corner-era run and Kasner-angle eigenvalue loop separately.

## Particle production sanity check and the two-ghost (sixth-order) theory (2026-10-06, 21:00+)

* `fields_light.npz`: the probe masses swept down to m ≪ H_b. Boson production saturates
  (n ~ 0.7, no adiabatic vacuum for m < 3H/2 — the sqrt warning); fermion production peaks at
  m ≈ H_b and falls on both sides (conformal protection at m → 0, adiabatic at m → ∞).
  Zel'dovich–Starobinsky 1971 reproduced; the code does particle production correctly.
  Dirac large-m curves in `fields_n.png` are under-sampled (400 points) and the scalar's
  lower-left fuzz is the rtol = 1e-10 floor — cosmetic, to fix with denser output.
* `bianchi9_six.py`: S = ∫√−g [R/2 − (k₂/4)C² − (k₄/4)C□C], k₂ = 1/M₁² + 1/M₂², k₄ = 1/M₁²M₂²
  (spin-2 poles at 0, M₁, M₂ with alternating Krein sign), Bianchi I. Derivation done and
  checked: Weyl traceless; C² = C□C = 0 isotropic; √−g(C□C + |∇C|²) is a total derivative on
  the ansatz (EL vanish — validates the frame covariant derivative). Structure found:
  α-equation is **fourth** order with coefficient −6|β′|²e^{3α}/M₁²M₂² (vanishes at zero shear;
  the one-ghost coefficient 6(M² + |β′|²)e^{3α}/M² never did); the β⁽⁶⁾ coefficient matrix is
  exactly the transverse projector |β′|²𝟙 − β′β′ᵀ, so the longitudinal shear mode is of lower
  order; Noether identity α′E_α + β′·E_β = d(Con)/dt verified. **The two-ghost minisuperspace
  is a constrained higher-order system.** Recipe for the next session, on the β₋ ≡ 0 branch
  (θ = π): α⁗ from E_α, β₊⁽⁵⁾ from Con, β₋⁽⁶⁾ from E_{β₋}; 15-d closed system; then the Kasner
  exponent table vs M₁t for M₂/M₁ = 3 and the Krein-collision question. Script runs to the
  singular-matrix point; everything before it is correct and reusable.

## Kasner-epoch onset (2026-10-06, night)

`bianchi9_kasner.py` (+ `bianchi9_kasner_check.py`, `kasner_onset.png/.npz`, `kasner_linear.pkl`):
Bianchi I limit of `L` (the e^{3α}-homogeneous part; GR part checked = 3e^{3α}(−α′²+β′²)),
vacuum Kasner α = ⅓ln t, β = ⅓ln t (cos θ, sin θ) verified to solve the full equations,
linearised 10-d system in x = Mt (one parameter: H/M = 1/3x), integrated from x = 200 to 0.01.

* Symbolic t → 0: the pure-Weyl indicial determinant is s³(s−2)³·Q(s,θ) with Q ≡ 0 (conformal
  invariance: the C² operator does not see the Kasner background in the ghost directions at
  leading order); exponents come from the next order and were obtained numerically.
* Frozen spectrum in Euler form xλ: ghost pair −0.5 ± ix for H ≪ M (adiabatic t^{−1/2}); Re
  crosses 0 at **H/M ≈ 0.2**; for H ≫ M the pair becomes 0.27 ± 0.56i (log-periodic,
  t^{0.27} cos(0.56 ln t)) and 0.45 ± 0.35i. The other frozen exponents (−4.8, −2.5) are not
  realised as amplitude growth (frozen eigenvalues mislead in a non-autonomous Euler system).
* Integrated exponents for H ≫ M (fit over x ∈ [0.02, 0.1], all θ): generic data δβ₊ ∝ t^{−1.0},
  δα ∝ t^{−1.0}, δβ₋ ∝ t^{−0.2…−0.3}. GR time-shift mode gives t^{−0.99999999999} (exact
  check). GR Kasner-angle mode (the physical GR perturbation) evolves in the quadratic theory
  as δβ ∝ t^{−0.30}, δα ∝ t^{−0.81}.
* Scaled ghost amplitude |δβ|√x is flat (adiabatic) for 0.01 < H/M < 0.1 within the
  preparation transient, with the knee at 0.2–0.3 for all four Kasner angles.

**Statement:** in a Kasner epoch the Weyl ghost is adiabatic for H ≲ 0.2M; beyond that the
Einstein branch is a repeller, but only by a power law — ghost t^{−1} against the physical
GR mode t^{−0.3}, relative t^{−0.7} — not an exponential runaway. Combined with the bounce
result (kick ∝ e^{−(π/2)M/H_b}, negligible for H_b < 0.15M): for H < 0.2M nothing happens;
above it the ghost drifts off the Einstein branch polynomially in each Kasner epoch, and any
exponential behaviour must come from bounces at H ≳ M, where κ ~ H ~ M and the kick is O(1),
and from the chaos re-entering there. The step-zero literature question — does the C² term
leave mixmaster chaotic as H/M → ∞ — is therefore the deciding fact, not a formality.

Cross-reference: the dip of the coupling |C(τ)| at the bounce (`bog_diag.png`) is the same
feature as the flattening of the M = 10, 20 curves at α ≈ −0.4 in `action_sweep.png`: the
ghost frame depends on β̇₊² and α̇, even under the time reflection about the turning point, so
its rate of change 2β̇₊β̈₊ vanishes there. Parity and derivative order of the coupling about
the turning point decide whether a bounce drives a heavy mode at first order; in Yang–Mills
the coupling w·∇V(q) is linear and largest at channel entry — opposite parity, opposite
answer.

## Work from the parallel branch of the 06 Oct session (merged 21:58)

A second branch of the same session ran in the same container and left four pieces, now
recorded here (its own `bianchi9_kasner.py` was overwritten by the exponent version above
and is restored as `bianchi9_kasner_transfer.py`).

* `bianchi9_kasner_transfer.py` (log `kasner_pi.log`): the transfer-matrix squeeze r along a
  pure Kasner epoch (Radau in log t, corrected frame), independent of the exponent method:
  r < 1e-4 for H/M ≤ 0.054, r = 0.12 at H/M = 0.136, r = 1.03 at 0.34. Same onset, ≈ 0.2,
  from a different observable. Ghost frozen frequency 0.9999 M at Mt = 100 with Re λ = 3H/2
  to 4 digits (adiabatic kinematics confirmed).
* `bianchi9_fields.py` (`fields.log`, `fields.npz`, `fields_n.png`): **healthy test fields on
  the same bounce** — a massive scalar (sees only the volume) and a massive Dirac fermion
  (sees the anisotropy through the square roots of the walls, axial coupling
  A = ⅛Σ a/bc), homogeneous, no back-reaction. Adiabatic particle number n just past the
  bounce (α = −0.5):

  | m/H_b | 0.099 | 0.071 | 0.050 | 0.035 | 0.025 |
  |---|---|---|---|---|---|
  | scalar n | 2.0e-5 | 2.5e-6 | 2.4e-7 | 3.3e-8 | 3.6e-9 |
  | Dirac n | 5.1e-5 | 8.9e-6 | 1.9e-6 | 1.6e-7 | 1.1e-7 |

  Both fall by ~×8–10 per factor 1.4 in m/H_b (exponential-type, as a smooth bounce gives).
  The ghost's bounce excitation (n_ghost = sinh²r, with r ∝ e^{−(π/2)M/H_b}) is at or below
  these numbers: **a mixmaster bounce excites the Weyl ghost no more than it excites an
  ordinary massive particle of the same mass.** This is the right control and the cleanest
  one-line statement of the bounce result.
* `bianchi9_onset_frozen.py` (`onset_frozen.npz`): frozen spectrum along the real Bianchi IX
  collapse for M = 2, 5, 10, 20 (Re λ_ghost, |Im λ_ghost|, max Re λ vs H/M). Data only: frozen
  eigenvalues mislead about growth once H ≳ M (see the Kasner section); kept for reference.
* `bianchi9_six.py` (`six.log`): sixth-order Weyl gravity R/2 − (k₂/4)C² − (k₄/4)C□C with two
  ghosts at M₁, M₂ (Krein signs +, −, + for 0, M₁, M₂), Bianchi I. Derivation and checks pass
  (C□C + |∇C|² a total derivative; Einstein manifold invariant; phase space 16-d). Numerics
  stopped: the α⁽⁴⁾ coefficient is −6|β′|²e^{3α}/(M₁²M₂²), which vanishes for isotropic
  motion, so the fourth-order regularity result does not carry to sixth order — a Kasner
  point with β′ = 0 is singular for the solve. Unfinished; the question (do the two ghosts
  Krein-collide along Kasner) is open.

## Nonlinear seeded runs through bounces at H_b ~ M, and closure (2026-10-06, 22:10)

`bianchi9_seed.py` (`seed.log`, `seed.npz`, `seed_growth.png`): GR data + a ghost seed
(δβ₊″ = 10⁻⁶ M H₀), full nonlinear equations, **Radau** (implicit: 3–15 s per run where DOP853
took an hour), M = 1, 2, 3, 5 so that the bounce sits at H_b/M = 2, 1, 0.66, 0.4; run to
H/M = 8. Deviation from the Einstein manifold in ghost units (δβ″/MH, δβ‴/M²H) relative to
the seed; unseeded reference stays on the manifold to 3e-7; geometry back-reaction stays
linear in the seed.

Result: all four curves fall on **one universal track, deviation ∝ (H/M)^{2.0}** (fits 2.03,
2.04, 2.06, 2.34), whatever the bounce. In these units that is exactly δβ ∝ 1/t — the Kasner
exponent of the section above, reproduced by the full nonlinear flow through a real bounce.
The bounce only lifts the trajectory onto the track with an O(1) prefactor: ×2.5 at
H_b/M = 0.4, ×6 at 0.66, ×13 at 1, ×50 at 2 (κ ~ H ~ M, kick unsuppressed, as predicted).
After the bounce every trajectory forgets it. Ignition visible at H/M ≈ 0.2–0.3 on every curve.

**Literature (step zero, resolved):** Querella, thesis, gr-qc/9902044 (open), with
Cotsakis–Demaret–De Rop–Querella 1993 and Barrow–Cotsakis 1989: in fourth-order gravity the
mixmaster chaos is structurally unstable near the singularity; the generic approach is a
monotonic, isotropising power law. So at H ≫ M there is no billiard; the departure from the
Einstein branch found here is the transition to the quadratic theory's own asymptotics.

**Closure (four independent methods agree: transfer matrix, Bogoliubov quadrature, Kasner
exponents, nonlinear seeds):**
1. H < 0.2M — the Weyl ghost is adiabatic on the Einstein branch; a mixmaster bounce excites
   it no more than it excites a healthy scalar of the same mass (kick ∝ e^{−(π/2)M/H_b};
   coupling dips at the turning point by parity).
2. H ≈ 0.2M — it wakes in a Kasner epoch: amplitude law changes from t^{−1/2} to log-periodic.
3. H > 0.2M — the Einstein branch is a repeller by a **power law**, δβ ∝ 1/t (physical GR mode
   t^{−0.3}); prefactor set by the first bounce at H ≳ M; no exponential runaway anywhere.
4. H ≫ M — chaos ends (literature); the Weyl ghost does not ride the mixmaster, it ends it.
Open: the sixth-order two-ghost case (β′ = 0 singularity), other wall types/Kasner angles
for the prefactor (cheap now with Radau), the n_ghost vs n_scalar comparison on one axis.

## Prefactor map: Kasner angle and polarisation (2026-10-06, 22:25)

`bianchi9_seed.py --theta ... --seed-pol p|m` (`seed_th*_*.npz`), M = 2 and 5, θ = π, 2.9,
2.6, 2.3, 2.0 (head-on to glancing; glancing bounces occur later, at larger H_b/M).
* Kick at matched H_b/M ≈ 1–1.3 (growth from α = −0.2 to −0.5): ×13 (θ = π), ×12 (2.9),
  ×9 (2.6). Mild angle dependence: **O(10) at H_b ≈ M** is a fair single number.
* Polarisation: head-on, the β₋ polarisation (odd under the bounce's β₋ → −β₋ symmetry) is
  barely excited (total ×9 vs ×800 for β₊ at M = 2; ×5 vs ×420 at M = 5) and drifts with the
  slow Kasner law t^{−0.3}, not t^{−1}. The universal (H/M)² track is that of the polarisation
  **aligned with the shear**; oblique bounces mix the two (θ = 2.0: the β₋ seed ends above
  the β₊ one), so in a real mixmaster sequence the aligned law governs after the first
  oblique bounce.

## Next

* Two-ghost theory: constrained reduction (recipe above), exponent table, Krein collision yes/no.
* Suite runs 4–5 (control + ghost on one trajectory, m = M resonance); light as a probe once understood.

* (Optional, low priority) superadiabatic series to optimal order to exhibit the e^{−(π/2)M/H_b}
  kick explicitly — a methods point, not a physics one.
* Direct tangent-matrix integration through Mt ~ 1 (Radau; scaled variables) with the
  transverse / on-manifold split done against physical directions (dH·v = 0, modulo the
  time-translation mode); Kasner-angle sweep; confirm "no exponent below 1" is angle-independent.
* Independent review of the two-day conclusion by a model that has not seen the work.
* Port the rhs to C (same pieces) for ensembles over bounce parameters (u, wall type).
* Literature check before anything is written: Cotsakis–Demaret–De Rop–Querella
  (1993) and Barrow–Cotsakis on mixmaster in fourth-order gravity — compare their
  minisuperspace Lagrangian with `L` in the pkl.
