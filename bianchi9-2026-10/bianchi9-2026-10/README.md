# bianchi9-2026-10 — Weyl ghost in Bianchi IX minisuperspace of quadratic gravity

Started 2026-10-05. Same question as the Lee–Wick Yang–Mills paper (root of this
repository), asked of the physical Weyl ghost: in the diagonal Bianchi IX
minisuperspace of

    S = ∫ d⁴x √−g [ R/2 − C_{abcd}C^{abcd}/(4M²) ]        (M̄_P = 1, no R² term yet)

the vacuum-GR mixmaster solutions form an invariant submanifold. Is it a repeller
for the ghost, and at what rate per bounce as a function of curvature/M²?

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

## Next

* Next-order adiabatic frame (first WKB correction to the frozen modes), or a wall-free
  control run at the same M to calibrate the O((H/M)²) frame artefact; then r(H_b/M) over three
  decades decides power law vs exponential.
* Kasner-epoch linear stability of the ghost vs shear/M (frozen coefficients, then Floquet
  along a full Kasner segment): where exactly does the instability switch on, and is it a
  Krein collision?
* Gain-per-bounce law in M/H_b over 3 decades once the action is exact; check e^{−cM/H_b}.
* Port the rhs to C (same pieces) for ensembles over bounce parameters (u, wall type).
* Literature check before anything is written: Cotsakis–Demaret–De Rop–Querella
  (1993) and Barrow–Cotsakis on mixmaster in fourth-order gravity — compare their
  minisuperspace Lagrangian with `L` in the pkl.
