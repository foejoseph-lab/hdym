# axion_su2/ — the isotropic SU(2) background of axion-SU(2) cosmology, with the Lee–Wick ghost

Plan item 5, decided 12 Oct 2026: one new background, the axial (axion-SU(2)) periodic one,
no lattice. The homogeneous isotropic SU(2) field A^a_i = δ^a_i a(t)Q(t) is the background of
axion-SU(2) (chromo-natural) inflation. After inflation, once the axion decouples, it oscillates
under its quartic self-coupling (Fujita, Murai, Nakayama, Yin, arXiv:2312.06889, whose radiation-era
equation is Q̈ + 3HQ̇ + m²Q + 2g²Q³ = 0; they assume, and do not show, that the field stays
isotropic). In this repository's units (g = M = 1) and in flat space it is the orbit x = y = z of
the 12-d Lee–Wick system of `ym/`, s̈ = −2s³. The purpose is the Yang–Mills paper's wording: which
of its features belong to chaotic channels and which survive on a periodic background.

Stages: (1) flat space, no axion: every linear sector on the isotropic orbit — this folder so far;
(2) expansion: Yang–Mills is conformal but the Lee–Wick term is not, so the ghost mass picks up
a(t) in conformal time — needs a derivation, not an assumption; (3) the axion stage (the φF F̃ drive
during reheating; Fujita et al., appendix A).

## Summary (read this first; written 2026-10-14 after the first stage-3 run)

**Question.** The same question as the Yang–Mills paper (`paper/`, `ym/`) and the Bianchi IX notes
(`bianchi9/`), asked of the one periodic background in cosmology that is an SU(2) Yang–Mills
configuration: the isotropic field A^a_i = δ^a_i Q(t) of axion-SU(2) inflation after inflation,
when Q oscillates under its quartic self-coupling and the axion φ oscillates at its mass m_φ and
drives it through φFF̃. Fujita et al. assume the field stays isotropic. With the Lee–Wick ghost of
mass M in the gauge sector, does it, and at what rate does anisotropy grow as a function of M?

**What closes item 5 (the exit criterion).** A growth-rate law for the homogeneous anisotropic sector
on this background as a function of M, measured and derived, with the two controls that bound it
(no ghost; ghost heavier than the drive), in the same shape as the Bianchi IX closure: a threshold
in M, the law on each side, the e-fold count through reheating, and what is scenario-dependent.
Hubble damping and the coupled axion (their A6) are listed as "what would change it" and are done
only if the flat-space result is borderline. It is not borderline.

**Answer (flat space, linear, homogeneous perturbations, M → ∞ background Q(t), one benchmark —
theirs; everything else is parameter-free in g, f, m_φ; `code/axion_aniso.py`, log
`results/axion_aniso.txt`, run in the cloud session 14 Oct, 374 s).**

1. **Threshold M = m_φ.** The anisotropic direction (the two permutation-odd modes, degenerate)
   obeys ε⁗/M² + ε″ + κ(t)Q(t)ε = 0 with κ = (g/f)φ̇. At second order in the drive its slow part
   obeys, in units g = M = 1 with ω = m_φ/M and X the slow amplitude of Q,
   s″ = + (ϕ²X²/2) · ω²/(ω² − 1) · s, ϕ = φ_amp/f.
   Two terms of the same order make it: the fast response of ε through the ghost pole,
   s″ = +c₁²/(2ω²(ω²−1)) s with c₁ = ωϕX (measured to 1e-4 on the frozen-Q problem, Part B3), and
   the jitter of Q itself, δQ = ϕX² sin(ωt)/ω, whose product with κ averages to −ϕ²X²/2 (Part D).
   So: **M > m_φ, stable** (the Kapitza stabilisation of the sign-changing coefficient, through the
   ghost, over-compensates the jitter term); **M → ∞, marginal** — the two terms cancel exactly, and
   the M → ∞ Hill equation on the same background shows drift only (the no-ghost control, Part D);
   **M < m_φ, tachyonic** with rate
   μ = gQ ϕ / √(2(1 − r²)), r = M/m_φ,
   not suppressed at small M, diverging at r → 1 (the ghost / free-branch combination resonance).
2. **Rate per Hubble time** at the end of inflation: μ/H_end = m_Q,end ϕ / √(2(1 − r²)), i.e.
   ≥ 3.5 ϕ for every M < m_φ at their m_Q,end = 5; 4.1 ϕ at r = ½, 8.1 ϕ at r = 0.9. Measured
   against the direct two-frequency integration at reduced scale separation (m_φ/gQ = 50 and 200
   instead of their 10⁵): growth per Q cycle 3.270 / 3.286 / 2.341 / 3.024 / 3.802 (ω = 2, 2, 2, 3,
   1.5; ϕ = 2, 2, 1, 2, 2) against the reduced slow equation 3.250 / 3.285 / 2.341 / 3.024 / 3.800
   (Part D, table). The frozen-Q (ghost-pole) term alone gives 1.80 / 1.81 / 1.28 / 0.91 / 2.29 there:
   the jitter term is the larger half.
3. **E-folds through matter-dominated reheating** (adiabatic estimate, not an FRW run): with
   Q ∝ 1/a and φ̇ ∝ a^{−3/2} the rate falls as a^{−5/2} and the total is N = (rate/H)_end, with the
   Q-cycle average ⟨|cn|⟩ = 0.599 and the turning-point factor 1.12 of the slow equation:
   N = 1.12 × 0.599 × m_Q ϕ/√(2(1−r²)) = 2.4 ϕ (r → 0), 2.7 ϕ (r = ½), 5.4 ϕ (r = 0.9), 9.8 ϕ (r = 0.97).
   ϕ (the axion displacement at the end of inflation, in units of f) is not stated by Fujita et al.;
   everything is linear in it.
4. **What the frozen picture gets wrong.** Frozen coefficients say the Yang–Mills branch
   ω² ≈ κQ is tachyonic half of every axion oscillation. It is not: κQ/m_φ² ~ 10⁻⁵ is five orders
   below the Kapitza threshold 0.454 (Mathieu b₁(q) = 0 at q = 0.908, reproduced by the 4-d system,
   Part B2). The instability is second order and goes through the ghost. The ghost's own parametric
   tongue at M = m_φ/2 (ghost-pair rate c₁/4, relative width qϕ ~ 10⁻⁵) lies below the band rate and
   is not the leading feature.

**Contrast.** Yang–Mills mechanics (`paper/`): chaotic channel visits kick the ghost with a power
law at every energy. Bianchi IX (`bianchi9/`): the smooth bounce is quiet, the ghost wakes at
H ≈ 0.2 M, and the Einstein branch is a repeller by a power law only. Here: a periodic background
with a fast drive, and the deciding fact is the sign of the ghost's response to the drive — above
its pole (m_φ > M) the response flips sign and turns the Kapitza stabilisation into a tachyon,
exponential, with a parameter-free e-fold count; below it nothing happens. Three backgrounds,
three mechanisms (a spectral tail; coupling parity at a turning point; a pole sign), one
pattern: the ghost is excited only where the background has a frequency that reaches M — the
channel chirp in Yang–Mills, H ≈ 0.2 M in Bianchi IX, m_φ > M here.

**Figures.** None yet; the numbers are in `results/axion_aniso.txt`.

**Status.** Working notes, 12–14 Oct 2026, in chronological order below. Flat space; the FRW
integration with the a^{−5/2} envelope, the coupled axion (their A6), the isotropic
(permutation-even) tangent, and the off-diagonal sector with the drive are not run. Each would
change a prefactor, not the threshold, unless it introduces a new scale between H and m_φ.

## Layout

```
code/
  iso_floquet.py     monodromy of every linear sector on x = y = z over E: diagonal Yang–Mills,
                     diagonal ghost, pure-YM off-diagonal blocks, and the off-diagonal Lee–Wick blocks
                     of ym/rung1_system.pkl (c_j = 1, the theory; c_j = 0 as a diagnostic).
                     Reads ../../ym/rung1_system.pkl.                       (float; added 2026-10-12)
  rung1_frw.py       ym/hdym_rung1.py with the Lee–Wick term / a(η)²: the 24-dim first-order system of
                     the off-diagonal components with a, a′, a″ as background symbols; checks that a = 1
                     reproduces ym/rung1_system.pkl exactly and that the gauge zero modes survive for
                     any a(η). Writes ../rung1_frw_system.pkl.                 (sympy; added 2026-10-14)
  frw_offdiag.py     the off-diagonal blocks of rung1_frw_system.pkl on the isotropic orbit: frozen-a
                     Floquet (checked against iso_floquet.py, extended to a = 0.25) and radiation-era
                     crossings a = η/η_s with c_j = 0 as the positive control.   (float; added 2026-10-14)
  axion_derive.py    stage 3, step 1: what φFF̃/(4f) adds to the homogeneous Lee–Wick system, seven
                     PASS/FAIL checks (reduction to (g/f)φ̇ det A, Gauss law, truncation, auxiliary
                     form, units, linearisation on x = y = z, FRW).            (sympy; added 2026-10-14)
  axion_aniso.py     stage 3, first run: the anisotropic sector ε⁗/M² + ε″ + κQε = 0 as the 4-d
                     first-order system. A: checks (frozen limit, Hamiltonian structure, M → ∞ against
                     Hill/Mathieu); B: the map of the frozen-Q periodic drive (pair analysis of the
                     monodromy, the Kapitza threshold, the ghost-pole law, the resonance near m_φ = M);
                     C: Fujita's numbers in g = M = 1 and the law per Hubble time; D: direct two-frequency
                     integration on the driven background with the no-ghost and M > m_φ controls.
                                                                              (float; added 2026-10-14)
results/
  iso_floquet.txt    its output (run in the cloud session of 12 Oct, Python 3.13; seconds per energy)
  frw_identity.txt   output of code/frw_identity.py
  frw_ghost.txt      output of code/frw_ghost.py
  rung1_frw.txt      output of code/rung1_frw.py (20 s)
  frw_offdiag.txt    output of code/frw_offdiag.py (200 s)
  axion_derive.txt   output of code/axion_derive.py (3 s), all checks PASS
  axion_aniso.txt    output of code/axion_aniso.py (374 s, cloud session 14 Oct, Python 3.13, scipy 1.18)
rung1_frw_system.pkl the FRW linearisation (K, Gauss matrix G, zero modes, L2; lambdify over
                     x,y,z,px,py,pz,c_j,a,a1,a2); regenerated by code/rung1_frw.py in 20 s
```
Also in `code/`: `frw_identity.py` (sympy: the Lee–Wick term on flat FRW in conformal time is the
flat one with M → M a(η)) and `frw_ghost.py` (the isotropic diagonal ghost crossing its band in a
radiation era). Both added 2026-10-12, run in the cloud session (seconds). The 14 Oct scripts were
run in the cloud session too (Python 3.13, sympy 1.14).

## Results so far (float, flat space, homogeneous perturbations only)

On x = y = z, s = A cn(√2 A t | k² = ½), A = (2E/3)^{1/4}, period 4K(½)/(√2 A):

- **Diagonal ghost.** Hess V = 2s² × (all-ones matrix). The two anisotropic ghost modes are free
  (w″ = −w). The isotropic mode obeys w″ = −(1 + 6s²)w, which in τ = √2 A t is the Lamé equation
  with n = 2, k² = ½, h = 3 + 1/(2A²). Its band edges are h = 1.268, 1.5, 3, 4.5, 4.732, so for
  h ≥ 3 the only instability is h ∈ (4.5, 4.732), i.e. **E ∈ (1/8, 1/6)** exactly. The scan finds
  one unstable interval, E ∈ [0.1263, 0.1661] on its grid, maximum rate 0.028.
- **Off-diagonal ghost (c_j = 1).** Floquet-stable at every E of the scan, 0.005–10. The three
  blocks agree to 1e-13, as isotropy requires. (With c_j = 0, the diagnostic, it is unstable for
  E ≳ 1.5.)
  Status (2026-10-14): extended by the frozen-a scan of `frw_offdiag.py` (flat space with M → a is
  the same problem at E = (3/2)/a⁴, checked to 1e-14 at five spot values): stable at every E of
  0.0024–384, 120 log-spaced points; with c_j = 0 the one unstable interval is E > 1.51
  (a < 0.998), max rate 0.21 per unit time at a = 1.
- **Yang–Mills perturbations, diagonal and off-diagonal.** Marginal (all multipliers 1, Jordan
  blocks): no exponential instability of the isotropic configuration among homogeneous
  perturbations at linear order. This is the assumption Fujita et al. make; at linear order and
  for homogeneous modes it holds marginally, not robustly. Inhomogeneous modes are not covered
  here (on a closed universe, Kumar, Lechtenfeld, Picanço Costa, arXiv:2102.08401, find bands of
  instability for spatial harmonics).

Compare: on the diagonal orbit x = y, z = 0 the paper finds several diagonal-ghost tongues
(sec:floquet); on the DR orbit the off-diagonal ghost is unstable over about half of E = 0.05–1
(`islands/README.md`). On the isotropic orbit there is one diagonal-ghost band and none for the
off-diagonal ghost.

## Stage 2: expansion (12 Oct, first pass)

- **Derivation.** On ds² = a²(dη² − dx²), ∇_μF^{μν} = a⁻⁴ × (flat divergence), √−g J·J =
  a⁻² × (flat J·J) and √−g F·F = flat F·F (`frw_identity.py`, sympy, all True). So in conformal time
  the Lee–Wick Lagrangian is the flat one with M → M a(η). In the auxiliary-field form
  (−J²/2M² = stat_w [w·J + M²|w|²/2]) M appears only algebraically, so the linear diagonal ghost
  obeys w″ = −(M²a² + Hess V)w on the unchanged conformal-time background.
- **Consequence on the isotropic orbit.** With comoving amplitude A, the frozen-a Lamé parameter is
  h = 3 + M²a²/(2A²): as the universe expands h rises, so the isotropic ghost crosses its single
  band once, while the physical energy falls through (1/6, 1/8), and is stable afterwards. The
  anisotropic ghost modes are never driven.
- **How much it grows in the crossing** (`frw_ghost.py`, radiation era a ∝ η): e-folds ≈ μ̄ × (time
  in the band), μ̄ = 0.041 per unit conformal time (A = 1), i.e. e-folds ≈ 1.6×10⁻³ M/H with H the
  Hubble rate at the band: 16.5 e-folds at H/M = 10⁻⁴, about 1 at H/M = 3×10⁻³, fractions of one
  for H/M ≳ 0.03. So a slow crossing (H ≪ M) amplifies the ghost a lot, linearly; a fast one
  does not. What H/M is at the crossing depends on whether the gauge field dominates the energy
  density, which this does not decide.
- Not covered: the off-diagonal blocks with M → M a(η) (ym/rung1_system.pkl is built with M = 1
  in fourth-order form; needs regenerating in the auxiliary form), backreaction, other eras.
  Status (2026-10-14): the off-diagonal blocks are done, in the fourth-order form (next bullet);
  backreaction and other eras remain open.
- **Off-diagonal blocks with M → M a(η)** (`rung1_frw.py`, `frw_offdiag.py`, 14 Oct). The
  fourth-order form was regenerated with the Lee–Wick term / a(η)², a(η) a free function: the
  u⁗ coefficients become −I/a², and a, a′, a″ enter the first-order matrix (a″ linearly, a′ up to
  a′², each over a power of a — the derivatives of 1/a²), the linearised Gauss matrix picks up a′.
  Checks (all PASS): a = 1, a′ = a″ = 0 reproduces `ym/rung1_system.pkl` symbolically, K and G;
  the gauge zero modes u = θ × A solve the system and the Gauss law for arbitrary a(η).
  Radiation era a = η/η_s, comoving amplitude 1, block (u12, u21), ln σ_max of the transfer matrix
  (raw state, no WKB normalisation; two background phases):
  - a from 0.3 to 3 (E_phys from 185 to 0.0185), c_j = 1: 7.5, 7.0, 7.1, 8.9, 9.8, 11.1, 12.0 for
    η_s = 1, 3, 10, 30, 100, 300, 1000 (H/M at a = 1 from 1 to 10⁻³; phase 0). From η_s = 30 to
    1000 that is +3.1 while ln(η span) grows by 3.5: growth ∝ (time)¹, the linear Jordan drift the
    flat Floquet already has (multiplier 1, Jordan block), not exponential. The O(1) base (7.5 at
    the fastest crossing) is the change of scale of the derivative components over a = 0.3 → 3.
  - the same with c_j = 0, the control: 7.4, 6.9, 7.2, 10.4, 18.1, 42.7, 127.6 — proportional to
    the conformal time spent in its frozen band a < 1 (η from 300 to 1000 at η_s = 1000):
    (127.6 − 7)/700 ≈ 0.17 per unit η against the frozen maximum 0.21 — the method sees an
    instability when one exists.
  - a from 1.4 to 2.3 (the diagonal ghost's band window): c_j = 1 gives 2.6 … 8.6, c_j = 0 gives
    2.6 … 7.7 over the same η_s, both ∝ (time)¹; neither has a band there.
  - the linearised Gauss charge is conserved along the FRW evolution to 1e-13 (relative).
  So on the isotropic orbit the expansion does not drive the off-diagonal ghost: no frozen band
  exists for c_j = 1 at any a in [0.25, 5], and the non-adiabatic terms (a′, a″) add nothing beyond
  the secular drift, at crossing rates from H ~ M to H ~ 10⁻³ M. The one ghost instability of this
  background under expansion remains the isotropic diagonal mode's single band crossing above.
  Not covered: a WKB-normalised measure for the fourth-order block (the raw ln σ_max carries the
  adiabatic O(1) part), other eras, backreaction.

## Stage 3: the axion stage (14 Oct, derivation only; no run yet)

What Fujita et al. appendix A is: a check that during matter-dominated reheating the inflaton's
perturbative decay to gauge bosons and stimulated emission are negligible (their A1–A5), and a
numerical run of the coupled isotropic system (their A6–A7, variables x = m_φ t, h = H/m_φ,
ϕ = φ/f, q = gQ/m_φ) showing the coupling barely moves q for their parameters (m_Q,end ≈ 5,
g = 10⁻¹²). Homogeneous isotropic mode only; nothing on anisotropic or inhomogeneous modes.

**Derivation** (`axion_derive.py`, sympy, all PASS; the Lee–Wick term is built from F alone and is
untouched; φ(t) is a given homogeneous function, its own equation — their (12) — not used):

- Homogeneous reduction: (1/4f) φ F·F̃ = ±(g/f) φ d(det A^a_i)/dt, so up to a total derivative
  L_axion = (g/f) φ̇ det A, with the sign fixed by matching their (13) at H = 0 (it is the choice of
  ε⁰¹²³; the oscillating stage does not depend on it). The term contains no A_0: the Gauss law is
  unchanged. The diagonal ansatz stays a consistent truncation, on which L_axion = (g/f) φ̇ xyz.
- Auxiliary form: the axion is the time-dependent potential −κ(t) xyz, κ = (g/f) φ̇, added to V
  in q = p − w, and nothing else; the Lee–Wick term keeps ∇V (not ∇V_tot). Verified against the
  4th-order form.
- Units g = M = 1 (time 1/M, fields M/g): the one new dimensionless quantity is κ/(gM) = φ̇/(fM)
  = 2ξ H/M with their ξ = φ̇/(2fH); their (13) at H = 0 reads u″ + 2u³ = (φ̇/(fM)) u².
- On FRW in conformal time the term is unchanged, (g/f) (dφ/dη) det A with comoving fields (the
  Chern–Simons density is conformal), so stage 2's M → M a(η) carries over.
- **The drive sources the ghost.** On shell w = (q̈ + ∇V)/M², and the driven background has
  q̈ + ∇V = κ ∇(xyz) ≠ 0: the ghost-free manifold w = W = 0 is not invariant once the axion is
  on, and the isotropic background is (Q, Q, Q; w₀, w₀, w₀) with w₀ = κQ²/M² + O(M⁻⁴). Its
  4th-order equation is printed in `results/axion_derive.txt`; as M → ∞ it is their (13).
- **Anisotropic sector on the isotropic background** (the two permutation-odd directions,
  degenerate, w₀ drops out exactly):
  ε⁗/M² + ε″ + κ(t) Q(t) ε = 0,
  a Pais–Uhlenbeck pair. Frozen coefficients: ω² = (M²/2)(1 ± √(1 − 4κQ/M²)); Yang–Mills branch
  ω² ≈ κQ, tachyonic whenever κQ < 0 (half of every axion oscillation); ghost branch
  ω² ≈ M² − κQ; the pair merges into complex ω² for κQ > M²/4. Without the axion this sector is
  exactly free (stage 1). In their variables κQ/M² = 2ξ m_Q (H/M)², and at the end of inflation
  ξ ≈ m_Q + 1/m_Q ≈ 5, so κQ > M²/4 there needs only H/M ≳ 1/√(8 ξ m_Q) ≈ 0.07.
  Frozen-coefficient statements are heuristics; the runs will say what the time dependence does.
- Isotropic tangent: printed; as M → ∞ it is ε″ + (6g²Q² − 2κQ) ε = 0.

**Conversion to Fujita et al.'s variables** (their (17), (A8): m_Q = gQ/H, ξ = φ̇/(2fH), x = m_φ t,
h = H/m_φ, ϕ = φ/f, q = gQ/m_φ; ours: g = M = 1, time in 1/M, fields in M/g): ω ≡ m_φ/M; Q_ours =
q m_φ/M; H/M = h m_φ/M; κ/(gM) = 2ξ H/M; m_Q = Q_ours/(H/M); the Hill coefficient of the anisotropic
equation κQ/M² = q (dϕ/dx) (m_φ/M)²; in the oscillating stage with V(φ) quadratic, κ = −g m_φ ϕ_amp
sin(m_φ t) and κQ/M² = (m_φ/M)² q ϕ_amp cos(·). Their (13) at H = 0 is Q″ + 2Q³ = κQ². Their Fig. 5
benchmark (used in `axion_aniso.py`): g = 10⁻¹², m_Q,end = 5, m_φ = 5×10¹⁰ GeV, f = 0.004 M_Pl,
so H_end = gM_Pl/m_Q² = 9.6×10⁴ GeV (their (33)), Q_end = M_Pl/m_Q (their (32)), gQ_end = 4.8×10⁵
GeV, q_end = 9.6×10⁻⁶, h_end = 1.9×10⁻⁶, ξ_end = m_Q + 1/m_Q = 5.2. Their m (the dark-matter mass)
is irrelevant at this epoch and set to 0; ϕ_amp is not stated by them and is taken as 1.

**First run (14 Oct, late): the anisotropic sector** (`axion_aniso.py`, `results/axion_aniso.txt`;
the Summary at the top quotes it). What was planned: one toy drive κQ = c₀ + c₁ cos ωt as positive
control, then the physical case with M as the scan parameter. What the scales dictated: with their
numbers m_φ : √(κQ) : gQ : H = 1 : 3×10⁻³ : 10⁻⁵ : 2×10⁻⁶, so the axion oscillates 10⁵ times per
oscillation of Q, and the sinusoidal drive with frozen Q *is* the physical fast problem; the slow Q
cycle sweeps its amplitude. c₀ ≠ 0 stays a control and is never quoted.

- *Method.* The 4-d monodromy over one drive period; multipliers come in reciprocal pairs
  {λ, 1/λ} (Hamiltonian), each read off x = λ + 1/λ from the characteristic polynomial, because the
  free Yang–Mills branch sits at a Jordan block (x = 2) where eigenvalues are √-sensitive to roundoff
  and x is not (Part A2: det = 1 to 1e-11, pairs to 1e-9). ω = 1/n exactly is avoided: there the ghost
  pair is also at x = 2 and the pair analysis is ε^{1/4}-sensitive. Part A3: as ω → 0 at fixed
  Mathieu (a, q) the 4-d rate agrees with the 2-d Hill equation (difference ∝ ω²) and with
  scipy's a_n, b_n verdicts at six (a, q) points.
- *The frozen-Q map* (Part B). Below ω = 1: stable except the ghost's 2/n tongues, O(c₁ⁿ), and, for
  c₁/ω² above the Kapitza threshold, the frozen-tachyon regime: threshold 0.4521 (ω = 0.1), 0.4539
  (ω = 0.03) against Mathieu's b₁(q) = 0 at q = 0.9080, c₁/ω² = 0.4540. Near ω = 1 the combination
  resonance (ghost at 1, free branch at 0): a Krein bubble below 1 and a peak above it, peak rate
  ∝ c₁^0.65, bubble width ∝ c₁^0.69 (c₁ = 0.001–0.03). Above ω = 1 the whole range is unstable by
  the ghost-pole law, measured rate / c₁(ω√(2(ω²−1)))⁻¹ = 1 to 9×10⁻⁵ for ω ≥ 1.5, c₁ ≤ 0.01
  (B3); below ω = 1 the same term is a frequency shift, x − 2 = −(νT)² with ν² = c₁²/(2ω²(1−ω²)),
  to three digits.
- *The physical point* (Part C). The frozen-Q rate along r = M/m_φ at the physical c₁ = 9.6×10⁻⁶ ω²:
  6.8×10⁻⁶ (r = 0.1), 7.1×10⁻⁶ (0.3), 9.6×10⁻⁶ (0.5, the tongue centre), 9.5×10⁻⁶ (0.7), 1.56×10⁻⁵
  (0.9), 2.79×10⁻⁵ (0.97) per unit M⁻¹, 0 for r = 1.03–2.5; linear in c₁ from 1000× down to 1× the
  physical value (rate/c₁ = 0.06671 at r = 0.3, 0.4852 at 0.7, the law 0.06671, 0.48517). The ghost
  tongue at r = ½: ghost-pair rate c₁/4 at the centre, falling to the band value 2.041 (per Hubble
  time) at ±0.75 qϕ — a separate pair, the larger counts, they do not add. Near r = 1 at the
  physical c₁: unstable for ω ∈ [0.9996, 1.0030], peak 2.36×10⁻⁴ per unit M⁻¹ at ω = 1.0002, i.e. 123
  per Hubble time, indicative only (the full law diverges there).
- *The direct check and the jitter term* (Part D). The two-frequency integration on the driven
  background Q″ + 2Q³ = κQ², at m_φ/gQ = 50 and 200, 30 Q periods, growth slope per cycle fitted over
  cycles 10–30: 3.270, 3.286, 2.341, 3.024, 3.802 for (ω, ϕ) = (2, 2), (2, 2), (2, 1), (3, 2),
  (1.5, 2) against the frozen-Q prediction 1.798, 1.813, 1.282, 0.906, 2.292 — short by a factor
  1.8–3.3 that does not depend on the scale separation. A synthetic test (single drive with the
  cn envelope, sign-changing or rectified; not in the log) agreed with the frozen-Q law to three
  digits, so the difference is the background: Q carries the Kapitza jitter δQ = ϕX² sin(ωt)/ω of its
  own response to κ, and κ δQ averages to −ϕ²X²/2, a direct tachyonic term of the same order as the
  ghost-pole term. The two add to s″ = (ϕ²X²/2)·ω²/(ω²−1)·s; the reduced slow equation with this
  μ(t) gives 3.250, 3.285, 2.341, 3.024, 3.800. Controls in the same table: ω = 0.8 and 0.6 (M > m_φ)
  slopes 0.009 and 0.014; the M → ∞ Hill equation on the same background 0.044 — the ln t drift of
  a free particle over cycles 10–30 is 0.055. The background's period agrees with the Kapitza-averaged
  quartic Q″ = −(2 + ϕ²)Q³ at the slow amplitude X = Q_amp(1 + 2ϕ²/(2 + ϕ²))^{1/4} to 0.1 % (the
  jitter velocity at t = 0 sets X), which is what Fujita's "the coupling barely moves q" looks like
  at ϕ = O(1): the shift is in the effective coupling, not in q.
- *Why there is no instability without the ghost.* As ω → 0 (M → ∞) the ghost-pole term becomes
  the ordinary Kapitza stabilisation −c₁²/(2ω²) = −ϕ²X²/2 and cancels the jitter term exactly:
  the anisotropic direction is marginal, as in stage 1 without the axion, and as the Hill control
  shows. This is the second-order linear form of the isotropisation Wolfson et al. find for the
  inflationary (fixed-sign φ̇) background.

Not run: the isotropic (permutation-even) tangent; the off-diagonal sector with the drive; the FRW
integration (the a^{−5/2} estimate above assumes the flat-space rates follow the envelope); the
coupled axion (their A6); ϕ_amp as a parameter (everything is linear in it). Not covered by the
method: the growth once ε is no longer small (nonlinear), and inhomogeneous modes.

**Literature pass (14 Oct), post-inflationary anisotropic stability of the isotropic SU(2)
background.** Nothing found that treats the anisotropic homogeneous perturbations while the axion
oscillates; the question above is open in the literature as far as this pass reached.
- During inflation (φ̇ of fixed sign): Wolfson, Maleknejad, Murata, Komatsu, Kobayashi,
  arXiv:2105.06259, Bianchi I with A^a_i = ψ e^{β}, β traceless, find numerically that the system
  always isotropises (a few e-folds in their examples), and note that the Chern–Simons term depends
  only on ψ = (det A)^{1/3}. That is the nonlinear form of CHECK 1 here (det A is the whole term);
  the linearised mass κQ of the Cartesian anisotropic modes is the second-order dependence of det A
  on the anisotropy, which their variables put into ψ. With κ > 0 the frozen Yang–Mills branch
  ω² = κQ is stable, consistent with their isotropisation. Earlier: Maleknejad, Erfani,
  arXiv:1311.3361 (same conclusion, chromo-natural model); Murata, Fujita, Kobayashi,
  arXiv:2211.09489 treat SU(N) ⊃ SU(2) embeddings, a different question.
- After inflation: Fujita et al. (2312.06889) assume isotropy (above). Khan, Pirzada, Mustafa,
  arXiv:2604.07044, decouple the axion first and decompose the k = 0 fluctuations of the massive
  (symmetry-broken) phase under the diagonal SO(3) into 1 ⊕ 3 ⊕ 5 with frequencies
  m² + 6g²Q², m² + 2g²Q², m² for trace, antisymmetric, traceless-symmetric: at m = 0 the quintet
  is gapless, which is the stage-1 result here (diagonal anisotropic modes free, isotropic mode
  6Q²). No Floquet or resonance analysis there. Daniel, Kamali, arXiv:2603.12339, keep the
  oscillating axion coupled in the background during the first e-fold of reheating, with
  Q″ + 2𝓗Q′ + (𝓗′ + 𝓗²)Q + 2a²g²Q³ = (gλ₁/f) a χ′ Q² in conformal time — with X = aQ this is the
  flat equation with κ = (gλ₁/f) dχ/dη, the same as CHECK 7 here — but study only tensor modes.
  Adshead, Giblin, Weiner, arXiv:1708.02944 (non-Abelian gauge preheating, lattice) use a dilaton
  coupling e^{φ/M} Tr F², not φFF̃, and no inflationary SU(2) background.
- The anisotropic equation ε⁗/M² + ε″ + κ(t)Q(t) ε = 0 is a Pais–Uhlenbeck oscillator with one
  time-dependent coefficient; Kuwabara, Yumibayashi, Harada, arXiv:1503.03657 note that the
  time-dependent PU oscillator does not decompose into two independent oscillators. So the run
  integrates the 4-d first-order system directly (as `frw_offdiag.py` does for the fourth-order
  off-diagonal blocks); nothing to borrow.

Not done: any run; Hubble damping and the decaying φ̇ amplitude;
the coupled axion (their A6); the off-diagonal sector with the drive (the Gauss law is unchanged,
so `rung1_frw.py` needs only the extra term (g/f) φ̇ det A).
Status (2026-10-14, late): the anisotropic-sector run is done (section "First run" above and the
Summary); Hubble damping is an adiabatic estimate there, not a run; the rest stands as listed.

## Status / open

- [ ] Stage 2: derivation done and the diagonal ghost's band crossing measured (above). Open:
      the off-diagonal blocks with M → M a(η); what H/M is at the crossing in a given scenario.
      Status (2026-10-14): the off-diagonal blocks are done (no instability from expansion, radiation
      era, a in [0.3, 3]); what H/M is at the crossing stays open (scenario-dependent).
- [ ] Stage 3: the axion stage. Status (2026-10-14): derivation done (`axion_derive.py`) and literature pass written
      (section above); the first anisotropic-sector run is next.
      Status (2026-10-14, late): the anisotropic-sector run is done (`axion_aniso.py`): threshold
      M = m_φ, the law on each side, the e-fold estimate through reheating, two controls; see the
      Summary. Open in this stage: the isotropic tangent, the off-diagonal sector with the drive, the
      FRW run, the coupled axion.
- [ ] Inhomogeneous (k ≠ 0) perturbations on the isotropic background: decide whether in scope.
