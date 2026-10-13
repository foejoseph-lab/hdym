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

## Layout

```
code/
  iso_floquet.py     monodromy of every linear sector on x = y = z over E: diagonal Yang–Mills,
                     diagonal ghost, pure-YM off-diagonal blocks, and the off-diagonal Lee–Wick blocks
                     of ym/rung1_system.pkl (c_j = 1, the theory; c_j = 0 as a diagnostic).
                     Reads ../../ym/rung1_system.pkl.                       (float; added 2026-10-12)
results/
  iso_floquet.txt    its output (run in the cloud session of 12 Oct, Python 3.13; seconds per energy)
  frw_identity.txt   output of code/frw_identity.py
  frw_ghost.txt      output of code/frw_ghost.py
```
Also in `code/`: `frw_identity.py` (sympy: the Lee–Wick term on flat FRW in conformal time is the
flat one with M → M a(η)) and `frw_ghost.py` (the isotropic diagonal ghost crossing its band in a
radiation era). Both added 2026-10-12, run in the cloud session (seconds).

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

## Status / open

- [ ] Stage 2: derivation done and the diagonal ghost's band crossing measured (above). Open:
      the off-diagonal blocks with M → M a(η); what H/M is at the crossing in a given scenario.
- [ ] Stage 3: the axion stage.
- [ ] Inhomogeneous (k ≠ 0) perturbations on the isotropic background: decide whether in scope.
