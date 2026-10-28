# axion_su2/ — the isotropic SU(2) background of axion-SU(2) cosmology, with the Lee–Wick ghost

Working notes 12–17 Oct 2026; rewritten 23 Oct 2026 from the state at that date. Status: in progress — the
coupled-axion run (Fujita et al.'s A6–A7 with the ghost) is the one open piece; the probe programme on this
background is complete, and its successor on the dilaton background (`dilaton_su2/README.md`) is closed.

## The question and where it stands

The homogeneous isotropic SU(2) field A^a_i = δ^a_i a(t)Q(t) is the background of axion-SU(2) (chromo-natural)
inflation. After inflation Q oscillates under its quartic self-coupling (Fujita, Murai, Nakayama, Yin,
arXiv:2312.06889; radiation-era equation Q̈ + 3HQ̇ + m²Q + 2g²Q³ = 0), driven through φFF̃ while the axion φ still
oscillates at its mass m_φ; Fujita et al. assume, and do not show, that the field stays isotropic. In this
repository's units (g = M = 1) and in flat space the undriven background is the orbit x = y = z of the 12-d Lee–Wick
system of `ym/`, s̈ = −2s³. The question of `paper/main.tex`, `ym/README.md` and `bianchi9/README.md`, asked of the
one periodic background in cosmology that is an SU(2) Yang–Mills configuration: with the Lee–Wick ghost of mass M in
the gauge sector, does the isotropic solution exist, is it stable, and does anisotropy grow, as a function of M?
Stages: (1) flat space, no axion; (2) expansion; (3) the axion drive during reheating; from 16 Oct, real particles
coupled to the ghost (the probe programme).

**Answer** (ω ≡ m_φ/M, r ≡ M/m_φ, ϕ the axion amplitude in units of f, Π(w) = 1/(w² − wM²) the full propagator),
for every homogeneous sector, on a background that is itself a solution of the Lee–Wick theory, with the controls
(no ghost, M → ∞; ghost heavier than the drive; the healthy twin):

1. **Anisotropic (and off-diagonal) sector: marginal at second order in the drive, for any heavy sector** — the
   two second-order terms are the same Π attached in the two orders and cancel identically; what remains is the ghost
   pair's parametric tongues at m_φ = 2M/n and the linear resonance at m_φ = M.
2. **Isotropic sector: a sharp threshold in the drive.** The drive renormalises the quartic coupling,
   λ_eff = 2 − ω²ϕ²Π(ω²); the isotropic solution is lost at λ_eff = 0, **ω²ϕ²Π = 2, i.e. ϕ² r²/(1 − r²) = 2**. Above
   it the orbit ceases to exist (no linear rate to measure); below it everything is marginal.
3. **M > m_φ: stable** in every sector for every drive tested (ϕ up to 3), as a ghost-pole statement; a second,
   pure-Yang–Mills loss exists on both sides of m_φ = M when the drive is slow relative to the Q oscillation
   (s_YM = 2ϕ gQ/m_φ ≈ 0.2–0.3). **M → ∞: marginal.**
4. **In Fujita et al.'s variables** both thresholds are conditions on m_φ/H_end (ghost-pole 300–630 for r = 0.7–0.9,
   slow-drive ≈ 130, at f = 0.004 M_Pl); at their benchmark m_φ/H_end = 5.2×10⁵ the ghost does nothing to the
   isotropic solution, for any M. Where either threshold is reachable the back-reaction on the axion is O(1), so the
   coupled system is required and the prescribed-axion result is a map, not a verdict.
5. **Real probes.** A probe at or above ≈ 3 gQ_amp is a spectator. The ghost-specific effects are one mechanism: a
   probe-induced anisotropic coefficient c = (g²/2)f² reaching the Lee–Wick pole-collision value c = M²/4
   (g f = M/√2) — directly from a strong oscillating scalar, by preheating of a doublet of mass ≈ M, or by a slow
   pump on an orbit with nonzero mean; the twin's pair never collides and Yang–Mills alone has nothing to collide.
   The pumped loss is a fraction, not a threshold or a rate (8 of 20 perturbed ghost runs within 60 T against 0 of
   10 controls). In chromo-natural inflation and its tail a massless doublet does nothing: the pump needs the axion
   drive within ~15 % of its own loss threshold, i.e. ϕ_end = O(1), which the benchmark does not give.

The healthy twin has the same anisotropic cancellation and an isotropic threshold where its own full propagator is
positive, 1/√2 < m_φ/M < 1. Against the chaotic and Bianchi IX backgrounds: the ghost pole does nothing here unless
m_φ > M, and then nothing at second order for the anisotropic field.

## Layout

Scripts in `code/`, each writing the log of the same name in `results/` (`probe_core.py` has none). All runs in a cloud session except `axion_probes_long.py` (run locally; its docstring), Python 3.13, scipy 1.18, sympy 1.14; float unless marked;
flat space, g = M = 1 unless stated. No figures.

| script (added) | run time; dependencies |
|---|---|
| `iso_floquet.py` (12 Oct) | seconds per energy; reads `../../ym/rung1_system.pkl` |
| `frw_identity.py` (12 Oct, sympy), `frw_ghost.py` (12 Oct) | seconds |
| `rung1_frw.py` (14 Oct, sympy) | 20 s; from `ym/hdym_rung1.py`; writes `../rung1_frw_system.pkl` |
| `frw_offdiag.py` (14 Oct) | 200 s; reads `rung1_frw_system.pkl` |
| `axion_derive.py` (14 Oct, sympy) | 3 s, all PASS (CHECK 2's first comparison line prints "False" by design: it records the ε⁰¹²³ sign choice; "ALL CHECKS PASSED") |
| `axion_aniso.py` (14 Oct) | 374 s; Part D withdrawn (History) |
| `axion_twins.py` (15 Oct) | 213 s; imports `axion_aniso` |
| `rung1_axion.py` (15 Oct, sympy) | 14 s, all PASS; from `ym/hdym_rung1.py`; writes `../rung1_axion_system.pkl` |
| `offdiag_axion.py` (15 Oct) | 23 s |
| `axion_dressed.py` (15 Oct) | 243 s (248 s also quoted); imports `axion_aniso`, `axion_twins` |
| `iso_lame_gap.py` (16 Oct) | ~2 min; imports `axion_frw_iso` |
| `axion_frw_iso.py` (16 Oct) | ~6 min |
| `probe_core.py` (16 Oct, numpy) | the one right-hand side with probes: (q, w) form, drives κ(t) and Θ(t) as slots, boson and fermion on or off, ghost or twin σ = ±1, `energy_parts`, `tangent_rhs`, adiabatic number; shared with `dilaton_su2/` — change it and rerun `probe_derive.py` |
| `probe_derive.py` (16 Oct, sympy) | 6 s, all PASS; checks `probe_core.py` against the Euler–Lagrange equations at random points for every switch, energy conserved; CHECK 5b added 16 Oct |
| `axion_probes.py` (16 Oct, revised) | 10 min (the first version's 7 min log was replaced); imports `probe_core` |
| `axion_probes_long.py` (16 Oct) | 10 min (first run 6 min); imports `axion_probes`, `probe_core`; local-run protocol (CHECK line, CONTROL-FAIL tags, DONE) |
| `axion_probes_spread.py` (16 Oct) | 17 min; imports `axion_probes` |
| `axion_probes_infl.py` (17 Oct) | 3.5 min; imports `axion_probes` |
| `axion_probes_fourier.py` (17 Oct) | 16 min; imports `axion_probes` |

Logs: `iso_floquet.txt`, `frw_identity.txt`, `frw_ghost.txt`, `rung1_frw.txt`, `frw_offdiag.txt`, `axion_derive.txt`,
`axion_aniso.txt`, `axion_twins.txt`, `rung1_axion.txt`, `offdiag_axion.txt`, `axion_dressed.txt`, `iso_lame_gap.txt`,
`axion_frw_iso.txt`, `probe_derive.txt`, `axion_probes.txt`, `axion_probes_long.txt`, `axion_probes_spread.txt`,
`axion_probes_infl.txt`, `axion_probes_fourier.txt`. Committed pickles in the folder root: `rung1_frw_system.pkl` (K,
Gauss matrix G, zero modes, L2; lambdify over x, y, z, px, py, pz, c_j, a, a1, a2) and `rung1_axion_system.pkl` (K, G,
zero modes, q⁗ of the diagonal background; lambdify over x, y, z, px, py, pz, x2, y2, z2, x3, y3, z3, c_j, kappa).

## Stage 1: flat space, no axion (12 Oct; homogeneous perturbations only)

On x = y = z, s = A cn(√2 A t | k² = ½), A = (2E/3)^{1/4}, period 4K(½)/(√2 A) (`iso_floquet.py`):

- **Diagonal ghost.** Hess V = 2s² × (all-ones matrix). The two anisotropic ghost modes are free (w″ = −w). The
  isotropic mode obeys w″ = −(1 + 6s²)w, in τ = √2 A t the Lamé equation with n = 2, k² = ½, h = 3 + 1/(2A²). Band
  edges h = 1.268, 1.5, 3, 4.5, 4.732, so for h ≥ 3 the only instability is h ∈ (4.5, 4.732), i.e. **E ∈ (1/8, 1/6)**
  exactly; the scan finds E ∈ [0.1263, 0.1661] on its grid, maximum rate 0.028.
- **Off-diagonal ghost (c_j = 1).** Floquet-stable at every E of the scan, 0.005–10; the three blocks agree to 1e-13.
  With c_j = 0 (diagnostic) unstable for E ≳ 1.5. Extended 14 Oct by the frozen-a scan of `frw_offdiag.py` (M → a is
  the same problem at E = (3/2)/a⁴, checked to 1e-14 at five spot values): stable at every E of 0.0024–384, 120
  log-spaced points; with c_j = 0 the one unstable interval is E > 1.51 (a < 0.998), max rate 0.21 per unit time at
  a = 1.
- **Yang–Mills perturbations, diagonal and off-diagonal.** Marginal (all multipliers 1, Jordan blocks): no
  exponential instability among homogeneous perturbations at linear order; Fujita et al.'s assumption holds
  marginally, not robustly. Inhomogeneous modes are not covered (on a closed universe Kumar, Lechtenfeld, Picanço
  Costa, arXiv:2102.08401, find bands of instability for spatial harmonics).

So one diagonal-ghost band and none for the off-diagonal ghost (compare the diagonal orbit's several tongues,
`paper/main.tex` sec:floquet, and the DR orbit's off-diagonal instability over about half of E = 0.05–1,
`islands/README.md`).

## Stage 2: expansion (12 and 14 Oct)

- **Derivation** (`frw_identity.py`, all True). On ds² = a²(dη² − dx²), ∇_μF^{μν} = a⁻⁴ × (flat divergence),
  √−g J·J = a⁻² × (flat J·J), √−g F·F = flat F·F: in conformal time the Lee–Wick Lagrangian is the flat one with
  M → M a(η). In the auxiliary form (−J²/2M² = stat_w [w·J + M²|w|²/2]) M appears only algebraically, so the linear
  diagonal ghost obeys w″ = −(M²a² + Hess V)w on the unchanged conformal-time background.
- **Isotropic diagonal ghost.** With comoving amplitude A the frozen-a Lamé parameter is h = 3 + M²a²/(2A²): the
  isotropic ghost crosses its single band once, while the physical energy falls through (1/6, 1/8), and is stable
  afterwards; the anisotropic ghost modes are never driven. Growth in the crossing (`frw_ghost.py`, radiation era
  a ∝ η): e-folds ≈ μ̄ × (time in the band), μ̄ = 0.041 per unit conformal time (A = 1), i.e. ≈ 1.6×10⁻³ M/H with H
  the Hubble rate at the band: 16.5 at H/M = 10⁻⁴, about 1 at 3×10⁻³, fractions of one for H/M ≳ 0.03. H/M at the
  crossing depends on whether the gauge field dominates the energy density (not decided here).
- **Off-diagonal blocks with M → M a(η)** (`rung1_frw.py`, `frw_offdiag.py`, 14 Oct). Fourth-order form with the
  Lee–Wick term / a(η)²: u⁗ coefficients −I/a²; a, a′, a″ enter the first-order matrix (a″ linearly, a′ up to a′²,
  each over a power of a); the linearised Gauss matrix picks up a′. Checks (all PASS): a = 1, a′ = a″ = 0 reproduces
  `ym/rung1_system.pkl` symbolically, K and G; the gauge zero modes u = θ × A solve the system and the Gauss law for
  arbitrary a(η). Radiation era a = η/η_s, comoving amplitude 1, block (u12, u21), ln σ_max of the transfer matrix
  (raw state, no WKB normalisation; two background phases); Gauss charge conserved to 1e-13. a from 0.3 to 3 (E_phys
  from 185 to 0.0185), c_j = 1: 7.5, 7.0, 7.1, 8.9, 9.8, 11.1, 12.0 for η_s = 1, 3, 10, 30, 100, 300, 1000 (H/M at
  a = 1 from 1 to 10⁻³; phase 0) — from η_s = 30 to 1000 that is +3.1 while ln(η span) grows by 3.5: growth ∝ (time)¹,
  the flat Jordan drift, not exponential (the O(1) base, 7.5, is the change of scale of the derivative components).
  c_j = 0, the control: 7.4, 6.9, 7.2, 10.4, 18.1, 42.7, 127.6, proportional to the conformal time in its frozen band
  a < 1 ((127.6 − 7)/700 ≈ 0.17 per unit η at η_s = 1000 against the frozen maximum 0.21) — the method sees an
  instability when one exists. a from 1.4 to 2.3 (the diagonal ghost's band window): c_j = 1 gives 2.6 … 8.6, c_j = 0
  gives 2.6 … 7.7, both ∝ (time)¹. So the expansion does not drive the off-diagonal ghost: no frozen band for c_j = 1
  at any a in [0.25, 5], and the non-adiabatic terms add nothing beyond the secular drift at crossing rates from
  H ~ M to H ~ 10⁻³ M.

## Stage 3: the axion drive

### Derivation (14 Oct, `axion_derive.py`, all PASS)

Fujita et al. appendix A: inflaton decay and stimulated emission are negligible in matter-dominated reheating
(A1–A5); the coupled isotropic system (A6–A7; x = m_φ t, h = H/m_φ, ϕ = φ/f, q = gQ/m_φ) barely moves q at their
parameters. Here the Lee–Wick term is built from F alone; φ(t) is a given function (its equation, their (12), unused).

- (1/4f) φ F·F̃ = ±(g/f) φ d(det A^a_i)/dt, so up to a total derivative L_axion = (g/f) φ̇ det A, sign fixed by matching
  their (13) at H = 0 (the ε⁰¹²³ choice; the oscillating stage does not depend on it). No A_0 in the term, so the Gauss
  law is unchanged; the diagonal ansatz stays a consistent truncation, L_axion = (g/f) φ̇ xyz.
- Auxiliary form: the time-dependent potential −κ(t) xyz, κ = (g/f) φ̇, added to V in q = p − w, nothing else; the
  Lee–Wick term keeps ∇V (not ∇V_tot). Verified against the fourth-order form.
- Units g = M = 1 (time 1/M, fields M/g): the one new dimensionless quantity is κ/(gM) = φ̇/(fM) = 2ξ H/M with their
  ξ = φ̇/(2fH); their (13) at H = 0 reads Q″ + 2Q³ = κQ².
- On FRW in conformal time the term is unchanged, (g/f)(dφ/dη) det A with comoving fields (the Chern–Simons density
  is conformal), so stage 2's M → M a(η) carries over (CHECK 7).
- **The drive sources the ghost.** On shell w = (q̈ + ∇V)/M² and the driven background has q̈ + ∇V = κ ∇(xyz) ≠ 0:
  the ghost-free manifold w = W = 0 is not invariant; the isotropic background is (Q, Q, Q; w₀, w₀, w₀),
  w₀ = κQ²/M² + O(M⁻⁴), its fourth-order equation printed in `results/axion_derive.txt` (their (13) as M → ∞).
- **Anisotropic sector** (the two permutation-odd directions, degenerate, w₀ drops out): ε⁗/M² + ε″ + κ(t)Q(t) ε = 0,
  a Pais–Uhlenbeck pair. Frozen coefficients ω² = (M²/2)(1 ± √(1 − 4κQ/M²)); Yang–Mills branch ω² ≈ κQ, ghost branch
  ω² ≈ M² − κQ; complex ω² for κQ > M²/4 (CHECK 6; in their variables κQ/M² = 2ξ m_Q (H/M)², ξ ≈ m_Q + 1/m_Q ≈ 5, so
  that would need H/M ≳ 1/√(8 ξ m_Q) ≈ 0.07, never reached at Fujita's parameters). Isotropic tangent as M → ∞:
  ε″ + (6g²Q² − 2κQ) ε = 0.

**Their variables** (their (17), (A8)): Q_ours = q m_φ/M; H/M = h m_φ/M; κ/(gM) = 2ξ H/M; m_Q = Q_ours/(H/M); the Hill
coefficient κQ/M² = q (dϕ/dx)(m_φ/M)²; with V(φ) quadratic, κ = −g m_φ ϕ_amp sin(m_φ t), κQ/M² = (m_φ/M)² q ϕ_amp cos(·).
Their Fig. 5 benchmark (used in `axion_aniso.py`): g = 10⁻¹², m_Q,end = 5, m_φ = 5×10¹⁰ GeV, f = 0.004 M_Pl, so
H_end = gM_Pl/m_Q² = 9.6×10⁴ GeV (their (33)), Q_end = M_Pl/m_Q (their (32)), gQ_end = 4.8×10⁵ GeV, q_end = 9.6×10⁻⁶,
h_end = 1.9×10⁻⁶, ξ_end = m_Q + 1/m_Q = 5.2; their dark-matter mass m set to 0. Scales m_φ : √(κQ) : gQ : H =
1 : 3×10⁻³ : 10⁻⁵ : 2×10⁻⁶: the sinusoidal drive with frozen Q is the physical fast problem, the slow Q cycle sweeping
its amplitude. The 14 Oct assumption ϕ_amp = 1 was replaced on 16 Oct by ϕ_end from their (A10)–(A11).

### Literature (14 Oct)

Nothing treats the anisotropic homogeneous perturbations while the axion oscillates. Inflation: Wolfson, Maleknejad,
Murata, Komatsu, Kobayashi, arXiv:2105.06259 (Bianchi I always isotropises; the Chern–Simons term depends only on
ψ = (det A)^{1/3}, the nonlinear form of CHECK 1); Maleknejad, Erfani, arXiv:1311.3361; Murata, Fujita, Kobayashi,
arXiv:2211.09489. After inflation: Khan, Pirzada, Mustafa, arXiv:2604.07044 (k = 0 fluctuations of the massive phase,
1 ⊕ 3 ⊕ 5 with frequencies m² + 6g²Q², m² + 2g²Q², m²; the gapless quintet at m = 0 is the stage-1 result; no Floquet
analysis); Daniel, Kamali, arXiv:2603.12339 (CHECK 7 with κ = (gλ₁/f) dχ/dη; tensor modes only); Adshead, Giblin,
Weiner, arXiv:1708.02944 (e^{φ/M} Tr F², not φFF̃). Kuwabara, Yumibayashi, Harada, arXiv:1503.03657: the
time-dependent Pais–Uhlenbeck oscillator does not decompose, so the 4-d system is integrated directly.

### The frozen-Q periodic problem (14–15 Oct; `axion_aniso.py` A–C, `axion_twins.py` A–C)

Statements about the periodic problem ε⁗/M² + ε″ + (c₀ + c₁ cos ωt)ε = 0 (pole terms only), which stand.

- *Method.* Monodromy over one drive period; multipliers in reciprocal pairs {λ, 1/λ}, read off x = λ + 1/λ from the
  characteristic polynomial (the free branch sits at a Jordan block, x = 2, where eigenvalues are √-sensitive to
  roundoff and x is not; A2: det = 1 to 1e-11, pairs to 1e-9); ω = 1/n exactly avoided. A3: as ω → 0 at fixed Mathieu
  (a, q) the 4-d rate agrees with the 2-d Hill equation (difference ∝ ω²) and with scipy's a_n, b_n at six points.
- *The frozen-Q map* (B). Below ω = 1: stable except the ghost's 2/n tongues, O(c₁ⁿ), and, for c₁/ω² above the
  Kapitza threshold, the frozen-tachyon regime: threshold 0.4521 (ω = 0.1), 0.4539 (ω = 0.03) against Mathieu's
  b₁(q) = 0 at q = 0.9080, c₁/ω² = 0.4540. Near ω = 1 the combination resonance (ghost at 1, free branch at 0): a Krein
  bubble below 1 and a peak above it, peak rate ∝ c₁^0.65, bubble width ∝ c₁^0.69 (c₁ = 0.001–0.03). Above ω = 1 the
  whole range is unstable by the ghost-pole law, measured rate / c₁(ω√(2(ω²−1)))⁻¹ = 1 to 9×10⁻⁵ for ω ≥ 1.5,
  c₁ ≤ 0.01 (B3); below ω = 1 the same term is a frequency shift, x − 2 = −(νT)² with ν² = c₁²/(2ω²(1−ω²)), to three
  digits. On the self-consistent background this pole law is cancelled by the jitter term.
- *The tongue at r = ½* (C3, physical c₁ = 9.6×10⁻⁶ ω²): ghost-pair rate c₁/4 at the centre (2.5 ϕ per Hubble time),
  falling to the band value 2.041 per Hubble time at ±0.75 qϕ, relative width qϕ ~ 10⁻⁵. Near r = 1 at the physical
  c₁: unstable for ω ∈ [0.9996, 1.0030], peak 2.36×10⁻⁴ per unit M⁻¹ at ω = 1.0002, i.e. 123 per Hubble time,
  indicative only (the law diverges there).
- *The twin suite* (`axion_twins.py`). L = ½s′² + Σσ_i(½h_i′² − ½M_i²h_i²) − ½c(s + Σh_i)²: ghost (reference),
  healthy twin (σ = +1), ghost + healthy at M_h = 2M (the sixth-order operator with poles 0, M, 2M and signs +, −, +),
  M_h = M (tuned), M_h = M/2 (reversed), each heavy mode alone (Mathieu). Fast responses f = c₁s cos(ωt)/ω²,
  h_i = −σ_i c₁ s cos(ωt)/(M_i² − ω²), D(w) = Σσ_i/(M_i² − w); the pole term is tachyonic iff D(m_φ²) > 0 — ghost:
  m_φ > M; healthy twin: m_φ < M; M_h > M: M < m_φ < M_h; M_h < M: m_φ outside [M_h, M]; M_h = M and no heavy mode:
  never. Checks (A, all PASS): the ghost case reproduces `axion_aniso.pairs` to 2×10⁻¹³; c₁ = 0 gives the frozen
  frequencies; the characteristic polynomial is palindromic to 10⁻¹¹; M_h = M gives rate ≤ 2×10⁻⁹ off resonance; a
  heavy mode alone of either sign has the Mathieu rate c₁/(4M) at ω = 2M and nothing (10⁻¹³) at 1.7M, 0.9M. Frozen-Q
  map (B, c₁ = 0.01): healthy twin unstable band ω ∈ [0.7075, 1.045] (the window 1/√2 < ω < 1 plus the resonance) and
  the ghost-free 2M tongue, stable above; M_h = 2: one band [0.955, 2.035], the 4M tongue; M_h = ½: bands
  [0.3675, 0.56] and [0.9525, 1.365]; M_h = M: isolated grid points at ω = 1/n and 2 only. Per Hubble time at the
  benchmark, ϕ = 1 (C): healthy twin 7.7 (r = 1.1), 3.2 (1.5), 1.5 (2.5), 0.9 (4), 0 for r < 1; two poles 8.8
  (r = 0.55), 6.1 (0.7), 8.4 (0.9), 0 outside ½ < r < 1. The sign rule is correct for the pole term alone.

### The self-consistent background and the cancellation (15 Oct; `rung1_axion.py`, `offdiag_axion.py`, `axion_dressed.py`)

- *The signal.* The direct runs of 14–15 Oct used the M → ∞ background Q″ = −2Q³ + κQ² (the massless pole's response
  only; for m_φ ~ M the ghost pole's share is O(1)). On it the gauge zero modes, exact only on a true solution, fail
  at O(κ): in `offdiag_axion.py` B the (u12, u21) block at ω = 0.8 grew 1.69 per cycle in its antisymmetric half,
  |ΦZ₀| reaching 4×10⁷ by cycle 12, Gauss drift 2×10⁵ (first version; its log was overwritten, no committed log).
- *The dressing.* The fourth-order isotropic equation gives J̃ = Q″ + 2Q³ − κQ² with J̃″ + (1 + 6Q²)J̃ = −(κQ²)″ + …,
  so the fast jitter of Q is δQ = κX²Π(ω²), Π(w) = 1/(w² − w) = −1/w + 1/(w − 1) (M = g = 1; the M → ∞ background has
  −1/w only); ⟨κδQ⟩ = ⟨κ²⟩X²Π. Measured: the background's period at ω = 0.8 is 198 (exact) against 305 (M → ∞) (no committed log), the
  drive on Q being κQ²/(1 − ω²).
- *Anisotropic sector: cancellation* (`axion_dressed.py` A1). With the DC offset c₀ = (c₁²/2)Π the massless pair's
  rate is 0 to roundoff for the ghost at ω = 3, 1.5, 2.5 ("5" listed until 28 Oct 2026: no such row in `axion_dressed.txt`); the residuals wherever |Π| is large scale as c₁² (4.6 c₁²
  for the ghost at ω = 0.8, 4.9 c₁² at 0.6 ("4.6 … at 0.8 and 0.6" until 28 Oct 2026); 0.66 c₁² for the twin at 1.5; 0.13 c₁² for the two-pole operator at 2.5), fourth
  order, 10⁻⁵ of the withdrawn law at the physical c₁. As ω → 0 (M → ∞) the pole term is the Kapitza stabilisation
  −c₁²/(2ω²) = −ϕ²X²/2 and cancels the jitter term exactly: marginal, as in stage 1 (the linear form of Wolfson et
  al.'s isotropisation).
- *Isotropic sector: no cancellation* (A2). Coefficients −2κQ against +κQ: with c₁_iso = 2c₁, c₀_iso = −(c₁_iso²/4)Π
  the model s″ = (3/4)c₁²Π s holds to 0.3 % (ghost ω = 3, 1.5, 2.5), 4 % (twin ω = 0.8), 0.5 % (two poles ω = 1.5;
  "0.5 %" for all three in the 15 Oct summary); zero wherever Π < 0.
- *The off-diagonal sector* (`rung1_axion.py`, all PASS; `offdiag_axion.py`). L₂ − L₂(κ=0) =
  κ(−x u₂₃u₃₂ − y u₁₃u₃₁ − z u₁₂u₂₁) and nothing else; u⁗ coefficients −I; κ = 0 reproduces `ym/rung1_system.pkl`
  symbolically; on the exact driven diagonal background (q, q′, q″, q‴ free, q⁗ from the diagonal equations) the
  gauge zero modes u = θ × A are exact for any κ(t). Quintet identity (A, 20 random background points, residuals
  ≤ 2×10⁻¹⁵): the symmetric subspace u_ab = u_ba of each 8-dim block is invariant under K, which restricted to it is
  the companion matrix of ε⁗ = −ε″ − κQε; the three blocks agree; the gauge zero mode is in ker G. The antisymmetric
  half is a Pais–Uhlenbeck pair whose light branch is the gauge zero mode itself (d ∝ Q(t), frequency √2 Q) with the
  one Gauss condition (the Yang–Mills mass 2Q² of the gauge triplet is what Gauss removes). Direct run (B; exact
  background, dressed data, ϕ = 1, m_φ/gQ = 100, 6 Q periods, ω = 0.8, 2, 3 and block 13 at 2): ln σ_max of the 8-dim
  block, its symmetric restriction, the 4-dim anisotropic sector and the block on the Gauss surface modulo the gauge
  direction agree to the printed precision in every row (e.g. 5.87, 6.56, 6.97, 7.25, 7.48, 7.66 at ω = 2: the O(1)
  transient, then ln t); |ΦZ₀| stays in [0.03, 1.5]; Gauss drift ≤ 1.1×10⁻¹⁴.
- *Threshold scan* (`axion_dressed.py` B; exact background, dressed data, m_φ/gQ = 100, 4 Q periods, period from
  smoothed zero crossings; table in the log). Isotropic threshold in ω²ϕ²Π(ω²): between 0.33 and 1.33 at ω = 2 (1.33
  slowly unstable), above 1.125 at ω = 3, above 1.8 at ω = 1.5, none for ω < 1 — "of order 1–2, not sharp"; the
  averaged estimate (3ω²ϕ²X²Π against 6⟨Q²⟩ = 2.74X²) put it at 2 (the ω = 2 rows sit on the ghost's tongue, below).
  Stable rows: the free drift 0.69, 0.41, 0.29 after the first-cycle transient in both sectors.
  Unstable rows: 2.5–6 per cycle in both at ω = 2, ϕ = 2–2.2 (the isotropic solution runs away in its own ghost
  content, the anisotropic sector carried along), |Q| reaching 1.2–2.6 Q_amp within four cycles or blowing up within a cycle (BLOW-UP rows,
  stopped at |Q| = 10 Q_amp); at (2, 2) it runs away by the twelfth cycle (no committed log; `axion_dressed.txt` B runs 4 cycles). On the M → ∞ model the background's period
  agrees with the Kapitza-averaged quartic Q″ = −(2 + ϕ²)Q³ at X = Q_amp(1 + 2ϕ²/(2 + ϕ²))^{1/4} to 0.1 % at ratio 200, 1.1 % at ratio 50 (`axion_aniso.txt` D, T_Q meas vs T_Q Kap: 60.57/61.24) (Fujita's
  "the coupling barely moves q": the shift is in the effective coupling).

### The isotropic threshold, exact (16 Oct; `iso_lame_gap.py`)

- *The Lamé-gap argument.* On the undriven orbit the isotropic tangent ε″ + 6Q²ε = 0 is Lamé with n = 2, k² = ½,
  h = 3: the lower edge of the band [3, 4.5], the Jordan tangent of the orbit family. A constant downward shift of h
  lands in the gap (1.5, 3), predicting instability for any positive tachyonic coefficient, against the scan.
  One-period monodromy: (a) constant shift −δ⟨Q²⟩: unstable for every δ > 0 (rate 0.35/T already at
  δ = 0.03); (b) a Q²-shaped shift on the unshifted orbit, (6 − δ)Q²: also unstable; (c) the Q²-shaped shift on the
  consistently shifted orbit, ε″ + 3λQ²ε = 0 on Q″ = −λQ³, λ = 2 − δ/3: multiplier exactly 1 for every λ > 0. The
  dressed isotropic equation (A2 with c₁_iso = 2κQ) is ε″ + (6 − 3κ²Π)Q²ε = 0 = ε″ + 3λ_eff Q² ε,
  λ_eff = 2 − κ_amp²Π = 2 − ω²ϕ²Π(ω²), on the averaged background Q″ = −λ_eff Q³: case (c) — the tachyonic term
  rescales the Q² coefficient, moving h and n(n+1)k² together along the band edge (the linearisation of an autonomous
  averaged Hamiltonian system). **Threshold λ_eff = 0, ω²ϕ²Π = 2 = ϕ² r²/(1 − r²)**, sharp; above it the averaged
  quartic potential is inverted and the loss takes 0.4–0.9 of a Kapitza period.
- *Against the exact background* ((Q, w) form, Q″ = w − 2Q³, w″ + (1 + 6Q²)w = κQ², checked against the pickle's q⁗
  to 5×10⁻¹⁶ and three rows of the 15 Oct scan in `axion_frw_iso.txt` Part 0; m_φ/gQ = 100, 8 periods): at ω = 1.8, 2.4, 3.0 steady at
  ω²ϕ²Π = 1.0, 1.33, 1.7 (amplitude 1.1 → 1.6 Q_amp, the (2/λ_eff)^{1/4} of the averaged orbit) and lost at 2.3
  (0.4–0.55 T); at exactly 2.0 lost at ω = 1.8, steady at 2.2–2.8 Q_amp at ω = 2.4, 3.0 (the next order holds it). At
  ω = 2.0 — the ghost pair's first tongue, which `axion_dressed.py` A skipped and B did not — lost already at 1.33
  (2.6 T) and 1.7. The m_φ = M resonance (Π singular) is the common limit of both thresholds, both → 0 there, and is
  not run separately.

### The FRW run and the physical window (16 Oct; `axion_frw_iso.py`)

Conformal time, comoving Q, M → M a(η), drive (g/f)(dφ/dη)Q³: Q″ = a²w − 2Q³, w″ + (a² + 6Q²)w = a κ_phys Q².
Reheating as in Fujita et al.: the axion dominates in its quadratic potential, H = H_end a^{−3/2}, ϕ(a) = ϕ_end a^{−3/2},
ϕ_end from (A10)–(A11); the gauge field's energy is not fed back into H. Parameters m_φ/H_end, r, m_Q,end = 5,
f = 0.004 M_Pl; a from 1 to 20; initial data Q = Q_end, physical Q̇ = 0, the dressed particular solution for w.

- Part A (25 points): for m_φ/H_end ≥ 300 the averaged prediction holds row by row — lost iff λ_eff(a = 1) < 0
  (r = 0.8, 0.9 at 300; r = 0.9 at 450; none at 700) — except the marginal row (300, r = 0.7), λ_eff(1) = −0.002, steady:
  ∀ (m, r) ∈ G_A, m ≥ 300, (m, r) ≠ (300, 0.7): LOST(m, r) ⇔ lam_eff(1) < 0 [grid]
  (G_A: the 25 rows of `results/axion_frw_iso.txt` Part A; columns "lam_eff(1)", "verdict"; the exception was not stated until 28 Oct 2026); the loss complete by a = 1.1–1.25, before the dilution crossing
  a_× = (ϕ_end²/(2(ω² − 1)))^{1/3} = 1.2–1.6: the expansion plays no part. Below threshold the orbit persists with the
  larger amplitude of small λ_eff (1.9–2.3 Q_end at λ_eff = 0.35–0.37, efold-0 max|Q|/Q_end 1.92 and 2.32; "2.0–2.3" until 28 Oct 2026) and relaxes as the drive dilutes.
- At m_φ/H_end = 100 every point is lost, including r = 1.2 (M > m_φ), and so is the no-ghost control (Part B): the
  slow-drive instability of the isotropic Hill equation, pure Yang–Mills' own, at s_YM = 2ϕ_end m_Q H_end/m_φ = 0.43.
  At 200, r = 1.2 is lost with the ghost (a = 2.4) and steady without (s_YM = 0.11): near the resonance the ghost
  lowers that threshold. Part B′ (flat space, m_φ/gQ = 40, 16 periods): pure-YM threshold between s_YM = 0.2 and 0.3 at
  ω = 0.5 and 0.833, the ghost one the same at ω = 0.5 and below s = 0.2 at ω = 0.833 (lost after 2.2 T); the dressing
  1/|1 − ω²| of the drive on Q (`axion_dressed.py` C) accounts for part of that, the combination resonance near ω = 1
  (not characterised) for the rest.
- Part C: both thresholds are quadratic in H_end/m_φ: ghost-pole m_φ/H_end < √(3/2)(M_Pl/f) r/√(1 − r²) (300 at
  r = 0.7, 630 at r = 0.9 for f = 0.004 M_Pl), slow-drive m_φ/H_end < √(2√3 m_Q (M_Pl/f)/0.25) ≈ 130 at the same f. At
  the benchmark (m_φ/H_end = 5.2×10⁵, ϕ_end = 8.3×10⁻⁴; quoted as 5×10⁵ and 8×10⁻⁴ in the 15–16 Oct summary)
  λ_eff(1) = 2 − 6.9×10⁻⁷ r²/(1 − r²) (ϕ_end² = (8.31×10⁻⁴)²; stated as 3.5×10⁻⁷ until 28 Oct 2026, which is the bound on 1 − r²): nothing, for any M, short of 1 − r² < 3.5×10⁻⁷. Where the loss happens it
  completes within Δa ≈ 0.1–0.4 of the end of inflation.
- *The coupled axion, sized (not run).* The φ equation from (g/f)φ̇Q³ is φ̈ + m_φ²φ = −(g/f) d(Q³)/dt; relative to
  m_φ²φ, with g²Q⁴ ~ ρ_Q ~ M_Pl²H² (their ε_B = 1 at the end of inflation) and ϕ_end from (A11), the back-reaction is
  of order √6 (M_Pl/f)(H_end/m_φ) = √2 ϕ_end: 1.2×10⁻³ at the benchmark, O(1) wherever either threshold is reachable.
  Their A6–A7 with the ghost (the (Q, w) form plus the φ equation, three second-order equations) decides the window.

## Probes: real particles coupled to the ghost (16–17 Oct)

Anything with SU(2) charge couples to the full gauge field q = p − w, ghost included, through its covariant
derivative (the healthy twin stays as a sign control). Derived in `probe_derive.py` (all PASS), implemented in
`probe_core.py`:

- **Boson:** a fundamental scalar doublet (the Higgs is one). With A₀ = 0 its Lagrangian depends on the gauge field
  only through Σ(A^a_i)², so its current is diagonal and the diagonal truncation survives; the Gauss law (ρ^a = 0,
  including the overall phase) leaves Φ = f(t)Φ₀, f real: L_f = ḟ² − (m_s² + (g²/4)|q|²)f², force −(g²/2) q f² on the
  gauge field. The free anisotropic directions (and the off-diagonal components) get a direct effective mass
  (g²/2)f², a second Hill drive at twice the scalar's frequency.
- **Fermion:** a classical Dirac doublet (8 complex components). The diagonal SO(3) (spin + colour) has a
  two-dimensional singlet sector invariant under H_D on every diagonal background: iζ′ = h(Q)ζ,
  h = [[m_D, 3gQ/2],[3gQ/2, −m_D]], Q = (x + y + z)/3, zero colour charge, current J^a_i = j δ^a_i, j = −½ ζ†σ_xζ,
  force g j on each gauge direction. It sees only the trace (no coupling to the quintet at linear order). ζ†ζ is
  conserved, so its energy is bounded by construction; the observable is n₋, its weight on the negative-energy level,
  meaningful only for m_p ≳ gQ_amp (at m_p = 0 the levels cross as Q passes through zero, n₋ ≈ ½ for every heavy
  sector).
- **Adjoint triplet** (record only): consistent only along an axis, then anisotropic by construction (mass
  g²(x² + y²)f² for n = ẑ); not in the core. **Light:** only through a carrier charged under both groups (a doublet
  with hypercharge); not done.
- **Twin:** σ = −1 flips the heavy sector to a healthy oscillator in the (p, w) variables, L_twin = L′ + |ẇ|² − M²|w|²
  (Θ = 1 only: with the dilaton on, a sign flip of the w sector makes the kinetic matrix degenerate at some Θ; a sign
  flip of the cross term alone leaves it indefinite, hence the added |ẇ|²). Its heavy coefficient along the isotropic
  direction is M² − Hess V (M² + Hess V for the ghost), negative for 6Q² > M², Q > 0.41 M: **the twin is a valid control
  only where Q ≪ M**, intrinsic to any Lagrangian twin that keeps the potential and matter coupled to q = p − w; at
  Q_amp = M the ghost-vs-M → ∞ comparison is the control (`probe_core.py` docstring).

Run setup: the exact (Q, w) background, probes at a fraction ε_p of its energy; the isotropic truncation of
`probe_core.rhs` checked against the full one at random points (4×10⁻¹⁶); hand-coded tangents checked against the
core's finite-difference linearisation at 20 random points and 20 along the orbit (errors 10⁻¹⁰). Mass labels
(corrected 17 Oct): Parts A and C of `axion_probes.py` run at gQ_amp = M/100 (Q_amp = 0.02 at ω = 2), Part B at
gQ_amp = M; the physical axis is m_p/gQ_amp; Part C's rows "m_p/gQ_amp = 50 (= M/2)" and "100 (= M)" are m_D = M and
2M, Part A's decoupled control "m_p = M (= 100 gQ_amp)" is m_p = ωM (2.4M, 3M).

### Results (`axion_probes.py`, revised 16 Oct)

- *A. The isotropic threshold with probe columns* (none / ghost + boson / ghost + fermion / ghost + both / twin +
  both; m_p/gQ_amp = 0, 1, 3, and the decoupled control; 8 cycles; m_φ/gQ = 100). Thresholds unchanged in every
  column: steady at ω²ϕ²Π = 1.0, 1.7, 2.0, lost at 2.3 (0.45 T at ω = 2.4, 0.55 T at ω = 3.0). The boson at m_p = 0
  trades 0.4–5× its energy with the ghost background, at gQ_amp 0.3–1.5×, at 3 gQ_amp 0.98–1.00 at ω²ϕ²Π = 1.0, 1.7 and 1.00–1.15 at 2.0 ("0.98–1.00" for all until 28 Oct 2026): **the probe decouples
  at m_p ≈ 3 gQ_amp**; at m_p = ωM every column equals the no-probe row to the printed digits (coupling ~10⁻⁴ of m_p²,
  force ~10⁻⁵ of 2Q³). At m_p = 0 with the ghost the background's amplitude swings between 0.3 and 2.3 Q_amp cycle to
  cycle; with the twin or at M → ∞ it stays at 1.33–1.6 Q_amp while the probe's energy swings within about ±1 of its
  start. The fermion's n₋ at m_p = gQ_amp is 0.1–0.6 with the ghost, 0.4–0.5 with the twin, ≈ 0.5 at M → ∞; at
  3 gQ_amp 0.001 (ghost), 0.07 (twin). The ghost sector's cycle-averaged energy is O(1–30) E_bg and negative at
  ω²ϕ²Π = 1.7 and O(10²) at 2.0, the twin's the same magnitudes positive, M → ∞ zero (`probe_core.energy_parts`; the
  pieces are not separately conserved with a drive on).
- *B. The anisotropic sector with the scalar's (g²/2)f²* (linear tangent on the undriven orbit, Q_amp = 1, no axion;
  the scalar a test field at ε_p, 30 periods, renormalised each period; ghost / twin / M → ∞). At ε_p ≤ 0.1 nothing
  grows in any column (rates ≤ 0.12 per period, the free drift; "≤ 0.1" until 28 Oct 2026: twin, m_s = 0, ε_p = 0.1 reads 0.1180). At ε_p = 0.5 **the scalar drives the ghost pair and
  not the twin**: per period ghost / twin / M → ∞ = 1.13 / −0.00 / −0.01 (m_s = 0), 0.68 / −0.01 / −0.01 (0.3),
  0.11 / 0.11 / 0.06 (0.5), 0.23 / 0.00 / −0.00 (0.7), nothing at m_s ≥ 0.8 away from the preheating rows m_s = 1.0,
  1.1, which are invalid as test-field rows: the scalar is itself pumped by the background (its own growth 0.36 and
  0.31 per period, its parametric resonance with Q² — the gauge preheating of a doublet), so f² grows without bound
  and all columns inherit it; L2 is the back-reacting case.
- *C. On the ghost's tongue* (ω = 2M, ω²ϕ²Π = 1.0). The boson's energy is untouched at m_p ≥ 3 gQ_amp with either
  heavy mode. The fermion's n₋ at m_D = M (level spacing Ω = 2m_D = 2M = the drive frequency): with the twin
  0.013 → 0.62 over eight cycles, with the ghost 0.000 → 0.015 (`axion_probes.txt` C, cycle 8; "0.009", the cycle-6 value, until 28 Oct 2026); at m_D = 2M nothing with either. Taken apart 17 Oct
  (`axion_probes_fourier.py`): the fermion couples to Q(t) alone, and the twin's dressed orbit has **14× the Fourier
  amplitude at Ω** (0.0302 vs 0.00215 Q_amp; 0.040 vs 0.0030 at the drive 2M; the same 0.0006 at M; dressed amplitudes
  1.56 vs 1.14 Q_amp). The fermion as a test field on each recorded waveform reproduces the coupled rows (0.50 vs
  0.017); a pure tone A cos Ωt at the twin's amplitude gives 0.50, at the ghost's 0.007; a Rabi estimate agrees
  ((3/2)g·A = 9×10⁻⁴ over 2100 time units, angle 1.9 rad, n₋ = sin² ≈ 0.66 against the measured 0.64; hand estimate, no committed log); the ghost's
  component sits at 2.000 (no committed log) against Ω = 2.0009, detuning-limited (scaled up 14× it gives only 0.035). The heavy pole's
  sign enters only through the orbit's dressed response: pair production by a strong oscillating field, not a ghost
  instability.

### The mechanism: pole collision at c = M²/4 (16 Oct; `axion_probes_long.py`, `axion_probes_spread.py`)

The ghost tangent of the anisotropic sector with a Hill coefficient c is ε⁗ + M²ε″ + M²c ε = 0, frequencies
ω² = [M² ± √(M⁴ − 4M²c)]/2, which leave the real axis for **c > M²/4** (the light and ghost poles collide at
ω² = M²/2; CHECK 6's criterion for κQ). The twin's pair, ω⁴ − (M² + 2c)ω² + M²c, has discriminant M⁴ + 4c² and never
collides; M → ∞ is second order with nothing to collide. With c = (g²/2)f² the collision is at **g f = M/√2**.

- *Part B matches it row by row:* at ε_p = 0.5 the scalar's c_max = 0.50, 0.45, 0.375, 0.30 (m_s = 0, 0.3, 0.5, 0.7; not in `axion_probes.txt`, no committed log)
  all exceed ¼ and the ghost column grows; at ε_p ≤ 0.1, c_max ≤ 0.1 and nothing grows. At c = 0.5, Im ω ≈ 0.32 above ¼
  for about half of each period: ≈ 0.8 per period against the measured 1.13.
- *L2, m_s ≈ M with back-reaction* (undriven orbit, Q_amp = 1, T = 5.244, scalar only, 100 periods; no twin rows,
  Q_amp = M). Controls (no scalar) steady, tangent drift 0.05 → 0.01 per period, ghost = M → ∞ exactly. At
  m_s = 0.8, 0.9, 1.1, 1.2 the ghost and M → ∞ columns agree: the scalar's energy within 0.15–16× of its start, the
  background within 3 %, the anisotropic tangent marginal, c_max ≤ 0.28. At m_s = 1.0 they do not: M → ∞ preheats the
  scalar to 20–70× its start (ε_p = 0.01), c_max reaching 0.52 harmlessly, the background giving up 5–30 % of its
  amplitude; **with the ghost (c_max 0.45, 0.55) the anisotropic tangent grows at 0.8 then 6.4 per period and the run
  ends between periods 30 and 60** (|Q| > 10 or non-finite); at ε_p = 0.1 the background's amplitude rises to
  3.3 Q_amp and stays (M → ∞: 0.7–0.97). One point in m_s; the window's width is not measured.
- *L1, m_p = 0 over 60 Kapitza periods* (both probes at 10 % of the background's energy; m_φ/gQ = 100). Controls
  steady: M → ∞ with or without the probes 1.32–1.37 Q_amp throughout, probe energy within −0.25 … 1.4 of its start.
  With the ghost at ω²ϕ²Π = 1.7 the isotropic solution is lost at 32.12 T (ω = 2.4) and 47.01 T (ω = 3.0), the
  amplitude wandering between 0.02 and 1.5 Q_amp first and the boson's energy reaching 2000× its start (18539× in the
  rerun's diagnostics); at 2.0 it survives the 60 cycles with the same wandering (0.1–3.2 Q_amp). The twin is steady in all four rows:
  ∀ (ω, ω²ϕ²Π) ∈ {2.4, 3.0} × {1.7, 2.0}: lost_twin = "-" ∧ |Q|max/Q_amp ∈ [1.62, 1.75] [grid]
  (`results/axion_probes_long.txt` L1, rows "twin both"; 1.65–1.66, 1.62–1.63 at 1.7, 1.66–1.70, 1.69–1.75 at 2.0).
  (Until 28 Oct 2026: "steady at 1.7 (1.40–1.49 Q_amp) and lost at 2.0 (ω = 2.4, 50 T; ω = 3.0 survives, amplitude rising
  to 2.0), the no-probe threshold itself" — not in the committed log; probably from the replaced first run.) **The loss coincides with the crossing:** c climbs slowly (0.05, 0.27 over the first two 20-cycle
  blocks) and first exceeds ¼ at cycle 32, 27, 47, 29 in the four lost rows, 6g²Q² (≤ 0.012 before) one cycle later
  (33, 28, 48, 29), in the loss cycle itself. With the twin or at M → ∞ c stays at 0. One-probe rows: boson only,
  c ≤ 0.04, the background wanders (0.1–1.4 Q_amp) but is not lost in 60 cycles; fermion only, flat (1.50 Q_amp).
- *The fermion is a constant force.* At m_D = 0, h = cgQσ_x for all t, so a state prepared in the positive eigenstate
  of h(Q_amp) stays a σ_x eigenstate and its current is constant (j = −6.912×10⁻⁷ at nine sampled times, to ten
  digits; no committed log): a fixed force g j₀ = −n/2 on each gauge direction, giving Q a nonzero mean. The `boson+F`
  column replaces the spinor by that force and reproduces the pump: lost at 27.08 T (ω = 2.4) and 28.52 T (ω = 3.0),
  the boson's energy reaching 15905× and 1032×, the first three cycles agreeing with the two-probe row digit for
  digit (no committed log). **A massless boson on an orbit with nonzero mean pumps the ghost background until the collision.** At m_D ≠ 0
  the fermion evolves (Part A's n₋ moves between 0.1 and 0.6 at gQ_amp) and is a participant.
- *The loss is realisation-dependent* (`axion_probes_spread.py`). `both` and `boson+F` integrate the same ODE and lost
  at 32.1 vs 27.1 T and 47.0 vs 28.5 T: roundoff amplified by a pumped, ultimately unstable system. On the `boson+F`
  row at ω²ϕ²Π = 1.7 with initial data perturbed by 10⁻⁶ (five seeds) at two tolerances: with the ghost **8 of 20 runs
  lose the background within 60 T** (ω = 2.4: 5 of 10, at 12–46 T; ω = 3.0: 3 of 10, at 30–48 T, all three at the
  tighter tolerance); c crossed ¼ in every lost run, within one cycle of the loss in six of the eight, 19 cycles
  before it in one. Perturbed controls: twin 0 of 5 and M → ∞ 0 of 5 lost, c ≤ 10⁻⁴, the background flat (1.66 and
  1.33 Q_amp). Quoted as 8/20 against 0/10 with the 60-cycle horizon.

### The inflationary run (17 Oct; `axion_probes_infl.py`)

"An orbit with nonzero mean" is the inflationary regime of this model (Adshead–Wyman arXiv:1203.2264 eqs. 1–6; their
−(λ/8f)𝒳FF̃ with F̃ = εF is this folder's (1/4f)φFF̃ with F̃ = εF/2 up to the ε⁰¹²³ sign of CHECK 2): ψ is held on the
attractor ψ_min = (μ⁴ sin(𝒳/f)/3g̃λH)^{1/3} by the axion's drive with φ̇ of fixed sign (constant κ gives
Q̈ = −2Q³ + κQ² the equilibrium Q₀ = κ/2), no fermion needed. Setup: constant κ, the isotropic fixed point Q₀
(w₀ = 2Q₀³, κ = 2Q₀(1 + 6Q₀²); checked to 10⁻¹⁷ for the core, the twin's and the M → ∞ systems with their own κ), the
massless doublet at ε_p = 0.01, 0.1, 0.5 of E_bg = (3/2)Q₀⁴, Q₀/M = 0.024, 0.1, 0.3, 60 periods of the isotropic mode;
then the tail: the axion starting to oscillate at m_φ = 0.8, 2.4, 3.0 M with the amplitude inflation left it
(κ₀ = 2Q₀). Columns ghost / twin / M → ∞, each with its no-scalar control.

- *Static inflation:* **nothing** with the ghost or at M → ∞: At Q₀ = 0.024 the ghost rows equal the M → ∞ rows to the third
  digit; |Q − Q₀|/Q₀ ≤ 0.56 (the scalar's static push at ε_p = 0.5), the scalar's energy within 0.5–1.0 of its start,
  c ≤ 0.054 against ¼, E_heavy at its fixed-point value (−0.002, −0.04, −0.36 E_bg at the three Q₀), no row lost in
  60 periods (∀ rows of `results/axion_probes_infl.txt` A with heavy ∈ {ghost, inf}: lost at = "-" [grid]); five
  10⁻⁶-perturbed seeds identical. The twin agrees at Q₀ = 0.024 and at 0.1 for ε_p ≤ 0.1; at Q₀ = 0.1, ε_p = 0.5 it is lost
  at 1.66 T (not stated until 28 Oct 2026; "nothing, in every column" until then), and its fixed point is unstable
  at 0.3 — the twin limit, not a result.
- *The tail:* **nothing**. Q oscillates at amplitude ≈ Q₀ in every column; the scalar's energy swings between 0 and
  1.1–1.6× its start in all three; c ≤ 0.001; no loss in 50 periods. The regime reached has ω²ϕ²Π = 4Q₀²Π ~ 10⁻³,
  three orders below the no-probe threshold.
- *Reading.* L1's pump needed a nonzero mean and a drive within ~15 % of the no-probe loss threshold (ω²ϕ²Π = 1.7,
  ϕ = O(1) in units of f); inflation proper supplies the mean and not the drive, its tail neither. So the pump needs
  ϕ_end = O(1), which the benchmark (ϕ_end ≈ 8×10⁻⁴) does not give. For the paper: the mechanism (a probe-induced
  coefficient reaching the collision, g f = M/√2 — not a tongue, not a resonance; the twin cannot do it; Yang–Mills
  alone cannot), scoped to that window; outside it real probes are spectators. For the dilaton's "threshold or pump"
  question: pump, then collision.

## Open

- The coupled axion: Fujita et al.'s A6–A7 with the ghost (the (Q, w) form plus the φ equation), required in the
  window ϕ_end = O(1) where the back-reaction √2 ϕ_end is O(1). The one open piece of the axion stage.
- The near-resonance (ω → 1) ghost lowering of the slow-drive threshold (the combination resonance), not characterised.
- Probes: the m_s window of the L2 preheating (one point measured); the loss fraction as a function of starting
  energy and ω²ϕ²Π; whether another source of a nonzero mean Q (the axion's own −κxyz, odd in Q) pumps without a
  constant force; the off-diagonal sector with the scalar's mass; light (a carrier doublet with hypercharge).
- Stage 2: H/M at the isotropic ghost's band crossing (scenario-dependent); backreaction; other eras; a
  WKB-normalised measure for the fourth-order block.
- Inhomogeneous (k ≠ 0) perturbations: decide whether in scope. Growth once ε is no longer small: not covered.

## History

- 12 Oct: folder opened as the one new (periodic, no lattice) background. Stage 1 (`iso_floquet.py`); stage 2
  first pass (`frw_identity.py`, `frw_ghost.py`).
- 14 Oct: off-diagonal blocks under expansion (`rung1_frw.py`, `frw_offdiag.py`); axion derivation
  (`axion_derive.py`); literature pass; first anisotropic-sector run (`axion_aniso.py`).
- 15 Oct: twin suite (`axion_twins.py`); off-diagonal generator and run (`rung1_axion.py`, `offdiag_axion.py`)
  exposed the inconsistent background; the correction (`axion_dressed.py`) the same evening. Withdrawn (15 Oct,
  evening): "threshold M = m_φ for the anisotropic sector, tachyonic for M < m_φ with rate μ = gQϕ/√(2(1 − r²)),
  e-folds through reheating 2.4 ϕ–9.8 ϕ" and the sign rule "D(m_φ²) > 0, healthy twin unstable on the other side" as
  physics — measured on the M → ∞ background (Part D of both scripts); on the self-consistent background the two
  second-order terms cancel. Parts A–C of both scripts stand as statements about the periodic problem.
- 16 Oct: threshold made exact (`iso_lame_gap.py`); FRW run (`axion_frw_iso.py`); ϕ_amp = 1 replaced by (A11); the
  "linear rate above threshold" item closed. Probe sector derived (`probe_derive.py`, `probe_core.py`); first probe
  run, revised the same night after an independent review: the Part-B twin tangent had a sign error (ζ″ = −ζ + cε,
  CHECK 5b) and the twin's dressed initial ẇ had the ghost's sign — withdrawn: "a strong scalar drives the heavy
  mode irrespective of Krein sign" and "the fermion loses 3 % of its energy with the ghost and 78 % with the twin";
  "spectator at m_p ≥ M" is a hierarchy statement. The long run (`axion_probes_long.py`, 6 min, rerun with two more
  diagnostics, 10 min) and the spread run (`axion_probes_spread.py`). Withdrawn the same night, in order: "a bounded
  three-way exchange and nothing in the thresholds"; "a massless probe lowers the loss threshold from ω²ϕ²Π = 2 to
  at most 1.7"; "a three-body pump needing both probes"; "the loss time as the number to report".
- 17 Oct: inflationary run (`axion_probes_infl.py`); Fourier diagnostic of Part C (`axion_probes_fourier.py`) with
  the mass-label corrections. Probe programme on this background complete; the same core was then run on the dilaton
  background (`dilaton_su2/README.md`, closed 17 Oct).
