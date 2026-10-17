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

## Summary (read this first; rewritten 2026-10-15 evening after the correction below; the 14 Oct version is kept under "Superseded summary")

**Question.** The same question as the Yang–Mills paper (`paper/`, `ym/`) and the Bianchi IX notes
(`bianchi9/`), asked of the one periodic background in cosmology that is an SU(2) Yang–Mills
configuration: the isotropic field A^a_i = δ^a_i Q(t) of axion-SU(2) inflation after inflation,
when Q oscillates under its quartic self-coupling and the axion φ oscillates at its mass m_φ and
drives it through φFF̃. Fujita et al. assume the field stays isotropic. With the Lee–Wick ghost of
mass M in the gauge sector, does the isotropic solution exist, is it stable, and does anisotropy
grow, as a function of M?

**What closes item 5 (the exit criterion).** For every homogeneous sector of the Lee–Wick field on
this background — the isotropic direction, the two diagonal anisotropic directions, the six
off-diagonal components — stable or unstable as a function of M/m_φ and the drive, measured on a
background that is itself a solution of the Lee–Wick theory, with the controls (no ghost; ghost
heavier than the drive; the healthy twin), in the shape of the Bianchi IX closure.

**The correction (15 Oct).** The first run (14 Oct) and the twin suite (15 Oct) measured the
perturbations on the M → ∞ background Q″ = −2Q³ + κQ², which is the gauge field's response to the
drive with its massless pole only. For m_φ comparable with M the ghost pole's share of that
response is O(1): the true background (the fourth-order isotropic equation, ghost included) jitters
by δQ = κX²Π(m_φ²) with the full propagator Π(w) = 1/(w² − wM²) = −1/w + 1/(w − M²), not −κX²/m_φ².
The inconsistency showed up as a gauge zero mode that failed by O(κ) in the off-diagonal run
(`offdiag_axion.py`, first attempt), because gauge zero modes are exact only on a true solution. On
the self-consistent background:

1. **The anisotropic sector is marginal at second order, for any heavy sector.** Its two
   second-order terms — the drive times the background's jitter, and the drive times the
   perturbation's own fast response — are the same response function Π attached in the two orders,
   and cancel identically. In the periodic 4-d problem: the DC offset c₀ = (c₁²/2)Π cancels the
   tachyon of ε⁗ + ε″ + (c₀ + c₁ cos ωt)ε = 0 to roundoff at ω = 3, 1.5, 2.5, 5 (`axion_dressed.py`
   A1; residuals where |Π| is large scale as c₁², a fourth-order term, 10⁻⁵ of the old law at the
   physical c₁). The sign rule of the twin suite, D(m_φ²) > 0, is a correct statement about the pole
   term alone and not about the physics. What remains in this sector: the ghost pair's parametric
   tongues at m_φ = 2M/n (rate c₁/4 at the first, relative width qϕ ~ 10⁻⁵; 2.5 ϕ per Hubble time at
   the centre at Fujita's benchmark — the one number of the 14 Oct run that stands) and the linear
   resonance at m_φ = M.
2. **The off-diagonal sector is the anisotropic sector.** On the isotropic background the symmetric
   combinations (u_ab + u_ba) obey exactly the anisotropic equation (the traceless-symmetric quintet of
   the diagonal SO(3) is one representation) and the antisymmetric ones are the gauge direction and
   its Gauss partner — proved at the algebra level for arbitrary background values
   (`offdiag_axion.py` A, residuals 10⁻¹⁶) and measured on the true background (B): the whole 8-dim
   block, its symmetric restriction and the anisotropic sector agree to the printed precision at
   ω = 0.8, 2, 3, drift only; the gauge mode stays O(1); Gauss conserved to 10⁻¹⁵. The self-coupling's
   role in this sector is the Yang–Mills mass 2Q² of the gauge triplet, which Gauss removes. The drive
   enters each block as the single potential −κ z u₁₂u₂₁ (`rung1_axion.py`, sympy, all checks PASS).
3. **The isotropic sector — the background's own direction — has a threshold in the drive.** Its
   second-order terms do not cancel (coefficient −2κQ against +κQ): on the free-particle model
   s″ = (3/4)c₁²Π s, confirmed to 0.5 % for the ghost, the healthy twin and the two-pole operator
   (`axion_dressed.py` A2). On the real background that direction is the tangent of the family of
   background orbits, with the quartic potential's own restoring curvature, and the tachyon has to
   beat it: measured threshold ω²ϕ²Π(ω²) between 0.33 and 1.33 at ω = 2 (1.33 slowly unstable),
   above 1.125 at ω = 3 and above 1.8 at ω = 1.5 over four cycles, none for ω < 1 (Part B, the
   (ϕ, ω) scan on the exact background) — of order 1–2, not sharp in this variable; an averaged
   estimate (3ω²ϕ²X²Π against 6⟨Q²⟩ = 2.74X²) puts it at 2. In physical variables: **ϕ² r²/(1 − r²)
   ≳ 1–2, r = M/m_φ**, i.e. r ≳ 0.7–0.8 at ϕ = 1, r ≳ 0.45–0.7 at ϕ = 2, r ≳ 0.9 at ϕ = 0.5, with ϕ the
   axion displacement at the end of inflation in units of f (not stated by Fujita et al.). Above it
   the isotropic solution runs away in its own ghost content by 2.5–6 e-folds per Q period at ω = 2,
   ϕ = 2–2.2 — the background is lost, and the anisotropic sector is carried along; the linear rate
   above threshold is not cleanly measured (the background's own growing mode is seeded by the
   approximate initial data). Below it, everything is marginal.
   Status (2026-10-16): the threshold is sharp and understood — it is λ_eff = 2 − ω²ϕ²Π(ω²) = 0, i.e.
   **ϕ² r²/(1 − r²) = 2** exactly in the averaged theory (the drive renormalises the quartic coupling; the
   Lamé-gap argument is resolved, section "The threshold, second look" below), and the exact background
   confirms it off the ω = 2M tongue: steady at ω²ϕ²Π = 1.0, 1.33, 1.7, lost at 2.0–2.3 (`iso_lame_gap.py`).
   The 1.33 row of Part B that was "slowly unstable" sits on that tongue. Above the threshold there is no
   linear rate to measure: the averaged quartic potential is inverted and the orbit no longer exists; the
   loss takes 0.4–0.9 of a Kapitza period. The open item "linear rate above threshold" is closed by that.
4. **M > m_φ: stable**, every sector, every drive tested (ϕ up to 3). **M → ∞: marginal** (the Hill
   control; pure Yang–Mills).
   Status (2026-10-16): with the amplitude of the drive relative to the Q oscillation allowed to grow (m_φ/gQ
   down to 40, not 100), a second, Yang–Mills-own loss of the isotropic solution appears on both sides of
   m_φ = M: the slow-drive (Kapitza-threshold) one of the isotropic Hill equation, at s_YM = 2ϕ gQ/m_φ
   between 0.2 and 0.3 without the ghost, the same with the ghost far from resonance and lower near it
   (`axion_frw_iso.py` B′). "Stable for M > m_φ" holds below that; above it pure Yang–Mills loses the
   isotropic solution too, so it is not a ghost statement.

**Contrast.** Yang–Mills mechanics (`paper/`): chaotic channel visits kick the ghost with a power
law at every energy. Bianchi IX (`bianchi9/`): the smooth bounce is quiet, the ghost wakes at
H ≈ 0.2 M, and the Einstein branch is a repeller by a power law only. Here: a periodic background
with a fast drive at m_φ; nothing at all unless m_φ > M, and then nothing at second order in the
drive for the anisotropic field, because the ghost dresses the background and the perturbation
alike; the one instability is of the isotropic solution itself, with a threshold set by the drive
amplitude. The 14 Oct sentence "three mechanisms, one pattern: the ghost is excited only where the
background has a frequency that reaches M" survives as the adiabatic half (m_φ > M is necessary);
the "pole sign" mechanism it named for this background does not. The healthy-twin control is
likewise redone: on its own dressed background the twin's anisotropic sector cancels the same way,
and its isotropic threshold sits where its own full propagator is positive, 1/√2 < m_φ/M < 1.
Status (2026-10-16): "nothing at all unless m_φ > M" is the ghost-pole statement; the slow-drive loss
(point 4) exists without the ghost. And in Fujita et al.'s own variables the picture closes: their
(A11) fixes the axion amplitude at the end of inflation, ϕ_end = √3 (H_end/m_φ)(M_Pl/f) (the 14 Oct
assumption ϕ_amp = 1 is replaced), so both thresholds are conditions on m_φ/H_end: the ghost-pole one
m_φ/H_end < √(3/2)(M_Pl/f) r/√(1 − r²) (300 at r = 0.7, 630 at r = 0.9 for f = 0.004 M_Pl), the slow-drive
one m_φ/H_end ≲ 130 (same f). At their benchmark m_φ/H_end = 5×10⁵ and ϕ_end = 8×10⁻⁴: the ghost does
nothing to the isotropic solution, for any M. Where either threshold is reachable the gauge field's
back-reaction on the axion is O(1) (relative size √2 ϕ_end, section below), so there the coupled
system (their A6) is required and the prescribed-axion result is a map, not a verdict.

**Figures.** None yet; the numbers are in `results/axion_dressed.txt` and `results/offdiag_axion.txt`.

**Status.** Working notes, 12–15 Oct 2026, in chronological order below, with the 14–15 Oct
sections kept as written and marked. Flat space; the FRW run, the coupled axion (their A6), the
linear rate above the isotropic threshold, and the m_φ = M resonance are open. Each would change
numbers, not the structure, unless it introduces a new scale between H and m_φ.
Status (2026-10-16): the FRW run is done (`axion_frw_iso.py`; the expansion plays no part in the loss,
which completes within Δa ≈ 0.1–0.4 of the end of inflation, and below threshold only dilutes the drive);
the linear rate is closed (no such rate: the orbit ceases to exist); the m_φ = M resonance is the common
limit of both thresholds (Π singular; both → 0 there) and is not run separately; the coupled axion is
now the one open piece, with a known reason (back-reaction √2 ϕ_end, O(1) in the window).

### Superseded summary (14–15 Oct; kept for the record, see "The correction" above)

The 14 Oct Summary stated: threshold M = m_φ for the anisotropic sector, tachyonic for M < m_φ with
rate μ = gQϕ/√(2(1 − r²)) from two second-order terms (ghost pole + jitter of Q), e-folds through
reheating 2.4 ϕ–9.8 ϕ, and (15 Oct, point 5) the sign rule D(m_φ²) > 0 with the healthy twin unstable
on the other side. **Status (2026-10-15, evening): withdrawn.** All of it was measured on the
M → ∞ background, inconsistent at O(1) for m_φ ~ M; on the self-consistent background the two
terms cancel (point 1 above). The runs themselves (`axion_aniso.py`, `axion_twins.py`) are correct
as statements about that model and are kept; their Parts A–B (checks, the frozen-Q map, the pole
law, the Kapitza threshold, the tongue) stand. (Only this subsection is superseded: the stage-1 and
stage-2 sections that follow the Layout are unaffected, they have no axion drive.)

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
  rung1_axion.py     ym/hdym_rung1.py with the axion drive κ(t) det A: the 24-dim first-order system of the
                     off-diagonal components on the exact driven diagonal background (q, q′, q″, q‴ free,
                     q⁗ from the diagonal equations); checks that the drive enters each block as the single
                     potential −κ z u₁₂u₂₁, the Gauss law is unchanged, κ = 0 reproduces ym/rung1_system.pkl,
                     the gauge zero modes are exact. Writes ../rung1_axion_system.pkl.  (sympy; added 2026-10-15)
  offdiag_axion.py   the off-diagonal blocks on the driven isotropic background: the quintet algebra (the
                     symmetric pair is the anisotropic equation exactly) and the direct run on the exact
                     fourth-order background with the 8-dim block, its symmetric restriction, the Gauss surface
                     and the 4-dim anisotropic sector side by side.                (float; added 2026-10-15)
  axion_dressed.py   the correction of 15 Oct: the background dressed by the ghost's own response (full
                     propagator Π); the periodic problem with the DC offset (anisotropic: cancellation;
                     isotropic: (3/4)c₁²Π) for every heavy sector; the (ϕ, ω) scan of both sectors on the exact
                     background (the isotropic threshold); Fujita's numbers. Imports axion_aniso, axion_twins.
                                                                              (float; added 2026-10-15)
  iso_lame_gap.py    the Lamé-gap argument for the isotropic direction tested three ways (constant shift,
                     Q²-shaped shift on the unshifted and on the consistently shifted orbit), then the
                     exact background off the ω = 2M tongue: threshold ω²ϕ²Π = 2. Imports axion_frw_iso.
                                                                              (float; added 2026-10-16)
  axion_frw_iso.py   the FRW run of the isotropic threshold through matter-era reheating in the (Q, w)
                     auxiliary form, ϕ_end from Fujita et al. (A11); the slow-drive (Kapitza-threshold)
                     loss with its no-ghost control (B, B′); both thresholds as conditions on m_φ/H_end and
                     the coupled-axion size estimate (C).                      (float; added 2026-10-16)
  probe_core.py      the one right-hand side for the homogeneous Lee–Wick SU(2) system with the real-particle
                     probes (16 Oct, see "Probes" below): (q, w) form, both drives (κ(t), Θ(t)) as slots, boson
                     (doublet modulus f) and fermion (singlet two-level ζ) on or off, ghost or twin (σ = ±1),
                     energy and its split (energy_parts), a finite-difference linearisation (tangent_rhs), the
                     fermion's adiabatic number. Shared with dilaton_su2/. Change it and rerun probe_derive.py.
                                                                              (numpy; added 2026-10-16)
  probe_derive.py    the probe sector derived with PASS/FAIL checks (scalar doublet: Gauss law fixes the phase,
                     truncation survives for any Φ; adjoint triplet for the record; Dirac doublet: the diagonal-SO(3)
                     singlet is invariant on every diagonal background, h = [[m, 3gQ/2],[3gQ/2, −m]], no linear
                     coupling to the quintet) and probe_core.py checked against the Euler–Lagrange equations at
                     random points for every switch combination, energy conserved.    (sympy; added 2026-10-16)
  axion_probes.py    the probe run on this background (section "Probes"; revised after the review of 16 Oct): the
                     isotropic threshold with the probe columns on the m_p/gQ_amp axis, boson energy, fermion n₋ and the
                     heavy sector's energy per cycle, twin and M → ∞ controls; the anisotropic sector with the scalar's
                     (g²/2)f² (ghost / twin / M → ∞, tangents checked against the core's linearisation, the test
                     scalar's own growth reported); the ghost's tongue at five probe masses. Imports probe_core.
                                                                              (float; added 2026-10-16, revised)
  axion_probes_long.py the two open items of axion_probes.py with back-reaction: m_p = 0 over 60 periods (the
                     threshold is lowered with the ghost) and m_s ≈ M on the undriven orbit (preheating feeds the
                     ghost's tongue); written to the local-run protocol (CHECK line, CONTROL-FAIL tags, DONE).
                     Imports axion_probes, probe_core.                        (float; added 2026-10-16)
  axion_twins.py     the twin suite: the massless anisotropic direction coupled by the same drive to heavy
                     modes (σ_i, M_i) of either Krein sign — healthy twin, sixth-order two-pole operator,
                     tuned and reversed mass orders, heavy mode alone; the general law D(ω²) = Σσ_i/(M_i²−ω²)
                     derived in the docstring, frozen-Q map (pair analysis for 2n modes), the law per Hubble
                     time, direct two-frequency integration. Imports axion_aniso.py.   (float; added 2026-10-15)
results/
  iso_floquet.txt    its output (run in the cloud session of 12 Oct, Python 3.13; seconds per energy)
  frw_identity.txt   output of code/frw_identity.py
  frw_ghost.txt      output of code/frw_ghost.py
  rung1_frw.txt      output of code/rung1_frw.py (20 s)
  frw_offdiag.txt    output of code/frw_offdiag.py (200 s)
  axion_derive.txt   output of code/axion_derive.py (3 s), all checks PASS
  axion_aniso.txt    output of code/axion_aniso.py (374 s, cloud session 14 Oct, Python 3.13, scipy 1.18)
  axion_twins.txt    output of code/axion_twins.py (213 s, cloud session 15 Oct, Python 3.13, scipy 1.18)
  rung1_axion.txt    output of code/rung1_axion.py (14 s, cloud session 15 Oct, sympy 1.14), all checks PASS
  offdiag_axion.txt  output of code/offdiag_axion.py (23 s, cloud session 15 Oct)
  axion_dressed.txt  output of code/axion_dressed.py (243 s, cloud session 15 Oct)
  iso_lame_gap.txt   output of code/iso_lame_gap.py (~2 min, cloud session 16 Oct, Python 3.13, scipy 1.18)
  axion_frw_iso.txt  output of code/axion_frw_iso.py (~6 min, cloud session 16 Oct, Python 3.13, scipy 1.18)
  probe_derive.txt   output of code/probe_derive.py (6 s, cloud session 16 Oct, sympy 1.14), all checks PASS (CHECK 5b added)
  axion_probes.txt   output of code/axion_probes.py, revised version (10 min, cloud session 16 Oct, Python 3.13, scipy 1.18)
  axion_probes_long.txt output of code/axion_probes_long.py (6 min, cloud session 16 Oct, Python 3.13, scipy 1.18)
rung1_frw_system.pkl the FRW linearisation (K, Gauss matrix G, zero modes, L2; lambdify over
                     x,y,z,px,py,pz,c_j,a,a1,a2); regenerated by code/rung1_frw.py in 20 s
rung1_axion_system.pkl the driven linearisation (K, G, zero modes, q⁗ of the diagonal background; lambdify over
                     x,y,z,px,py,pz,x2,y2,z2,x3,y3,z3,c_j,kappa); regenerated by code/rung1_axion.py in 14 s
```
Also in `code/`: `frw_identity.py` (sympy: the Lee–Wick term on flat FRW in conformal time is the
flat one with M → M a(η)) and `frw_ghost.py` (the isotropic diagonal ghost crossing its band in a
radiation era). Both added 2026-10-12, run in the cloud session (seconds). The 14 Oct scripts were
run in the cloud session too (Python 3.13, sympy 1.14).

## Stage 1 results (12 Oct; float, flat space, no axion, homogeneous perturbations only — these stand)

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
Status (2026-10-15, evening): **superseded** — the law, the direct runs of Part D and the e-fold
estimate were on the M → ∞ background, which misses the ghost's O(1) response to the drive when
m_φ ~ M; on the self-consistent background the two second-order terms cancel (the correction
section below). Parts A–B (checks, the frozen-Q map, the pole law, the Kapitza threshold) and C3
(the tongue at M = m_φ/2) stand.

**The twin suite (15 Oct): the sign rule** (`axion_twins.py`, `results/axion_twins.txt`; Summary
point 5 quotes it). The PU equation is L = ½s′² − (½h′² − ½M²h²) − ½c(t)(s + h)², ε = s + h, a
massless direction and the ghost seeing the drive through one combination. The suite keeps that and
varies the heavy sector, L = ½s′² + Σσ_i(½h_i′² − ½M_i²h_i²) − ½c(s + Σh_i)²: ghost (reference),
healthy twin (σ = +1), ghost + healthy at M_h = 2M (the sixth-order operator with poles 0, M, 2M and
signs +, −, +, the two-pole theory `bianchi9/` left at the β′ = 0 singularity), M_h = M (tuned),
M_h = M/2 (reversed order), and each heavy mode alone (Mathieu).

- *Law.* The fast responses are f = c₁s cos(ωt)/ω² for s and h_i = −σ_i c₁ s cos(ωt)/(M_i² − ω²)
  through each pole; the period average with the jitter term of Q (which cancels the 1/ω² Kapitza
  term exactly, as in the first run) is s″ = +(ϕ²X²ω²/2)·D(ω²)·s, D(w) = Σσ_i/(M_i² − w): tachyonic
  iff D(m_φ²) > 0, rate μ = gQϕ m_φ √(D(m_φ²)/2) in physical units. Ghost: m_φ > M. Healthy twin:
  m_φ < M. Two poles M_h > M: M < m_φ < M_h. M_h < M: m_φ outside [M_h, M]. M_h = M: never. No heavy
  mode: never (stage 1, and the Hill control of the first run).
- *Checks* (Part A, all PASS): the ghost case reproduces `axion_aniso.pairs` to 2×10⁻¹³; c₁ = 0 gives
  the frozen frequencies for every case; the characteristic polynomial is palindromic to 10⁻¹¹; M_h = M
  gives rate ≤ 2×10⁻⁹ off resonance; a heavy mode alone of either sign has the Mathieu rate c₁/(4M)
  at ω = 2M and nothing (10⁻¹³) at 1.7M, 0.9M.
- *Frozen-Q map* (Part B, c₁ = 0.01; pole terms only). Healthy twin: unstable band ω ∈ [0.7075, 1.045]
  (the pole-term prediction 1/√2 < ω < 1 plus the resonance) and the ghost-free 2M tongue; stable above.
  Two poles at M_h = 2: one band [0.955, 2.035] (the window plus both resonances), the 4M tongue;
  M_h = ½: bands [0.3675, 0.56] and [0.9525, 1.365]; M_h = M: isolated grid points at ω = 1/n and 2
  only. Spot rates agree with the pole law to three digits where the law is nonzero.
- *Per Hubble time* (Part C, their benchmark, ϕ = 1): healthy twin 7.7 (r = 1.1), 3.2 (1.5), 1.5 (2.5),
  0.9 (4); 0 for r < 1. Two poles (M_h = 2M): 8.8 (r = 0.55), 6.1 (0.7), 8.4 (0.9); 0 outside
  ½ < r < 1. The ghost column repeats the first run.
- *Direct* (Part D, m_φ/gQ = 100, ϕ = 2, 30 Q periods, slope per cycle over cycles 10–30 / the slow
  equation): ghost ω = 2: 3.283 / 3.278 (3.286 at ratio 200 in the first run). Healthy twin: ω = 0.8:
  3.783 / 3.769; ω = 0.6: 2.149 / 2.148; ω = 1.5: 0.009; ω = 2: 0.034 (drift level 0.053). Two poles
  M_h = 2: ω = 0.8: 0.003; ω = 1.5: 4.924 / 4.918; ω = 3: −0.002. M_h = M: ω = 0.8: 0.006; ω = 1.5:
  0.006. Every sign of D(ω²) is confirmed; where D > 0 the slope matches the slow equation to 0.4 %.

Not run: more than two heavy poles (the law is written for any number); the twin cases in the
physical e-fold estimate (same prefactors as the ghost, Summary point 3, with D in place of
1/(ω² − 1)); the twin cases on the off-diagonal sector.
Status (2026-10-15, evening): **superseded** — the direct runs of Part D, like those of the first
run, used the M → ∞ background; on the self-consistent background the sign rule governs only the
pole term and the anisotropic sector cancels (the correction section below). Parts A–C stand as
statements about the periodic problem.

**The off-diagonal sector with the drive, and the correction (15 Oct, evening)** (`rung1_axion.py`,
`offdiag_axion.py`, `axion_dressed.py`; logs in `results/`; the Summary quotes them).

- *Generator* (`rung1_axion.py`, sympy, 14 s, all checks PASS). `ym/hdym_rung1.py` with L_axion =
  κ(t) det A: det A has no A₀ (Gauss law unchanged), the diagonal ansatz stays a consistent truncation,
  L₂ − L₂(κ=0) = κ(−x u₂₃u₃₂ − y u₁₃u₃₁ − z u₁₂u₂₁) and nothing else, u⁗ coefficients −I, κ = 0
  reproduces `ym/rung1_system.pkl` symbolically. The background is the exact driven diagonal one:
  q, q′, q″, q‴ free, q⁗ from the diagonal fourth-order equations (ghost included) — the gauge zero
  modes u = θ × A are then exact for any κ(t). A first version substituted the M → ∞ background;
  on it the zero modes fail at O(κ), which is what started the correction.
- *The quintet identity* (`offdiag_axion.py` A, 20 random background points, residuals ≤ 2×10⁻¹⁵):
  the symmetric subspace u_ab = u_ba of each 8-dim block is invariant under K and K restricted to it
  is the companion matrix of ε⁗ = −ε″ − κQε; the three blocks agree; the gauge zero mode is in ker G.
  The antisymmetric half is a Pais–Uhlenbeck pair whose light branch is the gauge zero mode itself
  (d ∝ Q(t), frequency √2 Q, the background's own) with the one Gauss condition in it.
- *Direct run* (`offdiag_axion.py` B; exact background, dressed data, ϕ = 1, m_φ/gQ = 100, 6 Q
  periods, ω = 0.8, 2, 3 and block 13 at 2): ln σ_max of the 8-dim block, its symmetric restriction,
  the anisotropic 4-dim sector alongside and the block on the Gauss surface modulo the gauge
  direction agree to the printed precision in every row (e.g. 5.87, 6.56, 6.97, 7.25, 7.48, 7.66 at
  ω = 2: the O(1) transient, then ln t); |ΦZ₀| stays in [0.03, 1.5]; Gauss drift ≤ 1.1×10⁻¹⁴. On the
  M → ∞ background the same block at ω = 0.8 grew 1.69 per cycle in its antisymmetric half, with
  |ΦZ₀| reaching 4×10⁷ by cycle 12 and the Gauss drift 2×10⁵ — the gauge mode is not a solution on
  a background that is not one. That was the signal.
- *The dressing* (`axion_dressed.py`, docstring). The fourth-order isotropic equation gives
  J̃ = Q″ + 2Q³ − κQ² with J̃″ + (1 + 6Q²)J̃ = −(κQ²)″ + …, so the fast jitter of Q under the drive is
  δQ = κX²Π(ω²), Π(w) = 1/(w² − w) = −1/w + 1/(w − 1) (M = g = 1), where the M → ∞ background has
  −1/w only; its DC product with the drive, ⟨κδQ⟩ = ⟨κ²⟩X²Π, is what the 14 Oct "jitter term"
  computed with Π → −1/ω². Measured: the background's period at ω = 0.8 is 198 (exact) against 305
  (M → ∞), the drive on Q being κQ²/(1 − ω²).
- *Periodic problem* (`axion_dressed.py` A). A1: with c₀ = (c₁²/2)Π the massless pair's rate is 0 to
  roundoff for the ghost at ω = 3, 1.5, 2.5, 5, and the residuals of every case where |Π| is large
  scale as c₁² (4.6 c₁² for the ghost at ω = 0.8 and 0.6; 0.66 c₁² for the twin at 1.5; 0.13 c₁² for
  the two-pole operator at 2.5), i.e. fourth order. A2: the isotropic model with c₁_iso = 2c₁,
  c₀_iso = −(c₁_iso²/4)Π gives rate √((3/4)c₁_iso²Π) to 0.3 % (ghost ω = 3, 1.5, 2.5), 4 % (twin
  ω = 0.8), 0.5 % (two poles ω = 1.5); zero wherever Π < 0.
- *Threshold scan* (`axion_dressed.py` B; exact background, dressed data, m_φ/gQ = 100, 4 Q periods,
  period from the smoothed zero crossings): see the table in the log and Summary point 3. Stable
  rows show the free drift 0.69, 0.41, 0.29 after the first-cycle transient in both sectors; unstable
  rows grow by 2.5–6 per cycle in both and the background moves (|Q| reaching 1.2–2.6 Q_amp within
  four cycles) or blows up within a cycle (BLOW-UP rows, stopped at |Q| = 10 Q_amp). The run took
  248 s; an earlier version with a longer period-detection span hung for 50 min on the (2, 2) point,
  whose background runs away by the twelfth cycle — the threshold is slow there, not sharp.
- *What is open from this.* The linear rate above the isotropic threshold: it needs the background on
  its quasi-periodic solution to roundoff (the dressed initial data are O(10⁻³) off, seeding the
  growing mode); a Newton iteration on the Q-period map would do it. The m_φ = M resonance (Π
  singular; a first-order linear resonance of the background's ghost content with the drive). The
  Lamé-gap argument for the isotropic direction, which predicts instability for any positive
  tachyonic coefficient, is contradicted by the scan and is not understood.
  Status (2026-10-16): all three resolved or closed in the section "The threshold, second look" below:
  the Lamé-gap argument assumed a constant shift of h, the dressed term is Q²-shaped and keeps the tangent
  on the band edge, threshold ω²ϕ²Π = 2; no linear rate above it (the orbit is gone); the resonance is
  where both thresholds meet zero.

**The threshold, second look, and the FRW run (16 Oct)** (`iso_lame_gap.py`, `axion_frw_iso.py`; logs
`results/iso_lame_gap.txt`, `results/axion_frw_iso.txt`, cloud session 16 Oct, Python 3.13, scipy 1.18).

- *The Lamé-gap argument and why it fails.* On the undriven orbit the isotropic tangent ε″ + 6Q²ε = 0
  is Lamé with n = 2, k² = ½, h = 3: the lower edge of the band [3, 4.5], the Jordan tangent of the
  orbit family. A constant downward shift of h lands in the gap (1.5, 3), so "any positive tachyonic
  coefficient destabilises" — that was the argument. The one-period monodromy (`iso_lame_gap.py`):
  (a) constant shift −δ⟨Q²⟩: unstable for every δ > 0 (rate 0.35/T already at δ = 0.03); (b) a Q²-shaped
  shift on the unshifted orbit, (6 − δ)Q²: also unstable; (c) the Q²-shaped shift on the consistently
  shifted orbit, ε″ + 3λQ²ε = 0 on Q″ = −λQ³, λ = 2 − δ/3: multiplier exactly 1 for every λ > 0. The
  dressed second-order isotropic equation (A2, (3/4)c₁_iso²Π with c₁_iso = 2κQ) is ε″ + (6 − 3κ²Π)Q²ε = 0
  = ε″ + 3λ_eff Q² ε, λ_eff = 2 − κ_amp²Π = 2 − ω²ϕ²Π(ω²), and the averaged background is Q″ = −λ_eff Q³
  (the Kapitza-averaged quartic with Π in place of −1/ω²): case (c). The tachyonic term is not a shift of
  h but a rescaling of the Q² coefficient, which moves h and n(n+1)k² together along the band edge — as it
  must, since the averaged equation is the linearisation of an autonomous averaged Hamiltonian system,
  whose orbit family always has a Jordan tangent. Threshold: λ_eff = 0, **ω²ϕ²Π = 2 = ϕ² r²/(1 − r²)**,
  sharp. Above it the averaged quartic potential is inverted: the background runs away, not a linear mode
  on a periodic orbit (the orbit no longer exists). The README's averaged estimate of 15 Oct was this
  number; now with the reason.
- *Against the scan.* Exact background in the (Q, w) auxiliary form (Q″ = w − 2Q³, w″ + (1 + 6Q²)w = κQ²,
  checked against the pickle's q⁗ to 5×10⁻¹⁶ and against three Part B rows), m_φ/gQ = 100, 8 periods:
  at ω = 1.8, 2.4, 3.0 the background is steady at ω²ϕ²Π = 1.0, 1.33, 1.7 (amplitude 1.1 → 1.6 Q_amp,
  growing as λ_eff shrinks, the (2/λ_eff)^{1/4} of the averaged orbit) and lost at 2.3 (0.4–0.55 T); at
  exactly 2.0 it is lost at ω = 1.8 and steady at 2.2–2.8 Q_amp at ω = 2.4, 3.0 (the next order holds it).
  At ω = 2.0 — the centre of the ghost pair's first parametric tongue, m_φ = 2M, which Part A of
  `axion_dressed.py` skips and Part B did not — it is lost already at 1.33 (2.6 T) and 1.7. So Part B's
  "1–2, not sharp" was the tongue; off it the threshold is 2.
- *FRW run* (`axion_frw_iso.py`). Conformal time, comoving Q, M → M a(η), drive (g/f)(dφ/dη)Q³ (stage 2
  and `axion_derive.py` CHECK 7): Q″ = a²w − 2Q³, w″ + (a² + 6Q²)w = a κ_phys Q². Reheating as in Fujita
  et al.: the axion oscillates in its quadratic potential and dominates, H = H_end a^{−3/2}, ϕ(a) =
  ϕ_end a^{−3/2}, **ϕ_end = √3 (H_end/m_φ)(M_Pl/f)** from their (A10)–(A11) (the earlier "ϕ_amp not stated
  by them, taken as 1" missed these); the gauge field's energy is not fed back into H. Parameters
  m_φ/H_end, r = M/m_φ, m_Q,end = 5, f = 0.004 M_Pl; a from 1 to 20; initial data Q = Q_end, physical
  Q̇ = 0, the dressed particular solution for w. Results (Part A, 25 points):
  - for m_φ/H_end ≥ 300 the averaged prediction holds row by row: lost iff λ_eff(a = 1) < 0 (r = 0.8, 0.9
    at 300; r = 0.9 at 450; none at 700), the loss complete by a = 1.1–1.25, i.e. within a fraction of an
    e-fold and well before the dilution crossing a_× = (ϕ_end²/(2(ω² − 1)))^{1/3} = 1.2–1.6; the
    expansion plays no part. Below threshold the orbit persists with the larger amplitude of small λ_eff
    (2.0–2.3 Q_end at λ_eff = 0.35–0.37) and relaxes as the drive dilutes.
  - at m_φ/H_end = 100 every point is lost, including r = 1.2 (M > m_φ), and so is the no-ghost control
    (Part B): this is the slow-drive instability of the isotropic Hill equation, pure Yang–Mills' own, at
    s_YM = 2ϕ_end m_Q H_end/m_φ = 0.43. At 200, r = 1.2 is lost with the ghost (a = 2.4) and steady without
    (s_YM = 0.11): near the resonance the ghost lowers that threshold. Part B′ (flat space, m_φ/gQ = 40,
    16 periods) puts the pure-YM threshold between s_YM = 0.2 and 0.3 at ω = 0.5 and 0.833, the ghost
    one the same at ω = 0.5 and below s = 0.2 at ω = 0.833 (lost after 2.2 T). The dressing 1/|1 − ω²| of
    the drive on Q (Part C of `axion_dressed.py`) accounts for part of that, not all; the rest is the
    "combination resonance near ω = 1" of `axion_dressed.py`'s docstring, not characterised.
  - Both thresholds are quadratic in H_end/m_φ (Part C): ghost-pole m_φ/H_end < √(3/2)(M_Pl/f) r/√(1 − r²),
    slow-drive m_φ/H_end < √(2√3 m_Q (M_Pl/f)/0.25) ≈ 130 at f = 0.004 M_Pl. At the benchmark (m_φ/H_end = 5.2×10⁵, ϕ_end = 8.3×10⁻⁴) λ_eff(1) =
    2 − 3.5×10⁻⁷ r²/(1 − r²): nothing, for any M, short of 1 − r² < 3.5×10⁻⁷.
- *The coupled axion, sized* (not run). The Euler–Lagrange equation of (g/f)φ̇Q³ for φ is φ̈ + m_φ²φ =
  −(g/f) d(Q³)/dt. Relative to m_φ²φ, with g²Q⁴ ~ ρ_Q ~ M_Pl²H² (their ε_B = 1 at the end of inflation) and
  ϕ_end from (A11), the back-reaction is of order √6 (M_Pl/f)(H_end/m_φ) = **√2 ϕ_end**: 1.2×10⁻³ at the
  benchmark (consistent with their "the coupling barely moves q"), O(1) wherever either threshold above
  is reachable, since both need ϕ_end = O(1). So the window the prescribed-axion runs locate is exactly
  where the prescribed axion is not trustworthy: their A6–A7 with the ghost (the (Q, w) form plus the φ
  equation, three second-order equations) is the run that decides it. Not done.
- *What changes in the Summary* (Status lines there): point 3's threshold is sharp (2) and the rate item is
  closed; point 4's "stable for M > m_φ" is scoped by the slow-drive loss; the Contrast gets the
  m_φ/H_end form and the benchmark verdict; the coupled axion becomes the one open piece.

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

## Probes: real particles coupled to the ghost (16 Oct, derivation; no run yet)

Joseph, 16 Oct: the "healthy twin" of `axion_twins.py` was a misreading of what was wanted since Bianchi IX —
a representative boson, a fermion, and both, each actually coupled to the ghost, to see energy drain or feed.
The twin stays as a sign control (in the fourth-order operator the heavy pole's sign is forced, so a healthy
pole in its place is a code and mathematics check, not a model; a healthy heavy pole *alongside* the ghost is the
two-pole operator already in the suite). In this space, unlike the Yang–Mills paper's, real particles can be
coupled with nothing invented: anything with SU(2) charge couples to the full gauge field q = p − w, ghost
included, through its covariant derivative. Derived in `probe_derive.py` (all checks PASS):

- **Boson:** a fundamental scalar doublet (the Higgs is one). With A₀ = 0 its Lagrangian depends on the gauge
  field only through Σ(A^a_i)², so its current is diagonal and the diagonal truncation survives for any Φ; the
  Gauss law (the colour charge ρ^a = 0, which on a doublet includes the overall phase) leaves Φ = f(t)Φ₀ with f
  real. So the probe is one real field: L_f = ḟ² − (m_s² + (g²/4)|q|²)f², force −(g²/2) q f² on the gauge field.
  It gives the formerly free anisotropic directions a direct effective mass (g²/2)f² — a second Hill drive at
  twice the scalar's frequency — and the off-diagonal components the same.
- **Fermion:** a classical Dirac doublet (8 complex components). The diagonal SO(3) (spin + colour) has a
  two-dimensional singlet sector, invariant under H_D on *every* diagonal background (each α^iT^i preserves it
  separately): a two-level system iζ′ = h(Q)ζ, h = [[m_D, 3gQ/2],[3gQ/2, −m_D]], Q = (x + y + z)/3, with zero
  colour charge, current J^a_i = j δ^a_i, j = −½ ζ†σ_xζ, force g j on each gauge direction. It sees only the
  trace, so it does not couple to the anisotropic (quintet) sector at linear order; it couples to the isotropic
  sector and to the ghost through Q. Kept in phase space: no particle-number basis is needed (the adiabatic
  number is offered as a diagnostic only).
- **Adjoint triplet** (record only): consistent only with its direction along an axis, and then anisotropic by
  construction (mass g²(x² + y²)f² for n = ẑ). Not in the core.
- **Light:** a photon is not an SU(2) particle and its homogeneous mode is a free constant field; it enters only
  through a carrier charged under both groups (a doublet with hypercharge). Not done; the same limit Bianchi IX hit.
- **Twin:** σ = −1 flips the heavy sector to a healthy oscillator in the (p, w) variables, L_twin = L′ + |ẇ|² −
  M²|w|²; defined with Θ = 1 only (with the dilaton on, any sign flip of the w sector makes the kinetic matrix
  degenerate at some Θ), so the dilaton's twin control is done on the periodic-coefficient equations instead.

Plan (Joseph, 16 Oct): per background and sector, the columns no probe / ghost + boson / ghost + fermion /
ghost + both, and the same with the twin in place of the ghost; `axion_su2/` first (its existing results are the
no-probe column), then `dilaton_su2/` on the same core; Bianchi IX (the coupled "all three" run, item 4) only when
item 4 is reopened.

**First probe run on this background (16 Oct, late)** (`axion_probes.py`, `results/axion_probes.txt`, 7 min, cloud
session, Python 3.13, scipy 1.18; flat space, g = M = 1, m_φ/gQ = 100, the exact (Q, w) background with the probe
at 10 % of the background's energy; the isotropic truncation of `probe_core.rhs` is used for speed and checked
against it at random points, 4×10⁻¹⁶).

- *A. The isotropic threshold with the probes.* At probe mass m_p = M (and 3M) every column — none, ghost + boson,
  ghost + fermion, ghost + both, twin + both — is identical to the no-probe row to the printed digits: steady at
  ω²ϕ²Π = 1.0, 1.7, 2.0 and lost at 2.3 (0.45 T at ω = 2.4, 0.55 T at ω = 3.0), the probe's energy at 1.00 of its
  start through all eight cycles. A heavy probe is adiabatic on this background and neither drains nor feeds.
  At m_p = 0 (the probe's only mass is the one the background gives it) the picture changes: with the ghost the
  background's own amplitude swings between 0.3 and 2.3 Q_amp cycle to cycle and the probe's energy between
  −0.6 and 3.1 of its start (the classical Dirac energy is indefinite); with the twin, or with no heavy mode at
  all (M → ∞), the background stays flat (1.33–1.6 Q_amp) while the probe's energy swings within about ±1 of its
  start — its own bounded exchange with the quartic oscillator. So a massless probe and the ghost exchange energy
  through the background in a way neither does alone, bounded over 8 cycles; the loss threshold is unmoved (lost
  at 2.3, steady at 2.0, in every column). Not characterised further.
- *B. The anisotropic sector with the scalar's (g²/2)f².* Linear tangent on the undriven orbit (Q_amp = 1, no axion),
  the scalar a test field at ε_p of the background's energy, 30 periods, renormalised each period, ghost / twin /
  M → ∞ side by side (`aniso_scalar`). At ε_p ≤ 0.1 nothing grows in any column (rates ≤ 0.1 per period, the free
  drift). At ε_p = 0.5 the heavy mode grows, and the twin grows as much or more than the ghost: per period, ghost
  1.13 / twin 0.65 (m_s = 0), 0.68 / 1.32 (0.3), 0.11 / 1.39 (0.5), 0.23 / 0.44 (0.7), 0.00 / 0.21 (0.8),
  0.03 / 0.57 (1.2), nothing for m_s ≥ 1.5; the M → ∞ Hill control stays at ≤ 0.06 throughout. So a strong
  oscillating scalar is a parametric drive for the heavy mode irrespective of its Krein sign — not a ghost
  statement — and the resonance is not at the naive 2m_eff = 2M (its location is not characterised). The rows
  m_s = 1.0 and 1.1 are to be read differently: there the test scalar is itself pumped by the background (its own
  growth 0.36 and 0.31 per period — its parametric resonance with Q², the gauge preheating of a doublet), so f²
  grows without bound in the test-field setup and all three columns inherit it (ghost 30–80, twin ≈ 1, M → ∞
  ≈ 1 per period); the coupled run with back-reaction is what would decide that case.
- *C. On the ghost's tongue* (ω = 2M, ω²ϕ²Π = 1.0, m_p = M): the scalar's energy is untouched (1.00) with ghost or
  twin; the fermion loses 3 % of its energy over 8 cycles with the ghost and 78 % with the twin (both probes
  together: 1 % and 39 %) — a drain in both cases, larger for the healthy heavy mode. Not understood; recorded.

Status (2026-10-16, late night, after the review of `axion_probes.py`): the first run above is superseded by the
revised run below on four points — (i) its Part-B twin tangent had a sign error (ζ″ = −ζ **+** cε, derived in
`probe_derive.py` CHECK 5b) that faked the twin's growth: with it fixed the scalar drives the ghost pair and not
the twin, a ghost-specific result; (ii) m_p = M at m_φ/gQ = 100 is 100 gQ_amp and decoupled by construction
(coupling ~10⁻⁴ of m_p², force ~10⁻⁵ of 2Q³), so "spectator at m_p ≥ M" is a hierarchy statement — the axis is
m_p/gQ_amp and the decoupling is at ≈ 3 gQ_amp; (iii) the fermion's energy is bounded by construction (ζ†ζ is
conserved) and says nothing — the observable is its weight on the negative-energy level, n₋; (iv) the twin's
dressed initial ẇ had the ghost's sign, so the Part-C "78 % drain with the twin" was a transient and is withdrawn.
The energy split (massless / heavy / potential / probes) is exported from `probe_core.energy_parts`; with a
drive on, the pieces are not separately conserved, so the split is bookkeeping of where the drive's work went.

**Revised run (16 Oct, late night)** (`axion_probes.py` revised, `results/axion_probes.txt`, 10 min, cloud session;
regimes labelled: Parts A and C at gQ_amp = M/100, Part B at gQ_amp = M; the hand-coded tangents checked against a
finite-difference linearisation of the core at 20 random points and 20 along the orbit, errors 10⁻¹⁰):

- *A, the isotropic threshold on the m_p/gQ_amp axis* (0, 1, 3; 100 = M as the decoupled control). Eight-cycle
  thresholds unchanged in every column (steady at 2.0, lost at 2.3; the 60-cycle L1 statement stands). The boson at
  m_p = 0 trades 0.4–5× its energy with the ghost background, at gQ_amp 0.3–1.5×, at 3 gQ_amp 0.98–1.00: the probe
  decouples at m_p ≈ 3 gQ_amp. The fermion's n₋ at m_p = 0 sits near ½ for ghost, twin and M → ∞ alike — level
  crossing as Q passes through zero, not production, so n₋ is meaningful only for m_p ≳ gQ_amp; at m_p = gQ_amp it
  is 0.1–0.6 with the ghost, 0.4–0.5 with the twin, ≈ 0.5 at M → ∞; at 3 gQ_amp 0.001 (ghost), 0.07 (twin). The
  ghost sector's cycle-averaged energy is O(1–30) E_bg and negative at ω²ϕ²Π = 1.7 and O(10²) at 2.0 (the dressed
  amplitude), the twin's the same magnitudes positive, M → ∞ zero: the heavy sector carries many times the
  Yang–Mills energy on this background in either sign.
- *B, the anisotropic sector with the scalar, twin sign fixed.* At ε_p ≤ 0.1 nothing grows in any column. At
  ε_p = 0.5 **the scalar drives the ghost pair and not the twin**: per period ghost / twin / M → ∞ = 1.13 / −0.00 /
  −0.01 (m_s = 0), 0.68 / −0.01 / −0.01 (0.3), 0.11 / 0.11 / 0.06 (0.5), 0.23 / 0.00 / −0.00 (0.7), nothing at
  m_s ≥ 0.8 away from the preheating rows. A strong oscillating boson is a parametric drive for the ghost pair and
  not for a healthy pair of the same mass (the Krein sign enters the pair's equation through the sign of c in the
  ζ equation). The m_s = 1.0, 1.1 rows stay invalid as test-field rows (the scalar's own growth 0.36, 0.31; see
  L2 for the back-reacting case).
- *C, on the ghost's tongue, twin initial data fixed.* The boson's energy: untouched at m_p ≥ 3 gQ_amp with either
  heavy mode. The fermion's n₋ at m_p = M/2 (the probe's own tongue, 2m_p = M): with the twin 0.013 → 0.62 over
  eight cycles, with the ghost 0.000 → 0.009 at the same drive (the twin's dressed amplitude is larger, 1.56 vs
  1.14 Q_amp, a caveat); at m_p = M nothing with either. A healthy heavy mode driving a fermion between levels at
  2m_p = M is pair production by a strong oscillating field and expected for either sign; the asymmetry in rate
  is recorded, not understood.

What this says for item 6, so far: on this background a real probe at or above the Lee–Wick mass is a spectator;
the ghost-specific effects found are a bounded three-way energy exchange with a massless probe and nothing in the
thresholds. What is open: the m_p = 0 exchange over longer times and in the coupled (back-reacting) scalar case;
the heavy-mode resonance of Part B as a function of m_s with back-reaction; the twin's drain in Part C.
Status (2026-10-16, night): the first two of those are run (`axion_probes_long.py`, next paragraph) and the
"bounded" and "nothing in the thresholds" above are withdrawn for massless probes: the exchange is not bounded,
and the threshold is lowered.

**The long run (16 Oct, night)** (`axion_probes_long.py`, `results/axion_probes_long.txt`, 6 min, cloud session,
Python 3.13, scipy 1.18; same scope as above, all with back-reaction, all in phase space).

- *L1. m_p = 0 over 60 Kapitza periods* (both probes at 10 % of the background's energy; m_φ/gQ = 100). Controls
  steady: M → ∞ with no probe 1.32–1.37 Q_amp throughout; M → ∞ with both probes the same, the probe's energy
  within −0.25 … 1.4 of its start for all 60 cycles. **With the ghost the isotropic solution is lost below the
  no-probe threshold:** at ω²ϕ²Π = 1.7 it is lost at 32 T (ω = 2.4) and 47 T (ω = 3.0), with the background's
  amplitude wandering between 0.02 and 1.5 Q_amp first and the probe's energy reaching 2000× its start before the
  loss; at 2.0 it survives the 60 cycles with the same wandering (0.1–3.2 Q_amp). The twin is steady at 1.7
  (1.40–1.49 Q_amp, probe energy bounded) and lost at 2.0 (ω = 2.4, 50 T; ω = 3.0 survives with amplitude
  rising to 2.0), i.e. at the no-probe threshold itself, where everything is marginal. So: a massless real
  particle coupled to the ghost lowers the loss threshold of the isotropic solution from ω²ϕ²Π = 2 to at most
  1.7 on this time scale, with an unbounded energy exchange; the same particle coupled to a healthy heavy mode,
  or to Yang–Mills alone, does neither. Not characterised: the threshold's new value, the dependence on the probe's
  initial energy, which probe does it (both were on).
- *L2. m_s ≈ M with back-reaction, undriven orbit, 100 periods* (Q_amp = 1, T = 5.244, scalar only). Controls
  (no scalar) steady, tangent drift 0.05 → 0.01 per period, ghost = M → ∞ exactly. At m_s = 0.8, 0.9, 1.1, 1.2
  the ghost and M → ∞ columns agree: the scalar's energy moves within 0.15–16× of its start (the back-reacted
  version of the Part-B preheating at 1.1), the background within 3 %, the anisotropic tangent marginal in
  both. At m_s = 1.0 they do not: M → ∞ preheats the scalar to 20–70× its start (ε_p = 0.01) and the background
  gives up 5–30 % of its amplitude, tangent marginal; **with the ghost the anisotropic tangent grows at 0.8 then 6.4
  per period and the run ends between periods 30 and 60 (|Q| > 10 or non-finite)**, and at ε_p = 0.1 the
  background's amplitude rises to 3.3 Q_amp and stays there (M → ∞: 0.7–0.97). So the preheating of a doublet
  of mass ≈ M by the gauge oscillation, which saturates harmlessly without the ghost, feeds the ghost's
  anisotropic tongue with it. One point in m_s; the width of the window is not measured.
- *The twin on the Q_amp = 1 orbit is not a valid control* and its L2 rows are blank by construction: in the
  (p, w) twin the heavy sector's effective coefficient is M² − Hess V, negative for 6Q² > M², i.e. for Q > 0.41 M
  — the healthy "twin" is itself tachyonic on that orbit. It is a control only where Q ≪ M (L1, Part A). Noted
  in `probe_core.py`'s docstring.

For item 6, revised: on this background a probe at or above M is a spectator; a massless probe — or one whose
mass lets the background preheat it — is not: coupled to the ghost it lowers the isotropic loss threshold and
feeds the ghost's anisotropic tongue, and neither happens with a healthy heavy mode or with Yang–Mills alone.
That is the first ghost-specific statement of the probe programme.
Status (2026-10-16, late night): stands after the review; the revised run adds a second one — a strong oscillating
boson drives the ghost pair and not a healthy pair (Part B, twin sign fixed).

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
      Status (2026-10-15): the twin suite is done (`axion_twins.py`): the sign rule D(m_φ²) > 0, the
      healthy twin unstable on the other side, the sixth-order two-pole window; Summary point 5. The
      open list is unchanged.
      Status (2026-10-15, evening): the off-diagonal sector with the drive is done (`rung1_axion.py`,
      `offdiag_axion.py`) and exposed the inconsistent background; the correction (`axion_dressed.py`)
      withdraws the 14–15 Oct conclusions and replaces them (Summary, rewritten): anisotropic and
      off-diagonal sectors marginal at second order, the isotropic solution itself unstable above a
      drive threshold ϕ²r²/(1−r²) ≳ 1–2, everything stable for M > m_φ. Open: the linear rate
      above the threshold (needs the periodic background to roundoff), the m_φ = M resonance, the
      FRW run, the coupled axion.
      Status (2026-10-16): the threshold is ϕ²r²/(1−r²) = 2, sharp, with the Lamé-gap argument resolved
      (`iso_lame_gap.py`); the FRW run is done (`axion_frw_iso.py`): the expansion plays no part, the loss
      completes within Δa ≈ 0.1–0.4; a second, pure-Yang–Mills slow-drive loss at s_YM = 2ϕ gQ/m_φ ≈ 0.2–0.3
      on both sides of m_φ = M; ϕ_end fixed by their (A11), so both thresholds are conditions on m_φ/H_end
      (300–630 and ≈130 at f = 0.004 M_Pl), unreachable at their benchmark by many orders. The rate item
      is closed (no linear rate exists), the resonance is the common limit of both thresholds. Open: the
      coupled axion (their A6 with the ghost), required in the window because the back-reaction is
      √2 ϕ_end = O(1) there; the near-resonance ghost lowering of the slow-drive threshold (not
      characterised).
- [ ] Probes (16 Oct): derivation and core done (`probe_derive.py`, `probe_core.py`); the runs with the probe
      columns on this background are next (isotropic threshold table; anisotropic and off-diagonal sectors with
      the scalar's (g²/2)f² mass). Then the same on `dilaton_su2/`.
      Status (2026-10-16, night): the long run (`axion_probes_long.py`) finds the first ghost-specific probe
      results: a massless probe lowers the isotropic loss threshold (lost at ω²ϕ²Π = 1.7 within 32–47 T, steady
      with the twin and at M → ∞) and the m_s ≈ M preheated scalar feeds the ghost's anisotropic tongue. Open:
      the new threshold's value and its dependence on the probe energy; which probe does it; the m_s window.
      Status (2026-10-16, late): first run done (`axion_probes.py`, section "Probes"): heavy probes are spectators,
      thresholds unmoved; at m_p = 0 a bounded ghost–probe–background exchange; a strong scalar drives the heavy
      mode irrespective of sign. Open: the m_p = 0 exchange and the m_s ≈ M preheating case with back-reaction;
      the off-diagonal sector with the scalar's mass; then `dilaton_su2/`.
- [ ] Inhomogeneous (k ≠ 0) perturbations on the isotropic background: decide whether in scope.
