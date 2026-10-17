# dilaton_su2/ — the isotropic SU(2) background with a gauge-kinetic (dilaton) coupling, with the Lee–Wick ghost

Plan item 5b, decided 16 Oct 2026 (choice A of that day): the second periodic background under
item 5, the parity-even twin of `axion_su2/`'s φFF̃. The coupling is the gauge-kinetic function
Θ(φ) F^a_{μν}F^{aμν} of non-Abelian gauge preheating (Adshead, Giblin, Weiner, arXiv:1708.02944,
eq. (1), Θ = e^{φ/M}; the same kinetic function during inflation in Obata, Soda, arXiv:1602.06024,
eq. (3), written I(φ)²), on the homogeneous isotropic SU(2) field A^a_i = δ^a_i Q(t), no lattice,
with the dilaton φ(t) a prescribed homogeneous oscillation, as the axion was. Their coupling scale
M is written Λ here, because M is the Lee–Wick parameter; δ = φ/Λ, Θ = e^δ, kept exact (Adshead et
al. start reheating at φ ≈ 0.2 m_pl, so δ need not be small).

Template: `axion_su2/` (derivation script with PASS/FAIL checks first, then the (Q, w) auxiliary
form with the prescribed drive, the no-ghost and heavy-ghost controls, the averaged coupling, then
the threshold scan off any parametric tongue). Exit criterion as 5a's (every homogeneous sector,
stable or unstable as a function of M/m_φ and the drive, on a background that is itself a solution,
with the controls); cap: one closure paragraph in the Bianchi IX shape.

## Closure (17 Oct; in the shape of the Bianchi IX closing paragraph; the runs it rests on are in the next section)

On the isotropic SU(2) background with a prescribed dilaton Θ = e^{δ_a cos ωt} multiplying F² (column A, the cited
coupling; Joseph's choice of 17 Oct, column B as control), the Lee–Wick ghost's effective mass M²Θ(t) is modulated at
the dilaton's frequency, and the isotropic solution is lost — within a few ghost periods, in every realisation, with
no threshold — wherever ω/M sits in the ghost's own parametric tongue: at ω = 2M for every δ_a tested (0.1–1), over
1.75–2.3 M for δ_a = 1 (column B's tongue is wider, 1.3–2.4 M, and its rates twice A's). Off the tongue the isotropic
solution is steady for sixty Kapitza periods, perturbed seeds identical, the ghost sector carrying a bounded 0.1–2 of the
Yang–Mills energy; the anisotropic sector is independent of the background and has the same tongues and nothing else;
M → ∞ (Yang–Mills with g → gΘ) is steady everywhere. Real probes off the tongue: a massless or heavy doublet is a
spectator; a doublet of mass ≈ gQ_amp exchanges energy with the ghost sector by a factor ~10, periodically and bounded,
with nothing of the kind without the ghost, and the collision coefficient stays at 10⁻³ of its critical value. Contrast
with the axion (5a): there the drive vanished with Q, the ghost's O(1) response dressed background and perturbation
alike, and the one instability was a sharp threshold in the drive amplitude; here the drive multiplies the operator,
nothing cancels, and the one instability is a first-order parametric resonance of the ghost at m_φ = 2M/n — the
Yang–Mills paper's "the ghost is excited only where the background has a frequency that reaches M" in its most literal
form, and the dilaton's own benchmark (Adshead et al.'s φ_end ≈ 0.2 m_pl against an unstated Λ) is not read off here,
so where their model sits in (ω, δ_a) is left open. Item 5b closes with this paragraph (the cap of 16 Oct); open and
not worked on: Adshead et al.'s parameters in these units; the off-diagonal sector with Θ (`rung1_axion.py` with the
factor Θ on the Yang–Mills part); FRW (stage 2 carries over, not run); the ω = M linear resonance.

## Summary (16 Oct; the derivation only, no run yet)

**The sub-decision — does the dilaton multiply the Lee–Wick term too?** Derived in two columns of one
sympy pass (`code/dilaton_derive.py`, all 41 checks PASS, `results/dilaton_derive.txt`): column A,
Θ on F² only (the cited form); column B, Θ on the whole operator (string-frame-natural). What the
derivation settles:

- Both columns: the diagonal ansatz stays a consistent truncation; on it L_A = Θ(½q̇² − V) −
  |q̈ + ∇V|²/2M² and L_B = Θ(½q̇² − V − |q̈ + ∇V|²/2M²); in the (q, w) form of `ym/`, column A is the
  factor Θ(t) on the normal sector's kinetic term and potential and nothing else (L′_A = Θ(½q̇² − V) +
  ẇ·q̇ − w·∇V + M²w²/2), column B is L′_B = Θ L′ + Θ̇ w·q̇; on shell w = (q̈ + ∇V)/M² in both.
- The Gauss law. Column A: unchanged in content (L_YM has no Ȧ₀, so the A₀ equation is Θ times
  the Yang–Mills one plus the Lee–Wick one). Column B: the A₀ equation gains −Θ̇ ∂L_LW/∂Ȧ₀, zero on the
  diagonal ansatz but not off it — the ghost sector's constraint carries the dilaton, and the
  off-diagonal generator `ym/hdym_rung1.py` would have to be regenerated with Θ for column B
  (for column A it needs only the factor Θ on the Yang–Mills part, as `rung1_axion.py` needed only
  the extra potential).
- M → ∞, both columns: in the time τ with dτ = dt/Θ the system is Yang–Mills with g → gΘ(t),
  U_ττ = −Θ²∇V; on x = y = z, Q_ττ = −2Θ²Q³. For a fast drive the leading averaged coupling is
  ⟨Θ²⟩_τ = ∫e^δdt / ∫e^{−δ}dt, positive and exactly 1 for a δ symmetric about 0 (a cosine): the
  pure-Yang–Mills control never inverts the quartic potential and at leading order does not even
  renormalise it. Its only loss of the isotropic solution is the slow-drive parametric resonance of Q
  itself (`axion_su2/` B′).
- **The anisotropic sector is independent of the background**, exactly, in both columns (Q and w₀
  drop out of the linearised equation, by the same identity Σ_i ∂_i(v·Hess V·v) = 0 at x = y = z that
  removed w₀ for the axion). With u = ε′: column A, u″ + M²Θ(t) u = C — the ghost's effective mass
  M²Θ is modulated (never the mass M, which is fixed); column B, (Θu′)′ + M²Θu = CM², i.e. with u =
  e^{−δ/2}v, v″ + (M² − δ̈/2 − δ̇²/4)v = CM²e^{−δ/2} — the ghost's damping is modulated. So the axion's
  Π-cancellation question has nothing to act on here: the sector's only instabilities are the ghost
  pair's parametric tongues at m_φ = 2M/n, **first order in δ_a and not suppressed by Q** (the axion's
  first tongue had rate c₁/4 with c₁ = κQ). First-order Mathieu heuristics (the runs decide): at the
  first tongue ω = 2M the growth is δ_a M/4 (column A) and δ_a M/2 (column B); off the tongues the
  sector is stable with ε drifting linearly (the constant C).
- The isotropic background, column A (g = M = 1): Q″ = w − 2Q³, w″ + (Θ + 6Q²)w = −Θ̇Q′ (compare
  the axion's w″ + (1 + 6Q²)w = κQ²); its 4th-order form and the isotropic tangent are printed in the
  log, with the M → ∞ limits checked. The drive on the background is a friction, −Θ̇Q′, not a
  potential: whether the second-order DC term through the full propagator is a renormalised
  coupling (the axion's λ_eff = 2 − ω²ϕ²Π, with a threshold) or a secular energy pump (a rate) is the
  open question for the isotropic sector and the next script.
- FRW: √−g Θ F·F is independent of a(η), so the dilaton drive is conformal like the axion's; stage 2
  of `axion_su2/` (M → M a(η)) carries over to both columns.

Status (2026-10-17): decided by Joseph — **column A primary** (literature-backed), **column B as the control**.

Recommendation for the sub-decision (Joseph's call): column A as the primary — it is the cited coupling,
leaves the Gauss law and the ghost sector's constraint structure untouched, and the off-diagonal
machinery of `axion_su2/` carries over with one factor; column B as the twin control in the runs,
since its anisotropic equation is as cheap (a Hill equation in u) and it is the one a string-frame
reader will ask about. The two differ in which ghost coefficient the dilaton modulates (effective
mass vs damping) and by a factor 2 at the first tongue; the derivation does not single one out.

**The runs (17 Oct)** (`code/dilaton_iso.py`, `results/dilaton_iso.txt`, 11 min, cloud session, Python 3.13, scipy 1.18;
flat space, g = M = 1, gQ_amp = M/100, Θ = e^{δ_a cos ωt}; checks: column A's isotropic rhs against `probe_core.rhs` with
Θ on, 10⁻¹⁶; column B's (q, w) equations derived by sympy from L′_B at import, equal to A at Θ = 1 and energy-conserving
at constant Θ, 10⁻¹⁴; the anisotropic monodromy against the Mathieu heuristics, three digits).
- *A, the isotropic background, no probes* (ω = 0.5–3.0, δ_a = 0.1, 0.3, 1.0, 8 periods; then 60 at ω = 2.4, 3.0 with
  three seeds). Lost at |Q| = 10 Q_amp within 0.08–2.8 T: ω = 2.0 for all δ_a (A and B); ω = 2.4, δ_a = 1 (A and B); ω =
  1.5 and 3.0 at δ_a = 1 (B only). Every lost row sits in the anisotropic tongue of the same column (Part B); every row
  outside it is steady, |Q|max/Q_amp = 1.00 for 60 periods, seeds identical, E_heavy bounded (max 0.09, 0.07, 1.9 E_bg
  at (2.4, 0.3), (3.0, 0.3), (3.0, 1.0) in A). Mechanism: the isotropic w obeys w″ + (Θ + 6Q²)w = −Θ̇Q′, whose
  homogeneous part is the anisotropic Mathieu equation (6Q² ≈ 10⁻³), sourced by the drive and fed back through
  Q″ = w − 2Q³ — the ghost's tongue in the isotropic sector. No threshold: at ω = 2 the loss is as fast at δ_a = 0.1
  (2.8 T) as the tongue rate predicts (δ_a/4 = 0.025 per unit time, 2.8 T ≈ 600 time units ≈ 15 e-folds). M → ∞ steady
  in every row, as the derivation's ⟨Θ²⟩_τ = 1 requires. B is lost in more rows than A (its tongue is wider and twice
  as fast) and for the same reason.
- *B, the anisotropic map* (one-period monodromy, ω = 0.3–3.2, δ_a = 0.1, 0.3, 1.0). Column A: δ_a = 0.1 unstable at
  ω = 2.00–2.05 (max 0.0250 per unit time at 2.00; heuristic 0.0250); δ_a = 0.3 at 1.90–2.15 (max 0.0744; heuristic
  0.0750) plus a trace at ω = 1; δ_a = 1 at 1.75–2.3+ (max 0.252 at 2.2) and the second tongue 1.0–1.1 and 0.7. Column
  B: 1.95–2.10 (max 0.0500; heuristic 0.0500); 1.75–2.30 (max 0.151); 1.30–2.4+ (max 0.54). A healthy heavy mode has
  the same equation (a Hill equation with positive coefficient): not a ghost statement, the dilaton's anisotropic map.
- *C, the probes* (column A, off the tongue at (ω, δ_a) = (2.4, 0.3) and (3.0, 1.0), m_p/gQ_amp = 0, 1, 3, ghost and
  M → ∞, 8 periods; then the boson at m_p = gQ_amp for 60). m_p = 0 and 3 gQ_amp: ghost and M → ∞ agree (boson energy
  0.4–0.5 and 0.98; fermion n₋ ≈ ½ level crossing and 0.001). **m_p = gQ_amp: the boson's energy rises to 9–10× its
  start with the ghost and stays within 0.1–0.9 without it**, the background dipping to 0.62–0.75 Q_amp meanwhile;
  over 60 periods the swing repeats identically in each 20-cycle block (0.86–9.2, 0.87–9.0, 0.90–9.0 at (2.4, 0.3);
  0.92–9.1, 0.89–9.4, 0.88–9.5 at (3.0, 1.0)) — a periodic, bounded, ghost-mediated exchange, not a pump; c ≤ 0.0006
  against the collision value ¼ throughout; nothing lost. The fermion: identical with and without the ghost.
- *Not read off:* Adshead et al.'s φ_end and Λ in these units (δ_a) and m_φ/M; their 1708.02944 states φ₀ ≈ 0.19 m_pl
  but the excerpt read gives no Λ, so no row is called their benchmark.

**Not run yet** (16 Oct list; Status (2026-10-17): the first three are done in `dilaton_iso.py` above, the averaged
isotropic equation is moot — the instability is a tongue, not an averaged threshold — and the FRW run and Adshead
et al.'s numbers stay open, by the cap): the anisotropic Floquet scan of u″ + e^{δ_a cos ωt}u = 0 (and column B's)
over (ω, δ_a) with δ_a up to O(1); the isotropic background in the (Q, w) form with the controls; the averaged
isotropic equation; the FRW run; Adshead et al.'s numbers in g = M = 1 (their φ_end, Λ, g; not yet read off).

## Layout

```
code/
  dilaton_derive.py  step 1: what Θ(φ)F² adds to the homogeneous Lee–Wick system, seven PASS/FAIL checks
                     in two columns (Θ on F² only / on the whole operator): Gauss law, truncation, closed
                     forms, auxiliary forms, units, the M → ∞ time-rescaling identity, linearisation on
                     x = y = z (anisotropic sector background-free, first integrals, Mathieu heuristics),
                     FRW.                                                      (sympy; added 2026-10-16)
  dilaton_iso.py     the runs: the isotropic background under Θ(t) in columns A (from probe_core, Θ on) and B (sympy from
                     L′_B at import) against M → ∞, 8 and 60 periods with seeds; the anisotropic monodromy map for A and
                     B; the probes on column A with the collision diagnostic. Imports ../../axion_su2/code/probe_core.
                                                                              (float + sympy; added 2026-10-17)
results/
  dilaton_derive.txt output of code/dilaton_derive.py (73 s, cloud session 16 Oct, sympy 1.14), all PASS
  dilaton_iso.txt    output of code/dilaton_iso.py (11 min, cloud session 17 Oct, Python 3.13, scipy 1.18, sympy 1.14)
```

## Status / open

- [x] Derivation, both columns (16 Oct).
- [x] Sub-decision on the column: A primary, B control (Joseph, 17 Oct).
- [x] Anisotropic Floquet scan (17 Oct, `dilaton_iso.py` B): the tongues at 2/n, heuristics confirmed.
- [x] Isotropic background with the drive (17 Oct, `dilaton_iso.py` A): neither threshold nor pump — the ghost's tongue in
      the isotropic sector, lost within a few ghost periods on it, steady off it; M → ∞ steady everywhere.
- [x] Probes (17 Oct, `dilaton_iso.py` C): spectators at m_p = 0 and 3 gQ_amp; a bounded ×10 ghost-mediated exchange at gQ_amp.
- [x] Closure paragraph written (above); item 5b closed by the cap of 16 Oct.
- [ ] Open, not worked on: Adshead et al.'s parameters in g = M = 1; FRW; the off-diagonal sector with Θ; ω = M.
