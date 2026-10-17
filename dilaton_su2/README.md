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

**Not run yet:** the anisotropic Floquet scan of u″ + e^{δ_a cos ωt}u = 0 (and column B's) over (ω, δ_a)
with δ_a up to O(1); the isotropic background in the (Q, w) form with the controls; the averaged
isotropic equation; the FRW run; Adshead et al.'s numbers in g = M = 1 (their φ_end, Λ, g; not yet
read off).

## Layout

```
code/
  dilaton_derive.py  step 1: what Θ(φ)F² adds to the homogeneous Lee–Wick system, seven PASS/FAIL checks
                     in two columns (Θ on F² only / on the whole operator): Gauss law, truncation, closed
                     forms, auxiliary forms, units, the M → ∞ time-rescaling identity, linearisation on
                     x = y = z (anisotropic sector background-free, first integrals, Mathieu heuristics),
                     FRW.                                                      (sympy; added 2026-10-16)
results/
  dilaton_derive.txt output of code/dilaton_derive.py (73 s, cloud session 16 Oct, sympy 1.14), all PASS
```

## Status / open

- [x] Derivation, both columns (16 Oct).
- [x] Sub-decision on the column: A primary, B control (Joseph, 17 Oct).
- [ ] Anisotropic Floquet scan (a Hill equation in u; no background needed).
- [ ] Isotropic background with the drive: threshold or pump; no-ghost and M > m_φ controls; averaged equation.
- [ ] Adshead et al.'s parameters in g = M = 1 (their appendix and parameter tables not yet read).
- [ ] FRW; the off-diagonal sector (column A: `rung1_axion.py` with Θ on the Yang–Mills part).
