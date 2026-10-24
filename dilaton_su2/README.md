# dilaton_su2/ — the isotropic SU(2) background with a gauge-kinetic (dilaton) coupling and the Lee–Wick ghost

The parity-even twin of `axion_su2/`: the gauge-kinetic function Θ(φ) F^a_{μν}F^{aμν} of non-Abelian gauge preheating (Adshead, Giblin, Weiner, arXiv:1708.02944, eq. (1), Θ = e^{φ/M}; Obata, Soda, arXiv:1602.06024, eq. (3), I(φ)²) on the isotropic SU(2) field A^a_i = δ^a_i Q(t), no lattice, φ(t) a prescribed oscillation. Their M is written Λ (M is the Lee–Wick parameter); δ = φ/Λ, Θ = e^δ exact (Adshead et al. start reheating at φ ≈ 0.2 m_pl).

Working notes 16–17 Oct 2026; rewritten 23 Oct 2026 from the state at that date. Status: closed 17 Oct (the cap of 16 Oct, one day).

## Result

**Question.** Does the dilaton drive excite the Lee–Wick ghost, and how — a threshold in the drive amplitude (the axion's case) or something else? Column A: Θ on F² only (the cited form); column B: Θ on the whole operator (string-frame-natural). Decision of 17 Oct: **A primary** (literature-backed; Gauss law untouched), **B the control**; they differ in the modulated ghost coefficient (mass vs damping) and by a factor 2 at the first tongue.

**Answer (17 Oct).** With Θ = e^{δ_a cos ωt} on F², the ghost's effective mass M²Θ(t) is modulated and the isotropic solution is lost — within a few ghost periods, every realisation, no threshold — wherever ω/M sits in the ghost's parametric tongue: ω = 2M for every δ_a tested (0.1–1), 1.75–2.3 M for δ_a = 1 (B: wider, 1.3–2.4 M, rates twice A's). Off the tongue the isotropic solution is steady for sixty Kapitza periods, seeds identical, the ghost sector a bounded 0.1–2 of the Yang–Mills energy; the anisotropic sector is background-independent, same tongues, nothing else; M → ∞ is steady everywhere. Probes off the tongue: spectators at m_p = 0 and 3 gQ_amp; at m_p ≈ gQ_amp a bounded ×10 ghost-mediated exchange, collision coefficient 10⁻³ of critical. Unlike the axion's threshold (there the drive vanished with Q), here the drive multiplies the operator and the one instability is the ghost's first-order parametric resonance at m_φ = 2M/n — "the ghost is excited only where the background has a frequency that reaches M", literally.

## Layout

| path | what it is |
|---|---|
| `code/dilaton_derive.py` | sympy, 16 Oct: the derivation below, two columns, PASS/FAIL checks (Gauss law, truncation, closed and auxiliary forms, units, M → ∞ identity, linearisation, FRW); the count is given as "all 41 checks PASS" and as "seven PASS/FAIL checks". |
| `code/dilaton_iso.py` | float + sympy, 17 Oct: runs A, B, C below; column A from `probe_core` (Θ on), B from L′_B at import. Imports `../../axion_su2/code/probe_core.py`. |
| `results/dilaton_derive.txt` | output (73 s, sympy 1.14), all PASS. |
| `results/dilaton_iso.txt` | output (11 min, Python 3.13, scipy 1.18, sympy 1.14). |

Runs: flat space, g = M = 1, gQ_amp = M/100, Θ = e^{δ_a cos ωt}; run from `code/`.

## Derivation (16 Oct)

- The diagonal ansatz is a consistent truncation in both columns; L_A = Θ(½q̇² − V) − |q̈ + ∇V|²/2M², L_B = Θ(½q̇² − V − |q̈ + ∇V|²/2M²). In the (q, w) form of `ym/`: L′_A = Θ(½q̇² − V) + ẇ·q̇ − w·∇V + M²w²/2, L′_B = Θ L′ + Θ̇ w·q̇; on shell w = (q̈ + ∇V)/M².
- Gauss law. A: unchanged in content (L_YM has no Ȧ₀). B: the A₀ equation gains −Θ̇ ∂L_LW/∂Ȧ₀, zero on the diagonal ansatz, not off it — the ghost's constraint carries the dilaton.
- M → ∞, both columns: in τ, dτ = dt/Θ, Yang–Mills with g → gΘ(t), U_ττ = −Θ²∇V; on x = y = z, Q_ττ = −2Θ²Q³. Fast-drive average ⟨Θ²⟩_τ = ∫e^δdt / ∫e^{−δ}dt, positive, exactly 1 for δ symmetric about 0: the control never inverts the quartic nor renormalises it at leading order; its only loss is the slow-drive resonance of Q (`axion_su2/` B′).
- **The anisotropic sector is background-independent**, exactly, both columns (Q and w₀ drop out by Σ_i ∂_i(v·Hess V·v) = 0 at x = y = z). With u = ε′: A, u″ + M²Θ(t) u = C — the effective mass M²Θ is modulated, never M; B, (Θu′)′ + M²Θu = CM², i.e. u = e^{−δ/2}v, v″ + (M² − δ̈/2 − δ̇²/4)v = CM²e^{−δ/2} — the damping. The only instabilities: the ghost pair's tongues at m_φ = 2M/n, **first order in δ_a, not suppressed by Q** (axion: rate c₁/4, c₁ = κQ). Mathieu heuristic at ω = 2M: growth δ_a M/4 (A), δ_a M/2 (B); off the tongues ε drifts linearly (C). A healthy heavy mode obeys the same Hill equation: not a ghost statement.
- Isotropic background, A (g = M = 1): Q″ = w − 2Q³, w″ + (Θ + 6Q²)w = −Θ̇Q′ (axion: w″ + (1 + 6Q²)w = κQ²); 4th-order form and tangent in the log, M → ∞ limits checked. The drive is a friction, −Θ̇Q′.
- FRW: √−g Θ F·F is independent of a(η), the drive conformal; stage 2 of `axion_su2/` (M → M a(η)) carries over; not run.

## Runs (17 Oct, `results/dilaton_iso.txt`)

Checks: A's rhs against `probe_core.rhs`, Θ on, 10⁻¹⁶; B equals A at Θ = 1, conserves energy at constant Θ, 10⁻¹⁴; monodromy against the Mathieu heuristics, three digits.

**A — isotropic background** (ω = 0.5–3.0, δ_a = 0.1, 0.3, 1.0, 8 periods; then 60 at ω = 2.4, 3.0, three seeds). Lost at |Q| = 10 Q_amp within 0.08–2.8 T: ω = 2.0 for all δ_a (A and B); ω = 2.4, δ_a = 1 (A and B); ω = 1.5 and 3.0 at δ_a = 1 (B only). Every lost row sits in its column's tongue (map B); every other row is steady, |Q|max/Q_amp = 1.00 for 60 periods, seeds identical, E_heavy bounded (max 0.09, 0.07, 1.9 E_bg at (2.4, 0.3), (3.0, 0.3), (3.0, 1.0) in A). Mechanism: the w equation's homogeneous part is the anisotropic Mathieu equation (6Q² ≈ 10⁻³), sourced by the drive, fed back through Q″ = w − 2Q³. No threshold: at ω = 2 the loss at δ_a = 0.1 (2.8 T) is as fast as the tongue rate predicts (δ_a/4 = 0.025 per unit time, 2.8 T ≈ 600 time units ≈ 15 e-folds, as logged). M → ∞ steady in every row, as ⟨Θ²⟩_τ = 1 requires. B is lost in more rows, same reason; an averaged equation is moot.

**B — anisotropic map** (one-period monodromy, ω = 0.3–3.2; rates per unit time, heuristic bracketed):

| δ_a | A unstable at | A max rate | B unstable at | B max rate |
|---|---|---|---|---|
| 0.1 | ω = 2.00–2.05 | 0.0250 at 2.00 (0.0250) | 1.95–2.10 | 0.0500 (0.0500) |
| 0.3 | 1.90–2.15, a trace at ω = 1 | 0.0744 (0.0750) | 1.75–2.30 | 0.151 |
| 1.0 | 1.75–2.3+, second tongue 1.0–1.1 and 0.7 | 0.252 at 2.2 | 1.30–2.4+ | 0.54 |

(The Result's "1.75–2.3 M" for A, δ_a = 1 is this "1.75–2.3+"; isotropic A is lost at ω = 2.4.)

**C — probes** (A, off the tongue at (ω, δ_a) = (2.4, 0.3), (3.0, 1.0), m_p/gQ_amp = 0, 1, 3, ghost and M → ∞, 8 periods; the boson at m_p = gQ_amp for 60). m_p = 0, 3 gQ_amp: ghost and M → ∞ agree (boson energy 0.4–0.5, 0.98; fermion n₋ ≈ ½, 0.001). **m_p = gQ_amp: the boson's energy rises to 9–10× its start with the ghost, stays within 0.1–0.9 without**, the background dipping to 0.62–0.75 Q_amp; the swing repeats per 20-cycle block (0.86–9.2, 0.87–9.0, 0.90–9.0 at (2.4, 0.3); 0.92–9.1, 0.89–9.4, 0.88–9.5 at (3.0, 1.0)) — bounded, not a pump; c ≤ 0.0006 against the collision value ¼ (the Result's "10⁻³"); nothing lost. Fermion: identical with and without.

## Open, not worked on

- Adshead et al.'s parameters in g = M = 1 (φ_end, Λ, g; m_φ/M): 1708.02944 states φ₀ ≈ 0.19 m_pl, the excerpt read gives no Λ; no row is called their benchmark.
- The off-diagonal sector with Θ: A, `axion_su2/code/rung1_axion.py` with Θ on the Yang–Mills part; B, `ym/hdym_rung1.py` regenerated with Θ.
- FRW; the ω = M linear resonance.

## History

- 16 Oct: opened after `axion_su2/` with its exit criterion (every homogeneous sector, stable or unstable in M/m_φ and the drive, with controls), cap one paragraph; derivation run.
- 17 Oct: A primary, B control; runs; closing paragraph (Result) written; closed by the cap.

See `axion_su2/README.md` (template, axion contrast), `ym/README.md` ((q, w) form).
