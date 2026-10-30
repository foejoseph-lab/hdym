# Math in this repository: an index

Started 28 Oct 2026. A pointer list, not a home: each result stays described in its folder's README or note, and
files stay where they are. Listed here is what holds beyond the numbers of one run: theorems and certificates,
symbolic derivations, and methods that would carry over to another model. Tags as in the top-level README, "Claim
notation": [cert] certified by a checker in this repository, [exact] a symbolic identity checked by the script,
[method] a technique, not a claim. Numerical measurements of the model are not listed.

## Certified

| statement | tag | where |
|---|---|---|
| The Dahlqvist–Russberg orbit at E = 1/2 exists in the certified box, is elliptic with nonzero twist, and its period is the full period (lemma RET); with Moser's twist theorem (an import) the regular island exists | [cert] | `islands/README.md`; untrusted generator `islands/code/cap_dr_orbit.py`; trusted checker `islands/code/cert_lemmas.py`, `cert_check.py`, `cert_twist.py`; formula `islands/code/twist_formula.json` (generator `twist_symbolic2.py`); outputs `islands/results/cert16_*` |
| The diagonal ghost on that orbit is Floquet-unstable at E = 1/2 (real multiplier off the circle; ghost trace in [5.3469864, 5.3471379]) | [cert] | `islands/README.md`; same checker, `--ghost` |
| Reviews of the trusted layer | — | `islands/notes/REVIEW_cert_lemmas_2026-10-11.md`, `REVIEW_cert_twist_2026-10-11.md` |
| Invariant tori, cantori, the role of dimension and of a negative-energy ghost: background and expectations, with what was measured (the planar DR island's edge is sticky) | [est] background; [grid] for the measurement | `islands/notes/NOTE_tori.md`; `islands/code/cantori_residence.py`; `rabbit_holes/ym/notes/RH9_island_full.md` |
| Certified pair escape times (ball arithmetic: a floating-point centre and a proven radius) and the certified fractions of the uncertainty-exponent sweep | [cert] | `ym/hdym_ball.py`; `ym/README.md`, "Certified pair escape times"; `ym/notes/REVIEW_hdym_ball_2026-10-11.md` |

## Symbolic derivations

| statement | tag | where |
|---|---|---|
| Reduction from the 3+1 theory to the 12-d homogeneous system: Gauss law, off-diagonal equations, symmetric criticality, auxiliary-field Hamiltonian | [exact] | `ym/hdym_derive.py` |
| The full nonlinear nine-component system in temporal gauge: fourth-order equations, three Gauss constraints, the Ostrogradsky energy; generated C checked against sympy | [exact] | `ym/hdym_full_derive.py` |
| The linearisation about the diagonal background, a 24-dim system that splits exactly into three 8-dim off-diagonal blocks | [exact] | `ym/hdym_rung1.py` (`ym/rung1_system.pkl`) |
| The general-n law, the 4/(9n − 4) resonant share, the elliptic form of Z₂, the c_j = 0 off-diagonal law 0.0210 E^{5/4} | [exact] derivations, every constant checked | `ym/hdym_analytic.py`, `ym/data/analytic_output.txt` |
| The density of states Z = ∫δ(1 − H) d³q d³p at E = 1, by two independent quadratures, and ν₀ = 8π²/Z | quadrature (not certified) | `ym/hdym_Z.py` |
| Selection rule: on the reflection-symmetric DR orbit a single pair cannot fold at −1, so a collision at −1 forces a quadruplet | [exact] (theorem); numerical check of the half-period lemma | `rabbit_holes/ym/notes/RH7_selection_rule.md`, `RH7_literature.md`; `rabbit_holes/ym/code/rh7_selection_rule.py`; `islands/code/dr_offdiag_halfperiod.py` |
| The exact gauge reduction of the off-diagonal block 12 (moving gauge, colour-rotation zero mode and charge) | [exact] construction; numerical checks | `rabbit_holes/ym/code/rh_reduced_block.py`; `rabbit_holes/ym/notes/RH_reduced_block.md` |
| The canonical Ostrogradsky form of block 12 and the two-mode resonance normal form (parts A and B) | [exact] | `rabbit_holes/ym/code/rh27_averaging.py`; `rabbit_holes/ym/notes/RH27_averaging.md` |
| On x = y = z the diagonal ghost's isotropic mode is the Lamé equation with n = 2 (finite-gap) | [exact] reduction; band edges numerical | `rabbit_holes/ym/code/rh5_lame.py`; `rabbit_holes/ym/notes/RH5_nonplanar.md` |
| Bianchi IX minisuperspace of quadratic gravity: the equations, their regularity, the Ostrogradsky map and pulled-back Ω, the two-ghost (sixth-order) system | [exact] (checks 7–8 of the derivation are floating-point) | `bianchi9/bianchi9_derive.py`, `bianchi9_symp.py`, `bianchi9_six.py`; `bianchi9/README.md` |
| FRW identity, the axion-SU(2) model, its linearisations and the probe equations | [exact] | `axion_su2/code/frw_identity.py`, `rung1_frw.py`, `axion_derive.py`, `rung1_axion.py`, `probe_derive.py` |
| The dilaton-coupled model (Gauss law, truncation, closed and auxiliary forms, M → ∞ identity, linearisation, FRW) | [exact] | `dilaton_su2/code/dilaton_derive.py` |

## Methods

| method | tag | where |
|---|---|---|
| Untrusted generator plus a small trusted checker in outward-rounded dyadic rationals, each lemma elementary, designed for a Lean 4 port | [method] | `islands/README.md` |
| Precision tiers: ball arithmetic whose centre precision is raised when the radius says so | [method] | `ym/notes/NOTE_precision_tier.md`; `ym/hdym_diffpair.py` |
| Claim notation with evidence tags | [method] | top-level `README.md` |

See also hdym-su3's `MATH.md` (the SU(3) embedding and doublet block, the Lamé reduction on SU(3)'s isotropic
orbits, and the single-tangent bias).
