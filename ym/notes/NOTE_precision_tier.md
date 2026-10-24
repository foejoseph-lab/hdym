# Precision sized to the horizon: what `hdym_ball.py` is, why its lost pairs are arithmetic, and the tier that would recover them

13–14 Oct 2026 (two notes, merged 23 Oct 2026). Nothing here changed any code or result; nothing described as a
next step has been built. Reference data: `data/results_precision/2026-10-12/ball_sweep.txt` (E = 0.3, ε = 0.01, ghost
observable R = 10, M = 400); smoke test `tools/arb_smoke.py`.

## 1. The scheme, and what it is in standard terms

The scheme behind the roundoff study, stated in four steps: (1) use floating point as far as it will go; (2) bound it
above and below; (3) where it gets fuzzy, switch to exact integer / fixed point; (4) repeat.

Steps 1–2 are ball (interval) arithmetic: a centre computed in floating point and a radius that is a proven upper bound
on everything the centre has lost. Step 3 cannot mean "exact": a kick evaluates x·(y² + z²), inputs of b bits give a
product of about 3b, a Yoshida step composes three kicks, so one step takes b to roughly 27b and k steps to 27^k · b —
~10⁹ bits after 5 steps, ~10¹⁶ after 10, against the 25 000 steps (100 time units at dt = 0.004) the escape-ordering
claim needs. "Switch to exact" therefore means "switch to more bits, still rounded", and a rounding needs a bound on
what was thrown away, which is the radius of step 2. The scheme collapses, correctly, into ball arithmetic whose
centre precision is raised when the bound says so — multiword ball arithmetic with precision tiers; the questions
are how many tiers and when to switch.

## 2. The bit budget, and why the first tier is double-double

In a chaotic flow a radius grows like e^{λ_eff t}: the centre loses λ_eff t / ln 2 bits per unit time however
carefully each step is rounded, so the precision a certified statement needs is a budget spent linearly over the
horizon of the claim:

    bits available in the centre − bits that must survive at the end = λ_eff · t_horizon / ln 2.

From the committed sweep: for the pairs that get certified, the radius of the reference member is ~4 × 10⁻²⁶ at t = 10
and ~10⁻¹⁴ to 10⁻¹³ at escape (t_esc median ~100) — ~12 decades over ~90 time units, λ_eff ≈ 0.3 per time unit, ~0.43
bits per time unit. A pair is lost when the radius reaches R/2 = 5, another ~14 decades, which at this rate happens for
pairs with t_esc beyond roughly 200; the lost fraction in the sweep is 9–25 %, the tail of the escape-time distribution
(κ = 0.027) beyond that point. A single double (53 bits, start ~10⁻¹⁶ relative) would reach R/2 after ~120 time units —
at the median escape time, with nothing certified. Double-double (106 bits, start ~10⁻³²) buys the ~50 extra bits that
cover the horizon. That is why `hdym_ball.py` carries its centres as double-double: step 3 done once, at the precision
the horizon requires. Not adaptive, but sized by the same arithmetic.

## 3. What of the scheme exists

| step | where | status |
|---|---|---|
| 1. float as far as it holds | every production kernel | — |
| 2. bound above and below | `hdym_ball.py`: componentwise balls (v, rv) on difference numbers (d, rd), radii propagated with proven per-operation bounds (Rump-style 2⁻⁵² per rounding, `INFL`, `ETA`; double-word bounds of Joldes–Muller–Popescu) | rigorous; independently reviewed 11 Oct (`REVIEW_hdym_ball_2026-10-11.md`) |
| 3. more precision where the bound says so | `hdym_ball.py`: double-double centres, fixed tier; `hdym_precision.py`: the empirical version (mantissa sweep, and a fixed-point mode sc = 2^F that rounds every operation to a multiple of 1/sc) | fixed tier |
| 4. repeat | — | not implemented |

The one thing the four steps do not contain is the cure for the wrapping effect: a box enclosure of a chaotic flow
rotates as well as stretches, and a box that must contain a rotated thin box grows by a geometric factor per step on
top of e^{λt}. More bits do not help; with 12 components and λ ≈ 0.3 plain boxes are gone within a few time units.
`hdym_ball.py`'s ellipsoid (Lohner-type) frame for the tangent map (`ell_update`, `ebound`, `step_remainder`) is the
standard remedy and the reason the certificate reaches ~100 time units at all; it is also the part that was hardest to
review, which is where the rigour lives.

Side note from the drift census (`tools/drift_census.py`, 13 Oct): the three precision-study kernels
(`hdym_precision`, `hdym_diffpair`, `hdym_ball`) are built with `--fmad=false` on CUDA and `-ffp-contract=off` on the
CPU (the ball kernel with explicit `fma()` in `ddmul`), because their rounding bounds assume one rounding per written
operation; the five production kernels use NVRTC's default, which contracts. The census's `[CUDA fmad=false]` block is
the one that corresponds to the precision kernels, and there every value channel is bitwise with the CPU. The comment
in `hdym_diffpair.py` that `-ffp-contract=off` "keeps the v channel identical to the production kernel" is true of the
production kernel's C build, not its CUDA build. No certified statement depends on this.

## 4. Why the lost pairs are arithmetic (13 Oct)

The sweep's own diagnostic column raised the alternative: the QR containment radius (`frame radii/|Y| … QR`) reaches
7 × 10⁻² of |Y| at δ = 10⁻²⁰ and 10⁻³–10⁻² at several other δ — the ellipsoid frame losing shape, which inflates the
enclosure through geometry, not rounding, and which a wider centre would not touch (the 11 Oct review flagged the
same, F7). Test, without new arithmetic: the per-step rounding defects of the centre arithmetic (the three `def[...]`
lines of `bpair`) multiplied by 2¹⁰ and 2²⁰ in a scratch copy, the frame's constants (U2 tangent balls, step remainder,
QR containment) unchanged. Sweep of 200 pairs at δ = 10⁻⁸, 10⁻¹⁶, 10⁻²⁴, R = 10, E = 0.3, ε = 0.01, tmax 1500 (cloud CPU;
not committed): lost fraction 0.255 / 0.21 / 0.20 at full precision, 0.41 / 0.32 / 0.35 with 10 bits removed,
0.555 / 0.48 / 0.505 with 20 — the same factor at every δ, log(lost) linear in the bits removed, slope 0.039–0.046 per
bit. The losses are arithmetic. The slope implies λ_eff = κ ln 2 / slope ≈ 0.40–0.48: the paper's ensemble-weighted
separation rate (0.42), measured independently — the lost pairs are the fast-separating tail and spend their budget at
≈ 0.42, while §2's 0.3 is the rate of the certified pairs. The QR column rising (2 × 10⁻³ → 10⁻²) is the frame
enclosing a larger error set, the same lever a wider centre moves the other way. Extrapolated to +53 bits: lost
fraction 15–25 % → 2–3 % (an extrapolation over 53 bits from a 20-bit range; the rerun is the measurement).

## 5. The tier, if built

**Chosen per pair from t = 0, not switched on mid-flight.** The radius at time t is dominated by rounding errors
committed early and amplified by the flow since; promoting the centre to a wider type when rv/|v| crosses a threshold
shrinks only the roundings made after that moment, a vanishing share. So "repeat" means: every pair runs in
double-double, and a pair that is lost is rerun from t = 0 in the wider tier. Only the lost pairs (9–25 %) pay the cost.
(A first version of this section had a radius-triggered promotion rule; withdrawn the same day for this reason.)

**What the tier is.** Not a hand-written triple-double type (159 bits, Fabiano–Muller–Picot 2019, ~80 lines of new
arithmetic with their error constants) but **Arb** (F. Johansson, 2013–; part of FLINT 3 since 2023; Python binding
`python-flint`): a midpoint at a chosen precision plus a proven radius, every operation propagating both, the precision
a per-call parameter — the scheme's step 4 verbatim. What Arb does not supply is the enclosure strategy for a chaotic
flow: integrating the 12-d system in Arb balls naively, the wrapping effect eats the radius within a few time units at
any precision. The frame of `hdym_ball.py` stays. Concretely:

- the centre arithmetic of the member and difference balls in Arb at a chosen precision — in C, linking FLINT's C
  library into the kernel (MSYS2 `mingw-w64-x86_64-flint`; Linux `libflint-dev`), since the frame (`ell_update`, step
  remainder, QR containment) is C inside the per-step loop; `python-flint` for selftests and references, where it can
  replace `mpmath` by a certified reference;
- the frame unchanged, consuming centres as doubles (Arb midpoint rounded to double, one rounding of 2⁻⁵² relative
  charged to the radius) and radii as doubles; the `DDU`, `INFL`, `ETA` charges of the centre arithmetic disappear
  (Arb's radius already contains them), the U2 charges and the frame constants stay; escape bookkeeping, classification
  and output format unchanged;
- precision as a dial: 53 + 0.6 · t_horizon bits plus margin, from the budget of §2 at λ_eff ≈ 0.42 — ~250 bits for a pair
  escaping at t = 320, ~350 at t = 500;
- the lost pairs run 40 000–300 000 steps each at dt = 0.005 (not 25 000).

**Smoke test** (`tools/arb_smoke.py`, 14 Oct): one Yoshida step of the 12-d system at 200 bits in python-flint and in
FLINT's C library through ctypes, against a 120-digit mpmath reference. Windows (Python 3.12, python-flint 0.9.0,
MSYS2 gcc + FLINT): reference inside every ball on both paths, midpoints of the two paths bit-identical, one-step radius
1.6 × 10⁻⁶⁰, 24 µs per step in Python and 6.4 µs in C — ~0.7 min and ~0.2 min of centre arithmetic per 300 000-step pair,
the frame adding its own unchanged cost. Same results on a cloud CPU (43 and 11 µs). Console output not committed. No
obstacle found on either platform.

**What it would buy.** 53 more bits ≈ 16 more decades of radius ≈ 120 more time units of horizon at λ_eff ≈ 0.3: pairs
lost today from t_esc ≳ 200 would be lost from t_esc ≳ 320, cutting the lost fraction from ~15 % to well under 1 % with
κ = 0.027 — and since every uncertified pair in the sweep is a lost one, and the sweep says nothing below δ = 10⁻¹⁶
because the pairs that would be certified there are the long-lived ones, this is what would let the certificate speak
below 10⁻¹⁶, down to the 10⁻⁴⁰ the sweep already samples. (Estimate; the sweep after the change is the measurement.)

**Checks.** A wider centre can only tighten: every certified-different count ≥ the 2026-10-12 value and every
certified-same count ≥ it; f_lo never falls, f_hi never rises; one row the other way means a constant is wrong. Also
rerun a sample of already-certified pairs in the wider arithmetic (exit intervals must overlap, classification must
agree — the monotone check alone cannot fail when only lost pairs are rerun), and note that the sub-step count
n = ceil(…) of a wider centre can occasionally differ, i.e. a slightly different pair-dependent map is certified. Still
needed independently of all this: the Muller–Rideau (ACM TOMS 2022) check on the double-word constants of the existing
kernel, which certifies everything the sweep currently claims.

## 6. Blast radius

Changes: `ym/hdym_ball.py` only (the Arb centre path in C, a per-pair tier flag and the rerun-the-lost-pairs loop,
one selftest), after its migration onto `csrc/hdym_model.h` (phase 2; before/after on the committed 2026-10-12 sweep,
all certified fractions identical); `ym/README.md` (two lines, and `python-flint` / FLINT under Requirements); one new
dated folder `ym/data/results_precision/<date>/` (~5 min at M = 400 on the GPU). Unchanged: the paper, every other
kernel, every committed result, `islands/` (`cert_*.py` is exact Fraction arithmetic on a finite problem and needs no
tiers), `bianchi9/`, `axion_su2/`. Order if taken up: phase-2 migration; the Arb centre path with its selftest and an
independent read of the dozen interface lines; the sweep; the tightening checks.
