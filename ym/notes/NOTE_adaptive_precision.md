# Precision sized to the horizon: the scheme, what of it exists, and what a third word would change

2026-10-13. Written up from the author's description of the scheme he had in mind for the
roundoff question before `hdym_ball.py` was written, so that it is on record in his terms and
so that the next step, if taken, is built from it. Nothing here changes any code or result.

## 1. The scheme, as stated

1. Use floating point as far as it will get you.
2. Bound it above and below.
3. Where floating point gets fuzzy, switch to exact integer / fixed point.
4. Repeat.

## 2. What it is in standard terms

Steps 1–2 are ball (interval) arithmetic: a centre computed in floating point and a radius
that is a proven upper bound on everything the centre has lost. Step 3 is the one that needs
adjusting, and the adjustment is forced by the flow, not by taste.

Exact arithmetic has no floor to stand on here. A kick evaluates `x·(y² + z²)`: inputs of `b`
bits give a product of about `3b`, a Yoshida step composes three kicks, so one step takes `b`
to roughly `27b` and `k` steps to `27^k · b`. From 53 bits that is ~10⁹ bits after 5 steps,
~10¹⁶ after 10, and the escape-ordering claim needs 25 000 steps (100 time units at
`dt = 0.004`). No storage helps; the representation grows exponentially in the step count.
So "switch to exact" has to mean "switch to more bits, still rounded", and the moment one
rounds one needs a bound on what was thrown away — which is the radius of step 2. The scheme
therefore collapses, correctly, into: *ball arithmetic whose centre precision is raised when
the bound says so*. That is a precise and standard idea (multiword ball arithmetic with
precision tiers); the only question is how many tiers, and when to switch.

## 3. The bit budget, and why the first tier is double-double

In a chaotic flow a radius grows like `e^{λ_eff t}`, i.e. the centre loses `λ_eff t / ln 2`
bits per unit time regardless of how carefully each step is rounded. The precision a certified
statement needs is therefore not concentrated at a few "fuzzy" steps; it is a budget spent
linearly over the horizon of the claim:

    bits available in the centre  −  bits that must survive at the end  =  λ_eff · t_horizon / ln 2.

Measured, from the committed sweep (`data/results_precision/2026-10-12/ball_sweep.txt`, E = 0.3,
ε = 0.01, ghost observable R = 10, M = 400): the radius of the reference member is `~4e-26` at
`t = 10` and `~1e-14` to `1e-13` at escape (`t_esc` median ~100) for the pairs that get
certified; that is ~12 decades over ~90 time units, `λ_eff ≈ 0.3` per time unit, ~0.43 bits
per time unit. A pair is *lost* when the radius reaches `R/2 = 5`, another ~14 decades, which at
this rate happens for pairs with `t_esc` beyond roughly 200; the lost fraction in the sweep is
9–25 %, which is the tail of the escape-time distribution (`κ = 0.027`) beyond that point.

A single double (53 bits) starting at `~10⁻¹⁶` relative would reach `R/2` after ~120 time
units — at the median escape time, with nothing certified. Double-double (106 bits, start
`~10⁻³²`) buys the ~50 extra bits that cover the horizon of the claim. That is why
`hdym_ball.py` carries its centres as double-double: it is step 3 of the scheme, done once,
at the precision the horizon requires. Not adaptive, but sized by the same arithmetic.

## 4. What of the scheme exists already, and the one thing it did not contain

| step | where it is | status |
|---|---|---|
| 1. float as far as it holds | every production kernel | — |
| 2. bound above and below | `hdym_ball.py`: componentwise balls `(v, rv)` on difference numbers `(d, rd)`, radii propagated with proven per-operation bounds (Rump-style `2⁻⁵²` per rounding, `INFL`, `ETA`; double-word bounds of Joldes–Muller–Popescu) | rigorous, reviewed 2026-10-11 |
| 3. more precision where the bound says so | `hdym_ball.py`: double-double centres, fixed tier; `hdym_precision.py`: the *empirical* version of steps 1–3 (mantissa sweep, and a fixed-point mode `sc = 2^F` that rounds every operation to a multiple of `1/sc`) | fixed tier, not adaptive |
| 4. repeat | — | not implemented |

The one idea the four steps do not contain is the cure for the wrapping effect. A box
enclosure of a chaotic flow does not only stretch along the unstable direction; it rotates, and
a box that must contain a rotated thin box grows by a geometric factor per step *on top of*
`e^{λt}`. More bits do not help, because the loss is in the shape of the enclosure, not in the
arithmetic; with 12 components and `λ ≈ 0.3`, plain boxes are gone within a few time units.
`hdym_ball.py`'s ellipsoid (Lohner-type) frame for the tangent map (`ell_update`, `ebound`,
`step_remainder`) is the standard remedy and is the reason the certificate reaches ~100 time
units at all. It is also the part that was hardest to review, which is where the rigour lives.
So: the scheme is the right scheme; the code is the scheme plus the wrapping cure; what the
code lacks relative to the scheme is step 4.

Side note from the drift census (`tools/drift_census.py`, 2026-10-13): the three
precision-study kernels (`hdym_precision`, `hdym_diffpair`, `hdym_ball`) are built with
`--fmad=false` on CUDA and `-ffp-contract=off` on the CPU (the ball kernel with explicit
`fma()` calls in `ddmul`), because their rounding bounds assume one rounding per written
operation; the five production kernels are built with NVRTC's default, which contracts. The
census's `[CUDA fmad=false]` block is therefore the one that corresponds to the precision
kernels' production builds, and there every value channel is bitwise with the CPU. One
consequence worth knowing: the comment in `hdym_diffpair.py` that `-ffp-contract=off` "keeps
the v channel identical to the production kernel" is true of the production kernel's C build
and not of its CUDA build, which differs from both by the few ulps the census measures. No
certified statement depends on this (the certificate is about the integrator map with its own
rounding), but a float-vs-certified comparison on the GPU should expect it.

## 5. Step 4: a third word — chosen per pair from t = 0, not switched on mid-flight

*Revised 2026-10-13 after an independent read (two objections, both upheld; the first version
of this section had a promotion-by-radius rule and is superseded).*

Why a radius-triggered switch does not work here. The radius at time `t` is dominated by
rounding errors committed early and amplified by the flow at `e^{λt}` since. Promoting the
centre to a wider type when `rv/|v|` crosses a threshold shrinks only the roundings made after
that moment, a vanishing share of the enclosure; the error already inside keeps growing at the
same rate. A member promoted halfway through its budget gains almost nothing. The 53 extra bits
of a third word materialise only if the wider centre is used from the first step. That is §3's
own argument (precision sized to the horizon in advance) applied to the tier, and it rules
out "promote when the bound says so" for errors that are amplified exponentially; that rule is
right for errors made locally and wrong for these.

What survives of step 4 is adaptivity at the level of the *pair*, not the step: every pair runs
in double-double; a pair that is lost is rerun from `t = 0` in the wider tier. "Repeat" keeps
its meaning, the payoff is the same as always-triple, and only the lost pairs (9–25 % in the
sweep) pay the triple-word cost.

What the tier is:

- a triple-double type `td` (159 bits) with `add`, `sub`, `mul`, `scl` and proven relative
  error bounds. The algorithms and constants come from Fabiano, Muller and Picot, *Algorithms
  for triple-word arithmetic*, IEEE Trans. Computers 68 (2019); the per-operation charge
  replaces `DDU = 2⁻¹⁰¹` with the triple-word analogue (of order `2⁻¹⁵⁴`), to be taken from
  that paper, not derived here;
- a second kernel instantiation with `td` centres, selected per pair by a flag; the member's
  initial centre is the double-double one with a zero third word (exact);
- the nonlinear kick in `td` is one more instantiation of `HDYM_DEFINE_KICK_NL` in
  `csrc/hdym_model.h` (`ADD = tdadd`, …), not a new copy;
- the ellipsoid frame, the step remainder, the escape bookkeeping and the classification are
  unchanged: they consume centres as doubles and radii as doubles;
- a selftest: a lost pair rerun in `td` against a 60-digit `mpmath` orbit (the existing
  selftests 2, 2b, 3 give the pattern), checking that the certified exit steps are unchanged
  and the radius at escape is smaller;
- FMA policy unchanged: `--fmad=false` / `-ffp-contract=off`, explicit `fma()` only.

### 5a. First, find out what actually loses the pairs (diagnostics only, no new arithmetic)

The payoff estimate below assumes the lost pairs die of arithmetic: long-lived orbits that run
out of bits. The sweep's own diagnostic column says that may not be so: the QR containment
radius (`frame radii/|Y| … QR`) reaches `7e-02` of `|Y|` at `δ = 10⁻²⁰` and `10⁻³`–`10⁻²` at
several other `δ`. That is the ellipsoid frame losing *shape*, which inflates the enclosure
through geometry, not through rounding, and a wider centre does not touch it (the independent
review of 2026-10-11 flagged the same, F7).

So before any tier is built: for a handful of lost pairs, log which contribution dominates the
radius at the moment of loss — propagated early error, new rounding, step remainder, or QR
containment. `bpair` already has every one of these in hand at each step; this is a print,
not a change. If the answer is arithmetic, §5 applies. If it is QR containment, the remedy is
in the frame (re-orthogonalise more often, or carry the frame itself in double-double), and the
third word would be the wrong tool.

*Status (2026-10-13, later): 5a answered without the planned diagnostic, by making the
committed kernel worse instead of writing new arithmetic. The per-step rounding defects of the
centre arithmetic (the three `def[...]` lines of `bpair`) were multiplied by 2¹⁰ and 2²⁰ in a
scratch copy; the frame's own constants (`U2` tangent balls, step remainder, QR containment)
were unchanged. Sweep: 200 pairs, δ = 10⁻⁸, 10⁻¹⁶, 10⁻²⁴, ghost escape R = 10, E = 0.3,
ε = 0.01, tmax 1500 (cloud CPU, not committed). Lost fraction 0.255 / 0.21 / 0.20 at full
precision, 0.41 / 0.32 / 0.35 with 10 bits removed, 0.555 / 0.48 / 0.505 with 20: the same
factor at every δ, and log(lost) linear in the bits removed, slope 0.039–0.046 per bit. The
losses are arithmetic. The slope implies λ_eff = κ ln 2 / slope ≈ 0.40–0.48, the paper's
ensemble-weighted separation rate (0.42), measured independently: the lost pairs are the
fast-separating tail and spend their budget at ≈ 0.42; §3's 0.3 is the rate of the certified
pairs. The QR column rising (2e-3 → 1e-2) is the frame enclosing a larger error set, the same
lever a wider centre moves the other way. Extrapolating the line to +53 bits: lost fraction
15–25 % → 2–3 % (an extrapolation over 53 bits from a 20-bit range; the rerun is the
measurement). Two checks for any implementation: rerun a sample of already-certified pairs in
the wider arithmetic (exit intervals must overlap, classification must agree; §5b's monotone
check alone cannot fail when only lost pairs are rerun), and note that the sub-step count
n = ceil(...) of a wider centre can occasionally differ, i.e. a slightly different
pair-dependent map is certified.*

### 5b. What the tier would buy, if 5a says arithmetic

With the budget of §3: 53 more bits ≈ 16 more decades of radius ≈ 120 more time units of
horizon at `λ_eff ≈ 0.3`. Pairs are lost today from `t_esc ≳ 200`; with the third word from
`t_esc ≳ 320`, which with `κ = 0.027` would cut the lost fraction from ~15 % to well under
1 % — *if* the losses are arithmetic. Since every uncertified pair in the sweep is a lost one,
and the sweep says nothing below `δ = 10⁻¹⁶` because the pairs that would be certified there
are the long-lived ones, this is what would let the certificate speak below `10⁻¹⁶`, down to
the `10⁻⁴⁰` the sweep already samples, instead of stopping at `10⁻¹⁴`. (Estimate, conditional
on 5a; the sweep after the change is the measurement.)

Check that cannot fail if the tier is right: a wider centre can only tighten, so every
certified-different count must be ≥ the 2026-10-12 value and every certified-same count ≥ it
too; `f_lo` never falls, `f_hi` never rises. A single row going the other way means a constant
is wrong.

## 6. "What would this change?" — blast radius, ignoring time

Files that change:

- `ym/hdym_ball.py` — the only file with new arithmetic: the `td` type and its four
  operations (~80 lines of C), one macro instantiation for the kick, a per-pair tier flag and
  the rerun-the-lost-pairs loop on the Python side (~30 lines), one new constant with its
  citation, one selftest (~40 lines of Python). Before that, the 5a diagnostic: a print of
  the four radius contributions at loss (~15 lines, temporary). Nothing in `ell_update`,
  `ebound`, `step_remainder`, `estatus`, `classify` or the classification output format.
- `ym/csrc/hdym_model.h` — nothing, once `hdym_ball.py` has been migrated onto it (phase 2);
  the `td` kick is an instantiation of the macro already there.
- `ym/README.md` — one `Status (date):` line on the `hdym_ball.py` row and one on the
  `results_precision/` row.
- `ym/data/results_precision/<date>/` — one new dated folder from the rerun (`ball_sweep.npz`,
  `.png`, `.txt`, `selftest.txt`), ~5 minutes at `M = 400` on the GPU.

Files that do not change: the paper, every other kernel, every committed result, `islands/`,
`bianchi9/`, `axion_su2/`, the certificate checker (`islands/code/cert_*.py`, which is exact
`Fraction` arithmetic on a finite problem and needs no tiers).

So: one file of real work, two README lines, one five-minute rerun, and an independent read of
the ~80 lines of `td` arithmetic against the cited paper before the rerun is trusted; plus,
first, a diagnostic that may say the third word is the wrong tool and the frame is the right
one, which would be a different (and smaller) change in the same file. The author is not on
anything.

## Status

*(2026-10-13)* Written; nothing implemented. *(2026-10-13, later)* §5 revised after an
independent read: the radius-triggered promotion of the first version cannot recover error
already amplified, so the tier is chosen per pair from `t = 0` (lost pairs rerun wider); and
the sweep's QR containment column means the lost pairs may be a frame problem, not an
arithmetic one, so the 5a diagnostic comes first. Order of work if taken up: 5a diagnostic on
the committed kernel; migrate `hdym_ball.py` onto `csrc/hdym_model.h` as part of phase 2
(before/after on the committed `2026-10-12` sweep, all certified fractions identical); then
whichever of the tier or the frame change 5a points to, with its selftest; then the sweep;
then the tightening check of §5b.
