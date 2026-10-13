# Arbitrary-precision ball arithmetic as a library: the third word without writing it

2026-10-13, late. Companion to `NOTE_adaptive_precision.md`, whose §5 describes a hand-written
triple-double tier for the lost pairs of the ball sweep. This note records that the tier does
not have to be written: it can be imported. Nothing here changes any code or result.

## What the scheme was describing

`NOTE_adaptive_precision.md` §1 gives the author's scheme: float as far as it holds, bound it,
raise the precision where the bound says so, repeat. That is a description of **Arb**
(F. Johansson, 2013–; since 2023 part of FLINT 3; Python binding `python-flint`): a number is a
midpoint at a chosen precision plus a radius that is a proven upper bound on everything the
midpoint has lost, every operation propagates both, and the precision is a parameter you set
per call. Its documented way of working is the scheme's step 4 verbatim: compute at some
precision, inspect the radius, and if it is too wide, double the precision and compute again.

So the arithmetic half of the scheme is a standard library with a decade of use. What Arb does
**not** supply is the enclosure strategy for a chaotic flow: integrating the 12-d system in Arb
balls naively, the wrapping effect eats the radius within a few time units at any precision,
exactly as in doubles. The ellipsoid (Lohner) frame of `hdym_ball.py` — `ell_update`,
`ebound`, `step_remainder`, the QR containment — is the part that has to stay, and it is the
part that is already reviewed.

## The lost pairs as a CPU rerun in Arb

The 2026-10-13 degradation test (`NOTE_adaptive_precision.md` §5a Status) established that the
lost pairs of the sweep (9–25 % of 400) die of arithmetic, so a wider centre would recover
most of them. Instead of the hand-written `td` type of §5:

- rerun only the lost pairs, on the CPU, with the **centre arithmetic** of the member and
  difference balls done in `python-flint` `arb` at a chosen precision (200 bits is more than
  the 159 of triple-double; the cost of bits in Arb is mild);
- keep the **frame** exactly as it is: the tangent balls, the step remainder, `ell_update`
  and the containment radius consume centres as doubles (`arb` midpoint rounded to double, an
  exact operation to read, with the rounding charged to the radius) and radii as doubles;
- the per-operation charge `DDU` and the `INFL` / `ETA` inflation of the centre arithmetic
  disappear: Arb's radius already contains the rounding of every operation it performs. The
  `U2` charges of the tangent balls and the frame constants are untouched;
- escape bookkeeping, classification and output format unchanged, so the rerun's rows drop
  into the same `ball_sweep` format and the monotone check of §5b applies row by row.

Cost: ~60 pairs × 25 000 steps × a few hundred operations at 200 bits, in Python calling
Arb. Minutes to an hour on one core, once, as a check. Not a production path, and not meant
to be: if it tightens the way §5b predicts, that is the result, and the fast `td` kernel is
then an optimisation one can choose to build or not.

What it removes from the plan: the ~80 lines of new triple-word arithmetic and their error
constants, the Fabiano–Muller–Picot citation, and most of the review burden. The remaining
review item for this path is the interface: that the centre read out of an `arb` into the
double-valued frame is charged correctly (one rounding, `2⁻⁵²` relative, to the radius) and
that the frame's inputs are the ball's, not the midpoint's. A dozen lines.

What it does not remove: the Muller–Rideau (ACM TOMS 2022) check on the double-word constants
of the *existing* kernel, which certifies everything the sweep currently claims and is
unaffected by anything in this note.

## Precision as a dial, not a tier

With Arb the question "double-double or triple-double" dissolves. The budget of
`NOTE_adaptive_precision.md` §3 — bits available minus bits that must survive equals
`λ_eff t / ln 2`, at `λ_eff ≈ 0.42` for the lost pairs — becomes a formula for the precision
to request per pair: `53 + 0.6 · t_horizon` bits, plus margin. A pair that escapes at
`t = 320` wants ~250 bits; one that escapes at `t = 500` wants ~350. Arb takes either. The
step-4 loop "if lost, double the precision and rerun" is then literally what Johansson's
documentation recommends, applied to pairs.

## Status

*(2026-10-13)* Written; nothing implemented. If taken up, it goes after the phase-2 migration
of `hdym_ball.py` onto `csrc/hdym_model.h` (so the C and Python paths share the kick's text by
construction) and replaces §5's `td` tier as the first thing to try. Dependencies to add to
`ym/README.md` Requirements if so: `python-flint` (wheels for Python 3.12 and 3.13 on
Windows, Linux and macOS).
