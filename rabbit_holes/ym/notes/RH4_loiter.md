# Rabbit hole 4 — the loiter's statistics, from the section-E checkpoints

Session of 21 Oct 2026 (Claude, cloud). Script `code/rh4_loiter.py` (reads `islands/results/esc_E.npz`; no
integration), log `results/rh4_loiter.txt`, data `results/rh4_loiter.npz`, figure `results/rh4_loiter.png`.
Scope: the 6400 trajectories of the 20 Oct section-E run (100 energies × 64 seeds, δ = δ_diag = 10⁻⁴, tmax 3000,
checkpoints every 50 time units; the FP-mode caveat of `NOTE_esc_E_kernel_2026-10-20.md` applies and does not
affect statistics).

First, a correction to the list. The entry says checkpoint data exist for "6400 + 2 × 16384 trajectories". The
orient `.npz` files carry no `chk` array: `offdiag_escape.py`'s `run_orient` computes it and does not pass it to
`np.savez` (lines 465–466). Holes 1 and 12, which the list says share "the same data", have no checkpoint data
in the repository; only the section-E run does. (The orient maps do carry `tesc`, `sector`, `diag` and `basis`.)

## Definitions

t₁ = the time the diagonal ghost first reaches |w| = 1, interpolated log-linearly between checkpoints; when the
ghost crosses 1 and blows up inside one checkpoint interval (the common case: 2053 of 3778 escapes), t₁ is
extrapolated from the last two checkpoints with their local exponential rate, and if that lands after the
escape the loiter is clamped to 0 (982 trajectories, "unresolved"). Loiter L = t_esc − t₁. A trajectory that
reaches 1 and does not escape by tmax is a right-censored loiter (334, almost all below E = 0.2, where t₁ ≈
2300 leaves < 700 units to escape). 4112 trajectories reached |w| = 1 (the pooled fits count 4110: the lowest energy, 0.04999…, falls outside the bin edge).

## What the data say

1. **The loiter is short.** Median L = 20 time units (0.63 periods) over the 61 energies with ≥ 20 escaped
   loiters; a quarter of all escapes blow up within the checkpoint interval in which the ghost reached 1. The
   list's "median 10–170" is the upper part of this; at 50-unit resolution the median is not resolved.
2. **It scales with the period, not with absolute time** (the opposite of the list's "roughly E-independent
   in absolute terms"): Spearman of the per-energy median with E is −0.09 (p 0.48) in periods and −0.32
   (p 0.01) in time units. The loiter is a number of orbits.
3. **The tail is heavier than exponential, modestly.** On the resolved part (L ≥ 50 time units, or ≥ 2 periods;
   13–63 % of the trajectories depending on E) a Pareto tail beats an exponential by Δ log-likelihood 3–68 for
   one extra parameter in every E bin, with tail index γ ≈ 1.3 (E = 0.2–0.35), 2.6 (0.35–0.6), 2.9 (0.6–1) in
   time units and p90/median of the excess 3.6–6.3 against the exponential's 3.3. Below E = 0.2 the tail is
   censoring-dominated (322 of 616 censored) and the fitted γ = 0.27 is a tmax artefact. The Kaplan–Meier
   curves (figure, log-log panel) bend down rather than run straight, so this is "a stretched tail", not a
   clean power law: KAM stickiness is present but not dominant at these E.
4. **Long loiters end differently.** Escapes with L > 3 periods are led by an off-diagonal sector in 12 % of
   cases, against 4–5 % for L ≤ 1 period (`== 4.` in the log): the longer the ghost wanders, the more often the
   blow-up goes through block 12 or 23 rather than the diagonal ghost.
5. Where the diagonal ghost is Floquet-stable on the orbit but an off-diagonal band is open (E = 0.058,
   0.11, 0.13), the ghost is pumped to |w| = 1 by t ≈ 1600–2700 and then does **not** escape by 3000 in
   most cases (E = 0.058: 46 reached 1, 2 escaped). That is hole 8's loophole seen from the other side: order-one
   ghost amplitude with no escape over the remaining 300–1400 time units; tmax decides what is "benign".

## What this needs to become a measurement

The 50-unit checkpoint spacing is the limit. The right fix is three numbers per trajectory from the kernel:
the first times at which |w|, and each sector norm, cross 1 (and perhaps 0.1 and 10), recorded in
`one_trajectory` where `sector_norm` and `ghost_amp` are already evaluated at every accepted step — no extra
cost, no larger arrays. With that, the loiter distribution is resolved to the step size and the survival fit
is honest down to L = 0. A rerun of section E at 20 energies × 64 seeds would take ~1 h on 24 threads
(16893 s for 100 × 64); the 2 orient maps are 5.5 and 10.7 h each and are not needed for this.

Sandbox job (Joseph's machine; nothing here was run at size):
```
# after adding the crossing times to ym/csrc/hdym_full_escape.c and offdiag_escape.py's npz
python islands/code/offdiag_escape.py --section E --nE 20 --nseed 64 --tag esc_E_cross > islands/results/esc_E_cross.txt 2>&1
```

## Paper-worthiness

Not a paper. A paragraph for the fractal section if the paper wants it: the loiter is 0–3 orbits for 90 % of
escapes, scales with the period, and has a stretched tail of a few percent of trajectories that wander ≥ 10
orbits, more often ending through the off-diagonal sector. With the crossing times recorded it would be one
figure (survival in periods, by E bin) and could carry the "benign ghost as a distribution" sentence the list
suggests.
