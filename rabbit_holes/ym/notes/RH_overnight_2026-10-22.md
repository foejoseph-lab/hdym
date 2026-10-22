# The overnight runs of 21–22 Oct, first reading (holes 1, 4, 12)

Runs on Joseph's machine (24 threads), all CHECK lines PASS, no integrator failures: `esc_phase_0156` (12 phases ×
256, 2802 s after the restart with patch 13), `esc_orient_030_n64` (64², 13948 s), `esc_E_cross` (20 × 64, 3162 s).
Script `code/rh_overnight_2026-10-22.py`, log `results/rh_overnight_2026-10-22.txt`. The data files live in
`islands/results/` on Joseph's side and are not in this patch.

## Hole 1 — closed, negative

The phase ensemble (3072 independent cells, so a plain permutation null is valid, unlike the maps) shows no phase
effect once the amplitude is controlled: residual of log t_esc after log |c|, binned by the measured arg c, has an
across-bin range of z = +1.3 (12 bins over 2π) and +0.7 (6 bins mod π); the escape fraction z = +1.8 / +1.4. The raw
per-α table does show 96 % escape at α = π/2, 3π/2 against 90 % elsewhere, but that is |c|: the section's seed
direction Re(e^{iα}v) is normalised after construction, so |c| runs 0.56 → 1.21 with α (a design flaw of
`run_phase`, to be fixed by seeding at equal |c|; the regression absorbs it, slope −0.076 dex/dex, residual sd 0.084).
The maps' 4σ (RH1_phase.md) was spatial structure the tile null did not fully absorb, or map geometry correlated
with phase. Verdict: at the few-per-cent level this ensemble resolves, a Krein instability's seed phase does not
change the escape. One sentence in the paper at most, now a negative one.

## Hole 4 — a measurement

With the kernel's crossing times the loiter L = t_esc − t(|w| = 1) is resolved to the integrator step.
- Pooled over 20 energies: median 18.0 time units = 0.58 periods (min 6.6); per energy, 0.42–0.62 periods for
  E ≥ 0.33, 1.1–1.8 at 0.24–0.28.
- **Scales with the period, not absolute time**: Spearman of the per-energy median with E is −0.05 in periods and
  −0.46 in time units (11 energies with ≥ 20 escaped loiters).
- **The linear clock, verified**: t(|w| 0.1 → 1) equals ln(10)/rate_diag where the diagonal ghost is Floquet-unstable
  (68.7 vs 65.9 at E = 0.33; 37.2 vs 38.5 at 0.39; 13.6 vs 18.2 at 0.62; 10.8 vs 11.2 at 0.85; 13.0 vs 9.9 at 1.0),
  and is much faster than the ghost's own rate where an off-diagonal band leads (41.6 vs 183 at E = 0.24): the
  off-diagonal mode pumps the diagonal ghost through its last decade in a fraction of the ghost's own e-folding time.
- **Tail**: p90/median 4.3–5.7 against 3.3 for an exponential; Pareto beats exponential by Δll = 2–30 with γ ≈ 1.9
  (E 0.2–0.35), 8.9 (0.35–0.6: effectively exponential), 3.3 (0.6–1). Stretched, modest, as the checkpoint data said.
- **The low-E loophole, measured**: at E = 0.0585 and 0.110 (off-diagonal bands; diagonal ghost Floquet-stable) the
  ghost is pumped to |w| = 1 in 56 and 62 of 64 cells, and then 53 and 45 of them do **not** escape by t = 3000: loiters
  of 16–21 periods (median 847, 949 time units), censored. Order-one ghost amplitude with no runaway for the rest
  of the run — hole 8's loophole, and the "benign ghost" as a measured regime rather than a case.
- **Who crosses 1 first**: the diagonal ghost, in 71 % of cells at the Krein band (E = 0.156) and 89 % at the
  parametric band (E = 0.30), even though the off-diagonal mode is the linearly unstable one. The trigger pumps the
  executor past amplitude 1 before passing 1 itself.
- Loiter at the two bands: E = 0.156 median 8.3 periods (p90 27); E = 0.30 median 0.74 (p90 5.1). A 10× difference,
  but the parametric bands at 0.0585 and 0.110 are in the same 16–21-period class as 0.156, so this is the low-E
  regime, not a Krein signature.

## Hole 12 — still open

`chk` stores sector norms, not directions; "a quadruplet's runaway rotates while a real pair's runs along a
direction" needs state snapshots (a handful of full 36-component states per trajectory around t(|w| = 1)). Not in
these files; a kernel addition of the same size as `tcross` if wanted.

## For the paper

Hole 4's bullets are the paragraph: loiter ≈ 0.6 orbits median, scaling with the period, stretched tail, the
linear phase as the clock, the off-diagonal trigger pumping the diagonal executor, and a low-energy regime where
an order-one ghost sits for tens of orbits without escaping. Hole 1 becomes one negative sentence.

## The day run (22 Oct): the low-E regime at tmax = 12000 — a slow tail, not a stuck state

`esc_E_lowE_long` (20 energies in [0.05, 0.13] × 64, tmax 12000, all checks PASS, 6559 s), read by
`code/rh4_lowE_long.py` → `results/rh4_lowE_long.txt`.

- **In the off-diagonal-led bands the loiter is a Poisson process with a long time constant.** E = 0.0581: 64 cells
  reach |w| = 1 (median t = 2376), 55 escape by 12000, loiter median 3846 = 72 periods, censored-MLE τ = 106 periods;
  Kaplan–Meier survival 0.97 / 0.84 / 0.55 / 0.37 / 0.17 at L = 1000 / 2000 / 4000 / 6000 / 8000 against the exponential's
  0.84 / 0.70 / 0.49 / 0.35 / 0.24 — exponential after a dead time of ~1000. E = 0.1063 and 0.1118: all 64 escape,
  median 36–38 periods, τ = 43–45 periods, Kaplan–Meier on the exponential to ±0.1 everywhere. Against the
  diagonal-ghost tongues at E = 0.0914 and 0.1236 (loiter median 1.3–1.7 periods): the escape after the ghost reaches
  order one is memoryless in both classes, with a rate that differs by a factor 30–100 depending on which sector did
  the pumping. Last night's "16–21-period censored loiters" were the first 3000 units of these.
- **At Floquet-stable energies (between the bands) a few cells reach |w| = 1 by transient growth (hole 8) and then
  none escape in 12000**: 3–13 cells per energy, loiters ≥ 55–140 periods (median of the censored ≥ 150–230 periods).
  That is the genuinely stuck regime — an order-one ghost that sits for hundreds of orbits — and it is the
  linearly-stable one, not the slowly-unstable one.
- The paper sentence this supports: after the ghost reaches order one the runaway is a memoryless event whose rate
  is set by the linear mechanism that pumped it — of order one per orbit when the diagonal ghost is itself unstable,
  of order 10⁻² per orbit when an off-diagonal band did the pumping — and at linearly stable energies a transiently
  pumped ghost does not run away on any time scale reached here (≥ 200 orbits).
