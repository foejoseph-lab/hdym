# Rabbit hole 4 — the loiter: what happens between the ghost reaching order one and the blow-up

21–22 Oct 2026. Three stages: the 20 Oct section-E checkpoints (`code/rh4_loiter.py` → `results/rh4_loiter.{txt,npz,png}`,
50-unit resolution), then the kernel's crossing times (added 21 Oct: first times at which |w| reaches 0.1 and 1 and each
off-diagonal sector norm reaches 1, `tcross` (N, 5) in every section's npz) with the 20 × 64 run `islands/results/esc_E_cross.*`
(21–22 Oct night, 3162 s on 24 threads; `code/rh_overnight_2026-10-22.py` → `results/rh_overnight_2026-10-22.txt`), then the
low-energy run at tmax = 12000, `islands/results/esc_E_lowE_long.*` (20 energies in [0.05, 0.13] × 64, 22 Oct, 6559 s;
`code/rh4_lowE_long.py` → `results/rh4_lowE_long.txt`). All checks PASS, no integrator failures. The loiter is
L = t_esc − t(|w| = 1); the FP-mode caveat of `islands/results/NOTE_esc_E_kernel_2026-10-20.md` applies to the
checkpoint run and does not affect statistics.

## Result (crossing times; supersedes the checkpoint reading)

- **Short, and scaling with the period.** Pooled over 20 energies the median loiter is 18.0 time units = 0.58
  periods (minimum 6.6); per energy 0.42–0.62 periods for E ≥ 0.33 and 1.1–1.8 at 0.24–0.28. Spearman of the
  per-energy median with E is −0.05 in periods and −0.46 in time units (11 energies with ≥ 20 escaped loiters): the
  loiter is a number of orbits. (The checkpoint reading of the 6400 trajectories of `esc_E` — 4112 reached |w| = 1, 2053 of 3778 escapes inside one
  50-unit interval (no committed log: not printed in `rh4_loiter.txt`), 982 unresolved, 334 censored — had a median of 20 units = 0.63 periods and Spearman −0.09 (p 0.48)
  in periods vs −0.32 (p 0.01) in time units; the hole as posed had said "roughly E-independent in absolute terms",
  the opposite.)
- **The linear phase is the clock.** t(|w| 0.1 → 1) equals ln 10/rate_diag where the diagonal ghost is
  Floquet-unstable on the orbit (68.7 vs 65.9 at E = 0.33; 37.2 vs 38.5 at 0.39; 13.6 vs 18.2 at 0.62; 10.8 vs 11.2 at
  0.85; 13.0 vs 9.9 at 1.0) and is far faster than the ghost's own rate where an off-diagonal band leads (41.6 vs
  183 at E = 0.24): the off-diagonal mode pumps the diagonal ghost through its last decade in a fraction of the
  ghost's own e-folding time.
- **Who crosses 1 first.** The diagonal ghost, in 71 % of cells at the Krein band (E = 0.156) and 89 % at the
  parametric band (E = 0.30), even though the off-diagonal mode is the linearly unstable one: the trigger pumps the
  executor past amplitude 1 before passing 1 itself.
- **Tail.** p90/median 4.3–5.7 against 3.3 for an exponential; a Pareto tail beats an exponential by Δ log-likelihood
  2–30 with γ ≈ 1.9 (E 0.2–0.35), 8.9 (0.35–0.6, effectively exponential), 3.3 (0.6–1). Stretched, modestly, as the
  checkpoint data had said (there: Δll 3–68, γ ≈ 1.3 / 2.6 / 2.9 on the resolved part L ≥ 50, p90/median 3.6–6.3,
  Kaplan–Meier curves bending rather than straight — KAM stickiness present, not dominant — and escapes with L > 3
  periods led by an off-diagonal sector in 12 % of cases against 4–5 % for L ≤ 1 period).
- **Loiter at the two bands.** E = 0.156 median 8.3 periods (p90 27); E = 0.30 median 0.74 (p90 5.1). A 10× difference,
  but the parametric bands at E = 0.0585 and 0.110 are in the same class as 0.156, so this is the low-E regime, not a
  Krein signature.

## The low-energy regime (tmax = 12000): a slow Poisson tail, not a stuck state

At tmax = 3000 the off-diagonal bands at E = 0.0585 and 0.110 (diagonal ghost Floquet-stable) pumped the ghost to
|w| = 1 in 56 and 62 of 64 cells and 53 and 45 of them had not escaped: loiters of 16–21 periods, censored. At
tmax = 12000:

- **In the off-diagonal-led bands the escape is a Poisson process with a long time constant.** E = 0.0581: 64 cells
  reach |w| = 1 (median t = 2376), 55 escape by 12000, loiter median 3846 = 72 periods, censored-MLE τ = 106 periods;
  Kaplan–Meier survival 0.97 / 0.84 / 0.55 / 0.37 / 0.17 at L = 1000 / 2000 / 4000 / 6000 / 8000 against the
  exponential's 0.84 / 0.70 / 0.49 / 0.35 / 0.24 — exponential after a dead time of ~1000. E = 0.1063 and 0.1118: all
  64 escape, median 36–38 periods, τ = 43–45 periods, Kaplan–Meier on the exponential to ±0.1 everywhere. Against the
  diagonal-ghost tongues at E = 0.0914 and 0.1236 (loiter median 1.3–1.7 periods): the escape after the ghost reaches
  order one is memoryless in both classes, with a rate that differs by a factor 30–100 depending on which sector did
  the pumping.
- **At Floquet-stable energies between the bands a transiently pumped ghost does not run away**: at Floquet-stable energies between the bands (E = 0.050–0.087 excluding the two bands and E = 0.0526) 3–13 cells per energy reach |w| = 1 by transient growth (hole 8) and all but one sit censored at tmax = 12000 (exception, noted 28 Oct 2026: at E = 0.0526 all 64 cells reach |w| = 1 and 9 escape, τ = 875 periods, row 0.0526 of `rh4_lowE_long.txt`) — loiters ≥ 22–140 periods, medians ≥ 154–212 periods; the one escape is at E = 0.0786, after 9352 time units = 190 periods (`rh4_lowE_long.txt`). The genuinely stuck regime — an order-one ghost that sits for hundreds of orbits — is the
  linearly stable one, not the slowly unstable one.

## For the paper

One paragraph (`FOR_THE_REWRITE_2026-10-25.md`, item 7): after the ghost reaches order one the runaway is a memoryless
event whose rate is set by the linear mechanism that pumped it — of order one per orbit when the diagonal ghost is
itself unstable, of order 10⁻² per orbit when an off-diagonal band did the pumping — with a stretched tail, scaling
with the period; and at linearly stable energies a transiently pumped ghost does not run away on any time scale
reached here (≥ 200 orbits). Not: "the loiter is independent of energy".

## Checkpoint-era correction, kept for the record

The orient `.npz` files of 21 Oct (`esc_orient_0156`, `esc_orient_030`) carry no `chk` array: `run_orient` computed
it and did not pass it to `np.savez`. Fixed 21 Oct; `esc_orient_030_n64` (22 Oct) has it.
