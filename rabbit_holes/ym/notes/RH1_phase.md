# Rabbit hole 1 — does the seed's phase on a Krein quadruplet change the escape? (closed, negative)

21–22 Oct 2026. Scripts `code/rh1_phase.py` (the committed maps; two monodromies, ~1 min with the nulls) and
`code/rh_overnight_2026-10-22.py` (the phase ensemble). Logs `results/rh1_phase.txt`, `results/rh_overnight_2026-10-22.txt`;
figure `results/rh1_phase.png`. Data: `islands/results/esc_orient_0156.npz` (E = 0.156, block 12 Krein quadruplet,
128², 88.7 % escape by 3000, median 1952; Krein rate 0.0062 there), `esc_orient_030.npz` (E = 0.30, block 23 real pair, 100 % escape,
median 631), `esc_phase_0156.*` (the ensemble, 21–22 Oct night, 12 phases × 256 cells, 2802 s on 24 threads).

## The idea

A Krein quadruplet's unstable subspace is two-dimensional and the seed's projection on it is one complex number
c = w·s. Its modulus sets the linear phase (the nodal structure in `islands/code/offdiag_escape_nodal.py`); its
argument is the phase of the growing mode and decides at which point of the orbit the ghost reaches order one. For
a real pair c is real: the seed can choose a sign, nothing more. A dependence of the escape on arg c would be a
nonlinear signature available to a Krein instability and not to a parametric one. As posed on 21 Oct, from
scratch scripts: a 6–15 % swing of the non-escape fraction and 1870–2090 of the median with arg c, and a 0.57
correlation of the smoothed map with its antipode against a null of 0.00 ± 0.035.

## From the committed maps (21 Oct)

Krein map, arg c mod π (s → −s identified; 6 wedges of ~2700 cells): escape fraction by wedge 86.8–92.9 %
(non-escape 7–13 %), median t_esc 1900–2048, and after removing the linear predictor log|c| (corr −0.12, residual
sd 0.085 dex, i.e. 22 % loiter scatter per cell) the median residual still swings 0.039 dex (9 %) across wedges. Against a null that permutes 8 × 8
tiles of the map — phase wedges are contiguous regions of the sphere, so a cell-wise null overstates the
significance, which is what the 0.00 ± 0.035 did — z = +3.5 (escape fraction), +3.9 (median), +5.6 (residual) for
the mod-π binning, z ≈ +2 for 12 bins over [0, 2π). Controls: the real-pair map gives z ≈ 0 for sign c (median
634 vs 627); a placebo phase (the seed's projection on a stable complex pair of the same block) gives
z = −0.2 … +0.6 on all three statistics. So the method does not manufacture an effect, and on one map at one
energy the effect was a signature at ≈ 4σ, not a measurement.

## The ensemble (22 Oct): no effect

3072 independent cells (so a plain permutation null is valid), seeded on the quadruplet with a random admixture and
diagonal seed per cell, intended at fixed |c| (which it was not, see below). Residual of log t_esc after log|c|, binned by the measured arg c:
across-bin range z = +1.3 (12 bins over 2π) and +0.7 (6 bins mod π); escape fraction z = +1.8 / +1.4. The raw
per-phase table does show 96 % escape at α = π/2, 3π/2 against 90 % elsewhere, but that is |c|: the section's seed
direction Re(e^{iα}v) is normalised after construction, so |c| runs 0.56 → 1.21 with α (a flaw of `run_phase`;
seeding at equal |c| would remove it; the regression absorbs it, slope −0.076 dex/dex, residual sd 0.084).
The maps' 4σ was spatial structure the tile null did not fully absorb, or map geometry correlated with phase.

**Verdict.** At the few-per-cent level this ensemble resolves, a Krein instability's seed phase does not change
the escape. For the paper: one negative sentence at most (`FOR_THE_REWRITE_2026-10-25.md`, item 7).
