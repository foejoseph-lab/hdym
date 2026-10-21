# Rabbit hole 1 — phase-dependent loiter at a Krein band, redone from the committed maps, with the control

Session of 21 Oct 2026 (Claude, cloud). Script `code/rh1_phase.py` (two monodromies, ~1 min with the nulls), log
`results/rh1_phase.txt`, figure `results/rh1_phase.png`. Input: `islands/results/esc_orient_0156.npz` (E = 0.156,
block 12 Krein quadruplet, 128², 88.7 % escape by 3000, median 1952) and `esc_orient_030.npz` (E = 0.30, block 23
real pair, 100 % escape, median 631). No runs. Replaces the session scratch `phase.py` / `nodal2.py`.

## The idea, stated so it can be tested

A Krein quadruplet's unstable subspace is two-dimensional and the seed's projection on it is one complex number
c = w·s. Its modulus sets the linear phase (the nodal structure already in `offdiag_escape_nodal.py`); its
argument is the phase of the growing mode and so decides *at which point of the orbit the ghost reaches order
one*. For a real pair c is real: the seed can choose a sign, nothing more. So a dependence of the escape on
arg c is a nonlinear signature available to a Krein instability and not to a parametric one. The hole's claim
was a 6–15 % swing in the non-escape fraction and 1870–2090 in the median with arg c, and a 0.57 correlation of
the map with its antipode against a null of 0.00 ± 0.035.

## What the committed maps give

Krein map, arg c mod π (s → −s identified; 6 wedges of ~2700 cells):
- escape fraction by wedge 86.8–92.9 % (non-escape 7–13 %; the list's 6–15 % at 12 bins is reproduced);
- median t_esc by wedge 1900–2048 (the list's 1870–2090 reproduced);
- after removing the linear predictor log|c| (which explains little here: corr −0.12, residual sd 0.085 dex),
  the median residual still swings by 0.039 dex (9 %) across wedges.

Against a null that **permutes 8 × 8 tiles of the map** (keeping the loiter texture's spatial correlation —
phase wedges are contiguous regions of the sphere, so a cell-wise null overstates the significance, which is
what the list's "0.00 ± 0.035" did): z = +3.5 (escape fraction), +3.9 (median), +5.6 (residual) for the mod-π
binning; z ≈ +2 for the 12-bin [0, 2π) binning, which halves the cells per bin and whose extra structure (the
s → −s partner) is not there.

Control map (real pair), the same machinery: sign c gives median 634 vs 627 (z ≈ 0 against the tile null);
a placebo phase — the argument of the seed's projection on a *stable* complex pair of the same block — gives
z = −0.2 … +0.6 on all three statistics. The method does not manufacture a phase effect.

## Verdict

The effect is real at the ≈ 4σ level and **small**: a 4–7 % modulation of the escape fraction and of the
median time, period π in arg c, against a loiter scatter of 0.085 dex (22 %) per cell. It is what the hole said
it was — a nonlinear distinction between the two mechanisms that a parametric band cannot show — but on one
map at one energy it is a signature, not a measurement. The list's price (an ensemble per phase bin, one night
on 20 cores) is the right next step if it is wanted: 6 phase wedges × ~500 seeds × 2 diagonal-seed draws would
turn z ≈ 4 into a curve and would also answer whether the phase effect scales with the Krein rate (0.0062 here)
or with the splitting of the two pairs.

For the paper: not now. One sentence in the fractal/orient discussion at most ("the Krein map carries a weak
dependence on the phase of the seed's projection, period π, absent in the parametric map"), with the numbers
above. Standalone paper: no; as a section of the hole-2/7 paper after the ensemble run: possibly, as the
nonlinear counterpart of the typed classification.

## Sandbox job (if wanted)

```
# offdiag_escape.py needs a --section phase: seeds s = delta (cos a e_r + sin a e_i) with e_r, e_i the real and
# imaginary parts of the Krein right eigenvector direction (so |c| = const and arg c = a), a on a 6- or 12-point
# grid, nseed random orthogonal admixtures per a, both diagonal seeds; ~6 x 500 x 2 trajectories at E = 0.156,
# ~1.2 s each -> ~2 h on 24 threads.  Not written.
```
