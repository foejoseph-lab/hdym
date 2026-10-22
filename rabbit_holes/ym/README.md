# rabbit_holes/ym — the Yang–Mills spin-1 ghost holes, worked

Session of 21 Oct 2026 (Claude, cloud; Joseph's rules: nothing over 20 min run here, Joseph pushes). The list
being worked is `RABBITHOLES_ym_2026-10-21.md` (copied in verbatim). Layout mirrors the repo: `code/`,
`results/`, `notes/`. Everything here reads committed data from `islands/` and `ym/`; the three scans that
were run cost 41 s, 184 s and nothing (log parsing). Nothing here is paper work and nothing in `islands/` or
`ym/` was touched.

## Done tonight

| # | hole | outcome | where |
|---|---|---|---|
| 7 | Derivation of the selection rule | **Two of the three rules are now exact**: from M = N² with N real symplectic, a real negative multiplier of M off the circle is never isolated, so a single pair cannot fold at −1 and every single-pair instability is a +1 fold (two lines). Parity of the coupling verified to 10⁻¹¹ (self even / cross odd harmonics; swap symmetry in the ± sector basis). "Collisions only at −1" is exact in the decoupled symmetric limit (where Krein bands are the sector's period-doubling tongues) but the coupling is O(1) and that model sits a rung off; the real rule is "collisions at N = ±(0.44–0.68) π". Hole 5's "why block 12" falls out in one sentence (only block 12's components mix under the reflection). | `notes/RH7_selection_rule.md`, `notes/RH7_literature.md`, `code/rh7_selection_rule.py`, `results/rh7_*` |
| 2 | Resonance ladder | **A formula**: the block-12 band centres lie on T_c = 4.057 (n + 0.210), n = 7–13, odd n parametric, even n Krein, rms 0.11 time units; the ghost accumulates phase at a constant effective frequency 0.774 (23 % below the Lee–Wick mass) over a factor 12 in E. Out-of-sample: rung 14 (Krein) predicted at E = 0.0421, **found at 0.0427–0.0435**. Above E ≈ 1 the ladder breaks (the two pairs are a rung apart; rung 6 is parametric, not Krein). A canonical rotation number was tried and does not come out of the u12 phase. | `notes/RH2_ladder.md`, `code/rh2_ladder.py`, `code/rh2_predict.py`, `results/rh2_*` |
| 2+7 | First-order averaging (the sympy reduction) | **Reviewed and cut down** (see the Review section of the note): the exact parts stand (canonical form, normal form); the "mode-0 rungs predict the bands" claim fails a half-rung null (bands + 8 % cover 79 % of log E) and rung 11; the c_j family shows the c_j = 1 bands are the remnant of a strong cancellation (c_j = 0.25 is far more unstable); the gauge-reduced 4-dim block is the zeroth order to build next. Before review the claim read: the canonical Ostrogradsky form of block 12 (sympy, checked to 1e-16), the two-mode resonance normal form (sympy: sum resonances of signed frequencies unstable, half-width 2|g|), and the time-averaged block's ghost mode ω₀(E) = 1.05–1.33 reproduce all eight block-12 bands at ω₀T = nπ within 8 %, even n parametric, odd n Krein (period doubling of a single mode, forbidden by the lemma, realised as the collision). Widths and the odd-rung mechanism need second order — one sympy sitting, specified in the note. | `notes/RH27_averaging.md`, `code/rh27_averaging.py`, `results/rh27_*` |
| 2+7 (cont.) | The gauge-reduced 4-dim ghost block | **Built and exact**: Gauss row + moving gauge (constant gauge fixing is inadmissible on this orbit), then the colour-rotation symmetry and its charge; the 4-dim reduced monodromy reproduces the four ghost multipliers to 1e-11 at every E and c_j. Two modes: a light one soft at the turning points (the quasi-static line, and the transient hyperbolicity behind the low-c_j static instability) and a heavy one at 1–1.7 (the line near 1). First-order averaging of the exact block is still not Floquet-accurate. | `notes/RH_reduced_block.md`, `code/rh_reduced_block.py`, `results/rh_reduced_block_*` |
| 5 | δ(E), the splitting of the two pairs | Measured from B005 (`rh5_splitting.py`): Δ = (φ_A + φ_B) mod 2π is 0.07–0.25 π for E = 0.05–0.35 (near-degenerate pairs, the decoupled-symmetric structure), 0.5–0.8 π at 0.5–0.72, exactly 0 at E = 0.80 (the fake degeneracy), 0.95 π at E = 1 (antipodal). The ladder works where Δ is small. | `code/rh5_splitting.py`, `results/rh5_splitting.*`; note in RH27's Review |
| 1 | Phase-dependent loiter at a Krein band | **Closed, negative** (overnight ensemble, 3072 independent cells: residual phase effect z = +0.7–1.3 after controlling |c|; `notes/RH_overnight_2026-10-22.md`). Earlier, from the maps: reproduced from the committed maps in a committed script, with the control the list only asserted: the Krein map's escape fraction and median time swing 4–7 % with arg(w·s), period π, at ≈ 4σ against a tile-permuting null (the list's cell-wise null overstated it); the parametric map and a placebo phase give z ≈ 0. Real, small, a signature not a measurement; the ensemble-per-phase run stays the price. | `notes/RH1_phase.md`, `code/rh1_phase.py`, `results/rh1_*` |
| 6 | Which polarisation carries which sign | **Negative, with the reason**: no polarisation basis exists in the repo, and "the ghost's pair" needs the frequency sign, which is the rotation number hole 2 could not get. Measured instead: at the turning point both pairs are an exact 50/50 mix of the reflection sectors (symmetry), the component/sector/heavy-light contents reshuffle with E and along the orbit, no E-independent assignment. Reading: the two pairs are the same sector ghost twice, split by a phase δ(E) — tabulate δ, not a polarisation. | `notes/RH6_polarisation.md`, `code/rh6_polarisation.py`, `results/rh6_*` |
| 4 | Loiter statistics | **Now a measurement** (overnight section E with crossing times: median 0.58 periods, scales with the period, linear clock verified, tail stretched, and a low-E regime of 16–21-period loiters without escape; `notes/RH_overnight_2026-10-22.md`). Earlier, from checkpoints: median loiter 20 time units (0.6 periods); it scales with the period, not absolute time (the list says the opposite); a stretched, modestly heavy tail (Pareto beats exponential in every E bin, γ ≈ 1.3–2.9, not a clean power law); long loiters end through the off-diagonal sector 3× more often. The 50-unit checkpoints are the limit: a quarter of escapes are unresolved. **Correction to the list**: the orient `.npz` files have no `chk` arrays (`run_orient` never saves them), so holes 1 and 12 have no checkpoint data in the repo. | `notes/RH4_loiter.md`, `code/rh4_loiter.py`, `results/rh4_*` |

## Verdicts on "interesting in their own right"

(Books checked, `notes/RH7_literature.md`: the lemma is the Swift–Wiesenfeld 1984 mechanism; ours is the
symplectic multiplicity version and its identification with the Krein collision — cite, then state.)

- **Holes 2 + 7 together are a short paper's worth of material** (and after `RH27_averaging.md` they have their mechanism), or a strong section: a half-period
  symmetry forbids the ghost's period-doubling for a single pair and forces every −1 instability to be a
  two-pair Krein event (exact, general for any reflection-symmetric periodic background in a higher-derivative
  theory), and on the DR orbit the band centres follow a one-parameter ladder with an E-independent effective
  ghost frequency and a confirmed out-of-sample prediction. Title-shaped: "Krein collisions as dressed
  period doubling: a selection rule for Lee–Wick ghosts on symmetric periodic orbits". What it would still
  need: a two-harmonic
  Mathieu-type reduction that predicts the rung positions (the thing that would make it more than a
  measurement), and the same ladder on one more orbit (the x ↔ y-rotated copies are the same orbit; a
  non-planar one is hole 5's open question).
- Hole 4 is a paragraph, not a paper, and needs the kernel to record crossing times before it is a measurement.
- The rest of the list, in this light: 12 needs an orient run that saves `chk` (or crossing times) — 5.5–10.7 h,
  sandbox; 1 needs its ensemble-per-phase run (~2 h, `RH1_phase.md`); 3 needs the positive norm (half a day of thought, no runs); 8 is an hour on 20 cores; 9 a day of
  porting; 10 a few `hdym_island.py` runs; 11 an hour. None of them was touched.

## Sandbox jobs this session asks for (not run here; all over 20 min or on Joseph's kernel)

1. Hole 4 proper: add first-crossing times of |w| = 1 and of each sector norm = 1 to `one_trajectory` in
   `ym/csrc/hdym_full_escape.c` and to `offdiag_escape.py`'s npz; rerun section E at 20 energies (~1 h / 24 threads).
2. Holes 1 and 12: an orient rerun that saves `chk` (one line in `run_orient`) or the crossing times — 128² at
   E = 0.156 is 5.5 h.

## Smaller things found on the way

- `offdiag_escape.py` line 465: `chk` is computed per orient row and dropped; the list's "2 × 16384" checkpoint
  sets do not exist.
- B005's member labels (k0..k3) are not meaningful across the wide Krein band 0.28–0.49 (all four leave the
  circle; the tracker re-pairs on re-entry), so "which pair folded at 0.5694" is not a question the log can answer.
- The frozen block-12 K at the turning point (E = 0.3) has a real pair ±1.32 besides the two imaginary ones:
  the frozen picture is hyperbolic at the bounce, which bears on hole 3 (the growth is not an adiabatic
  frequency effect there).
