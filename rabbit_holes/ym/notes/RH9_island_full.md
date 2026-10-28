# Rabbit hole 9 — the E = 0.3 island in the full nine-component system (23 Oct 2026)

Script `islands/code/island_full.py` (22 Oct): `hdym_island.py`'s ball sampling around the DR point — 8 shells
r = 10⁻⁶ … 10⁻², 200 points each, plus a sea control — run through the nine-component kernel
`ym/csrc/hdym_full_escape.c` instead of the 12-d diagonal-sector integrator, with the six off-diagonal components
seeded at δ_od (Gauss-projected) and the diagonal ghost at |w| = √(2εE). Conversion of the 12-d state to the
36-component one checked to |E_full − E_12|/E ≤ 10⁻¹⁴ and |G| ≤ 4 × 10⁻¹³ on every shell. Runs of 23 Oct, 24 threads,
E = 0.3, ε = 10⁻⁸, tmax 5000, rtol 10⁻¹¹, no integrator failures: `islands/results/island_full_E0.3_eps1e-8.{txt,npz}`
(δ_od = 10⁻⁸, three decades under the diagonal ghost's 7.7 × 10⁻⁵; 2999 s) and `_od0.{txt,npz}` (δ_od = 0, the
exactly-diagonal control; 632 s). Queue: `code/hole9_queue.md`.

## The question

`islands/README.md` had the E = 0.3 island as an intact hole for the diagonal ghost at ε = 10⁻¹⁰–10⁻⁸ (hole 10) and
"not a refuge in the full theory" from the linear off-diagonal multiplier only: block (u23, u32) has the real
multiplier 1.5522 on the DR orbit at E = 0.3 (rate 0.0125) while the diagonal ghost is Floquet-stable there. Whether
the island's interior — the tori around the orbit, not just its centre — empties once all nine components are
integrated, through which sector, and when, had not been measured.

## Result

| shell r | δ_od = 10⁻⁸: escaped | t_esc median [p10, p90] | escaping sector | δ_od = 0: escaped | t_esc median |
|---|---|---|---|---|---|
| 10⁻⁶ | 200/200 | 1340 [1249, 1557] | diagonal 191, block 23: 9 | 0/200 | — |
| 10⁻⁵ | 200/200 | 1358 [1243, 1595] | diagonal 182, block 23: 18 | 0/200 | — |
| 5 × 10⁻⁵ | 200/200 | 1359 [1228, 1550] | diagonal 184, block 23: 16 | 0/200 | — |
| 10⁻⁴ | 200/200 | 1346 [1240, 1540] | diagonal 184, block 23: 16 | 0/200 | — |
| 2 × 10⁻⁴ | 200/200 | 1369 [1242, 1664] | diagonal 185, block 23: 15 | 10/200 (5 %) | 962 |
| 5 × 10⁻⁴ | 199/200 | 1276 [696, 1575] | diagonal 190, block 23: 9 | 96/200 (48 %) | 979 |
| 10⁻³ | 198/200 | 914 [511, 1463] | diagonal 193, block 23: 5 | 155/200 (78 %) | 825 |
| 10⁻² | 198/200 | 567 [384, 964] | diagonal 198 | 200/200 | 560 |
| sea control | 200/200 | 499 [382, 737] | diagonal 200 | 200/200 | 512 |

- **The control reproduces `hdym_island`.** With no off-diagonal seed the core shells r ≤ 10⁻⁴ do not escape in 5000,
  the r = 2 × 10⁻⁴ shell escapes 5 %, the halo 48–78 %, the sea 100 % — the hole-10 table at ε = 10⁻⁸, line for
  line. The nine-component kernel on an exactly diagonal state is the diagonal theory.
- **With the off-diagonal sector seeded at 10⁻⁸ the whole core empties**, every shell, 200 of 200, with a median
  escape time of 1340–1370 and p10 ≈ 1230–1250 — a narrow packet, the same at r = 10⁻⁶ and at 10⁻⁴: the time is set by
  the off-diagonal growth from the seed, not by where in the island the trajectory starts. The queue's estimate from
  the block-23 rate, ln(10⁸)/0.0125 ≈ 1500, is the right scale.
- **Trigger and executor, as everywhere else on this orbit.** The diagonal ghost is the escaping sector in 91–96 % of
  the core cells and block 23 in the rest; |w| reaches 1 at a median t ≈ 1313–1342 and a block norm reaches 1 at
  ≈ 1321–1352 — the diagonal ghost, pumped nonlinearly, passes amplitude 1 before the off-diagonal mode that drives
  it does, as in the loiter measurement (`RH4_loiter.md`: 71–89 % of cells). The halo shells, which escape on their
  own at 48–78 % in the control, escape at 99–100 % with the seed, and faster.

## Reading

The E = 0.3 island is not a refuge in the full homogeneous theory, and not only at its centre: at ε = 10⁻⁸, where the
diagonal ghost alone leaves it intact, an off-diagonal admixture three decades smaller empties its whole core in
≈ 38 orbits (1340/35.28, median t_esc over T at E = 0.3; stated as 25 until 28 Oct 2026), through the block-23 parametric instability of the orbit, with the diagonal ghost doing the escaping.
Together with hole 10 this closes the island question for this background on both sides: the diagonal ghost's
back-reaction erases it above ε ≈ 5 × 10⁻⁸, and below that the off-diagonal sector does. Not measured: whether an
off-diagonal seed small enough (the integrator's own truncation, ~10⁻¹⁰, is the floor the control sits at) or a
different energy (E = 1, where the diagonal ghost is itself unstable, is the paper's other island energy) changes the
time scale; the mechanism would not change.

## For the paper

Watch-list item 3 (`FOR_THE_REWRITE_2026-10-25.md`): the open item (iii) on the regular islands can now say that the
E = 0.3 island is emptied by the off-diagonal sector throughout its core at ε = 10⁻⁸ (this note) and by the diagonal
ghost's back-reaction above ε ≈ 5 × 10⁻⁸ (hole 10); the positive-measure question for the full theory remains open
only as a statement about other islands. Hole closed.
