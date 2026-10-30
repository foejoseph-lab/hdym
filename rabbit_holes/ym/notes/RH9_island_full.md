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
time scale; the mechanism would not change. (29 Oct 2026: E = 1 measured, next section; the small-seed question stays
open.)

## A second seed, and E = 1 (28–29 Oct 2026)

**Seed 8 at E = 0.3, ε = 10⁻⁸** (`islands/results/island_full_E0.3_eps1e-8_s8.{txt,npz}`; same command as seed 7 with
`--seed 8`, 4281 s on 6 threads). It reproduces seed 7: every core shell r = 10⁻⁶ … 2 × 10⁻⁴ escapes 200 of 200, median
t_esc 1356–1369 (seed 7: 1340–1369), the diagonal ghost the escaping sector in 181–189 of 200 cells (seed 7: 182–191),
block 23 in the rest (one cell block 12); halo 196–200 of 200; sea 200/200, median 511 (seed 7: 499). The packet does
not depend on the sample.

**E = 1, ε = 10⁻¹⁰, off-diagonal seed 10⁻⁸** (`islands/results/island_full_E1_eps1e-10.{txt,npz}`, seed 7, 3625 s):

- Every shell escapes 200 of 200, median t_esc 76–86 in the island (r = 10⁻⁶ … 10⁻²), against 137 in the sea. The
  diagonal ghost is the escaping sector in every cell. \|w\| reaches 1 at a median 61–65; a block norm reaches 1 in only
  15–26 of 200 cells per shell, after \|w\| (median 70–84).
- So at E = 1 the off-diagonal sector plays no part: the diagonal ghost's own Floquet instability on the orbit empties
  the island, faster than the sea. The scale fits a linear estimate: from √(2εE) = 1.41 × 10⁻⁵ at the smoke run's rate
  0.233 (islands README; no committed log), ln(1/1.41 × 10⁻⁵)/0.233 ≈ 48 to \|w\| = 1, against 61–65 measured [est].
- No `--delta-od 0` control was run at E = 1. With every escape diagonal it would be expected to match [est]; it is
  cheap if wanted.
- **The 12-d run agrees, and corrects the islands README.** The README described the 12-d `hdym_island.py` run at
  E = 1, ε = 10⁻¹⁰ as a fast plateau "with a sticky halo that escapes slower than the sea", from a log that was not
  committed. Rerun 29 Oct (`islands/results/island_E1_eps1e-10.log`, `islands/results/island_E1_eps1e-10_c0_N200.json`,
  `--N 200 --backend c`): every shell escapes 200 of 200, median t_esc 74.7–77.4 for r ≤ 10⁻³ and 84.9 at r = 10⁻²,
  against 136.1 in the sea; core FTLE at t = 26.1 is 0.0015–0.0027 against the sea's 0.4150. That matches the
  nine-component run to about 5 % (76–86 against 137).
- **Two statements, separated.** (a) The halo is sticky in the base flow's sense [grid]: the r = 10⁻² shell has an FTLE
  between island and sea (0.0814 at t = 26.1, against 0.0015–0.0027 and 0.4150) and the widest spread of escape times
  (log10 spread 0.41 against 0.25–0.28; p90 156). (b) Its ghost escape is not slower than the sea's at E = 1 [grid]:
  median 84.9 against 136.1, in both codes. The old sentence joined (a) to a version of (b) the log does not show.
  (29 Oct, from the direct test below: the r = 10⁻² shell's intermediate FTLE is proximity, a stay of about 10 periods,
  not trapping; the sticky layer proper is much closer to the island.)
- **A reading, the author's to keep or cut [est].** Stickiness near an island is usually attributed to cantori, the
  remnants of broken KAM tori that slow transport across the island boundary (literature from memory, to verify:
  Aubry–Mather; MacKay, Meiss & Percival 1984; Karney 1983; Chirikov & Shepelyansky 1984). Time held near the island
  then lengthens or shortens the ghost's escape depending on what the ghost does near the orbit. At E = 1 the orbit is
  in a diagonal-ghost tongue (Floquet rate 0.233, smoke run, no committed log) against a sea rate of 0.088
  (`ym/notes/NOTE_rung1.md`), so a sticky halo escapes faster than the sea. At E = 0.3, off the tongue, where the core
  is a non-escaping hole, the same kind of halo escapes slower than the sea (control `_od0`: 5–78 % of the halo shells
  against 100 % of the sea). Tested directly 29 Oct: next section.

## Cantori: the direct test (29 Oct 2026)

`islands/code/cantori_residence.py`, logs `islands/results/cantori_{edge,survive,shells,trace}.txt`, figures
`islands/results/cantori_survivors.png`, `cantori_tail.png` (`islands/code/plot_cantori.py`). Ghost-free planar
Yang–Mills at E = 1/2 (scale-free, so every E), section x = 0, vx < 0, where the DR orbit's period map P² has its fixed
point c = (0.781328800089, 0.481317450222) with rotation number 0.455360 (near 5/11). Distance from the island is the
area-normalised radius ρ of the linearised map's invariant ellipse (aspect 27). The criterion for a positive or negative
result was written into the script's docstring before the production run. All CHECKs pass (fixed point, area
preservation, the fixed point under the stepper, energy ≤ 5.3 × 10⁻⁸, the h/2 rerun within 0.18 of its tolerance).

- **The regular region is not round.** The radial scan finds intact curves to ρ = 2.15 × 10⁻⁴ and escapes from
  3.16 × 10⁻⁴ at 5 of 8 angles; the 28,164 starts (14.1 %) that never leave by t = 300,000 fill a region with three
  corners reaching ρ ≈ 6.5 × 10⁻⁴ (`cantori_survivors.png`, last panel).
- **200,000 chaotic-layer starts** (ρ₀ ∈ [2, 8] × 10⁻⁴, uniform in area), first exit past ρ_out = 1.5 × 10⁻³:
  - The bulk leaves exponentially: rate 0.0090 per unit time from t = 300 to 500, i.e. within about 30 periods.
    That exponential predicts 17 exiting starts left at t = 1,000 and 0.002 at t = 2,000. Observed: 203 and 168.
  - The rest leave on a long tail: exits per decade 101 in [10³, 10⁴), 76 in [10⁴, 10⁵), 26 in [10⁵, 3 × 10⁵).
    The local log-log slope of the survival is 0.27–0.37 from t ≈ 2 × 10³ to 2 × 10⁴ and 0.57 around 2–5 × 10⁴
    [grid]. It steepens near tmax, where counting only the starts that exit by tmax cuts the curve off.
  - The tail is 1.2 × 10⁻³ of the exiting starts, and it is the same at every ρ_out from 1.5 × 10⁻³ to 10⁻¹ (76 and 26
    in the last two decades at each): once a long-lived orbit leaves the layer it goes straight to the sea.
- **Where the long-lived orbits live.** The 12 longest (exits at t = 1.6–2.7 × 10⁵, about 5,300–8,600 periods) all have
  ρ percentiles 10/50/90 % = 3.06/3.85–3.90/6.36–6.51 × 10⁻⁴: one thin layer, the same for every one. The longest's
  17,236 section points trace the outer boundary of the regular region, corners included (`cantori_tail.png`, right).
  At t = 3,000 and 30,000 the surviving starts sit on that boundary (`cantori_survivors.png`).
- **The proximity control** (`cantori_shells.txt`, Euclidean discs in (y, vy), 4,000 starts each, to t = 30,000): from
  the r = 10⁻² disc, 3,998 of 4,000 leave ρ = 0.1 by t = 500 and 2 never do (regular). Starts that far out stay about
  10–15 periods and show no tail.

**Result, by the stated criterion: positive [grid].** Exits continue for more than two decades past the bulk, at
a roughly constant log-log slope, far above the bulk's exponential, from orbits confined to a thin layer hugging the
outermost KAM boundary. That is the standard stickiness signature of an island boundary, which the literature attributes
to cantori and the island-around-island hierarchy (from memory, to verify: MacKay, Meiss & Percival 1984; Karney 1983;
Chirikov & Shepelyansky 1984; Meiss 1992 review). Not separated here: cantori from satellite-island trapping (the
regular region's corners suggest a resonance chain, unresolved); the exponent (censored, and of a conditional survival,
so not directly the literature's recurrence exponent); 3 dof, where Arnold diffusion competes; and the ghost.
**Consequence for the halo statements above:** at r = 10⁻² the E = 1 halo is proximity. The cantori-type layer lies at
ρ ≈ 3–7 × 10⁻⁴, inside the shells `hdym_island.py` calls halo, and holds of order 10⁻³ of the chaotic starts near the
edge. The physics reading is the author's. Background (tori, cantori, dimension, the ghost): `islands/notes/NOTE_tori.md`.

## For the paper

Watch-list item 3 (`FOR_THE_REWRITE_2026-10-25.md`): the open item (iii) on the regular islands can now say that the
E = 0.3 island is emptied by the off-diagonal sector throughout its core at ε = 10⁻⁸ (this note) and by the diagonal
ghost's back-reaction above ε ≈ 5 × 10⁻⁸ (hole 10); the positive-measure question for the full theory remains open
only as a statement about other islands. Hole closed.
