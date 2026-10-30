# islands/ — the Dahlqvist–Russberg island, the ghost on it, and a computer-assisted proof

Work of 8–23 Oct 2026; rewritten 23 Oct 2026 from the state at that date. Everything here concerns the stable
periodic orbit of Yang–Mills mechanics found by Dahlqvist & Russberg (1990), with initial data from Contopoulos &
Harsoula (arXiv:2302.12071), and what the Lee–Wick ghost does on and around it — the paper's open item on the regular
islands. The side questions this raised are worked in `rabbit_holes/ym/`.

## Results

**The orbit and the island (proved).** The DR orbit at E = ½ exists and is elliptic, Hénon index
HI ∈ [−0.96096564804637046, −0.96088055206183126], non-resonant to order 4, and its Birkhoff twist is nonzero,
τ ∈ [−3.0474820853652777 × 10⁵, −1.2280792906437143 × 10⁵] — all by a certificate checked in exact rational arithmetic
(`results/cert16_check_output.txt`, `results/cert16_twist_output.txt`, both ending `PROVED: unique periodic orbit in X at
E = 1/2`). So the island exists, modulo Moser's twist theorem, which remains an import. The certified box contains
x₀ ∈ [3.14640122769836, 3.14640122769838], p_x0 ∈ [0.00179319341881, 0.00179319341882] (the Krawczyk image K, rounded
outward; the box X itself has width 2 × 10⁻¹¹ in x₀ and p_x0, exactly one periodic orbit in it), agreeing with the
published DR data (x₀ = 3.14640122769753, p_x0 = 0.0017931934191, `code/dr_orbit_ghost.py`) to about 10⁻¹² absolute (12 significant digits in x₀, 9 in p_x0), not beyond (stated as "to
every published digit" until 28 Oct 2026); the
certified period is the full period (lemma RET: p_y(T) ∈ [0.9999983922, 0.9999983923] > 0, not the reflected half-period
return), and T is the 9th crossing of the global section {y = 0, p_y > 0} (float, `code/dr_orbit_ghost.py`; not checked
by the certificate), so the twist's return map is the local Poincaré map of the orbit, as Moser's theorem needs.

**The diagonal ghost on it.** Floquet-unstable at E = ½ (proved: real multiplier off the circle; ghost trace
∈ [5.3469864, 5.3471379]) and, by float scan, for all E ≳ 0.315 plus narrow tongues below (0.09, 0.16–0.19, 0.24–0.26;
the "diag-ghost rate" column of `results/esc_E.txt` also shows thin tongues at E = 0.0531 (rate 0.0010) and 0.1239
(rate 0.0037)); Floquet-stable between tongues. In the full 12-d Lee–Wick system (`code/hdym_island.py`) with a ghost seed
ε = 10⁻¹⁰ the island is a non-escaping hole at E = 0.3 (no escape within tmax, 5000 at E = 0.3 by the default
min(5000, 40/(0.027 E^2.24)), `code/hdym_island.py`) and a fast smooth
plateau at E = 1 (inside the unstable region), with a halo (r = 10⁻²: FTLE between island and sea, the widest spread of escape times) whose ghost escape is
nevertheless faster than the sea's, median 84.9 against 136.1; that halo is proximity, not trapping, and the island's
genuinely sticky layer (cantori-type, an algebraic exit-time tail) is much thinner, just outside its last KAM curve
(`code/cantori_residence.py`, 29 Oct 2026) (`results/island_E1_eps1e-10.log`, rerun 29 Oct 2026;
stated as "a sticky halo that escapes slower than the sea" until then, from a log that was not committed; the
nine-component run `results/island_full_E1_eps1e-10.txt` agrees; RH9, last section, has a reading); at the paper's
ε = 0.01 the ghost's back-reaction (|w| = √(2εE) = 0.077 at E = 0.3, 0.1 at ½, 0.14 at 1, against an island
half-width ~2 × 10⁻⁴ in x at E = ½, scaled by (2E)^{1/4} elsewhere) erases the island. The crossover is at
|w| ≈ half-width: intact up to ε = 10⁻⁸, leaking at 10⁻⁷ (core escape 21–25 % with the core FTLE still ≈ 0: carried
out by the back-reaction, not by chaos), erased at 10⁻⁶, ε_c ≈ (half-width)²/2E ≈ 5 × 10⁻⁸ at E = 0.3
(`results/island_E0.3_eps1e-{9,8,7,6}.log`, 22 Oct; `rabbit_holes/ym/notes/RH10_eps_crossover.md`). The same
hole / plateau / erasure pattern holds in the exact 2-dof planar subsystem (`--planar` runs); this does not mean the
3-dof island is sealed — in 3-dof it leaks by Arnold diffusion (2.7 % at r = 2 × 10⁻⁴; no committed log), and the
z-direction Floquet trace −1.53 (float, `code/dr_orbit_ghost.py`; its output is not committed) is not certified.

**The off-diagonal ghost on it: not a refuge.** The six off-diagonal components (`ym/rung1_system.pkl`, c_j = 1)
on the DR orbit are Floquet-unstable over about half of E ∈ [0.05, 1] (13 intervals), including E ∈ [0.281, 0.494];
at E = 0.3 block (u23, u32) has the real multiplier 1.5522 (rate 0.0125, against the sea-averaged g₁ = 3.37 × 10⁻³ of
the paper's off-diagonal table) while the diagonal ghost is Floquet-stable there. Nonlinearly (`code/island_full.py`,
23 Oct): with the off-diagonal components seeded at 10⁻⁸ every core shell r = 10⁻⁶ … 10⁻⁴ of the E = 0.3 island
escapes, 200 of 200 per shell, median t_esc 1340–1360 (ln(10⁸)/0.0125 ≈ 1500), the diagonal ghost the escaping
sector in 91–96 % of cells and block 23 in the rest; with no off-diagonal seed (`_od0`) the island is intact, 0 of 200
in the four core shells, and the halo and sea match hole 10's table (`rabbit_holes/ym/notes/RH9_island_full.md`).

**The typed instabilities and the half-period lemma.** Every unstable multiplier is either a real pair leaving the
circle at +1 (parametric) or a complex quadruplet (a Krein collision of two unit-circle pairs of opposite
signature). Krein quadruplets occur in block 12 only, E ∈ [0.0761, 0.0773], [0.1519, 0.1613], [0.3469, 0.4937]
(max |μ| = 2.87, rate 0.032, at E = 0.412); every other unstable cell is parametric at +1; no pair ever leaves at −1
(`results/dr_offdiag_floquet_typed.txt`). The reason is exact: the orbit returns at T/2 with (x, p_x) restored and p_y
reversed, so the block monodromy is M = N² with N = S̃Φ_{T/2} real symplectic, and a real negative multiplier of M off
the circle forces a quadruplet in N — a single pair cannot fold at −1 (`code/dr_offdiag_halfperiod.py`; the derivation
and its literature in `rabbit_holes/ym/notes/RH7_selection_rule.md`, `RH7_literature.md`). Measured signatures
(`results/krein_b005.txt`, 400 energies): of eleven energies where the two block-12 pairs coincide in M, five are true
meetings in N (E = 0.0756, 0.0779, 0.1508, 0.3443, 0.4975) and six are fake degeneracies (N has λ and −λ: E = 0.054,
0.056, 0.094, 0.199, 0.283, 0.798) (stated as "six true, five fake" until 28 Oct 2026: miscount; the log's "meetings"
block); the six edges of the three Krein bands are the log's "transitions" entries marked COLLISION (E = 0.07613,
0.07787, 0.15190, 0.16252, 0.34693, 0.49745), and at every one the pairs carry opposite Krein sign (ψ₁⁺ = −ψ₂⁺ to 0.03–0.09 π) and meet
at N ≈ ±i; every parametric event is one pair's own fold at N = ±1. The frozen-channel background (x = D, y, z
transverse oscillators; `code/channel_offdiag_monodromy.py`) is Floquet-stable in every block within its validity
(ε = 4π²J/D³ ≤ 0.1, 188 of 576 cases); its one apparent weak Krein cell (D = 3.89, J = 0.0315, |μ| = 1.002) is the gauge
Jordan pair split by √ε (`rabbit_holes/ym/notes/RH11_weakcell.md`). The DR orbit itself has ε ≈ 0.2, so its
instabilities belong to the channel bounce, not the adiabatic interior; on the exact gauge-reduced block the unstable
mode of the strong Krein band grows e^{1.4–1.7} at each bounce and shrinks e^{−1} at each turning point
(`rabbit_holes/ym/notes/RH3_segments.md`).

*Scope and checks of the Floquet scans.* `dr_offdiag_floquet.py`: DOP853 rtol 10⁻¹², orbit closure ~3 × 10⁻¹³, det M = 1 to
6 × 10⁻¹⁰ (stated as 10⁻¹⁰ until 28 Oct 2026; worst row E = 0.75 of `results/dr_offdiag_floquet.txt`), multipliers in
(μ, 1/μ) pairs to 1.8 × 10⁻⁹, and with c_j = 0 the block reproduces the scalar Hill equations of the paper's
off-diagonal section. B001 (typed): 400 log-spaced energies in [0.05, 1]; 60 Krein cells, 137 parametric; every unstable
eigenvector satisfies the linearised Gauss law to 5 × 10⁻¹¹; the pure-YM off-diagonal block is stable at every E (so the
E = 0.3 block-23 instability is the ghost's); closure 4 × 10⁻¹³, det M − 1 ≤ 10⁻⁹, pairing 10⁻⁹, c_j = 0 Hill reduction
10⁻¹¹; in the Broucke (A, B) diagram block 12 rides the seam between the Krein parabola and the period-doubling line from
E ≈ 0.07 to 0.49 and always exits through the parabola, blocks 13 and 23 through the +1 line only. B002 (frozen channel):
D ∈ [0.5, 6] (24 log points), J ∈ [0.005, 0.5] (6), ψ ∈ {0, π/2}, c_j ∈ {1, 0}; the depth changes by ε = 4π²J/D³ over a
period, cells with ε > 0.1 marked outside the model, the c_j = 0 Hill check holding to O(ε) (9 × 10⁻⁶ at ε = 3 × 10⁻³); the
one in-model Krein cell had rate 0.001; every strong instability lies outside the model and most of it violates Gauss.
B005: MᵀΩM = Ω to 10⁻⁸–10⁻¹¹ in all three blocks at E = 0.3 and to 9.8 × 10⁻⁸ worst over the block-12 scan
(`results/krein_b005.txt`, "checks over the scan"), N² = M to 9 × 10⁻¹¹, NᵀΩN = Ω to 6 × 10⁻⁹; the "pass-through at
E ≈ 0.76" discussed on 19 Oct is the fake degeneracy at 0.798, where nothing meets.

**Nonlinear escape from the orbit** (`code/offdiag_escape.py`, kernel `ym/csrc/hdym_full_escape.c`: the full nonlinear
homogeneous SU(2) Lee–Wick system, all nine components, temporal gauge; model from `ym/hdym_full_derive.py`). Section E
(100 energies × 64 seeds, both seeds 10⁻⁴, tmax 3000, rtol 10⁻¹¹; `results/esc_E.*`, 20 Oct): the escape times lie
between ln(1/δ)/rate and ln(R/δ)/rate of the leading linear instability at every E. Where the diagonal ghost is
Floquet-unstable (E ≳ 0.31 and its tongues) it leads, and above E ≈ 0.33 its rate (0.04–0.23) exceeds every off-diagonal
rate, so block 12's Krein band at 0.35–0.49 (rate ≤ 0.032) never leads an escape there. In the gaps where the diagonal
ghost is stable the off-diagonal instabilities take over at their own linear rates — E = 0.153–0.158 (block 12 Krein,
median 2276 vs 1460–2900), 0.227–0.246 (block 12 parametric, 914 vs 575–1150), 0.28–0.31 (block 23 parametric) — and in
each the off-diagonal mode grows first and the diagonal ghost, pumped nonlinearly, overtakes and does the escaping.
Where nothing is linearly unstable nothing escapes by 3000 at this seed:
∀ E ∈ G_esc ∩ ([0.113, 0.149] ∪ [0.195, 0.221] ∪ {0.27221}) ∖ {0.12394}: rate_diag(E) = 0 ∧ escaped(E) = 0/64 [grid]
(G_esc: the 100 energies of `results/esc_E.txt`; columns "diag-ghost rate", to its printed 4 decimals, and "escaped"). The excluded row E = 0.12394 is a
thin diagonal tongue (rate_diag 0.0037, 32/64 escape with median t_esc 2622 near tmax 3000, so the half is censoring,
not a probability); the range 0.113–0.149 was stated without this exception until 28 Oct 2026. Model and kernel checks: A'''' coefficient −I, diagonal reduction, linearisation = `rung1_system.pkl`'s K to 10⁻¹⁵, Gauss
law and Ostrogradsky energy conserved to 10⁻¹⁸, generated C = sympy to 10⁻¹⁵; the integrator agrees with scipy DOP853 to
6 × 10⁻⁹ over two periods (logged: 5.9 × 10⁻¹⁰–2.3 × 10⁻⁹ in every `results/esc_*.txt` at δ = 10⁻⁴; at δ = 10⁻² this
CHECK and the escape-time CHECK fail, `esc_E_delta1e-2_b12.txt`). Smoke run (19 Oct, 6 energies × 8 seeds at 10⁻⁶,
tmax 1500; no committed log): at E = 1 the diagonal ghost
(Floquet rate 0.233) escapes at t ≈ 96 against ln(R/δ)/rate = 99; at E = 0.30 block 23's parametric pair drives the escape at
t ≈ 950 (linear 1100–1840); E = 0.05, 0.09 no escape by 1500. The
orient maps (seed-direction sphere at fixed E, 128²): at E = 0.412 structureless (median 183, 96 % of cells escape
through the diagonal ghost); at E = 0.30 (block 23's parametric tongue, 100 % escape, median 631, p10–p90 525–886) and
E = 0.156 (block 12's Krein band, 88.7 %, median 1952, 1527–2606) the maps carry the mechanism as geometry
(`code/offdiag_escape_nodal.py`; its `*_nodal.png` figures and its printed correlations are not committed, and the
lag-1 correlations and the quadruplet multiplier below have no committed log): a seed grows along an unstable Floquet mode with amplitude |w·s|, so
delayed escape sits on the nodal set of the unstable subspace — a great circle for a real pair (multiplier 1.552; corr(log t_esc,
log |proj|) = −0.44; lag-1 correlation of the map 0.28/0.20), two antipodal points for a quadruplet (−1.272 ± 0.181i;
corr −0.11, lag-1 0.09: the loiter dominates the texture, as in the paper's Figure 3). At escape the diagonal ghost is the largest sector in 90.5 % (E = 0.30: 14821 of 16380 cells of `esc_orient_030.npz`, `sector`; stated as 93 % until 28 Oct 2026, source not recorded) and 76 % (E = 0.156) of cells:
the off-diagonal mode triggers, the diagonal ghost executes. What follows the linear phase — the loiter, the low-energy
Poisson regime, the seed-size loophole at δ = 10⁻², the absence of a phase effect, the nonlinear indistinguishability of
the two mechanisms — is measured in `rabbit_holes/ym/` (holes 1, 4, 8, 12) from the runs listed under Layout.

## Layout

```
code/
  dr_orbit_ghost.py            DR orbit, Hénon index, coarse diagonal-ghost Floquet scan in E              (float)
  scan.py                      fine Floquet scan: python scan.py E_min E_max N                              (float)
  interior.py                  base and ghost Lyapunov exponents inside the island                          (float)
  hdym_island.py               escape times / FTLE around the DR orbit in the 12-d system (diagonal ghost; imports
                               hdym_escape_map.py, hdym_precision.py from ym/; pairs through hdym_precision.CPairs
                               with the SHARED sub-step schedule, a truncation-level difference from production)
  island_full.py               the same ball sampling through the nine-component kernel (off-diagonal seed
                               --delta-od, Gauss-projected; conversion checked to 1e-15)              (added 22 Oct)
  dr_offdiag_floquet.py        the six off-diagonal components (ym/rung1_system.pkl, c_j = 1) on the DR orbit:
                               block monodromies vs E                                                (float; 12 Oct)
  dr_offdiag_floquet_typed.py  B001: the same with the TYPE of each unstable multiplier (REAL = parametric,
                               CPLX = Krein), Gauss violation of unstable eigenvectors; imports the above      (19 Oct)
  channel_offdiag_monodromy.py B002: the blocks on a frozen-depth channel background, monodromy over one
                               transverse period vs (D, J), typed                                            (19 Oct)
  dr_offdiag_halfperiod.py     is M the square of a symplectic half-period map N = S~ Phi_{T/2}? (the lemma) (19 Oct)
  dr_offdiag_krein.py          B005 on N, block by block, 400 energies; plot_offdiag_krein.py its figure      (19 Oct)
  krein_b005.py                B005 for block 12: Omega from the Ostrogradsky structure of L2 (CHECK M^T Omega M =
                               Omega), N (CHECK N^2 = M), every meeting of the two pairs classified           (21 Oct)
  offdiag_escape.py            NONLINEAR escape times of the nine-component system seeded with off-diagonal (and
                               diagonal-ghost) perturbations: sections E (energy x seeds), orient (seed-direction
                               map at fixed E), phase (ensemble over the seed's phase on a Krein quadruplet);
                               records crossing times (tcross), growth curves (chk) and state snapshots (snap);
                               needs ym/hdym_full_derive.py run once                                 (19–22 Oct)
  offdiag_escape_nodal.py      overlays the linear prediction (nodal set of the unstable Floquet subspace) on an
                               orient map                                                                    (21 Oct)
  cap_dr_orbit.py              UNTRUSTED certificate generator (mpmath intervals, Lohner/Taylor)
  cert_lemmas.py               TRUSTED lemma layer: pure Python Fractions, one function per lemma (401 lines; 387 until 28 Oct 2026)
  cert_check.py                TRUSTED driver: parses the certificate, calls the lemmas, prints verdicts (214 lines; 208 until 28 Oct 2026)
  cert_twist.py                TRUSTED layer for the twist: VAR2/VAR3, FRAME2, TWIST (351 lines)
  test_checker.py              regression tests of the checker (seconds)
  twist_symbolic.py            exact sympy formula for the Birkhoff twist + float evaluation
  twist_symbolic2.py           writes twist_formula.json, the return-map formula cert_twist evaluates
  twist_formula.json           generated by twist_symbolic2.py; md5 64cc7b9c35a3096454605a4c3013c2d5 (trusted input)
  twist_estimate.py            independent twist estimate from rotation numbers
results/
  cert16.json                  certificate (cap_dr_orbit.py 16 1000 1e-11), md5 690878da67b071b6ce277aea0fb02e98
  cert16_check_output.txt      --ghost --E 1/2 --verbose, post-review checker with RET (12 Oct): PROVED
  cert16_twist_output.txt      --ghost --twist (12 Oct): PROVED, tau nonzero (predates the RET print line; RET is a
                               property of the certificate the twist layer reads unchanged)
  cert.json                    p = 24 certificate (generator default), md5 b2dc3efca31bd6b848048b34cc8527bc; no
                               checker output committed
  twist_sym.log                stdout of twist_symbolic.py (float twist tau = -2.140067e5)
  dr_offdiag_floquet.txt       12 Oct, cloud, python 3.13
  dr_offdiag_floquet_typed.txt 19 Oct, python 3.12.5 numpy 2.2.6 sympy 1.13.1, 20 cores; a cloud run agrees row for row
  channel_offdiag_monodromy.txt 19 Oct; cloud run agrees in every in-model cell
  dr_offdiag_halfperiod.txt, dr_offdiag_krein.{txt,npy}, offdiag_krein_b12.png   B005 on N, 19 Oct
  krein_b005.txt               21 Oct, cloud, 514 s
  hdym_full_derive.txt         the model's checks
  (every esc_* log prints "diagonal ghost not seeded" on line 3, a fixed string; the run parameters line
                               gives the actual delta_diag, 1e-4 except esc_E_delta1e-2_b12; the print is fixed in the code from 28 Oct 2026,
                               the committed logs keep the old text)
  esc_E.{txt,npz}              section E, 100 x 64, 20 Oct, 24 threads, 16893 s (see NOTE_esc_E_kernel_2026-10-20.md)
  esc_orient_0412.*            128^2 at E = 0.412, 32975 s (same kernel as esc_E)
  esc_orient_0156.*, esc_orient_030.*   128^2 at E = 0.156 (19987 s) and 0.30 (38688 s), 21 Oct, fixed kernel, one
                               diagonal seed per map; no chk arrays. (An E = 0.642 map started on the pre-fix kernel
                               with the diagonal ghost unseeded was stopped at 83 of 128 rows on 21 Oct and deleted on
                               23 Oct; every row had escaped 128/128 — above E ≈ 0.33 the diagonal ghost leads regardless.)
  esc_phase_0156.*             phase ensemble, 12 x 256, 21–22 Oct night, 2802 s         (rabbit hole 1)
  esc_orient_030_n64.*         64^2 at E = 0.30 with chk, 13948 s                       (rabbit hole 12, superseded by snap)
  esc_E_cross.*                section E with crossing times, 20 x 64, 3162 s           (rabbit hole 4)
  esc_E_lowE_long.*            20 energies in [0.05, 0.13] x 64, tmax 12000, 22 Oct, 6559 s   (rabbit hole 4)
  esc_E_delta1e-2_b12.*        48 energies in [0.045, 0.2] x 64, block-12 seed 1e-2 only, 7963 s   (rabbit hole 8;
                               two kernel CHECKs FAIL at this seed, log ends SOME CHECK FAILED)
  esc_E_snap_krein.*, esc_E_snap_param.*   4 x 128 in each block-12 band with state snapshots, 297 s / 822 s   (hole 12)
  island_E0.3_eps1e-{6,7,8,9}.{log,json}   hdym_island.py, 22 Oct                        (rabbit hole 10)
  island_E1_eps1e-10.log, island_E1_eps1e-10_c0_N200.json   hdym_island.py --E 1.0 --eps 1e-10 --N 200 --backend c,
                               29 Oct (the json was written to results_island/ at the repository root, the script's
                               old default, and moved here; the default is now islands/results)
  cantori_{edge,survive,shells,trace}.{txt,npz}, cantori_survive/, cantori_shells/ (per-chunk npz)
                               cantori_residence.py, 29 Oct, 22 workers; survive 7623 s (rabbit hole 9, last section)
  cantori_survivors.png, cantori_tail.png   plot_cantori.py from the npz above
  island_full_E0.3_eps1e-8.{txt,npz}, _od0.{txt,npz}   island_full.py, 23 Oct             (rabbit hole 9)
  island_full_E0.3_eps1e-8_s8.*, island_full_E1_eps1e-10.*   island_full.py, 28–29 Oct, 6 threads, 4281 s / 3625 s
                               (rabbit hole 9: second seed; E = 1)
  esc_dscan_{1e-4,3.1623e-4,1e-3,3.1623e-3}_b12.*   section E at E = 0.13239 and 0.19375 x 64, block-12 seed only,
                               28–29 Oct, 6 threads (rabbit hole 8, δ-crossover; the escape-time CHECK fails at
                               3.1623e-4 and 3.1623e-3, logs end SOME CHECK FAILED)
  NOTE_esc_E_kernel_2026-10-20.md   why esc_E and esc_orient_0412 stand without a rerun after the kernel's FP fix
```

Not committed: the `hdym_island.py` logs at E = 1 (until 29 Oct 2026, now `island_E1_eps1e-10.log`) and at E = 0.3 with ε = 10⁻¹⁰ and 0.01 (two planar), the outputs of
`scan.py` / `interior.py`, the figures of the 19 Oct Broucke-diagram brief, and the logs of the 10 Oct review's
corrected step (which passed all 1000 + 1000 steps of chains P and B), the 19 Oct smoke run of `offdiag_escape.py`,
the output of `dr_orbit_ghost.py` (and its `dr_ghost_floquet.npy`, written to the working directory) and of
`offdiag_escape_nodal.py` (its printed correlations), and the 3-dof Arnold-diffusion run behind the 2.7 %.

## Reproduce the proof

```
python3 code/test_checker.py                                                                       # seconds
python3 code/cap_dr_orbit.py 16 1000 1e-11 results/cert16.json                                     # ~5 min; untrusted
python3 code/cert_check.py results/cert16.json --ghost --E 1/2 --verbose > results/cert16_check_output.txt     # ~18–35 min
python3 code/cert_check.py results/cert16.json --ghost --twist code/twist_formula.json --E 1/2 --verbose > results/cert16_twist_output.txt   # ~57 min (362 s + 3046 s in the log; stated ~48 until 28 Oct 2026)
```

The checker verifies, step by step, with outward-rounded dyadic rationals: a Picard rough enclosure per step (PIC) and
of the variational matrix (PICV), the Lagrange remainder of the Taylor step (LAG), the mean-value Jacobian (MVT), a
Neumann-series inverse (NEU), set inclusion in the Lohner frame (step), the accumulated Jacobian in a moving frame
(FRAME), the full-period return p_y(T) > 0 (RET), Krawczyk's inclusion (KRA), and the trace/symplectic stability lemmas
(T, SYMP); for the twist, VAR2/VAR3, FRAME2, TWIST, plus the area-preservation identities (det L ∋ 1, the zetabar
coefficient ∋ 0, Re(conj λ c₁) ∋ 0) as a guard against a wrong formula. Lemma names are the function names. The
generator may be wrong in any way; a bad certificate fails to check. Trusted: the checker (`cert_lemmas.py`,
`cert_check.py`, `cert_twist.py`), Python integer arithmetic, `fractions.Fraction`, `math.isqrt` and `math.factorial`, the `json` parser and
`float.fromhex` (parsing hex strings to exact dyadics), and for the twist the formula
`code/twist_formula.json`, hence its generator `twist_symbolic2.py` and sympy (regenerated byte-identical on 11 Oct).
The sign of τ and the resonance checks are decided on exact Fractions; printed endpoints are outward-rounded.
Regenerating `cert16.json` on 10 Oct gave the same file (same md5). Design target: Lean — each lemma is an elementary
statement about polynomial ODE enclosures or rational linear algebra; Moser's twist theorem would remain an import.

## Review history of the checker

A review started on 10 Oct, two independent cold reviews on 11 Oct (`notes/REVIEW_cert_lemmas_2026-10-11.md`,
`notes/REVIEW_cert_twist_2026-10-11.md`) and one out of project context on 12 Oct found the lemmas sound and fixed,
the same day each: MVT's
variational remainder used X(ξ) = I (now a Picard enclosure, PICV; same fix in VAR2/VAR3); the step box did not have to
contain the centre (now the hull with c); `assert` soundness checks (now CheckFailed); ghost frames chosen per step,
which for an adversarial certificate could insert a non-identity factor into the ghost monodromy (now G at every step
or none; not triggered by cert16); the initial-data inclusion checked on the outward-rounded box (now exact); verdicts
printed before Krawczyk succeeded; endpoints rounded inward in the printout (two quoted τ endpoints had been rounded
inward from the printed doubles); E neither printed nor compared with ½ (now printed exactly, optional `--E`); lo ≤ hi unchecked at parse; the area-preservation
identities printed but not enforced; and — the 12 Oct finding — nothing checked p_y(T) > 0, although a float
integration shows the orbit reaches y = 0 at t = 15.52667 = T/2 with (x, p_x) restored and p_y = −0.999998, so a
certificate converging there would have "proved" verdicts about the half-period map (lemma RET; regression test 5).
Every certified interval was unchanged by the fixes (the earlier outputs had the same numbers, rounded to nearest;
float twist τ ≈ −2.140067 × 10⁵; cos θ ∈ [−0.960963, −0.960884]).
Not resolved by the twist review: an independent float fit of the return-map Taylor coefficients. `test_checker.py`
tests 2 and 3 fail on the checker of 10 Oct and pass since 11 Oct.

## Open

- [ ] Second certificate at E = 0.3 (ghost Floquet-STABLE on the orbit).
- [ ] 3-dof: the z-direction Floquet trace (−1.53, float) is not certified; the 3-dof island leaks by Arnold diffusion.
- [ ] Independent float fit of the twist's return-map coefficients (twist review, not resolved).
- [ ] Reduced-block items: `rabbit_holes/ym/README.md`, Open. (The δ-crossover scan was done 28–29 Oct 2026: RH8.)

## Conventions

`hdym_island.py` uses the repository's conventions: state (p, w, P, W), q = p − w, E = |P|²/2 + V(q), ghost
|w|² + |W|² = 2εE, p = q + w, P = P_q − W, Yoshida-4 with dt = 0.005 and channel sub-stepping, escape at any
|component| > 10⁴. The DR point in these coordinates: q* = λ(3.14640122769837, 0, 0), P* = λ²(0.00179319341881,
√(1 − p_x*²), 0), λ = (E/½)^{1/4}; copies 0–23 by cyclic permutation, 90° in-plane rotation, time reversal. The
escape kernel `hdym_full_escape.c` uses DOPRI5 with the standard err^(−1/5) controller (PI before 20 Oct), pins the
floating-point state (MXCSR 0x1F80 and the x87 control word) in every thread and computes the step factor with pure
arithmetic, so that escape times are bit-reproducible across threads and platforms; a 48-copy determinism CHECK and a
scipy DOP853 cross-check (state at 2T to 10⁻⁶, escape time to 10⁻³) run at every launch. Runs before that fix
(`esc_E`, `esc_orient_0412`) had about 1/24 of their cells integrated in a different FP mode on one thread (found
when 123 bit-identical cells of a θ = 0 row escaped at 2965.95 and 5 at 1818.84; a first FTZ/DAZ reading was wrong,
the 48-copy determinism check still failing with a t_esc spread of 0.087 at E = 0.549 until the x87 word was pinned); their
statistics stand (`results/NOTE_esc_E_kernel_2026-10-20.md`). Escape-time convergence in tolerance holds only for
seeds ≳ 10⁻⁵, hence the default δ = δ_diag = 10⁻⁴; an unseeded ghost-free orbit at E = 1 escapes at t ≈ 120 from the
integrator's 10⁻¹⁰ truncation seed, which is why both sectors are seeded at equal amplitude.

## History

8–9 Oct: orbit, Hénon index, ghost Floquet scan, certificate and checker, island runs. 9–10 Oct: twist enclosure,
first review and fixes, reruns. 11 Oct: two independent reviews, fixes, regression tests. 12 Oct: reruns with the
post-review checker; RET; the off-diagonal ghost found Floquet-unstable on the orbit (the island not a refuge in the
full theory). 19 Oct: B001 (typed), B002 (frozen channel), B005 on the half-period map. 19–21 Oct: the nonlinear
escape kernel, section E, the orient maps, the FP fix (a first FTZ reading was wrong; the cause was the x87 control
word used by MinGW's `pow()`), the nodal-set reading, B005 for block 12. 21–23 Oct: the rabbit-hole runs (Layout);
`island_full.py`; the ε-crossover and the nine-component island (holes 10 and 9).
