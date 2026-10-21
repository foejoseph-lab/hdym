# islands/ — the Dahlqvist–Russberg island, the ghost on it, and a computer-assisted proof

Session of 8–9 October 2026 (developed interactively with Claude; see the paper's
acknowledgement). Everything here concerns the stable periodic orbit of Yang–Mills mechanics
found by Dahlqvist & Russberg (1990), with initial data from Contopoulos & Harsoula,
arXiv:2302.12071, and what the Lee–Wick ghost does on and around it.

## Results in one paragraph

The DR orbit at E = 1/2 exists and is elliptic (Hénon index in [−0.96097, −0.96088]) and
non-resonant to order 4 — **proved**, by a certificate checked in exact rational arithmetic.
The linearised ghost on it is **Floquet-unstable** at E = 1/2 (proved, real multiplier off the
unit circle) and, by float scan, for all E ≳ 0.315 plus narrow tongues below; between tongues
it is Floquet-stable. Its Birkhoff twist is τ ≈ −2.14×10⁵ (exact formula, float evaluation;
rigorous enclosure pending), so Moser's theorem applies once τ is certified: the island exists.
*Status (2026-10-10): the enclosure was computed on 9–10 Oct, τ ∈ [−3.05×10⁵, −1.23×10⁵]; its
check output is not yet committed (see Status / open).*
*Status (2026-10-10, later): rerun with the fixed checker (see Status / open): τ ∈
[−3.0475×10⁵, −1.2280×10⁵], nonzero, output in `results/cert16_twist_output.txt`. The island
is now a theorem modulo Moser's twist theorem and the outstanding independent review of the
checker.*
In the full 12-d Lee–Wick system with a tiny ghost seed (ε = 10⁻¹⁰) the island is an
infinite-escape-time hole at E = 0.3 and a fast smooth plateau at E = 1, with a sticky halo
that escapes slower than the sea; at the paper's ε = 0.01 the ghost's back-reaction (|w| ≈ 0.1
≫ island width 2×10⁻⁴) erases the island entirely. The same holds in the exact 2-dof planar
subsystem, so nothing about the three-field lift changes the 2-dof analysis.
*Status (2026-10-12): five qualifications of this paragraph, from the 10 Oct audit and the 12 Oct
off-diagonal computation. (a) Widths: `hdym_island.py` prints the island HALF-width, ~2×10⁻⁴ in x
at E = ½ (scaled by λ = (2E)^{1/4} elsewhere); "island width 2×10⁻⁴" above should read half-width,
and the docstring's "(~1e-4)" is a rougher figure. The ghost amplitude is |w| = √(2εE): 0.077 at
E = 0.3, 0.1 at E = ½, 0.14 at E = 1 (the docstring's 0.07 and the 0.1 above are the same estimate
at different E); either way ≫ the island. (b) "Floquet tongue (e.g. 1.0)" in the docstring and
"Floquet-unstable for all E ≳ 0.315" here describe the same thing: E = 1 lies in the unstable
region of the diagonal ghost. (c) "Infinite-escape-time" means no escape within the run's tmax
(3000 in the usage line, default 5000 at E = 0.3; the logs are not committed). (d) "Nothing about
the three-field lift changes the 2-dof analysis" refers to the hole / plateau / erasure pattern,
which the planar runs reproduce; it does not mean the 3-dof island is sealed: in 3-dof it leaks by
Arnold diffusion (Status / open, 3-dof item), and the `--planar` sea control stays 3-dof.
(e) All of the above is the colour-diagonal ghost. The off-diagonal ghost is Floquet-unstable on
the DR orbit at E = 0.3 (Status / open, off-diagonal item), so in the full homogeneous theory the
E = 0.3 island is not a refuge, at least at its centre.*

## Layout

```
code/
  dr_orbit_ghost.py    DR orbit, Hénon index, coarse ghost Floquet scan in E       (float)
  dr_offdiag_floquet.py  the six off-diagonal ghost components (ym/rung1_system.pkl, c_j = 1)
                       on the DR orbit: block monodromies vs E                (float; added 2026-10-12)
  dr_offdiag_floquet_typed.py  B001: the same monodromies with the TYPE of each unstable multiplier
                       (REAL = parametric, CPLX = Krein quadruplet), all 8 printed where the off-circle
                       count changes, Gauss violation of unstable eigenvectors; imports dr_offdiag_floquet
                                                                                (float; added 2026-10-19)
  channel_offdiag_monodromy.py  B002: the blocks on a frozen-depth channel background (x = D,
                       y,z = sqrt(2J/D) cos), monodromy over one transverse period vs (D, J), typed
                                                                                (float; added 2026-10-19)
  offdiag_escape.py    NONLINEAR escape times of the full nine-component system seeded with off-diagonal
                       (and diagonal-ghost) perturbations of the DR orbit; sections E (energy x seeds) and
                       orient (2-d seed-direction map at fixed E); needs ym/hdym_full_derive.py run once;
                       kernel ym/csrc/hdym_full_escape.c (DOPRI5, OpenMP)             (float; added 2026-10-19)
  scan.py              fine Floquet scan: python scan.py E_min E_max N             (float)
  interior.py          base and ghost Lyapunov exponents inside the island         (float)
  cap_dr_orbit.py      UNTRUSTED certificate generator (mpmath intervals, Lohner/Taylor)
  cert_lemmas.py       TRUSTED lemma layer: pure Python Fractions, one function per lemma;
                       holds the rounding grid                       (split out 2026-10-10)
  cert_check.py        TRUSTED driver: parses the certificate, calls the lemmas, prints verdicts
  cert_twist.py        TRUSTED layer for the twist (VAR2/VAR3, FRAME2, TWIST); used by
                       cert_check.py --twist twist_formula.json            (added 2026-10-10)
  test_checker.py      fast regression tests of the checker (seconds): frame consistency,
                       exact initial inclusion, outward printing          (added 2026-10-11)
  twist_symbolic.py    exact sympy formula for the Birkhoff twist + float evaluation
  twist_symbolic2.py   writes twist_formula.json, the return-map formula cert_twist evaluates
  twist_formula.json   generated by twist_symbolic2.py
  twist_estimate.py    independent twist estimate from rotation numbers
  hdym_island.py       escape times / FTLE around the DR orbit in the 12-d system
                       (imports hdym_escape_map.py, hdym_precision.py from ym/)
results/
  cert16.json          certificate, md5 690878da67b071b6ce277aea0fb02e98
  cert16_check_output.txt   checker output (PROVED).  Status (2026-10-10): produced by the
                       checker before the MVT / centre-in-box fixes (see Status / open);
                       the recheck with the fixed checker is not yet committed.
                       Status (2026-10-10, later): replaced by the output of the fixed
                       checker (run 10 Oct): PROVED, every printed number equal to the
                       earlier output's
  cert16_twist_output.txt   output of `--ghost --twist` with the fixed checker (10 Oct):
                       PROVED, τ nonzero                                  (added 2026-10-10)
  cert.json            p = 24 certificate (generator default), md5 b2dc3efca31bd6b848048b34cc8527bc;
                       no checker output committed
  twist_sym.log
  dr_offdiag_floquet.txt        output of dr_offdiag_floquet.py (12 Oct, cloud session)
  dr_offdiag_floquet_typed.txt  output of dr_offdiag_floquet_typed.py --ncores 20 (19 Oct, Joseph's machine,
                       python 3.12.5 numpy 2.2.6 sympy 1.13.1); a cloud run the same day agrees row for row
  channel_offdiag_monodromy.txt output of channel_offdiag_monodromy.py --ncores 20 (19 Oct, same machine;
                       cloud run agrees in every in-model cell)
```

Not yet committed: the `results/island_E*_eps*.log`
outputs of `hdym_island.py` (E = 1, 0.3; ε = 1e-10, 0.01; two planar), the output of
`scan.py` / `interior.py`, and a `cert_check.py --twist` run.
*Status (2026-10-10): the `--twist` run is now committed (`results/cert16_twist_output.txt`).*

## Reproduce the proof

```
python3 code/cap_dr_orbit.py 16 1000 1e-11 results/cert16.json   # ~5 min; untrusted
python3 code/cert_check.py results/cert16.json --ghost > results/cert16_check_output.txt   # ~35 min; this is the proof
```
*Status (2026-10-11): after the review fixes the commands are*
```
python3 code/test_checker.py                                       # seconds
python3 code/cert_check.py results/cert16.json --ghost --E 1/2 --verbose > results/cert16_check_output.txt
python3 code/cert_check.py results/cert16.json --ghost --twist code/twist_formula.json --E 1/2 --verbose > results/cert16_twist_output.txt
```

The checker verifies, step by step, with outward-rounded dyadic rationals: a Picard rough
enclosure per step (lemma PIC), the Lagrange remainder of the Taylor step (LAG), the
mean-value Jacobian (MVT, with the variational matrix enclosed over the step by PICV), a
Neumann-series inverse (NEU), set inclusion in the Lohner frame
(step), the accumulated Jacobian in a moving frame (FRAME), Krawczyk's inclusion (KRA), and
the trace/symplectic stability lemmas (T, SYMP). Lemma names are the function names. The
generator may be wrong in any way; a bad certificate fails to check. What is trusted: the
checker (`cert_lemmas.py`, 370 lines, and the driver `cert_check.py`, 167 lines, plus
`cert_twist.py`, 338 lines, for the twist; see Status / open for the review), Python integer
arithmetic, and `math.isqrt`.
*Status (2026-10-11): after the review fixes: `cert_lemmas.py` 387 lines, `cert_check.py` 208,
`cert_twist.py` 351. For the twist the trusted base also includes the formula the checker
evaluates, `code/twist_formula.json` (md5 64cc7b9c35a3096454605a4c3013c2d5), hence its generator
`code/twist_symbolic2.py` and sympy; the review of 11 Oct regenerated it byte-identical and read
the derivation. The checker now also requires the area-preservation identities (det L ∋ 1, the
zetabar coefficient ∋ 0, Re(conj λ c₁) ∋ 0), a cheap guard against a wrong formula.*
*Which orbit: the certified box contains x₀ ∈ [3.14640122769836, 3.14640122769838],
p_x0 ∈ [0.00179319341881, 0.00179319341882] (Krawczyk enclosure in `results/cert16_check_output.txt`),
which agrees with the Dahlqvist–Russberg data quoted under Conventions at E = ½
(q* = 3.14640122769837, p_x* = 0.00179319341881) to every digit given there (the quoted p_x* is
rounded: the certified value is 0.0017931934188137–0.0017931934188142). Krawczyk proves there is
exactly one periodic orbit in the box (width 2e-11); it is the DR orbit in the sense that its
initial data agree with the published ones to every published digit.*

*Status (2026-10-12): the lemma list is now PIC, PICV, LAG, MVT, NEU, step, FRAME, RET, KRA, T,
SYMP (plus VAR2/VAR3, FRAME2, TWIST for the twist).*

Design target: Lean. The lemma shapes were chosen so each is an elementary statement about
polynomial ODE enclosures or rational linear algebra; Moser's twist theorem would remain an
import.

## Status / open

- [x] Twist enclosure: add 2nd/3rd-order variational series to generator and checker (VAR2, VAR3),
      evaluate the `twist_symbolic.py` formula in intervals. Then the island is a theorem modulo Moser.
      *Status (2026-10-10): done. Checker side: `cert_twist.py` (VAR2/VAR3, FRAME2, TWIST) with the
      formula `twist_formula.json` from `twist_symbolic2.py`. Run on Joseph's machine 9–10 Oct
      (`cert_check.py results/cert16.json --ghost --twist code/twist_formula.json`, ≈48 min): reported
      τ ∈ [−3.05×10⁵, −1.23×10⁵], nonzero. The output of that run is not in the repository yet; to
      close this, rerun with the output redirected to `results/cert16_twist_output.txt` and commit it.*
      *Status (2026-10-10, later): the τ interval above came from the checker before the fixes
      under "Independent review"; it needs the rerun with the fixed checker.*
      *Status (2026-10-10, latest): rerun with the fixed checker on Joseph's machine,
      output `results/cert16_twist_output.txt`: τ ∈ [−304748.21, −122807.93], NONZERO; cos θ ∈
      [−0.960962, −0.960884], non-resonant to order 4. Consistent with the 9–10 Oct interval.*
      *Status (2026-10-12): two of the quoted endpoints above were rounded inward from the
      printed doubles (found when the checker began printing outward-rounded endpoints, review
      of 11 Oct, C-F3). Outward: τ ∈ [−304748.21, −122807.92], cos θ ∈ [−0.960963, −0.960884].
      The sign of τ and the resonance checks were always decided on exact Fractions.*
- [ ] Second certificate at E = 0.3 (ghost Floquet-STABLE on the orbit).
- [ ] Independent review of `cert_check.py`.
      *Status (2026-10-10): a review was started on 10 Oct and stopped before its written report;
      from its working files, two gaps in the step lemmas and one robustness issue, all fixed the
      same day: (1) MVT's Lagrange remainder of the variational matrix X(h) used X(ξ) = I; now
      multiplied by a Picard enclosure of X on [0, h] (lemma PICV); the same gap in the twist
      layer's VAR2/VAR3 step remainder is fixed the same way (`_dir_rough`); (2) the box B of a
      step did not have to contain the centre c (r is untrusted and need not contain 0), which
      LAG and MVT require; now B is the hull with c; (3) soundness checks written as `assert`
      (dropped by `python -O`) now raise CheckFailed. The reviewer's corrected step passed all
      1000 + 1000 steps of chains P and B of `cert16.json` (base flow only; its logs, not
      committed, are in Joseph's `review_scratch_2026-10-10.zip`). The full recheck (--ghost,
      and --twist) with the fixed checker is pending. The review's per-lemma verdicts were lost;
      the review of the remaining lemmas and of `cert_twist.py` is still open.
      Regenerating `cert16.json` with `cap_dr_orbit.py 16 1000 1e-11` on 10 Oct gave the same
      file (same md5).*
      *Status (2026-10-10, later): full recheck with the fixed checker done on Joseph's machine:
      `--ghost` PROVED with every printed number equal to the earlier output
      (`results/cert16_check_output.txt`); `--ghost --twist` PROVED, τ nonzero
      (`results/cert16_twist_output.txt`). Still open: the review of the remaining lemmas and of
      `cert_twist.py`.*
      *Status (2026-10-11): reviewed cold by two independent sessions, findings in
      `notes/REVIEW_cert_lemmas_2026-10-11.md` and `notes/REVIEW_cert_twist_2026-10-11.md`.
      Lemma by lemma sound, except: (C-F1, breaks `--ghost` for an adversarial certificate) the
      ghost frames were chosen per step, so omitting G on one step inserted a non-identity factor
      into the ghost monodromy; now G at every step or at none (run_chain). (C-F2, gap) the
      initial-data inclusion was checked on the outward-rounded box; now exact. Neither is
      triggered by cert16 (G on all chain-B steps; steps[0] on the grid; both checked by the
      reviewer). Also fixed: verdicts are printed only after Krawczyk succeeds (C-F4 = B-F1);
      endpoints printed outward-rounded (C-F3); E printed exactly, optional --E (C-F6);
      lo ≤ hi checked at parse (C-F5); "elliptic" wording (C-F7); the area-preservation guards
      and the trusted-base note (B-F2, B-F3); docstring and silent-zero items (B-F4, B-F5).
      Regression tests in `code/test_checker.py`; tests 2 and 3 fail on the checker of 10 Oct
      and pass now. The recheck of cert16 with this checker (commands under "Reproduce the
      proof") is pending.
      *Status (2026-10-12): done on Joseph's machine with the post-review checker and `--E 1/2`:
      `results/cert16_check_output.txt` and `results/cert16_twist_output.txt` both end
      `PROVED: unique periodic orbit in X at E = 1/2`; the area-preservation guards passed; all
      certified intervals are the same as before the review fixes, now printed outward-rounded
      (e.g. HI ∈ [−0.96096564804637046, −0.96088055206183126], ghost trace ∈ [5.3469864, 5.3471379],
      τ ∈ [−3.0474820853652777×10⁵, −1.2280792906437143×10⁵]).*
      *Status (2026-10-12, later): an independent review (out of project context, 12 Oct) found
      that nothing checked p_y(T) > 0. F = 0 and energy conservation give only p_y(T) = ±p_y0, and
      this orbit really has the reflected zero: a float integration from the certified point
      reaches y = 0 at t = 15.52667 = T/2 with (x, p_x) restored and p_y = −0.999998. A
      certificate converging there would have "proved" verdicts about the half-period map. New
      lemma RET (`lemma_full_return`): the validated final enclosure Zf must have p_y > 0. For
      cert16, p_y on Zf ∈ [0.9999983922, 0.9999983923]: the certified T is the full period and the
      verdicts are about the monodromy; regression test 5 in `test_checker.py`. Also recorded:
      T is the 9th crossing of the global section {y = 0, p_y > 0}, so the "return map" of the
      twist is the local Poincaré map of the periodic orbit (time ≈ T), as Moser's theorem
      needs, not the first-return map of that global section. A rerun of the `--ghost` check
      with this checker (~18 min) would put RET into the committed output.* Not resolved by review B: an independent float fit of the return-map
      coefficients (described in its findings file).*
      *Status (2026-10-12, latest): the `--ghost` rerun is done on Joseph's machine (Spyder kernel
      interpreter; `cert_check.py results/cert16.json --ghost --E 1/2 --verbose`, checker with the
      RET print line): `results/cert16_check_output.txt` now shows `RET: p_y at time T in
      [9.9999839222565711e-1, 9.9999839222912168e-1] > 0` and ends `PROVED`; every other printed
      number is unchanged. `results/cert16_twist_output.txt` predates RET and was not rerun: RET is a
      property of the certificate (chain B's final enclosure Zf), which the twist layer reads
      unchanged, so the RET line above covers the twist verdict too.*
- [x] Off-diagonal ghost on the DR orbit (everything above is the colour-diagonal ghost).
      *Status (2026-10-12): `code/dr_offdiag_floquet.py`, output `results/dr_offdiag_floquet.txt`
      (float, DOP853 rtol 1e-12; run in the cloud session of 12 Oct, Python 3.13). Monodromy of the three
      8-dim off-diagonal blocks of ym/rung1_system.pkl (c_j = 1, the theory) over one period of the DR
      orbit scaled to energy E. Checks: orbit closure ~3e-13, det M = 1 to 1e-10, ghost multipliers
      in (mu, 1/mu) pairs, and with c_j = 0 the block reproduces the scalar Hill equations of the
      paper's off-diagonal section. Result: the off-diagonal ghost is Floquet-UNSTABLE on the orbit
      over about half of E = 0.05-1 (unstable intervals in the output), including E ∈ [0.281, 0.494]:
      at E = 0.3 block (u23, u32) has the real multiplier 1.5522 (rate 0.0125, against the sea-averaged
      off-diagonal g1 = 3.37e-3 of the paper's off-diagonal table), while the diagonal ghost is
      Floquet-stable there. So the E = 0.3 island, a non-escaping hole for the diagonal ghost at
      ε = 1e-10, is not a refuge in the full homogeneous theory, at least at its centre. The 12-d
      `hdym_island.py` runs have no off-diagonal components. Not certified; this is the central
      orbit only, not the surrounding tori.*
      *Status (2026-10-19): typed. `code/dr_offdiag_floquet_typed.py`, output
      `results/dr_offdiag_floquet_typed.txt` (B001 of the bridge work order), reproduces the 13 unstable
      intervals of `dr_offdiag_floquet.txt` exactly and gives each unstable multiplier its type: a real
      pair leaving the unit circle at +1 is parametric; a complex quadruplet is a Krein collision (two
      unit-circle pairs of opposite Krein signature meeting). Scope: DR orbit, off-diagonal ghost,
      c_j = 1, E ∈ [0.05, 1], 400 log-spaced points. Result: the type never changes inside a
      (block, interval). Krein quadruplets occur in block 12 only, E ∈ [0.0761, 0.0773], [0.1519, 0.1613]
      and [0.3469, 0.4937] (60 cells; max |mu| = 2.87, rate 0.032, at E = 0.412). Every other unstable
      cell (137) is parametric with the pair leaving at +1; no pair ever leaves at −1. Every unstable
      eigenvector satisfies the linearised Gauss law to 5e-11. The pure-YM off-diagonal block is stable
      at every E (so the E = 0.3 block-23 instability quoted above is the ghost's). Checks: closure 4e-13,
      det M − 1 ≤ 1e-9, (mu, 1/mu) pairing 1e-9, c_j = 0 Hill reduction 1e-11, no multiplier masked by
      the "drop the 4 nearest 1" rule. In the Broucke (A, B) diagram of the ghost pair, block 12 rides the
      seam between the Krein parabola and the period-doubling line from E ≈ 0.07 to 0.49 and always exits
      through the parabola; blocks 13 and 23 exit through the +1 line only (figures in the 19 Oct review
      brief, not committed). Open (B005, work order v3): the Krein signature of each pair, measured
      rather than inferred from pass/collide.*
      *Status (2026-10-19, B002): `code/channel_offdiag_monodromy.py`, output
      `results/channel_offdiag_monodromy.txt`: monodromy of each block over one transverse period
      T = 2π/D on the frozen-depth channel background x = D, px = 0, y = A cos Dt, z = A cos(Dt + ψ),
      A = √(2J/D), for D ∈ [0.5, 6] (24 log points), J ∈ [0.005, 0.5] (6), ψ ∈ {0, π/2}, c_j ∈ {1, 0}.
      The background is not an exact solution: the depth changes by ε = 4π²J/D³ (relative) over one
      period, cells with ε > 0.1 are marked OUTSIDE-MODEL, the gauge zero mode's Jordan pair splits by
      ~√ε there, so the verdict uses all 8 multipliers and the Gauss column tags growth in the broken
      constraint sector; the exported K has the background equations of motion substituted, so the
      c_j = 0 Hill check holds to O(ε) (9e-6 at ε = 3e-3). Result, scoped to ε ≤ 0.1 (188 of 576 cases):
      the c_j = 1 off-diagonal ghost is Floquet-stable in every block at both phases except one weak
      Krein cell (D = 3.89, J = 0.0315, |mu| = 1.002, rate 0.001, Gauss-respecting). Every strong
      instability lies outside the model and most of it violates Gauss (the F118 trap). The DR orbit
      itself has ε ≈ 0.2, so its off-diagonal instabilities, both types, belong to the channel bounce,
      not to the adiabatic channel interior. One scan, per the cap; not tuned.*
- [ ] Nonlinear off-diagonal escape map (`code/offdiag_escape.py`; the "all nine components" run that the
      off-diagonal item above says `hdym_island.py` lacks). *Status (2026-10-19): code written and checked,
      not yet run at size. The model is `ym/hdym_full_derive.py` (full nonlinear homogeneous SU(2) Lee-Wick,
      temporal gauge; its checks: A'''' coefficient −I, diagonal reduction, linearisation = rung1_system.pkl's
      K to 1e-15, Gauss law and Ostrogradsky energy conserved to 1e-18, generated C = sympy to 1e-15). The
      kernel's integrator agrees with scipy DOP853 to 6e-9 over two periods and in escape time to the stated
      tolerance (CHECK lines in the log). Smoke run in the cloud (6 energies x 8 seeds, both seeds 1e-6, tmax 1500):
      at E = 1 the diagonal ghost (Floquet rate 0.233 on the orbit) escapes at t ≈ 96 against ln(R/δ)/rate = 99;
      at E = 0.30, where the diagonal ghost is Floquet-stable, block 23's parametric pair drives the escape at
      t ≈ 950 (linear 1100–1840); E = 0.05, 0.09 no escape by 1500. An unseeded ghost-free orbit at E = 1 escapes
      at t ≈ 120 from the integrator's 1e-10 truncation seed, which is why both sectors are seeded at equal
      amplitude; and the escape time converges in tolerance only for seeds ≳ 1e-5 (convergence table in the
      docstring), so the defaults are δ = δ_diag = 1e-4, rtol 1e-11. Near a Krein band edge the two block-12
      modes beat and amplify the seed ~10³ while Floquet-stable (E = 0.05: 1e-4 → 0.23 and back, energy
      conserved), so the sector that leads a runaway is read as "largest before escape", not "first past a
      threshold". Committed logs: none yet; the sized runs are Joseph's (commands in the docstring).*
      *Status (2026-10-20): section E run on Joseph's machine (100 energies x 64 seeds, both seeds 1e-4, tmax 3000,
      rtol 1e-11, 24 threads, 16893 s; log `results/esc_E.txt`, data `results/esc_E.npz`; all checks PASS, no
      integrator failures). Reading: the nonlinear escape times lie between ln(1/δ)/rate and ln(R/δ)/rate of
      the leading linear instability at every E. Where the diagonal ghost is Floquet-unstable on the orbit
      (E ≳ 0.31 and its tongues at 0.09, 0.16–0.19, 0.24–0.26) it leads, and above E ≈ 0.33 its rate (0.04–0.23)
      exceeds every off-diagonal rate, so block 12's Krein band at 0.35–0.49 (rate ≤ 0.032) never leads an escape
      there. In the gaps where the diagonal ghost is stable the off-diagonal instabilities take over at their own
      linear rates — E = 0.153–0.158 (block 12 Krein, median 2276 vs 1460–2900), 0.227–0.246 (block 12
      parametric, 914 vs 575–1150), 0.28–0.31 (block 23 parametric) — and in each the off-diagonal mode grows
      first and the diagonal ghost, pumped nonlinearly, overtakes and does the escaping. Where nothing is
      linearly unstable (E ≈ 0.113–0.149, 0.195–0.22, 0.272) nothing escapes by 3000. The orient map at
      E = 0.412 (128², `results/esc_orient_0412.*`, 32975 s) is structureless (median 183, 96 % of cells escape
      through the diagonal ghost): at that E the diagonal ghost wins regardless of the off-diagonal seed, and
      the per-cell random diagonal seed decorrelates the map; the script now holds one diagonal seed per map
      (`--diag-orient fixed`) and the maps worth making are at E = 0.156 and 0.30, where the off-diagonal sector
      leads. Decomposing the escape time with the checkpoints: the diagonal ghost first reaches amplitude 1 at
      ln(1/δ)/rate within the 50-unit checkpoint resolution at every E from 0.09 to 0.9 (the linear rate is the
      clock), and what follows is a bounded nonlinear loiter of median 10–170 and 90th percentile 50–630 time
      units (zero to ~ten orbit periods) before the blow-up — negligible against a linear phase of thousands at
      low E, comparable to it above E ≈ 0.5, which is why the packets there spread upward. The upper line
      ln(R/δ)/rate, which assumes exponential growth to R, is not the right upper edge; the end is a loiter then
      a sharp blow-up.*
      *Status (2026-10-20, later): a thread-dependent floating-point mode was found and fixed. The 0156 orient run's
      θ = 0 row (128 bit-identical cells) gave two escape times, 2965.95 for 123 cells and 1818.84 for 5 (one thread's
      share of 24): on Windows numpy/MKL sets flush-to-zero on the main thread, OpenMP workers start with the default
      mode, a denormal in the O(δ²) z-sector of the planar orbit was flushed on one thread only, and the chaotic
      loiter made the last bit into 1100 time units. `hdym_full_escape.c` now pins MXCSR = 0x1F80 (FTZ/DAZ off, the
      Linux mode) in every thread and the log reports the modes found. Consequence for the committed section-E
      run: about 1/24 of its cells were integrated in the other mode; their escape times are valid trajectories of
      the same equations up to a last-bit difference in the seed, so the statistics and every statement above
      stand, but cell-level reproducibility between platforms holds only for runs after this fix. The 0156 and 030
      orient maps are rerun with the fix.*
- [ ] ε-crossover scan at E = 0.3: `hdym_island.py --eps 1e-9 … 1e-6`.
- [ ] 3-dof: the DR orbit's z-direction Floquet trace (−1.53, float) is not certified; the
      island is sealed in 2-dof and leaks by Arnold diffusion in 3-dof (seen: 2.7 % at r = 2e-4).

## Conventions

`hdym_island.py` uses the repo's conventions verbatim: state (p, w, P, W), q = p − w,
E = |P|²/2 + V(q), ghost |w|² + |W|² = 2εE, p = q + w, P = P_q − W, Yoshida-4 with dt = 0.005
and channel sub-stepping, escape at any |component| > 10⁴. The DR point in these
coordinates: q* = λ(3.14640122769837, 0, 0), P* = λ²(0.00179319341881, √(1 − p_x*²), 0),
λ = (E/½)^{1/4}; copies 0–23 by cyclic permutation, 90° in-plane rotation, time reversal.
*Status (2026-10-12): "verbatim" except the channel sub-stepping: `hdym_island.py` runs pairs
through `hdym_precision.CPairs` with the SHARED sub-step schedule (both members take the larger
n), not the production per-trajectory schedule of `hdym_escape_map`; its docstring says so. The
difference is at truncation level.*
