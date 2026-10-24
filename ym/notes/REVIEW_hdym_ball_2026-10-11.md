# Review of `ym/hdym_ball.py` (Prompt A), 2026-10-11

Record of an independent review (fresh session, no prior contact with the repository), 2026-10-11; kept as written.
F1–F6 were fixed on 11 Oct and the production sweep rerun on 12 Oct with identical certified fractions; F7 (QR
containment) was answered on 13 Oct by the degradation test of `NOTE_precision_tier.md` — `ym/README.md`,
"Certified pair escape times". Repository cloned read-only (shallow clone at HEAD, 2026-10-11 00:58 EDT). Budget 20 min wall clock.
Line numbers refer to `ym/hdym_ball.py` at that HEAD unless stated.

Selftest run here: `python hdym_ball.py --selftest` passed (selftest 0: all six polynomial
bounds ok, K2 P row with slack 8(Q0+Q1+Q2); selftest 1: worst dd error 0.029 x 2^-101;
selftest 2: three 60-digit orbits inside the bound, worst ratio ~1e-4 to 1e-5). sympy and mpmath
had to be installed in this environment (not a repository issue).

Overall: no error found that breaks the main guarantee. The ball rules, `step_remainder`, the
tangent balls and the ellipsoid update are sound as far as checked. The findings are one
edge-case logic flaw in the escape bookkeeping (F1), an untested path (F2), and overstated or
loose wording (F3-F6).

## Findings

### F1 — breaks the guarantee (edge case): a lost member b can still be recorded as certainly escaped
- Where: `bpair`, line 593 sets `out[3] = -2` on bit 4; line 599 `if (eb) out[3 * N + i] = t;`
  then overwrites it in the same step. Related: line 608, the return value of `ell_update` is
  ignored at the member switch.
- Must guarantee: a member flagged lost (bit 4: NaN, `E_k > R/2`, or `ell_update` returned -1)
  never contributes a certified exit time.
- What it does: in mode 2, if member b has bit 4 and bit 2 in the same step, `classify` sees a
  certified exit time for b. If bit 4 came only from `E_k > R/2`, E is still a valid bound and
  nothing is wrong. If it came from `ell_update` returning -1 (line 391: `L` already overwritten
  with the unchecked `L'`, and `*rho` not yet updated), E is not a proven bound, so the recorded
  exit time is not certified. NaN is safe (NaN comparisons cannot set bit 2). At line 608 a -1
  return leaves b's frame unverified with no flag at all.
- Reachability: very unlikely (a zero diagonal in R although M contains `sb I`, sb > 0).
- Fix: `if (eb && !lb) out[3 * N + i] = t;`; at line 608, `if (ell_update(...)) { mark b lost;
  break; }`.

### F2 — gap that needs an argument (test coverage): the difference-ball (rd) path has no end-to-end test
- Where: `selftest`, lines ~896-913.
- Must guarantee: b_true in (v + d) +- (rv + rd): the D rules in `badd/bsub/bmul/bscl` and
  `def[1]` hold in practice, not only on paper.
- What it does: selftest 2 compares only member a (output rows hold a's centre and bound). The rd
  propagation, which is the reason the code exists (certifying delta << 1e-16), is never
  compared with an independent orbit. With t = 40 and R = 1e4 nothing escapes, so the escape
  bookkeeping, the mode switch (2 -> 1 / 0) and the J = I fold never execute.
- I checked the D rules by hand (see below) and found them correct; this is a guard against
  future edits.
- Fix: also export b's final centre and bound, compare with a 60-digit orbit started from a + D
  (exact rational sum) at delta ~ 1e-30; add a short run with small R so the switch path runs.

### F3 — gap that needs an argument (wording): "rigorous bracket on the uncertainty exponent"
- Where: line 1007; the module docstring "f_lo = P(certified different) <= f_true <= ...".
- f_lo <= f_true <= f_hi holds rigorously for the empirical fraction over the M sampled pairs
  (given F1). It does not hold for the probability (Monte Carlo error ~ 1/sqrt(M)). Slopes fitted
  to f_lo(delta) and f_hi(delta) do not bound the slope of f_true: f_true can wander between the
  curves, so its fitted slope can fall outside [alpha(f_lo), alpha(f_hi)], and each alpha has a
  statistical error.
- Observed: the small sweep run here (`--M 32 --delta-max 1e-6 --delta-min 1e-10 --per-decade
  0.5 --tmax 300`, output in `ball_sweep_small.txt`) prints, under that heading, for |dt| > 10:
  alpha from f_lo 0.092 +- 0.034, from f_hi 0.048 +- 0.013, i.e. the "lower" curve gives the
  larger exponent. The two numbers are not an ordered bracket.
- Fix: "empirical fraction over the sampled pairs" in the docstring; print "alpha fitted to the
  lower / upper curve". A rigorous claim about alpha needs the uncertified fraction to be small
  compared with f itself at every delta used in the fit.

### F4 — cosmetic (semantics): the certified map depends on the partner
- Where: module docstring "What is certified"; line 553 (`n = max(na, nb)`).
- In mode 2 both members are stepped with the same `n = max(n_a, n_b)`, so b's sub-step schedule
  depends on a; after a switch b uses its own n. So what is certified is the escape time of a
  pair-dependent, non-autonomous step schedule, not of one fixed map applied to each initial
  condition separately. Both members always share the schedule, so the comparison is fair; only
  the wording ("the escape-time function of this symplectic map") is loose.

### F5 — cosmetic: "exact differences D = S_b - S_a"
- Where: `ym/hdym_diffpair.py:334` (`sample_dpairs` docstring).
- D is computed in floating point (sines, square roots): D = D_intended (1 + O(1e-16)). The
  certificate is for the pair (a, a + D) with D the stored doubles, which is fine; the separation
  is delta (1 + O(1e-16)), not exactly delta. Suggest "D computed to about 1 ulp of D".

### F6 — cosmetic: the default escape observable differs from the paper's
- `--esc ghost` (default) certifies max(|w|, |W|) > 10; the production landscape uses the max over
  all 12 components. The help text says so; worth one line wherever a bracket is compared with
  the paper's landscape.

### F7 — cosmetic (observation, not an error): the QR containment radius dominates in longer runs
- In selftest 2 (t = 40) the diagnostic rho'/|Y|_F is ~1e-12; in the small sweep (t <= 300) it
  reaches 3e-5 to 8e-4, while the tangent radius stays ~1e-14. So (|N|_2 - 1)|L'|_F, not the
  tangent arithmetic, is what inflates the enclosure, and it likely drives the "lost" fraction
  (0.19-0.34 here). Not a soundness issue; a hint for tightening (e.g. one step of iterative
  refinement of N, or bounding |N|_2 - 1 via |N N^T - I| instead of Gershgorin on N N^T).

## Checked and found sound

- Independent test of `badd`, `bsub`, `bmul`, `bscl` (scratch harness, not in the repo): the
  kernel functions compiled from `c_src()` with `-ffp-contract=off`, 200 000 random balls (values
  1e-3..1e3, differences 1e-40..1e-5, relative radii 1e-32..1e-14), true points at random ball
  corners, results compared in exact rational arithmetic: 0 violations, worst |error|/radius =
  1 - 1.4e-14 (tight, as expected for corners of linear operations).

- `badd`, `bsub` (line 152 ff.): AccurateDWPlusDW error <= 3u^2/(1-4u) |exact|, charged as
  2^-101 |computed| (large slack); `ddabs` >= |h + l| even after its own rounding (1 + 2.3e-16
  rounds to 1 + 2^-52, |l| <= u|h|); 3-term radius sums covered by INFL; underflow by ETA.
  D(a + b) = Da + Db.
- `bmul` (line 164): D(ab) = Da b + a Db + Da Db; all nine ball cross terms present; rounding of
  p1, p2, p3, s, d each charged; the ~14-term radius sum (gamma_15 ~ 1.7e-15) covered by INFL
  (1.4e-14). DWTimesDW (FMA) error <= 4u^2 < 2^-101.
- `bscl`: c an exact double. The stage constants fl(YC*h), fl(YD*h), h = fl(dt/n) are the same
  expressions in `yoshida_b`, `yoshida_t` and `step_remainder`.
- `step_remainder` (line 460): magnitude propagation matches the drift/kick updates (q = p - w
  moves by c(P + W)); FP, FW, K1, K2 confirmed by sympy (the P row is W row - 1 because the
  Hess V terms cancel under d/dw at fixed p); composition bound |DPhi(u) - DPhi(v)| <=
  G^2 sum a_i K2_i |u - v| holds since every l_j >= 1 and the stage bounds are boxes (convex), so
  segments between image points stay inside; the remainder M2 r0 (sh + r0/2) follows from the
  integral form; x~ and x + s e lie in S0; same n and stage order as `yoshida_b`; final rounding
  covered by INFL^3 + ETA.
- `tangent_step` and the `tb` rules (line 232 ff., 295): tangent of a product t_a b + a t_b with
  all six cross radius terms; exact inputs with zero radii; same constants and n. Encloses
  J(x~) col.
- `ell_update` (line 335): (1) |(JL - Y)u| <= |RY|_F; carry = sqrt(12) G rho valid since |J|_2 <=
  sqrt(12)|J|_inf <= sqrt(12) G with x~ in S0. (2) Minkowski: with w = (u/sa, b/sb), |w|^2 <=
  c/(1+c) + 1/(1+c) = 1; sa, sb rounded up by FNU^2; |D|_F <= U2 |fl(sa Y)|_F. (3) Nothing
  assumed of the QR. (4) M = L'N + E gives {Mu} subset {L'v : |v| <= |N|_2} + ball(|E|_F)
  subset {L'v : |v| <= 1} + ball((|N|_2 - 1)|L'|_F + |E|_F); Gershgorin on N N^T; gamma_13 =
  1.44e-15 < 2e-15 and gamma_24 = 2.66e-15 < 3e-15 for the sums actually formed (13 and 24
  terms); `sqrt_up` covers sums of <= 288 squares. The tr1 == 0 branch is right. Underflow of
  squares would only matter if all entries were below ~1e-154, impossible here (errors are at
  least ~1e-31 times the state).
- `ebound`: Cauchy-Schwarz row bound plus rho, rounded up.
- Error recursion in `bpair`: e_{s+1} = DPhi(x~) e_s + remainder + defect, with the defect from
  one-step balls started at the exact double-word centre (rv = rd = 0); member b's defect rv + rd
  and centre v + d; base point the hi part with shift 2U|x~|; Eb, emax, rho fed to
  `step_remainder` and `carry` are start-of-step values; the switch fold with J = I and carry =
  rho is right (apart from F1).
- `estatus` (line 511): possible (|c| + E > R) / certain (|c| - E > R) with E including the
  hi-part shift; exit step in [first possible, first certain] follows.
- `classify` (line 743) and the sweep: diff (blo - ahi > T or alo - bhi > T) and same (union span
  <= T) are correct implications from valid intervals; lost (-2) and unescaped (-1) pairs are
  uncertified; f_lo, f_hi use the ok-mask denominator as described.
- Build flags: `-ffp-contract=off -fno-fast-math` (C) and `--fmad=false` (CUDA) match the rounding
  assumptions; explicit `fma()` calls are intended.

## Not reached / not done
- The small sweep was run once (49 s; output `REVIEW_hdym_ball_2026-10-11_sweep_small.txt` beside this file), but it was
  not used to test the difference ball end to end (F2 describes that test). The default sweep
  (`--M 2000 --delta-min 1e-40 --tmax 3000`) was not run (far over budget).
- `nsub_b` and the CuPy path were read, not run (no GPU here).
- `hdym_escape_map` (`base_point`, `default_width`) and `hdym_uncertainty` (`sample_pairs`,
  `fit_alpha`) were not reviewed beyond the data flow.
