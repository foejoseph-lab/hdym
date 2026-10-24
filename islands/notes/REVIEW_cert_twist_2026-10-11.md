# Review of `islands/code/cert_twist.py` (Prompt B), 2026-10-11

Record of an independent cold review (no prior contact with the repository), 2026-10-11; kept as written. Findings
F1–F5 were fixed on 11 Oct; the independent float fit of the return-map coefficients is still open —
`islands/README.md`, "Review history of the checker". Repo cloned read-only at HEAD of 2026-10-11 01:39 EDT. Budget 20 min; used ~11 min. Status: COMPLETE for the items listed under "Reached".

Verdict: no finding that breaks the proof. The twist layer, the formula and the glue are sound
as far as I could check; F1-F3 are gaps in what the checker enforces/declares, F4-F5 cosmetic.

## Findings

### F1 — gap that needs an argument (output/glue): twist verdict printed before Krawczyk is known
`cert_check.py:122-159`. The twist enclosure is valid only for the TRUE fixed point (x*, px*)
and TRUE period T*; the env uses x_s = X[0], px_s = X[1], T in [T_lo, T_hi], i.e. it relies on
(x*, px*, T*) in X, which is the conclusion of KRA. The code computes `ok` at line 114 but does
not test it before running the twist block; the line "BIRKHOFF TWIST ... NONZERO (Moser twist
hypothesis certified)" is printed regardless, and only the last line says NOT PROVED. Same for
the "ELLIPTIC" verdict of lemma T. With the committed certificate `ok` is True, so the
committed result is not affected; but the printed claims are conditional.
Fix: `if not ok: raise CheckFailed(...)` right after `lemma_krawczyk`, before T / SYMP / TWIST.

### F2 — gap (robustness, not soundness): consistency checks are printed, not enforced
`cert_twist.py:335`, `cert_check.py:155-157`. At the true point det L = 1, the zetabar
coefficient G[(0,1)] = 0 and Re(conj(lambda) c1) = 0 (area preservation; cubic normal-form
coefficient is a conjugacy invariant). The checker prints the last two but does not require
the enclosures to contain these values. Soundness does not need them (each step is an interval
extension evaluated at the true data), but enforcing `0 in G01`, `0 in Re(conj(lambda) c1)`,
`1 in det L` costs nothing and would catch a wrong `twist_formula.json` or a wrong direction
convention (a cheap guard on the trust hole of finding F3).

### F3 — gap (trust boundary): `twist_formula.json` is trusted input
The checker evaluates whatever trees the JSON holds, so the formula (and its generator
`twist_symbolic2.py`, which uses sympy) is part of the trusted base; the README's list of what
is trusted (cert_lemmas, cert_check, cert_twist, Python ints, isqrt) omits it.
What I did: regenerated it with `twist_symbolic2.py` (sympy 1.x, 2026-10-11): byte-identical,
md5 64cc7b9c35a3096454605a4c3013c2d5; its float cross-check gives tau = -2.140067e5, as
`twist_symbolic.py`. I read the generator's derivation (see "checked" below) and found it
correct. Fix: list the JSON + generator + sympy in the trusted base in the README, record the
md5, and add the F2 guards.

### F4 — cosmetic: docstring on non-resonance
`cert_twist.py:17` says lambda^q != 1 for q = 1,2,3 "is enforced by the divisions". The
divisions exclude q = 1 (lambda - 1, lambda^2 - lambda) and q = 3 (lambda^2 - conj lambda);
q = 2 and q = 4 are excluded by the explicit cos theta checks at :295-296. The code is right,
the comment is not.

### F5 — cosmetic: silent zero in `A_series`
`cert_twist.py:40-41` (and `variational_coeffs` in cert_lemmas) uses `I(0)` when `k >= len(yy)`.
`taylor_coeffs(., q)` returns yy of length q, so A_q is silently 0. It is never used (the
coefficient n+1 needs A_0..A_n, n <= order-1), so no harm, but an explicit
`raise CheckFailed` would make the invariant checked rather than assumed.

## Checked and found sound

- Math, item 1. Moser's twist theorem needs: analytic area-preserving map, elliptic fixed
  point, lambda^q != 1 for q = 1..4, first Birkhoff coefficient != 0. cos theta not in
  {1,-1,-1/2,0} is exactly q <= 4 for |lambda| = 1. Section {y=0, py>0, H=E} with (x, px) is
  canonical (omega restricts to dx^dpx). tau != 0 does not depend on the section coordinates:
  the cubic normal-form coefficient is invariant under near-identity conjugacy and only scales
  by |alpha|^2 under linear changes commuting with the rotation.
- `lemma_twist` canonical frame (:297-304): v = (b, lambda - a) = u + i w is an eigenvector of
  L; S = [u w]/sqrt(area) has det 1; the sign branch is taken only when the area interval
  excludes 0. Each quantity is an interval extension at the true L, so the true S is enclosed.
  zeta' = conj-rotation of zeta; lambda is read from G[(1,0)], so the convention is consistent.
- Kuznetsov c1 (:328-333): term-by-term Kuznetsov's Neimark–Sacker formula with
  g_jk = j!k! [zeta^j zetabar^k]; for z -> lambda z e^{i tau |z|^2}, c1 = i lambda tau, so
  tau = Im(conj(lambda) c1). `cdiv` checks |b|^2 > 0 before dividing.
- Generator `twist_symbolic2.py`: lift delta = xi1 e_x + xi2 e_2 + (c2 xi2^2 + c3 xi2^3) e_py from
  the sqrt series (c1 = -px/py, c2 = -E/py^3, c3 = -px E/py^5 since px^2+py^2 = 2E); third-order
  Taylor of Phi_T in delta, keeping exactly the terms of order <= 3 (e_py appears only via
  O(xi^2), so only B(.,e_py) with linear partners and no C with e_py); return-time correction by
  3 fixed-point iterations of y(T+dt) = 0 (each gains one order); f2 = Df f, f3 = Df f2.
  Constant terms dropped use Phi_T*(z*) = z*, true at the true data, which are enclosed.
- A2 / A3 (:54-67, :105-112, :183-190): Hessian / third derivative of f = (px, py, -x y^2, -x^2 y)
  checked by hand; the c-equation c' = Ac + A2(m_i,b_jl) + A2(m_j,b_il) + A2(m_l,b_ij) + A3 is
  the standard third variation; TRIPS are sorted so the PAIR_IDX lookups exist.
- `dir_series`: coefficient n+1 = (rhs)_n/(n+1) for m, b_zero, b_full, c; product caches are
  extended before use and only read coefficients <= n.
- `lemma_dir_step`: Taylor part on B with m(0) = m_k, b_full(0) = b_k, b_zero(0) = c(0) = 0;
  remainder = (p+1)-th coefficient of the autonomous augmented polynomial ODE about xi, evaluated
  with z in Z (PIC, verified by lemma_step) and m, b_zero, b_full, c in their Picard enclosures
  over [0, h] (`_dir_rough`). Argument order into `dir_series` (b_full0 = bf_r, b_zero0 = bz_r)
  matches the signature.
- `_dir_picard_tools.picard`: linear inhomogeneous Picard lemma; returned vn = v0 + [0,d](A(Z) vt
  + F) is a valid (tighter) enclosure because the integral of a function with values in a box
  over [0, t] lies in t * box.
- FRAME2 (`DirState`, run_chain:344-349): b(t_{k+1}) = J_k b(t_k) + Delta_k (linearity in b(0));
  with J_k in A_{k+1} T_k A_k^{-1} (T_k = Ts[0] from lemma_step, J over B which contains the
  true z(t_k)), beta_{k+1} in T_k beta_k + Ainv Delta_k. m_k = A_k R_k m0 uses R_k before its
  update; A_0 = I is checked, so beta_0 = gamma_0 = 0 is right (D^2 phi_0 = 0).
- `lemma_dir_time_interval` + glue: z(t) in Zf on [T_lo, T_hi] by PIC (lemma_rough_time_interval),
  m in Jfinal m0, b/c by Picard from b_fin/c_fin. m0 = (e_x, e_2, e_py) with e_2 containing
  (0,0,1,-px*/py*) since px* in X[1]; env x_s, px_s, py_s, E enclose the true values (given F1).

## Targeted tests run (scratch only, repository untouched)

1. Regenerated `twist_formula.json` with `twist_symbolic2.py`: byte-identical (F3).
2. `lemma_twist` on a map with KNOWN twist: z -> lambda z (1 + i tau0 |z|^2), lambda = (3+4i)/5,
   tau0 = 7/3, put in a non-orthogonal symplectic frame S0 = [[2,1],[1,1]] and conjugated by a
   non-symplectic near-identity map h (quadratic + cubic terms), coefficients fed as exact
   rationals. Result tau = 7/3 exactly (and -7/3 for the reflected orientation, with lambda
   conjugated: the sign of tau follows the orientation, nonzero-ness does not); Re(conj(lambda)
   c1) encloses 0. Confirms the frame, the complexification, Kuznetsov's formula and the
   invariance argument.
3. One `lemma_dir_step` (h = 1/20, p = 12, point initial data at the DR orbit, Z = z0 +- 0.2
   checked by PIC) vs. a float DOP853 integration of the same augmented ODE (rtol 1e-13):
   all 20 b- and 16 c-components inside the enclosures (widths ~1e-16).
4. First 20 steps of chain B of `cert16.json` through `run_chain` with `DirState` (VAR2/VAR3 +
   FRAME2 + the real certificate frames), compared with a float integration from the box
   centre: all m, b, c components inside (widths b 1.7e-10, c 9e-11). 94 s.
5. Committed `cert16_twist_output.txt`: Re(conj(lambda) c1) in [-1.14e5, 1.15e5] contains 0,
   cos theta matches the monodromy HI, the float tau = -2.14e5 lies inside the certified
   interval, KRA True.

## Not resolved

- Independent float check of the return-map Taylor coefficients (P1, P2) by a least-squares fit
  of the numerically integrated Poincare map (solve_ivp with a y = 0 event) against the formula:
  my event detection returned no crossing twice, so I stopped this line (rule 3). For the
  author (a few minutes): sample P on a 9x9 grid of radius 2e-3 around (x_s, px_s), locate
  the return with dense output + root finding on y near t = T, fit monomials of degree <= 5,
  compare degrees 1-3 with the formula evaluated as in `twist_symbolic2.py` (redirect the output
  to a file). Test 2 above covers lemma_twist; this would cover the generator independently.

## Reached / not reached

Reached: items 1 (mathematics), 2 (provenance; regenerated and diffed, derivation read),
3 (dir_series, lemma_dir_step; tests 3-4), 4 (DirState/FRAME2, lemma_dir_time_interval;
test 4 covers FRAME2 over 20 steps only), 5 (lemma_twist; test 2), 6 (glue; F1).
Not reached: a full `--twist` run (~50 min, not allowed); FRAME2 over the whole chain and
`lemma_dir_time_interval` were checked by reading only, not by test. The base lemmas of
cert_lemmas.py (lemma_step, PIC, PICV, KRA, ...) were taken as given (Prompt C's scope).
To exercise the end-to-end path, the author can run:
`cd islands/code && python3 cert_check.py ../results/cert16.json --ghost --twist twist_formula.json > ../results/cert16_twist_recheck.txt`

## Addendum (written after the review above was complete): effect of the Prompt C findings

Source: `REVIEW_cert_lemmas_2026-10-11.md` (independent session, 01:24-01:30 EDT). I read it
only after finishing everything above; nothing above was changed.

Closing the dependency: this review took the base lemmas as given. C found these sound:
interval primitives, `isqrt`, `_div`, TAY, VAR, PIC, PICV, LAG, MVT, NEU, the base part of
`lemma_step`, FRAME for m = 0, `lemma_rough_time_interval`, KRA with F and DF. These are
exactly the pieces the twist layer uses: Z and B from `lemma_step`, Xt from PICV, Ts[0] and
Ainv from the base frame, Rs[0], Jfinal and Zf over [T_lo, T_hi], and KRA. With C's verdicts
the twist result rests on reviewed code throughout.

C's findings, one by one:

- C-F1 (ghost frames switch per step; breaks `--ghost` for an adversarial certificate):
  does NOT affect the twist. FRAME2 uses only the base frames (Ts[0], `ex["Ainv"]`, Rs[0],
  `final["A"]`), never G. The fix (all steps have "G" or none) does not touch the twist path.
- C-F2 (initial inclusion checked against the outward-rounded box): affects the twist in the
  same way as the base flow, because the directional systems ride on chain B. Not triggered
  by cert16 (C checked that steps[0] lies on the grid). C's fix closes it for the twist too.
- C-F3 (printed endpoints rounded to nearest): applies to the twist lines too (tau, cos theta,
  lambda). The decisions themselves are made on Fractions (`tau[1] < 0`, the resonance and
  ellipticity tests in `lemma_twist`), so only the quoted digits are affected:
  tau in [-304748.21, -122807.93] could be off in the last digit. Same fix.
- C-F4 = my F1, found independently (verdicts printed regardless of Krawczyk). C rated it
  cosmetic, I rated it a gap; the fix is the same. For the twist it matters slightly more,
  because the twist evaluates the formula at X and [T_lo, T_hi], which is meaningful only if
  the true orbit is in there.
- C-F5 (no lo <= hi check at parse): the twist path reads no new certificate fields (it uses
  the step data and Zfinal, which C covered). Not affected.
- C-F6 (E is read from the certificate, never compared with 1/2): applies to the twist too;
  E enters the env and the formula (c2, c3). cert16 has E = 1/2 (C checked). Same fix.
- C-F7 ("stable" for Lemma T): no effect. `lemma_twist` re-derives ellipticity from the
  return map L itself (|cos theta| < 1), so Moser's hypotheses do not depend on Lemma T.
- C's "not reached" (the certified orbit is the DR orbit): applies equally to the twist claim.

Net: no new finding against the twist layer; C-F2 and C-F6 carry over to it and are not
triggered by cert16; C-F1 does not reach it.
