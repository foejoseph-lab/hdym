# Review: `islands/code/cert_lemmas.py` and `islands/code/cert_check.py`

Reviewer: independent Claude session (no prior contact with the repository), 2026-10-11.
Clone: `https://github.com/foejoseph-lab/hdym`, shallow clone of the default branch at 01:24 EDT.
Budget: 20 minutes wall clock, started 01:24 EDT.

Status: COMPLETE (about 6 minutes of the 20-minute budget used).

## Findings

### F1. Ghost frames: the frame used can change between steps, so the accumulated ghost Jacobian need not telescope. Severity: breaks the proof (`--ghost` conclusion), for an adversarial certificate.

- Location: `cert_lemmas.py:238-243` (`lemma_step`), `cert_lemmas.py:334-339` and `:356` (`run_chain`).
- What the lemma must say (FRAME): M = F_N R_N with R_{k+1} = (F_{k+1}^{-1} J_k F_k) R_k encloses J_{N-1} ... J_0 **only if the same frame F_k appears on the left in step k-1 and on the right in step k**, and F_0 = I.
- What the code does: per step, it uses the G frames if **both** `step` and `nxt` carry `"G"`, and otherwise it falls back to the A frames. The choice is made separately for each step. If the certificate has `"G"` on some steps and not on others, the product contains a factor such as G_k A_k^{-1} or A_k G_k^{-1}, which is not the identity. The final frame (`:356`) is G_N if the last step has `"G"`, otherwise A_N, independently of what the previous steps used. The identity check at `:334-339` covers G_0 only when steps[0] has `"G"`.
- Effect: a certificate that omits `"G"` on one step can insert an arbitrary matrix into the ghost monodromy and so control the SYMP verdict. The base flow (m = 0) always uses the A frames and is unaffected.
- Is cert16 affected: **no**. chainB has `"G"` on all 1001 steps, G_0 = I, and chainP has none and is run with m_list = [0] only. Checked with a script that reads the certificate. So the committed `--ghost` verdict stands, but the checker would not catch a bad certificate of this kind.
- Demonstrated: I took steps 0-3 of chainB of cert16 and deleted `"G"`/`"Ginv"` from step 2 (scratch script `rv/t1.py`; the repository was not touched). `run_chain(..., m_list=[0,1])` accepted both the original and the altered slice. The two ghost Jacobian enclosures differ by up to 6e-2 while their widths are about 1e-11, and some entries are disjoint, so at least one accepted enclosure does not contain the true value.
- Suggested fix: decide once per chain, from steps[0]. Either every step has `"G"` (with the right length) or none does, and raise CheckFailed otherwise. An alternative is to drop the fallback entirely.

### F2. The initial-data inclusion is checked against the outward-rounded box, not the set c + A r. Severity: gap that needs an argument (it closes when c and r lie on the 2^-K grid).

- Location: `cert_lemmas.py:326-329`.
- What the lemma must say: the initial data lie in the set S_0 = c_0 + A_0 r_0, because `lemma_step` proves phi_h(S_0) ⊂ S_1 for points of the form c + A rho, rho in r, and nothing else.
- What the code does: `set_box` rounds c + A r outward to 2^-K. With A_0 = I (checked a few lines later), the rounded box equals the set exactly when every c_i, r_i endpoint is a multiple of 2^-K. Doubles with |v| < 2^-148 are not. A point in the rounding slack is not covered by the step lemma.
- cert16: every c and r endpoint of steps[0] in both chains lies on the 2^-200 grid (checked), so the gap does not affect the committed result.
- Suggested fix: after the identity check, test `initial_box[i][0] - c[i][0] >= r[i][0]` and `initial_box[i][1] - c[i][1] <= r[i][1]` in exact Fractions, with no rounding. An alternative is to reject certificates whose numbers are off the grid.

### F3. Printed numbers are rounded to nearest, not outward. Severity: cosmetic.

- Location: `cert_check.py:127, 132-133, 137` (the `float(...)` and `!r` of interval endpoints, and `float(s[0])/2`).
- These are the numbers quoted in the README, for example HI in [-0.96097, -0.96088]. Converting a Fraction endpoint to the nearest double can move it inward by up to half an ulp, so a printed interval may fail to contain the true value at the last digit. The README's 5-digit rounding is outward in practice, but nothing guarantees it.
- Suggested fix: print `floor`/`ceil` to a stated number of digits from the Fractions, or print the exact Fractions alongside.

### F4. Verdicts are printed even when Krawczyk fails, and the docstring describes the ghost test differently from the code. Severity: cosmetic.

- `cert_check.py:130-137` print ELLIPTIC and UNSTABLE before and regardless of `ok`. Without the existence of an orbit, these say nothing about any orbit. Suggested fix: print them only if `ok`, or label them conditional.
- The docstring item (3) says "if |tr| > 2 on the in-plane block". The code uses SYMP (`lemma_symplectic_4x4`), which is the correct test. The docstring should be updated.

### F5. Certificate intervals are not validated (lo <= hi) at parse time. Severity: cosmetic, as far as I can see.

- `cert_check.py:51` (`piv`). I traced the effect of an inverted interval. An inverted `r` behaves as its hull everywhere it is used, because it only enters through `mul`. An inverted `r'` or `Z` can never pass `inside`. An inverted box is rejected by `inside(xhat, X)`. So this is not exploitable, but an explicit check would make the argument unnecessary.

### F6. The printed conclusion does not state the energy E, which is read from the certificate. Severity: cosmetic, bordering on a gap.

- Location: `cert_check.py:77, 85-86`.
- The claim is "at E = 1/2", but E comes from the untrusted certificate and is never printed or compared with 1/2. A certificate at a different E would print the same RESULT line. cert16 has `E = 0x1.0p-1` = 1/2 (checked), so the committed claim holds.
- Suggested fix: print E as an exact Fraction in the header, or accept `--E 1/2` and raise CheckFailed on a mismatch.

### F7. The word "stable" for Lemma T. Severity: cosmetic.

- `cert_lemmas.py:297` and `cert_check.py:33-34` say "linearly stable (elliptic)". The monodromy of an autonomous periodic orbit has the eigenvalue 1 in a Jordan block (period-energy shear), so it is not linearly stable in the Lyapunov sense. What is proved is that the linearised Poincaré return map is elliptic (spectrally stable), which is what the Moser layer needs. A wording change is enough.

### Independent numerical sanity check (non-rigorous)

Step 0 of chainB of cert16 was checked with the checker's own lemmas (PIC, LAG, MVT with m = 0 and with m = 1). It was compared with a 60-digit `mpmath.odefun` integration of the flow and of the variational equation from the step centre c. phi_h(c) lies inside the LAG enclosure (width 1.7e-27), and X(h; c) lies inside the MVT enclosure J (width 4e-12) for both the base flow and the ghost. Script: `rv/t3.py`, `rv/t3g.py` (scratch).

## Lemmas checked and found sound (one line each)

- `rnd`, `add`, `sub`, `neg`, `mul`, `scal`, `div_int`, `hull`, `inside`: floor and ceil to 2^-K are exact on Fractions, `neg` and `hull` are exact, and `mul` takes the min and max of the four endpoint products. A random test of 3000 cases, checking endpoints and interior points against exact Fractions, found 0 violations.
- `isqrt`: `int(n)` is the floor for n >= 0. The lower end floor(sqrt(floor n)) equals floor(sqrt n). The upper end m+1 satisfies (m+1)^2 >= floor(n)+1 > n. Tested: inclusion holds and the width is at most 2 grid steps. The input is required to be >= 0.
- `_div`: it rejects a denominator interval containing 0, encloses 1/b by the reciprocals of the endpoints with outward rounding, then multiplies. Tested.
- TAY (`taylor_coeffs`): f = (px, py, -x y^2, -x^2 y) is H = |p|^2/2 + x^2 y^2/2. The Cauchy-product recurrences x_{k+1} = px_k/(k+1) and px_{k+1} = -(x*yy)_k/(k+1) are correct, and the indices are in range (yy and xx have length p, which is enough for coefficients up to p). Inclusion isotonicity gives enclosure over a box.
- VAR (`variational_coeffs`): A matches Hess V = [[y^2, 2xy], [2xy, x^2]], with m added on the diagonal for the ghost. X_{k+1} = (1/(k+1)) sum_{i<=k} A_i X_{k-i}. A_p is set to 0 where yy[p] is missing, but A_p is never used for X up to index p, and in MVT the series is built with p+1, so A_p is present there.
- PIC (`lemma_rough`): for a continuous path u in Z, z0 + int_0^t f(u) lies in B + [0,h]·hull f(Z) componentwise. So the Picard operator maps C([0,h], Z) into itself, Schauder gives existence, and the Lipschitz property gives uniqueness. The non-strict inclusion is enough because Z is compact and convex. This also holds for h < 0 as a backward step, so the sign of h is harmless.
- PICV (`lemma_var_rough`): the same argument for the linear system X' = A(z(t))X with X(0) in J. The returned J + [0,δ]A(Z)Jt is a valid, tighter enclosure. Inflating the entries only changes the candidate, which is then verified.
- LAG (`lemma_taylor_point`): the componentwise Lagrange remainder z_i^{(p+1)}(ξ_i)/(p+1)! is the (p+1)-th Taylor coefficient at z(ξ_i), which lies in Z by PIC, because c is in B, and B is now the hull with c.
- MVT (`lemma_jacobian`): X(ξ+s) = X̃(s; z(ξ)) X(ξ), so the remainder is X_{p+1}[Z]·Xt in this order, and the code uses `mmul(XZ[p+1], Xt)`. Xt is PICV over [0,h] from I. Each entry of the mean-value integral of Dphi along a segment in the convex box B lies in the entrywise hull. Correct.
- NEU (`lemma_inverse`): ||A^{-1} - B|| = ||((BA)^{-1} - I)B|| <= ||B||·e/(1-e). Entrywise bounds follow from the inf-norm, and the interval R covers every point A and B.
- `lemma_step` (base part): for z = c + Aρ, phi(z) = phi(c) + J̃Aρ with J̃ in J, so A'^{-1}(phi(z) - c') lies in w, and w ⊂ r' gives phi(z) in c' + A'r'. Using the enclosure of A'^{-1} twice independently is an over-approximation, so it is sound. (The ghost frame part is F1.)
- KRA (`lemma_krawczyk`) and the construction of F and DF (`cert_check.py:88-124`): F = (x(T) - x0, px(T) - px0, y(T)) with z0 = (x0, 0, px0, sqrt(2E - px0^2)) is on the energy surface. Fc comes from chain P at exactly T_center (checked both ways). DF is Jfinal·dz0/d(x0,px0) - I in the (x, px) rows, with column T = f(Zf) and rows reordered to (x, px, y). It covers all of X because chain B runs from the box to T_lo and then PIC/PICV run over [T_lo, T_hi], and dpy/dpx0 = -px0/py0 with py0 bounded away from 0 by `_div`. x̂ is in X (checked), and the strict inclusion is the standard Krawczyk uniqueness theorem.
- `lemma_rough_time_interval`: z(t) is in Zf by PIC from the box at T_lo. X(t) = X̃(t - T_lo)·X(T_lo) satisfies the linear ODE with initial set Jacc, so PICV applies. δ = T_hi - T_lo >= 0 because the box is validated.
- T (`lemma_trace_stability`): for an autonomous Hamiltonian, the monodromy has a double eigenvalue 1 (flow direction and energy), and the remaining pair satisfies λ + 1/λ = tr - 2. The interval tests are correct.
- SYMP (`lemma_symplectic_4x4`): the characteristic polynomial (λ^2 - s1λ + 1)(λ^2 - s2λ + 1) gives a = s1 + s2 and b = e2 = s1s2 + 2, so s^2 - as + (b - 2) = 0. A complex s or |s| > 2 forces a multiplier off the circle. The enclosures of the roots are sound.
- FRAME (`run_chain`) for m = 0: A_0 = I is checked exactly, and R telescopes because the same A_k is used on both sides. The product of enclosures contains the chain-rule product. (F1 and F2 cover the exceptions.)

- Driver glue: no `assert` remains in either file, the only `except` is the top-level CheckFailed handler, and floats appear only in printing. The README's HI interval [-0.96097, -0.96088] is an outward rounding of the printed endpoints (see F3 for why that is luck, not design).

## Summary

The base-flow proof (existence and uniqueness by Krawczyk, ellipticity by Lemma T) checks out lemma by lemma. I found one real soundness hole: F1, the per-step choice of ghost frames, which affects the `--ghost` verdict for an adversarial certificate. It does not affect cert16, and I demonstrated it on a 4-step slice. F2 is a small gap that does not affect cert16 either. F3-F7 are cosmetic.

## Not reached, or out of scope

- The `dirstate` / `cert_twist` path inside `run_chain` and the `--twist` glue in `cert_check.py:104-109, 138-159` belong to Prompt B and were not reviewed.
- I did not check that the orbit in the box is the Dahlqvist–Russberg orbit. The checker proves a unique periodic orbit in X at the certificate's E; comparing it with the published initial data is a one-line float comparison the author can add to the README.
- Full runs were not attempted. After fixing F1 and F2, the author should rerun:
  `python3 code/cert_check.py results/cert16.json --ghost > results/cert16_check_output_postreview.txt 2>&1` (~35 min)
- Suggested regression test for F1 (fast): run `run_chain` on chainB steps 0-3 with `"G"` deleted from one step and assert that CheckFailed is raised. See scratch `rv/t1.py` for the construction.

Session ended at about 01:30 EDT, within budget.
