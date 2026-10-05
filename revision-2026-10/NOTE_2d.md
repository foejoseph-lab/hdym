# The two-field model: derivation, CPU checks, and what the GPU run tests

*2026-10-01. Companion to `hdym_2d.py`, `hdym2d_nu0.py`, `hdym2d_spectrum.py`. Section 5 is a
result about the 3-d paper and should go into `main.tex` before arXiv.*

## 1. The model and why it needs a regulator

The z = 0 plane of SU(2) Yang–Mills mechanics, A^a_i = diag(x, y, 0), is invariant
(ż = ṗ_z = 0 is preserved exactly) and on it V = ½x²y². So the "two-field x²y² model" is not
a different toy; it is the same gauge theory on a sub-ansatz. The ghost of the 3-d model
restricted to this plane has Hessian

    Hess V = [[ y²,  2xy,   0   ],
              [ 2xy,  x²,   0   ],
              [ 0,    0,  x²+y² ]]

so (w_x, w_y) is the 2-d ghost and w_z decouples exactly, driven by K_zz = x² + y².

**Pure x²y² in 2-d has infinite density of states.** Z₂ = ∫δ(1−H) d²q d²p = 2π·area{x²y² < 2},
and the area is log-divergent because the channel cross-section falls as 1/|x| (in 3-d it
falls as 1/x², which is why Z = 557.1 is finite there). Dynamically: a channel visit with
transverse action J lasts 2√2/J, the entry measure is uniform in J, so the mean visit time
∫dJ/J diverges. Time spent per octave of depth is constant, a trajectory sinks deeper the
longer it is followed, and the fraction of time in any bounded region — including the
resonant depth D = E^{-1/4} — decays as 1/ln T. Both λ_YM and λ⊥ then drift to zero as
1/ln T: slowly enough to give a number at any finite T, fast enough (~25% between T = 10⁴
and 10⁵) to make that number meaningless and to defeat an a + b/T convergence test. This is
Simon's (1983) classical/quantum distinction for x²y², already cited in the paper, and the
reason the 2-d literature mostly uses the Pullen–Edmonds form.

**Regulator.** Close the channels with a mass:

    V = ½x²y² + ½μ²(x² + y²),      μ² = μ₁² E^{1/2}.

With μ² scaled this way, V_E(E^{1/4} q) = E·V_1(q) exactly, so every energy is the rescaled
E = 1 system with regulator μ₁, the Hessian scaling Hess V(t; E) = E^{1/2} h(E^{1/4} t) holds,
and the spectral argument goes through unchanged. The channel is cut off at depth
D_cut = √2/μ₁ (14.1 at μ₁ = 0.1, 47 at μ₁ = 0.03). The ghost mass becomes 1 + μ², a 1%
effect at most (included below). The regulator is a parameter of the toy, not of the
physics; the test of the derivation is that **λ depends on μ₁ only through Z₂(μ₁)**, and
running two values of μ₁ checks the 1/Z factor, which the 3-d model cannot do.

## 2. Entry flux: ν₀ = 4π/Z₂

Same argument as §5.3 of the paper with one transverse action instead of two. On the surface
|x| = A, with transverse action–angle variables (J, φ) and ∫dp_x p_x δ(E_x − p_x²/2) = 1 for
the along-channel flux, the microcanonical entry-rate density is dJ dφ/Z₂ per end, so
2π dJ/Z₂ per end and, with two ends,

    ν₀ = 4π/Z₂,      Z₂ = 2π ∫Θ(1 − V_1) d²q = 8π ∫₀^{√2/μ₁} √((2 − μ₁²x²)/(x² + μ₁²)) dx.

| μ₁ | D_cut | Z₂ | ν₀ |
|---|---|---|---|
| 0.1 | 14.1 | 189.73 | 0.06623 |
| 0.03 | 47.1 | 275.32 | 0.04564 |

Z₂ grows as 16√2π ln(1/μ₁) = 71.1 per e-fold.

The along-channel motion at E = 1: ½ẋ² + J|x| + ½μ²x² = 1 (transverse frequency
Ω = √(x²+μ²) ≈ |x|, adiabatic invariant J). Turning point: J·D_max + ½μ²D_max² = 1.

**Checked** (`hdym2d_nu0.py`, 256 trajectories, E = 1, 59 000–62 000 visits per μ₁):

- rate(any channel, D_max > D) / [2 J_max(D)] with J_max = (1 − ½μ²D²)/D gives ν₀ to
  0.3% at every D from 4 to 12 (μ₁ = 0.1) and 4 to 24 (μ₁ = 0.03); 1.5–3% high at D = 32–40
  for μ₁ = 0.03, where the regulator term is 20–40% of the energy and the adiabatic |x|
  approximation to Ω starts to show;
- rate(J_entry < j)/2j = ν₀ to 0.3% for j = 0.02–0.3: the action measure is uniform;
- D_max·J + ½μ²D_max² = 1.000 (median), IQR ±0.025.

So the flux formula holds in 2-d with the regulator and the Z-dependence (0.0662 vs 0.0456)
is reproduced.

## 3. Spectrum of the diagonal Hessian entry

### 3.1 The resonant piece (drives the ghost)

During an x-channel visit at depth D, y² = (J/Ω)(1 + cos 2φ): a chirp of amplitude J/D at
frequency 2D. With the chirp spectrum π⟨(a²/2) δ(u − Ω)⟩ and the passage time
2dD/√(2(1 − JD − ½μ²D²)),

    ŝ_xx^res(u) = (2π/u²) ⟨J² δ(u − 2|x|)⟩ = (2π/u²) ν₀ ∫₀^{J_max} J² dJ /√(2(1 − JD − ½μ²D²)),  D = u/2.

Substituting J = c·a/D with c = 1 − ½μ²D²:

    ŝ_xx^res(u) = 16π I₂ ν₀ c^{5/2} / u⁵,      I₂ = ∫₀¹ a² da/√(2(1−a)) = B(3,½)/√2 = 16/(15√2),

    A_res ≡ u⁵ ŝ_xx^res / c^{5/2} = (256π/(15√2)) ν₀ = 37.93 ν₀    (2.511 at μ₁ = 0.1, 1.730 at 0.03).

[3-d check of the same algebra: I = ∫_{a+b<1}(a²+b²)/√(2(1−a−b)) = (2/3√2)B(4,½) = 64/(105√2),
the paper's value.]

### 3.2 The detuned piece (does not drive the ghost) — and the measured 5/4

`hdym2d_spectrum.py` (128 trajectories, Welch/Hann, segments of 164 time units) measures
u⁵ŝ_xx/(A_res c^{5/2}) = **1.25 ± 0.02, flat over u = 8–28 at μ₁ = 0.1 and u = 8–48 at
μ₁ = 0.03**. The u⁻⁵ shape and the c^{5/2} regulator factor are exactly right; the constant is
5/4 of A_res. The extra quarter is the y-channel:

In the y-channel, y is the along-channel coordinate and x the transverse oscillator,
x² = (J/|y|)(1 + cos 2φ). The equation of motion ÿ = −x²y − μ²y therefore contains the
oscillating force −J sgn(y) cos 2φ at frequency 2|y|; the slow coordinate responds with
δy = J sgn(y) cos 2φ/(4y²), and

    K_xx = y² = y_slow² + 2y δy = y_slow² + (J/2D) cos 2φ.

Amplitude J/(2D) at frequency 2|y|, against the resonant chirp's J/D at frequency 2|x|: one
quarter of the power, with the same J-statistics and the same passage integral. Total
u⁵ŝ_xx = (5/4) A_res c^{5/2}. The y-channel piece is detuned from the parametric resonance
of w_x (natural frequency √(1 + y²) ≈ |y| + 1/2|y|) by 1/D, against a resonance half-width
J/(4D²) ≤ 1/(4D³); it drives nothing for D ≳ 2, exactly the paper's argument for the 3-d
detuned power. Hence in 2-d the **resonant share of the diagonal-entry spectrum is 4/5**, and
the ratio of the measured λ_diag to (1/8)E^{3/4}ŝ_meas(2E^{-1/4}) applied to the *full*
measured spectrum should be 0.80 (the 3-d analogue of the 0.28 test in §5.2 of the paper).

### 3.3 The z-component's driver

K_zz = x² + y² + μ². In the x-channel: the y² chirp (J/D) plus the second-order response
2xδx = (J/2D) cos 2φ of the along-channel coordinate, coherent (same φ): amplitude 3J/(2D),
power 9/4 A_res; the y-channel gives the same. Prediction u⁵ŝ_zz = (9/2) A_res c^{5/2},
**all of it detuned** (w_z has frequency √(1 + x² + y²) ≈ D in either channel). Measured:
4.3–4.4 over u = 8–28 (2–4% low; second-order corrections to the response are O(1/D²)).

## 4. The ghost law in 2-d

With λ = S_K(2ω₀)/(8ω₀²), ω₀² = 1 + μ², and the resonant piece only,

    λ_diag(E) = (π/(15√2)) ν₀ E² · c^{5/2} (1+μ²)^{-7/2}
              = (4π²/(15√2)) E²/Z₂(μ₁) · (1 − μ₁²/(2√E))^{5/2} (1 + μ₁²√E)^{-7/2}
              = 1.861 E²/Z₂(μ₁) · [regulator factors].

Exponent 2 = ¾ (kinematic) + ¼ (one transverse action: P(J < j) ∝ j, exceedance D⁻¹) +
1 (driving power ∝ (J/D)²). The 3-d decomposition 9/4 = ¾ + 2/4 + 1 differs only in the
middle term, which is the count of transverse actions. The regulator factor matters:
0.88 at E = 0.01, μ₁ = 0.1. Validity: E ≫ μ₁⁴/4 (= 2.5×10⁻⁵ at μ₁ = 0.1).

| μ₁ | E | u = 2E^{-1/4} | regulator factor | predicted λ_diag |
|---|---|---|---|---|
| 0.1 | 0.1 | 3.56 | 0.950 | 9.32×10⁻⁵ |
| 0.1 | 0.03 | 4.81 | 0.923 | 8.15×10⁻⁶ |
| 0.1 | 0.01 | 6.32 | 0.877 | 8.60×10⁻⁷ |
| 0.03 | 0.03 | 4.81 | 0.994 | 6.04×10⁻⁶ |
| 0.03 | 0.01 | 6.32 | 0.989 | 6.68×10⁻⁷ |

**CPU results** (`hdym_2d.py --backend c`, N = 8192, μ₁ = 0.1, two cores; `results_2d_cpu/`):

| E | T | λ_diag (measured) | a + b/T extrap. | predicted | ratio | λ_z | λ_z/λ_diag |
|---|---|---|---|---|---|---|---|
| 0.1 | 12 000 | 8.19×10⁻⁵ ± 1.2% | 7.93×10⁻⁵ | 9.32×10⁻⁵ | 0.87 | 1.84×10⁻⁵ | 0.22 |
| 0.03 | 60 000 | 7.96×10⁻⁶ ± 1.8% | 8.10×10⁻⁶ | 8.15×10⁻⁶ | 0.98 (0.99 extrap.) | 1.5 ± 1.3×10⁻⁷ | 0.02 |

**GPU results** (RTX 4080 SUPER, 2026-10-01, N = 10⁵ at E = 0.1 and 2×10⁵ otherwise,
`--target-rel 0.02`; 0.03–0.05 ms/step, 46 min in all):

| E | μ₁ | u | T | λ_diag | ± | extrap. | predicted | ratio (extrap.) | λ_z | λ_z/λ_diag |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.1 | 3.56 | 6 000 | 8.077×10⁻⁵ | 0.5% | 8.21×10⁻⁵ | 9.322×10⁻⁵ | 0.866 (0.881) | 1.90×10⁻⁵ | 0.235 |
| 0.03 | 0.1 | 4.81 | 60 000 | 7.759×10⁻⁶ | 0.4% | 7.91×10⁻⁶ | 8.155×10⁻⁶ | 0.951 (0.970) | 1.33×10⁻⁷ | 0.017 |
| 0.01 | 0.1 | 6.32 | 340 000 | 8.847×10⁻⁷ | 0.4% | 8.79×10⁻⁷ | 8.598×10⁻⁷ | 1.029 (1.022) | 4.0 ± 2.6×10⁻⁹ | 0.005 |
| 0.03 | 0.03 | 4.81 | 100 000 | 5.747×10⁻⁶ | 0.4% | 5.75×10⁻⁶ | 6.041×10⁻⁶ | 0.951 (0.952) | 8.6×10⁻⁸ | 0.015 |

- **Exponent:** local slope 0.03 → 0.01 is 1.98 (predicted 2.05 with the regulator factors);
  the 3-d value 2.25 is excluded at ~10σ. With n = 2 in the paper and n = 1 here, the
  decomposition ¾ + n/4 + 1 is tested at two values of n.
- **1/Z:** μ₁ = 0.1 → 0.03 changes Z₂ by 45% and the regulator correction by 7%; the ratio to
  prediction is 0.951 in both cases. The 5% shortfall at u = 4.8 is independent of Z and of
  the regulator; it is the same approach from below that the 3-d data show at the same u
  (0.971 at E = 0.03, 0.995 at u = 6.7). The law holds to 5% over 0.01 ≤ E ≤ 0.03 and to 3%
  at the lowest point.
- **Detuned driver:** λ_z/λ_diag = 0.235, 0.017, 0.005 at E = 0.1, 0.03, 0.01 — falling
  faster than any power of E, as detuned driving plus exponentially decaying interior power
  should. Direct confirmation of the paper's inertness claim.
- **Resonant share:** with the measured spectrum (`spec_2d_mu0.1.npz`) at the resonant u,
  λ_diag/[(1/8)E^{3/4}ŝ_meas(2ω₀E^{-1/4})/ω₀²] = 0.83 at both E = 0.03 and 0.01, against
  the derived 4/5 (§3.2). Second system in which the ghost sees exactly the resonant
  fraction of the tail (2/7 in 3-d).
- Note `ymerr = 2.3×10⁻⁵` on the μ₁ = 0.03 run (channels to depth 47 at the sub-stepping
  threshold); harmless here, use `--dt0 0.005` for μ₁ ≤ 0.01.

Earlier CPU runs (N = 8192, two cores) gave 0.87 and 0.98 at E = 0.1 and 0.03, consistent.

## 5. A result for the 3-d paper: the detuned share is exactly 5/7

The same second-order bookkeeping in 3-d derives the number §5.2 of the paper measures but
does not explain. For K_xx = y² + z² at E = 1:

- x-channel (resonant for w_x): (J_y/D) cos 2φ_y + (J_z/D) cos 2φ_z; averaged over the
  relative phase, power ∝ (J_y² + J_z²)/D². [The paper's formula.]
- y-channel (detuned): the direct chirp z² → (J_z/D) cos 2φ_z, plus the second-order
  response of y to both transverse oscillators, 2yδy = (J_x cos 2φ_x + J_z cos 2φ_z)/(2D),
  coherent in φ_z: K_xx^osc = (3J_z/2D) cos 2φ_z + (J_x/2D) cos 2φ_x, power
  (9J_z² + J_x²)/(4D²).
- z-channel: likewise (9J_y² + J_x²)/(4D²).

The uniform measure on the two transverse actions makes all ⟨J_a²⟩ equal and the three
channels equally visited, so resonant : detuned : total = 2 : 5 : 7.

    A/A_res = 7/2  →  A = 21.5     (measured 21.6 ± 0.3, Eq. (A) of the paper)
    resonant share = 2/7 = 0.2857   (measured 0.283, 0.285, 0.290, 0.276 at E = 0.008–0.03)

Suggested change to §5.2: replace "the rest being the detuned y- and z-channel chirps and
the interior" by the derivation above and "the detuned share is 5/2 of the resonant one,
so A = (7/2)A_res = 21.5 against the measured 21.6 ± 0.3, and the ghost sees 2/7 of the
tail power, against the measured 0.276–0.290". The remark about the interior can go: the
tail is entirely channel power. For K_zz-type (the lowest-eigenvalue) spectra the same
method gives the constant too, if wanted.

## 6. What the GPU run tests, run by run

1. **E = 0.1, 0.03, 0.01 at μ₁ = 0.1.** Exponent: local slopes should be ≈ 2 (1.97 allowing
   the regulator factor's E-dependence; the 3-d value 2.25 is excluded at the 10σ level by
   the E = 0.01/0.03 pair alone). Prefactor: ratio to the table ≈ 1 at E ≤ 0.03 to the
   same ~5% as 3-d. This is the test of the exponent decomposition.
2. **E = 0.03 (and 0.01 if time) at μ₁ = 0.03.** Z₂ changes by 45%; the ratio should stay 1.
   This is the test of the 1/Z factor.
3. **λ_z/λ_diag at each E.** Should fall with E (0.22 at E = 0.1 on CPU). It is the
   direct test of "detuned power is inert"; a value that does not fall would mean the
   inertness argument is incomplete and the 2/7 of §5 is a coincidence.
4. Resonant-share test: λ_diag / [(1/8)E^{3/4}ŝ_meas(2E^{-1/4})] = 0.80 using the measured
   spectrum (`spec_2d_mu0.1.npz`); the 2-d analogue of the 0.28 test.
5. (Optional, cheap) `--mu1 0` at E = 0.03, long run, no convergence target: watch λ_diag(T)
   fall as 1/ln T. One figure's worth of footnote.

## 7. Files

- `hdym_2d.py` — kernel (CuPy / C / numpy), per-component estimator, prediction printed with
  every run, CSV with `ratio` and `lam_z`.
- `hdym2d_nu0.py` — flux check (§2). `hdym2d_spectrum.py` — spectrum check (§3).
- `nu0_2d_mu0.1.npz`, `nu0_2d_mu0.03.npz`, `spec_2d_mu0.1.npz`, `spec_2d_mu0.03.npz` — data.

## 8. The general-n statement (`hdym_analytic.py`, sympy, all checks pass)

For a diagonal potential V = ½Σ_{i<j} q_i²q_j² with N = n + 1 degrees of freedom (n
transverse actions per channel), the two inputs — microcanonical entry flux
ν₀ = 2(2π)ⁿ/Z_N and adiabatic channel passage — give, exactly:

    I_n   = ∫_{Σa<1} (Σa_i²) /√(2(1−Σa)) dⁿa = √(2π) n / Γ(n + 5/2)
    A_res = u^{n+4} ŝ_res(u) = 2π·2^{n+2} I_n ν₀
    λ_diag = (π/16) I_n ν₀ E^{(n+7)/4} = C_n E^{(n+7)/4} / Z_N,   C_n = 2^{n−5/2} π^{n+3/2} n / Γ(n + 5/2)

    n = 1: C₁ = 2√2π²/15 = 1.861,  exponent 2;      n = 2: C₂ = 16√2π³/105 = 6.682,  exponent 9/4.

Second-order response of the along-channel coordinate: in channel c the oscillator a adds
2YδY = (J_a/2D) cos 2φ_a to Y²; so in K_xx (channel c ≠ x) the coefficient of cos 2φ_b is
3/2 for b ∉ {x, c} and 1/2 for b = x, and

    total tail / resonant tail = (9n − 4)/4,     resonant share = 4/(9n − 4):   4/5 (n=1),  2/7 (n=2).

Regulated 2-d density of states in closed form, b² = μ⁴/2, m = 1/(1+b²):

    Z₂(μ) = 8√2π √(1+b²) [K(m) − E(m)]  →  16√2π ln(1/μ) + 8√2π(ln 4√2 − 1) + O(μ⁴ ln μ)
          = 189.733 (μ₁ = 0.1), 275.316 (0.03); the asymptotic form is already exact to 4 digits.

The only non-closed-form number left in the 3-d law is Z₃ = 557.1, a single quadrature.
Everything else is exact *given ergodicity of the Yang–Mills flow on the energy shell*,
with corrections O(E^{1/2}) from adiabaticity and from the weak-driving expansion. That is
how the paper should state it.
