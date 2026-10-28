# ym/ — Lee–Wick Yang–Mills mechanics: code and data for the paper

Everything behind `paper/main.tex`. Rewritten 23 Oct 2026 from the state at that date. The system is the spatially
homogeneous, colour-diagonal sector of SU(2) Lee–Wick Yang–Mills on a three-torus: 12-dimensional Hamiltonian, units
g = M = 1, state (p, w, P, W) with q = p − w, V = ½(x²y² + y²z² + z²x²),
H = ½|P|² − ½|W|² + V(q) + w·∇V(q) − ½|w|². The ghost-free manifold w = W = 0 is invariant and carries pure
Yang–Mills mechanics; the paper measures and derives its transverse Lyapunov exponent λ⊥(E), extends the derivation to
the two-field model, and linearises the full homogeneous theory about the diagonal background (the off-diagonal sector).

All scripts are in this one folder so that imports between them resolve without setting `PYTHONPATH`. Committed
results are in `data/`; notes in `notes/`; the shared C model and the drift census in `csrc/` and `tools/`; the
first-draft scripts in `early/`. The off-diagonal sector on a periodic orbit, and the nonlinear escape kernel built on
`hdym_full_derive.py`, live in `islands/` and `rabbit_holes/ym/`.

## Requirements

Python ≥ 3.10, numpy, scipy, matplotlib; sympy for the derivation scripts (`hdym_derive.py`, `hdym_analytic.py`,
`hdym_rung1.py`, `hdym_full_derive.py`). The production runs use CuPy with an NVIDIA GPU (`--backend cupy`); the same
C kernels compile with gcc for a multicore CPU (`--backend c`, built in a temporary directory), and a pure-numpy
reference backend exists for checking. All compiler calls except `tools/arb_smoke.py`'s (which repeats the PATH fix
itself) go through `cbuild.py`, which puts the compiler's own bin
folder first on the compiler's PATH and prints the compiler's output and exit code on failure (added 12 Oct after
gcc's `cc1.exe` failed silently when started from a Spyder 6 console whose runtime folder ships older `zlib1.dll`,
`libzstd.dll`, `libwinpthread-1.dll`, `libgcc_s_seh-1.dll` ahead of `C:\msys64\mingw64\bin`). Optional, for
`tools/arb_smoke.py` only: `python-flint` and FLINT's C library.

## Scripts

Section, figure and table names refer to `paper/main.tex`.

| file | what it does | paper |
|---|---|---|
| `hdym_core.py` | equations of motion, Yoshida-4 symplectic integrator with per-trajectory channel sub-stepping, CPU Lyapunov estimators, adaptive DOP853 confirmation of blow-up | The model; Numerical methods |
| `hdym_derive.py` | symbolic (sympy) reduction from the 3+1 theory to the 12-d system: Gauss law, off-diagonal equations, symmetric criticality, auxiliary-field Hamiltonian | ref. note:derive |
| `hdym_full_derive.py` | the FULL nonlinear homogeneous system, all nine A^a_i, temporal gauge: fourth-order equations solved for A'''', the three Gauss constraints, the Ostrogradsky energy; checks against `rung1_system.pkl` (linearisation, 10⁻¹⁵; the script's tolerance is 10⁻⁹, `hdym_full_derive.py:248`; no committed log in `ym/`) and the diagonal closed form; writes `full_system.pkl` and the generated `csrc/hdym_full_rhs.h` for `csrc/hdym_full_escape.c` (19 Oct) | — (nonlinear escape runs, `islands/`) |
| `hdym_model.py` | the Python side of `csrc/hdym_model.h`: pastes the header into C and CUDA kernel sources, compiles `csrc/hdym_model.c` for ctypes (`csrc/README.md`) | — |
| `hdym_gpu.py` | transverse Lyapunov exponent λ⊥(E), whole integrator in one kernel; `--ghost vec` (physical ghost) or `--ghost diag` (three independent scalar oscillators, per-component estimator) | fig:lambda, tab:low |
| `hdym_gpu_landscape.py` | `hdym_gpu.py` plus the (D, J) histogram output, results bit-identical | (landscape figure, not in the paper) |
| `hdym_spectrum.py` | multitaper spectrum of Hess V along pure-YM trajectories at E = 1; the u⁻⁶ tail | Spectral formula and the u⁻⁶ tail |
| `hdym_tailA.py` | reads a saved `hdym_spectrum.py` run (`spectrum.npz`) and prints the band slopes and the tail constant A = ⟨u⁶ s(u)⟩ over u = 8–45 (12 Oct; `hdym_spectrum.py` itself does not print A) | eq:A |
| `hdym_channel_check.py` | channel-passage checks: adiabatic invariant J, parabolic x(t), chirp at 2\|x\|, P(J) ∝ J | Channel statistics |
| `hdym_nu0.py` | channel entry rate ν₀ at E = 1, measured (the predicted 8π²/Z comes from `hdym_Z.py`) | ν₀ from the density of states, fig:nu0 |
| `hdym_analytic.py` | sympy derivations, every constant checked: general-n law, 4/(9n−4) share, Z₂ elliptic form, c_j = 0 law 0.0210 E^{5/4} | Derivation; Off-diagonal sector |
| `hdym_Z.py` | the density of states Z = ∫δ(1−H) d³q d³p at E = 1, two independent quadratures (radial/angular and Cartesian), then ν₀ = 8π²/Z and the λ_diag and c_j = 0 coefficients; seconds; assumes no value of Z (12 Oct) | ν₀ from the density of states; eq:final |
| `hdym_2d.py`, `hdym2d_nu0.py`, `hdym2d_spectrum.py` | two-field model: kernel, entry flux, spectrum (`notes/NOTE_2d.md`) | Test on the two-field model, tab:twofield |
| `hdym_rung1.py` → `rung1_system.pkl` | linearisation of the full homogeneous theory about the diagonal background (`notes/NOTE_rung1.md`); also the input of the Floquet work in `islands/` | Off-diagonal sector |
| `hdym_rung1_lyap.py` | CPU Lyapunov spectra of the off-diagonal blocks (reads `rung1_system.pkl` from the working directory) | Off-diagonal sector |
| `hdym_rung1_gpu.py` | kernel for the off-diagonal ghost exponents: `--cj`, growth-location split, (D, J) landscape | Off-diagonal sector, tab:offdiag |
| `hdym_spectrum_xy.py` | spectrum of the xy Hessian entry at E = 1 | Off-diagonal sector |
| `hdym_escape_map.py` | escape-time landscape on 2-d sections of initial conditions | fig:maps (a, b) |
| `hdym_uncertainty.py` | uncertainty exponent α of the escape landscape | fig:maps (c) |
| `hdym_precision.py`, `hdym_diffpair.py`, `hdym_ball.py` | the roundoff study (below): mantissa-precision sweep; pair integration in "difference numbers"; certified pair escape times in ball arithmetic (rigour status in the docstring of `hdym_ball.py`); `hdym_precision.py` is also imported by `islands/code/hdym_island.py` | The escape-time landscape |
| `make_fig1.py`, `make_figs.py`, `make_fig_Q.py` | regenerate fig:lambda, fig:maps, fig:nu0 and fig:Q from `data/` into `paper/` | |
| `make_fig_landscape.py` | (D, J) landscape figure; its inputs `results_gpu/landscape_E0.1_diag.npz`, `results_rung1/landscape_E0.1_cj{0,1}.npz` are not committed | not in the paper |

`paper/fig_floquet.pdf` (fig:floquet): its data come from `early/floquet.py` (committed output `data/early/floquet.npz`);
the plotting script is not in the repository.

### Known README-vs-code mismatches (found 12 Oct, held)

Found by a README-against-code pass on 12 Oct (repository at 8b866214) and held, not started, by decision of that day;
when one is taken up, fix it in its home and strike it here. (The bianchi9 list of the same pass is in `bianchi9/README.md`,
verified against its scripts on 23 Oct; these `ym/` items were not re-verified.) Code that a paper claim or note depends on:
`hdym_uncertainty.py` default thresholds 2/5/10 vs the committed 5/10/20 (a README rerun would plot threshold 5 as
"> 10"; feeds α); `hdym_rung1_lyap.py` `load_K` leaves c_j free and cannot run (`NOTE_rung1` CPU spectra; the paper's
table is from `hdym_rung1_gpu.py`, own K); `make_fig1.py` / `make_fig_Q.py` hard-code points and run lengths although
the README says figures regenerate from `data/`. Wording: the `hdym_spectrum.py` row claimed A and the `hdym_nu0.py`
row claimed 8π²/Z (both rows reworded on 23 Oct to say what the scripts compute; struck); the numpy reference backend renormalises by the vector norm in both ghost
modes; `islands/code/dr_orbit_ghost.py`, `scan.py`, `interior.py` label digits ~10⁻¹² outside the certified box as
the DR orbit; stale checker line counts (now in `islands/README.md`). Also: `hdym2d_spectrum.py`'s docstring and
printout label A_res/ν₀ as 37.93 where the closed form gives 37.913 (`notes/NOTE_2d.md`).

### `csrc/` and `tools/` — the shared model and the drift census (13 Oct)

Every kernel script above carries its own copy of V, ∇V, Hess V·w and the Yoshida-4 coefficients, in Python, in C
and in CUDA (as strings inside the `.py` files): fifteen copies. `csrc/hdym_model.h` is the one copy they are meant
to converge on, written in the same operation order as the kernels so that a kernel built on it reproduces the
committed data bit for bit; `csrc/README.md` says what is in it and how a kernel uses it. No script has been moved
onto it yet (phase 2, pending): the scripts and their data are as they were. `tools/drift_census.py` evaluates every
existing copy at the same random points and compares it with the header (bitwise / roundoff / DIFFERENT; run from
`ym/`: `python tools/drift_census.py --out data/drift_census_<date>.txt`). Result (`data/drift_census_2026-10-13_cloud.txt`,
cloud CPU): 83 quantities bitwise, 22 at roundoff (≤ 7 ulp), none different, CUDA skipped; on an RTX 4080 SUPER
(`data/drift_census_2026-10-13_win.txt`) the same plus the CUDA copies at roundoff (≤ 8 ulp of their terms, because
NVRTC contracts `a*b + c` into fused multiply-adds and gcc does not), bitwise again with `--fmad=false`. The CPU
roundoff rows are the copies that compute the Yoshida coefficients from the formula instead of the literals
(`hdym_core.py`, `early/`; 1, 3, 1, 0 ulp for YC0/YC1/YD0/YD1, census rows "Yoshida C" 3 ulp, "Yoshida D" 1 ulp; stated
as "one ulp in three of the four" until 28 Oct 2026), the two operation orders of Hess V·w (linearised vs
nonlinear kernels), `hdym2d_nu0.grad`, which writes `x*y*y + mu2*x` where the kernel has `x*(y2 + mu2)`,
`hdym_diffpair` "dV vs V(q+dq)−V(q)" (7 ulp, the maximum) and `hdym_ball.py [C]` `kick_b` (4 ulp, double-word centre).
`tools/arb_smoke.py` (14 Oct) is the smoke test for Arb as the centre arithmetic of `hdym_ball.py`
(`notes/NOTE_precision_tier.md`). `csrc/` also holds the nonlinear escape kernel `hdym_full_escape.c` and its
generated right-hand side, documented in `csrc/README.md` and `islands/README.md`.

### `early/` — the first-draft scripts (restored)

Plain-numpy scripts from before the GPU code, removed from the repository on 1 Oct (as `other_code.7z`) and
restored unchanged except for two import paths. Several numbers in the paper still rest on them. Each writes its
output to the working directory.

| file | what it does | output |
|---|---|---|
| `transverse.py` | linearised ghost on pure-YM trajectories: λ⊥ for the vector ghost and the diagonal ghost, tachyonic fraction; `python transverse.py E T N`. Its diagonal ghost is renormalised with the **vector norm**, not per component | `transverse_E*.npz` (committed in `data/`); at E = 1 also `spectrum_E1.npz` |
| `lyap_full.py` | Benettin maximal Lyapunov exponent of the full 12-d system at E = 0.3, ε = 0.01, with escape times | `lyap_E0.3_eps*.npz` |
| `escape.py` | escape times at other ε (e.g. 10⁻⁴) and the blow-up mechanism; needs `lyap_full.py` beside it | `escape_eps*.npz` |
| `floquet.py` | Floquet growth of the ghost on the diagonal periodic orbit x = y, z = 0; unstable tongues | `floquet.npz` (committed in `data/early/`) |
| `spectrum.py` | spectrum of K_xx = y² + z² at E = 1, conditioned on channel occupancy | `spectrum_E1_full.npz` |
| `spec.py` | spectrum of Hess V along pure-YM trajectories at E = 1 (uses `../hdym_core.py`) | `spec.npz` |
| `channels.py` | channel-visit statistics at E = 1, independent recount for fig:nu0 | `nu0_E1_rederived.npz` |
| `integrals.py` | search for polynomial first integrals of the 12-d system, with a separable control | printed only |

## Notes (`notes/`)

| file | what it is |
|---|---|
| `NOTE_2d.md` | the two-field model: derivation, checks, and the test of the ghost law at n = 1 (1–12 Oct) |
| `NOTE_rung1.md` | the off-diagonal linearisation, the neutrality result, the off-diagonal ghost exponents (1–12 Oct) |
| `NOTE_precision_tier.md` | the precision scheme behind `hdym_ball.py`, why its lost pairs are arithmetic, and the Arb tier that would recover them (13–14 Oct; nothing built) |
| `REVIEW_hdym_ball_2026-10-11.md` (+ `_sweep_small.txt`) | record of the independent review of `hdym_ball.py`; all findings acted on (below) |

## Data (`data/`)

| file / folder | contents |
|---|---|
| `transverse_vec_percomp.csv` | physical ghost, E = 0.008, 0.015, 0.02, 0.03, N = 2×10⁵ |
| `transverse_diag_percomp.csv` | diagonal ghost, E = 0.008, 0.015, 0.02, N = 2×10⁵ (per-component estimator) |
| `transverse_diag_percomp_hi.csv` | diagonal ghost, E = 0.05 … 1, N = 10⁵ (per-component estimator) |
| `results_percomp.7z` | the run folders behind the three CSVs (`lam_i_E*.npy`, `history_E*.npy`, `transverse_gpu.csv`), including the diagonal E = 0.03 run (N = 10⁵), which has no row in the CSV. Unpack to `data/results_percomp/` (ignored by git) for `make_fig_Q.py`: `7z x results_percomp.7z -oresults_percomp`, or with `py7zr` |
| `transverse_E*.npz` | independent numpy re-run at E = 0.215, 0.3, 0.43, 0.6, 1 (open symbols in fig:lambda), written by `early/transverse.py` (`lam_d` = vector-norm diagonal estimator) |
| `nu0_E1.npz` | channel visits (depth, entry action) at E = 1 |
| `Z_density_of_states.txt` | output of `hdym_Z.py`: Z = 557.088806 (the two methods agree to 10⁻¹³), ν₀ = 0.141731, λ_diag = 0.011994 E^{9/4}, c_j = 0 law 0.020990 E^{5/4} (12 Oct) |
| `analytic_output.txt` | output of `hdym_analytic.py` (12 Oct): every symbolic check `True`; 3-d law 0.01199 E^{9/4} (with Z = 557.1); predicted A = (7/2) A_res = 21.49; Z₂ = 189.7330 (μ = 0.1), 275.3160 (μ = 0.03); c_j = 0 law 0.0210 E^{5/4} |
| `results_spectrum/` | `spectrum.npz` and `prediction.csv` from `hdym_spectrum.py --N 2048 --T 20000 --t-burn 500 --workers 12` (12 Oct; log `spectrum_log.txt`) |
| `tailA_output.txt` | output of `hdym_tailA.py results_spectrum`: 98 304 segments, band slopes −6.14, −6.05, −6.00, −5.98, −5.89, A = 21.65 ± 0.13 (paper eq:A: 21.6 ± 0.3; predicted (7/2) A_res = 21.49) |
| `map256.npz`, `ghost256.npz`, `unc_ym.npz`, `unc_ghost.npz` | escape maps and uncertainty-exponent data at E = 0.3, ε = 0.01 |
| `nu0_2d_mu*.npz`, `spec_2d_mu*.npz` | two-field entry flux and spectrum at μ₁ = 0.1, 0.03 |
| `results_2d/` | two-field GPU runs: summaries `transverse_2d_mu*.csv` (`hdym_2d.py` appends to `<out>/transverse_2d.csv`; one file per μ₁ here), per-trajectory exponents `lam_i_E*_mu*.npy`, `history_*.npy` |
| `results_2d_cpu/` | two-field CPU runs, N = 8192, with their CSVs and histories |
| `results_rung1/` | off-diagonal GPU runs: summary `rung1_offdiag.csv`, per-trajectory exponents (c_j = 1 and 0), `history_*.npy` |
| `rung1_lyap_E*.npy`, `rung1_ym_E*.npy` | off-diagonal CPU Lyapunov spectra (ghost blocks, YM blocks) |
| `spec_xy_E1.npz` | xy-entry spectrum at E = 1 |
| `drift_census_2026-10-13_{cloud,win}.txt` | the drift census (above) |
| `results_precision/` | the roundoff study (below): the extended-precision difference-pair run at E = 0.3, ε = 0.01; `2026-10-03/` the precision study of that day (`results/`, `figures/`; its `ball_sweep.npz` predates the rigour fixes); `2026-10-10/` the certified ball sweep with the rigorous step remainder; `2026-10-12/` the production sweep with the post-review kernel (`ball_sweep.npz`, `.png`, log `ball_sweep.txt`, `selftest.txt`) |

### Certified pair escape times (`hdym_ball.py`, `results_precision/2026-10-12/`)

What the sweep certifies, after an independent review out of project context on 12 Oct: the paper's α and its "no
flattening down to 10⁻¹⁴" come from the float sweeps (`hdym_uncertainty.py`, observable max |component| > 10⁴). The
ball sweep (ghost escape |ghost| > 10, E = 0.3, ε = 0.01, M = 400, 20 δ × 3 thresholds) certifies a lower bound only:
certified-different pairs are a rigorous lower count, f_true ≥ f_lo for the sampled pairs, and a one-sided 95 %
Clopper–Pearson bound turns that into a statement about the section — at δ = 10⁻¹⁴ × width, 9 of 400 pairs certified
different (f_lo = 0.0225), section fraction ≥ 0.012 (|Δt| > 10; at 10⁻¹²: 36/400, ≥ 0.068). So the fine structure at 10⁻¹⁴ is not an
artefact of floating-point roundoff, certified for the integrator map (the exact Yoshida map with its pair-dependent
step schedule) and the ghost observable. That the map inherits it from the continuous flow is plausible (the
time-discretisation error is smooth in the initial data) but not certified. It does not certify α (slopes of f_lo,
f_hi do not bracket it), and from 10⁻¹⁶ down to the 10⁻⁴⁰ the sweep samples it says nothing either way: no pair is certified different, and the 9–18 %
uncertified pairs are exactly the lost (long-lived) ones. The `.png` draws f_lo = 0 clipped to 10⁻⁴ on the log axis.

History of the kernel: the step remainder made rigorous on 10 Oct; the two frame-arithmetic constants computed rather
than asserted the same day (a 100-pair, 11-δ comparison gave identical certified fractions); reviewed cold on 11 Oct
(`notes/REVIEW_hdym_ball_2026-10-11.md`: nothing breaks the main guarantee; fixed the same day — a lost member b could
still be recorded as certainly escaped, F1; new selftests 2b (member b's difference ball against a 60-digit orbit at δ = 10⁻³⁰) and 3 (escape bookkeeping
and both member-switch paths against 60-digit exit steps), F2; the f_lo / f_hi bracket is for the sampled pairs and fitted exponents are not a bracket on α, F3, wording and
printout; pair-dependent step schedule, D to 1 ulp, escape observable, F4–F6, wording); the production sweep rerun on
12 Oct with that kernel (Windows, gcc 15.2 — no committed log records the compiler; the Windows census of the same
night, `data/drift_census_2026-10-13_win.txt`, shows gcc 16.2.0): every certified fraction identical to `2026-10-10/`, and the selftest
numbers identical to a Linux run — the rigour fixes changed bounds, not classifications. The review's F7 (QR
containment dominates the lost pairs) was answered on 13 Oct: the losses are arithmetic (`notes/NOTE_precision_tier.md`).

The production sweep, its flags as recorded in the `args` of `2026-10-12/ball_sweep.npz` (working directory inferred
from the log's last line `wrote 2026-10-12\ball_sweep.png`):

```
python ../../hdym_ball.py --E 0.3 --eps 0.01 --section ym --M 400 --delta-max 1e-2 --delta-min 1e-40 --per-decade 0.5 --tmax 1500 --esc ghost --R 10 --backend c --out 2026-10-12 --tag ball_sweep
```

## Reproducing the main result

Commands write into `--out` relative to the working directory. The four low-energy diagonal-ghost points of tab:low
(the parameter-free comparison with 0.0120 E^{9/4}):

```
python ym/hdym_gpu.py --E 0.03 --N 100000 --ghost diag --backend cupy --t-burn-ym 200 --target-abs 3e-7 --t-block 5000 --tmax 300000 --out results_diag_percomp
python ym/hdym_gpu.py --E 0.02 0.015 0.008 --N 200000 --ghost diag --backend cupy --t-burn-ym 200 --target-abs 3e-8 --t-block 20000 --tmax 1500000 --out results_diag_percomp
```

The physical ghost at the same energies, and the diagonal ghost at higher energies:

```
python ym/hdym_gpu.py --E 0.03 0.02 0.015 0.008 --N 200000 --ghost vec --backend cupy --t-burn-ym 200 --target-abs 3e-8 --t-block 20000 --tmax 1500000 --out results_vec_percomp
python ym/hdym_gpu.py --E 0.05 0.07 0.1 0.15 0.2 0.3 0.5 1.0 --N 100000 --ghost diag --backend cupy --t-burn-ym 200 --target-abs 3e-7 --out results_diag_percomp_hi
```

On an RTX 4080 SUPER these take roughly 0.5 h, 4 h, 4 h and 10 min respectively. Each run writes
`<out>/transverse_gpu.csv` plus per-energy `history_E*.npy` and `lam_i_E*.npy`; the three CSVs in `data/` are those
files renamed.

Spectrum, channel statistics and ν₀ (CPU, pure Yang–Mills, no ghost):

```
python ym/hdym_spectrum.py --N 2048 --T 20000 --t-burn 500 --workers 12
python ym/hdym_channel_check.py --N 2048 --T 4000 --A-in 3
python ym/hdym_nu0.py
```

Escape landscape and uncertainty exponent at E = 0.3, ε = 0.01 (the committed `map256.npz`, `ghost256.npz` were made
with `--backend cupy`; `--backend c` reproduces them only to roundoff, `csrc/README.md`):

```
python ym/hdym_escape_map.py --E 0.3 --eps 0.01 --section ym    --n 256 --backend c
python ym/hdym_escape_map.py --E 0.3 --eps 0.01 --section ghost --n 256 --backend c
python ym/hdym_uncertainty.py --E 0.3 --eps 0.01 --section ym    --M 20000 --backend cupy --delta-min 1e-14 --delta-max 1e-4 --n-delta 11
python ym/hdym_uncertainty.py --E 0.3 --eps 0.01 --section ghost --M 20000 --backend cupy --delta-min 1e-14 --delta-max 1e-4 --n-delta 11
```

Two-field and off-diagonal runs: commands and GPU summary lines are in `notes/NOTE_2d.md` and `notes/NOTE_rung1.md`.

## Figures and paper

```
python ym/make_fig1.py && python ym/make_figs.py && python ym/make_fig_Q.py   # from any directory
cd paper && latexmk -pdf main.tex
```

## A note on the diagonal-ghost estimator

With the off-diagonal Hessian entries dropped the three ghost components are independent scalar oscillators.
Renormalising the norm of the three-vector returns the *largest* of three finite-time exponents and overestimates the
single-mode rate by a factor 1.5–1.7 at E ≤ 0.03 (1.1 at E = 1). `hdym_gpu.py --ghost diag` renormalises each
component separately and averages the three log-growths; the vector ghost (`--ghost vec`) uses the six-vector norm,
which is correct for the coupled system. See the footnote in the Numerical methods section of the paper.
