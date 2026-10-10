# hdym — the Lee–Wick ghost of Yang–Mills mechanics

Code and data for *The Lee–Wick ghost of Yang–Mills mechanics is malicious at every energy:
a parameter-free E^{9/4} instability law and an area-filling escape landscape* (J. Foe, 2026).

The system is the spatially homogeneous, colour-diagonal sector of SU(2) Lee–Wick Yang–Mills
on a three-torus: 12-dimensional Hamiltonian, units g = M = 1, state (p, w, P, W) with
q = p − w, V = ½(x²y² + y²z² + z²x²),
H = ½|P|² − ½|W|² + V(q) + w·∇V(q) − ½|w|².  The ghost-free manifold w = W = 0 is invariant and
carries pure Yang–Mills mechanics; the paper measures and derives its transverse Lyapunov
exponent λ⊥(E).

The code and the derivations were developed in interactive sessions with Claude (Anthropic);
the numerical experiments and judgement calls are the author's.

## Requirements

Python ≥ 3.10, numpy, scipy, matplotlib.  The production runs use CuPy with an NVIDIA GPU
(`--backend cupy`); the same C kernel compiles with gcc for a multicore CPU (`--backend c`),
and a pure-numpy reference backend exists for checking.  No other dependencies.

## Scripts

| file | what it does | paper |
|---|---|---|
| `hdym_core.py` | equations of motion, Yoshida-4 symplectic integrator with per-trajectory channel sub-stepping, CPU Lyapunov estimators, adaptive DOP853 confirmation of blow-up | Sec. 2–3 |
| `hdym_derive.py` | symbolic (sympy) reduction from the 3+1 theory to the 12-d system: Gauss law, off-diagonal equations, symmetric criticality, auxiliary-field Hamiltonian | Ref. [14] |
| `hdym_gpu.py` | transverse Lyapunov exponent λ⊥(E) with the whole integrator inside one kernel; `--ghost vec` (physical ghost) or `--ghost diag` (three independent scalar oscillators, per-component estimator) | Fig. 1, Table 1 |
| `hdym_spectrum.py` | multitaper spectrum of Hess V along pure-YM trajectories at E = 1; the u⁻⁶ tail and the tail constant A | Sec. 5.1 |
| `hdym_channel_check.py` | the channel-passage checks: adiabatic invariant J, parabolic x(t), chirp at 2|x|, P(J) ∝ J | Sec. 5.2 |
| `hdym_nu0.py` | channel entry rate ν₀ at E = 1 against 8π²/Z | Fig. 2 |
| `hdym_escape_map.py` | escape-time landscape on 2-d sections of initial conditions | Fig. 3a,b |
| `hdym_uncertainty.py` | uncertainty exponent α of the escape landscape | Fig. 3c |
| `make_fig1.py`, `make_figs.py` | regenerate the figures from the data files below | |
| `main.tex`, `fig_*.pdf` | the paper source and figures | |

## Data

| folder / file | contents |
|---|---|
| `transverse_vec_percomp.csv` | physical ghost, E = 0.008, 0.015, 0.02, 0.03, N = 2×10⁵ |
| `transverse_diag_percomp.csv` | diagonal ghost, E = 0.008, 0.015, 0.02, N = 2×10⁵ (per-component estimator) |
| `transverse_diag_percomp_hi.csv` | diagonal ghost, E = 0.05 … 1, N = 10⁵ (per-component estimator) |
| `transverse_E*.npz` | independent numpy re-run at E = 0.215, 0.3, 0.43, 0.6, 1 (open symbols in Fig. 1) |
| `nu0_E1.npz` | channel visits (depth, entry action) at E = 1 |
| `map256.npz`, `ghost256.npz`, `unc_*.npz` | escape maps and uncertainty-exponent data at E = 0.3, ε = 0.01 |

## Reproducing the main result

The four low-energy diagonal-ghost points of Table 1 (the parameter-free comparison with
0.0120 E^{9/4}):

```
python hdym_gpu.py --E 0.03 --N 100000 --ghost diag --backend cupy --t-burn-ym 200 --target-abs 3e-7 --t-block 5000 --tmax 300000 --out results_diag_percomp
python hdym_gpu.py --E 0.02 0.015 0.008 --N 200000 --ghost diag --backend cupy --t-burn-ym 200 --target-abs 3e-8 --t-block 20000 --tmax 1500000 --out results_diag_percomp
```

The physical ghost at the same energies, and the diagonal ghost at higher energies:

```
python hdym_gpu.py --E 0.03 0.02 0.015 0.008 --N 200000 --ghost vec --backend cupy --t-burn-ym 200 --target-abs 3e-8 --t-block 20000 --tmax 1500000 --out results_vec_percomp
python hdym_gpu.py --E 0.05 0.07 0.1 0.15 0.2 0.3 0.5 1.0 --N 100000 --ghost diag --backend cupy --t-burn-ym 200 --target-abs 3e-7 --out results_diag_percomp_hi
```

On an RTX 4080 SUPER these take roughly 0.5 h, 4 h, 4 h and 10 min respectively.  Each run
writes `<out>/transverse_gpu.csv` plus per-energy `history_E*.npy` (checkpoint history) and
`lam_i_E*.npy` (per-trajectory finite-time exponents); the three CSVs in this repository are
those files renamed.

Spectrum, channel statistics and ν₀ (CPU, pure Yang–Mills, no ghost):

```
python hdym_spectrum.py --N 2048 --T 20000 --t-burn 500 --workers 12
python hdym_channel_check.py --N 2048 --T 4000 --A-in 3
python hdym_nu0.py
```

Escape landscape and uncertainty exponent at E = 0.3, ε = 0.01:

```
python hdym_escape_map.py --E 0.3 --eps 0.01 --section ym    --n 256 --backend c
python hdym_escape_map.py --E 0.3 --eps 0.01 --section ghost --n 256 --backend c
python hdym_uncertainty.py --E 0.3 --eps 0.01 --section ym    --M 20000 --backend cupy --delta-min 1e-14 --delta-max 1e-4 --n-delta 11
python hdym_uncertainty.py --E 0.3 --eps 0.01 --section ghost --M 20000 --backend cupy --delta-min 1e-14 --delta-max 1e-4 --n-delta 11
```

Figures and paper, from the repository root: `python make_fig1.py && python make_figs.py && latexmk -pdf main.tex`.
(`make_figs.py` expects the Fig. 2–3 data files `nu0_E1.npz`, `map256.npz`, `ghost256.npz`, `unc_ym.npz`, `unc_ghost.npz` beside it.)

## A note on the diagonal-ghost estimator

With the off-diagonal Hessian entries dropped the three ghost components are independent
scalar oscillators.  Renormalising the norm of the three-vector returns the *largest* of three
finite-time exponents and overestimates the single-mode rate by a factor 1.5–1.7 at
E ≤ 0.03 (1.1 at E = 1).  `hdym_gpu.py --ghost diag` renormalises each component separately and
averages the three log-growths; the vector ghost (`--ghost vec`) uses the six-vector norm, which
is correct for the coupled system.  See the footnote in Sec. 3 of the paper.
