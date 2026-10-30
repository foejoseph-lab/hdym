# Where the ghost grows, and the spectral formula rechecked (29 Oct 2026)

Two scripts, two logs:
- `hdym_jeans_sea.py` → `data/jeans_sea.txt`: 528 trajectories, τ = 20,000 at E = 1 scaled, 22 workers, 5125 s, run on
  the author's machine. Each sector has its own tangent: the vector ghost and the three c_j = 0 off-diagonal blocks.
  Growth is booked by channel depth.
- `hdym_random_oscillator.py` → `data/random_oscillator.txt` and `.npz`: 48 trajectories, run on Claude's side. These
  are spectra of individual Hessian entries, used in the paper's formula.

Naive ghost operator throughout, SU(2) colour-diagonal sector, hdym's GPU initial recipe (q, p standard normal, scaled
to the energy; burn-in 300 in scaled time).

**Origin.** Session notes of 27–29 Oct, not committed, parked a "Boltzmann–Jeans" idea: the ghost grows in the sea,
not in the channels, and its rate follows the spectrum of the drive restricted to sea time. If that spectrum falls
exponentially, the low-E law would be exponentially suppressed rather than a power. This note tests that idea.

**A correction to the session (29 Oct).** Asked why stochastic methods did not appear in the project, Claude said
they had been used only in that parked note. That was wrong. The paper's derivation is the randomly modulated
oscillator (eq:spectral, `paper/main.tex`: Arnold, Papanicolaou & Wihstutz 1986; Jarzynski 1993), with the u⁻⁶ tail
of the diagonal Hessian entry giving E^{9/4}. The u⁻² tail of xy gives E^{5/4} for the c_j = 0 blocks
(eq:fivequarters).

## 1. The sea-time idea is refuted

- **The sea spectrum does fall exponentially.** Chunks whose deepest excursion stays below 3 (41.5 % of the time) have
  S_sea ~ exp(−ω/0.415) over ω = 4–14: log-fit rms 0.169, against 1.02 for a power law. The full spectrum is a power
  law: mean entry ω^(−2.01), rms 0.012.
- **The off-diagonal blocks grow in the channels.** At E = 0.1, λ booked to channel time is 2.35 × 10⁻³ against
  4.8 × 10⁻⁴ in the sea; at E = 0.05, 9.7 × 10⁻⁴ against 5.6 × 10⁻⁵.
- **The vector ghost grows mostly in the sea at low E, but does not follow the sea spectrum either.**
  - About 64 % of its log-growth at E = 0.1 and 61 % at E = 0.05 comes from the 41.5 % of the time spent in the sea.
    That is the parked note's finding, which held for this sector only.
  - Its ratio to the sea-spectrum prediction still varies by a factor of 3.6 between E = 1 and 0.05.
  - The three-vector norm takes the largest of three finite-time exponents (`../README.md`), so it has no clean scalar
    formula.
- **The rate does not track the sea spectrum.** The ratio of λ to S_sea(2Ω)/(8Ω²) runs from 0.9 to 19.8 for the blocks
  between E = 1 and 0.05. Against the full spectrum it stays at about 2 (next section).

So the low-E laws are not set by the sea. The paper's remark that a bounded forcing of fixed bandwidth would give
λ ∝ e^(−cE^(−1/4)) describes what the sea alone would do. The channels' power-law tail is what prevents it.

## 2. eq:spectral with a fresh spectrum

**Tail exponents** (`data/random_oscillator.txt`, part 1; the law exponent is (p + 3)/4):
- H_xx = y² + z²: p = 6.08 (ω 6–14) and 6.09 (8–20), giving 2.27. The paper's multitaper slopes are 6.14 … 5.89.
- z²: p = 6.11 and 6.03.
- xy and 2xy: p = 2.01–2.09, giving 1.25–1.27.

**Diagonal ghost** (part 2): measured rate over eq:spectral with the full spectrum = 0.267, 0.254, 0.266 at
E = 0.02, 0.015, 0.008, and 0.267 at E = 0.1.
- The paper predicts the resonant share 2/7 = 0.286 and measures 0.276–0.290 at E = 0.03–0.008. This independent,
  smaller estimate (48 trajectories, Kaiser window, 20-unit chunks) is about 7 % lower.
- Between E = 0.15 and 0.05 the ratio dips to 0.17–0.30. That is below the u⁻⁶ regime (u < 4.3), as the paper notes
  for E = 0.05.
- Dropping the y- and z-channel chunks (the detuned power) moves the ratio toward 1: 0.47–0.83 at the lowest E. A
  chunk-level split is too coarse to reproduce the paper's exact 2/7 accounting.

**c_j = 0 blocks** (part 3), λ_inc of `data/jeans_sea.txt` over eq:spectral with S_h = S_z² + S_xy, no fitted
parameter:

| E | 1 | 0.5 | 0.3 | 0.2 | 0.1 | 0.07 | 0.05 |
|---|---|---|---|---|---|---|---|
| ratio | 0.63 | 0.70 | 0.87 | 0.93 | 1.07 | 1.02 | 1.06 |
| 0.0210 E^{5/4} / measured | 0.60 | 0.68 | 0.90 | 0.98 | 0.78 | 0.85 | 0.87 |

- With the measured spectrum itself, eq:spectral gives the block exponent to within 2–7 % at E ≤ 0.1.
- The closed form is 13–22 % low there. This confirms the paper's attribution of that excess to the knee of the xy
  spectrum at u ≤ 3.6 [grid], and extends the check to E = 0.07 and 0.05.
- Above E ≈ 0.3 the weak-driving approximation fails, as expected.

**Checks against existing numbers.** Blocks at E = 0.1: 1.49–1.55 × 10⁻³ against hdym's GPU 1.543 × 10⁻³ (offset-free
slope, t = 20,000). Vector ghost at E = 0.3: 0.016471 against 0.016506 (`data/transverse_E0.3.npz`). Energy drift
≤ 1.9 × 10⁻⁴.

## 3. Reading

Nothing here changes the paper's account; it confirms it with an independent spectrum and adds the c_j = 0 blocks to
the sectors where eq:spectral works quantitatively with the measured spectrum. The parked sea-time idea is closed:
negative. The vector ghost's preference for sea time at low E is a measured fact with no derivation here. For the
rewrite: possibly the blocks' 2–7 % (the author's call; nothing added to the watch-list).

```
python ym/hdym_jeans_sea.py --workers 22 --per 24 --tau 20000 --out jeans_out > ym/data/jeans_sea.txt 2>&1
python ym/hdym_random_oscillator.py --workers 22 > ym/data/random_oscillator.txt
```
