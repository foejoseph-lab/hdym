# The gauge-reduced 4-dim ghost block of block 12 — built and verified on the Floquet side (22 Oct 2026)

Script `code/rh_reduced_block.py` (numerical; the Ostrogradsky form, the Gauss row and their orbit derivatives are
exact sympy objects, the reduction and integration are float). Logs `results/rh_reduced_block_cj1.txt` (E = 0.1,
0.3, 0.6; the E = 0.3 row rerun with the corrected Ω in `_cj1_E03.txt`), `_cj025.txt`, `_cj05.txt`; frozen spectra in
`results/rh_reduced_block_E*_cj*.npz`. This is the zeroth order `RH27_averaging.md` ended by asking for.

## Construction

Block 12 (8-dim state x = (u12, u21, u′, u″, u‴)) carries one linearised Gauss row g(t)·x = 0 (first class: a
random Gauss-clean initial condition keeps g·x at 10⁻¹⁵ over a period). With Ω(t) the Ostrogradsky symplectic form
of the block (c_j kept symbolic — the c_j = 1 form gives a wrong gauge direction at c_j ≠ 1, found by a wrong first result at c_j ≠ 1), the
gauge direction is v_G = Ω⁻¹g.

**Constant gauge fixing is not admissible on the DR orbit.** The unit gauge direction v_G(t) sweeps the
whole 8-space over a period (singular values of its samples 0.74, 0.39, 0.36, 0.30, 0.18, 0.16, 0.09, 0.09 at E = 0.1),
and the best fixed functional f has min |f·v_G| ≈ 10⁻³ somewhere on the orbit, so the gauge-fixed flow x′ = Kx + λv_G
with λ = −f·Kx/f·v_G blows up there (a multiplier 4 × 10⁸ in the first attempt). The moving gauge f(t) = v_G(t) is
used instead, with the connection term from the exact orbit derivatives: λ = −(v̇_G·x + v_G·Kx)/|v_G|², v̇_G =
Ω⁻¹(ġ − Ω̇v_G). Off-surface components g·x and v_G·x are then conserved (10⁻¹¹ over a period) and the surface
S(t) = {g·x = 0, v_G·x = 0} is invariant (leak 10⁻¹⁴).

**Second reduction.** The two leftover dimensions are the global colour-rotation zero mode z(t) = (−y, x, −p_y, p_x,
…) (an exact solution, Gauss-clean to 10⁻¹⁸) and its conserved charge c(t) = Ωz. Quotienting by z (made orthogonal to
v_G) and restricting to c·x = 0, with a second moving gauge term λ₂z, gives the 4-dim block.

## Checks (every E, every c_j in {0.25, 0.5, 1})

- CHECK A (6-dim): the four ghost multipliers of the full block are reduced multipliers to ≤ 8 × 10⁻¹¹; the other two
  are the unipotent pair at 1.
- CHECK A4 (4-dim): the reduced monodromy's spectrum *is* the four ghost multipliers, to ≤ 8 × 10⁻¹¹, including the
  c_j = 0.25, E = 0.6 case with |μ| = 34 and 16.5 (the double +1 fold).
- The reduced 4 × 4 Ω is nondegenerate and the reduced monodromy is symplectic for it to 3 × 10⁻¹³.

So the 4-dim reduced block is exact (not a model): it is the ghost dynamics with the constraint and symmetry sectors
removed, and anything computed from it is a statement about the theory.

## First outputs

**Frozen spectrum along the orbit** (eigenvalues of the restricted generator at fixed t; the gauge terms make it not
exactly Hamiltonian, real parts ≤ 0.1, so read the imaginary parts as instantaneous frequencies), c_j = 1:

| E | turning point t = 0 | bounce t = T/4 | t = T/2 |
|---|---|---|---|
| 0.1 | 0.71, 1.70 | 0.81, 1.44 | 0.71, 1.70 |
| 0.3 | 0.34, 1.29 | 0.62, 1.57 | 0.33, 1.29 |
| 0.6 | 0.13, 1.02 | 0.73, 1.64 | real pair ±0.17 and 1.02 |

Two modes: a light one that goes soft at the turning points (0.13–0.7, hyperbolic at E = 0.6, t = T/2) and a heavy
one at 1.0–1.7. This is the sweep of `rh5_instfreq.png` seen from the Hamiltonian side, with one disagreement with
that reading, not resolved: the frozen spectrum puts the light mode's softening at the **turning points** (x extremal,
y = 0), not at the bounce (x = 0), where both modes are stiffest, while the instantaneous-frequency estimator of
`rh5_instfreq.py` puts the dips at the bounce (`RH27_averaging.md`, item 10). The quasi-static spectral line is the light mode; the main line near 1 is the heavy mode dressed
down by the sweep.

**c_j family** (frozen, E = 0.6, t = T/2): c_j = 1 real pair ±0.17; c_j = 0.5 real ±0.06/0.23; c_j = 0.25 real
±0.12/0.29 — the light mode is hyperbolic at the turning points for every c_j, and the static instability of the
c_j = 0.25 family is **not** a mean-stiffness effect: the averaged generator at c_j = 0.25 has purely imaginary
spectrum (0.59, 0.66 at E = 0.6; 0.82, 0.85 at 0.3; 0.98, 1.07 at 0.1) while the exact multipliers are 34 and 16.5.
Of the two alternatives — negative mean stiffness, or transient-negative stiffness (a Kapitza-type effect the averaged
block cannot see) — the second is the one that holds: the +1 instability at low c_j is the light mode being hyperbolic for a fraction of each period,
and the physical c_j = 1 is where that fraction is smallest. At c_j = 1 the frozen spectrum is hyperbolic at E = 0.6
only; at c_j = 0.25 at every E (fraction of the period with a real eigenvalue: 1.0 in all logs, because the gauge
terms give small real parts everywhere — the hyperbolic *light-mode* fraction is the quantity to extract next).

**Averaged block (parallel-transported frame).** Not Floquet-accurate: at c_j = 1 the averaged
frequencies (1.16/1.61, 1.43/1.31, 1.12/1.50 at E = 0.1/0.3/0.6) give ωT mod 2π at 1.14/1.82, 0.05/0.73, 0.62/0.14 π
against the measured 0.10/0.26, 0.40/0.56, 0.00/0.77 π. The modulation is too strong for first-order averaging even on
the exact reduced block; the frozen spectrum and the bounce transfer matrix are the objects to use.

## What this settles and what it opens

- Settled: the 4-dim block exists, is exact, and carries two modes — light (swept, soft at the turning points) and
  heavy (≈ 1–1.7). The "near 1" line of the spectra is the heavy mode; the quasi-static line the light one.
- Settled: the low-c_j static instability is transient hyperbolicity of the light mode, not negative mean stiffness.
- Open, now well-posed: (i) the hyperbolic fraction of the light mode per period as a function of (E, c_j) — one
  number per case, from the saved frozen spectra; (ii) the bounce/turning-point transfer matrix: integrate the 4-dim
  reduced flow across a window around t = 0 (where the light mode is soft) and read its sector-mixing entries
  (`RH27_averaging.md`, part E), with the soft point located at the turning point; (iii) the Krein width at
  [0.151, 0.163] across c_j as the one first-order test; (iv) a proper Hamiltonian form H_red(t) = −Ω_red K_red in a
  symplectic frame (the current frame is orthonormal, which is why the frozen generator is not exactly Hamiltonian).
