# Rabbit hole 12 — nonlinear character of a Krein vs a parametric runaway (22 Oct 2026)

Runs `islands/results/esc_E_snap_krein.*` (block 12's Krein band, E = 0.152–0.161, 4 × 128, all checks PASS, 297 s on 24 threads)
and `esc_E_snap_param.*` (its parametric band, E = 0.226–0.262, 4 × 128, 822 s), both with the kernel's new state
snapshots at the crossing events (patch 22). Script `code/rh12_direction.py`, log `results/rh12_direction.txt`,
per-energy arrays `results/rh12_esc_E_snap_*_E*.npz`. 32–48 trajectories per energy propagated linearly.

**Question.** Once the runaway is nonlinear, does a quadruplet's rotate in the off-diagonal state space while a real
pair's runs along a fixed direction? Measured between the two snapshots |w| = 0.1 and |w| = 1 (the diagonal ghost's
amplitude; 1–3 periods apart at these energies, 16 at the band edge E = 0.152): the block-12 direction d at each,
the angle between them, the same angle for the linear Floquet propagation of the first snapshot over the same
interval, the angle between the nonlinear and the linear end states, and the angle of the nonlinear end state to
the unstable Floquet subspace at that phase of the orbit (Euclidean angles in (u, u′, u″, u‴); not canonical, but
the same convention for both bands).

**Corrected 28 Oct 2026** (rerun, `results/rh12_direction.txt`, wall 2030 s, cloud). For a Krein quadruplet the script
built the unstable subspace from both conjugate eigenvectors, whose real and imaginary parts span only a 2-plane; the QR
step then returned two further, essentially arbitrary columns, so the "unstable plane" was four-dimensional and the
angle to it came out too small. The script now keeps only independent directions (SVD), and prints two nulls: the
exact median for an isotropic direction (cos²θ ~ Beta(k/2, (8 − k)/2), Euclidean metric assumed) and a phase-shifted
null (the same end state against the subspace at a wrong phase, uniform in 0.15–0.85 of the period; no metric
assumption). The angles to the linear end state are unaffected by the bug; their small differences from the table
below come from the trajectory count (now up to `--nmax` 64).

| band | E | type | trajectories | angle(nl end, lin end), median [p10, p90] | angle(nl end, unstable subspace), median [p10, p90] | phase-shifted null | isotropic null |
|---|---|---|---|---|---|---|---|
| Krein | 0.152 | quadruplet, \|μ\| 1.07 | 45 | 76° [42, 88] | **52°** [40, 68] | 59° [43, 77] | 63° (k = 2) |
| Krein | 0.155 | quadruplet, 1.27 | 64 | 71° [44, 88] | **53°** [31, 71] | 53° [34, 64] | 63° |
| Krein | 0.158 | quadruplet, 1.28 | 63 | 72° [41, 87] | **48°** [36, 72] | 54° [34, 75] | 63° |
| Krein | 0.161 | quadruplet, 1.14 | 64 | 74° [45, 88] | **52°** [30, 76] | 53° [36, 71] | 63° |
| param. | 0.226 | real pair, 1.48 | 64 | 73° [41, 88] | **67°** [44, 85] | 66° [42, 86] | 75° (k = 1) |
| param. | 0.237 | two real pairs, 1.80 | 64 | 79° [54, 89] | **56°** [37, 78] | 54° [32, 75] | 63° (k = 2) |
| param. | 0.249 | real pair, 1.77 | 64 | 72° [49, 88] | **71°** [45, 87] | 63° [42, 85] | 75° (k = 1) |
| param. | 0.262 | real pair, 1.38 | 64 | 78° [43, 88] | **72°** [47, 85] | 66° [38, 83] | 75° (k = 1) |

As of 22 Oct (subspace angles wrong for the quadruplets, references not computed; kept as the record):

| band | E | type | angle(nl end, lin end), median [p10, p90] | angle(nl end, unstable subspace), median | random-vector reference |
|---|---|---|---|---|---|
| Krein | 0.152 | quadruplet, |μ| 1.07 | 75° [27, 89] | **37°** | 60° (2-plane in 8 dims; not computed in the script) |
| Krein | 0.155 | quadruplet, 1.27 | 76° [44, 89] | **36°** | 60° |
| Krein | 0.158 | quadruplet, 1.28 | 74° [44, 88] | **38°** | 60° |
| Krein | 0.161 | quadruplet, 1.14 | 75° [45, 87] | **32°** | 60° |
| param. | 0.226 | real pair, 1.48 | 73° [36, 87] | 66° | 69° (line in 8 dims; not computed in the script) |
| param. | 0.237 | two real pairs, 1.80 | 78° [47, 89] | 53° | 60° |
| param. | 0.249 | real pair, 1.77 | 71° [46, 89] | 71° | 69° |
| param. | 0.262 | real pair, 1.38 | 75° [40, 88] | 72° | 69° |

**Result.** By the time the diagonal ghost is at 0.1 the off-diagonal sector no longer follows its own linear mode in
either band: the nonlinear end state is 70–78° from the linearly propagated one, with the same spread, for
quadruplets and real pairs alike. "Rotating vs directional" is therefore not observable at the nonlinear stage — the
diagonal ghost has taken the sector over (consistent with `RH4_loiter.md`: the diagonal ghost crosses
amplitude 1 before the off-diagonal mode does in 71–89 % of cells). The residual reported on 22 Oct (the Krein end
state 32–38° from its unstable plane against 60° for a random direction) was an artefact of the four-dimensional
"plane" (see the correction above). Corrected, the Krein end states sit 48–53° from the plane against 53–59° at a wrong
phase, and the parametric ones 67–72° from the line against 63–66°: within the spread of the nulls in both bands
[grid]. Neither mechanism's unstable subspace is remembered after the handover.

**Reading.** The distinction between the two mechanisms is a property of the linear stage (nodal points vs nodal
circle, `offdiag_escape_nodal.py`; the lemma; the segment localisation of RH3) and of the trigger (which sector
pumps, and at what rate, `RH4_loiter.md`); once the diagonal ghost is at order one the runaway is the diagonal ghost's,
and it does not remember which mechanism started it (the 22 Oct "residual alignment with a plane" withdrawn 28 Oct
2026). A referee asking
"do the two mechanisms differ once nonlinear" gets: in their trigger yes, in their execution no. Hole closed.
