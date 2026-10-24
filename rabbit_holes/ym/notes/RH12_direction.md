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

| band | E | type | angle(nl end, lin end), median [p10, p90] | angle(nl end, unstable subspace), median | random-vector reference |
|---|---|---|---|---|---|
| Krein | 0.152 | quadruplet, |μ| 1.07 | 75° [27, 89] | **37°** | 60° (2-plane in 8 dims) |
| Krein | 0.155 | quadruplet, 1.27 | 76° [44, 89] | **36°** | 60° |
| Krein | 0.158 | quadruplet, 1.28 | 74° [44, 88] | **38°** | 60° |
| Krein | 0.161 | quadruplet, 1.14 | 75° [45, 87] | **32°** | 60° |
| param. | 0.226 | real pair, 1.48 | 73° [36, 87] | 66° | 69° (line in 8 dims) |
| param. | 0.237 | two real pairs, 1.80 | 78° [47, 89] | 53° | 60° |
| param. | 0.249 | real pair, 1.77 | 71° [46, 89] | 71° | 69° |
| param. | 0.262 | real pair, 1.38 | 75° [40, 88] | 72° | 69° |

**Result.** By the time the diagonal ghost is at 0.1 the off-diagonal sector no longer follows its own linear mode in
either band: the nonlinear end state is 70–78° from the linearly propagated one, with the same spread, for
quadruplets and real pairs alike. "Rotating vs directional" is therefore not observable at the nonlinear stage — the
diagonal ghost has taken the sector over (consistent with `RH4_loiter.md`: the diagonal ghost crosses
amplitude 1 before the off-diagonal mode does in 71–89 % of cells). The one residual of the type is weak but
systematic: the Krein end state keeps a memory of its unstable **plane** (32–38° from it against 60° for a random
direction, at all four energies), while the parametric end state is at random with respect to its unstable **line**
(66–72° against 69°). A quadruplet's two-dimensional unstable subspace survives the handover as a preferred plane; a
real pair's one-dimensional direction does not survive at all.

**Reading.** The distinction between the two mechanisms is a property of the linear stage (nodal points vs nodal
circle, `offdiag_escape_nodal.py`; the lemma; the segment localisation of RH3) and of the trigger (which sector
pumps, and at what rate, `RH4_loiter.md`); once the diagonal ghost is at order one the runaway is the diagonal ghost's,
and it does not remember which mechanism started it beyond a residual alignment with a plane. A referee asking
"do the two mechanisms differ once nonlinear" gets: in their trigger yes, in their execution no. Hole closed.
