# Hole 9 — the island in the full nine-component system (run Friday 23 Oct while at work; both resumable)

Two runs, ~10–20 min each on 24 threads (1800 trajectories, tmax 5000):

```
!"{sys.executable}" -B -u islands\code\island_full.py --E 0.3 --eps 1e-8 --tag island_full_E0.3_eps1e-8 --resume >> islands\results\island_full_E0.3_eps1e-8.txt 2>&1
!"{sys.executable}" -B -u islands\code\island_full.py --E 0.3 --eps 1e-8 --delta-od 0 --tag island_full_E0.3_eps1e-8_od0 --resume >> islands\results\island_full_E0.3_eps1e-8_od0.txt 2>&1
```

The first seeds the six off-diagonal components at 10⁻⁸ (Gauss-projected; three decades under the diagonal ghost's
amplitude 7.7 × 10⁻⁵); the second leaves them to the integrator's truncation (the "exactly diagonal" control).
hdym_island's result at the same E, ε (rabbit hole 10): core shells r ≤ 10⁻⁴ escape 0 of 800, sea control 100 %
with median ≈ 500. The question is whether the core shells now escape, through which sector, and when (the block-23
parametric rate on the orbit is 0.0125: from 10⁻⁸ to order one is ln(10⁸)/0.0125 ≈ 1500 time units, so an emptied
island should show t_esc ≈ 1500–3000 in the core shells with block 23 the first sector to reach norm 1).
Optional third, if the first two finish early: `--eps 1e-10` (the paper's other island energy is E = 1; not needed).

Then send the two .txt logs (and the .npz if convenient); the reading and the closing note follow.
