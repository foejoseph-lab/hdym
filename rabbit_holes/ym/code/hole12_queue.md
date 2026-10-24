# Hole 12 — the snapshot runs (queue of 22 Oct; run and closed the same day)

Both runs done, logs committed 23 Oct (`islands/results/esc_E_snap_krein.*`, `esc_E_snap_param.*`); analysis `rh12_direction.py`;
hole closed in `../notes/RH12_direction.md`. Kept as the record. (Operational note that applied at the time: the kernel is rebuilt
at every launch, so a run with the snapshot kernel must not start while a job on the previous kernel is in flight in the same checkout.)

The kernel now snapshots the full 36-component state at each crossing event (|w| = 0.1, |w| = 1, each block norm = 1),
stored as `snap` (N, 5, 36) float32 in every section's npz.  Hole 12's question -- does a Krein quadruplet's runaway
ROTATE in the off-diagonal state space while a real pair's runs along a fixed direction -- is answered by the angle
between the block-12 off-diagonal direction at the |u_12| = 0.1-level... (not recorded) -- so the comparison is between
the snapshots at |w| = 0.1 and |w| = 1 and at |u_12| = 1, against the linear prediction from the Floquet mode propagated
over the same interval.  Two runs, both --block 12 so the sector is clean:

```
!"{sys.executable}" -B -u islands\code\offdiag_escape.py --section E --Emin 0.152 --Emax 0.161 --nE 4 --nseed 128 --block 12 --tag esc_E_snap_krein --resume >> islands\results\esc_E_snap_krein.txt 2>&1
!"{sys.executable}" -B -u islands\code\offdiag_escape.py --section E --Emin 0.226 --Emax 0.262 --nE 4 --nseed 128 --block 12 --tag esc_E_snap_param --resume >> islands\results\esc_E_snap_param.txt 2>&1
```

(Krein band of block 12 at 0.152-0.161, parametric band at 0.223-0.265; 4 energies x 128 seeds each; ~1 h each on 24
threads at these low rates.)  Analysis (`rh12_direction.py`): for each trajectory, the unit
off-diagonal direction d(t) of block 12 at the three snapshots; the angle between d(|w|=0.1) and d(|w|=1); and the
same angle for the linear Floquet solution started from the seed over the same time interval.  Rotating runaway:
the nonlinear angle tracks the linear one (which winds at the Krein angle per period); directional: both near 0 or pi.
