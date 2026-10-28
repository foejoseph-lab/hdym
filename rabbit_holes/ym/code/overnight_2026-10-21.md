# Overnight queue, 21 Oct 2026 (24 threads) — holes 1, 4, 12

Run on the night of 21–22 Oct; all CHECK lines PASS; logs and `.npz` committed in `islands/results/`
(`esc_phase_0156`, `esc_orient_030_n64`, `esc_E_cross`); read in `../notes/RH1_phase.md` and `../notes/RH4_loiter.md`.
Hole 12 needed the state snapshots added afterwards (`hole12_queue.md`). The analysis scripts named at the end
became `rh_overnight_2026-10-22.py`. Followed on 22 Oct by `day_2026-10-22.bat` (`esc_E_lowE_long`, read by `rh4_lowE_long.py`).
Kept as the record of what was run and why.

Three runs, sequential, ~5 h in all. Each writes its npz after every row, so a killed run keeps what it has and
`--resume` continues it. **For an unattended night use `overnight_2026-10-21.bat`** (same three commands, every one
with `--resume` and `>>` on its log): after a reboot, launch the .bat again and each run continues from its last
completed row — exposure is one row, 2.5–5 min. Registering the .bat in Task Scheduler "at log on" makes the
resume automatic. Check `islands/results/overnight_2026-10-21.log` (the .bat's own log; not committed) for the start/done stamps and exit codes.

A note on the kernel edit and floats: the crossing-time code reads the accepted state and writes a separate array;
it touches no state, no step size and no control flow, so trajectories are bit-identical to the previous kernel's
(the 48-copy determinism CHECK and the scipy cross-check still run at every launch and would say otherwise). The only
new arithmetic is `sqrt` on non-negative arguments, already used by `escaped()` and the checkpoints. From the repo root (Spyder: prefix each with `!"{sys.executable}" -B -u ` as in
`offdiag_escape.py`'s docstring; `ym/hdym_full_derive.py` must have been run once in this checkout).

```
# 1. hole 1: the phase ensemble at the Krein band (12 phases x 256 cells, |c| fixed (it was not: 0.56-1.21 with alpha, ../notes/RH1_phase.md), random admixture and
#    diagonal seed per cell).  ~1 h.
python islands/code/offdiag_escape.py --section phase --E 0.156 --block 12 --nalpha 12 --nseed 256 --tag esc_phase_0156 > islands/results/esc_phase_0156.txt 2>&1

# 2. hole 12: a parametric-band map that KEEPS its growth curves (chk, float32) and crossing times, for the
#    rotating-vs-directional runaway comparison against run 1's Krein trajectories.  64^2 cells, ~2.7 h.
python islands/code/offdiag_escape.py --section orient --E 0.30 --n 64 --tag esc_orient_030_n64 > islands/results/esc_orient_030_n64.txt 2>&1

# 3. hole 4: section E with first-crossing times (|w| = 0.1, 1; each sector norm = 1), 20 energies x 64 seeds.  ~1 h.
python islands/code/offdiag_escape.py --section E --nE 20 --nseed 64 --tag esc_E_cross > islands/results/esc_E_cross.txt 2>&1
```

What to look for in the logs: `ALL CHECKS PASS` at the end of each (the kernel's scipy cross-check, determinism
check and orbit-closure check are unchanged by the kernel edit; they passed in the cloud smoke test of 21 Oct).
Sizes: esc_phase_0156.npz ≈ 3072 x 60 x 6 x 4 B ≈ 4.4 MB; esc_orient_030_n64.npz ≈ 6 MB; esc_E_cross.npz ≈ 4 MB.

Analysis scripts to run on them afterwards: `rabbit_holes/ym/code/rh1_phase.py` can be pointed at the phase
run (to be adapted: it currently reads the orient maps); `rh4_loiter.py` reads `chk` and will read `tcross`
when present (to be adapted); hole 12's script does not exist yet (written 22 Oct: rh12_direction.py).
