# Overnight queue, 21 Oct 2026 (Joseph's machine, 24 threads) — holes 1, 4, 12

Three runs, sequential, ~5 h in all. Each writes its npz after every row, so a killed run keeps what it has and
`--resume` continues it. From the repo root (Spyder: prefix each with `!"{sys.executable}" -B -u ` as in
`offdiag_escape.py`'s docstring; `ym/hdym_full_derive.py` must have been run once in this checkout).

```
# 1. hole 1: the phase ensemble at the Krein band (12 phases x 256 cells, |c| fixed, random admixture and
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
when present (to be adapted); hole 12's script does not exist yet.
