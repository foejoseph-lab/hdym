# Held to-do lists from the code map of 12 Oct 2026 (repository then at 8b866214)

Moved out of the context file on 16 Oct (reviewer's edit 4): repository facts, held and not started by
Joseph's decision of 12 Oct. Each line is a mismatch between a README and the code it describes, found by the
"README against code" pass of the code map. When one is taken up, fix it in its home (the folder README or the
script) and strike it here with a dated Status line.

## ym/ (plan item 3)

code that a paper claim or note depends on — `hdym_uncertainty.py` default thresholds
   2/5/10 vs the committed 5/10/20 (a README rerun would plot threshold 5 as "> 10"; feeds α);
   `hdym_rung1_lyap.py` `load_K` leaves c_j free and cannot run (NOTE_rung1 CPU spectra; the paper's
   table is from `hdym_rung1_gpu.py`, own K); `make_fig1.py` / `make_fig_Q.py` hard-code points and
   run lengths although the README says figures regenerate from `data/`. Wording — the
   `hdym_spectrum.py` row claims A; the `hdym_nu0.py` row claims 8π²/Z; numpy reference backend
   renormalises by the vector norm in both ghost modes; `dr_orbit_ghost.py`, `scan.py`,
   `interior.py` label digits ~1e-12 outside the certified box as the DR orbit; stale checker line
   counts. Held, not started (Joseph, 12 Oct).

## bianchi9/ (plan item 4)

six figures in `results/` with no writer; `bianchi9_bogoliubov.py`
   default `--out bogoliubov.npz` vs `bog.npz`, and no endpoint term removed vs the README's
   estimator; `bianchi9_symp.py` checks Ω without the constraint projection; Dirac axial coefficient
   ¼ (README, docstring) vs Σ/8 stored; `bianchi9_kasner.py` integrates x = 60→0.02 vs README
   200→0.01; `bianchi9_kasner_check.py` fixes θ = π vs README "all θ"; the Frobenius σ table file
   was overwritten (`bianchi9_kasner_transfer.py` still calls itself `bianchi9_kasner.py`);
   `bianchi9_six.py` has no Noether check and no β⁽⁶⁾ printout; `kasner_linear.pkl` is written by
   `bianchi9_kasner.py` and read by nothing.
