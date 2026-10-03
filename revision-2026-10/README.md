# revision-2026-10 — material added to the paper after the 2026-10-01 draft

Everything the revised manuscript (`main_merged.tex`, the version to be posted) uses beyond the
files in the repository root. The root files are the 2026-10-01 draft and its data; nothing there
was changed. Start with `HANDOFF_2026-10-02.md`.

| what | files |
|---|---|
| derivations (sympy, every constant checked) | `hdym_analytic.py` — general-n law, 4/(9n−4) share, Z₂ elliptic form, c_j = 0 law 0.0210 E^{5/4} |
| two-field model (Sec. "test on the two-field model") | `hdym_2d.py` (kernel), `hdym2d_nu0.py`, `hdym2d_spectrum.py`, `nu0_2d_*.npz`, `spec_2d_*.npz`, `results_2d/` (GPU per-trajectory exponents), `results_2d_cpu/`, `NOTE_2d.md` |
| off-diagonal sector (Sec. "The off-diagonal ghost sector") | `hdym_rung1.py` → `rung1_system.pkl` (linearisation), `hdym_rung1_lyap.py` (CPU spectra, `rung1_*.npy`), `hdym_rung1_gpu.py` (kernel: `--cj`, growth split, (D,J) landscape), `hdym_spectrum_xy.py` + `spec_xy_E1.npz`, `results_rung1/` (GPU per-trajectory exponents, c_j = 1 and 0), `NOTE_rung1.md` |
| figure 5 | `make_fig_Q.py`, `fig_Q.pdf` (needs `results_percomp.7z` from the root unpacked to `results_percomp/`, and `results_rung1/`) |
| landscape figure (not in the paper) | `make_fig_landscape.py`; `hdym_gpu_landscape.py` is the root `hdym_gpu.py` plus the (D,J) histogram output, results bit-identical |
| paper | `additions.tex` (the blocks, with insertion points), `main_merged.tex`/`.pdf` (merged, reviewed) |

GPU summary lines for every run are quoted in the two NOTE files; the CSVs (`transverse_2d.csv`,
`rung1_offdiag.csv`) are in the author's `results_2d/` and `results_rung1/` folders.
