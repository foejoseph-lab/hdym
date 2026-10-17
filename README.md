# hdym — the Lee–Wick ghost of Yang–Mills mechanics

Code, data and paper source for *The Lee–Wick ghost of Yang–Mills mechanics is unstable at every
energy: parameter-free instability laws and an area-filling escape landscape* (J. Foe, 2026),
and for two follow-on investigations of open items listed in its Discussion.

The code and the derivations were developed in interactive sessions with Claude (Anthropic);
the numerical experiments and judgement calls are the author's.

## Layout

| folder | what it is |
|---|---|
| `paper/` | the current paper: `main.tex`, `main.pdf`, and every figure it includes |
| `ym/` | all Yang–Mills code behind the paper (`ym/README.md`: scripts, run commands, data index); results in `ym/data/`; working notes in `ym/notes/` |
| `islands/` | the regular (Dahlqvist–Russberg) island and the ghost on it, including a computer-assisted proof (`islands/README.md`) — the paper's open item on the regular islands |
| `bianchi9/` | the Weyl ghost in Bianchi IX minisuperspace of quadratic gravity (`bianchi9/README.md`) — the paper's open item on Bianchi IX; code at the top, run outputs in `bianchi9/results/` |
| `axion_su2/` | the isotropic SU(2) background of axion-SU(2) cosmology with the Lee–Wick ghost (`axion_su2/README.md`) — plan item 5, the new background; started 2026-10-12 |
| `dilaton_su2/` | the same isotropic SU(2) background with a gauge-kinetic (dilaton) coupling Θ(φ)F² instead of the axion's φFF̃ (`dilaton_su2/README.md`) — plan item 5b, the parity-even twin; started 2026-10-16 |

`TODO_held.md` holds the code-map to-do lists of 12 Oct (README-against-code mismatches in `ym/` and `bianchi9/`), held by decision.

`islands/`, `bianchi9/`, `axion_su2/` and `dilaton_su2/` are dated working notes, written after the paper text; the paper does
not yet incorporate them.

Conventions that keep this consistent: each result is described in one place (its folder's README
or note) and other files point there rather than restating it; when a later run settles or replaces
something, the earlier text keeps its wording and gets a dated `Status (YYYY-MM-DD):` line saying
what holds now; files a document relies on but that are not committed are listed as such.

## History

The 2026-10-01 draft (old title "… is malicious at every energy: a parameter-free E^{9/4}
instability law …") and the separate `revision-2026-10/` folder are superseded by `paper/main.tex`.
The repository as it stood before this layout is tagged `draft-2026-10-01`
(`git checkout draft-2026-10-01`).

## Quick start

Python ≥ 3.10 with numpy, scipy, matplotlib (sympy for derivation scripts; CuPy for GPU runs).

```
python ym/make_fig1.py && python ym/make_figs.py     # figures from committed data into paper/
cd paper && latexmk -pdf main.tex
```

`ym/make_fig_Q.py` additionally needs `ym/data/results_percomp.7z` unpacked to
`ym/data/results_percomp/`.
