# hdym — the Lee–Wick ghost of Yang–Mills mechanics

Code, data and paper source for *The Lee–Wick ghost of Yang–Mills mechanics is unstable at every energy:
parameter-free instability laws and an area-filling escape landscape* (J. Foe, 2026), and for the follow-on
investigations of the open items in its Discussion. Rewritten 23 Oct 2026 from the state at that date.

The code and the derivations were developed in interactive sessions with Claude (Anthropic); the numerical
experiments and judgement calls are the author's. The paper text in `paper/main.tex` is that of 10 Oct 2026; a
rewrite is planned (25 Oct or the first free weekend after) and the paper has not been submitted anywhere.

## Layout

| folder | what it is | status |
|---|---|---|
| `paper/` | the paper: `main.tex`, `main.pdf`, and every figure it includes | text of 10 Oct; rewrite pending |
| `ym/` | all Yang–Mills code behind the paper (`ym/README.md`: scripts, run commands, data index, the roundoff study); results in `ym/data/`, notes in `ym/notes/`, the shared C model and the nonlinear escape kernel in `ym/csrc/` | complete for the paper as written |
| `islands/` | the regular (Dahlqvist–Russberg) island: a computer-assisted proof that the orbit is elliptic with nonzero twist, the diagonal and off-diagonal ghost on it, and the nonlinear escape from it (`islands/README.md`) — the paper's open item on the regular islands | proof done 12 Oct; the island is not a refuge in the full theory (12 and 23 Oct) |
| `rabbit_holes/ym/` | twelve side questions of the off-diagonal ghost on the DR orbit, worked 21–23 Oct (`rabbit_holes/ym/README.md`): a selection rule that is a theorem, a ladder that is a measurement, the exact gauge-reduced ghost block, and the nonlinear loiter; `notes/FOR_THE_REWRITE_2026-10-25.md` is the watch-list for the paper rewrite | eleven closed, one run and unread, four open items |
| `bianchi9/` | the Weyl ghost in Bianchi IX minisuperspace of quadratic gravity (`bianchi9/README.md`) — the paper's open item on Bianchi IX | worked 5–10 Oct; closure written, two-ghost case unfinished |
| `axion_su2/` | the isotropic SU(2) background of axion-SU(2) cosmology with the Lee–Wick ghost (`axion_su2/README.md`), a new background | worked 12–17 Oct; coupled-axion stage open |
| `dilaton_su2/` | the same isotropic SU(2) background with a gauge-kinetic (dilaton) coupling Θ(φ)F², the parity-even twin (`dilaton_su2/README.md`) | 16–17 Oct, closed |

`islands/`, `rabbit_holes/`, `bianchi9/`, `axion_su2/` and `dilaton_su2/` were written after the paper text; the paper
does not yet incorporate them. What of them belongs in the rewrite, and at what strength, is in
`rabbit_holes/ym/notes/FOR_THE_REWRITE_2026-10-25.md` for the Yang–Mills work; the other three backgrounds have their
own closing paragraphs.

Conventions: each result is described in one place (its folder's README or note) and other files point there;
every number in a README traces to a committed log, script or figure named beside it; files a document relies on
but that are not committed are listed as such; dates are kept where they say when a result was obtained. Held
README-vs-code mismatches found on 12 Oct are listed in `ym/README.md` and `bianchi9/README.md` under that heading.

## History

The 1 Oct 2026 draft (old title "… is malicious at every energy: a parameter-free E^{9/4} instability law …") and
a separate `revision-2026-10/` folder are superseded by `paper/main.tex`; the repository as it stood before the
present layout is tagged `draft-2026-10-01` (`git checkout draft-2026-10-01`). The layout `paper/`, `ym/`,
`bianchi9/`, `islands/` dates from 10 Oct; `axion_su2/` from 12 Oct, `dilaton_su2/` from 16 Oct, `rabbit_holes/`
from 21 Oct. All markdown files were rewritten on 23 Oct from the state at that date (the previous files had grown
as dated notebooks with appended status lines); the history of every file is in git.

## Quick start

Python ≥ 3.10 with numpy, scipy, matplotlib (sympy for derivation scripts; CuPy for GPU runs).

```
python ym/make_fig1.py && python ym/make_figs.py     # figures from committed data into paper/
cd paper && latexmk -pdf main.tex
```

`ym/make_fig_Q.py` additionally needs `ym/data/results_percomp.7z` unpacked to `ym/data/results_percomp/`.
