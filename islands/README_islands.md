# islands/ — the Dahlqvist–Russberg island, the ghost on it, and a computer-assisted proof

Session of 8–9 October 2026 (developed interactively with Claude; see the paper's
acknowledgement). Everything here concerns the stable periodic orbit of Yang–Mills mechanics
found by Dahlqvist & Russberg (1990), with initial data from Contopoulos & Harsoula,
arXiv:2302.12071, and what the Lee–Wick ghost does on and around it.

## Results in one paragraph

The DR orbit at E = 1/2 exists and is elliptic (Hénon index in [−0.96097, −0.96088]) and
non-resonant to order 4 — **proved**, by a certificate checked in exact rational arithmetic.
The linearised ghost on it is **Floquet-unstable** at E = 1/2 (proved, real multiplier off the
unit circle) and, by float scan, for all E ≳ 0.315 plus narrow tongues below; between tongues
it is Floquet-stable. Its Birkhoff twist is τ ≈ −2.14×10⁵ (exact formula, float evaluation;
rigorous enclosure pending), so Moser's theorem applies once τ is certified: the island exists.
In the full 12-d Lee–Wick system with a tiny ghost seed (ε = 10⁻¹⁰) the island is an
infinite-escape-time hole at E = 0.3 and a fast smooth plateau at E = 1, with a sticky halo
that escapes slower than the sea; at the paper's ε = 0.01 the ghost's back-reaction (|w| ≈ 0.1
≫ island width 2×10⁻⁴) erases the island entirely. The same holds in the exact 2-dof planar
subsystem, so nothing about the three-field lift changes the 2-dof analysis.

## Layout

```
code/
  dr_orbit_ghost.py    DR orbit, Hénon index, coarse ghost Floquet scan in E       (float)
  scan.py              fine Floquet scan: python scan.py E_min E_max N             (float)
  interior.py          base and ghost Lyapunov exponents inside the island         (float)
  cap_dr_orbit.py      UNTRUSTED certificate generator (mpmath intervals, Lohner/Taylor)
  cert_check.py        TRUSTED checker: pure Python Fractions, one function per lemma
  twist_symbolic.py    exact sympy formula for the Birkhoff twist + float evaluation
  twist_estimate.py    independent twist estimate from rotation numbers
  hdym_island.py       escape times / FTLE around the DR orbit in the 12-d system
                       (imports hdym_escape_map.py, hdym_precision.py from the repo root)
  run_island.sh        the four production runs
results/
  cert16.json          certificate, md5 690878da67b071b6ce277aea0fb02e98
  cert16_check_output.txt   checker output (PROVED)
  twist_sym.log
  island_E*_eps*.log   six hdym_island runs (E = 1, 0.3; ε = 1e-10, 0.01; two planar)
```

## Reproduce the proof

```
python3 code/cap_dr_orbit.py 16 1000 1e-11 results/cert16.json   # ~5 min; untrusted
python3 code/cert_check.py results/cert16.json --ghost            # ~35 min; this is the proof
```

The checker verifies, step by step, with outward-rounded dyadic rationals: a Picard rough
enclosure per step (lemma PIC), the Lagrange remainder of the Taylor step (LAG), the
mean-value Jacobian (MVT), a Neumann-series inverse (NEU), set inclusion in the Lohner frame
(step), the accumulated Jacobian in a moving frame (FRAME), Krawczyk's inclusion (KRA), and
the trace/symplectic stability lemmas (T, SYMP). Lemma names are the function names. The
generator may be wrong in any way; a bad certificate fails to check. What is trusted: the
checker (≈400 lines, unreviewed as of 9 Oct), Python integer arithmetic, and `math.isqrt`.

Design target: Lean. The lemma shapes were chosen so each is an elementary statement about
polynomial ODE enclosures or rational linear algebra; Moser's twist theorem would remain an
import.

## Status / open

- [ ] Twist enclosure: add 2nd/3rd-order variational series to generator and checker (VAR2, VAR3),
      evaluate the `twist_symbolic.py` formula in intervals. Then the island is a theorem modulo Moser.
- [ ] Second certificate at E = 0.3 (ghost Floquet-STABLE on the orbit).
- [ ] Independent review of `cert_check.py`.
- [ ] ε-crossover scan at E = 0.3: `hdym_island.py --eps 1e-9 … 1e-6`.
- [ ] 3-dof: the DR orbit's z-direction Floquet trace (−1.53, float) is not certified; the
      island is sealed in 2-dof and leaks by Arnold diffusion in 3-dof (seen: 2.7 % at r = 2e-4).

## Conventions

`hdym_island.py` uses the repo's conventions verbatim: state (p, w, P, W), q = p − w,
E = |P|²/2 + V(q), ghost |w|² + |W|² = 2εE, p = q + w, P = P_q − W, Yoshida-4 with dt = 0.005
and channel sub-stepping, escape at any |component| > 10⁴. The DR point in these
coordinates: q* = λ(3.14640122769837, 0, 0), P* = λ²(0.00179319341881, √(1 − p_x*²), 0),
λ = (E/½)^{1/4}; copies 0–23 by cyclic permutation, 90° in-plane rotation, time reversal.
