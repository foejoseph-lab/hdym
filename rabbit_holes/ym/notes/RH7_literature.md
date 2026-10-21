# Rabbit hole 7 — the books (21 Oct 2026, half an hour of searching)

Verdict: **the Lemma is the Swift–Wiesenfeld mechanism (1984), well known in the dissipative setting; what is
ours is the symplectic multiplicity statement and its identification with the ghost's Krein collision.** The
paper should cite, not claim.

## What exists

- **J. W. Swift and K. Wiesenfeld, "Suppression of period doubling in symmetric systems", Phys. Rev. Lett. 52,
  705 (1984).** The original. A symmetry S with the orbit satisfying x(t + T/2) = S x(t) makes the period map
  the square of the half-period map composed with S, P = P̃². A multiplier −1 of P needs ±i in P̃, which for a
  single real eigenvalue is impossible; period doubling is therefore suppressed for symmetric orbits except in
  "extraordinary cases" (a complex pair of P̃ at ±i), which they show cannot occur in a large class of dissipative
  systems (driven damped pendulum, Lorenz). Their even/odd-harmonic signature (symmetric response has only odd
  harmonics; even ones appear when the symmetry breaks) is the parity statement of our step 1.
- **J. S. W. Lamb and J. A. G. Roberts, "Time-reversal symmetry in dynamical systems: a survey", Physica D 112,
  1–39 (1998)**, and **J. S. W. Lamb, "k-symmetry and return maps of space-time symmetric flows", Nonlinearity 11
  (1998)**: the general framework — a space-time (spatio-temporal) symmetry of order k forces the return map to
  be a k-th root of the monodromy, M = Gᵏ, with the classification of which local bifurcations survive. Our S̃
  with the T/2 shift is the k = 2 case.
- **U. Frauenfelder and A. Moreno, "On doubly symmetric periodic orbits", Celest. Mech. Dyn. Astron. 135, 20
  (2023), arXiv:2301.01803.** The symplectic cousin: for orbits symmetric under two commuting anti-symplectic
  involutions, Theorem A says a doubly symmetric orbit cannot be negative hyperbolic in dimension four, Corollary
  D that it cannot period-double there, and they note **Corollary D fails in dimension six** (a numerical
  planar-to-spatial period doubling). Their mechanism is the B-sign of the reduced 2 × 2 monodromy, not the
  square root, and their symmetry is reversing (anti-symplectic) where ours is a canonical symmetry composed with
  a half-period shift — different hypotheses, same conclusion in the two-dimensional reduced case, and the same
  escape hatch in higher dimension: once there are two pairs, the "extraordinary case" (N at ±i) is available.
- A current example of the mechanism used exactly as we use it, with the parity selection written out:
  A. C. Marti and E. D. Leonel, "Antiperiodicity in the Duffing–Holmes oscillator: symmetry origin, parity
  selection, and transitions", arXiv:2607.05604 (2026) — dissipative, so their ±i case is excluded by det < 1.

## What is ours, after reading

1. The **symplectic multiplicity form**: with N real symplectic, N's ±i pair comes as a quadruplet
   {ir, −ir, i/r, −i/r}, so a real negative multiplier of M off the circle has multiplicity ≥ 2 — the Swift–
   Wiesenfeld "extraordinary case" is not excluded in a Hamiltonian system (no det < 1) but is forced to be a
   *two-pair* event. In the dissipative literature the extraordinary case is a curiosity; here it is the
   generic −1 instability.
2. Its **identification with the Krein collision** of the ghost sector: the two-pair −1 event of M is a Krein
   collision of N at ±i between pairs of opposite signature (B005 measured the signs), and in the decoupled
   symmetric limit it is the sector's period-doubling tongue doubled. "Krein bands are dressed period
   doubling" is, as far as this search went, not in the literature.
3. The application to Ostrogradsky / higher-derivative ghost dynamics, where the sector has two pairs by
   construction (a massive ghost plus the gauge/rotation modes), so the Frauenfelder–Moreno dimension-six
   escape is the normal case rather than the exception.

## For the paper

Cite Swift–Wiesenfeld for the suppression, Lamb–Roberts / Lamb for the k-th-root structure, Frauenfelder–Moreno
for the symplectic statement and its failure beyond two degrees of freedom; state (1)–(2) as a short lemma with
proof, in the off-diagonal section. Not a standalone paper on its own; with the ladder of hole 2 and a reduction
that predicts the rungs it would be one.

Sources: [Frauenfelder & Moreno 2023 (Springer)](https://link.springer.com/article/10.1007/s10569-023-10135-6),
[arXiv:2301.01803](https://ar5iv.labs.arxiv.org/html/2301.01803),
[Swift & Wiesenfeld 1984 (NAU record)](https://experts.nau.edu/en/publications/suppression-of-period-doubling-in-symmetric-systems/),
[Marti & Leonel 2026, arXiv:2607.05604](https://arxiv.org/html/2607.05604),
[Lamb, k-symmetry papers (Imperial)](https://www.ma.imperial.ac.uk/~jswlamb/papers/kbif.abs),
[Rucklidge & Silber 1998 (M = Gⁿ for spatio-temporal symmetries)](https://eprints.whiterose.ac.uk/986).
