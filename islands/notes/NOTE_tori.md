# Invariant tori, cantori and the ghost: what the words mean here and what is expected (29 Oct 2026)

A record, so the question does not have to be reopened from memory. Mostly standard theory, stated for this project;
what this repository measured is marked and pointed to. Literature is from memory, to verify.

## The words

- **Torus** is meant topologically: an n-dimensional surface in a 2n-dimensional phase space with the topology of
  n circles multiplied together, deformed and not a round doughnut. In the planar diagonal system (2 degrees of freedom,
  4-d phase space, 3-d energy shell) a torus is a 2-d surface. On a Poincaré section it shows as a closed curve.
- **Why "angular" keeps coming up.** Near a regular orbit the motion can be written in action-angle variables. The
  actions J_i = ∮ p dq are conserved and label the torus; the angles θ_i wind around it at fixed frequencies.
  - The angles are angles on the torus, and the actions are generalised momenta.
  - Neither is the spatial angular momentum L_i = ε_ijk A_j^a P_k^a of the gauge field.
  - The colour-diagonal sector always has L = 0, so everything in this note is at L = 0.
  - L enters elsewhere (hdym-su3, the generic-L ensembles) as a separate conserved quantity.
- **KAM tori.** For an orbit that is elliptic (all Floquet multipliers on the unit circle) and has twist (frequency
  changing with amplitude), the KAM theorem (Kolmogorov, Arnold, Moser) says most tori of the nearby integrable
  approximation survive: those with sufficiently irrational frequency ratios. Together they fill a set of positive
  measure, the island. The certificate in this folder is of that kind: elliptic, with nonzero twist, then Moser's
  theorem (`islands/README.md`).
- **Cantori.** When a torus breaks it leaves a Cantor-set remnant, a "cantorus" (Aubry–Mather theory; Percival). That is
  an invariant set full of gaps, which still slows transport through it. This is the standard explanation of
  stickiness: chaotic orbits linger near an island's edge, and their exit times have an algebraic tail rather than an
  exponential one (MacKay, Meiss & Percival 1984; Karney 1983; Chirikov & Shepelyansky 1984). Satellite island chains
  near the edge (island-around-island hierarchy) add to it.

## What dimension does

- **2 degrees of freedom:** a 2-d torus divides the 3-d energy shell, so KAM tori are walls. An island is sealed, and
  cantori are partial walls.
- **3 or more degrees of freedom:** a torus no longer divides the energy shell, so orbits can drift around tori
  (Arnold diffusion). Islands leak, slowly. The 3-dof diagonal DR island does (2.7 % at r = 2 × 10⁻⁴; no committed
  log, `code/hdym_island.py` docstring). Cantori have no clean analogue, and stickiness comes from partial barriers.
- **Infinitely many (the field theory):** KAM results exist for model nonlinear PDEs (Kuksin, Wayne, Craig, Pöschel,
  Bourgain): finite-dimensional tori survive. None is known for Yang–Mills, and none is expected to confine.

## What the ghost does (expected [est]; measured where marked)

- **Negative energy does not forbid tori.** KAM is local and needs ellipticity and twist, not a positive Hamiltonian.
- **The ghost does change the linearisation.** The orbit is a candidate for full-dimensional tori only if every
  direction is elliptic, ghost directions included. On the DR orbit:
  - at E = 0.3 the off-diagonal block 23 is parametrically unstable (rate 0.0125; measured, `islands/README.md`);
  - at E = 1 the diagonal ghost is unstable (rate 0.233; smoke run, no committed log).

  With a hyperbolic direction, no full-dimensional tori surround the orbit. What generically survives is
  lower-dimensional "whiskered" tori, elliptic in some directions and hyperbolic in the rest (Graff, Zehnder, Eliasson).
  They are remnants, but they do not confine. That fits rabbit hole 9 (measured): a 10⁻⁸ off-diagonal seed empties the
  E = 0.3 island's core.
- **Long-time stability is weaker too.** With modes of both energy signs, Nekhoroshev-type bounds, which rely on
  convexity-like conditions, are not guaranteed. Diffusion can then be fast, not exponentially slow.
- **Summary:** some remnant of the torus structure is expected with any ghost, in any dimension. Confinement is expected
  only where the full linearisation is elliptic and there are just two degrees of freedom.

## What this repository has measured

- The planar ghost-free island and its sticky edge: an algebraic exit-time tail over more than two decades, from
  orbits hugging the outermost KAM boundary (`code/cantori_residence.py`, `rabbit_holes/ym/notes/RH9_island_full.md`,
  "Cantori: the direct test", 29 Oct 2026). Positive by the criterion stated before the run. Not separated: cantori
  from satellite-island trapping.
- Not measured: stickiness at 3 dof, or with the ghost.
