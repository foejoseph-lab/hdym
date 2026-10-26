# Rabbit hole 5 — are Krein collisions a planar-orbit phenomenon? Answered on the x = y = z orbit (22 Oct 2026)

Script `code/rh5_nonplanar.py`, log `results/rh5_nonplanar.txt` (300 energies in [0.02, 3]), data `rh5_nonplanar.npz`.
The "why block 12" half of the hole was answered in `RH7_selection_rule.md` (only block 12's components mix under
the planar orbit's reflection; blocks 13 and 23 are exactly block-diagonal there). The open half was whether, off
the plane, all three blocks couple and collide. The non-planar periodic orbit used is the exact x = y = z orbit of
Yang–Mills mechanics (each component the quartic oscillator s″ = −2s³; z ≠ 0 throughout; T = 15.4 → 4.4 over the
range, a third of the DR orbit's; closure 10⁻¹¹–10⁻¹²). On it the three blocks are related by the cyclic
permutation and each carries a component cross-coupling of relative size 0.19 (the DR orbit's block 12: O(1)).

## Result

**No off-diagonal instability of any type in any block at any of the 300 energies** — every ghost multiplier on
the unit circle to 10⁻⁴, in all three blocks (which are equivalent by the permutation symmetry, and come out
identical). Not for lack of opportunity:

- the two ghost pairs of each block **cross each other with opposite Krein signature** at E ≈ 0.07 (angles 0.394π
  with sign +1 and 0.409π with sign −1, the ordering reversed on the other side) and again at E ≈ 0.37 (0.329π/−1,
  0.366π/+1): the "dangerous" configuration of Krein's criterion, which on the DR orbit always collides — here they
  pass through;
- a pair passes through +1 at E ≈ 0.046 (angle 0.013π), ≈ 0.13 (0.036π) and ≈ 0.70 (0.023π) and through −1 at
  E ≈ 0.13 (0.929π) and 0.03 (0.945π) without folding;
- the diagonal ghost on the same orbit *does* have a parametric tongue (E = 1/8 to 1/6 exactly, peak |μ| = 1.306 at
  E ≈ 0.144; first quoted here as "|μ| = 1.22 at E = 0.13–0.16" without a log — see the Lamé section below), so the orbit is
  not trivially non-resonant;
- at E ≈ 2–2.4 the pairs' near-meeting is same-sign (safe), as expected.

So on this orbit the couplings that would open a Krein band at an opposite-sign crossing, or a parametric band at
±1, vanish — a selection rule of the orbit's larger symmetry group (cyclic permutation of the three components and
s → −s), stricter than the DR orbit's single reflection. (The candidate half-period involution of the DR analysis,
S̃ = diag on derivative parity, is not a symmetry here — residual 2 — so the lemma of RH7 is not what forbids the
folds; the actual symmetry's representation on the block is not worked out.)

## Reading for hole 5

Krein collisions are not "block 12's": they belong to the planar DR orbit, whose reflection symmetry mixes exactly
one block's components and whose cross-coupling is O(1). On the most symmetric non-planar orbit all blocks mix
weakly and none goes unstable at all. Off the plane in general (the chaotic sea), the paper's Lyapunov rates
g₁ ≥ g₂ show the off-diagonal ghost *is* unstable, so instability depends on the orbit through its symmetry; the
right generalisation of "Krein collisions live in one block" is "on a periodic orbit, which blocks collide is set
by which components the orbit's symmetry mixes, and a symmetry large enough forbids the collision altogether".
The chaotic-sea typing (do the sea's growth directions rotate like a quadruplet's or run like a pair's) is not a
Floquet question and stays with hole 12's state-snapshot requirement. Hole closed at the level asked.

## Addendum (26 Oct 2026): the diagonal ghost on x = y = z is Lamé n = 2 — one tongue, exactly

Script `code/rh5_lame.py`, log `results/rh5_lame.txt`. On x = y = z = s(t), s = A cn(√2 A t | k² = ½), E = 3A⁴/2, and
Hess V = 2s²J (J the all-ones matrix): eigenvalues 6s² on (1, 1, 1) and 0, 0. The diagonal ghost ẅ = −(1 + Hess V)w
therefore splits into two anisotropic modes ẅ + w = 0 (always stable) and the isotropic mode ẅ + (1 + 6s²)w = 0
(the same equation as `axion_su2/`'s isotropic w″ + (1 + 6Q²)w), which in u = √2 A t is the **Lamé equation with
n(n+1)k² = 3, k² = ½, so n = 2**, spectral parameter h = 3 + 1/(2A²). Lamé n = 2 is a finite-gap potential (the
finite-gap solutions of inverse scattering): its five band edges are the Lamé-polynomial eigenvalues 1 + k², 1 + 4k²,
4 + k², 2(1 + k²) ∓ 2√(1 − k² + k⁴) = 1.268, 1.5, 3, 4.5, 4.732, with gaps (1.5, 3) and (4.5, 3 + √3) and nothing
above. Since h > 3 for any positive ghost mass, only the second gap is reachable, so the diagonal ghost on this orbit
has **exactly one instability tongue in all of E, E ∈ (1/8, 1/6)**, and is stable at every other energy, in particular
at every lower one. (The unreachable gap ends at h = 3, the ghost-free limit, where the isotropic mode is the orbit's
own tangent.) Checked: closed-form edges against a numerical Lamé discriminant (1.6 × 10⁻¹²); a direct Floquet scan
of the full 3 × 3 diagonal ghost along the integrated orbit over the RH5 range E ∈ [0.02, 3] finds unstable points only
inside (1/8, 1/6); the tongue's edges by bisection agree with 1/8 and 1/6 to 3 × 10⁻¹³; peak |μ| = 1.306 at E ≈ 0.144.

Reading: the finite number of tongues is a property of the background, not of the symmetry: an elliptic (integrable)
background gives a finite-gap Hill operator and finitely many tongues; the DR orbit is not elliptic and its tongues
form the open-ended ladder of hole 2. This complements the symmetry argument above (which couplings vanish) with a
second one (why so few resonances exist at all). The off-diagonal blocks on this orbit are coupled two-component
systems; whether they are finite-gap matrix potentials is not checked.
