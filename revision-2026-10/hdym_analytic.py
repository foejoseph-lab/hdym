"""hdym_analytic.py -- symbolic derivation of the channel-tail constants for a diagonal
potential with N = n + 1 degrees of freedom,

    V = 1/2 sum_{i<j} q_i^2 q_j^2   (+ optional regulator 1/2 mu^2 |q|^2),

in the single-mode (diagonal-ghost) approximation of the paper.  Every number quoted in
main.tex Sec. 5 and in NOTE_2d.md Sec. 3-5 is reproduced here from the two inputs
  (i) the microcanonical entry flux  nu0 = 2 (2 pi)^n / Z_N   (per channel, both ends),
  (ii) adiabatic channel passage     1/2 xdot^2 + J |x| = 1,  q_a = sqrt(2 J_a / |x|) cos phi_a,
and checked against mpmath.  Run:  python hdym_analytic.py

Outputs
  1. passage integral  I_n = int_{simplex} (sum a_i^2) / sqrt(2 (1 - sum a)) d^n a = sqrt(2 pi) n / Gamma(n + 5/2)
  2. resonant tail      u^{n+4} s_res(u) = 2 pi 2^{n+2} I_n nu0
  3. ghost law          lambda_diag = (pi/16) I_n nu0 E^{(n+7)/4} = C_n E^{(n+7)/4} / Z_N
  4. second-order response of the along-channel coordinate and the detuned/resonant split:
     total tail / resonant tail = (9n - 4) / 4,  resonant share 4 / (9n - 4)   (5/4 and 7/2)
  5. regulated 2-d density of states  Z_2(mu) = 8 sqrt2 pi sqrt(1+b^2) [K(m) - E(m)],  b^2 = mu^4/2,
     m = 1/(1+b^2),  and its small-mu asymptotics  Z_2 ~ 16 sqrt2 pi ln(1/mu) + 8 sqrt2 pi (ln(4 sqrt2) - 1)
"""
import sympy as sp
import mpmath as mp

n = sp.symbols('n', positive=True, integer=True)
u, E, Z, D, J, Om, Y, t = sp.symbols('u E Z D J Omega Y t', positive=True)

print("=" * 78)
print("1. Passage integral over the n-simplex (Dirichlet integral)")
# int_{sum a < 1} a_1^2 (1 - sum a)^{-1/2} d^n a = Gamma(3) Gamma(1)^{n-1} Gamma(1/2) / Gamma(3 + (n-1) + 1/2)
one_term = sp.gamma(3) * sp.gamma(sp.Rational(1, 2)) / sp.gamma(n + sp.Rational(5, 2))
I_n = sp.simplify(n * one_term / sp.sqrt(2))
print("   I_n =", I_n)
for nn, ref in ((1, sp.Rational(16, 15) / sp.sqrt(2)), (2, sp.Rational(64, 105) / sp.sqrt(2))):
    val = sp.simplify(I_n.subs(n, nn))
    # brute-force check by direct integration
    if nn == 1:
        a = sp.symbols('a', positive=True)
        direct = sp.integrate(a ** 2 / sp.sqrt(2 * (1 - a)), (a, 0, 1))
    else:
        a, b = sp.symbols('a b', positive=True)
        direct = sp.integrate(sp.integrate((a ** 2 + b ** 2) / sp.sqrt(2 * (1 - a - b)), (b, 0, 1 - a)), (a, 0, 1))
    print(f"   n={nn}: I_n = {val} = {float(val):.6f}; direct integration {sp.nsimplify(sp.simplify(direct))}; "
          f"paper/note value {ref}: match = {sp.simplify(val - ref) == 0 and sp.simplify(direct - ref) == 0}")

print("=" * 78)
print("2. Resonant tail of the diagonal entry K_xx = sum_{b != x} q_b^2 in the x-channel")
# chirp amplitude of q_b^2 is J_b / D at frequency 2 D; chirp spectrum pi <(a^2/2) delta(u - 2|x|)>
# time at depth D per visit: 2 dD / sqrt(2 (1 - J D)); entries at rate nu0 d^n J (uniform)
nu0 = sp.symbols('nu0', positive=True)
# <sum J_b^2 delta(u - 2|x|)> = nu0 int d^n J (sum J^2) / sqrt(2 (1 - J D)),  J = 2 a / u  ->  nu0 (2/u)^{n+2} I_n
s_res = 2 * sp.pi / u ** 2 * nu0 * (2 / u) ** (n + 2) * I_n
A_res = sp.simplify(s_res * u ** (n + 4))
print("   u^{n+4} s_res(u) = A_res =", A_res)
print("   n=2:", sp.simplify(A_res.subs(n, 2)), "= 32 pi I nu0 :", sp.simplify(A_res.subs(n, 2) - 32 * sp.pi * I_n.subs(n, 2) * nu0) == 0)
print("   n=1:", sp.simplify(A_res.subs(n, 1)), "= 256 pi/(15 sqrt2) nu0 :", sp.simplify(A_res.subs(n, 1) - 256 * sp.pi / (15 * sp.sqrt(2)) * nu0) == 0)

print("=" * 78)
print("3. Ghost law: lambda = (1/8) E^{3/4} s(2 E^{-1/4}) with the resonant piece only")
lam = sp.Rational(1, 8) * E ** sp.Rational(3, 4) * s_res.subs(u, 2 * E ** sp.Rational(-1, 4))
lam = sp.powsimp(sp.simplify(lam), force=True)
print("   lambda_diag =", lam)
coef = sp.simplify(lam / E ** ((n + 7) / 4) / nu0)
print("   = (pi/16) I_n nu0 E^{(n+7)/4}:", sp.simplify(coef - sp.pi / 16 * I_n) == 0)
# flux: (2 pi)^n per channel end, two ends
nu0_expr = 2 * (2 * sp.pi) ** n / Z
C_n = sp.simplify(coef * nu0_expr * Z)        # coef is per unit nu0
print("   C_n (so that lambda = C_n E^{(n+7)/4} / Z_N) =", C_n)
for nn, ref in ((2, 32 * sp.pi ** 3 / (105 * sp.sqrt(2))), (1, 4 * sp.pi ** 2 / (15 * sp.sqrt(2)))):
    v = sp.simplify(C_n.subs(n, nn))
    print(f"   n={nn}: C = {v} = {float(v):.5f}  (reference {ref} = {float(ref):.5f}, match {sp.simplify(v - ref) == 0});"
          f"  exponent (n+7)/4 = {sp.Rational(nn + 7, 4)}")
print(f"   3-d with Z = 557.1: {float(C_n.subs(n, 2)) / 557.1:.5f} E^(9/4)   [paper: 0.0120]")
print(f"   2-d with Z2(0.1) = 189.73: {float(C_n.subs(n, 1)) / 189.73:.5e} E^2   [note: 9.809e-3]")

print("=" * 78)
print("4. Second-order response and the detuned/resonant split")
# In channel c (along-coordinate Y, |Y| = D), transverse oscillators q_a = sqrt(2 J_a/Om) cos(phi_a), Om = |Y|.
# Equation of motion of Y:  Ydd = -Y sum_a q_a^2 = -sgn(Y) sum_a J_a (1 + cos 2 phi_a).
# Oscillating force F = -sgn(Y) J_a cos(2 phi_a) at frequency 2 Om  ->  response dY = -F/(2 Om)^2:
phi = sp.symbols('phi', real=True)
F = -J * sp.cos(2 * phi)                                   # per oscillator, sgn(Y) = +1
dY = sp.simplify(-F / (2 * Om) ** 2)                      # dY'' = F with phi = Om t
# check: d^2/dt^2 [dY(phi = Om t)] == F
chk = sp.simplify(sp.diff(dY.subs(phi, Om * t), t, 2) - F.subs(phi, Om * t))
print("   dY = ", dY, "  (solves dY'' = F:", chk == 0, ")")
Y_sq_osc = sp.simplify(2 * Om * dY)                        # oscillating part of Y^2 = (Y + dY)^2
print("   oscillating part of Y^2 from oscillator a:  2 Y dY =", Y_sq_osc, "   [amplitude J_a / (2 D)]")
direct_chirp = J / Om                                       # amplitude of q_b^2 itself, for b transverse in channel c
print("   direct chirp amplitude of q_b^2:", direct_chirp, "              [J_b / D]")
# K_xx in channel c != x:  Y^2 + sum_{b not in {x,c}} q_b^2.  Coefficient of cos 2 phi_b:
coef_b_other = sp.Rational(1, 2) + 1       # b not in {x, c}: appears in Y^2 response AND directly
coef_b_x = sp.Rational(1, 2)               # b = x: only through the Y^2 response
# powers (phases independent across oscillators), all <J_b^2> equal:
detuned_per_channel = (n - 1) * coef_b_other ** 2 + coef_b_x ** 2
resonant = n                               # x-channel: sum_{a != x} J_a^2, n terms, amplitude 1
total_over_res = sp.simplify((resonant + n * detuned_per_channel) / resonant)
print("   coefficient of cos 2phi_b in K_xx, channel c != x:  b not in {x,c}:", coef_b_other, "; b = x:", coef_b_x)
print("   (total tail)/(resonant tail) =", total_over_res, " ;  resonant share =", sp.simplify(1 / total_over_res))
for nn in (1, 2, 3):
    print(f"   n={nn}: total/res = {total_over_res.subs(n, nn)},  share = {sp.simplify(1 / total_over_res).subs(n, nn)}"
          f" = {float(1 / total_over_res.subs(n, nn)):.4f}")
print("   3-d: A = 7/2 A_res = 7/2 * 32 pi I nu0 with nu0 = 8 pi^2 / 557.1 :",
      f"{float(sp.Rational(7, 2) * A_res.subs(n, 2).subs(nu0, 8 * sp.pi ** 2 / 557.1)):.2f}   [measured 21.6 +- 0.3]")
print("   2-d: total/res = 5/4   [measured 1.25 +- 0.02 at mu1 = 0.1 and 0.03]")
# decoupled z-component driver in 2-d: K_zz = x^2 + y^2 in the x-channel = Y^2 (response, 1/2) + y^2 (direct, 1), coherent
print("   2-d K_zz = x^2 + y^2 in either channel: coefficient 1/2 + 1 = 3/2, power 9/4 per channel, total 9/2 A_res   [measured 4.3-4.4]")

print("=" * 78)
print("5. Regulated 2-d density of states Z2(mu) = 2 pi * area{ x^2 y^2 / 2 + mu^2 (x^2 + y^2)/2 < 1 }")
mu, s, b = sp.symbols('mu s b', positive=True)
# area = 4 int_0^{sqrt2/mu} sqrt((2 - mu^2 x^2)/(x^2 + mu^2)) dx ;  x = sqrt2 s / mu  ->  8 int_0^1 sqrt(1 - s^2)/sqrt(2 s^2 + mu^4) ds
# so Z2 = 16 pi int_0^1 sqrt(1-s^2)/sqrt(2 s^2 + mu^4) ds = 8 sqrt2 pi G(b), b^2 = mu^4/2,
#    G(b) = int_0^1 sqrt((1 - s^2)/(s^2 + b^2)) ds
G_int = sp.Integral(sp.sqrt((1 - s ** 2) / (s ** 2 + b ** 2)), (s, 0, 1))
# closed form by s = cos(theta):  G(b) = (1/sqrt(1+b^2)) int_0^{pi/2} sin^2 th / sqrt(1 - m sin^2 th) dth
#                                       = sqrt(1+b^2) [K(m) - E(m)],   m = 1/(1+b^2)   (sympy/mpmath parameter convention)
def G_num(bv):
    return mp.quad(lambda x: mp.sqrt((1 - x * x) / (x * x + bv * bv)), [0, 1])
def G_closed(bv):
    m = 1 / (1 + bv * bv)
    return mp.sqrt(1 + bv * bv) * (mp.ellipk(m) - mp.ellipe(m))
ok = True
for bv in (mp.mpf('1.3'), mp.mpf('0.5'), mp.mpf('0.1'), mp.mpf('0.00707'), mp.mpf('0.000636')):
    gn, gc = G_num(bv), G_closed(bv)
    ok &= abs(gn - gc) < 1e-9
    print(f"   b={float(bv):.3g}: G numeric {float(gn):.10f}  closed form {float(gc):.10f}")
print("   G(b) = sqrt(1+b^2) [K(m) - E(m)],  m = 1/(1+b^2):", ok)
print("   Z2(mu) = 8 sqrt2 pi sqrt(1+b^2) [K(m) - E(m)],  b^2 = mu^4/2,  m = 1/(1+b^2)")
for muv in (0.1, 0.03):
    bv = mp.mpf(muv) ** 2 / mp.sqrt(2)
    print(f"   mu1={muv}: Z2 = {float(8 * mp.sqrt(2) * mp.pi * G_closed(bv)):.4f}   [quadrature in hdym_2d.py: {189.73 if muv == 0.1 else 275.32}]")
# small-b asymptotics: E(m->1) = 1, K(m) = ln(4/sqrt(1-m)) + ..., sqrt(1-m) = b/sqrt(1+b^2)
# G(b) = int_0^1 sqrt(1-s^2)/sqrt(s^2+b^2) ds ;  split at s = 1: G = int_0^1 [sqrt(1-s^2) - 1]/s ds + int_0^1 ds/sqrt(s^2+b^2) + O(b^2 ln b)
part1 = sp.log(2) - 1          # int_0^1 (sqrt(1-s^2)-1)/s ds
part2 = sp.asinh(1 / b)
Gasym = sp.simplify(part1 + part2)
print("   small-b:  G(b) = asinh(1/b) + int_0^1 (sqrt(1-s^2)-1)/s ds + O(b^2 ln b) =", Gasym)
print("            = ln(2/b) + ln 2 - 1 + O(b^2 ln b) = ln(4/b) - 1;  with b = mu^2/sqrt2:")
Z2asym = sp.expand(8 * sp.sqrt(2) * sp.pi * (sp.log(4 * sp.sqrt(2) / mu ** 2) - 1))
print("   Z2(mu) ~ 8 sqrt2 pi [ 2 ln(1/mu) + ln(4 sqrt2) - 1 ] = 16 sqrt2 pi ln(1/mu) + ", float(8 * sp.sqrt(2) * sp.pi * (sp.log(4 * sp.sqrt(2)) - 1)))
for muv in (0.1, 0.03):
    print(f"   mu1={muv}: asymptotic Z2 = {float(Z2asym.subs(mu, muv)):.2f}")
print("   pure x^2 y^2 (mu -> 0): Z2 -> infinity logarithmically; nu0 = 4 pi / Z2 -> 0 (NOTE_2d.md Sec. 1)")
print("=" * 78)
print("Summary:  lambda_diag = C_n E^{(n+7)/4} / Z_N,   C_n = (pi/8)(2 pi)^n sqrt(2 pi) n / Gamma(n + 5/2)")
print("          resonant share of the tail = 4/(9n - 4);   n = 1: 4/5,  n = 2: 2/7")

print("=" * 78)
print("6. Off-diagonal ghost block with c_j = 0 (hdym_rung1): the xy-entry resonance, lambda ~ E^{5/4}")
# m = z^2 I - xy sigma_x has fixed eigenvectors (1, +-1): w_+- obey  w'' = -(1 + z^2 -+ xy) w.
# In the x-channel, xy = x * sqrt(2 J_y / |x|) cos(phi_y) = sqrt(2 J_y D) cos(phi_y): chirp at frequency D
# (not 2D), amplitude^2 / 2 = J_y D.  Resonance of the unit oscillator at driving frequency 2 <-> D = 2.
# Both channels resonate (the ghost frequency stays ~1 in both for c_j = 0).
# <J_y delta(u - |x|)>_time = nu0 int_{J<1/u} d^2J J_y * 2 / sqrt(2 (1 - J u))  =  2 nu0 u^-3 I1',
a, b = sp.symbols('a b', positive=True)
I1p = sp.integrate(sp.integrate(a / sp.sqrt(2 * (1 - a - b)), (b, 0, 1 - a)), (a, 0, 1))
print("   I1' = int_{a+b<1} a / sqrt(2(1-a-b)) =", sp.nsimplify(sp.simplify(I1p)))
s_xy_one_channel = sp.pi * u * 2 * nu0 * I1p / u ** 3          # pi <(a^2/2) delta(u - Omega)>, a^2/2 = J_y u
s_xy = 2 * s_xy_one_channel                                     # x- and y-channels
A_xy = sp.simplify(s_xy * u ** 2)
print("   u^2 s_xy(u) =", A_xy, "=", float(A_xy.subs(nu0, 8 * sp.pi ** 2 / 557.1)), "at Z = 557.1   [paper: 'u^-2 tail']")
lam_pm = sp.simplify(sp.Rational(1, 8) * E ** sp.Rational(3, 4) * s_xy.subs(u, 2 * E ** sp.Rational(-1, 4)))
C_pm = sp.simplify(lam_pm / E ** sp.Rational(5, 4) / nu0 * (8 * sp.pi ** 2))     # per 1/Z
print("   lambda_+- =", lam_pm, "=", sp.simplify(C_pm), "* E^{5/4} / Z")
val = float(C_pm) / 557.1
print(f"   = {val:.4f} E^(5/4)   [c_j = 0 GPU: E = 0.3 -> {val * 0.3 ** 1.25:.3e} predicted, 5.07e-3 measured (mean of g1, g2); "
      f"E = 0.1 -> {val * 0.1 ** 1.25:.3e} predicted]")
print("   (The z^2 part of K_+- adds the diagonal-type E^{9/4} term with half the paper's prefactor; negligible here.)")
