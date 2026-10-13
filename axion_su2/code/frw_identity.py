"""frw_identity.py -- the Lee-Wick term on spatially flat FRW, in conformal time (sympy check).

Claim.  On ds^2 = a(eta)^2 (deta^2 - dx^2), for any antisymmetric F^{mu nu} (and so for the
gauge-covariant divergence too, since the gauge part carries no metric),
    J^nu = nabla_mu F^{mu nu} = (1/sqrt(-g)) d_mu( sqrt(-g) F^{mu nu} ),
and with F_{mu nu} the same tensor as in flat space (Yang-Mills is conformal: F^{mu nu} =
a^-4 F_flat^{mu nu}) this gives J^nu = a^-4 J_flat^nu, hence
    sqrt(-g) g_{nu rho} J^nu J^rho = a^-2 (J_flat . J_flat)   and   sqrt(-g) F^2 = F_flat^2.
So the Lee-Wick Lagrangian  -F^2/4 + J^2 / (2 M^2)  in conformal time is the flat-space one with
M replaced by M a(eta).  In the auxiliary-field form used by ym/hdym_core.py the mass enters only
algebraically (as -(M^2/2)|w|^2 after the field redefinition), so a time-dependent M(eta) adds no
derivative terms: the linear diagonal ghost obeys  w'' = -(M^2 a(eta)^2 + Hess V(q)) w  on the
unchanged conformal-time background q'' = -grad V(q).

This script checks the two metric identities by brute force (Christoffel symbols) for a
homogeneous antisymmetric F^{mu nu}(eta) with all six components free, which is the case the
homogeneous ansatz needs.  The auxiliary-field step is checked in frw_ghost.py's docstring by
derivation, not here.
Usage: python axion_su2/code/frw_identity.py > axion_su2/results/frw_identity.txt
Added 2026-10-12.
"""
import sympy as sp

eta = sp.Symbol('eta')
X = sp.symbols('x1 x2 x3')
co = [eta, *X]
a = sp.Function('a')(eta)
g = sp.diag(a ** 2, -a ** 2, -a ** 2, -a ** 2)
gi = g.inv()
sqrtg = a ** 4

# Christoffel symbols Gamma^l_{mn}
Gam = [[[sum(gi[l, k] * (sp.diff(g[k, m], co[n]) + sp.diff(g[k, n], co[m]) - sp.diff(g[m, n], co[k]))
              for k in range(4)) / 2 for n in range(4)] for m in range(4)] for l in range(4)]

# generic homogeneous antisymmetric F^{mu nu}(eta) (upper indices), from flat-space components
f = {}
for m in range(4):
    for n in range(m + 1, 4):
        f[(m, n)] = sp.Function(f'G{m}{n}')(eta)
Fflat = sp.zeros(4, 4)
for (m, n), v in f.items():
    Fflat[m, n], Fflat[n, m] = v, -v
Fup = Fflat / a ** 4            # conformal: F_{mu nu} unchanged, raising with g^-1 gives a^-4

# covariant divergence by Christoffels:  nabla_m F^{m n} = d_m F^{mn} + G^m_{mk} F^{kn} + G^n_{mk} F^{mk}
Jcov = [sp.simplify(sum(sp.diff(Fup[m, n], co[m]) for m in range(4))
                    + sum(Gam[m][m][k] * Fup[k, n] for m in range(4) for k in range(4))
                    + sum(Gam[n][m][k] * Fup[m, k] for m in range(4) for k in range(4))) for n in range(4)]
Jflat = [sum(sp.diff(Fflat[m, n], co[m]) for m in range(4)) for n in range(4)]

ok1 = all(sp.simplify(Jcov[n] - Jflat[n] / a ** 4) == 0 for n in range(4))
eta_flat = sp.diag(1, -1, -1, -1)
J2curved = sp.simplify(sqrtg * sum(g[m, n] * Jcov[m] * Jcov[n] for m in range(4) for n in range(4)))
J2flat = sum(eta_flat[m, n] * Jflat[m] * Jflat[n] for m in range(4) for n in range(4))
ok2 = sp.simplify(J2curved - J2flat / a ** 2) == 0
Fdown = g * Fup * g
F2curved = sp.simplify(sqrtg * sum(Fdown[m, n] * Fup[m, n] for m in range(4) for n in range(4)))
F2flat = sum((eta_flat * Fflat * eta_flat)[m, n] * Fflat[m, n] for m in range(4) for n in range(4))
ok3 = sp.simplify(F2curved - F2flat) == 0
print("nabla_mu F^{mu nu} = a^-4 (flat divergence), all four nu:", ok1)
print("sqrt(-g) J.J = a^-2 (flat J.J):", ok2)
print("sqrt(-g) F.F = flat F.F (conformal invariance of the Yang-Mills term):", ok3)
print("=> -F^2/4 + J^2/(2 M^2) on FRW (conformal time) = flat-space Lagrangian with M -> M a(eta):",
      ok1 and ok2 and ok3)
