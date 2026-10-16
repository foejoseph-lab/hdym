"""probe_core.py -- the one right-hand side for the homogeneous Lee-Wick SU(2) system with the
representative probes (plan item 5, 16 Oct: real particles coupled to the ghost), shared by
axion_su2/ and dilaton_su2/.  Everything is an initial-value problem in phase space; no Green's
function, no quantisation.  Derived and checked against sympy in probe_derive.py (random points,
every switch combination, energy conservation); change this file and rerun that one.

System (units g = M = 1 by default; g and M are parameters so the checks can keep them):

  gauge field, diagonal ansatz A^a_i = delta^a_i q_i, in the (q, w) auxiliary form of ym/
      (hdym_derive.py CHECK 4; w = (qddot + grad V)/M^2 on shell for the ghost):
        L' = Theta(t) (1/2|qdot|^2 - V) + qdot.wdot - w.grad V + M^2/2 |w|^2 + kappa(t) x y z
             - (g^2/4)|q|^2 f^2 + fdot^2 - ms^2 f^2 + L_D
  drives (either, both, or none; prescribed functions of t):
        axion    kappa(t) = (g/f_a) phidot:  the potential -kappa x y z      (axion_derive.py)
        dilaton  Theta(t) = e^{phi/Lambda}:  on the Yang-Mills part only, column A of dilaton_derive.py
  probes:
        boson    a fundamental scalar doublet Phi; the Gauss law forces Phi = f(t) Phi_0 with f real
                 (probe_derive.py CHECK 1), so the probe is one real field f with
                 L_f = fdot^2 - (ms^2 + (g^2/4)(x^2 + y^2 + z^2)) f^2  (f = |Phi|, hence the 1 not 1/2);
                 it couples to the full q = p - w, ghost included, and pushes back with -(g^2/2) q f^2.
        fermion  a classical Dirac doublet psi (8 complex components) in the diagonal-SO(3) singlet
                 sector, which is invariant on every diagonal background (CHECK 3) and is a
                 two-level system zeta in C^2:  i zeta' = h zeta,  h = [[mD, c g Q], [c g Q, -mD]],
                 c = 3/2 (the sign is a basis convention), Q = (x + y + z)/3; it feeds g j into
                 each gauge equation with j = -(c/3) zeta^dagger sigma_x zeta (= dL_D/dq_i).  It sees only the
                 trace x + y + z, so it does not couple to the anisotropic (quintet) perturbations at
                 linear order (CHECK 3).
  ghost or twin:
        ghost (sigma = +1) is the Lee-Wick theory above.  twin (sigma = -1) flips the heavy sector to
        a healthy oscillator in the (p, w) variables, L_twin = L' + |wdot|^2 - M^2|w|^2, i.e.
        1/2|pdot|^2 + 1/2|wdot|^2 - V(q) - w.grad V(q) - M^2/2 |w|^2 with p = q + w: a mathematics /
        code control, not a model.  It is defined only with Theta = 1 (with the dilaton on, any sign
        flip of the w sector makes the kinetic matrix degenerate at some Theta; the dilaton's twin
        control is done on the periodic-coefficient equations instead).

State vector s (18 reals):  q(3), qdot(3), w(3), wdot(3), f, fdot, Re zeta(2), Im zeta(2).
Probes that are off keep their entries at zero and contribute nothing.
"""
import numpy as np

C_FERMION = 1.5      # the singlet coupling constant, from probe_derive.py CHECK 3


class Params:
    def __init__(self, g=1.0, M=1.0, ms=1.0, mD=1.0, sigma=+1, boson=False, fermion=False,
                 kappa=None, Theta=None, Thetadot=None):
        self.g, self.M, self.ms, self.mD = g, M, ms, mD
        self.sigma = sigma                      # +1 ghost, -1 twin
        self.boson, self.fermion = boson, fermion
        self.kappa = kappa if kappa is not None else (lambda t: 0.0)
        self.Theta = Theta if Theta is not None else (lambda t: 1.0)
        self.Thetadot = Thetadot if Thetadot is not None else (lambda t: 0.0)
        if sigma == -1 and (Theta is not None):
            raise ValueError("the twin is defined only with Theta = 1 (see docstring)")


def V(q, g):
    x, y, z = q
    return g**2 / 2 * (x * x * y * y + y * y * z * z + z * z * x * x)


def gradV(q, g):
    x, y, z = q
    return g**2 * np.array([x * (y * y + z * z), y * (z * z + x * x), z * (x * x + y * y)])


def hessV(q, g):
    x, y, z = q
    return g**2 * np.array([[y * y + z * z, 2 * x * y, 2 * x * z],
                            [2 * x * y, z * z + x * x, 2 * y * z],
                            [2 * x * z, 2 * y * z, x * x + y * y]])


def grad_xyz(q):
    x, y, z = q
    return np.array([y * z, z * x, x * y])


def fermion_h(Q, P):
    return np.array([[P.mD, C_FERMION * P.g * Q], [C_FERMION * P.g * Q, -P.mD]])


def fermion_current(zeta):
    """j = -(c/3) zeta^dagger sigma_x zeta, the same for all three gauge directions on the singlet."""
    return -(C_FERMION / 3.0) * 2.0 * np.real(np.conj(zeta[0]) * zeta[1])


def unpack(s):
    q, v, w, u = s[0:3], s[3:6], s[6:9], s[9:12]
    f, fd = s[12], s[13]
    zeta = s[14:16] + 1j * s[16:18]
    return q, v, w, u, f, fd, zeta


def rhs(t, s, P):
    q, v, w, u, f, fd, zeta = unpack(s)
    g, M = P.g, P.M
    Th, Thd, ka = P.Theta(t), P.Thetadot(t), P.kappa(t)
    gV, H = gradV(q, g), hessV(q, g)
    # the extra force on q from drives and probes (enters the q-equation, not the w-equation)
    F = ka * grad_xyz(q)
    if P.boson:
        F = F - (g**2 / 2) * q * f * f
    if P.fermion:
        F = F + g * fermion_current(zeta)
    if P.sigma == +1:
        qdd = M**2 * w - gV
        wdd = -H @ w - Th * M**2 * w - Thd * v + F
    else:
        # L_twin = L' + |wdot|^2 - M^2 |w|^2, Theta = 1
        wdd = H @ w - M**2 * w - F
        qdd = -gV - M**2 * w - 2.0 * wdd
    out = np.zeros_like(s)
    out[0:3], out[3:6], out[6:9], out[9:12] = v, qdd, u, wdd
    if P.boson:
        out[12] = fd
        out[13] = -(P.ms**2 + (g**2 / 4) * (q @ q)) * f
    if P.fermion:
        Q = (q[0] + q[1] + q[2]) / 3.0
        zd = -1j * (fermion_h(Q, P) @ zeta)
        out[14:16], out[16:18] = zd.real, zd.imag
    return out


def energy(t, s, P):
    """Total energy of L' (Legendre transform), the probes included.  Conserved when kappa and Theta
    are constant.  The fermion's part is zeta^dagger h zeta (indefinite, as for any classical Dirac field)."""
    q, v, w, u, f, fd, zeta = unpack(s)
    g, M = P.g, P.M
    Th, ka = P.Theta(t), P.kappa(t)
    if P.sigma == +1:
        # L' = Theta (v^2/2 - V) + v.u - w.gradV + M^2 w^2/2 + ...: E = Theta (v^2/2 + V) + v.u + w.gradV - M^2 w^2/2
        E = Th * ((v @ v) / 2 + V(q, g)) + v @ u + w @ gradV(q, g) - M**2 / 2 * (w @ w)
    else:
        E = 0.5 * (v @ v) + v @ u + (u @ u) + V(q, g) + w @ gradV(q, g) + M**2 / 2 * (w @ w)
    E -= ka * q[0] * q[1] * q[2]
    if P.boson:
        E += fd * fd + (P.ms**2 + (g**2 / 4) * (q @ q)) * f * f
    if P.fermion:
        Q = (q[0] + q[1] + q[2]) / 3.0
        E += np.real(np.conj(zeta) @ (fermion_h(Q, P) @ zeta))
    return E


def fermion_adiabatic_number(Q, zeta, P):
    """Diagnostic only: the weight of zeta on the instantaneous negative-energy eigenvector of h
    (a positive-energy state that leaks there has 'produced a pair' in the first-quantised picture)."""
    h = fermion_h(Q, P)
    ev, U = np.linalg.eigh(h)
    neg = U[:, 0]
    return abs(np.conj(neg) @ zeta)**2 / max(np.real(np.conj(zeta) @ zeta), 1e-300)
