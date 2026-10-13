"""frw_ghost.py -- the isotropic diagonal ghost crossing its Lame band in a radiation-dominated FRW.

Derivation (frw_identity.py checks the metric part).  In conformal time eta the Lee-Wick Lagrangian
on spatially flat FRW is the flat one with M -> M a(eta).  Write the higher-derivative term with
the auxiliary field of ym/hdym_core.py:  -J^2/(2 M^2) = stat_w [ w.J + (M^2/2)|w|^2 ],  J = q'' + grad V,
integrate w.q'' by parts (no M involved), and vary:  q'' = -grad V - M^2 w,  w'' = -(M^2 + Hess V) w.
M enters only algebraically, so M -> M a(eta) needs no extra terms.  On the ghost-free background
(w = 0) q(eta) is the flat-space orbit in conformal time (Yang-Mills is conformal), and the linear
diagonal ghost obeys
    w'' = -( M^2 a(eta)^2 + Hess V(q(eta)) ) w.
Units M = 1 (physical mass), comoving isotropic amplitude A = 1 (any A rescales to this).

On x = y = z = s, s = cn(sqrt2 eta | 1/2): isotropic ghost mode w'' = -(a^2 + 6 s^2) w; frozen-a Lame
parameter h = 3 + a^2/2, unstable for h in (4.5, 4.732), i.e. a in (sqrt3, (2 sqrt3)^(1/2)) =
(1.7321, 1.8612): the physical energy E_phys = (3/2)/a^4 (units M^4/g^2) falls through (1/6, 1/8).
Anisotropic ghost modes w'' = -a^2 w: no driving.

Radiation era: a = eta / eta_s.  The physical Hubble rate at the band, H/M = a'/a^2 = 1/(eta_s a^2),
is about 0.31/eta_s.  For each eta_s, integrate the isotropic mode from a = 1.4 to a = 2.3 and report
ln sigma_max of the transfer matrix in WKB-normalised variables (sqrt(Om) w, w'/sqrt(Om)),
Om^2 = a^2 + 3 (a^2 + the period average of 6 s^2): the e-folds gained in crossing.  The same with
a frozen at the band centre gives the Floquet rate mu for the comparison  e-folds ~ mu * (eta spent
in the band).

Float only.  Usage: python axion_su2/code/frw_ghost.py > axion_su2/results/frw_ghost.txt
Added 2026-10-12.
"""
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipk

RT = dict(method='DOP853', rtol=1e-10, atol=1e-12)
A_LO, A_HI = np.sqrt(3.0), np.sqrt(2.0 * np.sqrt(3.0))
T_ORB = 4.0 * ellipk(0.5) / np.sqrt(2.0)          # period of s at A = 1


def crossing(eta_s, a0=1.4, a1=2.3):
    e0, e1 = a0 * eta_s, a1 * eta_s

    def rhs(e, y):
        s, sd = y[0], y[1]
        a = e / eta_s
        om2 = a * a + 6.0 * s * s
        return [sd, -2.0 * s ** 3, y[3], -om2 * y[2], y[5], -om2 * y[4]]
    # background phase at e0: start at a turning point (s = 1); the phase is irrelevant for the
    # e-fold count at slow crossing and is varied below at fast crossing
    out = []
    for ph in (0.0, 0.25, 0.5):
        # advance the flat orbit by ph periods to set the phase
        s0 = solve_ivp(lambda t, y: [y[1], -2 * y[0] ** 3], (0, ph * T_ORB), [1.0, 0.0], **RT).y[:, -1] if ph else [1.0, 0.0]
        om0 = np.sqrt(a0 * a0 + 3.0)
        y0 = [s0[0], s0[1], 1.0, 0.0, 0.0, 1.0]
        y = solve_ivp(rhs, (e0, e1), y0, **RT).y[:, -1]
        om1 = np.sqrt(a1 * a1 + 3.0)
        Mt = np.array([[y[2], y[4]], [y[3], y[5]]])
        S0 = np.diag([np.sqrt(om0), 1 / np.sqrt(om0)])
        S1 = np.diag([np.sqrt(om1), 1 / np.sqrt(om1)])
        out.append(np.log(np.linalg.svd(S1 @ Mt @ np.linalg.inv(S0), compute_uv=False)[0]))
    return out


def frozen_rate(a):
    def rhs(t, y):
        s = y[0]
        om2 = a * a + 6.0 * s * s
        return [y[1], -2.0 * s ** 3, y[3], -om2 * y[2], y[5], -om2 * y[4]]
    y = solve_ivp(rhs, (0, T_ORB), [1.0, 0.0, 1.0, 0.0, 0.0, 1.0], **RT).y[:, -1]
    mu = np.max(np.abs(np.linalg.eigvals(np.array([[y[2], y[4]], [y[3], y[5]]]))))
    return np.log(mu) / T_ORB


def main():
    t0 = time.time()
    print(f"frw_ghost.py  python {sys.version.split()[0]}  numpy {np.__version__}  ({sys.executable})")
    print(f"band (frozen a): a in ({A_LO:.4f}, {A_HI:.4f}); flat-orbit period {T_ORB:.4f}")
    aa = np.linspace(A_LO, A_HI, 41)[1:-1]
    rates = np.array([frozen_rate(a) for a in aa])
    mu_mean = rates.mean()
    print(f"frozen-a Floquet rate across the band: max {rates.max():.4f}, mean {mu_mean:.4f} (per unit eta)")
    print(f"\n{'eta_s':>8} {'H/M at band':>12} {'eta in band':>12} {'mu_mean*eta':>12} | e-folds gained (3 background phases)")
    for eta_s in (1, 3, 10, 30, 100, 300, 1000, 3000):
        efs = crossing(eta_s)
        Hm = 1.0 / (eta_s * 1.8 ** 2)
        dn = (A_HI - A_LO) * eta_s
        print(f"{eta_s:8d} {Hm:12.2e} {dn:12.2f} {mu_mean * dn:12.2f} | " + "  ".join(f"{e:8.3f}" for e in efs), flush=True)
    print(f"\n({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
