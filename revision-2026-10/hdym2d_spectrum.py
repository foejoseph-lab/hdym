"""hdym2d_spectrum.py -- two-sided power spectrum of the diagonal Hessian entry K_xx = y^2
(and of K_zz = x^2 + y^2, the decoupled z-component's driver) along trajectories of the
regulated two-field model at E = 1.

Prediction (NOTE_2d.md, Sec. 3): the tail of s_xx is purely the resonant x-channel chirp,
    u^5 s_xx(u) = A_res = (256 pi / (15 sqrt2)) nu0 = 37.93 nu0,   nu0 = 4 pi / Z2(mu1),
up to the regulator cutoff at u ~ 2 sqrt2/mu1 (and the regulator correction factor
(1 - mu^2 D^2/2)^{5/2} with D = u/2, which is applied below).  In 2-d the diagonal entries
receive no detuned power from the other channel (in the y-channel K_xx = y^2 is slow), so
the measured tail should equal A_res -- unlike 3-d, where the resonant share was 28%.

Welch / Hann estimate on segments of `seg` time units, sampled every `every` steps.
Usage: python hdym2d_spectrum.py --mu1 0.1 --N 128 --T 3000
"""
import argparse
import numpy as np
from hdym2d_nu0 import Z2, V, step


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mu1', type=float, default=0.1)
    ap.add_argument('--N', type=int, default=128)
    ap.add_argument('--T', type=float, default=3000.0)
    ap.add_argument('--tburn', type=float, default=500.0)
    ap.add_argument('--dt', type=float, default=0.01)
    ap.add_argument('--every', type=int, default=2)
    ap.add_argument('--seg', type=float, default=163.84)
    ap.add_argument('--seed', type=int, default=11)
    a = ap.parse_args()
    mu2 = a.mu1 ** 2
    rng = np.random.default_rng(a.seed)
    q = rng.standard_normal((2, a.N))
    p = rng.standard_normal((2, a.N))
    s = np.ones(a.N)
    for _ in range(60):
        E = 0.5 * s ** 4 * np.sum(p * p, 0) + V(s * q, mu2)
        s *= (1.0 / E) ** 0.25 * 0.5 + 0.5
    q *= s
    p *= s ** 2
    nburn = int(a.tburn / a.dt)
    for _ in range(nburn):
        q, p = step(q, p, a.dt, mu2)
    nrec = int(a.T / a.dt / a.every)
    Kxx = np.empty((a.N, nrec))
    Kzz = np.empty((a.N, nrec))
    for k in range(nrec):
        for _ in range(a.every):
            q, p = step(q, p, a.dt, mu2)
        Kxx[:, k] = q[1] ** 2 + mu2
        Kzz[:, k] = q[0] ** 2 + q[1] ** 2 + mu2
    dts = a.dt * a.every
    M = int(round(a.seg / dts))
    win = np.hanning(M)
    Teff = dts * np.sum(win ** 2)
    nseg = nrec // M
    u = 2 * np.pi * np.fft.rfftfreq(M, dts)

    def spec(K):
        S = np.zeros(len(u))
        for i in range(a.N):
            for j in range(nseg):
                seg = K[i, j * M:(j + 1) * M]
                seg = seg - seg.mean()
                X = dts * np.fft.rfft(win * seg)
                S += np.abs(X) ** 2 / Teff
        return S / (a.N * nseg)

    Sxx = spec(Kxx)
    Szz = spec(Kzz)
    Z = Z2(a.mu1)
    nu0 = 4 * np.pi / Z
    Ares = 256 * np.pi / (15 * np.sqrt(2)) * nu0
    print(f"mu1={a.mu1}  Z2={Z:.2f}  nu0={nu0:.4f}  predicted A_res = u^5 s_xx = {Ares:.3f}")
    print(f"{a.N} trajectories x {nseg} segments of {M * dts:.1f} time units; du = {u[1]:.4f}")
    print("   u-band      u^5 s_xx    /A_res   /(A_res*corr)    u^5 s_zz / A_res")
    ucut = 2 * np.sqrt(2) / a.mu1
    for lo, hi in ((3, 4), (4, 5), (5, 6), (6, 8), (8, 10), (10, 12), (12, 16), (16, 20), (20, 24), (24, 28), (28, 36), (36, 48)):
        if lo >= ucut:
            break
        m = (u >= lo) & (u < hi)
        um = u[m]
        D = um / 2
        corr = np.clip(1 - 0.5 * mu2 * D * D, 0, None) ** 2.5
        v = np.mean(um ** 5 * Sxx[m])
        vc = np.mean(um ** 5 * Sxx[m] / np.where(corr > 0, corr, np.nan))
        vz = np.mean(um ** 5 * Szz[m])
        print(f"  {lo:3d}-{hi:3d}    {v:9.3f}   {v / Ares:6.3f}      {vc / Ares:6.3f}          {vz / Ares:7.3f}")
    np.savez(f"spec_2d_mu{a.mu1}.npz", u=u, Sxx=Sxx, Szz=Szz, mu1=a.mu1, Z2=Z, nu0=nu0)


if __name__ == '__main__':
    main()
