import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9})
tr = np.trapezoid

# Fig 1: lambda_perp(E)   (PLACEHOLDER POINTS: replace with full transverse_gpu.csv)
E = np.array([1.0, 0.215, 0.05, 0.03, 0.015])
vec = np.array([8.98e-2, 6.69e-3, 4.13e-5, 9.96e-6, 2.08e-6])
dg = np.array([4.0e-3, 4.42e-4, 1.95e-5, 6.21e-6, 1.38e-6])
Ee = np.logspace(-2.2, 0.1, 100)
fig, ax = plt.subplots(figsize=(4.6, 3.6))
ax.loglog(E, vec, 'o', label='vector ghost (measured)')
ax.loglog(E, dg, 's', label='diagonal ghost (measured)')
ax.loglog(Ee, 0.027 * Ee**2.24, '-', lw=1, color='C0', alpha=.6, label=r'fit $0.027\,E^{2.24}$')
ax.loglog(Ee, 0.0120 * Ee**2.25, '--', lw=1.2, color='k', label=r'derived $(32\pi^3/105\sqrt{2}Z)\,E^{9/4}$')
ax.loglog(Ee, 0.38 * Ee**0.25, ':', color='gray', lw=1, label=r'$\lambda_{\rm YM}=0.38E^{1/4}$')
ax.set_xlabel(r'$E$  (units of $M^4/g^2$)'); ax.set_ylabel(r'$\lambda_\perp$  (units of $M$)')
ax.legend(fontsize=7, loc='upper left'); ax.set_ylim(5e-7, 2)
plt.tight_layout(); plt.savefig('fig_lambda.pdf'); plt.savefig('fig_lambda.png', dpi=150)

# Fig 3: maps + uncertainty
fig, ax = plt.subplots(1, 3, figsize=(9.6, 3.1))
for k, (f, lab, xl, yl) in enumerate([('../map256.npz', '(a) YM section', '$u$', '$v$'),
                                      ('../ghost256.npz', '(b) ghost section', r'$\theta$', r'$\phi$')]):
    d = np.load(f, allow_pickle=True); t = d['t_esc']; x1, x2 = d['ax1'], d['ax2']
    im = ax[k].imshow(np.log10(t), origin='lower', cmap='magma', extent=[x1[0], x1[-1], x2[0], x2[-1]],
                      aspect='auto', interpolation='nearest')
    ax[k].set_title(lab); ax[k].set_xlabel(xl); ax[k].set_ylabel(yl)
    plt.colorbar(im, ax=ax[k], label=r'$\log_{10} t_{\rm esc}$')
for name, mk in (('ym', 'o'), ('ghost', 's')):
    d = np.load('../unc_%s.npz' % name, allow_pickle=True); dl = d['deltas'] / float(d['width'])
    j = 1; a, se = d['alphas'][j]
    ax[2].loglog(dl, d['F'][j], mk + '-', ms=3, label=r'%s: $\alpha=%.3f\pm%.3f$' % (name, a, se))
ax[2].set_xlabel(r'$\delta$ / section width'); ax[2].set_ylabel(r'$f(\delta)$, $|\Delta t_{\rm esc}|>10$')
ax[2].set_title('(c) uncertainty exponent'); ax[2].legend(fontsize=7); ax[2].set_ylim(0.1, 1.1)
plt.tight_layout(); plt.savefig('fig_maps.pdf'); plt.savefig('fig_maps.png', dpi=150)

# Fig 2: nu0
d = np.load('../nu0_E1.npz'); D = d['D']; J = d['J']; T = float(d['Ttot'])
Ds = np.geomspace(3.5, 70, 25); r = np.array([np.sum(D > x) / T for x in Ds])
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
ax[0].loglog(Ds, r, 'o', ms=3, label='measured')
ax[0].loglog(Ds, 3 * 0.1417 / (2 * Ds**2), 'k--', label=r'$3\nu_0/2D^2$, $\nu_0=8\pi^2/Z$')
ax[0].set_xlabel(r'depth $D$'); ax[0].set_ylabel(r'rate of visits with $D_{\max}>D$'); ax[0].legend(fontsize=7); ax[0].set_title('(a)')
g = J > 0
ax[1].hist(D[g] * J[g], bins=60, range=(0.5, 1.5)); ax[1].set_xlabel(r'$D_{\max}\, J$'); ax[1].set_title('(b) adiabatic relation')
plt.tight_layout(); plt.savefig('fig_nu0.pdf'); plt.savefig('fig_nu0.png', dpi=150)
print('figures written')
