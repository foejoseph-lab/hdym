# Fig 1: lambda_perp(E).  Run from build/.
#   diagonal ghost : per-component estimator, all energies, from
#                    ../results_diag_percomp/transverse_gpu.csv (E <= 0.03, N = 2e5) and
#                    ../results_diag_percomp_hi/transverse_gpu.csv (E >= 0.05, N = 1e5)
#   vector ghost   : ../results_vec_percomp/transverse_gpu.csv (E <= 0.03, N = 2e5) plus the
#                    earlier production scan at E >= 0.05 (values below)
#   re-run         : ../transverse_E*.npz if present (open symbols), else skipped
import csv, glob, os
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9})


def read_csv(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    E = np.array([float(r['E']) for r in rows]); lam = np.array([float(r['lam_mean']) for r in rows])
    o = np.argsort(E); return E[o], lam[o]


Ev_new, vec_new = read_csv('../results_vec_percomp/transverse_gpu.csv')
Ed_lo, dg_lo = read_csv('../results_diag_percomp/transverse_gpu.csv')
Ed_hi, dg_hi = read_csv('../results_diag_percomp_hi/transverse_gpu.csv')
if not np.any(np.isclose(Ed_lo, 0.03)):      # the E = 0.03 run (N = 1e5) preceded the overnight CSV
    Ed_lo = np.r_[Ed_lo, 0.03]; dg_lo = np.r_[dg_lo, 4.364e-6]
Ed = np.r_[Ed_lo, Ed_hi]; dg = np.r_[dg_lo, dg_hi]
E_hi = np.array([0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0])
vec_hi = np.array([4.53e-5, 1.17e-4, 4.31e-4, 1.94e-3, 5.22e-3, 1.65e-2, 3.03e-2, 8.79e-2])
Ev = np.r_[Ev_new, E_hi]; vec = np.r_[vec_new, vec_hi]
files = glob.glob('../transverse_E*.npz')
rr = {float(np.load(f)['E']): (float(np.load(f)['lam_v']), float(np.load(f)['lam_d'])) for f in files}
Er = np.array(sorted(rr)); rv = np.array([rr[e][0] for e in Er]); rd = np.array([rr[e][1] for e in Er])
if not files:
    print('note: no ../transverse_E*.npz found; re-run (open) symbols omitted')

Ee = np.logspace(-2.4, 0.1, 100)
fig, ax = plt.subplots(figsize=(4.6, 3.6))
ax.loglog(Ev, vec, 'o', color='C0', ms=5, label='vector ghost')
ax.loglog(Ed, dg, 's', color='C1', ms=5, label='diagonal ghost')
if files:
    ax.loglog(Er, rv, 'o', mfc='none', mec='C0', ms=9, label='vector, independent re-run')
    ax.loglog(Er, rd, 's', mfc='none', mec='C1', ms=9, label='diagonal, independent re-run')
ax.loglog(Ee, 0.0120 * Ee**2.25, '--', lw=1.2, color='k', label=r'derived $(32\pi^3/105\sqrt{2}Z)\,E^{9/4}$')
ax.loglog(Ee, 0.38 * Ee**0.25, ':', color='gray', lw=1, label=r'$\lambda_{\rm YM}=0.38E^{1/4}$')
ax.set_xlabel(r'$E$  (units of $M^4/g^2$)'); ax.set_ylabel(r'$\lambda_\perp$  (units of $M$)')
ax.legend(fontsize=6.5, loc='upper left'); ax.set_ylim(1e-7, 2)
plt.tight_layout(); plt.savefig('fig_lambda.pdf'); print('fig_lambda.pdf written')
for e, d in zip(Ed, dg):
    v = vec[np.argmin(np.abs(Ev - e))]
    print(f'E={e:<6g} vec={v:.3e} diag={d:.3e} vec/diag={v/d:5.2f} diag/pred={d/(0.012*e**2.25):.3f}')
