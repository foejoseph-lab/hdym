"""make_fig_Q.py -- Q = sigma^2 T / lambda from the per-trajectory finite-time exponents (lam_i_*.npy).
For growth that arrives in discrete events Q is the typical log-gain per event and is independent
of T; this reading requires lambda*T >> 1 (otherwise sigma is dominated by bounded O(1) wobble that
scales as 1/T, not diffusively), so only runs with lambda*T >= 1.5 e-folds are plotted.  Errors are
bootstrap over trajectories.  Run from the repository root with the result folders as laid out below.
"""
import glob, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9})
rng = np.random.default_rng(0)
PC = 'results_percomp'      # unpacked results_percomp.7z (folders results_diag_percomp, results_diag_percomp_hi)
R1 = 'results_rung1'        # off-diagonal runs (lam_i_E*.npy, lam_i_E*_cj0.npy)
LAMT_MIN = 1.5

def Q(lam, T, nb=300):
    lam = np.asarray(lam, float); n = lam.size
    q = lam.var(ddof=1) * T / lam.mean()
    bs = [(lambda s: s.var(ddof=1) * T / s.mean())(lam[rng.integers(0, n, n)]) for _ in range(nb)]
    return q, float(np.std(bs)), lam.mean() * T

def series(files_T, comp=None):
    out = []
    for f, E, T in files_T:
        L = np.load(f); L = L[comp] if comp is not None else L
        q, e, lt = Q(L, T)
        if lt >= LAMT_MIN: out.append((E, q, e))
    return np.array(sorted(out))

diag = []
for d in ('results_diag_percomp', 'results_diag_percomp_hi'):
    for f in glob.glob(f'{PC}/{d}/lam_i_E*.npy'):
        E = float(f.split('_E')[1][:-4]); T = np.load(f.replace('lam_i', 'history'))[-1, 0]; diag.append((f, E, T))
def r1(tag):
    fs = glob.glob(f'{R1}/lam_i_{tag}.npy'); return fs[0] if fs else None
cj1 = [(r1(t), E, T) for t, E, T in (('E1', 1, 12000), ('E0.3', 0.3, 12000), ('E0.1', 0.1, 60000), ('E0.03', 0.03, 120000), ('E0.01', 0.01, 240000)) if r1(t)]
cj0 = [(r1(t), E, T) for t, E, T in (('E0.3_cj0', 0.3, 30000), ('E0.1_cj0', 0.1, 30000)) if r1(t)]
S = {'diagonal ghost': (series(diag), 's'),
     r'off-diagonal ghost, $(D^\mu F_{\mu0})^2$ off': (series(cj0, 0), 'v'),
     'off-diagonal ghost, the theory': (series(cj1, 0), '^')}
fig, ax = plt.subplots(figsize=(4.6, 3.3))
for lab, (p, mk) in S.items():
    if len(p): ax.errorbar(p[:, 0], p[:, 1], yerr=p[:, 2], fmt=mk + '-', ms=4, lw=0.8, capsize=2, label=lab)
Es = np.array([0.01, 0.4]); ax.loglog(Es, 1.9 * Es ** -0.5, 'k:', lw=0.8, label=r'$1.9\,E^{-1/2}$')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('$E$'); ax.set_ylabel(r'$Q=\sigma^2T/\lambda$  (log-gain per growth event)')
ax.legend(fontsize=7, loc='upper right'); ax.set_ylim(0.5, 40); ax.set_xlim(0.007, 1.5)
plt.tight_layout(); plt.savefig('fig_Q.pdf'); plt.savefig('fig_Q.png', dpi=150)
for lab, (p, _) in S.items(): print(lab, [(float(e), round(float(q), 2), round(float(s), 2)) for e, q, s in p])
