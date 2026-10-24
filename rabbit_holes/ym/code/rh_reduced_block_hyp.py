"""rh_reduced_block_hyp.py -- where along the orbit the light reduced mode is soft or hyperbolic, from the committed
frozen spectra of rh_reduced_block.py (results/rh_reduced_block_E*_cj*.npz: 'inst' = the four eigenvalues of the
restricted 4-dim generator at 400 times over one period, 'te' = those times), and how that compares with the
instantaneous phase rate of the Floquet solutions (results/rh5_instfreq.txt).  Answers open item (i) of
RH_reduced_block.md and settles the turning-point / bounce question of RH27_averaging.md item 10.  No integration.

A frozen eigenvalue counts as hyperbolic when |Im| < 1e-3 and |Re| > 0.1 (the gauge terms give real parts <= 0.1
where the spectrum is elliptic).  Printed per (E, c_j): the light frequency at t = 0 and its minimum, the fraction of
the period that is hyperbolic, the integral of the hyperbolic rate over the period (what the frozen picture would
give as log-growth per period if it were followed adiabatically), and the hyperbolic windows in t/T.  Then the
overlay at the 25 times of rh5_instfreq.txt.  Usage (from rabbit_holes/ym/): python code/rh_reduced_block_hyp.py
"""
import glob, os, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')


def windows(mask, tf):
    out, inwin = [], False
    for i, h in enumerate(mask):
        if h and not inwin:
            start, inwin = tf[i], True
        if not h and inwin:
            out.append((start, tf[i])); inwin = False
    if inwin:
        out.append((start, tf[-1]))
    # merge windows separated by a single sample (eigenvalue-pairing flicker)
    merged = []
    for a, b in out:
        if merged and a - merged[-1][1] < 2.5 / len(tf):
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    return merged


def main():
    print(__doc__.split('\n')[0])
    print("\nfrozen spectrum of the reduced 4-dim generator along the orbit; hyperbolic = |Im| < 1e-3 and |Re| > 0.1")
    print(f" {'E':>4} {'c_j':>5} | light f(t=0)   min | heavy f(t=0) | hyperbolic fraction | int Re dt (hyp) | windows t/T")
    files = sorted(glob.glob(os.path.join(RES, 'rh_reduced_block_E*_cj*.npz')),
                   key=lambda f: tuple(float(x) for x in re.search(r'E([\d.]+)_cj([\d.]+)\.npz', f).groups()))
    spectra = {}
    for f in files:
        E, cj = (float(x) for x in re.search(r'E([\d.]+)_cj([\d.]+)\.npz', f).groups())
        d = np.load(f); ev = d['inst']; te = d['te']; tf = te / te.max()
        light = np.array([np.sort(np.abs(l.imag))[0] for l in ev])
        heavy = np.array([np.sort(np.abs(l.imag))[-1] for l in ev])
        hre = np.array([max([abs(l.real) for l in evv if abs(l.imag) < 1e-3], default=0.0) for evv in ev])
        hyp = hre > 0.1
        integ = np.trapezoid(np.where(hyp, hre, 0.0), te)
        ws = ' '.join(f'[{a:.2f},{b:.2f}]' for a, b in windows(hyp, tf))
        print(f" {E:4.1f} {cj:5.2f} |   {light[0]:5.3f}   {light.min():5.3f} |    {heavy[0]:5.3f}     |        {hyp.mean():4.2f}         |     {integ:6.2f}      | {ws}")
        if cj == 1.0:
            spectra[E] = (tf, light, heavy, hre)

    inst_path = os.path.join(RES, 'rh5_instfreq.txt')
    if not os.path.exists(inst_path):
        return
    txt = open(inst_path, encoding='utf-8').read()
    for b in re.split(r'\nE=', txt)[1:]:
        E = float(b.split(',')[0])
        if E not in spectra:
            continue
        tf, light, heavy, hre = spectra[E]
        print(f"\nE = {E}: frozen light / heavy frequency (and hyperbolic |Re|) of the c_j = 1 generator against the phase rate of the two Floquet pairs")
        print("  t/T  |   x       y   | frozen light  heavy  |Re| | pair 1   pair 2")
        for line in b.splitlines()[1:]:
            m = re.match(r'\s+([\d.]+) \| ([+\-\d.]+) \| ([+\-\d.]+) \|\s+([+\-\d.]+) \|\s+([+\-\d.]+)', line)
            if not m:
                continue
            tt, x, y, p1, p2 = (float(v) for v in m.groups())
            i = np.argmin(np.abs(tf - tt))
            print(f"  {tt:5.3f} | {x:+6.3f} {y:+6.3f} |    {light[i]:5.3f}  {heavy[i]:5.3f}  {hre[i]:5.2f} | {p1:+6.3f}  {p2:+6.3f}")


if __name__ == '__main__':
    main()
