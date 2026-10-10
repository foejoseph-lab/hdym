#!/usr/bin/env python3
"""
cert_check.py -- exact-arithmetic checker for a periodic-orbit certificate of
Yang-Mills mechanics  V = x^2 y^2 / 2.  Driver: parses the certificate, calls the lemmas
in cert_lemmas.py (and cert_twist.py for --twist), prints the verdicts.

The trusted layer is cert_lemmas.py + cert_twist.py + this driver.  It uses only Python integers and fractions.Fraction;
no floats enter any computation (the certificate's numbers are parsed from hex-float
strings into exact dyadic rationals).  Every function below is the shadow of one
lemma a proof assistant would have to state; the names say which.

Certificate (JSON) layout
  {
    "E": "<hex>",                      energy
    "p": 24,                            Taylor order used for the point and Jacobian series
    "K": 200,                           outward-rounding grid: all intervals rounded to 2^-K
    "box":  {"x0": [lo,hi], "px0": [lo,hi], "T": [lo,hi]},       Krawczyk box X
    "center": {"x0": v, "px0": v, "T": v},                        the point x^ (inside X)
    "Y": [[..3x3 hex..]],                                         preconditioner (any matrix)
    "chainP": {"h": "<hex>", "steps": [ {"c":[4], "A":[[4x4]], "r":[4 x [lo,hi]], "Z":[4 x [lo,hi]]}, ... ]},
    "chainB": {"h": "<hex>", "steps": [ ... same ... ], "Zfinal": [4 x [lo,hi]]}
  }
  chainP: sets along the orbit from the CENTER point for time T_center = h * len(steps).
  chainB: sets along the orbit from the BOX  for time T_lo     = h * len(steps),
          then Zfinal must be a rough enclosure over [T_lo, T_hi].
  steps[k] is the set at time k*h.  Step k maps steps[k] -> steps[k+1].  steps[0] must
  contain the initial data (checked).  The checker recomputes every enclosure itself; the
  certificate only supplies the *frames* (c, A, r, Z) the generator chose.

What is proved if every check passes
  (1) there is exactly one (x0, px0, T) in X with  x(T)=x0, px(T)=px0, y(T)=0  for the
      initial condition (x0, 0, px0, sqrt(2E - px0^2))  -- a periodic orbit;
  (2) tr M of the monodromy over the orbit lies in the printed interval; if it is in (0,4)
      the orbit is linearly stable (elliptic) [Lemma T];
  (3) with --ghost, the same for the ghost monodromy  w'' = -(1 + Hess V) w:  if |tr| > 2
      on the in-plane block, the ghost is Floquet-unstable on this orbit.

Usage:  python3 cert_check.py cert.json [--ghost] [--verbose]
"""
import json, sys
from fractions import Fraction as Fr
from cert_lemmas import *          # noqa: F401,F403  (the lemma layer; one module, one grid)
import cert_lemmas as cl

# ============================================================== parsing
def px(s):
    """hex-float string ('0x1.9p+1') or rational string ('num/den') -> exact Fraction."""
    if "/" in s:
        n, d = s.split("/"); return Fr(int(n), int(d))
    return Fr(float.fromhex(s))
def piv(pair):   return (px(pair[0]), px(pair[1]))
def pvec(v):     return [piv(a) for a in v]
def pmat(M):     return [[I(px(a)) for a in row] for row in M]
def pmat_iv(M):  return [[piv(a) for a in row] for row in M]

def parse_chain(ch):
    steps = []
    for s in ch["steps"]:
        rec = {"c": [I(px(a)) for a in s["c"]], "A": pmat(s["A"]),
               "Binv": pmat(s["Binv"]), "r": pvec(s["r"]), "Z": pvec(s["Z"]) if "Z" in s else None}
        if "G" in s:                                   # extra Jacobian frames (ghost etc.)
            rec["G"] = [pmat(G) for G in s["G"]]; rec["Ginv"] = [pmat(G) for G in s["Ginv"]]
        steps.append(rec)
    out = {"h": px(ch["h"]), "steps": steps}
    if "Zfinal" in ch: out["Zfinal"] = pvec(ch["Zfinal"])
    return out

def main():
    path = sys.argv[1]
    ghost = "--ghost" in sys.argv
    verbose = "--verbose" in sys.argv
    twist_path = None
    if "--twist" in sys.argv: twist_path = sys.argv[sys.argv.index("--twist") + 1]
    cert = json.load(open(path))
    set_grid(int(cert.get("K", 200)))
    p = int(cert["p"])
    E = I(px(cert["E"]))
    bx = cert["box"]; ct = cert["center"]
    X = [piv(bx["x0"]), piv(bx["px0"]), piv(bx["T"])]
    xhat = [I(px(ct["x0"])), I(px(ct["px0"])), I(px(ct["T"]))]
    for i in range(3):
        if not inside(xhat[i], X[i]): raise CheckFailed("center not inside box")
    Y = pmat(cert["Y"])
    m_list = [0, 1] if ghost else [0]
    print(f"certificate: p={p}, K={cl.K}, box widths "
          f"x0 {float(width(X[0])):.1e}  px0 {float(width(X[1])):.1e}  T {float(width(X[2])):.1e}")

    # ---- chain P: from the center point, time T_center
    chP = parse_chain(cert["chainP"])
    py0c = isqrt(sub(scal(Fr(2), E), mul(xhat[1], xhat[1])))
    z0c = [xhat[0], I(0), xhat[1], py0c]
    print("chain P (center) ...")
    zT, _, Tp = run_chain(chP, p, [0], "P", verbose, z0c)
    if not inside(xhat[2], I(Tp)) or not inside(I(Tp), xhat[2]):
        raise CheckFailed(f"chain P total time {float(Tp)} != center T {float(xhat[2][0])}")
    Fc = [sub(zT[0], xhat[0]), sub(zT[2], xhat[1]), zT[1]]
    print("  F(center) =", [f"[{float(a[0]):.3e}, {float(a[1]):.3e}]" for a in Fc])

    # ---- chain B: from the box, time T_lo, then rough step over [T_lo, T_hi]
    chB = parse_chain(cert["chainB"])
    py0 = isqrt(sub(scal(Fr(2), E), mul(X[1], X[1])))
    z0b = [X[0], I(0), X[1], py0]
    print("chain B (box) ...")
    dirstate = None
    if twist_path:
        import cert_twist as ctw
        c1v = neg(cl._div(X[1], py0))                                  # dpy/dpx at the fixed point
        m0 = [[I(1), I(0), I(0), I(0)], [I(0), I(0), I(1), c1v], [I(0), I(0), I(0), I(1)]]   # e_x, e_2, e_py
        dirstate = ctw.DirState(m0)
    zTlo, Jaccs, Tlo = run_chain(chB, p, m_list, "B", verbose, z0b, dirstate)
    if not Tlo == X[2][0]:
        raise CheckFailed(f"chain B total time {float(Tlo)} != T_lo {float(X[2][0])}")
    delta = X[2][1] - X[2][0]
    Zf = chB["Zfinal"]
    Jfinal, fZ = lemma_rough_time_interval(zTlo, delta, Zf, Jaccs[0], 0)
    # DF over the box:  d z0 / d(x0, px0) = [[1,0],[0,0],[0,1],[0, -px0/py0]]
    dpy = neg(cl._div(X[1], py0))
    dz0 = [[I(1), I(0)], [I(0), I(0)], [I(0), I(1)], [I(0), dpy]]
    JZ = mmul(Jfinal, dz0)
    rows = [0, 2, 1]
    DF = [[JZ[r][0], JZ[r][1], fZ[r]] for r in rows]
    DF[0][0] = sub(DF[0][0], I(1)); DF[1][1] = sub(DF[1][1], I(1))

    ok, Kv = lemma_krawczyk(xhat, X, Fc, DF, Y)
    print("\nKrawczyk:")
    for name, k, x in zip(["x0 ", "px0", "T  "], Kv, X):
        print(f"  {name}: K=[{float(k[0])!r}, {float(k[1])!r}]  ratio {float(width(k)/width(x)):.3f}")
    print("  K strictly inside X :", ok)

    trM = cl._trace(Jfinal)
    verdict, s = lemma_trace_stability(trM)
    print(f"\nbase monodromy: tr M in [{float(trM[0])!r}, {float(trM[1])!r}]  "
          f"-> HI = (tr-2)/2 in [{float(s[0])/2!r}, {float(s[1])/2!r}]  {verdict}")
    if ghost:
        Jg, _ = lemma_rough_time_interval(zTlo, delta, Zf, Jaccs[1], 1)
        verdict_g, trG, disc = lemma_symplectic_4x4(Jg)
        print(f"ghost monodromy: tr in [{float(trG[0])!r}, {float(trG[1])!r}], disc in [{float(disc[0]):.4g}, {float(disc[1]):.4g}]  -> {verdict_g}")
    if twist_path:
        import cert_twist as ctw
        m_fin, b_fin, c_fin = run_chain.dir_final
        mt, bt, ctt = ctw.lemma_dir_time_interval(Zf, delta, Jfinal, dirstate.m0, b_fin, c_fin)
        env = {'x_s': X[0], 'px_s': X[1], 'py_s': py0, 'E': E}
        names_m = ['Mx', 'M2', 'Mpy']; names_b = ['Bxx', 'Bx2', 'B22', 'Bxpy', 'B2py']; names_c = ['Cxxx', 'Cxx2', 'Cx22', 'C222']
        for nm, v in zip(names_m, mt):
            for i in range(4): env[f'{nm}_{i}'] = v[i]
        for nm, v in zip(names_b, bt):
            for i in range(4): env[f'{nm}_{i}'] = v[i]
        for nm, v in zip(names_c, ctt):
            for i in range(4): env[f'{nm}_{i}'] = v[i]
        print("\ndirectional variations at T: max widths  m %.1e  b %.1e  c %.1e" % (
            max(float(width(a)) for v in mt for a in v), max(float(width(a)) for v in bt for a in v),
            max(float(width(a)) for v in ctt for a in v)))
        formula = ctw.load_formula(twist_path)
        tau, lam, ct, re_part, zb = ctw.lemma_twist(formula, env)
        print(f"return map: cos theta in [{float(ct[0])!r}, {float(ct[1])!r}]  (non-resonant to order 4: checked)")
        print(f"            lambda = ({float(lam[0][0]):.6f}..{float(lam[0][1]):.6f}) + i({float(lam[1][0]):.6f}..{float(lam[1][1]):.6f});  zetabar coeff width {float(max(width(zb[0]), width(zb[1]))):.1e}")
        print(f"            Re(conj(lambda) c1) in [{float(re_part[0]):.3e}, {float(re_part[1]):.3e}]  (0 for area-preserving)")
        print(f"BIRKHOFF TWIST tau in [{float(tau[0])!r}, {float(tau[1])!r}]  ->  "
              f"{'NONZERO (Moser twist hypothesis certified)' if tau[1] < 0 or tau[0] > 0 else 'UNDECIDED (enclosure contains 0)'}")
    print("\nRESULT:", "PROVED: unique periodic orbit in X" if ok else "NOT PROVED (Krawczyk inclusion failed)")

if __name__ == "__main__":
    try:
        main()
    except CheckFailed as e:
        print("\nCHECK FAILED:", e)
        sys.exit(1)
