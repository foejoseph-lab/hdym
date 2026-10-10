"""Search for polynomial first integrals F of degree <= d in the 12 phase-space variables (p,w,P,W)
of H = 1/2|P|^2 - 1/2|W|^2 + V(p-w) + w.gradV(p-w) - 1/2|w|^2, by requiring {H,F}=0 at many random points.
Null space of the resulting linear map = space of integrals (contains 1 and H, and H^2 if d>=8).
Positive control: separable V_sep = (x^4+y^4+z^4)/4 -> three partial energies are integrals."""
import numpy as np, itertools, sys
d = int(sys.argv[1]) if len(sys.argv) > 1 else 4
rng = np.random.default_rng(0)
# monomials in 12 variables of total degree <= d
mons = [e for k in range(d+1) for e in itertools.product(range(k+1), repeat=12) if sum(e) == k]
mons = np.array(sorted(set(map(tuple, mons)))); M = len(mons)
def make_H(sep):
    def V(q):
        x, y, z = q.T
        return (x**4+y**4+z**4)/4 if sep else 0.5*(x*x*y*y + y*y*z*z + z*z*x*x)
    def gradV(q):
        x, y, z = q.T
        return np.stack([x**3, y**3, z**3], 1) if sep else np.stack([x*(y*y+z*z), y*(x*x+z*z), z*(x*x+y*y)], 1)
    def hessw(q, w):
        x, y, z = q.T; wx, wy, wz = w.T
        if sep: return np.stack([3*x*x*wx, 3*y*y*wy, 3*z*z*wz], 1)
        return np.stack([(y*y+z*z)*wx + 2*x*y*wy + 2*x*z*wz, 2*x*y*wx + (x*x+z*z)*wy + 2*y*z*wz, 2*x*z*wx + 2*y*z*wy + (x*x+y*y)*wz], 1)
    def dH(X):  # returns dH/dp, dH/dw, dH/dP, dH/dW  (each n x 3)
        p, w, P, W = X[:, 0:3], X[:, 3:6], X[:, 6:9], X[:, 9:12]; q = p - w
        Hp = gradV(q) + hessw(q, w); Hw = -(hessw(q, w) + w)
        return Hp, Hw, P, -W
    return dH
def run(sep, npts=6000):
    dH = make_H(sep)
    X = rng.normal(size=(npts, 12))
    Hp, Hw, HP, HW = dH(X)
    # {H,F} = sum_i (dH/dp_i dF/dP_i - dH/dP_i dF/dp_i) + (dH/dw_i dF/dW_i - dH/dW_i dF/dw_i)
    # build for each monomial: dF/dx_j = e_j x^(e - 1_j)
    logX = np.log(np.abs(X)); sgn = np.sign(X)
    A = np.empty((npts, M))
    for m, e in enumerate(mons):
        val = np.zeros(npts)
        for j in range(12):
            if e[j] == 0: continue
            ee = e.copy(); ee[j] -= 1
            dm = e[j] * np.prod(X**ee, axis=1)
            if j < 3:   val -= HP[:, j] * dm          # -dH/dP_i dF/dp_i
            elif j < 6: val -= HW[:, j-3] * dm        # -dH/dW_i dF/dw_i
            elif j < 9: val += Hp[:, j-6] * dm        # +dH/dp_i dF/dP_i
            else:       val += Hw[:, j-9] * dm        # +dH/dw_i dF/dW_i
        A[:, m] = val
    # column scaling for conditioning
    sc = np.linalg.norm(A, axis=0); sc[sc == 0] = 1; A /= sc
    s = np.linalg.svd(A, compute_uv=False)
    return s
for sep in (True, False):
    s = run(sep)
    tail = s[-8:] / s[0]
    nnull = int(np.sum(s/s[0] < 1e-9))
    print(f"{'separable control' if sep else 'Lee-Wick YM'}: {M} monomials up to degree {d}; smallest singular values / largest:",
          " ".join(f"{v:.1e}" for v in tail), f" -> null-space dimension {nnull}")
