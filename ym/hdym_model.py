#!/usr/bin/env python3
"""
hdym_model.py -- Python side of csrc/hdym_model.h, the one copy of the model.

Three things live here:

  header()          the text of csrc/hdym_model.h, for prepending to a CuPy RawKernel source
                    (together with cuda_prelude()) so the GPU kernels compile the same
                    characters as the CPU ones.
  c_prelude(),      the lines a C / CUDA source needs before `#include`-ing or pasting the
  cuda_prelude()    header (math.h and the DEVFN definition).
  lib()             the compiled csrc/hdym_model.c loaded through ctypes, with the vectorised
                    entry points wrapped as numpy functions: V, gradV, hessV, hessV_dot_lin,
                    kick_lin, yoshida_lin, kick_nl, yoshida_nl, yoshida_ym, V2d, gradV2d,
                    kick2d, yoshida2d, offdiag_Kv, offdiag_mass, y4.

Build: gcc through cbuild.check_call into csrc/build/ (created on first use; add it to
.gitignore).  Rebuilt whenever hdym_model.h or hdym_model.c is newer than the library.
Set CC to pick a compiler that is not called gcc.

Array convention for the wrappers: component-major, shape (ncomp, N), C-contiguous, as in the
kernels' state layout.  Functions that act in place on a state array return the same array.

This module is deliberately free of numpy *implementations* of the model: the numpy reference
copies in hdym_core.py are the Python home of the same formulas, and tools/drift_census.py
checks the two against each other.
"""
import ctypes
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
CSRC = os.path.join(_HERE, 'csrc')
HEADER = os.path.join(CSRC, 'hdym_model.h')
SOURCE = os.path.join(CSRC, 'hdym_model.c')
BUILD = os.path.join(CSRC, 'build')


def header():
    with open(HEADER, encoding='utf-8') as fh:
        return fh.read()


def c_prelude():
    return "#include <math.h>\n#define DEVFN static inline\n"


def cuda_prelude():
    return "#define DEVFN __device__ __forceinline__\n"


def c_source_with_header(body):
    """C translation unit: prelude + header text + body.  The header is pasted rather than
    #include'd so that the C and CUDA paths are built the same way and no include path has
    to be passed to the compiler."""
    return c_prelude() + header() + "\n" + body


def cuda_source_with_header(body):
    return cuda_prelude() + header() + "\n" + body


def _libname():
    return 'hdym_model' + ('.dll' if sys.platform == 'win32' else '.so')


def build(force=False, verbose=True):
    """Compile csrc/hdym_model.c into csrc/build/ if needed; return the library path."""
    sys.path.insert(0, _HERE)
    import cbuild
    os.makedirs(BUILD, exist_ok=True)
    lib = os.path.join(BUILD, _libname())
    stale = (force or not os.path.exists(lib)
             or os.path.getmtime(lib) < max(os.path.getmtime(HEADER), os.path.getmtime(SOURCE)))
    if stale:
        cc = os.environ.get('CC', 'gcc')
        cmd = [cc, '-O3', '-shared', '-fPIC', '-I', CSRC, '-o', lib, SOURCE, '-lm']
        cbuild.check_call(cmd)
        if verbose:
            print(f"hdym_model: built {lib} with {cc}", flush=True)
    return lib


_lib = None


class _Lib:
    def __init__(self, path):
        self.path = path
        L = ctypes.CDLL(path)
        dp = ctypes.POINTER(ctypes.c_double)
        d, i = ctypes.c_double, ctypes.c_int
        sig = {
            'hdym_y4': [dp],
            'hdym_V_n': [dp, dp, i],
            'hdym_gradV_n': [dp, dp, i],
            'hdym_hessV_n': [dp, dp, i],
            'hdym_hessV_dot_lin_n': [dp, dp, dp, d, i],
            'hdym_kick_lin_n': [dp, d, d, i],
            'hdym_yoshida_lin_n': [dp, d, d, i, i],
            'hdym_kick_nl_n': [dp, d, i],
            'hdym_yoshida_nl_n': [dp, d, i, i],
            'hdym_yoshida_ym_n': [dp, d, i, i],
            'hdym2d_V_n': [dp, dp, d, i],
            'hdym2d_gradV_n': [dp, dp, d, i],
            'hdym2d_kick_lin_n': [dp, d, d, d, i],
            'hdym2d_yoshida_lin_n': [dp, d, d, d, i, i],
            'hdym_offdiag_Kv_n': [dp, dp, dp, dp, d, i],
            'hdym_offdiag_mass_n': [dp, dp, i, i],
        }
        for name, argtypes in sig.items():
            f = getattr(L, name)
            f.argtypes = argtypes
            f.restype = None
        self.L = L
        self._dp = dp

    # -- helpers
    def _arr(self, a, ncomp=None):
        a = np.ascontiguousarray(a, dtype=np.float64)
        if ncomp is not None and (a.ndim != 2 or a.shape[0] != ncomp):
            raise ValueError(f"expected shape ({ncomp}, N), got {a.shape}")
        return a

    def _p(self, a):
        return a.ctypes.data_as(self._dp)

    # -- wrappers
    def y4(self):
        out = np.empty(7)
        self.L.hdym_y4(self._p(out))
        return out[:4].copy(), out[4:].copy()

    def V(self, q):
        q = self._arr(q, 3); out = np.empty(q.shape[1])
        self.L.hdym_V_n(self._p(q), self._p(out), q.shape[1]); return out

    def gradV(self, q):
        q = self._arr(q, 3); out = np.empty_like(q)
        self.L.hdym_gradV_n(self._p(q), self._p(out), q.shape[1]); return out

    def hessV(self, q):
        q = self._arr(q, 3); out = np.empty((9, q.shape[1]))
        self.L.hdym_hessV_n(self._p(q), self._p(out), q.shape[1]); return out

    def hessV_dot_lin(self, q, w, od=1.0):
        q = self._arr(q, 3); w = self._arr(w, 3); out = np.empty_like(q)
        self.L.hdym_hessV_dot_lin_n(self._p(q), self._p(w), self._p(out), float(od), q.shape[1])
        return out

    def kick_lin(self, S, h, od=1.0):
        S = self._arr(S, 12); self.L.hdym_kick_lin_n(self._p(S), float(h), float(od), S.shape[1]); return S

    def yoshida_lin(self, S, h, od=1.0, nsteps=1):
        S = self._arr(S, 12)
        self.L.hdym_yoshida_lin_n(self._p(S), float(h), float(od), int(nsteps), S.shape[1]); return S

    def kick_nl(self, S, h):
        S = self._arr(S, 12); self.L.hdym_kick_nl_n(self._p(S), float(h), S.shape[1]); return S

    def yoshida_nl(self, S, h, nsteps=1):
        S = self._arr(S, 12); self.L.hdym_yoshida_nl_n(self._p(S), float(h), int(nsteps), S.shape[1]); return S

    def yoshida_ym(self, S, h, nsteps=1):
        S = self._arr(S, 6); self.L.hdym_yoshida_ym_n(self._p(S), float(h), int(nsteps), S.shape[1]); return S

    def V2d(self, q, mu2):
        q = self._arr(q, 2); out = np.empty(q.shape[1])
        self.L.hdym2d_V_n(self._p(q), self._p(out), float(mu2), q.shape[1]); return out

    def gradV2d(self, q, mu2):
        q = self._arr(q, 2); out = np.empty_like(q)
        self.L.hdym2d_gradV_n(self._p(q), self._p(out), float(mu2), q.shape[1]); return out

    def kick2d(self, S, h, od, mu2):
        S = self._arr(S, 10)
        self.L.hdym2d_kick_lin_n(self._p(S), float(h), float(od), float(mu2), S.shape[1]); return S

    def yoshida2d(self, S, h, od, mu2, nsteps=1):
        S = self._arr(S, 10)
        self.L.hdym2d_yoshida_lin_n(self._p(S), float(h), float(od), float(mu2), int(nsteps), S.shape[1])
        return S

    def offdiag_Kv(self, q, p, v, cj):
        q = self._arr(q, 3); p = self._arr(p, 3); v = self._arr(v, 8); out = np.empty_like(v)
        self.L.hdym_offdiag_Kv_n(self._p(q), self._p(p), self._p(v), self._p(out), float(cj), q.shape[1])
        return out

    def offdiag_mass(self, q, name):
        q = self._arr(q, 3); out = np.empty((4, q.shape[1]))
        self.L.hdym_offdiag_mass_n(self._p(q), self._p(out), int(name), q.shape[1]); return out


def lib(force_build=False):
    """The compiled model, built on first call."""
    global _lib
    if _lib is None or force_build:
        _lib = _Lib(build(force=force_build))
    return _lib


if __name__ == '__main__':
    m = lib()
    YC, YD = m.y4()
    print("Yoshida-4 as compiled:", YC, YD, " sums:", YC.sum(), YD.sum())
    rng = np.random.default_rng(0)
    q = rng.standard_normal((3, 4))
    print("V(q)      ", m.V(q))
    print("grad V(q) ", m.gradV(q))
    print("library   ", m.path)
