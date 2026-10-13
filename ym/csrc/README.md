# csrc/ — the one copy of the model

C sources shared by the CPU (gcc → ctypes) and CUDA (CuPy RawKernel) kernels of `ym/`: the
Yang–Mills potential of the colour-diagonal sector, its gradient and Hessian, the Yoshida-4
coefficients, the kick and drift stages of the symplectic integrators, the two-field
(μ-regularised) variant, and the generated off-diagonal block K(q, p). Nothing here is
cryptographic or security-related; "model" means the equations of motion.

| file | what it is |
|---|---|
| `hdym_model.h` | the model. Plain C, header-only (`DEVFN` functions), no dependencies beyond `<math.h>`. Written in the **same operation order** as the kernels it replaces, so a kernel built on it reproduces the committed data bit for bit on the same compiler and flags. |
| `hdym_model.c` | exported, vectorised entry points over the header for ctypes (`hdym_V_n`, `hdym_kick_lin_n`, ...). Used by `../hdym_model.py` and `../tools/drift_census.py`; the production kernels do not call through it. |
| `build/` | where `../hdym_model.py` compiles `hdym_model.c` (git-ignored) |

## How a kernel uses it

The header is pasted into the kernel source as text, not `#include`d, so that the C and CUDA
paths are built the same way and no include path has to be passed to any compiler:

```python
import hdym_model
C_SRC    = hdym_model.c_source_with_header(BODY)      # "#include <math.h>\n#define DEVFN static inline\n" + header + BODY
CUDA_SRC = hdym_model.cuda_source_with_header(BODY)   # "#define DEVFN __device__ __forceinline__\n" + header + BODY
```

`BODY` is what is left of a script's `STEP_SRC` once its copy of the model is removed: the
per-trajectory loop, renormalisation, histogramming, escape test. It calls `hdym_kick_lin`,
`hdym_yoshida_lin`, `hdym_kick_nl`, `hdym_yoshida_nl`, `hdym_yoshida_ym`, `hdym2d_*`,
`hdym_offdiag_Kv` and the `YC*`/`YD*` macros by name.

## Two orders of Hess V · w, both kept

The linearised-ghost kernels (`hdym_gpu`, `hdym_gpu_landscape`) and the nonlinear kernels
(`hdym_escape_map` and its precision / diffpair / ball twins) sum the three terms of each
component in a different order, and the two differ at roundoff. The header keeps both:
`hdym_kick_lin` / `hdym_hessV_dot_lin` for the first, the macro `HDYM_DEFINE_KICK_NL` for the
second. The macro takes the arithmetic as parameters (`ADD`, `SUB`, `MUL`, `SCL`) because
`hdym_precision.kick_r`, `hdym_diffpair.kick_d` and `hdym_ball.kick_b` / `kick_t` are that
same text in rounded, dual and ball arithmetic; instantiating the macro makes "same
operations, same order as `kick_nl`" true by construction instead of by inspection.

## Yoshida-4 literals

`YC0 … YD2` are the literals every committed kernel was built with. The closed form
(c = 2^{1/3}, w1 = 1/(2−c), w0 = −c/(2−c)) evaluated in double differs from them by one ulp in
three of the four; `hdym_core.py` and the `early/` scripts compute it that way. The literals
stay because they produced the data. The census (`../tools/drift_census.py`) reports the
one-ulp rows as `roundoff`, which is what they are.

## Status

*(2026-10-13)* Header and census written; the census passes on the committed copies (83
bitwise, 22 roundoff ≤ 7 ulp, 0 different; CUDA copies to be run on a machine with a GPU).
No script has been moved onto the header yet; that is the next step, one script at a time,
each checked against a short fixed-seed run of the committed version before and after.
