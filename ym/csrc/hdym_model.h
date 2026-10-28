/* hdym_model.h -- the one copy of the model shared by every kernel in ym/.
 *
 * The Yang-Mills potential of the colour-diagonal sector, its gradient, the action of its
 * Hessian, the Yoshida-4 coefficients, the kick/drift stages of the symplectic integrators,
 * the two-field (mu-regularised) variant, and the generated off-diagonal block K(q, p).
 * Plain C.  Compiled by gcc for the CPU kernels (through cbuild.py) and prepended, as text,
 * to the CuPy RawKernel sources for the GPU kernels, so both backends compile the same
 * characters.  The header has no dependencies beyond <math.h> on the CPU side; the CUDA side
 * gets fabs/sqrt/pow/ceil from the device runtime.
 *
 * Conventions (paper/main.tex, Sec. 2):
 *   q = (x, y, z)           Yang-Mills coordinates,  V = 1/2 (x^2 y^2 + y^2 z^2 + z^2 x^2)
 *   linearised ghost        p = YM background, w = ghost, P = pdot, Pi = wdot,
 *                           wddot = -(1 + Hess V(p)) w             (hdym_gpu, hdym_gpu_landscape)
 *   full nonlinear system   state (p, w, P, W), q = p - w, pdot = P, wdot = -W,
 *                           Pdot = -(grad V(q) + Hess V(q) w),  Wdot = +(Hess V(q) w + w)
 *                                                               (hdym_escape_map and its
 *                                                                precision/diffpair/ball twins)
 *   two-field model         V = 1/2 x^2 y^2 + 1/2 mu2 (x^2 + y^2), third ghost component
 *                           decoupled with K_zz = x^2 + y^2 + mu2  (hdym_2d)
 *
 * Operation order is part of the interface.  Every expression below is written in the same
 * order as the kernel it replaces, so that a kernel built on this header reproduces the
 * committed data bit for bit on the same compiler and flags.  Two orders exist for
 * Hess V . w and both are kept on purpose:
 *   - hdym_hessV_dot_lin / hdym_kick_lin : the linearised-ghost order with the `od` switch
 *     (hx = (y2+z2)*wx + od*(2xy*wy + 2xz*wz)), as in hdym_gpu.py;
 *   - HDYM_DEFINE_KICK_NL                 : the nonlinear order (hx = ((yz*wx) + (txy*wy)) +
 *     (txz*wz)), as in hdym_escape_map.py, which hdym_precision, hdym_diffpair and hdym_ball
 *     reproduce in rounded, dual and ball arithmetic.  The macro takes the arithmetic as
 *     parameters so those three instantiate the same text instead of a copy.
 * They differ by roundoff at od = 1.  Do not "simplify" one into the other.
 *
 * DEVFN is defined by the including source: `static inline` for C, `__device__
 * __forceinline__` for CUDA.  A default is supplied for standalone compilation.
 */
#ifndef HDYM_MODEL_H
#define HDYM_MODEL_H

#ifndef DEVFN
#define DEVFN static inline
#endif

/* ---------------------------------------------------------------- Yoshida-4 coefficients
 * c = 2^(1/3), w1 = 1/(2 - c), w0 = -c/(2 - c):
 *   YC = (w1/2, (w0+w1)/2, (w0+w1)/2, w1/2),   YD = (w1, w0, w1).
 * The literals are the ones every committed kernel was built with.  Evaluating the formula
 * in double gives 0.6756035959798289, -0.17560359597982877, 1.3512071919596578 -- 1, 3 and 1
 * ulp off, three of the four (YD1 equal; hdym_core.py computes them that way).  The literals are kept
 * because they produced the data; the difference is at the 1e-16 level and only changes
 * roundoff.                                                                              */
#define YC0 0.6756035959798288
#define YC1 -0.17560359597982886
#define YC2 -0.17560359597982886
#define YC3 0.6756035959798288
#define YD0 1.3512071919596576
#define YD1 -1.7024143839193153
#define YD2 1.3512071919596576

/* ---------------------------------------------------------------- the potential */
DEVFN double hdym_V(const double* q)
{
    double x = q[0], y = q[1], z = q[2];
    return 0.5 * (x * x * y * y + y * y * z * z + z * z * x * x);
}

/* grad V, written as in hdym_gpu.kick / hdym_rung1_gpu.gradV */
DEVFN void hdym_gradV(const double* q, double* g)
{
    double x = q[0], y = q[1], z = q[2];
    double x2 = x * x, y2 = y * y, z2 = z * z;
    g[0] = x * (y2 + z2);  g[1] = y * (z2 + x2);  g[2] = z * (x2 + y2);
}

/* the 3x3 Hessian, row-major (hdym_spectrum.hess, hdym_gpu.hessV_min_eig) */
DEVFN void hdym_hessV(const double* q, double* H)
{
    double x = q[0], y = q[1], z = q[2];
    H[0] = y * y + z * z;  H[1] = 2.0 * x * y;      H[2] = 2.0 * x * z;
    H[3] = 2.0 * x * y;    H[4] = z * z + x * x;    H[5] = 2.0 * y * z;
    H[6] = 2.0 * x * z;    H[7] = 2.0 * y * z;      H[8] = x * x + y * y;
}

/* Hess V(q) . w in the linearised-ghost order.  od = 1: full Hessian (physical ghost);
   od = 0: off-diagonal entries dropped (three scalar oscillators driven by K_xx, K_yy, K_zz).
   Order as in hdym_gpu.py / hdym_gpu_landscape.py.                                       */
DEVFN void hdym_hessV_dot_lin(const double* q, const double* w, double* hw, double od)
{
    double x = q[0], y = q[1], z = q[2];
    double x2 = x * x, y2 = y * y, z2 = z * z;
    double wx = w[0], wy = w[1], wz = w[2];
    hw[0] = (y2 + z2) * wx + od * (2.0 * x * y * wy + 2.0 * x * z * wz);
    hw[1] = od * (2.0 * x * y * wx + 2.0 * y * z * wz) + (z2 + x2) * wy;
    hw[2] = od * (2.0 * x * z * wx + 2.0 * y * z * wy) + (x2 + y2) * wz;
}

/* ---------------------------------------------------------------- linearised ghost
   kick: P -= h grad V(p),  Pi -= h (w + Hess V(p) w).   Text of hdym_gpu.kick.           */
DEVFN void hdym_kick_lin(const double* p, const double* w, double* P, double* Pi, double h, double od)
{
    double x = p[0], y = p[1], z = p[2];
    double x2 = x * x, y2 = y * y, z2 = z * z;
    double gx = x * (y2 + z2), gy = y * (z2 + x2), gz = z * (x2 + y2);
    double wx = w[0], wy = w[1], wz = w[2];
    double hx = (y2 + z2) * wx + od * (2.0 * x * y * wy + 2.0 * x * z * wz);
    double hy = od * (2.0 * x * y * wx + 2.0 * y * z * wz) + (z2 + x2) * wy;
    double hz = od * (2.0 * x * z * wx + 2.0 * y * z * wy) + (x2 + y2) * wz;
    P[0] -= h * gx;  P[1] -= h * gy;  P[2] -= h * gz;
    Pi[0] -= h * (wx + hx);  Pi[1] -= h * (wy + hy);  Pi[2] -= h * (wz + hz);
}

DEVFN void hdym_yoshida_lin(double* p, double* w, double* P, double* Pi, double h, double od)
{
    int k;
    for (k = 0; k < 3; k++) { p[k] += YC0 * h * P[k]; w[k] += YC0 * h * Pi[k]; }
    hdym_kick_lin(p, w, P, Pi, YD0 * h, od);
    for (k = 0; k < 3; k++) { p[k] += YC1 * h * P[k]; w[k] += YC1 * h * Pi[k]; }
    hdym_kick_lin(p, w, P, Pi, YD1 * h, od);
    for (k = 0; k < 3; k++) { p[k] += YC2 * h * P[k]; w[k] += YC2 * h * Pi[k]; }
    hdym_kick_lin(p, w, P, Pi, YD2 * h, od);
    for (k = 0; k < 3; k++) { p[k] += YC3 * h * P[k]; w[k] += YC3 * h * Pi[k]; }
}

/* pure Yang-Mills step (hdym_rung1_gpu.yoshida): q += c h p, p -= d h grad V(q) */
DEVFN void hdym_yoshida_ym(double* q, double* p, double h)
{
    double g[3]; int k;
    for (k = 0; k < 3; k++) q[k] += YC0 * h * p[k];
    hdym_gradV(q, g); for (k = 0; k < 3; k++) p[k] -= YD0 * h * g[k];
    for (k = 0; k < 3; k++) q[k] += YC1 * h * p[k];
    hdym_gradV(q, g); for (k = 0; k < 3; k++) p[k] -= YD1 * h * g[k];
    for (k = 0; k < 3; k++) q[k] += YC2 * h * p[k];
    hdym_gradV(q, g); for (k = 0; k < 3; k++) p[k] -= YD2 * h * g[k];
    for (k = 0; k < 3; k++) q[k] += YC3 * h * p[k];
}

/* channel sub-step count: n = ceil(om dt / omega_dt_max), at least 1, at most cap (<= 0: no cap) */
DEVFN int hdym_nsub(double om, double dt, double omega_dt_max, int cap)
{
    int n = (int)ceil(om * dt / omega_dt_max);
    if (n < 1) n = 1;
    if (cap > 0 && n > cap) n = cap;
    return n;
}

/* ---------------------------------------------------------------- full nonlinear system
 * The kick of hdym_escape_map.kick_nl, as a macro over the arithmetic so that the rounded
 * (hdym_precision.kick_r), dual (hdym_diffpair.kick_d) and ball (hdym_ball.kick_b, kick_t)
 * versions are the same text.  Parameters:
 *   NAME          the function to define
 *   T             the number type (double, du, bd, tb, ...)
 *   ADD(a, b)     a + b           SUB(a, b)  a - b          MUL(a, b)  a * b
 *   SCL(c, a)     c * a with c an exact double (the 2.0 and the step h)
 * Plain double: HDYM_DEFINE_KICK_NL(hdym_kick_nl, double, HDYM_ADD, HDYM_SUB, HDYM_MUL, HDYM_SCL)
 * Rounded:      ... (kick_r, double, RA, RS, RM, RM)   (RM(c, a) = rnd(c * a))
 * Dual:         ... (kick_d, du, dadd, dsub, dmul, dscl)
 * Ball:         ... (kick_b, bd, badd, bsub, bmul, bscl);  tangent: (kick_t, tb, tadd, tsub, tmul, tscl)
 * The rounded version's macros read two extra parameters (d, sc) from the enclosing scope;
 * HDYM_DEFINE_KICK_NL_X(NAME, T, ADD, SUB, MUL, SCL, int d, double sc) adds them to the
 * signature.                                                                             */
#define HDYM_KICK_NL_BODY(T, ADD, SUB, MUL, SCL)                                              \
    T x = SUB(p[0], w[0]), y = SUB(p[1], w[1]), z = SUB(p[2], w[2]);                          \
    T x2 = MUL(x, x), y2 = MUL(y, y), z2 = MUL(z, z);                                         \
    T yz = ADD(y2, z2), zx = ADD(z2, x2), xy = ADD(x2, y2);                                   \
    T gx = MUL(x, yz), gy = MUL(y, zx), gz = MUL(z, xy);                                      \
    T wx = w[0], wy = w[1], wz = w[2];                                                        \
    T txy = MUL(SCL(2.0, x), y), txz = MUL(SCL(2.0, x), z), tyz = MUL(SCL(2.0, y), z);        \
    T hx = ADD(ADD(MUL(yz, wx), MUL(txy, wy)), MUL(txz, wz));                                 \
    T hy = ADD(ADD(MUL(txy, wx), MUL(zx, wy)), MUL(tyz, wz));                                 \
    T hz = ADD(ADD(MUL(txz, wx), MUL(tyz, wy)), MUL(xy, wz));                                 \
    P[0] = SUB(P[0], SCL(h, ADD(gx, hx)));                                                    \
    P[1] = SUB(P[1], SCL(h, ADD(gy, hy)));                                                    \
    P[2] = SUB(P[2], SCL(h, ADD(gz, hz)));                                                    \
    W[0] = ADD(W[0], SCL(h, ADD(hx, wx)));                                                    \
    W[1] = ADD(W[1], SCL(h, ADD(hy, wy)));                                                    \
    W[2] = ADD(W[2], SCL(h, ADD(hz, wz)));

#define HDYM_DEFINE_KICK_NL(NAME, T, ADD, SUB, MUL, SCL)                                      \
    DEVFN void NAME(const T* p, const T* w, T* P, T* W, double h)                             \
    { HDYM_KICK_NL_BODY(T, ADD, SUB, MUL, SCL) }

#define HDYM_DEFINE_KICK_NL_X(NAME, T, ADD, SUB, MUL, SCL, EXTRA1, EXTRA2)                    \
    DEVFN void NAME(const T* p, const T* w, T* P, T* W, double h, EXTRA1, EXTRA2)             \
    { HDYM_KICK_NL_BODY(T, ADD, SUB, MUL, SCL) }

#define HDYM_ADD(a, b) ((a) + (b))
#define HDYM_SUB(a, b) ((a) - (b))
#define HDYM_MUL(a, b) ((a) * (b))
#define HDYM_SCL(c, a) ((c) * (a))

/* the plain-double instance: hdym_escape_map.kick_nl */
HDYM_DEFINE_KICK_NL(hdym_kick_nl, double, HDYM_ADD, HDYM_SUB, HDYM_MUL, HDYM_SCL)

/* drift of the nonlinear system: p += ch P, w -= ch W  (hdym_escape_map.yoshida_nl) */
DEVFN void hdym_drift_nl(double* p, double* w, const double* P, const double* W, double ch)
{
    int k;
    for (k = 0; k < 3; k++) { p[k] += ch * P[k]; w[k] -= ch * W[k]; }
}

DEVFN void hdym_yoshida_nl(double* p, double* w, double* P, double* W, double h)
{
    hdym_drift_nl(p, w, P, W, YC0 * h);
    hdym_kick_nl(p, w, P, W, YD0 * h);
    hdym_drift_nl(p, w, P, W, YC1 * h);
    hdym_kick_nl(p, w, P, W, YD1 * h);
    hdym_drift_nl(p, w, P, W, YC2 * h);
    hdym_kick_nl(p, w, P, W, YD2 * h);
    hdym_drift_nl(p, w, P, W, YC3 * h);
}

/* ---------------------------------------------------------------- two-field model
   V = 1/2 x^2 y^2 + 1/2 mu2 (x^2 + y^2); p = (x, y), w = (wx, wy, wz) with wz the decoupled
   z-ghost of the 3-d model on z = 0.  Text of hdym_2d.kick / hdym_2d.yoshida.           */
DEVFN double hdym2d_V(const double* q, double mu2)
{
    double x = q[0], y = q[1];
    return 0.5 * x * x * y * y + 0.5 * mu2 * (x * x + y * y);
}

DEVFN void hdym2d_gradV(const double* q, double* g, double mu2)
{
    double x = q[0], y = q[1];
    double x2 = x * x, y2 = y * y;
    g[0] = x * (y2 + mu2);  g[1] = y * (x2 + mu2);
}

DEVFN void hdym2d_kick_lin(const double* p, const double* w, double* P, double* Pi, double h,
                           double od, double mu2)
{
    double x = p[0], y = p[1];
    double x2 = x * x, y2 = y * y;
    double gx = x * (y2 + mu2), gy = y * (x2 + mu2);
    double wx = w[0], wy = w[1], wz = w[2];
    double hx = (y2 + mu2) * wx + od * 2.0 * x * y * wy;
    double hy = od * 2.0 * x * y * wx + (x2 + mu2) * wy;
    double hz = (x2 + y2 + mu2) * wz;
    P[0] -= h * gx;  P[1] -= h * gy;
    Pi[0] -= h * (wx + hx);  Pi[1] -= h * (wy + hy);  Pi[2] -= h * (wz + hz);
}

DEVFN void hdym2d_yoshida_lin(double* p, double* w, double* P, double* Pi, double h, double od, double mu2)
{
    int k;
    for (k = 0; k < 2; k++) p[k] += YC0 * h * P[k];
    for (k = 0; k < 3; k++) w[k] += YC0 * h * Pi[k];
    hdym2d_kick_lin(p, w, P, Pi, YD0 * h, od, mu2);
    for (k = 0; k < 2; k++) p[k] += YC1 * h * P[k];
    for (k = 0; k < 3; k++) w[k] += YC1 * h * Pi[k];
    hdym2d_kick_lin(p, w, P, Pi, YD1 * h, od, mu2);
    for (k = 0; k < 2; k++) p[k] += YC2 * h * P[k];
    for (k = 0; k < 3; k++) w[k] += YC2 * h * Pi[k];
    hdym2d_kick_lin(p, w, P, Pi, YD2 * h, od, mu2);
    for (k = 0; k < 2; k++) p[k] += YC3 * h * P[k];
    for (k = 0; k < 3; k++) w[k] += YC3 * h * Pi[k];
}

/* ---------------------------------------------------------------- off-diagonal block
 * K(q, p) v for the (u12, u21) block of the linearised full homogeneous theory
 * (paper Sec. 8, Eq. (L2)); v = (u12, u21, u12', u21', u12'', u21'', u12''', u21''').
 * cj = 1 is the theory, cj = 0 switches the (D^mu F_mu0)^2 term off.
 * GENERATED: the body between the markers is sympy's cse output as pasted into
 * hdym_rung1_gpu.py (KBLOCK); it is reproduced here verbatim and is not hand-edited.
 * Its provenance is rung1_system.pkl / hdym_rung1.py.                                   */
DEVFN void hdym_offdiag_Kv(const double* q, const double* p, const double* v, double* out, double cj)
{
    double x = q[0], y = q[1], z = q[2], px = p[0], py = p[1], pz = p[2];
    double v0 = v[0], v1 = v[1], v2 = v[2], v3 = v[3], v4 = v[4], v5 = v[5], v6 = v[6], v7 = v[7];
    double a6, a7;
    /* --- begin generated (hdym_rung1_gpu.KBLOCK) --- */
    double x0 = pow(z, 2);
    double x1 = v0*x0;
    double x2 = pow(z, 4);
    double x3 = 4*pz*z;
    double x4 = 2*pow(pz, 2);
    double x5 = 2*x0;
    double x6 = pow(x, 2);
    double x7 = cj*x6;
    double x8 = pow(y, 3);
    double x9 = pow(x, 3);
    double x10 = x9*y;
    double x11 = 2*cj;
    double x12 = px*py;
    double x13 = pow(y, 2);
    double x14 = x13*x6;
    double x15 = x13*x7;
    double x16 = v1*x0;
    a6 = 2*cj*pow(px, 2)*v0 + 2*cj*px*v3*y + cj*v1*x*x0*y + cj*v1*x9*y + cj*v5*x*y + 2*px*py*v1 - px*v2*x*x11 + 2*px*v3*y + 2*py*v3*x + 2*v0*x0*x13 + 2*v0*x0*x6 - v0*x14 - v0*x15 - v0*x2 - v0*x4 - v1*x*x8 + v1*x*y - v1*x10 - v1*x11*x12 - v2*x3 - v4*x5 - v4*x7 - v4 + 2*v5*x*y - x1*x7 - x1;
    a7 = 2*cj*pow(py, 2)*v1 + 2*cj*py*v2*x + cj*v0*x*x0*y + cj*v0*x*x8 + cj*v4*x*y - cj*v5*x13 - cj*x13*x16 + 2*px*py*v0 + 2*px*v2*y + 2*py*v2*x - py*v3*x11*y - v0*x*x8 + v0*x*y - v0*x10 - v0*x11*x12 + 2*v1*x0*x13 + 2*v1*x0*x6 - v1*x14 - v1*x15 - v1*x2 - v1*x4 - v3*x3 + 2*v4*x*y - v5*x5 - v5 - x16;
    /* --- end generated --- */
    out[0] = v2; out[1] = v3; out[2] = v4; out[3] = v5; out[4] = v6; out[5] = v7;
    out[6] = a6; out[7] = a7;
}

/* the off-diagonal Yang-Mills mass block m (paper Eq. (L2)), for block '12', '13', '23':
   name = 0: [[z2, -xy], [-xy, z2]],  1: [[y2, -xz], [-xz, y2]],  2: [[x2, -yz], [-yz, x2]]
   (hdym_rung1_lyap.ym_mass).  Row-major 2x2.                                           */
DEVFN void hdym_offdiag_mass(const double* q, double* m, int name)
{
    double x = q[0], y = q[1], z = q[2];
    if (name == 0)      { m[0] = z * z; m[1] = -x * y; m[2] = -x * y; m[3] = z * z; }
    else if (name == 1) { m[0] = y * y; m[1] = -x * z; m[2] = -x * z; m[3] = y * y; }
    else                { m[0] = x * x; m[1] = -y * z; m[2] = -y * z; m[3] = x * x; }
}

#endif /* HDYM_MODEL_H */
