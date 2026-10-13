/* hdym_model.c -- exported, vectorised entry points over hdym_model.h for ctypes.
 *
 * Every function takes arrays laid out (ncomp, N) C-contiguous, i.e. component-major as in
 * the kernels' state layout, and loops over N.  These are for checking and for Python-side
 * use (the drift census, hdym_model.py); the production kernels include hdym_model.h
 * directly and never call through here.
 *
 * Build (done by hdym_model.py through cbuild.py):
 *   gcc -O3 -shared -fPIC -o hdym_model.so hdym_model.c -lm
 */
#include <math.h>
#define DEVFN static inline
#include "hdym_model.h"

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

/* Yoshida coefficients as the compiler sees them: out[0..3] = YC, out[4..6] = YD */
EXPORT void hdym_y4(double* out)
{
    out[0] = YC0; out[1] = YC1; out[2] = YC2; out[3] = YC3;
    out[4] = YD0; out[5] = YD1; out[6] = YD2;
}

/* q (3, N) -> V (N) */
EXPORT void hdym_V_n(const double* q, double* out, int N)
{
    int i;
    for (i = 0; i < N; i++) {
        double qq[3] = { q[i], q[N + i], q[2 * N + i] };
        out[i] = hdym_V(qq);
    }
}

/* q (3, N) -> g (3, N) */
EXPORT void hdym_gradV_n(const double* q, double* out, int N)
{
    int i, k;
    for (i = 0; i < N; i++) {
        double qq[3] = { q[i], q[N + i], q[2 * N + i] }, g[3];
        hdym_gradV(qq, g);
        for (k = 0; k < 3; k++) out[k * N + i] = g[k];
    }
}

/* q (3, N) -> H (9, N), row-major 3x3 per point */
EXPORT void hdym_hessV_n(const double* q, double* out, int N)
{
    int i, k;
    for (i = 0; i < N; i++) {
        double qq[3] = { q[i], q[N + i], q[2 * N + i] }, H[9];
        hdym_hessV(qq, H);
        for (k = 0; k < 9; k++) out[k * N + i] = H[k];
    }
}

/* q (3, N), w (3, N) -> Hess V . w (3, N), linearised order, with od */
EXPORT void hdym_hessV_dot_lin_n(const double* q, const double* w, double* out, double od, int N)
{
    int i, k;
    for (i = 0; i < N; i++) {
        double qq[3] = { q[i], q[N + i], q[2 * N + i] };
        double ww[3] = { w[i], w[N + i], w[2 * N + i] }, hw[3];
        hdym_hessV_dot_lin(qq, ww, hw, od);
        for (k = 0; k < 3; k++) out[k * N + i] = hw[k];
    }
}

/* S (12, N) = (p, w, P, Pi) in place: one linearised kick of size h */
EXPORT void hdym_kick_lin_n(double* S, double h, double od, int N)
{
    int i, k;
    for (i = 0; i < N; i++) {
        double p[3], w[3], P[3], Pi[3];
        for (k = 0; k < 3; k++) { p[k] = S[k * N + i]; w[k] = S[(3 + k) * N + i];
                                  P[k] = S[(6 + k) * N + i]; Pi[k] = S[(9 + k) * N + i]; }
        hdym_kick_lin(p, w, P, Pi, h, od);
        for (k = 0; k < 3; k++) { S[(6 + k) * N + i] = P[k]; S[(9 + k) * N + i] = Pi[k]; }
    }
}

/* S (12, N) in place: nsteps linearised Yoshida steps of size h (no sub-stepping) */
EXPORT void hdym_yoshida_lin_n(double* S, double h, double od, int nsteps, int N)
{
    int i, k, s;
    for (i = 0; i < N; i++) {
        double p[3], w[3], P[3], Pi[3];
        for (k = 0; k < 3; k++) { p[k] = S[k * N + i]; w[k] = S[(3 + k) * N + i];
                                  P[k] = S[(6 + k) * N + i]; Pi[k] = S[(9 + k) * N + i]; }
        for (s = 0; s < nsteps; s++) hdym_yoshida_lin(p, w, P, Pi, h, od);
        for (k = 0; k < 3; k++) { S[k * N + i] = p[k]; S[(3 + k) * N + i] = w[k];
                                  S[(6 + k) * N + i] = P[k]; S[(9 + k) * N + i] = Pi[k]; }
    }
}

/* S (12, N) = (p, w, P, W) in place: one nonlinear kick of size h */
EXPORT void hdym_kick_nl_n(double* S, double h, int N)
{
    int i, k;
    for (i = 0; i < N; i++) {
        double p[3], w[3], P[3], W[3];
        for (k = 0; k < 3; k++) { p[k] = S[k * N + i]; w[k] = S[(3 + k) * N + i];
                                  P[k] = S[(6 + k) * N + i]; W[k] = S[(9 + k) * N + i]; }
        hdym_kick_nl(p, w, P, W, h);
        for (k = 0; k < 3; k++) { S[(6 + k) * N + i] = P[k]; S[(9 + k) * N + i] = W[k]; }
    }
}

/* S (12, N) in place: nsteps nonlinear Yoshida steps of size h (no sub-stepping) */
EXPORT void hdym_yoshida_nl_n(double* S, double h, int nsteps, int N)
{
    int i, k, s;
    for (i = 0; i < N; i++) {
        double p[3], w[3], P[3], W[3];
        for (k = 0; k < 3; k++) { p[k] = S[k * N + i]; w[k] = S[(3 + k) * N + i];
                                  P[k] = S[(6 + k) * N + i]; W[k] = S[(9 + k) * N + i]; }
        for (s = 0; s < nsteps; s++) hdym_yoshida_nl(p, w, P, W, h);
        for (k = 0; k < 3; k++) { S[k * N + i] = p[k]; S[(3 + k) * N + i] = w[k];
                                  S[(6 + k) * N + i] = P[k]; S[(9 + k) * N + i] = W[k]; }
    }
}

/* S (6, N) = (q, p) in place: nsteps pure-YM Yoshida steps of size h */
EXPORT void hdym_yoshida_ym_n(double* S, double h, int nsteps, int N)
{
    int i, k, s;
    for (i = 0; i < N; i++) {
        double q[3], p[3];
        for (k = 0; k < 3; k++) { q[k] = S[k * N + i]; p[k] = S[(3 + k) * N + i]; }
        for (s = 0; s < nsteps; s++) hdym_yoshida_ym(q, p, h);
        for (k = 0; k < 3; k++) { S[k * N + i] = q[k]; S[(3 + k) * N + i] = p[k]; }
    }
}

/* two-field: q (2, N) -> V (N); g (2, N) */
EXPORT void hdym2d_V_n(const double* q, double* out, double mu2, int N)
{
    int i;
    for (i = 0; i < N; i++) { double qq[2] = { q[i], q[N + i] }; out[i] = hdym2d_V(qq, mu2); }
}

EXPORT void hdym2d_gradV_n(const double* q, double* out, double mu2, int N)
{
    int i;
    for (i = 0; i < N; i++) {
        double qq[2] = { q[i], q[N + i] }, g[2];
        hdym2d_gradV(qq, g, mu2);
        out[i] = g[0]; out[N + i] = g[1];
    }
}

/* S (10, N) = (p[2], w[3], P[2], Pi[3]) in place: one two-field linearised kick */
EXPORT void hdym2d_kick_lin_n(double* S, double h, double od, double mu2, int N)
{
    int i, k;
    for (i = 0; i < N; i++) {
        double p[2], w[3], P[2], Pi[3];
        for (k = 0; k < 2; k++) { p[k] = S[k * N + i]; P[k] = S[(5 + k) * N + i]; }
        for (k = 0; k < 3; k++) { w[k] = S[(2 + k) * N + i]; Pi[k] = S[(7 + k) * N + i]; }
        hdym2d_kick_lin(p, w, P, Pi, h, od, mu2);
        for (k = 0; k < 2; k++) S[(5 + k) * N + i] = P[k];
        for (k = 0; k < 3; k++) S[(7 + k) * N + i] = Pi[k];
    }
}

/* S (10, N) in place: nsteps two-field Yoshida steps */
EXPORT void hdym2d_yoshida_lin_n(double* S, double h, double od, double mu2, int nsteps, int N)
{
    int i, k, s;
    for (i = 0; i < N; i++) {
        double p[2], w[3], P[2], Pi[3];
        for (k = 0; k < 2; k++) { p[k] = S[k * N + i]; P[k] = S[(5 + k) * N + i]; }
        for (k = 0; k < 3; k++) { w[k] = S[(2 + k) * N + i]; Pi[k] = S[(7 + k) * N + i]; }
        for (s = 0; s < nsteps; s++) hdym2d_yoshida_lin(p, w, P, Pi, h, od, mu2);
        for (k = 0; k < 2; k++) { S[k * N + i] = p[k]; S[(5 + k) * N + i] = P[k]; }
        for (k = 0; k < 3; k++) { S[(2 + k) * N + i] = w[k]; S[(7 + k) * N + i] = Pi[k]; }
    }
}

/* q (3, N), p (3, N), v (8, N) -> K v (8, N) */
EXPORT void hdym_offdiag_Kv_n(const double* q, const double* p, const double* v, double* out,
                              double cj, int N)
{
    int i, k;
    for (i = 0; i < N; i++) {
        double qq[3] = { q[i], q[N + i], q[2 * N + i] };
        double pp[3] = { p[i], p[N + i], p[2 * N + i] };
        double vv[8], o[8];
        for (k = 0; k < 8; k++) vv[k] = v[k * N + i];
        hdym_offdiag_Kv(qq, pp, vv, o, cj);
        for (k = 0; k < 8; k++) out[k * N + i] = o[k];
    }
}

/* q (3, N) -> m (4, N) row-major 2x2 */
EXPORT void hdym_offdiag_mass_n(const double* q, double* out, int name, int N)
{
    int i, k;
    for (i = 0; i < N; i++) {
        double qq[3] = { q[i], q[N + i], q[2 * N + i] }, m[4];
        hdym_offdiag_mass(qq, m, name);
        for (k = 0; k < 4; k++) out[k * N + i] = m[k];
    }
}
