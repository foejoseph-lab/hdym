#!/usr/bin/env python3
"""
hdym_ball.py -- certified pair escape times with ball arithmetic on difference numbers.

Each quantity of a pair (a, b = a + D) is carried as four doubles
      (v, rv, d, rd):    a_true in [v - rv, v + rv],   D_true in [d - rd, d + rd]
(componentwise balls).  The centres v, d are computed with exactly the
operations of hdym_diffpair (bit-identical; see --selftest).  The radii are
upper bounds, propagated through every add / multiply with the standard ball
rules plus a bound on each rounding error (Rump-style: 2^-52 |result| per
rounding, every radius sum inflated by 1 + 2^-46 to cover the rounding of the
radius computation itself, plus a denormal floor).  No reliance on directed
rounding modes, so the same source runs on CPU and GPU.

Why on difference numbers.  A ball around each member separately starts at
~1e-16 and grows at the Lyapunov rate, so a pair at delta < 1e-15 could never
be certified to differ: the floor returns.  In the difference algebra every
term of the D update carries a factor of D, so the radius of D grows like
(relative error of the reference) x |D|: a separation of 1e-40 is certified
with a radius of ~1e-40 x (rv / |v|), as long as the reference ball stays small.

What is certified.  The exact orbit, in real arithmetic, of the integrator
map actually used: Yoshida-type compositions with the stored double
coefficients, step h = fl(dt/n), sub-step count n chosen from the computed
centres.  Each such step is an exactly symplectic map (a composition of exact
drifts and kicks), consistent with the flow to O(dt^4) up to ~1e-16 in the
coefficients.  Roundoff is fully accounted for; the time-discretisation error
relative to the continuous Hamiltonian flow is NOT (that needs Taylor-model
enclosures, e.g. CAPD).  So: "the escape-time function of this symplectic map
has these pair differences", rigorously, at any delta the radii allow.

Escape intervals.  At each step: 'possibly escaped' if some |c| + r > R,
'certainly escaped' if some |c| - r > R.  The true exit step lies in
[first possible, first certain].  A pair is CERTIFIED DIFFERENT at threshold
T if the two intervals are separated by more than T, CERTIFIED SAME if the
union of the intervals spans <= T.  Rigorous bracket:
      f_lo = P(certified different)  <=  f_true  <=  1 - P(certified same) = f_hi.

Examples
  python hdym_ball.py --selftest
  python hdym_ball.py --E 0.3 --eps 0.01 --M 2000 --delta-min 1e-40
  python hdym_ball.py --E 0.3 --eps 0.01 --M 20000 --backend cupy
"""
import argparse
import ctypes
import json
import os
import subprocess
import sys
import tempfile
import time

import math
import numpy as np


if hasattr(math, 'fma'):                      # Python >= 3.13
    _fma = math.fma
else:
    def _two_prod(a, b):
        """Dekker/Veltkamp exact product: a*b = p + e exactly (no fma needed, no overflow assumed)."""
        p = a * b
        c = 134217729.0 * a                   # 2^27 + 1
        ah = c - (c - a); al = a - ah
        c = 134217729.0 * b
        bh = c - (c - b); bl = b - bh
        e = ((ah * bh - p) + ah * bl + al * bh) + al * bl
        return p, e

    def _fma(a, b, c):
        """a*b + c with the product exact and the sum rounded in double-word arithmetic: within
        one ulp of a true fma, which keeps the double-word selftest error far below 2^-101."""
        p, e = _two_prod(a, b)
        s = p + c
        bb = s - p
        err = (p - (s - bb)) + (c - bb)        # two_sum error of p + c
        return s + (err + e)

import hdym_diffpair as dpm
import hdym_escape_map as em
import hdym_precision as hp
import hdym_uncertainty as unc

BALL_SRC = r"""
#define YC0 0.6756035959798288
#define YC1 -0.17560359597982886
#define YC2 -0.17560359597982886
#define YC3 0.6756035959798288
#define YD0 1.3512071919596576
#define YD1 -1.7024143839193153
#define YD2 1.3512071919596576
#define MAXSUB 16384
#define U2   2.220446049250313e-16          /* 2^-52: bound on |fl(x) - x| / |fl(x)| */
#define INFL 1.0000000000000142             /* 1 + 2^-46: covers radius-arithmetic rounding */
#define ETA  3.16e-322                      /* 64 x smallest denormal: underflow cover */

/* ---------------- double-word ("double-double") arithmetic ----------------
   Algorithms and relative error bounds from Joldes, Muller & Popescu, ACM TOMS 44 (2017):
   AccurateDWPlusDW (Alg. 6) <= 3u^2/(1-4u),  DWTimesDW with FMA (Alg. 12) <= 4u^2,
   u = 2^-53.  We charge DDU = 2^-101 (> 6u^2) per operation, relative to the result.     */
#define DDU 3.944304526105059e-31          /* 2^-101 */
typedef struct { double h, l; } dd;
DEVFN dd mkd(double h, double l) { dd r; r.h = h; r.l = l; return r; }
DEVFN dd fast2sum(double a, double b) { double s = a + b; return mkd(s, b - (s - a)); }
DEVFN dd two_sum(double a, double b)
{
    double s = a + b, bb = s - a;
    return mkd(s, (a - (s - bb)) + (b - bb));
}
DEVFN dd ddadd(dd x, dd y)
{
    dd s = two_sum(x.h, y.h), t = two_sum(x.l, y.l);
    double c = s.l + t.h;
    dd v = fast2sum(s.h, c);
    double w = t.l + v.l;
    return fast2sum(v.h, w);
}
DEVFN dd ddneg(dd x) { return mkd(-x.h, -x.l); }
DEVFN dd ddsub(dd x, dd y) { return ddadd(x, ddneg(y)); }
DEVFN dd ddmul(dd x, dd y)
{
    double ch = x.h * y.h, cl1 = fma(x.h, y.h, -ch);
    double tl0 = x.l * y.l, tl1 = fma(x.h, y.l, tl0), cl2 = fma(x.l, y.h, tl1);
    return fast2sum(ch, cl1 + cl2);
}
DEVFN double ddabs(dd x) { return fabs(x.h) * (1.0 + 2.3e-16); }   /* >= |h + l| */

/* ball on difference numbers with double-word centres:
   a_true in v +- rv,  D_true in d +- rd                                                    */
typedef struct { dd v; double rv; dd d; double rd; } bd;

DEVFN bd mkb(dd v, double rv, dd d, double rd)
{ bd r; r.v = v; r.rv = rv; r.d = d; r.rd = rd; return r; }

DEVFN bd badd(bd a, bd b)
{
    dd v = ddadd(a.v, b.v), d = ddadd(a.d, b.d);
    return mkb(v, (a.rv + b.rv + DDU * ddabs(v)) * INFL + ETA,
               d, (a.rd + b.rd + DDU * ddabs(d)) * INFL + ETA);
}
DEVFN bd bsub(bd a, bd b)
{
    dd v = ddsub(a.v, b.v), d = ddsub(a.d, b.d);
    return mkb(v, (a.rv + b.rv + DDU * ddabs(v)) * INFL + ETA,
               d, (a.rd + b.rd + DDU * ddabs(d)) * INFL + ETA);
}
DEVFN bd bmul(bd a, bd b)
{
    dd v = ddmul(a.v, b.v);
    dd p1 = ddmul(a.d, b.v), p2 = ddmul(a.v, b.d), p3 = ddmul(a.d, b.d);
    dd s = ddadd(p1, p2), d = ddadd(s, p3);
    double av = ddabs(a.v), bv = ddabs(b.v), ad = ddabs(a.d), bdd = ddabs(b.d);
    double rv = (av * b.rv + bv * a.rv + a.rv * b.rv + DDU * ddabs(v)) * INFL + ETA;
    double rd = (ad * b.rv + bv * a.rd + a.rd * b.rv
               + av * b.rd + bdd * a.rv + a.rv * b.rd
               + ad * b.rd + bdd * a.rd + a.rd * b.rd
               + DDU * (ddabs(p1) + ddabs(p2) + ddabs(p3) + ddabs(s) + ddabs(d))) * INFL + ETA;
    return mkb(v, rv, d, rd);
}
DEVFN bd bscl(double c, bd a)
{
    dd cc = mkd(c, 0.0);
    dd v = ddmul(cc, a.v), d = ddmul(cc, a.d);
    return mkb(v, (fabs(c) * a.rv + DDU * ddabs(v)) * INFL + ETA,
               d, (fabs(c) * a.rd + DDU * ddabs(d)) * INFL + ETA);
}

DEVFN void kick_b(const bd* p, const bd* w, bd* P, bd* W, double h)
{
    bd x = bsub(p[0], w[0]), y = bsub(p[1], w[1]), z = bsub(p[2], w[2]);
    bd x2 = bmul(x, x), y2 = bmul(y, y), z2 = bmul(z, z);
    bd yz = badd(y2, z2), zx = badd(z2, x2), xy = badd(x2, y2);
    bd gx = bmul(x, yz), gy = bmul(y, zx), gz = bmul(z, xy);
    bd wx = w[0], wy = w[1], wz = w[2];
    bd txy = bmul(bscl(2.0, x), y), txz = bmul(bscl(2.0, x), z), tyz = bmul(bscl(2.0, y), z);
    bd hx = badd(badd(bmul(yz, wx), bmul(txy, wy)), bmul(txz, wz));
    bd hy = badd(badd(bmul(txy, wx), bmul(zx, wy)), bmul(tyz, wz));
    bd hz = badd(badd(bmul(txz, wx), bmul(tyz, wy)), bmul(xy, wz));
    P[0] = bsub(P[0], bscl(h, badd(gx, hx)));
    P[1] = bsub(P[1], bscl(h, badd(gy, hy)));
    P[2] = bsub(P[2], bscl(h, badd(gz, hz)));
    W[0] = badd(W[0], bscl(h, badd(hx, wx)));
    W[1] = badd(W[1], bscl(h, badd(hy, wy)));
    W[2] = badd(W[2], bscl(h, badd(hz, wz)));
}
DEVFN void drift_b(bd* p, bd* w, const bd* P, const bd* W, double ch)
{
    int k;
    for (k = 0; k < 3; k++) { p[k] = badd(p[k], bscl(ch, P[k])); w[k] = bsub(w[k], bscl(ch, W[k])); }
}
DEVFN void yoshida_b(bd* p, bd* w, bd* P, bd* W, double h)
{
    drift_b(p, w, P, W, YC0 * h);
    kick_b(p, w, P, W, YD0 * h);
    drift_b(p, w, P, W, YC1 * h);
    kick_b(p, w, P, W, YD1 * h);
    drift_b(p, w, P, W, YC2 * h);
    kick_b(p, w, P, W, YD2 * h);
    drift_b(p, w, P, W, YC3 * h);
}
DEVFN int nsub_b(const bd* p, const bd* w, double dt, double omd, int sel)
{
    double om = 0.0, a;
    int k, n;
    for (k = 0; k < 3; k++) {
        a = sel ? fabs((p[k].v.h + p[k].d.h) - (w[k].v.h + w[k].d.h)) : fabs(p[k].v.h - w[k].v.h);
        if (a > om) om = a;
    }
    n = (int)ceil(om * dt / omd);
    if (n < 1) n = 1;
    if (n > MAXSUB) n = MAXSUB;
    return n;
}

/* ---------------- tangent-linear Yoshida step (for the Lohner frames) ---------------- */
typedef struct { double v, t; } tg;
DEVFN tg mkt(double v, double t) { tg r; r.v = v; r.t = t; return r; }
DEVFN tg tadd(tg a, tg b) { return mkt(a.v + b.v, a.t + b.t); }
DEVFN tg tsub(tg a, tg b) { return mkt(a.v - b.v, a.t - b.t); }
DEVFN tg tmul(tg a, tg b) { return mkt(a.v * b.v, a.t * b.v + a.v * b.t); }
DEVFN tg tscl(double c, tg a) { return mkt(c * a.v, c * a.t); }
DEVFN void kick_t(const tg* p, const tg* w, tg* P, tg* W, double h)
{
    tg x = tsub(p[0], w[0]), y = tsub(p[1], w[1]), z = tsub(p[2], w[2]);
    tg x2 = tmul(x, x), y2 = tmul(y, y), z2 = tmul(z, z);
    tg yz = tadd(y2, z2), zx = tadd(z2, x2), xy = tadd(x2, y2);
    tg gx = tmul(x, yz), gy = tmul(y, zx), gz = tmul(z, xy);
    tg wx = w[0], wy = w[1], wz = w[2];
    tg txy = tmul(tscl(2.0, x), y), txz = tmul(tscl(2.0, x), z), tyz = tmul(tscl(2.0, y), z);
    tg hx = tadd(tadd(tmul(yz, wx), tmul(txy, wy)), tmul(txz, wz));
    tg hy = tadd(tadd(tmul(txy, wx), tmul(zx, wy)), tmul(tyz, wz));
    tg hz = tadd(tadd(tmul(txz, wx), tmul(tyz, wy)), tmul(xy, wz));
    P[0] = tsub(P[0], tscl(h, tadd(gx, hx)));
    P[1] = tsub(P[1], tscl(h, tadd(gy, hy)));
    P[2] = tsub(P[2], tscl(h, tadd(gz, hz)));
    W[0] = tadd(W[0], tscl(h, tadd(hx, wx)));
    W[1] = tadd(W[1], tscl(h, tadd(hy, wy)));
    W[2] = tadd(W[2], tscl(h, tadd(hz, wz)));
}
DEVFN void drift_t(tg* p, tg* w, const tg* P, const tg* W, double ch)
{
    int k;
    for (k = 0; k < 3; k++) { p[k] = tadd(p[k], tscl(ch, P[k])); w[k] = tsub(w[k], tscl(ch, W[k])); }
}
DEVFN void yoshida_t(tg* p, tg* w, tg* P, tg* W, double h)
{
    drift_t(p, w, P, W, YC0 * h); kick_t(p, w, P, W, YD0 * h);
    drift_t(p, w, P, W, YC1 * h); kick_t(p, w, P, W, YD1 * h);
    drift_t(p, w, P, W, YC2 * h); kick_t(p, w, P, W, YD2 * h);
    drift_t(p, w, P, W, YC3 * h);
}
/* out = J(x) col, J the Jacobian of n sub-steps of size dt/n at base point x */
DEVFN void tangent_step(const double* x, const double* col, double* out, double dt, int n)
{
    tg p[3], w[3], P[3], W[3];
    int k, u;
    for (k = 0; k < 3; k++) {
        p[k] = mkt(x[k], col[k]);         w[k] = mkt(x[3 + k], col[3 + k]);
        P[k] = mkt(x[6 + k], col[6 + k]); W[k] = mkt(x[9 + k], col[9 + k]);
    }
    { double h = dt / n; for (u = 0; u < n; u++) yoshida_t(p, w, P, W, h); }
    for (k = 0; k < 3; k++) {
        out[k] = p[k].t; out[3 + k] = w[k].t; out[6 + k] = P[k].t; out[9 + k] = W[k].t;
    }
}

/* Ellipsoidal enclosure, square-root form (as in a square-root Kalman filter).
   Error set of a member: e in { L u : |u|_2 <= 1 },  L 12x12 (column-major L[c*12 + k]).
   Step:  e' = J e + (defect + remainder),  |defect + remainder|_2 <= beta.
   Linear image J L is exact (no wrapping); the Minkowski sum with the ball is enclosed by
       Sigma' = (1 + 1/c) Y Y^T + (1 + c) beta^2 I,   Y = J L,  c = |Y|_F / (beta sqrt 12)
   and L' with L' L'^T = Sigma' comes from a Householder QR of [sqrt(1+1/c) Y | sqrt(1+c) beta I]^T.
   Frame arithmetic (Y, QR) is floating point; its error, relative to |Y|, is covered by
   +1e-13 |Y|_F in beta and a 1 + 1e-12 inflation.                                          */
DEVFN void ell_update(double* L, const double* x, double dt, int n, const double* defect,
                      double rem, int usej)
{
    double Y[144], A[288], beta = 0.0, tr1 = 0.0;
    int c, k, j, r;
    for (c = 0; c < 12; c++) {
        if (usej) tangent_step(x, &L[c * 12], &Y[c * 12], dt, n);
        else for (k = 0; k < 12; k++) Y[c * 12 + k] = L[c * 12 + k];
    }
    for (k = 0; k < 144; k++) tr1 += Y[k] * Y[k];
    for (k = 0; k < 12; k++) beta += (defect[k] + rem) * (defect[k] + rem);
    beta = sqrt(beta) * INFL + 1e-13 * sqrt(tr1) + ETA;
    if (tr1 == 0.0) {
        for (k = 0; k < 144; k++) L[k] = 0.0;
        for (k = 0; k < 12; k++) L[k * 12 + k] = beta;
        return;
    }
    double cc = sqrt(tr1) / (beta * sqrt(12.0));
    double sa = sqrt(1.0 + 1.0 / cc), sb = sqrt(1.0 + cc) * beta;
    /* A = M^T, 24 x 12, row-major A[row*12 + col]:  rows 0..11 = sa * Y^T, rows 12..23 = sb I */
    for (r = 0; r < 12; r++) for (k = 0; k < 12; k++) A[r * 12 + k] = sa * Y[r * 12 + k];
    for (r = 0; r < 12; r++) for (k = 0; k < 12; k++) A[(12 + r) * 12 + k] = (r == k) ? sb : 0.0;
    /* Householder QR of A (24 x 12) in place; R in the upper triangle */
    for (j = 0; j < 12; j++) {
        double nrm = 0.0;
        for (r = j; r < 24; r++) nrm += A[r * 12 + j] * A[r * 12 + j];
        nrm = sqrt(nrm);
        if (nrm == 0.0) continue;
        double alpha = A[j * 12 + j] > 0 ? -nrm : nrm;
        double vj = A[j * 12 + j] - alpha;
        double vnorm2 = vj * vj;
        for (r = j + 1; r < 24; r++) vnorm2 += A[r * 12 + j] * A[r * 12 + j];
        if (vnorm2 == 0.0) continue;
        for (c = j + 1; c < 12; c++) {
            double s = vj * A[j * 12 + c];
            for (r = j + 1; r < 24; r++) s += A[r * 12 + j] * A[r * 12 + c];
            s = 2.0 * s / vnorm2;
            A[j * 12 + c] -= s * vj;
            for (r = j + 1; r < 24; r++) A[r * 12 + c] -= s * A[r * 12 + j];
        }
        A[j * 12 + j] = alpha;
    }
    /* L' = R^T :  column c of L' = row c of R  */
    for (c = 0; c < 12; c++) for (k = 0; k < 12; k++)
        L[c * 12 + k] = (k >= c ? A[c * 12 + k] : 0.0) * (1.0 + 1e-12);
}
/* componentwise bound E_k = |row k of L|_2 ; returns max_k E_k */
DEVFN double ebound(const double* L, double* E)
{
    int c, k;
    double m = 0.0;
    for (k = 0; k < 12; k++) {
        double s = 0.0;
        for (c = 0; c < 12; c++) s += L[c * 12 + k] * L[c * 12 + k];
        E[k] = sqrt(s) * (1.0 + 1e-15);
        if (E[k] > m) m = E[k];
    }
    return m;
}
/* Taylor remainder of one step, |Phi(x+e) - Phi(x) - J e| per component, for |e| <= r0 and
   all intermediate states bounded by X.  A kick adds h_k F(p, w), F cubic; for a monomial
   c a b c' with increments <= eps the remainder is <= |c| (3 X eps^2 + eps^3).  Increments of
   q = p - w are <= 2 r; the coefficient sum per force component is <= 8; sum |h_k| = 4.41 dt
   over the step; the error inside one step grows by < 1.1 (|J - 1| = O(dt) per stage).       */
DEVFN double taylor_rem(double X, double r0, double dt)
{
    double e = 2.0 * 1.1 * r0;
    return 1.1 * 4.41 * dt * 8.0 * (3.0 * X * e * e + e * e * e);
}

/* status of one member: centre c_k, error bound E_k  ->  bit0 possible, bit1 certain, bit2 lost */
/* GHOST_ESC: escape = ghost amplitude max(|w|, |W|) > R (components 3-5, 9-11); otherwise
   max over all 12 components, as in the production kernels.                              */
DEVFN int estatus(const double* c, const double* E, double R)
{
    double hi = 0.0, lo = 0.0;
    int k, st = 0;
    for (k = 0; k < 12; k++) {
        if (!(c[k] == c[k]) || !(E[k] == E[k])) st |= 4;
        if (E[k] > 0.5 * R) st |= 4;
        if (GHOST_ESC && !((k >= 3 && k < 6) || k >= 9)) continue;
        if (fabs(c[k]) + E[k] > hi) hi = fabs(c[k]) + E[k];
        if (fabs(c[k]) - E[k] > lo) lo = fabs(c[k]) - E[k];
    }
    if (hi > R) st |= 1;
    if (lo > R) st |= 2;
    return st;
}

/* S (24, N): reference state then difference, as hdym_diffpair.
   out rows: 0 ta_lo, 1 ta_hi, 2 tb_lo, 3 tb_hi  (hi: -1 not escaped by nmax, -2 enclosure lost)
             4 max error bound of a at its certain escape (absolute)
             then per checkpoint j: 5+2j  max error bound of a,  6+2j  of b  (-1: not alive)  */
DEVFN void bpair(const double* S, double* out, int N, int i, double dt, long long nmax, double R,
                 double omd, const long long* chk, int nchk)
{
    bd p[3], w[3], P[3], W[3];
    double Q[2][144], Ev[2][12], x[2][12], c[2][12], E[2][12], def[2][12], emax[2] = {0, 0};
    int k, j = 0, n, m, mode = 2;               /* 2: v = a, v+d = b;  0: only a in v;  1: only b in v */
    long long s;
    for (k = 0; k < 3; k++) {
        p[k] = mkb(mkd(S[(0 + k) * N + i], 0), 0.0, mkd(S[(12 + k) * N + i], 0), 0.0);
        w[k] = mkb(mkd(S[(3 + k) * N + i], 0), 0.0, mkd(S[(15 + k) * N + i], 0), 0.0);
        P[k] = mkb(mkd(S[(6 + k) * N + i], 0), 0.0, mkd(S[(18 + k) * N + i], 0), 0.0);
        W[k] = mkb(mkd(S[(9 + k) * N + i], 0), 0.0, mkd(S[(21 + k) * N + i], 0), 0.0);
    }
    for (m = 0; m < 2; m++) for (k = 0; k < 144; k++) Q[m][k] = 0.0;   /* exact initial data */
    for (k = 0; k < 5 + 2 * nchk; k++) out[k * N + i] = -1.0;
    for (s = 1; s <= nmax; s++) {
        bd* g[4] = {p, w, P, W};
        int alive[2] = {mode != 1, mode != 0};
        if (mode == 2) {
            int na = nsub_b(p, w, dt, omd, 0), nb = nsub_b(p, w, dt, omd, 1);
            n = na > nb ? na : nb;
        } else n = nsub_b(p, w, dt, omd, 0);
        /* base points of the pseudo-orbits at the start of the step */
        double X = 0.0;
        for (m = 0; m < 4; m++) for (k = 0; k < 3; k++) {
            double va = g[m][k].v.h, vb = ddadd(g[m][k].v, g[m][k].d).h;
            x[0][3 * m + k] = va; x[1][3 * m + k] = (mode == 2) ? vb : va;
            if (fabs(va) > X) X = fabs(va);
            if (fabs(vb) > X) X = fabs(vb);
            g[m][k].rv = 0.0; g[m][k].rd = 0.0;   /* one-step balls: defects */
        }
        { int u; double h = dt / n; for (u = 0; u < n; u++) yoshida_b(p, w, P, W, h); }
        for (m = 0; m < 4; m++) for (k = 0; k < 3; k++) {
            double cv = g[m][k].v.h, cb = ddadd(g[m][k].v, g[m][k].d).h;
            def[0][3 * m + k] = g[m][k].rv;
            def[1][3 * m + k] = (mode == 2) ? g[m][k].rv + g[m][k].rd : g[m][k].rv;
            c[0][3 * m + k] = cv; c[1][3 * m + k] = (mode == 2) ? cb : cv;
            if (fabs(cv) > X) X = fabs(cv);
            if (fabs(cb) > X) X = fabs(cb);
        }
        X = 1.5 * X + 1.0;
        int st[2] = {0, 0}, mem;
        for (mem = 0; mem < 2; mem++) {
            if (!alive[mem]) continue;
            int f = (mode == 2) ? mem : 0;       /* single member keeps its own frame */
            int fr = (mode == 2) ? mem : (mode == 0 ? 0 : 1);
            double shift = 0.0;                  /* base point = hi part, |true - hi| <= 2^-52 |hi| */
            for (k = 0; k < 12; k++) shift = fmax(shift, 2.0 * U2 * fabs(x[f][k]));
            /* tangent is evaluated at the hi part x~ of the base point: [J(x) - J(x~)] e is
               bounded by the cross term of the quadratic remainder, linear in e            */
            double rem = taylor_rem(X, emax[fr], dt)
                       + 1.1 * 4.41 * dt * 8.0 * 6.0 * X * (2.2 * emax[fr]) * (2.2 * shift);
            ell_update(Q[fr], x[f], dt, n, def[f], rem, 1);
            emax[fr] = ebound(Q[fr], E[fr]);
            for (k = 0; k < 12; k++) Ev[fr][k] = E[fr][k];
            for (k = 0; k < 12; k++) E[fr][k] += 2.0 * U2 * fabs(c[f][k]) * INFL;   /* centre = hi part */
            st[mem] = estatus(c[f], E[fr], R);
        }
        double t = s * dt;
        int done = 0;
        for (mem = 0; mem < 2; mem++) {
            if (!alive[mem]) continue;
            int lo = mem ? 2 : 0, hi = mem ? 3 : 1;
            if ((st[mem] & 1) && out[lo * N + i] < 0) out[lo * N + i] = t;
            if (st[mem] & 4) { out[hi * N + i] = -2.0; if (mode == 2 && mem == 0) out[3 * N + i] = -2.0; }
        }
        if ((st[0] & 4 && alive[0] && mode == 2) || (mode != 2 && (st[0] | st[1]) & 4)) break;
        if (mode == 2) {
            int ea = (st[0] & 2) != 0, eb = (st[1] & 2) != 0, lb = (st[1] & 4) != 0;
            if (ea) { out[1 * N + i] = t; out[4 * N + i] = emax[0]; }
            if (eb) out[3 * N + i] = t;
            int bdone = eb || lb;
            if (ea && bdone) break;
            if (ea) {                            /* continue b alone: v <- fl(v + d) */
                for (m = 0; m < 4; m++) for (k = 0; k < 3; k++) {
                    dd cc = ddadd(g[m][k].v, g[m][k].d);
                    g[m][k].v = cc; g[m][k].d = mkd(0, 0);
                    def[1][3 * m + k] = DDU * ddabs(cc);
                }
                ell_update(Q[1], x[1], dt, 1, def[1], 0.0, 0);   /* fold in the switch rounding */
                mode = 1;
            } else if (bdone) {
                for (m = 0; m < 4; m++) for (k = 0; k < 3; k++) g[m][k].d = mkd(0, 0);
                mode = 0;
            }
        } else {
            int mem1 = (mode == 0) ? 0 : 1;
            if (st[mem1] & 2) {
                out[(mem1 ? 3 : 1) * N + i] = t;
                if (!mem1) out[4 * N + i] = emax[0];
                done = 1;
            }
        }
        if (done) break;
        while (j < nchk && chk[j] <= s) {
            if (chk[j] == s && mode == 2) {
                out[(5 + 2 * j) * N + i] = emax[0];
                out[(6 + 2 * j) * N + i] = emax[1];
            }
            j++;
        }
    }
    /* validation rows: final centre of a (double-word, hi + lo) and its error bound */
    for (m = 0; m < 4; m++) for (k = 0; k < 3; k++) {
        out[(5 + 2 * nchk + 3 * m + k) * N + i] = (m == 0 ? p : m == 1 ? w : m == 2 ? P : W)[k].v.h;
        out[(17 + 2 * nchk + 3 * m + k) * N + i] = (m == 0 ? p : m == 1 ? w : m == 2 ? P : W)[k].v.l;
        out[(29 + 2 * nchk + 3 * m + k) * N + i] = Ev[0][3 * m + k];
    }
}
"""


GHOST = [1]


def _dll_search_path(cc):
    """Windows: Python >= 3.8 does not search PATH for a DLL's dependencies (libgomp-1.dll,
    libwinpthread-1.dll from MinGW), so add the compiler's bin directory explicitly."""
    if sys.platform == 'win32':
        import shutil
        exe = shutil.which(cc)
        if exe:
            try:
                os.add_dll_directory(os.path.dirname(os.path.abspath(exe)))
            except (AttributeError, OSError):
                pass


def c_src():
    return (f"#include <math.h>\n#define GHOST_ESC {GHOST[0]}\n#define DEVFN static inline\n"
            + BALL_SRC) + r"""
void bpairs(const double* S, double* out, int N, double dt, long long nmax, double R, double omd,
            const long long* chk, int nchk)
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
    for (i = 0; i < N; i++) bpair(S, out, N, i, dt, nmax, R, omd, chk, nchk);
}
"""


def cuda_src():
    return (f"#define GHOST_ESC {GHOST[0]}\n#define DEVFN __device__ __forceinline__\n"
            + BALL_SRC) + r"""
extern "C" __global__
void bpairs(const double* S, double* out, int N, double dt, long long nmax, double R,
            double omd, const long long* chk, int nchk)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) bpair(S, out, N, i, dt, nmax, R, omd, (const long long*)chk, nchk);
}
"""


class CBall:
    def __init__(self):
        d = tempfile.mkdtemp()
        src = os.path.join(d, 'hdym_ball.c')
        lib = os.path.join(d, 'hdym_ball' + ('.dll' if sys.platform == 'win32' else '.so'))
        with open(src, 'w') as fh:
            fh.write(c_src())
        cc = os.environ.get('CC', 'gcc')
        # -ffp-contract=off: an FMA rounds once, and the bounds assume each op rounds
        base = [cc, '-O3', '-march=native', '-ffp-contract=off', '-fno-fast-math', '-shared', '-fPIC', '-o', lib,
                src, '-lm']
        try:
            subprocess.check_call(base[:5] + ['-fopenmp'] + base[5:])
        except subprocess.CalledProcessError:
            subprocess.check_call(base)
        _dll_search_path(cc)
        self.lib = ctypes.CDLL(lib)
        dp, lp = ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_longlong)
        self.lib.bpairs.argtypes = [dp, dp, ctypes.c_int, ctypes.c_double, ctypes.c_longlong,
                                    ctypes.c_double, ctypes.c_double, lp, ctypes.c_int]
        print("C ball/difference-number kernel", flush=True)

    def run(self, S, dt, nmax, R, omd, chk):
        S = np.ascontiguousarray(S, dtype=np.float64)
        N = S.shape[1]
        chk = np.ascontiguousarray(chk, dtype=np.int64)
        out = np.empty((41 + 2 * len(chk), N))
        dp, lp = ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_longlong)
        self.lib.bpairs(S.ctypes.data_as(dp), out.ctypes.data_as(dp), N, dt, int(nmax), R, omd,
                        chk.ctypes.data_as(lp), len(chk))
        return out


class CupyBall:
    def __init__(self):
        import cupy as cp
        self.cp = cp
        self.kern = cp.RawKernel(cuda_src(), 'bpairs', options=('--fmad=false',))
        print("CuPy ball/difference-number kernel, --fmad=false", flush=True)

    def run(self, S, dt, nmax, R, omd, chk):
        cp = self.cp
        Sd = cp.asarray(np.ascontiguousarray(S, dtype=np.float64))
        N = Sd.shape[1]
        out = cp.empty((41 + 2 * len(chk), N), dtype=cp.float64)
        threads = 128
        self.kern(((N + threads - 1) // threads,), (threads,),
                  (Sd, out, np.int32(N), np.float64(dt), np.int64(nmax), np.float64(R),
                   np.float64(omd), cp.asarray(np.asarray(chk, dtype=np.int64)),
                   np.int32(len(chk))))
        cp.cuda.Stream.null.synchronize()
        return cp.asnumpy(out)


def classify(out, ok, th):
    """Return (cert_diff, cert_same, uncertified) boolean arrays."""
    alo, ahi, blo, bhi = out[0], out[1], out[2], out[3]
    both = ok & (ahi > 0) & (bhi > 0)
    diff = both & ((blo - ahi > th) | (alo - bhi > th))
    same = both & (np.maximum(ahi, bhi) - np.minimum(alo, blo) <= th)
    return diff, same, ok & ~diff & ~same


# ----------------------------------------------------------------- self-test
def _py_dd():
    """The kernel's double-word algorithms in Python (IEEE doubles, math.fma)."""
    import math

    def fast2sum(a, b):
        s = a + b
        return s, b - (s - a)

    def two_sum(a, b):
        s = a + b
        bb = s - a
        return s, (a - (s - bb)) + (b - bb)

    def add(x, y):
        sh, sl = two_sum(x[0], y[0])
        th, tl = two_sum(x[1], y[1])
        vh, vl = fast2sum(sh, sl + th)
        return fast2sum(vh, tl + vl)

    def mul(x, y):
        ch = x[0] * y[0]
        cl1 = _fma(x[0], y[0], -ch)
        tl1 = _fma(x[0], y[1], x[1] * y[1])
        cl2 = _fma(x[1], y[0], tl1)
        return fast2sum(ch, cl1 + cl2)
    return add, mul


def _mp_orbit(x0, tmax, dt=0.005, omd=0.25, dps=60):
    """Member a's orbit of the SAME map in 60-digit arithmetic: identical double
    coefficients, h = fl(dt/n), drift/kick constants fl(YC*h), fl(YD*h)."""
    import mpmath as mp
    mp.mp.dps = dps
    YC = (0.6756035959798288, -0.17560359597982886, -0.17560359597982886, 0.6756035959798288)
    YD = (1.3512071919596576, -1.7024143839193153, 1.3512071919596576)
    p = [mp.mpf(float(v)) for v in x0[0:3]]
    w = [mp.mpf(float(v)) for v in x0[3:6]]
    P = [mp.mpf(float(v)) for v in x0[6:9]]
    W = [mp.mpf(float(v)) for v in x0[9:12]]

    def kick(h):
        x, y, z = p[0] - w[0], p[1] - w[1], p[2] - w[2]
        x2, y2, z2 = x * x, y * y, z * z
        yz, zx, xy = y2 + z2, z2 + x2, x2 + y2
        g = (x * yz, y * zx, z * xy)
        txy, txz, tyz = 2 * x * y, 2 * x * z, 2 * y * z
        hv = (yz * w[0] + txy * w[1] + txz * w[2],
              txy * w[0] + zx * w[1] + tyz * w[2],
              txz * w[0] + tyz * w[1] + xy * w[2])
        hh = mp.mpf(h)
        for k in range(3):
            P[k] -= hh * (g[k] + hv[k])
            W[k] += hh * (hv[k] + w[k])

    def drift(ch):
        cc = mp.mpf(ch)
        for k in range(3):
            p[k] += cc * P[k]
            w[k] -= cc * W[k]
    nsteps = int(round(tmax / dt))
    for _ in range(nsteps):
        om = max(abs(float(p[k] - w[k])) for k in range(3))
        n = max(1, int(np.ceil(om * dt / omd)))
        h = dt / n
        for _u in range(n):
            drift(YC[0] * h); kick(YD[0] * h)
            drift(YC[1] * h); kick(YD[1] * h)
            drift(YC[2] * h); kick(YD[2] * h)
            drift(YC[3] * h)
    return p + w + P + W


def selftest(n_traj=3, tmax=40.0):
    from fractions import Fraction as Fr
    print("selftest 1: double-word add/mul error <= 2^-101 |result| (vs exact rationals)")
    add, mul = _py_dd()
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(20000):
        x = [float(v) for v in rng.standard_normal(1) * 10.0 ** rng.integers(-8, 8)]
        x = (x[0], x[0] * float(rng.standard_normal()) * 2 ** -60)
        y = [float(v) for v in rng.standard_normal(1) * 10.0 ** rng.integers(-8, 8)]
        y = (y[0], y[0] * float(rng.standard_normal()) * 2 ** -60)
        for op, ex in ((add, Fr(x[0]) + Fr(x[1]) + Fr(y[0]) + Fr(y[1])),
                       (mul, (Fr(x[0]) + Fr(x[1])) * (Fr(y[0]) + Fr(y[1])))):
            r = op(x, y)
            got = Fr(r[0]) + Fr(r[1])
            if ex != 0:
                worst = max(worst, float(abs(got - ex) / abs(got)) / 2 ** -101)
    print(f"   worst error / (2^-101 |result|) = {worst:.3f}  (must be < 1)")
    assert worst < 1.0

    print(f"selftest 2: an independent 60-digit orbit of the same map lies inside the "
          f"certified ellipsoid ({n_traj} trajectories to t = {tmax:g})")
    E, eps, width = 0.3, 0.01, em.default_width('ym', 0.3)
    Sa, D, ok = dpm.sample_dpairs('ym', E, eps, 1, 'xy', [0.0, 0.0], width, n_traj,
                                  1e-20 * width, hp.pair_rng(11, 1e-20))
    dt = 0.005
    out = CBall().run(np.concatenate([Sa, D], 0), dt, int(round(tmax / dt)), 1e4, 0.25, [])
    import mpmath as mp
    for i in range(n_traj):
        ch, cl, Eb = out[5:17, i], out[17:29, i], out[29:41, i]
        t0 = time.time()
        orb = _mp_orbit(Sa[:, i], tmax, dt)
        err = [abs(orb[k] - (mp.mpf(float(ch[k])) + mp.mpf(float(cl[k])))) for k in range(12)]
        ratio = max(float(err[k]) / Eb[k] for k in range(12))
        print(f"   traj {i}: max |true - centre| = {float(max(err)):.2e}, bound {Eb.max():.2e}, "
              f"worst ratio {ratio:.3f}  ({time.time() - t0:.0f}s)")
        assert ratio <= 1.0, "enclosure violated"
    print("selftest passed")


# ----------------------------------------------------------------- sweep
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--E', type=float, default=0.3)
    ap.add_argument('--eps', type=float, default=0.01)
    ap.add_argument('--section', choices=['ghost', 'ym'], default='ym')
    ap.add_argument('--plane', default='xy')
    ap.add_argument('--center', type=float, nargs=2, default=None)
    ap.add_argument('--width', type=float, default=None)
    ap.add_argument('--M', type=int, default=2000)
    ap.add_argument('--delta-max', type=float, default=1e-2)
    ap.add_argument('--delta-min', type=float, default=1e-40)
    ap.add_argument('--per-decade', type=float, default=0.5)
    ap.add_argument('--t-sep', type=float, nargs='+', default=[10.0, 50.0])
    ap.add_argument('--thresholds', type=float, nargs='+', default=[2.0, 5.0, 10.0])
    ap.add_argument('--dt', type=float, default=0.005)
    ap.add_argument('--tmax', type=float, default=3000.0)
    ap.add_argument('--esc', choices=['ghost', 'all'], default='ghost',
                    help="ghost: escape when max(|w|,|W|) > R (default R = 10); "
                         "all: max over all 12 components, as production (certifies poorly: "
                         "the final blow-up is where timing sensitivity diverges)")
    ap.add_argument('--R', type=float, default=None)
    ap.add_argument('--omega-dt-max', type=float, default=0.25)
    ap.add_argument('--backend', choices=['c', 'cupy'], default='c')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--pair-seed', type=int, default=12345)
    ap.add_argument('--out', default='results_precision')
    ap.add_argument('--tag', default=None)
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    GHOST[0] = 1 if a.esc == 'ghost' else 0
    if a.R is None:
        a.R = 10.0 if a.esc == 'ghost' else 1e4
    if a.center is None:
        a.center = [np.pi, np.pi] if a.section == 'ghost' else [0.0, 0.0]
    if a.width is None:
        a.width = em.default_width(a.section, a.E)
    tag = a.tag or f"ball_{a.esc}R{a.R:g}_{a.section}_E{a.E:g}_eps{a.eps:g}_M{a.M}"
    os.makedirs(a.out, exist_ok=True)
    be = CupyBall() if a.backend == 'cupy' else CBall()
    nd = int(round(a.per_decade * np.log10(a.delta_max / a.delta_min))) + 1
    rel = np.logspace(np.log10(a.delta_max), np.log10(a.delta_min), nd)
    chk = [int(round(t / a.dt)) for t in a.t_sep]
    nmax = int(np.ceil(a.tmax / a.dt))
    nT = len(a.thresholds)
    f_lo, f_hi, f_c = np.full((nT, nd), np.nan), np.full((nT, nd), np.nan), np.full((nT, nd), np.nan)
    unc_frac = np.full((nT, nd), np.nan)
    lost = np.zeros(nd)
    rv_esc = np.full(nd, np.nan)
    rv_chk = np.full((len(chk), nd), np.nan)
    rd_chk = np.full((len(chk), nd), np.nan)
    iw = np.full(nd, np.nan)
    t0 = time.time()
    for i, r in enumerate(rel):
        d = r * a.width
        Sa, D, ok = dpm.sample_dpairs(a.section, a.E, a.eps, a.seed, a.plane, a.center, a.width,
                                      a.M, d, hp.pair_rng(a.pair_seed, r))
        out = be.run(np.concatenate([Sa, D], 0), a.dt, nmax, a.R, a.omega_dt_max, chk)
        n_ok = ok.sum()
        lost[i] = np.mean(((out[1] == -2) | (out[3] == -2))[ok])
        for k, th in enumerate(a.thresholds):
            cd, cs, un = classify(out, ok, th)
            f_lo[k, i] = cd.sum() / n_ok
            f_hi[k, i] = 1.0 - cs.sum() / n_ok
            unc_frac[k, i] = un.sum() / n_ok
            mid = ok & (out[1] > 0) & (out[3] > 0)
            ca, cb = 0.5 * (out[0] + out[1]), 0.5 * (out[2] + out[3])
            f_c[k, i] = np.mean((np.abs(ca - cb) > th)[mid]) if mid.any() else np.nan
        e = out[4][ok & (out[4] >= 0)]
        rv_esc[i] = np.median(e) if e.size else np.nan
        for j in range(len(chk)):
            x = out[5 + 2 * j][ok & (out[5 + 2 * j] >= 0)]
            y = out[6 + 2 * j][ok & (out[6 + 2 * j] >= 0)]
            rv_chk[j, i] = np.median(x) if x.size else np.nan
            rd_chk[j, i] = np.median(y) if y.size else np.nan
        wa = (out[1] - out[0])[ok & (out[1] > 0)]
        iw[i] = np.median(wa) if wa.size else np.nan
        k = nT - 1
        print(f"  delta/w={r:.1e}  f(>{a.thresholds[k]:g}) in [{f_lo[k, i]:.4f}, {f_hi[k, i]:.4f}]"
              f" (centres {f_c[k, i]:.4f}, uncertified {unc_frac[k, i]:.4f})  lost {lost[i]:.3f}  "
              f"bound(a) t={a.t_sep[0]:g}: {rv_chk[0, i]:.1e} at esc: {rv_esc[i]:.1e}  "
              f"({time.time() - t0:.0f}s)", flush=True)
    print("\nrigorous bracket on the uncertainty exponent (fit to f_lo and f_hi separately):")
    al = []
    for k, th in enumerate(a.thresholds):
        alo = unc.fit_alpha(rel, f_lo[k])[:2]
        ahi = unc.fit_alpha(rel, f_hi[k])[:2]
        al.append((alo, ahi))
        print(f"  |dt|>{th:g}: alpha from f_lo {alo[0]:.3f}+-{alo[1]:.3f}, "
              f"from f_hi {ahi[0]:.3f}+-{ahi[1]:.3f}; max uncertified fraction "
              f"{np.nanmax(unc_frac[k]):.4f}")
    np.savez(os.path.join(a.out, tag + '.npz'), rel=rel, f_lo=f_lo, f_hi=f_hi, f_c=f_c,
             unc=unc_frac, lost=lost, rv_esc=rv_esc, rv_chk=rv_chk, rd_chk=rd_chk, iw=iw,
             args=json.dumps(vars(a)))
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.5))
        ax = axs[0]
        k = nT - 1
        ax.fill_between(rel, np.maximum(f_lo[k], 1e-4), f_hi[k], color='C0', alpha=0.25,
                        label='certified bracket [f_lo, f_hi]')
        ax.loglog(rel, f_lo[k], '-', color='C0', lw=1)
        ax.loglog(rel, f_hi[k], '-', color='C0', lw=1)
        ax.loglog(rel, f_c[k], 'o', ms=3, color='k', label='ball centres (float result)')
        ax.axvline(2e-15, color='0.5', ls=':', lw=1, label='fp64 floor, separate pairs')
        ax.set_xlabel('δ / section width'); ax.set_ylabel(f'f(δ), |Δt_esc| > {a.thresholds[k]:g}')
        ax.set_title('certified uncertainty fraction', fontsize=10)
        ax.legend(fontsize=7)
        ax = axs[1]
        unc_k = unc_frac[nT - 1]
        ax.semilogx(rel, 1 - unc_k, 'o-', ms=3, label='certified (different or same)')
        ax.semilogx(rel, lost, 's-', ms=3, label='enclosure lost (long-lived orbits)')
        ax.set_ylim(0, 1)
        ax.set_xlabel('δ / section width'); ax.set_ylabel('fraction of pairs')
        ax.set_title('how many pairs the guarantees cover', fontsize=10)
        ax.legend(fontsize=7)
        plt.suptitle(f"certified pairs (double-word ellipsoids on difference numbers), escape "
                     f"{'|ghost|' if a.esc == 'ghost' else 'max'} > {a.R:g}: E={a.E:g} eps={a.eps:g} "
                     f"{a.section} section, M={a.M}", fontsize=10)
        plt.tight_layout()
        png = os.path.join(a.out, tag + '.png')
        plt.savefig(png, dpi=150)
        print("wrote", png)
    except Exception as e:  # noqa
        print("plot skipped:", e)


if __name__ == '__main__':
    main()
