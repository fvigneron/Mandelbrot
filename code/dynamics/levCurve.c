//
//  levCurve.c
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the reference below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
//
//  Copyright 2019 - 2021 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the Lesser GNU Public Licence v3.0
//

#include <mpfr.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpc.h"
#include "levCurve.h"
#include "mandel.h"
#include "memFile.h"

// MARK: Creation, computation and destruction

/// @brief Creates an empty levelCurve of period @c per, angles of @c 2^{-tpow},
/// precision @c prec, start angle @c sta, end angle @c ena and radius @c r.
///
/// @param per the period of the level curve
/// @param tpow the power of two to obtain the smalles dyadic angle
/// @param sta the start angle
/// @param ena the end angle
/// @param prec the precision
/// @param r the radius
///
/// @return the level curve
static levc levc_alloc(uint per, int tpow, uint prec, ulong sta, ulong ena, mpfr_t r) {
    if(tpow == 0 || tpow > LEVC_MAX_2POW || prec < LEVC_MIN_PREC || r == NULL
       || ! mpfr_regular_p(r) || mpfr_cmp_si(r, per == 1 ? 0 : 2) <= 0
       || ena < sta || ena >= (1L << tpow)) {
        return NULL;
    }
    
    levc lc = malloc(sizeof(levelCurve));
    lc->period = per;
    lc->ang2pow = tpow;
    lc->prec = prec;
    lc->startAngle = sta;
    lc->endAngle = ena;
    
    mpfr_init2(lc->radius, prec);
    mpfr_set(lc->radius, r, MPFR_RNDN);
    
    mpfr_init2(lc->eps, LEVC_PREC_EPS);
    mpfr_set_zero(lc->eps, 1);
    
    lc->guard = LEVC_MIN_GUARD;
    
    mpv_init(lc->points, prec, 2 * levc_count(lc));
    
    return lc;
}

/// @brief Creates an empty levelCurve representing the upper half circle of period @c per, angles of @c 2^{-tpow},
/// precision @c prec and radius @c r.
///
/// @param per the period of the level curve
/// @param tpow the power of two to obtain the smalles dyadic angle
/// @param prec the precision
/// @param r the radius
///
/// @return the level curve
static levc levc_alloc_upper(uint per, int tpow, uint prec, mpfr_t r) {
    return levc_alloc(per, tpow, prec, 0, 1L << (tpow - 1), r);
}

/// @brief Computes in place the level curve which has the endpoints initialized , using  fp80.h.
///
/// Checks if the previous points are computed with max error @c 1E-17 and if all convergences
/// are not worse than @c maxErr.
///
/// @param lc the level curve to compute
/// @param r the radius or level of the curve
/// @param maxErr the upper bound for the error
///
/// @return @ref true if the computation succeeded, @ref false otherwise
static bool levc_fill_inl(levc lc, ldbl r, ldbl maxErr) {
    // compute how many targets are needed w.r.t. lc->radius and lc->guard
    ldbl a = 2 * PI * lc->guard * r / (r - 2);
    int tpowc = 4; // two power of the circle of targets
    while((1L << tpowc) <= a && tpowc < 30) {
        tpowc ++;
    }
    
    if(tpowc == 30) {
        return false;
    }
    
    int per = lc->period;
    int a2p = lc->ang2pow;
    int lctp = a2p - per + 1;
    tpowc = lctp > tpowc ? lctp : tpowc;     // if points on the curve are even denser
    
    int tpow = tpowc + per - 1;        // the power of two for the angles, >= lc->ang2pow
    // for any angle in the scale tpow, the mask gives the index of the target
    ulong tpowMask = (1L << tpowc) - 1;
    ldbl mula = ldexpl(PI, -tpowc + 1); // coefficient to get the actual angle
    // for every point on the level curve, there will be 2^extp points computed
    int extp = tpow - a2p;
    ulong extpMask = (1L << extp) - 1; // mask to check when to store a point
    
    // compute intermediary points
    ulong sta = lc->startAngle << extp;
    ulong ena = lc->endAngle << extp;
    
    bool ok = true;
    fp80 c, t;
    ulong am, i;
    ldbl err, eps = 0;
    
    mpv pts = lc->points;
    mpv_getcl(c, pts, 0); // start point on the right
    for (ulong ula = sta; ula <= ena && ok; ula++) {
        am = ula & tpowMask;
        a = am * mula;
        
        t->x = r * cosl(a);
        t->y = r * sinl(a);
        
        ok = mandel_sol_refl(c, &err, c, t, per, LEVC_MAX_ITER, 1E-17, 2);
        ok = ok && err <= maxErr;
        
        if((ula & extpMask) == 0 && ok) { // store the new value
            i = (ula - sta) >> extp;
            
            if(ula == sta || ula == ena) {
                mpv_getcl(t, pts, i);
                
                ok = fp80_dist(c, t) <= 1E-17;
            }
            
            mpv_setcl(pts, i, c);
            eps = err > eps ? err : eps;
        }
    }
    
    if(ok) {
        mpfr_set_ld(lc->eps, 2 * eps, MPFR_RNDU); // allow for rounding errors
    }
    
    return ok;
}

static bool levc_refine_and_check(levc lc, ulong ind, mpc t, mpfr_t eps, mpc c,
                                   mpc nt, mpfr_t m, bool check);

/// @brief Computes the next point @c ep from the starting point @c sp with target @c t and period @c per.
///
/// It uses fp80 as long as the modulus of the Newton term is amller than LEVC_MAX_FP80_STEP, then refines the
/// result up to the precision of @c ep.
///
/// @param ep the result
/// @param sp the start point
/// @param t the target
/// @param per the period
///
/// @return @ref true if the computation succeeded, @ref false otherwise
static bool levc_next(mpc ep, mpc sp, mpc t, int per, uint prec, mpfr_t err) {
    fp80 c8, t8, nt8, v8, d8;
    mpc_get80(c8, sp);
    mpc_get80(t8, t);
    
    if(prec == LEVC_MIN_PREC) { // try all in fp80
        if(mandel_sol_refl(c8, NULL, c8, t8, per, LEVC_MAX_ITER, 1E-17, 2)) {
            mpc_set80(ep, c8);
            
            return true;
        } else { // no point in trying with mp if the precision is minimal
            return false;
        }
    }
    
    mandel_val_derl(v8, d8, c8, per);
    fp80_sub(v8, v8, t8);
    ldbl md = fp80_mod(d8);
    ldbl mnt = fp80_mod(v8) / md;
    
    if(mnt < LEVC_MIN_FP80_STEP || md * LEVC_MIN_FP80_STEP > 1) {
        return mandel_sol_ref(ep, err, sp, t, per, LEVC_MAX_ITER, 2);
    }
    
    int iter = 0;
    do {
        fp80_div(nt8, v8, d8);
        fp80_sub(c8, c8, nt8);
        iter ++;
        
        if(iter < LEVC_MAX_ITER) {
            mandel_val_derl(v8, d8, c8, per);
            fp80_sub(v8, v8, t8);
            md = fp80_mod(d8);
            mnt = fp80_mod(v8) / md;
        } else {
            return false; // no more iterates to refine
        }
    } while(mnt >= LEVC_MIN_FP80_STEP && md * LEVC_MIN_FP80_STEP <= 1);
    
    mpc_set80(sp, c8);
    
    return mandel_sol_ref(ep, err, sp, t, per, LEVC_MAX_ITER - iter, 2);
}

/// @brief Sets the value of the point with the given index @c ind in the level curve @c lc to @c p and updates
/// the error of the points in the curve to @c err.
///
/// @param lc the level curve
/// @param ind the index of the point
/// @param p the point
/// @param err the error of p
static void levc_set_point(levc lc, ulong ind, mpc p, mpfr_t err) {
    mpv_setc(lc->points, ind, p);
    
    if(mpfr_cmp(err, lc->eps) > 0) {
        mpfr_set(lc->eps, err, MPFR_RNDU);
    }
}

/// @brief Computes in place the level curve which has endpoints initialized with max error
/// @c eps, using mpc.h or fp80.h if the precision requested is LEVC_MIN_PREC.
///
/// @param lc the level curve to compute
/// @param eps the max error of the existing endpoints
///
/// @return @ref true if the computation succeeded and the endpoints error was bounded cy @c eps, @ref false otherwise
static bool levc_fill_in(levc lc) {
    if(levc_segs(lc) == 1) { // only two points, no fill in
        return true;
    }
    
    ldbl r = mpfr_get_ld(lc->radius, MPFR_RNDN);
    
    if(lc->prec == LEVC_MIN_PREC) { // low prec requested
        return levc_fill_inl(lc, r, mpfr_get_ld(lc->eps, MPFR_RNDD));
    }
    
    // compute how many targets are needed w.r.t. lc->radius and lc->guard
    ldbl ra = 2 * PI * lc->guard * r / (r - 2);
    int tpowc = 4; // two power of the circle of targets
    while((1L << tpowc) <= ra && tpowc < 30) {
        tpowc ++;
    }
    
    if(tpowc == 30) {
        return false;
    }
    
    int per = lc->period;
    int a2p = lc->ang2pow;
    int lctp = a2p - per + 1;
    tpowc = lctp > tpowc ? lctp : tpowc;     // if points on the curve are even denser
    
    int tpow = tpowc + per - 1;        // the power of two for the angles, >= lc->ang2pow
    // for any angle in the scale tpow, the mask gives the index of the target
    ulong tpowMask = (1L << tpowc) - 1;
    // for every point on the level curve, there will be 2^extp points computed
    int extp = tpow - a2p;
    ulong extpMask = (1L << extp) - 1; // mask to check when to store a point
    
    // compute intermediary points
    ulong sta = lc->startAngle << extp;
    ulong ena = lc->endAngle << extp;
    
    bool ok = true;
    
    uint prec = lc->prec;
    levc tc = levc_circle(tpowc, prec, lc->radius);
    mpv ts = tc->points;
    
    mpc c, t;
    mpc_init(c, prec);
    mpc_init(t, prec);
    
    mpfr_t err, d;
    mpfr_init2(err, prec);
    mpfr_init2(d, prec);
    
    mpv pts = lc->points;
    mpv_getc(c, pts, 0); // start point on the right
    
    // use the targets to compute intermediary points, store the nodes of lc
    for (ulong a = sta; a <= ena && ok; a++) {
        mpv_getc(t, ts, a & tpowMask);
        
        ok = levc_next(c, c, t, per, prec, err);
        
        if((a & extpMask) == 0 && ok) { // store the new value
            if(a == sta || a == ena) {
                mpv_getc(t, pts, (a - sta) >> extp);
                mpc_dist(d, c, t);
                
                ok = mpfr_cmp(d, lc->eps) <= 0;
            }
            
            levc_set_point(lc, (a - sta) >> extp, c, err);
        }
    }
    
    mpfr_clear(err);
    mpfr_clear(d);
    
    mpc_clear(t);
    mpc_clear(c);
    
    levc_free(tc);
    
    return ok;
}

static inline ldbl psi_rpnl(int n, ldbl c) {
    ldbl v = c;
    for (int i = 1; i < n; i++) {
        v *= v;
        v += c;
    }
    
    return v;
}

static bool psi_n_real_l(int n, fp80 c, ldbl r) {
    if(n < 1 || c == NULL || (r <= 2 && r >= -2)) {
        return false;
    }
    
    if(n == 1) {
        c->x = r;
        c->y = 0;
        
        return true;
    }
    
    ldbl x, y, m;
    if(r > 0) {
        x = 0.25;
        y = psi_rpnl(n, 1) >= r ? 1 : powl(r, ldexpl(1, 1 - n));
        m = (x + y) / 2;
        
        ldbl om;
        do {
            om = m;
            
            if(psi_rpnl(n, m) > r) {
                y = m;
            } else {
                x = m;
            }
            
            m = (x + y) / 2;
        } while(om != m);
    } else {
        ldbl l = -r;
        x = -2 - ldexpl(l, 2 * (1 - n));
        y = -2;
        m = (x + y) / 2;
        
        ldbl om;
        do {
            om = m;
            
            if(psi_rpnl(n, m) > l) {
                x = m;
            } else {
                y = m;
            }
            
            m = (x + y) / 2;
        } while(om != m);
    }
    
    c->x = m;
    c->y = 0;
    
    return true;
}

static inline void psi_rpn(long n, mpfr_t v, mpfr_t c) {
    mpfr_set(v, c, MPFR_RNDN);
    for (long i = 1; i < n; i++) {
        mpfr_sqr(v, v, MPFR_RNDN);
        mpfr_add(v, v, c, MPFR_RNDN);
    }
}

static bool psi_n_real(long n, mpc c, mpfr_t r) {
    if(n < 1 || c == NULL || r == NULL || (mpfr_cmp_si(r, 2) <= 0 && mpfr_cmp_si(r, -2) >= 0)) {
        return false;
    }
    
    if(n == 1) {
        mpc_setr(c, r);
        
        return true;
    }
    
    mpfr_t x, y, m, om, b;
    long prec = mpc_prec(c);
    mpfr_init2(x, prec);
    mpfr_init2(y, prec);
    mpfr_init2(m, prec);
    mpfr_init2(om, prec);
    mpfr_init2(b, prec);
    
    if(mpfr_cmp_si(r, 0) > 0) {
        mpfr_set_d(x, 0.25, MPFR_RNDN);
        
        mpfr_set_d(m, 1, MPFR_RNDN);
        psi_rpn(n, b, m);
        if(mpfr_cmp(b, r) >= 0) {
            mpfr_set_d(y, 1, MPFR_RNDN);
        } else {
            mpfr_div_2ui(m, m, n - 1, MPFR_RNDN);
            mpfr_pow(y, r, m, MPFR_RNDN);
        }
        
        mpfr_add(m, x, y, MPFR_RNDN);
        mpfr_div_2ui(m, m, 1, MPFR_RNDN);
        
        do {
            mpfr_set(om, m, MPFR_RNDN);
            
            psi_rpn(n, b, m);
            if(mpfr_cmp(b, r) > 0) {
                mpfr_set(y, m, MPFR_RNDN);
            } else {
                mpfr_set(x, m, MPFR_RNDN);
            }
            
            mpfr_add(m, x, y, MPFR_RNDN);
            mpfr_div_2ui(m, m, 1, MPFR_RNDN);
        } while(mpfr_cmp(om, m) != 0);
    } else {
        mpfr_t l;
        mpfr_init2(l, prec);
        mpfr_neg(l, r, MPFR_RNDN);
        
        mpfr_set_si(y, -2, MPFR_RNDN);
        mpfr_div_2ui(m, l, 2 * (n - 1), MPFR_RNDN);
        mpfr_sub(x, y, m, MPFR_RNDN);
        
        mpfr_add(m, x, y, MPFR_RNDN);
        mpfr_div_2ui(m, m, 1, MPFR_RNDN);
        
        do {
            mpfr_set(om, m, MPFR_RNDN);
            
            psi_rpn(n, b, m);
            if(mpfr_cmp(b, l) > 0) {
                mpfr_set(x, m, MPFR_RNDN);
            } else {
                mpfr_set(y, m, MPFR_RNDN);
            }
            
            mpfr_add(m, x, y, MPFR_RNDN);
            mpfr_div_2ui(m, m, 1, MPFR_RNDN);
        } while(mpfr_cmp(om, m) != 0);
    }
    
    mpc_setr(c, m);
    
    mpfr_clear(om);
    mpfr_clear(m);
    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(b);
    
    return true;
}


/// @brief Computes in place the level curve of the upper half circle using mpc.h.
///
/// @param lc the level curve to compute
static bool levc_compute_direct(levc lc) {
    int per = lc->period;
    bool ok = true;
    
    // init start and end points
    if(lc->prec == LEVC_MIN_PREC) {
        ldbl ep = 1E-17;
        
        ldbl r = mpfr_get_ld(lc->radius, MPFR_RNDN);
        
        fp80 v;
        ok = psi_n_real_l(per, v, r);
        mpv_setcl(lc->points, 0, v);
        
        ok = ok && psi_n_real_l(per, v, -r);
        mpv_setcl(lc->points, levc_segs(lc), v);
        
        mpfr_set_ld(lc->eps, ep, MPFR_RNDU);
    } else {
        long prec = lc->prec;
        mpfr_set_si_2exp(lc->eps, 1, LEVC_EXTRA_PREC - prec, MPFR_RNDN);
        
        mpc v;
        mpc_init(v, prec);
        
        ok = psi_n_real(per, v, lc->radius);
        mpv_setc(lc->points, 0, v);
        
        mpfr_neg(lc->radius, lc->radius, MPFR_RNDN);
        ok = ok && psi_n_real(per, v, lc->radius);
        mpv_setc(lc->points, levc_segs(lc), v);
        mpfr_neg(lc->radius, lc->radius, MPFR_RNDN);
        
        mpc_clear(v);
    }
        
    ok = ok && levc_fill_in(lc);
    
    return ok;
}

levc levc_new(uint per, uint tpow, uint prec, mpfr_t r, double guard) {
    if(guard < LEVC_MIN_GUARD || (prec == LEVC_MIN_PREC && per > LEVC_MAX_FP80_PREC) ||
       per > LEVC_MAX_PERIOD || tpow > LEVC_MAX_2POW) {
        return NULL;
    }
    
    if(per == 1) {
        return levc_circle(tpow, prec, r);
    }
    
    levc lc = levc_alloc_upper(per, tpow, prec, r);
    if(lc == NULL) {
        return NULL;
    }
    
    lc->guard = guard;
    
    if(levc_compute_direct(lc)) {
        return lc;
    }
    
    levc_free(lc);
    
    return NULL;
}

levc levc_circle(int tpow, uint prec, mpfr_t r) {
    levc c = levc_alloc(1, tpow, prec, 0, (1L << tpow) - 1, r);
    if(c == NULL) {
        return NULL;
    }
    
    mpfr_mul_2si(c->eps, c->radius, -prec, MPFR_RNDU);
    
    mpc u;
    mpc_init(u, prec);
    mpfr_t th;
    mpfr_init2(th, prec);
    for (ulong i = 0; i <= c->endAngle; i++) {
        mpfr_set_ui_2exp(th, i, -tpow, MPFR_RNDN);
        mpc_exp_2Pi_i(u, th);
        mpc_mulr(u, u, r);
        
        mpv_setc(c->points, i, u);
    }
    
    mpfr_clear(th);
    mpc_clear(u);
    
    return c;
}

void levc_free(levc lc) {
    if(lc == NULL || lc->points->vals == NULL || lc->points->sexp == NULL) {
        return;   // probably already cleared
    }
    
    mpv_clear(lc->points);
    mpfr_clear(lc->eps);
    mpfr_clear(lc->radius);
    
    free(lc);
}

// MARK: Pointwise operations

bool levc_angle(mpfr_t angle, levc lc, ulong ind) {
    if(angle == NULL || lc == NULL || ind >= levc_count(lc)) {
        return false;
    }
    
    mpfr_set_ui(angle, lc->startAngle + ind, MPFR_RNDN);
    mpfr_mul_2si(angle, angle, -((long) lc->ang2pow), MPFR_RNDN);
    
    return true;
}

ldbl levc_anglel(levc lc, ulong ind) {
    if(lc == NULL || ind >= levc_count(lc)) {
        return false;
    }
    
    ldbl angle = lc->startAngle + ind;
    angle = ldexpl(angle, -((int) lc->ang2pow));
    
    return angle;
}

bool levc_point(mpc p, levc lc, ulong ind) {
    if(p == NULL || lc == NULL || ind >= levc_count(lc)) {
        return false;
    }
    
    return mpv_getc(p, lc->points, ind);
}

bool levc_pointl(fp80 p, levc lc, ulong ind) {
    if(p == NULL || lc == NULL || ind >= levc_count(lc)) {
        return false;
    }
    
    return mpv_getcl(p, lc->points, ind);
}

// MARK: Refinements and derivations

/// @brief Return the base two logarithm of @c n if it is a positive power of @c 2, @c -1 otherwise.
///
/// @param n the number
///
/// @return @c k if @c n==2^k and @c k>=0, @c -1 otherwise
static inline int levc_log2(ulong n) {
    if(n == 0) {
        return -1;
    }
    
    int c = 0;
    while((n & 1) == 0) {
        n >>= 1;
        c ++;
    }
    
    return n == 1 ? c : -1;
}

/// @brief Creates a new curve with the same structure as @c src with the given precision @c prec.
///
/// @param src the source curve
/// @param prec the precision of the new curve
/// @param initPoints @ref true to initialize the vector of points, @ref false otherwise
/// @param copyPoints @ref true to copy the vector of points, @ref false otherwise
///
/// @return the new curve
static levc levc_copy(levc src, uint prec, bool initPoints, bool copyPoints) {
    levc lc = malloc(sizeof(levelCurve));
    
    lc->period = src->period;
    lc->ang2pow = src->ang2pow;
    lc->prec = prec;
    lc->startAngle = src->startAngle;
    lc->endAngle = src->endAngle;
    
    mpfr_init2(lc->radius, lc->prec);
    mpfr_set(lc->radius, src->radius, MPFR_RNDN);
    
    mpfr_init2(lc->eps, LEVC_PREC_EPS);
    mpfr_set(lc->eps, src->eps, MPFR_RNDN);
    
    lc->guard = src->guard;
    
    if(initPoints) {
        mpv_init(lc->points, lc->prec, src->points->count);
        
        if(copyPoints) {
            mpv_copy(lc->points, 0, src->points);
        }
    }
    
    return lc;
}

levc levc_sub_curve(levc src, ulong start, ulong step, ulong count) {
    int l2 = levc_log2(step);
    
    if(src == NULL || l2 < 0 || start % step != 0) {
        return NULL;
    }
    
    mpv pts = mpv_sub_vectorc(src->points, start, step, count);
    if(pts == NULL) {
        return NULL;
    }
    
    levc lc = levc_copy(src, src->prec, false, false);
    
    lc->ang2pow = src->ang2pow - l2;
    lc->startAngle = (src->startAngle + start) >> l2;
    lc->endAngle = lc->startAngle + count - 1;
    
    *lc->points = *pts;
    free(pts);              // ATTENTION: not mpv_free() !
    
    return lc;
}

/// @brief Refines the value of the point with index @c ind in the level curve @c lc. Checks if the error is at most @c eps.
///
/// @param lc the level curve
/// @param ind the index in @c lc
/// @param t the target, the value of @c p_n in the new found point
/// @param eps the max error of the existing result
/// @param c buffer for computations
/// @param nt buffer for computations
/// @param m buffer for computations
/// @param check @ref true to check the error, @ref false otherwise
///
/// @return @ref true if the computation succeeded and the error was bounded cy @c eps, @ref false otherwise
static bool levc_refine_and_check(levc lc, ulong ind, mpc t, mpfr_t eps, mpc c, mpc nt, mpfr_t m, bool check) {
    mpv_getc(c, lc->points, ind);
    
    bool conv = 0, div = 0, per = lc->period, prec = lc->prec;
    long exp;
    for (int j = 0; j < LEVC_MAX_ITER && ! conv && ! div; j++) {
        mandel_nt_sol(nt, c, t, per);
        mpc_sub(c, c, nt);
        
        exp = mpc_2exp(nt);
        div = exp >= 0;
        conv = -exp > prec * 3 / 4;
    }
    
    if(! conv || div) {
        return false;
    }
    
    mandel_nt_sol(nt, c, t, per);
    mpc_sub(c, c, nt);
        
    // update max error in the curve
    mpc_mod(m, nt);
    
    if(! check) {
        levc_set_point(lc, ind, c, m); // the refined value
        
        return true;
    }
    
    // check if error from previous point is less than eps
    mpv_getc(nt, lc->points, ind);
    mpc_sub(nt, nt, c);
    mpc_mod(m, nt);
    
    if(mpfr_cmp(m, eps) <= 0) {
        levc_set_point(lc, ind, c, m); // the refined value
        
        return true;
    }
    
    return false;
}

/// @brief Checks if the existing points are precise up to error @c eps and refines their precision to @c prec.
///
/// It performs @b approximative computations.
///
/// It is assumed that @c prec is larger than the previous precision of the points, thus it cannot be the minimal precision
/// and so all computations are performed with tools from mpc.h.
///
/// @param lc the level curve to refine in place
/// @param eps the max error of the existing points
///
/// @return @ref true if successfull, @ref false otherwise
static bool levc_refine_prec_direct(levc lc, mpfr_t eps) {
    mpfr_set_zero(lc->eps, 1); // will be recomputed for each point
    uint prec = lc->prec;
    uint tpow = lc->ang2pow;
    mpfr_ptr r = lc->radius;
    ulong count = levc_count(lc);
    
    mpc c, t, nt;
    mpc_init(c, prec);
    mpc_init(t, prec);
    mpc_init(nt, prec);
    bool ok = true;
    
    mpfr_t m;
    mpfr_init2(m, LEVC_PREC_EPS);
    if(lc->period > tpow) { // all values of p_n will be r
        mpc_setr(t, r);
        
        for (ulong i = 0; i < count && ok; i++) {
            ok = levc_refine_and_check(lc, i, t, eps, c, nt, m, true);
        }
    } else { // need a circle for targets, faster than recomputing each time
        int lcpm1 = lc->period - 1;
        levc cr = levc_circle(tpow - lcpm1, prec, r);
        
        for (ulong i = 0; i < count && ok; i++) {
            mpv_getc(t, cr->points, i >> lcpm1);
            
            ok = levc_refine_and_check(lc, i, t, eps, c, nt, m, true);
        }
        
        levc_free(cr);
    }
    
    mpfr_clear(m);
    
    mpc_clear(nt);
    mpc_clear(t);
    mpc_clear(c);
    
    return ok;
}

levc levc_refine(levc src, uint tpow, uint prec) {
    if(src == NULL || tpow > LEVC_MAX_2POW || tpow == 0 || prec < LEVC_MIN_PREC || src->period < 2) {
        return NULL;
    }
    
    if(tpow == src->ang2pow) {
        levc lc = levc_copy(src, prec, true, true);
        
        if(prec == src->prec) {
            return lc;
        } else if(prec < src->prec) {
            mpfr_t u;
            mpfr_init2(u, prec);
            
            if(mpv_ulp(u, lc->points)) {
                mpfr_max(lc->eps, lc->eps, u, MPFR_RNDU);
            }
            
            mpfr_clear(u);
            
            return lc;
        }
        
        // more precision is needed
        mpfr_div_2si(lc->eps, src->eps, prec - src->prec, MPFR_RNDN);
        if(levc_refine_prec_direct(lc, src->eps)) {
            return lc;
        }
        
        levc_free(lc);
        
        return NULL;
    } else if(tpow < src->ang2pow) {
        int l2 = src->ang2pow - tpow;
        levc lc = levc_sub_curve(src, 0, 1L << l2, (levc_segs(src) >> l2) + 1);
        
        if(prec == lc->prec) {
            return lc;
        }
        
        levc nc = levc_refine(lc, lc->ang2pow, prec);
        levc_free(lc);
        
        return nc;
    }
        
    // more points are needed
    int l2 = tpow - src->ang2pow;
    ulong sta = src->startAngle << l2;
    ulong ena = src->endAngle << l2;
    levc lc = levc_alloc(src->period, tpow, prec, sta, ena, src->radius);
    mpfr_div_2si(lc->eps, src->eps, prec - src->prec, MPFR_RNDN);
    
    mpv_get_setc(lc->points, 0, src->points, 0);
    mpv_get_setc(lc->points, levc_segs(lc), src->points, levc_segs(src));
    
    if(levc_fill_in(lc)) {
        return lc;
    }
    
    levc_free(lc);
    
    return NULL;
}

/// @brief Replaces the point fo index @c ind of the curve @c lc (whic has target @c pt of modulus @c pr) with the point
/// computed with target @c t.
///
/// If the @c lc->guard condition is not satisfied, it divides accordingly the segment @c [pt, t] and computes in several steps.
///
/// @warning All parameters are potentially used as buffers, the user should safely keep a copy of their values.
///
/// @param lc the level curve
/// @param ind the index in @c lc
/// @param pt previous target
/// @param pr previous radius, or @c |pt|
/// @param t the new target
/// @param c buffer
/// @param bt buffer
/// @param m buffer
/// @param b1 buffer
///
/// @return @ref true if successfull, @ref false otherwise
static bool levc_move_point(levc lc, ulong ind, mpc pt, mpfr_t pr, mpc t, mpc c,
                              mpc bt, mpfr_t m, mpfr_t b1, mpfr_t b2, bool epsUpdate) {
    mpv_getc(c, lc->points, ind);
    uint prec = lc->prec;
    
    do {
        mpc_sub(bt, t, pt);                       // bt = t - pt
        mpc_mod(m, bt);                           // m = |bt| = |t - pt|
        mpfr_sub_si(b1, pr, 2, MPFR_RNDN);        // b1 = pr - 2
        mpfr_div_d(b1, b1, lc->guard, MPFR_RNDN); // b1 = (pr - 2) / guard
        
        if(mpfr_cmp(b1, m) > 0) {
            if(levc_next(c, c, t, lc->period, prec, b2)) {
                if(epsUpdate) {
                    levc_set_point(lc, ind, c, b2);
                } else {
                    mpv_setc(lc->points, ind, c);
                }
                
                return true;
            } else {
                return false;
            }
        }
        
        mpfr_div(b1, b1, m, MPFR_RNDN);
        mpc_mulr(bt, bt, b1);
        mpc_add(bt, bt, pt);    // here bt = pt + (t - pt) * (pr - 2) / (lc->guard * |t - pt|)
        if(! levc_next(c, c, bt, lc->period, prec, b2)) {
            return false;
        }
        
        mpc_set(pt, bt);
        mpc_mod(pr, pt);
    } while(true);
}

/// @brief Lifts the curve @c lc to level @c r. It performs @b approximative computations.
///
/// @param lc the level curve to refine in place
/// @param r the new level of the curve
///
/// @return @ref true if successfull, @ref false otherwise
static bool levc_level_direct(levc lc, mpfr_t r) {
    mpfr_set_zero(lc->eps, 1); // will be recomputed for each point
    
    uint prec = lc->prec;
    uint tpow = lc->ang2pow;
    ulong count = levc_count(lc);
    
    mpc c, t, b, pt; // buffer for c, target, buffer, prev target
    mpc_init(c, prec);
    mpc_init(t, prec);
    mpc_init(b, prec);
    mpc_init(pt, prec);
    
    mpfr_t m, b1, pr, b2;
    mpfr_init2(m, prec);
    mpfr_init2(b1, prec);
    mpfr_init2(pr, prec);
    mpfr_init2(b2, prec);
    
    bool ok = true;
    if(lc->period > tpow) { // all values of p_n are r
        for (ulong i = 0; i < count && ok; i++) {
            mpc_setr(t, r);
            mpfr_set(pr, lc->radius, MPFR_RNDN);
            mpc_setr(pt, pr);
            
            ok = levc_move_point(lc, i, pt, pr, t, c, b, m, b1, b2, true);
        }
    } else { // need a circle for targets, faster than recomputing each time
        mpfr_t mr;
        mpfr_init2(mr, prec);
        mpfr_set(mr, lc->radius, MPFR_RNDN);
        mpfr_div(mr, mr, r, MPFR_RNDN);
        
        int lcpm1 = lc->period - 1;
        levc cr = levc_circle(tpow - lcpm1, prec, r);
        ulong mask = (1L << (tpow - lcpm1)) - 1;
        
        for (ulong i = 0; i < count && ok; i++) {
            mpv_getc(t, cr->points, i & mask);
            mpc_mulr(pt, t, mr);
            mpfr_set(pr, lc->radius, MPFR_RNDN);
            
            ok = levc_move_point(lc, i, pt, pr, t, c, b, m, b1, b2, true);
        }
        
        levc_free(cr);
        mpfr_clear(mr);
    }
    
    mpfr_clear(b2);
    mpfr_clear(pr);
    mpfr_clear(b1);
    mpfr_clear(m);
    
    mpc_clear(pt);
    mpc_clear(b);
    mpc_clear(t);
    mpc_clear(c);
    
    if(ok) {
        mpfr_set(lc->radius, r, MPFR_RNDN);
    }
    
    return ok;
}

levc levc_level(levc src, mpfr_t r) {
    if(src == NULL || r == NULL || src->period < 2 ||
       ! mpfr_regular_p(r) || mpfr_cmp_si(r, 2) <= 0) {
        return NULL;
    }
    
    levc lc = levc_copy(src, src->prec, true, true);
    
    if(mpfr_cmp(r, src->radius) == 0) {
        return lc;
    }

    if(levc_level_direct(lc, r)) {
        return lc;
    }
    
    levc_free(lc);
    
    return NULL;
}

/// @brief Increases the period of the curve @c lc by @c 1. It performs @b approximative computations.
///
/// @param lc the level curve to refine in place
///
/// @return @ref true if successfull, @ref false otherwise
static bool levc_next_period_direct(levc lc, bool squareLevel) {
    mpfr_set_zero(lc->eps, 1); // will be recomputed for each point
    
    // FIXME: tpow >= per; targets; start_angle
    
    uint prec = lc->prec;
    uint tpow = lc->ang2pow;
    ulong count = levc_count(lc);
    
    defs_mpc(prec, c, t, b, pt, nt); // buffer for c, target, buffer, prev target
    defs_mpfr(prec, m, b1, pr, b2, pr2);
    
    mpfr_ptr r = lc->radius;
    mpfr_set(pr2, r, MPFR_RNDN);
    mpfr_sqr(pr2, pr2, MPFR_RNDN);
    
    lc->period ++;
    int per = lc->period;
    mpv pts = lc->points;
    
    bool ok = true;
    if(per > tpow) { // all values of p_n are r
        mpc_setr(t, pr2);
        mpc_setr(nt, r);
        
        for (ulong i = 0; i < count && ok; i++) {
            mpv_getc(c, pts, i);
            mpc_add(pt, t, c);
            mpc_mod(pr, pt);
            
            ok = levc_move_point(lc, i, pt, pr, t, c, b, m, b1, b2, squareLevel);
            if(! squareLevel) {
                ok = ok && levc_move_point(lc, i, t, pr2, nt, c, b, m, b1, b2, true);
            }
        }
    } else { // need a circle for targets, faster than recomputing each time
        int tp = tpow - per + 1;
        levc cr = levc_circle(tp, prec, r);
        ulong mask = (1L << tp) - 1;
        
        for (ulong i = 0; i < count && ok; i++) {
            mpv_getc(nt, cr->points, i & mask);
            mpc_mulr(t, nt, r);
            
            mpv_getc(c, pts, i);
            mpc_add(pt, t, c);
            mpc_mod(pr, pt);
            
            ok = levc_move_point(lc, i, pt, pr, t, c, b, m, b1, b2, squareLevel);
            if(! squareLevel) {
                ok = ok && levc_move_point(lc, i, t, pr2, nt, c, b, m, b1, b2, true);
            }
        }
        
        levc_free(cr);
    }
    
    if(ok && squareLevel) {
        mpfr_set(lc->radius, pr2, MPFR_RNDN);
    }
        
    return ok;
}

levc levc_next_period(levc src, bool squareLevel) {
    if(src == NULL || src->period <= 1 || src->period == LEVC_MAX_PERIOD) {
        return NULL;
    }
    
    levc lc = levc_copy(src, src->prec, true, true);

    if(levc_next_period_direct(lc, squareLevel)) {
        return lc;
    }
    
    levc_free(lc);
    
    return NULL;
}

// MARK: IO functions

bool levc_write_csv(levc lc, char *fileName, uint digits) {
    if(lc == NULL || fileName == NULL || digits < 10 || digits > 1500) {
        return false;
    }
    
    FILE *f = fopen(fileName, "w");
    if(f == NULL) {
        return false;
    }
    
    if(fprintf(f, "%u %u %u\n", lc->period, lc->ang2pow, lc->prec) < 6) {
        fclose(f);
        
        return false;
    }
    
    if(fprintf(f, "%lu %lu\n", lc->startAngle, lc->endAngle) < 4) {
        fclose(f);
        
        return false;
    }
    
    char fmt[30];
    snprintf(fmt, 29, "%%.%dRf %%.%dRf\n", digits, digits);
    int len = 2 * digits + 200;
    char line[len + 1];
    int minLen = 2 * digits + 3;
    
    int n = mpfr_snprintf(line, len, fmt, lc->radius, lc->eps);
    if(n < minLen || n >= len || fprintf(f, "%s\n", line) < n) {
        fclose(f);
        
        return false;
    }
    
    if(fprintf(f, "%.6lf\n", lc->guard) < 10) {
        fclose(f);
        
        return false;
    }
    
    snprintf(fmt, 29, "%%.%dRf, %%.%dRf\n", digits, digits);
    ulong count = levc_count(lc);
    mpc c;
    mpc_init(c, lc->prec);
    for (ulong i = 0; i < count; i ++) {
        levc_point(c, lc, i);
        int n = mpfr_sprintf(line, fmt, c->x, c->y);
        if(n < minLen || n >= len || fprintf(f, "%s\n", line) < n) {
            fclose(f);
            mpc_clear(c);
            
            return false;
        }
    }
    
    fclose(f);
    mpc_clear(c);
    
    return true;
}

levc levc_read_csv(char *fileName) {
    if(fileName == NULL) {
        return NULL;
    }

    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return NULL;
    }

    int p, a, pr;
    if(fscanf(f, "%u %u %u", &p, &a, &pr) < 3 || p < 2 || a < 1 || pr < LEVC_MIN_PREC) {
        fclose(f);
        
        return NULL;
    }
    
    ulong sa, ea;
    if(fscanf(f, "%lu %lu", &sa, &ea) < 2 || ea < sa || ea > (1L << a)) {
        fclose(f);
        
        return NULL;
    }
    
    char sx[1600], sy[1600];
    mpfr_t r, eps;
    mpfr_init2(r, pr);
    mpfr_init2(eps, pr);
    if(fscanf(f, "%s%s", sx, sy) || mpfr_set_str(r, sx, 10, MPFR_RNDN) != 0
       || mpfr_set_str(eps, sy, 10, MPFR_RNDN) != 0) {
        fclose(f);
        mpfr_clears(r, eps, NULL);
        
        return NULL;
    }
    
    levc lc = levc_alloc(p, a, pr, sa, ea, r);
    mpfr_clear(r);
    
    if(lc == NULL) {
        levc_free(lc);
        fclose(f);
        mpfr_clear(eps);
        
        return NULL;
    }
    
    mpfr_set(lc->eps, eps, MPFR_RNDU);
    mpfr_clear(eps);
    
    int proof;
    if(fscanf(f, "%u %lf", &proof, &lc->guard) < 2 || (proof == 0 && lc->guard < LEVC_MIN_GUARD)) {
        levc_free(lc);
        fclose(f);
        
        return NULL;
    }
    
    bool ok = true;
    mpc c;
    mpc_init(c, pr);
    ulong count = levc_count(lc);
    for (ulong i = 0; i < count && ok; i ++) {
        if(fscanf(f, "%s%s", sx, sy) == 2) {
            long lx = strlen(sx);
            if(lx > 1 && sx[lx - 1] == ',') {
                sx[lx - 1] = 0;
            }
            
            bool ok = 0 == mpfr_set_str(c->x, sx, 10, MPFR_RNDN);
            ok = ok && 0 == mpfr_set_str(c->y, sy, 10, MPFR_RNDN);
            ok = ok && mpv_setc(lc->points, i, c);
        } else {
            ok = false;
        }
    }

    mpc_clear(c);
    fclose(f);

    if(! ok) {
        levc_free(lc);
        
        return NULL;
    }
    
    return lc;
}

bool levc_write_to(levc lc, char *fileName) {
    if(fileName == NULL) {
        return false;
    }

    FILE *f = fopen(fileName, "w");
    ulong wr = levc_write(lc, f, 0);
    fclose(f);
        
    return wr > 0 && wr == levc_fileLen(lc);
}

levc levc_read_from(char *fileName) {
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return false;
    }
    
    levc lc = levc_read(f, 0);
    
    fclose(f);
    
    return lc;
}

/// @brief Writes the file header of the curve @c lc to the file @c f, starting at absolute position @c pos if @c pos>=0,
/// at the current position otherwise.
///
/// @param lc the level curve
/// @param f the file
/// @param pos the absolute position if @c pos>=0, at the current position otherwise
///
/// @return @ref true if the operation completed successfully, @ref false otherwise
static bool levc_write_header(levc lc, FILE *f, long pos) {
    if(f == NULL) { // v already checked
        return false;
    }
        
    ulong flen = levc_fileLen(lc);
    mfile h = mfile_header(LEVC_FILE_ID, flen, levc_headerLen(lc));
    if(h == NULL) {
        return false;
    }
    
    mfile_puti(h, lc->period);
    mfile_puti(h, lc->ang2pow);
    mfile_puti(h, lc->prec);
    
    mfile_putl(h, lc->startAngle);
    mfile_putl(h, lc->endAngle);
    
    mwrite_mpfr(lc->radius, h);
    mwrite_mpfr(lc->eps, h);
    
    mfile_putb(h, 0);
    mfile_putd(h, lc->guard);
    
    // compute MD5 checksum
    MD5_CTX md5;
    byte md5sum[16];
    
    MD5_Init(&md5);
    bool ok = MD5_Update(&md5, h->data, (uint) h->len);
    
    ok = ok && mpv_update_md5(lc->points, &md5);
    
    ok = ok && MD5_Final(md5sum, &md5);
    mfile_putbs(h, md5sum, 16);
    
    ok = ok && h->len == levc_headerLen(lc) && mfile_write_to(f, pos, h, 0, h->len);
    
    mfile_free(h);
    
    return ok;
}

ulong levc_write(levc lc, FILE *f, long pos) {
    if(lc == NULL) {
        return 0;
    }
    
    if(! levc_write_header(lc, f, pos)) {
        return 0;
    }
    
    return levc_headerLen(lc) + mpv_write_points(lc->points, f, -1);
}

/// @brief Reads the header into a buffer levelCurve and computes the partial @c MD5 sum.
///
/// The returned object is not fully formed, do @b not use levc_free() to destroy it.
///
/// @param f the file
/// @param pos the absolute position in the file, unless @c pos<0, in which case the current position is used
/// @param md5 the @c MD5 checksum of the header
/// @param fmd5sum the @c MD5 checksum stored in the file
///
/// @return a levelCurve that represents the header of the sored level curve
static levc levc_read_header_md5(FILE *f, long pos, MD5_CTX *md5, byte fmd5sum[]) {
    mfile h = mfile_read_header(f, pos, LEVC_FILE_ID, LEVC_MIN_HEADER_LEN);
    if(h == NULL) {
        return NULL;
    }
    
    levc buf = malloc(sizeof(levelCurve));
    
    h->pos = 20; // skip de file ID, length of the file and of the header
    buf->period = mfile_geti(h);
    buf->ang2pow = mfile_geti(h);
    buf->prec = mfile_geti(h);
    buf->startAngle = mfile_getl(h);
    buf->endAngle = mfile_getl(h);
    
    if(! mread_mpfr(buf->radius, h)) {
        mfile_free(h);
        free(buf);
        
        return NULL;
    }
    
    if(! mread_mpfr(buf->eps, h)) {
        mpfr_clear(buf->radius);
        mfile_free(h);
        free(buf);
        
        return NULL;
    }
    
    mfile_getb(h);
    buf->guard = mfile_getd(h);
    
    mfile_getbs(h, fmd5sum, 16);
    
    long mpos = h->pos;
    bool ok = MD5_Update(md5, h->data, (uint) mpos - 16); // except the MD5
    
    mfile_free(h);
    
    if(! ok || mpos != levc_headerLen(buf)) {
        mpfr_clear(buf->radius);
        mpfr_clear(buf->eps);
        free(buf);
        
        return NULL;
    }
    
    return buf;
}

/// @brief Reads the header into a buffer levelCurve.
///
/// The returned object is not fully formed, do @b not use levc_free() to destroy it.
///
/// @param f the file
/// @param pos the absolute position in the file, unless @c pos<0, in which case the current position is used
///
/// @return a levelCurve that represents the header of the sored level curve
static levc levc_read_header(FILE *f, long pos) {
    mfile h = mfile_read_header(f, pos, LEVC_FILE_ID, LEVC_MIN_HEADER_LEN);
    if(h == NULL) {
        return NULL;
    }
    
    levc buf = malloc(sizeof(levelCurve));
    
    h->pos = 20; // skip de file ID, length of the file and of the header
    buf->period = mfile_geti(h);
    buf->ang2pow = mfile_geti(h);
    buf->prec = mfile_geti(h);
    buf->startAngle = mfile_getl(h);
    buf->endAngle = mfile_getl(h);
    
    if(! mread_mpfr(buf->radius, h)) {
        mfile_free(h);
        free(buf);
        
        return NULL;
    }
    
    if(! mread_mpfr(buf->eps, h)) {
        mpfr_clear(buf->radius);
        mfile_free(h);
        free(buf);
        
        return NULL;
    }
    
    mfile_getb(h);
    buf->guard = mfile_getd(h);
    
    long mpos = h->pos;
    
    mfile_free(h);
    
    if(mpos != levc_headerLen(buf) - 16) {
        mpfr_clear(buf->radius);
        mpfr_clear(buf->eps);
        free(buf);
        
        return NULL;
    }
    
    return buf;
}

levc levc_read(FILE *f, long pos) {
    byte fmd5sum[16];
    
    // compute MD5 checksum
    MD5_CTX md5;
    byte md5sum[16];
    
    MD5_Init(&md5);
    levc buf = levc_read_header_md5(f, pos, &md5, fmd5sum);
    
    if(buf == NULL) {
        return NULL;
    }
            
    // create the new curve
    levc lc = levc_alloc(buf->period, buf->ang2pow, buf->prec, buf->startAngle,
                         buf->endAngle, buf->radius);
    mpfr_set(lc->eps, buf->eps, MPFR_RNDN);
    lc->guard = buf->guard;

    mpfr_clear(buf->radius);
    mpfr_clear(buf->eps);
    free(buf);
    
    if(lc == NULL || ! mpv_read_points(lc->points, 0, lc->points->count, f, -1)) {
        levc_free(lc);
        
        return NULL;
    }
    
    bool ok = mpv_update_md5(lc->points, &md5);
    ok = ok && MD5_Final(md5sum, &md5);
    
    for (int i = 0; i < 16 && ok; i++) {
        ok = ok && fmd5sum[i] == md5sum[i];
    }
    
    if(! ok) {
        levc_free(lc);
        
        return NULL;
    }
    
    return lc;
}

levc levc_read_partial(FILE *f, long pos, ulong start, ulong step, ulong count) {
    levc buf = levc_read_header(f, pos);
            
    // check the parameters
    int l2 = levc_log2(step);
    ulong end = buf->startAngle + start + step * (count - 1);
    if(l2 < 0 || count < 2 || start % step != 0 || end > buf->endAngle) {
        free(buf);
        
        return NULL;
    }
    
    // create the new curve
    levc lc = levc_alloc(buf->period, buf->ang2pow - l2, buf->prec, buf->startAngle >> l2,
                         buf->endAngle >> l2, buf->radius);
    mpfr_set(lc->eps, buf->eps, MPFR_RNDN);
    lc->guard = buf->guard;

    mpfr_clear(buf->radius);
    mpfr_clear(buf->eps);
    free(buf);
            
    ulong tc = levc_count(buf);
    if(lc == NULL || ! mpv_read_points_partial(lc->points, tc, f, -1, start, step, count, true)) {
        levc_free(lc);
        
        return NULL;
    }
    
    return lc;
}
