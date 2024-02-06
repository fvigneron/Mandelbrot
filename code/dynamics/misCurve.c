//
//  misCurve.c
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

#include "misCurve.h"

// MARK: Creation, computation and destruction

/// @brief Creates an empty misCurve of pre-period @c pp, period @c per, angles of @c 2^{-tpow},
/// precision @c prec, start angle @c sta, end angle @c ena and radius @c r.
///
/// @param pp the pre-period of the level curve
/// @param per the period of the level curve
/// @param tpow the power of two to obtain the smalles dyadic angle
/// @param sta the start angle
/// @param ena the end angle
/// @param prec the precision
/// @param r the radius
///
/// @return the level curve
static levm levm_alloc(uint pp, uint per, int tpow, uint prec, ulong sta, ulong ena, mpfr_t r) {
    if(pp < 2 || per < 1 || tpow == 0 || tpow > LEVM_MAX_2POW || prec < LEVM_MIN_PREC || r == NULL
       || ! mpfr_regular_p(r) || mpfr_cmp_si(r, per == 1 ? 0 : 4) <= 0
       || ena < sta || ena >= (1L << tpow)) {
        return NULL;
    }
    
    levm lm = malloc(sizeof(misCurve));
    lm->prePeriod = pp;
    lm->period = per;
    lm->ang2pow = tpow;
    lm->prec = prec;
    lm->startAngle = sta;
    lm->endAngle = ena;
    
    mpfr_init2(lm->radius, prec);
    mpfr_set(lm->radius, r, MPFR_RNDN);
    
    lm->guard = LEVM_MIN_GUARD;
    
    mpv_init(lm->points, prec, 2 * levm_count(lm));
    
    return lm;
}

/// @brief Creates an empty misCurve of pre-period @c pp, period @c per representing the upper half circle, with angles
/// of @c 2^{-tpow}, precision @c prec and radius @c r.
///
/// @param pp the pre-period of the level curve
/// @param per the period of the level curve
/// @param tpow the power of two to obtain the smalles dyadic angle
/// @param prec the precision
/// @param r the radius
///
/// @return the level curve
static levm levm_alloc_upper(uint pp, uint per, int tpow, uint prec, mpfr_t r) {
    return levm_alloc(pp, per, tpow, prec, 0, 1L << (tpow - 1), r);
}

/// @brief Computes in place the misCurve which has the endpoints initialized , using  fp80.h.
///
/// @param lm the level curve to compute
/// @param r the radius or level of the curve
///
/// @return @ref true if the computation succeeded, @ref false otherwise
static bool levm_fill_inl(levm lm, ldbl r) {
    // compute how many targets are needed w.r.t. lm->radius and lm->guard
    ldbl a = 2 * PI * lm->guard * r / (r - 2);
    int tpowc = 4; // two power of the circle of targets
    while((1L << tpowc) <= a && tpowc < 30) {
        tpowc ++;
    }
    
    if(tpowc == 30) {
        return false;
    }
    
    int pp = lm->prePeriod;
    int per = lm->period;
    int n = pp + per;
    int a2p = lm->ang2pow;
    int lctp = a2p - n + 1;
    tpowc = lctp > tpowc ? lctp : tpowc;     // if points on the curve are even denser
    
    int tpow = tpowc + n - 1;        // the power of two for the angles, >= lm->ang2pow
    // for any angle in the scale tpow, the mask gives the index of the target
    ulong tpowMask = (1L << tpowc) - 1;
    ldbl mula = ldexpl(PI, -tpowc + 1); // coefficient to get the actual angle
    // for every point on the level curve, there will be 2^extp points computed
    int extp = tpow - a2p;
    ulong extpMask = (1L << extp) - 1; // mask to check when to store a point
    
    // compute intermediary points
    ulong sta = lm->startAngle << extp;
    ulong ena = lm->endAngle << extp;
    
    bool ok = true;
    fp80 c, t;
    ulong am, i;
    
    mpv pts = lm->points;
    mpv_getcl(c, pts, 0); // start point on the right
    for (ulong ula = sta; ula <= ena && ok; ula++) {
        am = ula & tpowMask;
        a = am * mula;
        
        t->x = r * cosl(a);
        t->y = r * sinl(a);
        
        ok = mandel_mis_root_refl(c, c, t, pp, per, LEVM_MAX_ITER, 1E-17, 2, 1);
        
        if((ula & extpMask) == 0 && ok) { // store the new value
            i = (ula - sta) >> extp;
            
            if(ula == sta || ula == ena) {
                mpv_getcl(t, pts, i);
                
                ok = fp80_dist(c, t) <= 1E-17;
            }
            
            mpv_setcl(pts, i, c);
        }
    }
    
    return ok;
}

/// @brief Computes in place the misCurve which has the endpoints initialized , using  fp80.h.
///
/// @param lm the level curve to compute
/// @param r the radius or level of the curve
///
/// @return @ref true if the computation succeeded, @ref false otherwise
static bool levm_fill_simple_inl(levm lm, ldbl r) {
    // compute how many targets are needed w.r.t. lm->radius and lm->guard
    ldbl a = 2 * PI * lm->guard * r / (r - 2);
    int tpowc = 4; // two power of the circle of targets
    while((1L << tpowc) <= a && tpowc < 30) {
        tpowc ++;
    }
    
    if(tpowc == 30) {
        return false;
    }
    
    int pp = lm->prePeriod;
    int per = lm->period;
    int n = pp + per - 1;
    int a2p = lm->ang2pow;
    int lctp = a2p - n + 1;
    tpowc = lctp > tpowc ? lctp : tpowc;     // if points on the curve are even denser
    
    int tpow = tpowc + n - 1;        // the power of two for the angles, >= lm->ang2pow
    // for any angle in the scale tpow, the mask gives the index of the target
    ulong tpowMask = (1L << tpowc) - 1;
    ldbl mula = ldexpl(PI, -tpowc + 1); // coefficient to get the actual angle
    // for every point on the level curve, there will be 2^extp points computed
    int extp = tpow - a2p;
    ulong extpMask = (1L << extp) - 1; // mask to check when to store a point
    
    // compute intermediary points
    ulong sta = lm->startAngle << extp;
    ulong ena = lm->endAngle << extp;
    
    bool ok = true;
    fp80 c, t;
    ulong am, i;
    
    mpv pts = lm->points;
    mpv_getcl(c, pts, 0); // start point on the right
    for (ulong ula = sta; ula <= ena && ok; ula++) {
        am = ula & tpowMask;
        a = am * mula;
        
        t->x = r * cosl(a);
        t->y = r * sinl(a);
        
        ok = mandel_miss_root_refl(c, c, t, pp, per, LEVM_MAX_ITER, 1E-17, 2, 1);
        
        if((ula & extpMask) == 0 && ok) { // store the new value
            i = (ula - sta) >> extp;
            
            if(ula == sta || ula == ena) {
                mpv_getcl(t, pts, i);
                
                ok = fp80_dist(c, t) <= 1E-17;
            }
            
            mpv_setcl(pts, i, c);
        }
    }
    
    return ok;
}

/// @brief Computes the next point @c ep from the starting point @c sp with target @c t and period @c per.
///
/// It uses fp80 as long as the modulus of the Newton term is amller than LEVC_MAX_FP80_STEP, then refines the
/// result up to the precision of @c ep.
///
/// @param ep the result
/// @param sp the start point
/// @param t the target
/// @param pp the pre-period
/// @param per the period
/// @param prec the precision
///
/// @return @ref true if the computation succeeded, @ref false otherwise
static bool levm_next(mpc ep, mpc sp, mpc t, int pp, int per, uint prec) {
    fp80 c8, t8, nt8, v8, d8;
    mpc_get80(c8, sp);
    mpc_get80(t8, t);
    
    if(prec == LEVM_MIN_PREC) { // try all in fp80
        if(mandel_mis_root_refl(c8, c8, t8, pp, per, LEVC_MAX_ITER, 1E-17, 2, 1)) {
            mpc_set80(ep, c8);
            
            return true;
        } else { // no point in trying with mp if the precision is minimal
            return false;
        }
    }
    
    mandel_mis_val_derl(v8, d8, c8, pp, per);
    fp80_sub(v8, v8, t8);
    ldbl md = fp80_mod(d8);
    ldbl mnt = fp80_mod(v8) / md;
    
    if(mnt < LEVM_MIN_FP80_STEP || md * LEVM_MIN_FP80_STEP > 1) {
        return mandel_mis_root_ref(ep, sp, t, pp, per, LEVM_MAX_ITER, 2);
    }
    
    int iter = 0;
    do {
        fp80_div(nt8, v8, d8);
        fp80_sub(c8, c8, nt8);
        iter ++;
        
        if(iter < LEVM_MAX_ITER) {
            mandel_mis_val_derl(v8, d8, c8, pp, per);
            fp80_sub(v8, v8, t8);
            md = fp80_mod(d8);
            mnt = fp80_mod(v8) / md;
        } else {
            return false; // no more iterates to refine
        }
    } while(mnt >= LEVM_MIN_FP80_STEP && md * LEVM_MIN_FP80_STEP <= 1);
    
    mpc_set80(sp, c8);
    
    return mandel_mis_root_ref(ep, sp, t, pp, per, LEVM_MAX_ITER - iter, 2);
}

/// @brief Computes in place the misCurve which has endpoints initialized, using by default mpc.h,
/// or fp80.h if the precision requested is LEVM_MIN_PREC.
///
/// @param lm the misCurve to compute
///
/// @return @ref true if the computation succeeded and the endpoints error was bounded cy @c eps, @ref false otherwise
static bool levm_fill_in(levm lm) {
    if(levm_segs(lm) == 1) { // only two points, no fill in
        return true;
    }
    
    ldbl r = mpfr_get_ld(lm->radius, MPFR_RNDN);
    
    if(lm->prec == LEVM_MIN_PREC) { // low prec requested
        return levm_fill_inl(lm, r);
    }
    
    // compute how many targets are needed w.r.t. lm->radius and lm->guard
    ldbl ra = 2 * PI * lm->guard * r / (r - 2);
    int tpowc = 4; // two power of the circle of targets
    while((1L << tpowc) <= ra && tpowc < 30) {
        tpowc ++;
    }
    
    if(tpowc == 30) {
        return false;
    }
    
    int pp = lm->prePeriod;
    int per = lm->period;
    int n = pp + per;
    int a2p = lm->ang2pow;
    int lctp = a2p - n + 1;
    tpowc = lctp > tpowc ? lctp : tpowc;     // if points on the curve are even denser
    
    int tpow = tpowc + n - 1;        // the power of two for the angles, >= lm->ang2pow
    // for any angle in the scale tpow, the mask gives the index of the target
    ulong tpowMask = (1L << tpowc) - 1;
    // for every point on the level curve, there will be 2^extp points computed
    int extp = tpow - a2p;
    ulong extpMask = (1L << extp) - 1; // mask to check when to store a point
    
    // compute intermediary points
    ulong sta = lm->startAngle << extp;
    ulong ena = lm->endAngle << extp;
    
    bool ok = true;
    
    uint prec = lm->prec;
    levc tc = levc_circle(tpowc, prec, lm->radius);
    mpv ts = tc->points;
    
    mpc c, t;
    mpc_init(c, prec);
    mpc_init(t, prec);
    
    mpfr_t err, d;
    mpfr_init2(err, prec);
    mpfr_init2(d, prec);
    
    mpv pts = lm->points;
    mpv_getc(c, pts, 0); // start point on the right
    
    // use the targets to compute intermediary points, store the nodes of lm
    for (ulong a = sta; a <= ena && ok; a++) {
        mpv_getc(t, ts, a & tpowMask);
        
        ok = levm_next(c, c, t, pp, per, prec);
        
        if((a & extpMask) == 0 && ok) { // store the new value
            if(a == sta || a == ena) {
                mpv_getc(t, pts, (a - sta) >> extp);
                mpc_dist(d, c, t);
            }
            
            mpv_setc(lm->points, (a - sta) >> extp, c);
        }
    }
    
    mpfr_clear(err);
    mpfr_clear(d);
    
    mpc_clear(t);
    mpc_clear(c);
    
    levc_free(tc);
    
    return ok;
}

/// @brief Computes in place the simple misCurve which has endpoints initialized, using by default mpc.h,
/// or fp80.h if the precision requested is LEVM_MIN_PREC.
///
/// @param lm the misCurve to compute
///
/// @return @ref true if the computation succeeded and the endpoints error was bounded cy @c eps, @ref false otherwise
static bool levm_fill_simple_in(levm lm) {
    if(levm_segs(lm) == 1) { // only two points, no fill in
        return true;
    }
    
    ldbl r = mpfr_get_ld(lm->radius, MPFR_RNDN);
    
    if(lm->prec == LEVM_MIN_PREC) { // low prec requested
        return levm_fill_simple_inl(lm, r);
    }
    
    // compute how many targets are needed w.r.t. lm->radius and lm->guard
    ldbl ra = 2 * PI * lm->guard * r / (r - 2);
    int tpowc = 4; // two power of the circle of targets
    while((1L << tpowc) <= ra && tpowc < 30) {
        tpowc ++;
    }
    
    if(tpowc == 30) {
        return false;
    }
    
    int pp = lm->prePeriod;
    int per = lm->period;
    int n = pp + per - 1;
    int a2p = lm->ang2pow;
    int lctp = a2p - n + 1;
    tpowc = lctp > tpowc ? lctp : tpowc;     // if points on the curve are even denser
    
    int tpow = tpowc + n - 1;        // the power of two for the angles, >= lm->ang2pow
    // for any angle in the scale tpow, the mask gives the index of the target
    ulong tpowMask = (1L << tpowc) - 1;
    // for every point on the level curve, there will be 2^extp points computed
    int extp = tpow - a2p;
    ulong extpMask = (1L << extp) - 1; // mask to check when to store a point
    
    // compute intermediary points
    ulong sta = lm->startAngle << extp;
    ulong ena = lm->endAngle << extp;
    
    bool ok = true;
    
    uint prec = lm->prec;
    levc tc = levc_circle(tpowc, prec, lm->radius);
    mpv ts = tc->points;
    
    mpc c, t;
    mpc_init(c, prec);
    mpc_init(t, prec);
    
    mpfr_t err, d;
    mpfr_init2(err, prec);
    mpfr_init2(d, prec);
    
    mpv pts = lm->points;
    mpv_getc(c, pts, 0); // start point on the right
    
    // use the targets to compute intermediary points, store the nodes of lm
    for (ulong a = sta; a <= ena && ok; a++) {
        mpv_getc(t, ts, a & tpowMask);
        
        // no optimization in fp80
        ok = mandel_miss_root_ref(c, c, t, pp, per, LEVM_MAX_ITER, 2);
        
        if((a & extpMask) == 0 && ok) { // store the new value
            if(a == sta || a == ena) {
                mpv_getc(t, pts, (a - sta) >> extp);
                mpc_dist(d, c, t);
            }
            
            mpv_setc(lm->points, (a - sta) >> extp, c);
        }
    }
    
    mpfr_clear(err);
    mpfr_clear(d);
    
    mpc_clear(t);
    mpc_clear(c);
    
    levc_free(tc);
    
    return ok;
}

/// @brief Computes in place the level curve of the upper half circle using mpc.h.
///
/// @param lm the level curve to compute
static bool levm_compute_direct(levm lm) {
    int pp = lm->prePeriod;
    int per = lm->period;
    int n = pp + per;
    
    // init start and end points with low precision
    ldbl ep = 1E-16;
    
    ldbl r = mpfr_get_ld(lm->radius, MPFR_RNDN);
    
    fp80 sp = {-2, 0}, t = {r, 0}, left, right;
    bool ok = mandel_sol_refl(sp, NULL, sp, t, n, 1000, ep, 1);
    ok = ok && mandel_mis_root_refl(left, sp, t, pp, per, 1000, ep, 1, 1);
    
    sp->x = 0.25 + 8 / pow(n, 1.9);
    ok = ok && mandel_sol_refl(sp, NULL, sp, t, n, 1000, ep, 1);
    ok = ok && mandel_mis_root_refl(right, sp, t, pp, per, 1000, ep, 1, 1);
    
    if(! ok) {
        return false;
    }
    
    mpv_setcl(lm->points, 0, right);
    mpv_setcl(lm->points, levm_segs(lm), left);
    
    long prec = lm->prec;
    if(prec > LEVM_MIN_PREC) {
        mpc c, vs, ve, t;
        mpc_init(c, prec);
        mpc_init(vs, prec);
        mpc_init(ve, prec);
        mpc_init(t, prec);
        
        mpv_getc(c, lm->points, 0);
        mpc_setr(t, lm->radius);
        ok = ok && mandel_mis_root_ref(vs, c, t, pp, per, 20, 3);
        
        mpv_getc(c, lm->points, levm_segs(lm));
        ok = ok && mandel_mis_root_ref(ve, c, t, pp, per, 100, 3);
        
        if(ok) {
            mpv_setc(lm->points, 0, vs);
            mpv_setc(lm->points, levm_segs(lm), ve);
        }
        
        mpc_clears(c, vs, ve, t, NULL);
    }
        
    ok = ok && levm_fill_in(lm);
    
    return ok;
}

/// @brief Computes in place the level curve of the upper half circle using mpc.h.
///
/// @param lm the level curve to compute
static bool levm_compute_simple_direct(levm lm) {
    int pp = lm->prePeriod;
    int per = lm->period;
    int n = pp + per - 1;
    
    // init start and end points with low precision
    ldbl ep = 1E-16;
    
    ldbl r = mpfr_get_ld(lm->radius, MPFR_RNDN);
    
    fp80 sp = {-2, 0}, t = {r, 0}, left, right;
    bool ok = mandel_sol_refl(sp, NULL, sp, t, n, 1000, ep, 1);
    ok = ok && mandel_miss_root_refl(left, sp, t, pp, per, 1000, ep, 1, 1);
    
    sp->x = 0.25 + 8 / pow(n, 1.9);
    ok = ok && mandel_sol_refl(sp, NULL, sp, t, n, 1000, ep, 1);
    ok = ok && mandel_miss_root_refl(right, sp, t, pp, per, 1000, ep, 1, 1);
    
    if(! ok) {
        return false;
    }
    
    mpv_setcl(lm->points, 0, right);
    mpv_setcl(lm->points, levm_segs(lm), left);
    
    long prec = lm->prec;
    if(prec > LEVM_MIN_PREC) {
        mpc c, vs, ve, t;
        mpc_init(c, prec);
        mpc_init(vs, prec);
        mpc_init(ve, prec);
        mpc_init(t, prec);
        
        mpv_getc(c, lm->points, 0);
        mpc_setr(t, lm->radius);
        ok = ok && mandel_miss_root_ref(vs, c, t, pp, per, 20, 3);
        
        mpv_getc(c, lm->points, levm_segs(lm));
        ok = ok && mandel_miss_root_ref(ve, c, t, pp, per, 100, 3);
        
        if(ok) {
            mpv_setc(lm->points, 0, vs);
            mpv_setc(lm->points, levm_segs(lm), ve);
        }
        
        mpc_clears(c, vs, ve, t, NULL);
    }
        
    ok = ok && levm_fill_simple_in(lm);
    
    return ok;
}

levm levm_new(uint pp, uint per, uint tpow, uint prec, mpfr_t r, double guard) {
    if(guard < LEVM_MIN_GUARD || (prec == LEVM_MIN_PREC && pp + per > LEVM_MAX_FP80_PREC)) {
        return NULL;
    }
    
    levm lm = levm_alloc_upper(pp, per, tpow, prec, r);
    if(lm == NULL) {
        return NULL;
    }
    
    lm->guard = guard;
    
    if(levm_compute_direct(lm)) {
        return lm;
    }
    
    levm_free(lm);
    
    return NULL;
}

levm levm_simple_new(uint pp, uint per, uint tpow, uint prec, mpfr_t r, double guard) {
    if(guard < LEVM_MIN_GUARD || (prec == LEVM_MIN_PREC && pp + per > LEVM_MAX_FP80_PREC)) {
        return NULL;
    }
    
    levm lm = levm_alloc_upper(pp, per, tpow, prec, r);
    if(lm == NULL) {
        return NULL;
    }
    
    lm->guard = guard;
    
    if(levm_compute_simple_direct(lm)) {
        return lm;
    }
    
    levm_free(lm);
    
    return NULL;
}

void levm_free(levm lm) {
    if(lm == NULL || lm->points->vals == NULL || lm->points->sexp == NULL) {
        return;   // probably already cleared
    }
    
    mpv_clear(lm->points);
    mpfr_clear(lm->radius);
    
    free(lm);
}

// MARK: Pointwise operations

bool levm_angle(mpfr_t angle, levm lm, ulong ind) {
    if(angle == NULL || lm == NULL || ind >= levm_count(lm)) {
        return false;
    }
    
    mpfr_set_ui(angle, lm->startAngle + ind, MPFR_RNDN);
    mpfr_mul_2si(angle, angle, -lm->ang2pow, MPFR_RNDN);
    
    return true;
}

ldbl levm_anglel(levm lm, ulong ind) {
    if(lm == NULL || ind >= levm_count(lm)) {
        return false;
    }
    
    ldbl angle = lm->startAngle + ind;
    angle = ldexpl(angle, -lm->ang2pow);
    
    return angle;
}

bool levm_point(mpc p, levm lm, ulong ind) {
    if(p == NULL || lm == NULL || ind >= levm_count(lm)) {
        return false;
    }
    
    return mpv_getc(p, lm->points, ind);
}

bool levm_pointl(fp80 p, levm lm, ulong ind) {
    if(p == NULL || lm == NULL || ind >= levm_count(lm)) {
        return false;
    }
    
    return mpv_getcl(p, lm->points, ind);
}

// MARK: Refinements and derivations

/// @brief Return the base two logarithm of @c n if it is a positive power of @c 2, @c -1 otherwise.
///
/// @param n the number
///
/// @return @c k if @c n==2^k and @c k>=0, @c -1 otherwise
static inline int levm_log2(ulong n) {
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
static levm levm_copy(levm src, uint prec, bool initPoints, bool copyPoints) {
    levm lm = malloc(sizeof(misCurve));
    
    lm->prePeriod = src->prePeriod;
    lm->period = src->period;
    lm->ang2pow = src->ang2pow;
    lm->prec = prec;
    lm->startAngle = src->startAngle;
    lm->endAngle = src->endAngle;
    
    mpfr_init2(lm->radius, lm->prec);
    mpfr_set(lm->radius, src->radius, MPFR_RNDN);
    lm->guard = src->guard;
    
    if(initPoints) {
        mpv_init(lm->points, lm->prec, src->points->count);
        
        if(copyPoints) {
            mpv_copy(lm->points, 0, src->points);
        }
    }
    
    return lm;
}

levm levm_sub_curve(levm src, ulong start, ulong step, ulong count) {
    int l2 = levm_log2(step);
    
    if(src == NULL || l2 < 0 || start % step != 0) {
        return NULL;
    }
    
    mpv pts = mpv_sub_vectorc(src->points, start, step, count);
    if(pts == NULL) {
        return NULL;
    }
    
    levm lm = levm_copy(src, src->prec, false, false);
    
    lm->ang2pow = src->ang2pow - l2;
    lm->startAngle = (src->startAngle + start) >> l2;
    lm->endAngle = lm->startAngle + count - 1;
    
    *lm->points = *pts;
    free(pts);              // ATTENTION: not mpv_free() !
    
    return lm;
}

/// @brief Refines the value of the point with index @c ind in the misCurve @c lm.
///
/// @param lm the misCurve
/// @param ind the index in @c lm
/// @param t the target, the value of @c p_n in the new found point
/// @param c buffer for computations
/// @param nt buffer for computations
/// @param m buffer for computations
///
/// @return @ref true if the computation succeeded and the error was bounded cy @c eps, @ref false otherwise
static bool levm_refine_point(levm lm, ulong ind, mpc t, mpc c, mpc nt, mpfr_t m) {
    mpv_getc(c, lm->points, ind);
    
    bool conv = 0, div = 0, per = lm->period, prec = lm->prec;
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
    
    mpv_setc(lm->points, ind, c); // the refined value
    
    return true;
}

/// @brief Refines their precision to @c prec.
///
/// It is assumed that @c prec is larger than the previous precision of the points, thus it cannot be the minimal precision
/// and so all computations are performed with tools from mpc.h.
///
/// @param lm the misCurve to refine in place
///
/// @return @ref true if successfull, @ref false otherwise
static bool levm_refine_prec_direct(levm lm) {
    uint prec = lm->prec;
    uint tpow = lm->ang2pow;
    mpfr_ptr r = lm->radius;
    ulong count = levm_count(lm);
    
    mpc c, t, nt;
    mpc_init(c, prec);
    mpc_init(t, prec);
    mpc_init(nt, prec);
    bool ok = true;
    
    mpfr_t m;
    mpfr_init2(m, LEVC_PREC_EPS);
    if(lm->period > tpow) { // all values of p_n will be r
        mpc_setr(t, r);
        
        for (ulong i = 0; i < count && ok; i++) {
            ok = levm_refine_point(lm, i, t, c, nt, m);
        }
    } else { // need a circle for targets, faster than recomputing each time
        int lcpm1 = lm->period - 1;
        ulong mask = (1L << lcpm1) - 1;
        levc cr = levc_circle(tpow - lcpm1, prec, r);
        
        for (ulong i = 0; i < count && ok; i++) {
            mpv_getc(t, cr->points, i & mask);
            
            ok = levm_refine_point(lm, i, t, c, nt, m);
        }
        
        levc_free(cr);
    }
    
    mpfr_clear(m);
    mpc_clear(nt);
    mpc_clear(t);
    mpc_clear(c);
    
    return ok;
}

levm levm_refine(levm src, uint tpow, uint prec) {
    if(src == NULL || tpow > LEVM_MAX_2POW || tpow == 0 || prec < LEVM_MIN_PREC) {
        return NULL;
    }
    
    if(tpow == src->ang2pow) {
        levm lm = levm_copy(src, prec, true, true);
        
        if(prec <= src->prec) {
            return lm;
        }
        
        // more precision is needed
        if(levm_refine_prec_direct(lm)) {
            return lm;
        }
        
        levm_free(lm);
        
        return NULL;
    } else if(tpow < src->ang2pow) {
        int l2 = src->ang2pow - tpow;
        levm lm = levm_sub_curve(src, 0, 1L << l2, (levm_segs(src) >> l2) + 1);
        
        if(prec == lm->prec) {
            return lm;
        }
        
        levm nc = levm_refine(lm, lm->ang2pow, prec);
        levm_free(lm);
        
        return nc;
    }
        
    // more points are needed
    int l2 = tpow - src->ang2pow;
    ulong sta = src->startAngle << l2;
    ulong ena = src->endAngle << l2;
    levm lm = levm_alloc(src->prePeriod, src->period, tpow, prec, sta, ena, src->radius);
    
    mpv_get_setc(lm->points, 0, src->points, 0);
    mpv_get_setc(lm->points, levm_segs(lm), src->points, levm_segs(src));
    
    if(levm_fill_in(lm)) {
        return lm;
    }
    
    levm_free(lm);
    
    return NULL;
}

/// @brief Replaces the point fo index @c ind of the curve @c lm (whic has target @c pt of modulus @c pr) with the point
/// computet with target @c t.
///
/// If the @c lm->guard condition is not satisfied, it divides accordingly the segment @c [pt, t] and computes in several steps.
///
/// @warning All parameters are potentially used as buffers, the user should safely keep a copy of their values.
///
/// @param lm the level curve
/// @param ind the index in @c lm
/// @param pt previous target
/// @param pr previous radius, or @c |pt|
/// @param t the new target
/// @param c buffer
/// @param bt buffer
/// @param m buffer
/// @param b1 buffer
/// @param b2 buffer
///
/// @return @ref true if successfull, @ref false otherwise
static bool levm_move_point(levm lm, ulong ind, mpc pt, mpfr_t pr, mpc t, mpc c,
                              mpc bt, mpfr_t m, mpfr_t b1, mpfr_t b2) {
    mpv_getc(c, lm->points, ind);
    uint prec = lm->prec;
    
    do {
        mpc_sub(bt, t, pt);                       // bt = t - pt
        mpc_mod(m, bt);                           // m = |bt| = |t - pt|
        mpfr_sub_si(b1, pr, 4, MPFR_RNDN);        // b1 = pr - 4
        mpfr_div_d(b1, b1, lm->guard, MPFR_RNDN); // b1 = (pr - 4) / guard
        
        if(mpfr_cmp(b1, m) > 0) {
            if(levm_next(c, c, t, lm->prePeriod, lm->period, prec)) {
                mpv_setc(lm->points, ind, c);
                
                return true;
            } else {
                return false;
            }
        }
        
        mpfr_div(b1, b1, m, MPFR_RNDN);
        mpc_mulr(bt, bt, b1);
        mpc_add(bt, bt, pt);    // here bt = pt + (t - pt) * (pr - 4) / (lm->guard * |t - pt|)
        if(! levm_next(c, c, bt, lm->prePeriod, lm->period, prec)) {
            return false;
        }
        
        mpc_set(pt, bt);
        mpc_mod(pr, pt);
    } while(true);
}

/// @brief Replaces the point fo index @c ind of the curve @c lm (whic has target @c pt of modulus @c pr) with the point
/// computet with target @c t. Uses the simplified Misiurewicz polynomial.
///
/// If the @c lm->guard condition is not satisfied, it divides accordingly the segment @c [pt, t] and computes in several steps.
/// Does not attempt to optimize using fp80.
///
/// @warning All parameters are potentially used as buffers, the user should safely keep a copy of their values.
///
/// @param lm the level curve
/// @param ind the index in @c lm
/// @param pt previous target
/// @param pr previous radius, or @c |pt|
/// @param t the new target
/// @param c buffer
/// @param bt buffer
/// @param m buffer
/// @param b1 buffer
/// @param b2 buffer
///
/// @return @ref true if successfull, @ref false otherwise
static bool levm_move_simple_point(levm lm, ulong ind, mpc pt, mpfr_t pr, mpc t, mpc c,
                              mpc bt, mpfr_t m, mpfr_t b1, mpfr_t b2) {
    mpv_getc(c, lm->points, ind);
    
    do {
        mpc_sub(bt, t, pt);                       // bt = t - pt
        mpc_mod(m, bt);                           // m = |bt| = |t - pt|
        mpfr_sub_si(b1, pr, 4, MPFR_RNDN);        // b1 = pr - 4
        mpfr_div_d(b1, b1, lm->guard, MPFR_RNDN); // b1 = (pr - 4) / guard
        
        if(mpfr_cmp(b1, m) > 0) {
            if(mandel_miss_root_ref(c, c, t, lm->prePeriod, lm->period, LEVM_MAX_ITER, 2)) {
                mpv_setc(lm->points, ind, c);
                
                return true;
            } else {
                return false;
            }
        }
        
        mpfr_div(b1, b1, m, MPFR_RNDN);
        mpc_mulr(bt, bt, b1);
        mpc_add(bt, bt, pt);    // here bt = pt + (t - pt) * (pr - 4) / (lm->guard * |t - pt|)
        if(! mandel_miss_root_ref(c, c, bt, lm->prePeriod, lm->period, LEVM_MAX_ITER, 2)) {
            return false;
        }
        
        mpc_set(pt, bt);
        mpc_mod(pr, pt);
    } while(true);
}

/// @brief Lifts the curve @c lm to level @c r.
///
/// @param lm the level curve to refine in place
/// @param r the new level of the curve
///
/// @return @ref true if successfull, @ref false otherwise
static bool levm_level_direct(levm lm, mpfr_t r) {
    uint prec = lm->prec;
    uint tpow = lm->ang2pow;
    ulong count = levm_count(lm);
    
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
    
    int ord = lm->prePeriod + lm->period;
    bool ok = true;
    if(ord > tpow) { // all values of p_n are r
        mpc_setr(t, r);
        mpfr_set(pr, lm->radius, MPFR_RNDN);
        mpc_setr(pt, pr);
        
        for (ulong i = 0; i < count && ok; i++) {
            ok = levm_move_point(lm, i, pt, pr, t, c, b, m, b1, b2);
        }
    } else { // need a circle for targets, faster than recomputing each time
        mpfr_t mr;
        mpfr_init2(mr, prec);
        mpfr_set(mr, lm->radius, MPFR_RNDN);
        mpfr_div(mr, mr, r, MPFR_RNDN);
        
        int lcpm1 = ord - 1;
        levc cr = levc_circle(tpow - lcpm1, prec, r);
        ulong mask = (1L << (tpow - lcpm1)) - 1;
        
        for (ulong i = 0; i < count && ok; i++) {
            mpv_getc(t, cr->points, i & mask);
            mpc_mulr(pt, t, mr);
            mpfr_set(pr, lm->radius, MPFR_RNDN);
            
            ok = levm_move_point(lm, i, pt, pr, t, c, b, m, b1, b2);
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
        mpfr_set(lm->radius, r, MPFR_RNDN);
    }
    
    return ok;
}

levm levm_level(levm src, mpfr_t r) {
    if(src == NULL || r == NULL ||
       ! mpfr_regular_p(r) || mpfr_cmp_si(r, 4) <= 0) {
        return NULL;
    }
    
    levm lm = levm_copy(src, src->prec, true, true);

    if(levm_level_direct(lm, r)) {
        return lm;
    }
    
    levm_free(lm);
    
    return NULL;
}

/// @brief Moves the points of the curve @c lm from hyperbolic (of period @c prePeriod+period) to
/// the @ref misCurve described by the parameters of @c lm.
///
/// @param lm the level curve to refine in place
///
/// @return @ref true if successfull, @ref false otherwise
static bool levm_from_hyp_direct(levm lm) {
    uint prec = lm->prec;
    uint tpow = lm->ang2pow;
    ulong count = levm_count(lm);
    
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
    
    mpfr_ptr r = lm->radius;
    
    int pp = lm->prePeriod;
    int per = lm->period;
    int n = pp + per;
    mpv pts = lm->points;
    
    bool ok = true;
    if(n > tpow) { // all values of p_n are r
        for (ulong i = 0; i < count && ok; i++) {
            mpc_setr(t, r);
            mpv_getc(c, pts, i);
            mandel_mis_val(pt, c, pp, per);
            mpc_mod(pr, pt);
            
            ok = levm_move_point(lm, i, pt, pr, t, c, b, m, b1, b2);
        }
    } else { // need a circle for targets, faster than recomputing each time
        int lcpm1 = n - 1;
        ulong mask = (1L << lcpm1) - 1;
        levc cr = levc_circle(tpow - lcpm1, prec, r);
        
        for (ulong i = 0; i < count && ok; i++) {
            mpv_getc(t, cr->points, i & mask);
            
            mpv_getc(c, pts, i);;
            mandel_mis_val(pt, c, pp, per);
            mpc_mod(pr, pt);
            
            ok = levm_move_point(lm, i, pt, pr, t, c, b, m, b1, b2);
        }
        
        levc_free(cr);
    }
    
    mpfr_clear(b2);
    mpfr_clear(pr);
    mpfr_clear(b1);
    mpfr_clear(m);
    
    mpc_clear(pt);
    mpc_clear(b);
    mpc_clear(t);
    mpc_clear(c);
    
    return ok;
}

/// @brief Moves the points of the simple curve @c lm from hyperbolic (of period @c prePeriod+period-1) to
/// the @ref misCurve described by the parameters of @c lm.
///
/// Simple curves correspond to the simplified Misiurewicz polynomial @c p_{pp,per}/p_{pp-1,per}.
///
/// @param lm the level curve to refine in place
///
/// @return @ref true if successfull, @ref false otherwise
static bool levm_simple_from_hyp_direct(levm lm) {
    uint prec = lm->prec;
    uint tpow = lm->ang2pow;
    ulong count = levm_count(lm);
    
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
    
    mpfr_ptr r = lm->radius;
    
    int pp = lm->prePeriod;
    int per = lm->period;
    int n = pp + per - 1;
    mpv pts = lm->points;
    
    bool ok = true;
    if(n > tpow) { // all values of p_n are r
        for (ulong i = 0; i < count && ok; i++) {
            mpc_setr(t, r);
            mpv_getc(c, pts, i);
            
            mandel_mis_val(pt, c, pp, per);
            mandel_mis_val(b, c, pp - 1, per);
            mpc_div(pt, pt, b);
            
            mpc_mod(pr, pt);
            
            ok = levm_move_simple_point(lm, i, pt, pr, t, c, b, m, b1, b2);
        }
    } else { // need a circle for targets, faster than recomputing each time
        int lcpm1 = n - 1;
        levc cr = levc_circle(tpow - lcpm1, prec, r);
        
        for (ulong i = 0; i < count && ok; i++) {
            mpv_getc(t, cr->points, i >> lcpm1);
            
            mpv_getc(c, pts, i);
            
            mandel_mis_val(pt, c, pp, per);
            mandel_mis_val(b, c, pp - 1, per);
            mpc_div(pt, pt, b);
            
            mpc_mod(pr, pt);
            
            ok = levm_move_simple_point(lm, i, pt, pr, t, c, b, m, b1, b2);
        }
        
        levc_free(cr);
    }
    
    mpfr_clear(b2);
    mpfr_clear(pr);
    mpfr_clear(b1);
    mpfr_clear(m);
    
    mpc_clear(pt);
    mpc_clear(b);
    mpc_clear(t);
    mpc_clear(c);
    
    return ok;
}

levm levm_from_hyp(levc src, uint pp, uint per, double guard) {
    int n = pp + per;
    if(src == NULL || pp < 2 || per < 1 || src->period != n || n > LEVM_MAX_PERIOD) {
        return NULL;
    }
    
    levm lm = levm_alloc(pp, per, src->ang2pow, src->prec, src->startAngle,
                         src->endAngle, src->radius);
    lm->guard = guard;
    mpv_copy(lm->points, 0, src->points);

    if(levm_from_hyp_direct(lm)) {
        return lm;
    }
    
    levm_free(lm);
    
    return NULL;
}

levm levm_simple_from_hyp(levc src, uint pp, uint per, double guard) {
    int n = pp + per - 1;
    if(src == NULL || pp < 2 || per < 1 || n > LEVM_MAX_PERIOD) {
        return NULL;
    }
    
    levm lm = levm_alloc(pp, per, src->ang2pow, src->prec, src->startAngle,
                         src->endAngle, src->radius);
    lm->guard = guard;
    mpv_copy(lm->points, 0, src->points);

    if(levm_simple_from_hyp_direct(lm)) {
        return lm;
    }
    
    levm_free(lm);
    
    return NULL;
}

// MARK: IO functions

bool levm_write_csv(levm lm, char *fileName, uint digits) {
    if(lm == NULL || fileName == NULL || digits < 10 || digits > 1500) {
        return false;
    }
    
    FILE *f = fopen(fileName, "w");
    if(f == NULL) {
        return false;
    }
    
    if(fprintf(f, "%u %u %u %u\n", lm->prePeriod, lm->period, lm->ang2pow, lm->prec) < 6) {
        fclose(f);
        
        return false;
    }
    
    if(fprintf(f, "%lu %lu\n", lm->startAngle, lm->endAngle) < 4) {
        fclose(f);
        
        return false;
    }
    
    char fmt[30];
    snprintf(fmt, 29, "%%.%dRf\n", digits);
    int len = 2 * digits + 200;
    char line[len + 1];
    int minLen = 2 * digits + 3;
    
    int n = mpfr_snprintf(line, len, fmt, lm->radius);
    if(n < minLen || n >= len || fprintf(f, "%s\n", line) < n) {
        fclose(f);
        
        return false;
    }
    
    if(fprintf(f, "%.6lf\n", lm->guard) < 10) {
        fclose(f);
        
        return false;
    }
    
    snprintf(fmt, 29, "%%.%dRf, %%.%dRf\n", digits, digits);
    ulong count = levm_count(lm);
    mpc c;
    mpc_init(c, lm->prec);
    for (ulong i = 0; i < count; i ++) {
        levm_point(c, lm, i);
        int n = mpfr_sprintf(line, fmt, c->x, c->y);
        if(n < minLen || n >= len || fprintf(f, "%s\n", line) < n) {
            fclose(f);
            
            return false;
        }
    }
    
    fclose(f);
    
    return true;
}

levm levm_read_csv(char *fileName) {
    if(fileName == NULL) {
        return NULL;
    }

    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return NULL;
    }

    int pp, p, a, pr;
    if(fscanf(f, "%u %u %u %u", &pp, &p, &a, &pr) < 4 || pp < 2 || p < 1 || a < 1 || pr < LEVM_MIN_PREC) {
        fclose(f);
        
        return NULL;
    }
    
    ulong sa, ea;
    if(fscanf(f, "%lu %lu", &sa, &ea) < 2 || ea < sa || ea > (1L << a)) {
        fclose(f);
        
        return NULL;
    }
    
    char sx[1600];
    mpfr_t r;
    mpfr_init2(r, pr);
    if(fscanf(f, "%s", sx) || mpfr_set_str(r, sx, 10, MPFR_RNDN) != 0) {
        fclose(f);
        
        return NULL;
    }
    
    levm lm = levm_alloc(pp, p, a, pr, sa, ea, r);
    mpfr_clear(r);
    
    if(lm == NULL) {
        levm_free(lm);
        fclose(f);
        
        return NULL;
    }
    
    if(fscanf(f, "%lf", &lm->guard) < 1 || lm->guard < LEVM_MIN_GUARD) {
        levm_free(lm);
        fclose(f);
        
        return NULL;
    }
    
    char sy[1600];
    bool ok = true;
    mpc c;
    mpc_init(c, pr);
    ulong count = levm_count(lm);
    for (ulong i = 0; i < count && ok; i ++) {
        if(fscanf(f, "%s%s", sx, sy) == 2) {
            long lx = strlen(sx);
            if(lx > 1 && sx[lx - 1] == ',') {
                sx[lx - 1] = 0;
            }
            
            bool ok = 0 == mpfr_set_str(c->x, sx, 10, MPFR_RNDN);
            ok = ok && 0 == mpfr_set_str(c->y, sy, 10, MPFR_RNDN);
            ok = ok && mpv_setc(lm->points, i, c);
        } else {
            ok = false;
        }
    }

    mpc_clear(c);
    fclose(f);

    if(! ok) {
        levm_free(lm);
        
        return NULL;
    }
    
    return lm;
}

bool levm_write_to(levm lm, char *fileName) {
    if(fileName == NULL) {
        return false;
    }

    FILE *f = fopen(fileName, "w");
    ulong wr = levm_write(lm, f, 0);
    fclose(f);
        
    return wr > 0 && wr == levm_fileLen(lm);
}

levm levm_read_from(char *fileName) {
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return false;
    }
    
    levm lm = levm_read(f, 0);
    
    fclose(f);
    
    return lm;
}

/// @brief Writes the file header of the curve @c lm to the file @c f, starting at absolute position @c pos if @c pos>=0,
/// at the current position otherwise.
///
/// @param lm the level curve
/// @param f the file
/// @param pos the absolute position if @c pos>=0, at the current position otherwise
///
/// @return @ref true if the operation completed successfully, @ref false otherwise
static bool levm_write_header(levm lm, FILE *f, long pos) {
    if(f == NULL) { // v already checked
        return false;
    }
        
    ulong flen = levm_fileLen(lm);
    mfile h = mfile_header(LEVM_FILE_ID, flen, levm_headerLen(lm));
    if(h == NULL) {
        return false;
    }
    
    mfile_puti(h, lm->prePeriod);
    mfile_puti(h, lm->period);
    mfile_puti(h, lm->ang2pow);
    mfile_puti(h, lm->prec);
    
    mfile_putl(h, lm->startAngle);
    mfile_putl(h, lm->endAngle);
    
    mwrite_mpfr(lm->radius, h);
    mfile_putd(h, lm->guard);
    
    bool ok = h->len == levm_headerLen(lm) && mfile_write_to(f, pos, h, 0, h->len);
    
    mfile_free(h);
    
    return ok;
}

ulong levm_write(levm lm, FILE *f, long pos) {
    if(lm == NULL) {
        return 0;
    }
    
    if(! levm_write_header(lm, f, pos)) {
        return 0;
    }
    
    return levm_headerLen(lm) + mpv_write_points(lm->points, f, -1);
}

/// @brief Reads the header into a buffer misCurve.
///
/// The returned object is not fully formed, do @b not use levm_free() to destroy it.
///
/// @param f the file
/// @param pos the absolute position in the file, unless @c pos<0, in which case the current position is used
///
/// @return a misCurve that represents the header of the sored level curve
static levm levm_read_header(FILE *f, long pos) {
    mfile h = mfile_read_header(f, pos, LEVM_FILE_ID, LEVM_MIN_HEADER_LEN);
    if(h == NULL) {
        return NULL;
    }
    
    levm buf = malloc(sizeof(misCurve));
    
    h->pos = HEADER_MIN_LEN; // skip de file ID, length of the file and of the header
    
    buf->prePeriod = mfile_geti(h);
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
    
    buf->guard = mfile_getd(h);
    
    long mpos = h->pos;
    mfile_free(h);
    
    if(mpos != levm_headerLen(buf)) {
        mpfr_clear(buf->radius);
        free(buf);
        
        return NULL;
    }
    
    return buf;
}

levm levm_read(FILE *f, long pos) {
    levm buf = levm_read_header(f, pos);
    
    if(buf == NULL) {
        return NULL;
    }
            
    // create the new curve
    levm lm = levm_alloc(buf->prePeriod, buf->period, buf->ang2pow, buf->prec, buf->startAngle,
                         buf->endAngle, buf->radius);
    lm->guard = buf->guard;

    mpfr_clear(buf->radius);
    free(buf);
    
    if(lm == NULL || ! mpv_read_points(lm->points, 0, lm->points->count, f, -1)) {
        levm_free(lm);
        
        return NULL;
    }
    
    return lm;
}

levm levm_read_partial(FILE *f, long pos, ulong start, ulong step, ulong count) {
    levm buf = levm_read_header(f, pos);
            
    // check the parameters
    int l2 = levm_log2(step);
    ulong end = buf->startAngle + start + step * (count - 1);
    if(l2 < 0 || count < 2 || start % step != 0 || end > buf->endAngle) {
        free(buf);
        
        return NULL;
    }
    
    // create the new curve
    levm lm = levm_alloc(buf->prePeriod, buf->period, buf->ang2pow - l2, buf->prec, buf->startAngle >> l2,
                         buf->endAngle >> l2, buf->radius);
    lm->guard = buf->guard;

    mpfr_clear(buf->radius);
    free(buf);
            
    ulong tc = levm_count(buf);
    if(lm == NULL || ! mpv_read_points_partial(lm->points, tc, f, -1, start, step, count, true)) {
        levm_free(lm);
        
        return NULL;
    }
    
    return lm;
}

