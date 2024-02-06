//
//  mpDisk.c
//  Mandelbrot
//
//  x, y are considered exact, their eventual errors are contained in r (considered exact at each time, too)
//  m is an approximation of the modulus, care should be exerced when used; computed lazily
//
//  the two or three arguments of functions are not necessarily distinct; this imposes the use of only one
//  set of buffers and the attribution of the final values only when the operands are not used anymore
//
//  Created by MIHALACHE Nicolae on 12/20/19.
//  Copyright © 2019 UPEC. All rights reserved.
//

#include <math.h>
#include <mpfr.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "mpc.h"
#include "fp80d.h"

#include "mpd.h"

bool mpd_init(mpd c, long prec, long precR) {
    if(c == NULL || prec < MP_MIN_PREC || precR < MP_MIN_PREC) {
        return false;
    }
    
    mpfr_init2(c->x, prec);
    mpfr_init2(c->y, prec);
    mpfr_init2(c->mu, prec);
    mpfr_init2(c->md, prec);
        
    mpfr_init2(c->r, precR);
    
    return true;
}

bool mpd_set_prec(mpd c, long prec, long precR) {
    if(c == NULL || prec < MP_MIN_PREC || precR < MP_MIN_PREC) {
        return false;
    }
    
    if(prec != mpc_prec(c)) {
        mpfr_set_prec(c->x, prec);
        mpfr_set_prec(c->y, prec);
        mpfr_set_prec(c->mu, prec);
        mpfr_set_prec(c->md, prec);
    }
        
    if(precR != mpfr_get_prec(c->r)) {
        mpfr_set_prec(c->r, precR);
    }
    
    return true;
}

bool mpd_iniz(mpd d, long prec, long precR, ulong *limbs) {
    bool ok = mpfr_initz(prec, limbs, d->x, d->y, d->mu, d->md, NULL) == 4;
    ok = ok && mpfr_iniz(d->r, precR, limbs + mpfr_limbs(prec) * 4) == 1;
    
    return ok;
}

int mpd_initz(long prec, long precR, ulong *limbs, mpd d, ...) {
    va_list ap;
    va_start(ap, d);
    
    bool ok = mpd_iniz(d, prec, precR, limbs);
    int count = 0;
    long size = mpd_limbs(prec, precR);
    
    mpd_struct *a;
    while(ok) {
        count ++;
        
        a = va_arg(ap, mpd_struct *);
        ok = mpd_iniz(a, prec, precR, limbs + count * size);
    }
    
    va_end(ap);
    
    return count;
}

bool mpd_clear(mpd c) {
    if(c == NULL) {
        return false;
    }
    
    mpfr_clear(c->x);
    mpfr_clear(c->y);
    mpfr_clear(c->mu);
    mpfr_clear(c->md);
    
    mpfr_clear(c->r);
    
    return true;
}

int mpd_clears(mpd d, ...) {
    va_list ap;
    va_start(ap, d);
    
    bool ok = mpd_clear(d);
    int count = 0;
    
    mpd_struct *a;
    while(ok) {
        count ++;
        
        a = va_arg(ap, mpd_struct *);
        ok = mpd_clear(a);
    }
    
    va_end(ap);
    
    return count;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Assignment from s (source, various types) into d (destination)
// precision (higher or lower) is imposed by the destination (round to nearest if necessary)
// modulus is indicative only (~buffer) and is set to NaN
// //////////////////////////////////////////////////////////////////////////////////////////

void mpd_set(mpd d, mpd s) {
    mpfr_set(d->x, s->x, MPFR_RNDN);
    mpfr_set(d->y, s->y, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    mpfr_set(d->r, s->r, MPFR_RNDU);
}

void mpd_set80(mpd d, fp80d s) {
    mpfr_set_ld(d->x, s->x, MPFR_RNDN);
    mpfr_set_ld(d->y, s->y, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    mpfr_set_ld(d->r, s->r, MPFR_RNDU);
}

void mpd_get80(fp80d d, mpd s) {
    d->x = mpfr_get_ld(s->x, MPFR_RNDN);
    d->y = mpfr_get_ld(s->y, MPFR_RNDN);
    d->r = mpfr_get_ld(s->r, MPFR_RNDN);
    d->m = hypotl(d->x, d->y);
}

void mpd_setr(mpd d, mpc c, mpfr_t r) {
    mpfr_set(d->x, c->x, MPFR_RNDN);
    mpfr_set(d->y, c->y, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    mpfr_set(d->r, r, MPFR_RNDU);
}

void mpd_setrl(mpd d, mpc c, long double r) {
    mpfr_set(d->x, c->x, MPFR_RNDN);
    mpfr_set(d->y, c->y, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    mpfr_set_ld(d->r, r, MPFR_RNDU);
}

void mpd_set80r(mpd d, fp80 c, long double r) {
    mpfr_set_ld(d->x, c->x, MPFR_RNDN);
    mpfr_set_ld(d->y, c->y, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    mpfr_set_ld(d->r, r, MPFR_RNDU);
}

void mpd_setd(mpd d, double x, double y, double r) {
    mpfr_set_d(d->x, x, MPFR_RNDN);
    mpfr_set_d(d->y, y, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    mpfr_set_d(d->r, r, MPFR_RNDU);
}

void mpd_setld(mpd d, long double x, long double y, long double r) {
    mpfr_set_ld(d->x, x, MPFR_RNDN);
    mpfr_set_ld(d->y, y, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    mpfr_set_ld(d->r, r, MPFR_RNDU);
}

void mpd_set_exact(mpd d, mpfr_t x, mpfr_t y) { // unspecified radius defaults to zero
    mpfr_set(d->x, x, MPFR_RNDN);
    mpfr_set(d->y, y, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    mpfr_set_ui(d->r, 0, MPFR_RNDN);
}

void mpd_set_center(mpd d, mpd a) { // equivalent to mpDisk_setExact(d, a->x, a->y)
    mpfr_set(d->x, a->x, MPFR_RNDN);
    mpfr_set(d->y, a->y, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    mpfr_set_ui(d->r, 0, MPFR_RNDN);
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: ULP = Unit Last Position
// For x != 0, one defines ulp(x) == 2^( exponent(x) - precision(x)  )
// It is the least amount by which x can be changed when the least significant bit is changed
// Our convention when x == 0 is :
//     * if no previous knowledge on x, then ulp(0) == 2^( - x->_mpfr_prec )
//     * if x is the result of a computation with prevUlp, then ulp(x) == prevUlp
// This convention conflicts with MPFR handling of zero: 0->_mpfr_exp == -(2^63 - 1)
// We therefore need to correct the exponent of zero when necessary.
// //////////////////////////////////////////////////////////////////////////////////////////

void mpd_ulp(mpfr_t u, mpfr_t x, long prevLogUlp) { // u = ulp(x) with special zero handling
    // usage: prevLogUlp is the log2 of the ulp immediately preceeding the operation that led to x == 0
    if(mpfr_zero_p(x)) {
        mpfr_set_si_2exp(u, 1, prevLogUlp, MPFR_RNDU);
        return;
    }
    
    mpfr_set_si_2exp(u, 1, x->_mpfr_exp - x->_mpfr_prec, MPFR_RNDU);
}

void mpd_hulp(mpfr_t u, mpfr_t x, long prevLogUlp) { // half ulp : u = mpd_ulp(u, x, prevExp) / 2
    if(mpfr_zero_p(x)) {
        mpfr_set_si_2exp(u, 1, prevLogUlp - 1, MPFR_RNDU);
        return;
    }
    
    mpfr_set_si_2exp(u, 1, x->_mpfr_exp - x->_mpfr_prec - 1, MPFR_RNDU);
}

static inline long mpd_ulpp(mpfr_t x) { // tool: log2( ulp(x) ) == exponent(x) - precision(x)
                                        // corrected such that exponent(0) = precision(0) instead of -(2^63 - 1)
                                        // The last p stands for position i.e. log2
    return (mpfr_zero_p(x) ? 0 : x->_mpfr_exp) - x->_mpfr_prec;
}

static long mpd_ulp_pos(mpfr_t x, mpfr_t y) { // tool: minimal log2 of ulp
    // if one parameter is zero, then the ulp of the other one is used
    // if both are zero, then returns - precision(y)
    if(mpfr_zero_p(x))
        return mpd_ulpp(y);  // x == 0 : corrected log2(ulp(y))
    
    if(mpfr_zero_p(y))
        return x->_mpfr_exp - x->_mpfr_prec; // x != 0 && y == 0 : log2(ulp(x))
    
    long ex = x->_mpfr_exp - x->_mpfr_prec;
    long ey = y->_mpfr_exp - y->_mpfr_prec;
    
    return ex <= ey ? ex : ey; // min{ log2(ulp(x)), log2(ulp(y)) }
}

void mpd_set_ulp(mpd d, mpfr_t x, mpfr_t y) { // unspecified radius defaults to rectangular ULP
    defs_mpfr(mpd_prec(d) + MPD_EXTRA_PREC, br1, br2);
    
    mpfr_set(d->x, x, MPFR_RNDN);
    mpfr_set(d->y, y, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    
    mpd_ulp(br1, d->x, mpd_ulpp(x));
    mpd_ulp(br2, d->y, mpd_ulpp(y));
    
    mpfr_hypot(d->r, br1, br2, MPFR_RNDU);
}

void mpd_set_half_ulp(mpd d, mpfr_t x, mpfr_t y) { // unspecified radius defaults to half rectangular ULP
    defs_mpfr(mpd_prec(d) + MPD_EXTRA_PREC, br1, br2);
    
    mpfr_set(d->x, x, MPFR_RNDN);
    mpfr_set(d->y, y, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    
    mpd_hulp(br1, d->x, mpd_ulpp(x));
    mpd_hulp(br2, d->y, mpd_ulpp(y));
    
    mpfr_hypot(d->r, br1, br2, MPFR_RNDU);
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Interval arithmetic tools
// Specific tools to compute real and imaginary parts in complex multiplications
// //////////////////////////////////////////////////////////////////////////////////////////

/// @brief Computes @c a*b-c*d and writes the result in @c res and a bound for all errors in @c err.
///
/// @warning the buffers @c b1, @c b2, @c br1, @c br2 , @c pb1 and @c pb2 of @c buff are changed !
///
/// @param res the result
/// @param err the bound for rounding errors
/// @param a first operand
/// @param b second operand
/// @param c third operand
/// @param d fourth operand
static void mpd_mms(mpfr_t res, mpfr_t err, mpfr_t a, mpfr_t b, mpfr_t c, mpfr_t d,
                    ulong prec, ulong rprec, bool highPrec) {
    // Interval arithmetic tool:    mms = multiply, multiply, substract
    // Compute res <- ab - cd and guaranties that the exact value lies in [res-err, res+err].
    // Uses b1, b2, br1 and br2 as buffers; in the case of precision loss, it uses the precise buffers pb1 and pb2

    // check if one of the terms is zero, with several advantages
    // 1. better error bounds in this case
    // 2. much faster when working with real or purely imaginary numbers
    int za = mpfr_zero_p(a);
    int zb = mpfr_zero_p(b);
    int zc = mpfr_zero_p(c);
    int zd = mpfr_zero_p(d);
    
    if(za || zb) { // ab == 0
        if(zc || zd) { // cd == 0
            mpfr_set_zero(res, 1);
            mpfr_set_zero(err, 1);
            
            return;
        }
        
        mpfr_mul(res, c, d, MPFR_RNDN);
        mpfr_neg(res, res, MPFR_RNDN);
        mpd_hulp(err, res, mpd_ulpp(res));
        
        return;
    }
    
    if(zc || zd) {
        mpfr_mul(res, a, b, MPFR_RNDN);
        mpd_hulp(err, res, mpd_ulpp(res));
        
        return;
    }
    
    defs_mpfr(prec + MPD_EXTRA_PREC, b1, b2);
    defs2_mpfr(rprec + MPD_EXTRA_PREC, br1, br2);
    
    mpfr_mul(b1, a, b, MPFR_RNDN);
    mpd_hulp(br1, b1, mpd_ulpp(b1));    // ab in [b1 - br1, b1 + br1]
    
    mpfr_mul(b2, c, d, MPFR_RNDN);
    mpd_hulp(br2, b2, mpd_ulpp(b2));    // cd in [b2 - br2, b2 + br2]
    
    long logE = mpd_ulp_pos(b1, b2);
    if(mpfr_sub(res, b1, b2, MPFR_RNDN)) {       // compute res == ab - cd
        // if res is not exact, collect all possible errors
        mpfr_add(br1, br1, br2, MPFR_RNDU); // br1 <- br1 + br2
        mpd_hulp(br2, res, logE); // error from rounding up res to final precision
    }
    mpfr_add(err, br1, br2, MPFR_RNDU);          // ab - cd in [res - err, res + err]
    
    if(! highPrec) {
        return;
    }
    
    // If the loss of exponent is larger than MPD_EXTRA_PREC in the substraction
    // then we could do better by using high precision buffers
    // For the products, exact results are expected if precH >= 2 * prec.
    long abE = b1->_mpfr_exp;
    long cdE = b2->_mpfr_exp;
    long resE = res->_mpfr_exp + MPD_EXTRA_PREC;
    long l1 = abE - resE, l2 = cdE - resE, loss = l1 > l2 ? l1 : l2;
    
    if(loss > 0 || mpfr_zero_p(res)) {
        defs_mpfr(prec + loss + MPD_EXTRA_PREC, pb1, pb2);
        
        if(! mpfr_mul(pb1, a, b, MPFR_RNDN)) { // exact result a * b
            mpfr_set_zero(br1, 0);
        } else { // rounding error in a * b product :  a * b in [pb1-br1, pb1+br1]
            mpd_hulp(br1, pb1, mpd_ulpp(pb1));
        }
        
        if(! mpfr_mul(pb2, c, d, MPFR_RNDN)) { // exact result c * d
            mpfr_set_zero(br2, 0);
        } else { // rounding error in c * d product : c * d in [pb2-br2, pb2+br2]
            mpd_hulp(br2, pb2, mpd_ulpp(pb2));
        }
        
        logE = mpd_ulp_pos(pb1, pb2);
        if(mpfr_sub(res, pb1, pb2, MPFR_RNDN)) { // compute res == ab - cd
            // if res is not exact, collect all possible errors
            mpfr_add(br1, br1, br2, MPFR_RNDU); // br1 <- br1+br2
            mpd_hulp(br2, res, logE); // error from rounding up res to final precision
        }
        mpfr_add(err, br1, br2, MPFR_RNDU);     // ab-cd in [res-err, res+err]
    }
}

/// @brief Computes @c a*b+c*d and writes the result in @c res and a bound for all errors in @c err.
///
/// @warning the buffers @c b1, @c b2, @c br1, @c br2 , @c pb1 and @c pb2 of @c buff are changed !
///
/// @param res the result
/// @param err the bound for rounding errors
/// @param a first operand
/// @param b second operand
/// @param c third operand
/// @param d fourth operand
static void mpd_mma(mpfr_t res, mpfr_t err, mpfr_t a, mpfr_t b, mpfr_t c, mpfr_t d,
                    ulong prec, ulong rprec, bool highPrec) {
    // Interval arithmetic tool:    mma = multiply, multiply, add
    // Compute res <- ab + cd and guaranties that the exact value lies in [res-err, res+err].
    // Uses b1, b2, br1 and br2 as buffers; in the case of precision loss, it uses the precise buffers pb1 and pb2
    
    // check if one of the terms is zero, with several advantages
    // 1. better error bounds in this case
    // 2. much faster when working with real or purely imaginary numbers
    int za = mpfr_zero_p(a);
    int zb = mpfr_zero_p(b);
    int zc = mpfr_zero_p(c);
    int zd = mpfr_zero_p(d);
    
    if(za || zb) { // ab == 0
        if(zc || zd) { // cd == 0
            mpfr_set_zero(res, 1);
            mpfr_set_zero(err, 1);
            
            return;
        }
        
        mpfr_mul(res, c, d, MPFR_RNDN);
        mpd_hulp(err, res, mpd_ulpp(res));
        
        return;
    }
    
    if(zc || zd) {
        mpfr_mul(res, a, b, MPFR_RNDN);
        mpd_hulp(err, res, mpd_ulpp(res));
        
        return;
    }
    
    defs_mpfr(prec + MPD_EXTRA_PREC, b1, b2);
    defs2_mpfr(rprec + MPD_EXTRA_PREC, br1, br2);
    
    mpfr_mul(b1, a, b, MPFR_RNDN);
    mpd_hulp(br1, b1, mpd_ulpp(b1));    // a * b in [b1-br1, b1+br1]
    
    mpfr_mul(b2, c, d, MPFR_RNDN);
    mpd_hulp(br2, b2, mpd_ulpp(b2));    // c * d in [b2-br2, b2+br2]
    
    long logE = mpd_ulp_pos(b1, b2);
    if(mpfr_add(res, b1, b2, MPFR_RNDN)) {       // compute res == ab + cd
        // if res is not exact, collect all possible errors
        mpfr_add(br1, br1, br2, MPFR_RNDU); // br1 <- br1+br2
        mpd_hulp(br2, res, logE); // error from rounding up res to final precision
    }
    mpfr_add(err, br1, br2, MPFR_RNDU);          // ab+cd in [res-err, res+err]
    
    if(! highPrec) {
        return;
    }
    
    // If the loss of exponent is larger than MP_DISK_EXTRA_PREC in the substraction
    // then we could do better by using the high precision buffers, if availlable.
    // For the products, exact results are normally expected if precH >= 2 * prec.
    long abE = b1->_mpfr_exp;
    long cdE = b2->_mpfr_exp;
    long resE = res->_mpfr_exp + MPD_EXTRA_PREC;
    long l1 = abE - resE, l2 = cdE - resE, loss = l1 > l2 ? l1 : l2;
    
    if(loss > 0 || mpfr_zero_p(res)) {
        defs_mpfr(prec + loss + MPD_EXTRA_PREC, pb1, pb2);
        
        // br1, br2 are reused, old values are discarded
        if(! mpfr_mul(pb1, a, b, MPFR_RNDN)) { // exact result a * b
            mpfr_set_zero(br1, 0);
        } else { // rounding error in a * b product :  a * b in [pb1-br1, pb1+br1]
            mpd_hulp(br1, pb1, mpd_ulpp(pb1));
        }
        
        if(! mpfr_mul(pb2, c, d, MPFR_RNDN)) { // exact result c * d
            mpfr_set_zero(br2, 0);
        } else { // rounding error in c * d product : c * d in [pb2-br2, pb2+br2]
            mpd_hulp(br2, pb2, mpd_ulpp(pb2));
        }
        
        logE = mpd_ulp_pos(pb1, pb2);
        if(mpfr_add(res, pb1, pb2, MPFR_RNDN)) { // compute res == ab + cd
            // if res is not exact, collect all possible errors
            mpfr_add(br1, br1, br2, MPFR_RNDU); // br1 <- br1+br2
            mpd_hulp(br2, res, logE); // error from rounding up res to final precision
        }
        mpfr_add(err, br1, br2, MPFR_RNDU);     // ab+cd in [res-err, res+err]
    }
}

static void mpd_sqs(mpfr_t res, mpfr_t err, mpfr_t a, mpfr_t b,
                    ulong prec, ulong rprec, bool highPrec) {
    // Interval arithmetic tool:    sqs = square, substract
    // Compute res <- a^2 - b^2 and guaranties that the exact value lies in [res-err, res+err].
    // Equivalent to mpd_mms(res, err, a, a, b, b, buff) but optimized.
    // Uses b1, b2, br1 and br2 as buffers; in the case of precision loss, it uses the precise buffers pb1 and pb2

    // check if one of the terms is zero, with several advantages
    // 1. better error bounds in this case
    // 2. much faster when working with real or purely imaginary numbers
    int za = mpfr_zero_p(a);
    int zb = mpfr_zero_p(b);
    
    if(za) { // a == 0
        if(zb) { // b == 0
            mpfr_set_zero(res, 1);
            mpfr_set_zero(err, 1);
            
            return;
        }
        
        mpfr_sqr(res, b, MPFR_RNDN);
        mpfr_neg(res, res, MPFR_RNDN);
        mpd_hulp(err, res, mpd_ulpp(res));
        
        return;
    }
    
    if(zb) {
        mpfr_sqr(res, a, MPFR_RNDN);
        mpd_hulp(err, res, mpd_ulpp(res));
        
        return;
    }
    
    // checks if a == b or a == -b, in which case returns 0 with error 0
    if(! mpfr_cmp_abs(a, b)) {
        mpfr_set_zero(res, 1);
        mpfr_set_zero(err, 1);
        
        return;
    }
    
    defs_mpfr(prec + MPD_EXTRA_PREC, b1, b2);
    defs2_mpfr(rprec + MPD_EXTRA_PREC, br1, br2);
    
    mpfr_sqr(b1, a, MPFR_RNDN);
    mpd_hulp(br1, b1, mpd_ulpp(b1));     // a^2 in [b1 - br1, b1 + br1]
    
    mpfr_sqr(b2, b, MPFR_RNDN);
    mpd_hulp(br2, b2, mpd_ulpp(b2));     // b^2 in [b2 - br2, b2 + br2]
    
    long logE = mpd_ulp_pos(b1, b2);
    if(mpfr_sub(res, b1, b2, MPFR_RNDN)) {        // compute res == a^2 - b^2
        // if res is not exact, collect all possible errors
        mpfr_add(br1, br1, br2, MPFR_RNDU); // br1 <- br1 + br2
        mpd_hulp(br2, res, logE); // error from rounding up res to final precision
    }
    mpfr_add(err, br1, br2, MPFR_RNDU);           // a^2 - b^2 in [res - err, res + err]
    
    if(! highPrec) {
        return;
    }
    
    // If the loss of exponent is larger than MP_DISK_EXTRA_PREC in the substraction
    // then we could do better by using the high precision buffers, if availlable.
    // For the products, exact results are normally expected if precH >= 2 * prec.
    long abE = b1->_mpfr_exp;
    long cdE = b2->_mpfr_exp;
    long resE = res->_mpfr_exp + MPD_EXTRA_PREC;
    long l1 = abE - resE, l2 = cdE - resE, loss = l1 > l2 ? l1 : l2;
    
    if(loss > 0 || mpfr_zero_p(res)) {
        defs_mpfr(prec + loss + MPD_EXTRA_PREC, pb1, pb2);
        
        if(! mpfr_sqr(pb1, a, MPFR_RNDN)) { // exact result a^2
            mpfr_set_zero(br1, 0);
        } else { // rounding error in a^2 :  a^2 in [pb1-br1, pb1+br1]
            mpd_hulp(br1, pb1, mpd_ulpp(pb1));
        }
        
        if(! mpfr_sqr(pb2, b, MPFR_RNDN)) { // exact result b^2
            mpfr_set_zero(br2, 0);
        } else { // rounding error in b^2 :  b^2 in [pb2-br2, pb2+br2]
            mpd_hulp(br2, pb2, mpd_ulpp(pb2));
        }
        
        logE = mpd_ulp_pos(pb1, pb2);
        if(mpfr_sub(res, pb1, pb2, MPFR_RNDN)) { // compute res == a^2 - b^2
            // if not exact, collect all possible errors
            mpfr_add(br1, br1, br2, MPFR_RNDU);   // br1 <- br1 + br2
            mpd_hulp(br2, res, logE); // error from rounding up res to final precision
        }
        mpfr_add(err, br1, br2, MPFR_RNDU);       // a^2-b^2 in [res-err, res+err]
    }
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Basic arithmetic operations with mpDisk(s)
// A zero return code indicates of mpfr_add,... indicates that the computation is exact.
// If the arithmetic that computes the new center is not exact, then one needs to adjust the radius up
// so that one will be allowed to assume that the new center is again exact. The adjustment needed
// depends on the nature of the operation.
// As we always round x and y to nearest, an adjustement of 1/2 of (rectangular) ulp is always mandatory.
// //////////////////////////////////////////////////////////////////////////////////////////

bool mpd_valid(mpd a) {
    // A valid disk has numeric coordinates, a positive (or zero) radius
    // and, if it is defined, a positive (or zero) modulus.
    // It does not check that the modulus is correct, which has lazy initialization.
    // Operations that change coordinates must either compute the correct modulus again, or set it to NaN.
    
    if(mpfr_number_p(a->r) && mpfr_number_p(a->x) && mpfr_number_p(a->y)) {
        // numeric coordinates and radius (not NaN, not infinity)
        if(! mpfr_zero_p(a->r) && a->r->_mpfr_sign < 0)
            return false; // but strictly negative radius
        
        // NaN modulus is ok, -0 is ok, strictly negative, not ok
        bool muok = mpfr_nan_p(a->mu) || mpfr_zero_p(a->mu) || a->mu->_mpfr_sign >= 0;
        bool mdok = mpfr_nan_p(a->md) || mpfr_zero_p(a->md) || a->md->_mpfr_sign >= 0;
        
        return muok && mdok; // modulus is positive (or NaN)
    }
    
    return false; // coordinate or radius is NaN or infinity
}

void mpd_add(mpd d, mpd a, mpd b) { // d = a + b
    // New radius is the sum of radii of a and b, rounded towards +infinity to include the radius rounding error
    // If the arithmetic that computes the new center is not exact, then one needs to adjust the radius up again
    // by a half-ulp diagonal, so that one will be allowed to assume that the new center is again exact.
    // If the new center happens to have a zero cordinate, then the ulp that predates the computation is used.
    // As MPFR guarantees the result up to the last bit: p + q == 0 <=> p == -q and thus ulp(p) = ulp(q).

    long logUlpx = mpd_ulp_pos(a->x, b->x);
    int xex = ! mpfr_add(d->x, a->x, b->x, MPFR_RNDN); // round nearest
    
    long logUlpy = mpd_ulp_pos(a->y, b->y);
    int yex = ! mpfr_add(d->y, a->y, b->y, MPFR_RNDN); // round nearest
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    
    mpfr_add(d->r, a->r, b->r, MPFR_RNDU); // round up
    
    if(xex && yex) { // both operations were exact, no need for radius adjusting
        return;
    }
    
    defs_mpfr(mpd_rad_prec(d) + MPD_EXTRA_PREC, br1, br2);
    
    if(! xex && ! yex) { // both operations were inexact, we need to add a half-ulp diagonal
        // one can use the ulp after the operation because MPFR guarantees the result up to the last bit
        mpd_hulp(br1, d->x, logUlpx);
        mpd_hulp(br2, d->y, logUlpy);
        mpfr_hypot(br1, br1, br2, MPFR_RNDU);
        
        mpfr_add(d->r, d->r, br1, MPFR_RNDU);
        
        return;
    }
    
    if(! xex) { // here only the sum in y coordiante was exact, we add a half-ulp along x
        mpd_hulp(br1, d->x, logUlpx);
        mpfr_add(d->r, d->r, br1, MPFR_RNDU);
        
        return;
    }
    
    // here only the sum in x coordinate was exact, we add a half-ulp along y
    mpd_hulp(br1, d->y, logUlpy);
    mpfr_add(d->r, d->r, br1, MPFR_RNDU);
}

void mpd_addc(mpd d, mpd a, mpc b) {
    long logUlpx = mpd_ulp_pos(a->x, b->x);
    int xex = ! mpfr_add(d->x, a->x, b->x, MPFR_RNDN); // round nearest
    
    long logUlpy = mpd_ulp_pos(a->y, b->y);
    int yex = ! mpfr_add(d->y, a->y, b->y, MPFR_RNDN); // round nearest
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    
    mpfr_set(d->r, a->r, MPFR_RNDU); // round up
    
    if(xex && yex) { // both operations were exact, no need for radius adjusting
        return;
    }
    
    defs_mpfr(mpd_rad_prec(d) + MPD_EXTRA_PREC, br1, br2);
    
    if(! xex && ! yex) { // both operations were inexact, we need to add a half-ulp diagonal
        // one can use the ulp after the operation because MPFR guarantees the result up to the last bit
        mpd_hulp(br1, d->x, logUlpx);
        mpd_hulp(br2, d->y, logUlpy);
        mpfr_hypot(br1, br1, br2, MPFR_RNDU);
        
        mpfr_add(d->r, d->r, br1, MPFR_RNDU);
        
        return;
    }
    
    if(! xex) { // here only the sum in y coordiante was exact, we add a half-ulp along x
        mpd_hulp(br1, d->x, logUlpx);
        mpfr_add(d->r, d->r, br1, MPFR_RNDU);
        
        return;
    }
    
    // here only the sum in x coordinate was exact, we add a half-ulp along y
    mpd_hulp(br1, d->y, logUlpy);
    mpfr_add(d->r, d->r, br1, MPFR_RNDU);
}

void mpd_sub(mpd d, mpd a, mpd b) { // d = a - b
    // New radius is the sum of radii of a and b, rounded towards +infinity to include the radius rounding error
    // If the arithmetic that computes the new center is not exact, then one needs to adjust the radius up again
    // by a half-ulp diagonal, so that one will be allowed to assume that the new center is again exact.
    // If the new center happens to have a zero cordinate, then the ulp that predates the computation is used.
    // As MPFR guarantees the result up to the last bit: p + q == 0 <=> p == -q and thus ulp(p) = ulp(q).

    long logUlpx = mpd_ulp_pos(a->x, b->x);
    int xex = ! mpfr_sub(d->x, a->x, b->x, MPFR_RNDN);  // round nearest
    
    long logUlpy = mpd_ulp_pos(a->y, b->y);
    int yex = ! mpfr_sub(d->y, a->y, b->y, MPFR_RNDN);  // round nearest
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    
    mpfr_add(d->r, a->r, b->r, MPFR_RNDU);  // round up
    
    if(xex && yex) { // both operations were exact, no need for radius adjusting
        return;
    }
    
    defs_mpfr(mpd_rad_prec(d) + MPD_EXTRA_PREC, br1, br2);
    
    if(! xex && ! yex) { // both operations were inexact, we need to add a half-ulp diagonal
        // one can use the ulp after the operation because MPFR guarantees the result up to the last bit
        mpd_hulp(br1, d->x, logUlpx);
        mpd_hulp(br2, d->y, logUlpy);
        mpfr_hypot(br1, br1, br2, MPFR_RNDU);
        
        mpfr_add(d->r, d->r, br1, MPFR_RNDU);
        
        return;
    }
    
    if(! xex) { // here only the sum in y coordiante was exact, we add a half-ulp along x
        mpd_hulp(br1, d->x, logUlpx);
        mpfr_add(d->r, d->r, br1, MPFR_RNDU);
        
        return;
    }
    
    // here only the sum in x coordinate was exact, we add a half-ulp along y
    mpd_hulp(br1, d->y, logUlpy);
    mpfr_add(d->r, d->r, br1, MPFR_RNDU);
}

/// @brief If it is not already computed, initializes @c d->mu with an upper bound for the modulus of the center of @c d.
///
/// It does not change any internal buffer of @c d.
///
/// @param d the disk
///
/// @return a pointer to @c d->mu which stores the upper bound of the center of @c d.
static __mpfr_struct* mpd_mod_up(mpd d) {
    if(mpfr_nan_p(d->mu)) {
        int xz = mpfr_zero_p(d->x);
        int yz = mpfr_zero_p(d->y);
        
        if(xz) {
            if(yz) {
                mpfr_set_zero(d->mu, 1);
            } else {
                mpfr_abs(d->mu, d->y, MPFR_RNDU);
            }
        } if(yz) {
            mpfr_abs(d->mu, d->x, MPFR_RNDU);
        } else {
            mpfr_hypot(d->mu, d->x, d->y, MPFR_RNDU);
        }
    }
    
    return d->mu;
}

/// @brief If it is not already computed, initializes @c d->md with an lower bound for the modulus of the center of @c d.
///
/// It does not change any internal buffer of @c d.
///
/// @param d the disk
///
/// @return a pointer to @c d->md which stores the lower bound of the center of @c d.
static __mpfr_struct* mpd_mod_down(mpd d) {
    if(mpfr_nan_p(d->md)) {
        int xz = mpfr_zero_p(d->x);
        int yz = mpfr_zero_p(d->y);
        
        if(xz) {
            if(yz) {
                mpfr_set_zero(d->md, 1);
            } else {
                mpfr_abs(d->md, d->y, MPFR_RNDD);
            }
        } if(yz) {
            mpfr_abs(d->md, d->x, MPFR_RNDD);
        } else {
            mpfr_hypot(d->md, d->x, d->y, MPFR_RNDD);
        }
    }
    
    return d->md;    
}

void mpd_mul(mpd d, mpd a, mpd b) { // d = a * b
    // New radius is the sum of
    //     (radius of a) * (radius of b) + (modulus of a) * (radius of b)
    //   + (modulus of b) * (radius of a) + rounding errors for the center.
    // Each intermidiary computation of the radius is rounded towards +infinity.
    // This guaranties that one may assume that the new center is again exact.
    
    // allocate 2 coords buffers on the stack
    ulong cp = mpd_prec(d);
    ulong rp = mpd_rad_prec(d);
    defs_mpfr(cp + MPD_EXTRA_PREC, bx, by);
    defs2_mpfr(rp + MPD_EXTRA_PREC, br1, br2, br3);
    
    // compute the real part in bx, rounding errors in br3
    mpd_mms(bx, br3, a->x, b->x, a->y, b->y, cp, rp, true);
    
    // compute the imaginary part in by, rounding errors in br2
    // no conflict with br2 used as buffer by mpd_mma(..., d) !!!
    mpd_mma(by, br2, a->x, b->y, a->y, b->x, cp, rp, true);
    
    // total rounding errors for the center in br3 (and br1, br2 are now free)
    mpfr_hypot(br3, br3, br2, MPFR_RNDU);
    
    mpfr_mul(br1, mpd_mod_up(a), b->r, MPFR_RNDU); // upper bound of (modulus of a) * (radius of b)
    mpfr_mul(br2, mpd_mod_up(b), a->r, MPFR_RNDU); // upper bound of (modulus of b) * (radius of a)
    
    // compute the new radius
    mpfr_mul(d->r, a->r, b->r, MPFR_RNDU);   //  (radius of a) * (radius of b)
    mpfr_add(d->r, d->r, br1, MPFR_RNDU); // (modulus of a) * (radius of b)
    mpfr_add(d->r, d->r, br2, MPFR_RNDU); // (modulus of b) * (radius of a)
    mpfr_add(d->r, d->r, br3, MPFR_RNDU); // rounding errors for the center // now d->r is correct
    
    // must be here, otherwise we may rewrite the centers of a or b
    mpfr_set(d->x, bx, MPFR_RNDN);
    mpfr_set(d->y, by, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
}

bool mpd_is_inf(mpd a) {
    // radius or one of the coordinates is infinite
    return mpfr_inf_p(a->r) || mpfr_inf_p(a->x) || mpfr_inf_p(a->y);
}

static inline void mpd_mul_2p(mpfr_t x, int tpow) {// tool: x = x * 2^tpow, rounded up
    mpfr_mul_2si(x, x, tpow, MPFR_RNDN); // MPFR guaranties that only the exponent is changed
}

//static inline void mpd_abs(mpfr_t x) { // tool: absolute value (faster than mpfr_abs)
//    x->_mpfr_sign = 1;
//}

bool mpd_div(mpd d, mpd a, mpd b) { // d = a / b == a * conj(b) / |b|^2
    if(mpd_contains_0(b, 0)) { // inverse of zero is the whole complex plane
        mpfr_set_zero(d->x, 1);
        mpfr_set_zero(d->y, 1);
        mpfr_set_inf(d->mu, 1);
        mpfr_set_zero(d->md, 1);
        mpfr_set_inf(d->r, 1);
        
        return false; // abnormal termination indicates division by zero
    }
    
    ulong cp = mpd_prec(d);
    ulong rp = mpd_rad_prec(d);
    defs_mpfr(cp + MPD_EXTRA_PREC, bx, by, b1, b2);
    defs2_mpfr(rp + MPD_EXTRA_PREC, br1, br2, br3, br4);
    
    // Step 1 : compute a * conj(b)
    
    // compute the real part in b4, rounding errors in br3
    mpd_mma(bx, br3, a->x, b->x, a->y, b->y, cp, rp, true);
    
    // compute the imaginary part in b5, rounding errors in br2
    mpd_mms(by, br2, a->y, b->x, a->x, b->y, cp, rp, true);
    
    // total rounding errors for the center in br3 (and br1, br2 are now free)
    mpfr_add(br3, br3, br2, MPFR_RNDU);   // hypot too slow, should not lose much here
    
    // Step 2 : compute the new radius; r = a->r, t = b->r
    mpfr_mul(br1, mpd_mod_up(a), b->r, MPFR_RNDU);   // upper bound of |a| * t
    mpfr_mul(br2, mpd_mod_up(b), a->r, MPFR_RNDU);   // upper bound of |b| * r
    mpfr_add(br1, br1, br2, MPFR_RNDU);              // upper bound of |a| * t + |b| * r
    
    mpfr_sub(br2, b->md, b->r, MPFR_RNDD);           // lower bound of |b| - t
    mpfr_mul(br2, br2, b->md, MPFR_RNDD);            // lower bound of |b| * (|b| - t)
    
    mpfr_div(br1, br1, br2, MPFR_RNDU);              // upper bound of (|a| * t + |b| * r) / (|b| * (|b| - t))
            
    // Step 3 : divide by |b|^2 and correct the radius
    
    // compute |b|^2 in b1
    mpfr_sqr(b1, b->x, MPFR_RNDN);
    mpfr_sqr(b2, b->y, MPFR_RNDN);
    mpfr_add(b1, b1, b2, MPFR_RNDN);
    
    // compute (bx, by) = a / b
    mpfr_div(bx, bx, b1, MPFR_RNDN);
    mpfr_div(by, by, b1, MPFR_RNDN);
    
    // br = lower bound of |b|^2
    mpfr_sqr(br2, b->x, MPFR_RNDD);
    mpfr_sqr(br4, b->y, MPFR_RNDD);
    mpfr_add(br2, br2, br4, MPFR_RNDD);
    
    // br3 = max error from rounding (bx, by) = a / b
    mpfr_div(br3, br3, br2, MPFR_RNDU);
    
    // upper bound of the new radius and final result
    mpfr_add(d->r, br1, br3, MPFR_RNDU);
    mpfr_set(d->x, bx, MPFR_RNDN);
    mpfr_set(d->y, by, MPFR_RNDN);
    
    // account for the rounding d = (bx, by)
    if(rp <= cp) { // should be true
        mpfr_nextabove(d->r);
    } else {
        mpfr_abs(bx, bx, MPFR_RNDU);
        mpfr_set(b1, bx, MPFR_RNDD);
        mpfr_nextabove(bx);
        mpfr_sub(br1, bx, b1, MPFR_RNDU);
        
        mpfr_abs(by, by, MPFR_RNDU);
        mpfr_set(b1, by, MPFR_RNDD);
        mpfr_nextabove(by);
        mpfr_sub(br2, by, b1, MPFR_RNDU);
        
        mpfr_add(br1, br1, br2, MPFR_RNDU);
        mpfr_add(d->r, d->r, br1, MPFR_RNDU);
    }
    
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    
    return true;
}

void mpd_sqr(mpd d, mpd a) { // similar to mpDisk_mul(d,a,a) but optimized
    // allocate 3 coords buffers on the stack
    ulong cp = mpd_prec(d);
    ulong rp = mpd_rad_prec(d);
    defs_mpfr(cp + MPD_EXTRA_PREC, bx, by);
    defs2_mpfr(rp + MPD_EXTRA_PREC, br1, br2, br3);
    
    // compute the real part in b4, rounding errors in br1
    mpd_sqs(bx, br1, a->x, a->y, cp, rp, true);
    
    // compute the imaginary part and rounding error in br2
    mpfr_mul(by, a->x, a->y, MPFR_RNDN); // b5 stores y
    mpd_mul_2p(by, 1);
    mpd_hulp(br2, by, mpd_ulpp(by));
    
    // total rounding errors for the center in br1
    mpfr_hypot(br1, br1, br2, MPFR_RNDU);
    
    mpd_mul_2p(mpd_mod_up(a), 1);
    mpfr_add(br2, a->mu, a->r, MPFR_RNDU); // second use of a->mu is guarateed to be initialized
    mpd_mul_2p(a->mu, -1);
    
    mpfr_mul(d->r, a->r, br2, MPFR_RNDU);
    
    // correct radius with rounding errors
    mpfr_add(d->r, d->r, br1, MPFR_RNDU);
    
    // must be here, otherwise we may rewrite the center of a
    mpfr_set(d->x, bx, MPFR_RNDN);
    mpfr_set(d->y, by, MPFR_RNDN);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
}

void mpd_add_si(mpd d, mpd a, long b) { // d = a + b
    int nex = mpfr_add_si(d->x, a->x, b, MPFR_RNDN);
    int ney = mpfr_set(d->y, a->y, MPFR_RNDN);
    mpfr_set(d->r, a->r, MPFR_RNDU);
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
    
    // most of the time the addition is exact, no need for corrections
    if(! nex && ! ney) {
        return;
    }
    
    defs_mpfr(mpd_rad_prec(d) + MPD_EXTRA_PREC, br1, br2);
    
    if( nex && ! ney ) { // x alone is not exact
        mpd_hulp(br1, d->x, mpd_ulpp(d->x));
        mpfr_add(d->r, d->r, br1, MPFR_RNDU);
        
        return;
    }
    
    if(! nex && ney) { // y alone is not exact
        mpd_hulp(br2, d->y, mpd_ulpp(d->y));
        mpfr_add(d->r, d->r, br2, MPFR_RNDU);
        
        return;
    }
    
     // both are not exact
    mpd_hulp(br1, d->x, mpd_ulpp(d->x));
    mpd_hulp(br2, d->y, mpd_ulpp(d->y));
    mpfr_hypot(br1, br1, br2, MPFR_RNDU);
    
    mpfr_add(d->r, d->r, br1, MPFR_RNDU);
}

void mpd_mull(mpd d, mpd a, long b) { // d = a * b
    defs_mpfr(mpd_rad_prec(d) + MPD_EXTRA_PREC, br1, br2);
    
    long abs_b = b < 0 ? -b : b;
    
    mpfr_mul_si(d->x, a->x, b, MPFR_RNDN);
    mpfr_mul_si(d->y, a->y, b, MPFR_RNDN);
    
    mpfr_mul_si(d->r, a->r, abs_b, MPFR_RNDU);
    mpd_hulp(br1, d->x, mpd_ulpp(d->x));
    mpd_hulp(br2, d->y, mpd_ulpp(d->y));
    mpfr_hypot(br1, br1, br2, MPFR_RNDU);
    
    mpfr_add(d->r, d->r, br1, MPFR_RNDU);
    
    mpfr_set_nan(d->mu);
    mpfr_set_nan(d->md);
}

void mpd_scale(mpd d, mpd a, int twoPow) { // d = a * 2^twoPow
    if(d == a) { // self-substitution, no rounding
        mpd_mul_2p(d->x, twoPow);
        mpd_mul_2p(d->y, twoPow);
        if(! mpfr_nan_p(d->mu))
            mpd_mul_2p(d->mu, twoPow);
        if(! mpfr_nan_p(d->md))
            mpd_mul_2p(d->md, twoPow);
        mpd_mul_2p(d->r, twoPow);

        return;
    }
    
    if(! mpfr_nan_p(a->mu))
        mpfr_mul_2si(d->mu, a->mu, twoPow, MPFR_RNDU);
    else
        mpfr_set_nan(d->mu);
    
    if(! mpfr_nan_p(a->md)) {
        mpfr_mul_2si(d->md, a->md, twoPow, MPFR_RNDD);
    } else {
        mpfr_set_nan(d->md);
    }
    
    mpfr_mul_2si(d->r, a->r, twoPow, MPFR_RNDU);
    
    defs_mpfr(mpd_rad_prec(d) + MPD_EXTRA_PREC, br1, br2);
    
    // target is different from source
    int nex = mpfr_mul_2si(d->x, a->x, twoPow, MPFR_RNDN);
    int ney = mpfr_mul_2si(d->y, a->y, twoPow, MPFR_RNDN);
    
    if(! nex) {
        if(! ney) {
            return;
        }
        
        mpd_hulp(br2, d->y, mpd_ulpp(d->y));
        mpfr_add(d->r, d->r, br2, MPFR_RNDU);
        
        return;
    }
    
    if(! ney) {
        mpd_hulp(br1, d->x, mpd_ulpp(d->x));
        mpfr_add(d->r, d->r, br1, MPFR_RNDU);
        
        return;
    }
    
    mpd_hulp(br1, d->x, mpd_ulpp(d->x));
    mpd_hulp(br2, d->y, mpd_ulpp(d->y));
    mpfr_hypot(br1, br1, br2, MPFR_RNDU);
    
    mpfr_add(d->r, d->r, br1, MPFR_RNDU);
}

bool mpd_contains_0(mpd a, bool guarantee) {
    if(guarantee) { // positive result is guaranteed, open disk
        return mpfr_cmp(mpd_mod_up(a), a->r) < 0;
    }
    
    // negative result is guaranteed, closed disk
    return mpfr_cmp(mpd_mod_down(a), a->r) <= 0;
}

bool mpd_contains_c(mpd d, mpc a, bool guarantee) {
    defs_mpfr(mpd_rad_prec(d) + MPD_EXTRA_PREC, br1, br2);
    
    if(guarantee) { // positive result is guaranteed, open disk
        mpfr_sub(br1, d->x, a->x, MPFR_RNDA);
        mpfr_sub(br2, d->y, a->y, MPFR_RNDA);
        mpfr_hypot(br1, br1, br2, MPFR_RNDU); // upper bound for the distance between "centers"
        
        return mpfr_cmp(d->r, br1) > 0;
    }
    
    // negative result is guaranteed, closed disk
    mpfr_sub(br1, d->x, a->x, MPFR_RNDZ);
    mpfr_sub(br2, d->y, a->y, MPFR_RNDZ);
    mpfr_hypot(br1, br1, br2, MPFR_RNDD); // lower bound for the distance between "centers"
    
    return mpfr_cmp(d->r, br1) >= 0;
}

bool mpd_contains_c80(mpd d, fp80 a, bool guarantee) {
    defs_mpfr(mpd_rad_prec(d) + MPD_EXTRA_PREC, br1, br2);
    
    mpfr_set_ld(br1, a->x, MPFR_RNDN); // should be exact
    mpfr_set_ld(br2, a->y, MPFR_RNDN); // should be exact
    
    if(guarantee) { // positive result is guaranteed, open disk
        mpfr_sub(br1, d->x, br1, MPFR_RNDA);
        mpfr_sub(br2, d->y, br2, MPFR_RNDA);
        mpfr_hypot(br1, br1, br2, MPFR_RNDU); // upper bound for the distance between "centers"
        
        return mpfr_cmp(d->r, br1) > 0;
    }
    
    // negative result is guaranteed, closed disk
    mpfr_sub(br1, d->x, br1, MPFR_RNDZ);
    mpfr_sub(br2, d->y, br2, MPFR_RNDZ);
    mpfr_hypot(br1, br1, br2, MPFR_RNDD); // lower bound for the distance between "centers"
    
    return mpfr_cmp(d->r, br1) >= 0;
}

bool mpd_contains(mpd d, mpd a, bool guarantee) {
    defs_mpfr(mpd_rad_prec(d) + MPD_EXTRA_PREC, br1, br2);
    
    if(guarantee) { // positive result is guaranteed, d open disk, a closed disk
        if(mpfr_cmp(d->r, a->r) <= 0)
            return 0;
        
        mpfr_sub(br1, d->x, a->x, MPFR_RNDA);
        mpfr_sub(br2, d->y, a->y, MPFR_RNDA);
        mpfr_hypot(br1, br1, br2, MPFR_RNDU); // upper bound for the distance between centers
        mpfr_add(br1, br1, a->r, MPFR_RNDU);
        
        return mpfr_cmp(d->r, br1) > 0;
    }
    
    // negative result is guaranteed, d closed disk, a open disk
    if(mpfr_cmp(d->r, a->r) < 0)
        return 0;
    
    mpfr_sub(br1, d->x, a->x, MPFR_RNDZ);
    mpfr_sub(br2, d->y, a->y, MPFR_RNDZ);
    mpfr_hypot(br1, br1, br2, MPFR_RNDD); // lower bound for the distance between centers
    mpfr_add(br1, br1, a->r, MPFR_RNDD);
    
    return mpfr_cmp(d->r, br1) >= 0;
}

bool mpd_intersect(mpd a, mpd b, bool guarantee) {
    defs_mpfr(mpd_rad_prec(a) + MPD_EXTRA_PREC, br1, br2);
    
    if(guarantee) { // positive result is guaranteed, open disks
        mpfr_sub(br1, a->x, b->x, MPFR_RNDA);
        mpfr_sub(br2, a->y, b->y, MPFR_RNDA);
        mpfr_hypot(br1, br1, br2, MPFR_RNDU); // upper bound for the distance between centers
        
        mpfr_add(br2, a->r, b->r, MPFR_RNDD); // lower bound for the sum of radii
        
        return mpfr_cmp(br2, br1) > 0;
    }
    
    // negative result is guaranteed, closed disks
    mpfr_sub(br1, a->x, b->x, MPFR_RNDZ);
    mpfr_sub(br2, a->y, b->y, MPFR_RNDZ);
    mpfr_hypot(br1, br1, br2, MPFR_RNDD); // lower bound for the distance between centers
    
    mpfr_add(br2, a->r, b->r, MPFR_RNDU); // upper bound for the sum of radii
    
    return mpfr_cmp(br2, br1) >= 0;
}

void mpd_min_mod(mpfr_t d, mpd s) {
    if(mpfr_cmp(mpd_mod_down(s), s->r) <= 0) {
        mpfr_set_zero(d, 1);
        
        return;
    }
    
    mpfr_sub(d, mpd_mod_down(s), s->r, MPFR_RNDD);
}

void mpd_max_mod(mpfr_t d, mpd s) {
    mpfr_add(d, mpd_mod_up(s), s->r, MPFR_RNDU);
}
