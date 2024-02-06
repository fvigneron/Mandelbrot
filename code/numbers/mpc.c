//
//  mpc.c
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2021.
//
//  Copyright 2019 - 2021 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the GNU Public Licence v3.0
//

#include <mpfr.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "mpc.h"
#include "fp80.h"

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Constructors
// //////////////////////////////////////////////////////////////////////////////////////////

mpc_struct *mpc_new(long prec) {
    if(prec < MP_MIN_PREC) {
        return NULL;
    }
    
    mpc_struct *c = malloc(sizeof(mpc_struct));
    mpc_init(c, prec);
    
    return c;
}

void mpc_free(mpc c) {
    if(c == NULL) {
        return;
    }
    
    mpc_clear(c);
    free(c);
}

bool mpc_init(mpc c, long prec) {
    if(prec < MP_MIN_PREC || c == NULL) {
        return false;
    }
    
    // initialize c: precision is set to prec, value is set to NaN, malloc
    mpfr_inits2(prec, c->x, c->y, NULL);
    
    return true;
}

int mpc_inits(long prec, mpc c, ...) {
    va_list ap;
    va_start(ap, c);
    
    bool ok = mpc_init(c, prec);
    int count = 0;
    
    mpc_struct *a;
    while(ok) {
        count ++;
        
        a = va_arg(ap, mpc_struct *);
        ok = mpc_init(a, prec);
    }
    
    va_end(ap);
    
    return count;
}

inline bool mpfr_iniz(mpfr_ptr x, long prec, ulong *limbs) {
    if(prec < MP_MIN_PREC || limbs == NULL || x == NULL) {
        return false;
    }
    
    x->_mpfr_exp = __MPFR_EXP_ZERO;
    x->_mpfr_sign = 1;
    x->_mpfr_prec = prec;
    x->_mpfr_d = limbs;
    
    return true;
}

int mpfr_initz(long prec, ulong *limbs, mpfr_ptr x, ...) {
    va_list ap;
    va_start(ap, x);
    
    bool ok = mpfr_iniz(x, prec, limbs);
    int count = 0;
    long size = mpfr_limbs(prec);
    
    mpfr_ptr a;
    while(ok) {
        count ++;
        
        a = va_arg(ap, mpfr_ptr);
        ok = mpfr_iniz(a, prec, limbs + count * size);
    }
    
    va_end(ap);
    
    return count;
}

bool mpc_iniz(mpc c, long prec, ulong *limbs) {
    return mpfr_initz(prec, limbs, c->x, c->y, NULL) == 2;
}

int mpc_initz(long prec, ulong *limbs, mpc c, ...) {
    va_list ap;
    va_start(ap, c);
    
    bool ok = mpfr_initz(prec, limbs, c->x, c->y, NULL) == 2;
    int count = 0;
    long size = mpc_limbs(prec);
    
    mpc_struct *a;
    while(ok) {
        count ++;
        
        a = va_arg(ap, mpc_struct *);
        ok = mpfr_initz(prec, limbs + count * size, a->x, a->y, NULL) == 2;
    }
    
    va_end(ap);
    
    return count;
}

bool mpc_set_prec(mpc c, long prec) {
    if(prec < MP_MIN_PREC || c == NULL) {
        return false;
    }
    
    if(prec == mpc_prec(c)) {
        return true;
    }
    
    // initialize c: precision is set to prec, value is set to NaN, malloc
    mpfr_set_prec(c->x, prec);
    mpfr_set_prec(c->y, prec);
    
    return true;
}

inline bool mpc_clear(mpc c) {
    if(c == NULL) {
        return false;
    }
    
    // free the space used by the mantissa (instead of costlier free + reinitialize)
    mpfr_clears(c->x, c->y, NULL);
    
    return true;
}

int mpc_clears(mpc c, ...) {
    va_list ap;
    va_start(ap, c);
    
    bool ok = mpc_clear(c);
    int count = 0;
    
    mpc_struct *a;
    while(ok) {
        count ++;
        
        a = va_arg(ap, mpc_struct *);
        ok = mpc_clear(a);
    }
    
    va_end(ap);
    
    return count;
}

void mpc_init_set(mpc d, mpc s) {
    mpc_init(d, mpc_prec(s));
    mpc_set(d, s);
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Assignment from s (source, various type) into d (destination)
// precision (higher or lower) is imposed by the destination (round to nearest if necessary)
// //////////////////////////////////////////////////////////////////////////////////////////

inline void mpc_set(mpc d, mpc s) {
    mpfr_set(d->x, s->x, MPFR_RNDN);
    mpfr_set(d->y, s->y, MPFR_RNDN);
}

inline void mpc_set0(mpc d) {
    mpfr_set_zero(d->x, 1);
    mpfr_set_zero(d->y, 1);
}

inline void mpc_setmp(mpc d, mpfr_t re, mpfr_t im) {
    mpfr_set(d->x, re, MPFR_RNDN);
    mpfr_set(d->y, im, MPFR_RNDN);
}

inline void mpc_setr(mpc d, mpfr_t re) {
    mpfr_set(d->x, re, MPFR_RNDN);
    mpfr_set_zero(d->y, 1);
}

inline void mpc_set80(mpc d, fp80 s) {
    mpfr_set_ld(d->x, s->x, MPFR_RNDN);
    mpfr_set_ld(d->y, s->y, MPFR_RNDN);
}

inline void mpc_setl(mpc d, long double x, long double y) {
    mpfr_set_ld(d->x, x, MPFR_RNDN);
    mpfr_set_ld(d->y, y, MPFR_RNDN);
}

inline void mpc_setd(mpc d, double x, double y) {
    mpfr_set_d(d->x, x, MPFR_RNDN);
    mpfr_set_d(d->y, y, MPFR_RNDN);
}

inline void mpc_seti(mpc d, long x, long y) {
    mpfr_set_si(d->x, x, MPFR_RNDN);
    mpfr_set_si(d->y, y, MPFR_RNDN);
}

inline bool mpc_set_str(mpc d, char *re, char *im) {
    bool ok = mpfr_set_str(d->x, re, 10, MPFR_RNDN) == 0;
    ok = ok && mpfr_set_str(d->y, im, 10, MPFR_RNDN) == 0;
    
    return ok;
}

inline void mpc_get80(fp80 d, mpc s) {
    d->x = mpfr_get_ld(s->x, MPFR_RNDN);
    d->y = mpfr_get_ld(s->y, MPFR_RNDN);
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Additions
// //////////////////////////////////////////////////////////////////////////////////////////

inline void mpc_add(mpc d, mpc a, mpc b) {
    mpfr_add(d->x, a->x, b->x, MPFR_RNDN);
    mpfr_add(d->y, a->y, b->y, MPFR_RNDN);
}

inline void mpc_sub(mpc d, mpc a, mpc b) {
    mpfr_sub(d->x, a->x, b->x, MPFR_RNDN);
    mpfr_sub(d->y, a->y, b->y, MPFR_RNDN);
}

inline void mpc_addi(mpc d, mpc a, long b) {// add b (long int) to the real part of a
    mpfr_add_si(d->x, a->x, b, MPFR_RNDN);
    mpfr_set(d->y, a->y, MPFR_RNDN);
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Multiplications
// //////////////////////////////////////////////////////////////////////////////////////////

void mpc_mul(mpc d, mpc a, mpc b) {
    defs_mpfr(mpc_prec(d) + MPC_EXTRA_PREC, b1, b2, b3, b4, b5);
    
    if(mpfr_zero_p(a->y)) {
        mpfr_mul(b1, a->x, b->x, MPFR_RNDN);
        mpfr_mul(d->y, a->x, b->y, MPFR_RNDN);
        mpfr_set(d->x, b1, MPFR_RNDN);

        return;
    }

    if(mpfr_zero_p(b->y)) {
        mpfr_mul(b1, a->x, b->x, MPFR_RNDN);
        mpfr_mul(d->y, a->y, b->x, MPFR_RNDN);
        mpfr_set(d->x, b1, MPFR_RNDN);

        return;
    }
    
    mpfr_mul(b1, a->x, b->x, MPFR_RNDN);
    mpfr_mul(b2, a->y, b->y, MPFR_RNDN);

    // for low precision, mul is not much slower than add
    if(mpc_prec(a) <= 128 || mpc_prec(b) <= 128) {
        mpfr_mul(b3, a->x, b->y, MPFR_RNDN);
        mpfr_mul(b4, a->y, b->x, MPFR_RNDN);

        mpfr_sub(d->x, b1, b2, MPFR_RNDN); // Re(a*b)
        mpfr_add(d->y, b3, b4, MPFR_RNDN); // Im(a*b)
    } else {
        mpfr_add(b3, a->x, a->y, MPFR_RNDN);
        mpfr_add(b4, b->x, b->y, MPFR_RNDN);

        mpfr_sub(d->x, b1, b2, MPFR_RNDN); // Re(a*b)
        mpfr_mul(b5, b3, b4, MPFR_RNDN);
        mpfr_sub(b5, b5, b1, MPFR_RNDN);
        mpfr_sub(d->y, b5, b2, MPFR_RNDN);
    }
}

inline void mpc_mulr(mpc d, mpc a, mpfr_t b) {
    mpfr_mul(d->x, a->x, b, MPFR_RNDN);
    mpfr_mul(d->y, a->y, b, MPFR_RNDN);
}

void mpc_sqr(mpc d, mpc a) { // one instruction shorter than mpComplex_mul(d,a,a)
    defs_mpfr(mpc_prec(d) + MPC_EXTRA_PREC, b1, b2);
    
    if(mpfr_zero_p(a->y)) { // much quicker if real
        mpfr_sqr(d->x, a->x, MPFR_RNDN);
        mpfr_set_zero(d->y, 1);
        
        return;
    }
    
    mpfr_sqr(b1, a->x, MPFR_RNDN);
    mpfr_sqr(b2, a->y, MPFR_RNDN);

    mpfr_mul(d->y, a->x, a->y, MPFR_RNDN);
    mpfr_mul_2ui(d->y, d->y, 1, MPFR_RNDN); // dy = 2*ax*ay

    // assignment to d->x last, to avoid corruption of a->x if d == a
    mpfr_sub(d->x, b1, b2, MPFR_RNDN); // dx = ax*ax - ay*ay
}

void mpc_sqrt(mpc d, mpc a) {
    ulong dp = mpfr_get_prec(d->x);
    ulong prec = dp + MPC_EXTRA_PREC;
    
    // check for loss of precision when |y| << |x|
    long de = mpfr_get_exp(a->x) - mpfr_get_exp(a->y);
    if(mpfr_regular_p(a->x) && mpfr_regular_p(a->x) && de > MPC_EXTRA_PREC) {
        prec = dp + (de > dp ? dp : de) + MPC_EXTRA_PREC;
    }
    
    defs_mpfr(prec, b1, b2, b3);
    mpc_mod(b1, a);
    
    mpfr_add(b2, b1, a->x, MPFR_RNDN);
    mpfr_div_2ui(b2, b2, 1, MPFR_RNDN);
    mpfr_sub(b3, b1, a->x, MPFR_RNDN);
    mpfr_div_2ui(b3, b3, 1, MPFR_RNDN);
    
    bool neg = mpfr_cmp_ui(a->y, 0) < 0;
    
    mpfr_sqrt(d->x, b2, MPFR_RNDN);
    mpfr_sqrt(d->y, b3, MPFR_RNDN);
    
    if(neg) {
        mpfr_neg(d->y, d->y, MPFR_RNDN);
    }
}

inline void mpc_muli(mpc d, mpc a, long b) { // scale by b (long int): d = b * a
    mpfr_mul_si(d->x, a->x, b, MPFR_RNDN);
    mpfr_mul_si(d->y, a->y, b, MPFR_RNDN);
}

inline void mpc_muld(mpc d, mpc a, double b) { // scale by b (long int): d = b * a
    mpfr_mul_d(d->x, a->x, b, MPFR_RNDN);
    mpfr_mul_d(d->y, a->y, b, MPFR_RNDN);
}

inline void mpc_scale(mpc d, mpc a, long twoPow) { // scale by power of 2 : d = a * 2^{twoPow}
    mpfr_mul_2si(d->x, a->x, twoPow, MPFR_RNDN);
    mpfr_mul_2si(d->y, a->y, twoPow, MPFR_RNDN);
}

inline void mpc_neg(mpc d, mpc a) {
    mpfr_neg(d->x, a->x, MPFR_RNDN);
    mpfr_neg(d->y, a->y, MPFR_RNDN);
}

inline void mpc_conj(mpc d, mpc a) {
    mpfr_set(d->x, a->x, MPFR_RNDN);
    mpfr_neg(d->y, a->y, MPFR_RNDN);
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Divisions
// //////////////////////////////////////////////////////////////////////////////////////////

bool mpc_div(mpc d, mpc a, mpc b) {
    if(mpfr_nan_p(a->x) || mpfr_nan_p(a->y) || mpfr_nan_p(b->x) || mpfr_nan_p(b->y)) {
        return false;
    }
    
    defs_mpfr(mpc_prec(d) + MPC_EXTRA_PREC, b1, b2, b3, b4, b5);
    
    if(mpfr_zero_p(b->y)) { // much quicker
        if(mpfr_zero_p(b->x)) {
            return false;
        }
        
        mpfr_div(b1, a->x, b->x, MPFR_RNDN);
        mpfr_div(d->y, a->y, b->x, MPFR_RNDN);
        mpfr_set(d->x, b1, MPFR_RNDN);
        
        return true;
    }
    
    mpfr_sqr(b1, b->x, MPFR_RNDN);
    mpfr_sqr(b2, b->y, MPFR_RNDN);
    mpfr_add(b5, b1, b2, MPFR_RNDN); // squared modulus of b
    
    mpfr_mul(b1, a->x, b->x, MPFR_RNDN);
    mpfr_mul(b2, a->y, b->y, MPFR_RNDN);
    mpfr_mul(b3, a->x, b->y, MPFR_RNDN);
    mpfr_mul(b4, a->y, b->x, MPFR_RNDN);
    
    mpfr_add(d->x, b1, b2, MPFR_RNDN); // real part of a * conj(b)
    mpfr_sub(d->y, b4, b3, MPFR_RNDN); // im part of a * conj(b)
    
    mpfr_div(d->x, d->x, b5, MPFR_RNDN); // Re(a/b)
    mpfr_div(d->y, d->y, b5, MPFR_RNDN); // Im(a/b)
    
    return true;
}

bool mpc_divr(mpc d, mpc a, mpfr_t b) {
    if(mpfr_zero_p(b)) {
        return false;
    }
    
    mpfr_div(d->x, a->x, b, MPFR_RNDN);
    mpfr_div(d->y, a->y, b, MPFR_RNDN);
    
    return true;
}

int mpc_divi(mpc d, mpc a, long b) {
    if(b == 0) {
        return 0;
    }
    
    mpfr_div_si(d->x, a->x, b, MPFR_RNDN);
    mpfr_div_si(d->y, a->y, b, MPFR_RNDN);
    
    return 1;
}

bool mpc_inv(mpc d, mpc a) {
    if(mpfr_zero_p(a->y)) { // much quicker
        if(mpfr_zero_p(a->x)) {
            return false;
        }
        
        mpfr_ui_div(d->x, 1, a->x, MPFR_RNDN);
        mpfr_set_zero(d->y, 1);
        
        return true;
    }
    
    defs_mpfr(mpc_prec(d) + MPC_EXTRA_PREC, b1, b2, b3);
    
    mpfr_sqr(b1, a->x, MPFR_RNDN);
    mpfr_sqr(b2, a->y, MPFR_RNDN);
    mpfr_add(b3, b1, b2, MPFR_RNDN); // squared modulus of a
    
    mpfr_set(d->x, a->x, MPFR_RNDN);
    mpfr_neg(d->y, a->y, MPFR_RNDN);
    
    mpfr_div(d->x, d->x, b3, MPFR_RNDN);
    mpfr_div(d->y, d->y, b3, MPFR_RNDN);
    
    return true;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Distance related
// //////////////////////////////////////////////////////////////////////////////////////////

inline void mpc_mod(mpfr_t d, mpc s) { // modulus of s : d = sqrt( sx*sx + sy*sy )
    if(mpfr_zero_p(s->y)) { // much quicker if real
        mpfr_abs(d, s->x, MPFR_RNDU);
        
        return;
    }
    
    mpfr_hypot(d, s->x, s->y, MPFR_RNDU);
}

ldbl mpc_modl(mpc s) {
    fp80 sl;
    mpc_get80(sl, s);
    
    return fp80_mod(sl);
}

void mpc_mod2(mpfr_t d, mpc s) {
    mpfr_sqr(d, s->x, MPFR_RNDU);
    if(mpfr_zero_p(s->y)) {
        return;
    }
    
    def_mpfr(d->_mpfr_prec, buf);
    mpfr_sqr(buf, s->y, MPFR_RNDU);
    
    mpfr_add(d, d, buf, MPFR_RNDU);
}

ldbl mpc_mod2l(mpc s) {
    fp80 sl;
    mpc_get80(sl, s);
    
    return fp80_mod2(sl);
}

void mpc_dist(mpfr_t d, mpc a, mpc b) { // modulus of s : d = sqrt( sx*sx + sy*sy )
    defs_mpfr(mpfr_get_prec(d) + MPC_EXTRA_PREC, b1, b2);
    
    mpfr_sub(b1, a->x, b->x, MPFR_RNDA);
    mpfr_sub(b2, a->y, b->y, MPFR_RNDA);
    
    mpfr_hypot(d, b1, b2, MPFR_RNDU);
}

void mpc_dist2(mpfr_t d, mpc a, mpc b) { // square of the modulus of s : d = sx*sx + sy*sy
    defs_mpfr(mpfr_get_prec(d) + MPC_EXTRA_PREC, b1, b2);
    
    mpfr_sub(b1, a->x, b->x, MPFR_RNDA);
    mpfr_sqr(b1, b1, MPFR_RNDA);
    mpfr_sub(b2, a->y, b->y, MPFR_RNDA);
    mpfr_sqr(b2, b2, MPFR_RNDA);
    
    mpfr_add(d, b1, b2, MPFR_RNDU);
}


bool mpc_close(mpc a, mpc b, mpfr_t error) {
    defs_mpfr(mpc_prec(a) + MPC_EXTRA_PREC, b1, b2);
    
    mpfr_sub(b1, a->x, b->x, MPFR_RNDA);
    mpfr_sqr(b1, b1, MPFR_RNDA);
    mpfr_sub(b2, a->y, b->y, MPFR_RNDA);
    mpfr_sqr(b2, b2, MPFR_RNDA);
    
    mpfr_add(b1, b1, b2, MPFR_RNDU);
    
    return mpfr_cmp(b1, error) < 0;
}


long double mpc_distl(mpc a, mpc b) { // modulus of s : d = sqrt( sx*sx + sy*sy )
    def_mpfr(mpc_prec(a) + MPC_EXTRA_PREC, b1);
    
    mpfr_sub(b1, a->x, b->x, MPFR_RNDA);
    ldbl x = mpfr_get_ld(b1, MPFR_RNDN);
    
    mpfr_sub(b1, a->y, b->y, MPFR_RNDA);
    ldbl y = mpfr_get_ld(b1, MPFR_RNDN);
    
    return hypotl(x, y);
}

ldbl mpc_dist2l(mpc a, mpc b) {
    def_mpfr(mpc_prec(a) + MPC_EXTRA_PREC, b1);
    
    mpfr_sub(b1, a->x, b->x, MPFR_RNDA);
    ldbl x = mpfr_get_ld(b1, MPFR_RNDN);
    
    mpfr_sub(b1, a->y, b->y, MPFR_RNDA);
    ldbl y = mpfr_get_ld(b1, MPFR_RNDN);
    
    return x * x + y * y;
}

inline long mpc_2exp(mpc d) {  // highest power of 2 to radix point
                                        // about [log2(max(|dx|,|dy|))] + 1
    long sx = d->x->_mpfr_exp;
    long sy = d->y->_mpfr_exp;
    
    return sx >= sy ? sx : sy;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Trigonometric
// //////////////////////////////////////////////////////////////////////////////////////////

void mpc_exp_2Pi_i(mpc c, mpfr_t theta) {
    defs_mpfr(mpc_prec(c) + MPC_EXTRA_PREC, b1, b2);
    
    mpfr_frac(b1, theta, MPFR_RNDN);
    if(mpfr_zero_p(b1)) {
        mpfr_set_ui(c->x, 1, MPFR_RNDN);
        mpfr_set_zero(c->y, 1);
        
        return;
    }
    
    if(mpc_prec(c) > 128) {
        mpfr_mul_2si(b2, b1, 3, MPFR_RNDN);
        if(mpfr_integer_p(b2)) {
            ldbl th = mpfr_get_ld(b1, MPFR_RNDN);
            th = th < 0 ? th + 1 : th;
            
            if(th == 0.5) {
                mpfr_set_si(c->x, -1, MPFR_RNDN);
                mpfr_set_zero(c->y, 1);
                
                return;
            }
            
            if(th == 0.25) {
                mpfr_set_zero(c->x, 1);
                mpfr_set_si(c->y, 1, MPFR_RNDN);
                
                return;
            }
            
            if(th == 0.75) {
                mpfr_set_zero(c->x, 1);
                mpfr_set_si(c->y, -1, MPFR_RNDN);
                
                return;
            }
            
            mpfr_sqrt_ui(c->x, 2, MPFR_RNDN);
            mpfr_div_2si(c->x, c->x, 1, MPFR_RNDN);
            mpfr_set(c->y, c->x, MPFR_RNDN);
            
            if(th == 0.375 || th == 0.625) {
                mpfr_neg(c->x, c->x, MPFR_RNDN);
            }
            
            if(th == 0.625 || th == 0.875) {
                mpfr_neg(c->y, c->y, MPFR_RNDN);
            }
            
            return;
        }
    }
    
    mpfr_const_pi(b2, MPFR_RNDN);
    mpfr_mul_2si(b2, b2, 1, MPFR_RNDN);
    
    mpfr_mul(b1, b1, b2, MPFR_RNDN);
    mpfr_sin_cos(c->y, c->x, b1, MPFR_RNDN);
}

bool mpc_eq(mpc a, mpc b) {
    if(a == NULL || b == NULL) {
        return false;
    }
    
    return mpfr_cmp(a->x, b->x) == 0 && mpfr_cmp(a->y, b->y) == 0;
}

bool mpc_is_exact(mpc c, mpfr_t err) {
    if(c == NULL || err == NULL || mpfr_cmp_si(err, 0) < 0) {
        return false;
    }
    
    if(mpfr_inf_p(err)) {
        return true;
    }
    
    bool ez = mpfr_zero_p(err);
    if(! ez && ! mpfr_regular_p(err)) {
        return false;
    }
    
    bool xz = mpfr_zero_p(c->x);
    bool yz = mpfr_zero_p(c->y);
    
    if(xz && yz) {
        return true;
    }
    
    long ee = mpfr_get_exp(err);
    long prec = mpfr_get_prec(c->y);
    long least = ee + prec + 1;
    
    if(xz) {
        if(mpfr_regular_p(c->y)) {
            return mpfr_get_exp(c->y) > least;
        }
        
        return false;
    }
    
    if(yz) {
        if(mpfr_regular_p(c->x)) {
            return mpfr_get_exp(c->x) > least;
        }
        
        return false;
    }
    
    if(! mpfr_regular_p(c->x) || ! mpfr_regular_p(c->y)) {
        return false;
    }
    
    return mpfr_get_exp(c->y) > least || mpfr_get_exp(c->y) > least;
}

bool mpc_ulp(mpfr_t ulp, mpc c) {
    if(c == NULL || ulp == NULL) {
        return false;
    }
    
    defs_mpfr(mpc_prec(c) + MPC_EXTRA_PREC, b1, b2);
    
    mpfr_abs(b1, c->x, MPFR_RNDU);
    mpfr_abs(b2, c->y, MPFR_RNDU);
    mpfr_add(ulp, b1, b2, MPFR_RNDU);
    
    mpfr_div_2ui(ulp, ulp, mpc_prec(c), MPFR_RNDU);
    
    return true;
}

bool mpc_is_number(mpc z) {
    return z != NULL && mpfr_number_p(z->x) && mpfr_number_p(z->y);
}

bool mpc_print(mpc c, int digits) {
    if(c == NULL || digits < 2) {
        return false;
    }
    
    char fmt[20];
    snprintf(fmt, 20, "(%%.%dRg, %%.%dRg)", digits, digits);
    
    mpfr_printf(fmt, c->x, c->y);
    
    return true;
}

bool mpc_snprint(char *str, int len, mpc c, int digits) {
    if(c == NULL || digits < 2 || str == NULL || len < 2 * digits + 12) {
        return false;
    }
    
    char fmt[20];
    snprintf(fmt, 20, "(%%.%dRg, %%.%dRg)", digits, digits);
    
    mpfr_snprintf(str, len, fmt, c->x, c->y);
    
    return true;
}

