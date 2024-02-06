//
//  u128c.c
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

#include "u128c.h"

// MARK: u128 comparisons

inline bool u128_sless(u128 a, u128 b) {
    return a->x < b->x || (a->x == b->x && a->y < b->y);
}

inline bool u128_smore(u128 a, u128 b) {
    return a->x > b->x || (a->x == b->x && a->y > b->y);
}

inline bool u128_leq(u128 a, u128 b) {
    return a->x < b->x || (a->x == b->x && a->y <= b->y);
}

inline bool nset_geq(u128 a, u128 b) {
    return a->x > b->x || (a->x == b->x && a->y <= b->y);
}

inline bool u128_eq(u128 a, u128 b, unsigned long eps) {
    return ! (a->x > b->x + eps || b->x > a->x + eps ||
              a->y > b->y + eps || b->y > a->y + eps);
}

// MARK: u128 conversions

bool u128_get(mpc c, u128 p) {
    return u128_getr(c->x, c->y, p);
}

bool u128_set(u128 p, mpc c) {
    return u128_setr(p, c->x, c->y);
}

bool u128_getr(mpfr_t x, mpfr_t y, u128 p) {
    mpfr_set_ui(x, p->xh, MPFR_RNDN);
    mpfr_mul_2ui(x, x, 64, MPFR_RNDN);
    mpfr_add_ui(x, x, p->xl, MPFR_RNDN);
    mpfr_mul_2si(x, x, -126, MPFR_RNDN);
    mpfr_sub_ui(x, x, 2, MPFR_RNDN);
    
    mpfr_set_ui(y, p->yh, MPFR_RNDN);
    mpfr_mul_2ui(y, y, 64, MPFR_RNDN);
    mpfr_add_ui(y, y, p->yl, MPFR_RNDN);
    mpfr_mul_2si(y, y, -126, MPFR_RNDN);
    
    return true;
}

bool u128_setr(u128 p, mpfr_t x, mpfr_t y) {
    if(mpfr_cmp_si(x, -2) < 0 || mpfr_cmp_ui(x, 2) >= 0 ||
       mpfr_cmp_si(y, -4) <= 0 || mpfr_cmp_ui(y, 4) >= 0) {
        return false;
    }
    
    // allocate a buffer on the stack
    def_mpfr(mpfr_get_prec(x) + MPC_EXTRA_PREC, buf);
    
    mpfr_add_ui(buf, x, 2, MPFR_RNDN);
    mpfr_mul_2ui(buf, buf, 62, MPFR_RNDN);
    p->xh = mpfr_get_ui(buf, MPFR_RNDZ);
    mpfr_sub_ui(buf, buf, p->xh, MPFR_RNDN);
    mpfr_mul_2ui(buf, buf, 64, MPFR_RNDN);
    p->xl = mpfr_get_ui(buf, MPFR_RNDN);
    
    mpfr_mul_2ui(buf, y, 62, MPFR_RNDN);
    mpfr_setsign(buf, buf, 0, MPFR_RNDN);
    
    p->yh = mpfr_get_ui(buf, MPFR_RNDZ);
    mpfr_sub_ui(buf, buf, p->yh, MPFR_RNDN);
    mpfr_mul_2ui(buf, buf, 64, MPFR_RNDN);
    p->yl = mpfr_get_ui(buf, MPFR_RNDN);
    
    return true;
}

bool u128_setl(u128 p, fp80 z) {
    if(z->x < -2 || z->y < 0 || z->x >= 2 || z->y >= 4) {
        return false;
    }
    
    ulong f = 1L << 62;
    p->xh = (z->x + 2) * f;
    p->xl = 0;
    p->yh = z->y * f;
    p->yl = 0;
    
    return true;
}

bool u128_getl(fp80 z, u128 p) {
    z->x = p->xh;
    z->x = ldexpl(z->x, -62);
    z->x -= 2;
    
    z->y = p->yh;
    z->y = ldexpl(z->y, -62);
    
    return true;
}
