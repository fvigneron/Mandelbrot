//
//  mpi.c
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

#include <stdlib.h>

#include "mpi.h"
#include "mpv.h"

mpi mpi_new(long prec) {
    if(prec < MPV_MIN_PREC) {
        return NULL;
    }
    
    mpi x = malloc(sizeof(mpi_struct));
    mpi_init(x, prec);
    
    return x;
}

void mpi_init(mpi x, long prec) {
    if(x == NULL) {
        return;
    }
    
    mpfr_init2(x->a, prec);
    mpfr_init2(x->b, prec);
}

void mpi_free(mpi x) {
    if(x == NULL) {
        return;
    }
    
    mpi_clear(x);
    free(x);
}

void mpi_clear(mpi x) {
    if(x == NULL) {
        return;
    }
    
    mpfr_clear(x->a);
    mpfr_clear(x->b);
}

// MARK: quick operations, no checks

inline void mpi_neg(mpi res, mpi x) {
    if(res != x) {
        mpfr_neg(res->a, x->b, MPFR_RNDD);
        mpfr_neg(res->b, x->a, MPFR_RNDU);
        
        return;
    }
    
    mpfr_t b;
    *b = *x->a;
    *x->a = *x->b;
    *x->b = *b;
    
    mpfr_neg(x->a, x->a, MPFR_RNDD);
    mpfr_neg(x->b, x->b, MPFR_RNDU);
}

void mpi_mul(mpi res, mpi x, mpi y) {
    // allocate a buffer on the stack
    ulong prec = mpi_prec(res);
    mpfr_t b1, b2;
    ulong limbs[mpfr_limbs(prec) * 2];
    
    mpfr_initz(prec, limbs, b1, b2, NULL);
    int sxa = mpfr_sgn(x->a);
    if(sxa >= 0) {
        int sya = mpfr_sgn(y->a);
        if(sya >= 0) {
            mpfr_mul(res->a, x->a, y->a, MPFR_RNDD);
            mpfr_mul(res->b, x->b, y->b, MPFR_RNDU);
        } else {
            int syb = mpfr_sgn(y->b);
            if(syb <= 0) {
                mpfr_mul(b1, x->b, y->a, MPFR_RNDD);
                mpfr_mul(res->b, x->a, y->b, MPFR_RNDU);
                mpfr_set(res->a, b1, MPFR_RNDD);
            } else {
                mpfr_mul(res->a, x->b, y->a, MPFR_RNDD);
                mpfr_mul(res->b, x->b, y->b, MPFR_RNDU);
            }
        }
    } else {
        int sxb = mpfr_sgn(x->b);
        if(sxb <= 0) {
            int sya = mpfr_sgn(y->a);
            if(sya >= 0) {
                mpfr_mul(b1, x->a, y->b, MPFR_RNDD);
                mpfr_mul(res->b, x->b, y->a, MPFR_RNDU);
                mpfr_set(res->a, b1, MPFR_RNDD);
            } else {
                int syb = mpfr_sgn(y->b);
                if(syb <= 0) {
                    mpfr_mul(b1, x->b, y->b, MPFR_RNDD);
                    mpfr_mul(res->b, x->a, y->a, MPFR_RNDU);
                    mpfr_set(res->a, b1, MPFR_RNDD);
                } else {
                    mpfr_mul(b1, x->a, y->b, MPFR_RNDD);
                    mpfr_mul(res->b, x->a, y->a, MPFR_RNDU);
                    mpfr_set(res->a, b1, MPFR_RNDD);
                }
            }
        } else { // here 0 is in the interior of x
            int sya = mpfr_sgn(y->a);
            if(sya >= 0) {
                mpfr_mul(res->a, x->a, y->b, MPFR_RNDD);
                mpfr_mul(res->b, x->b, y->b, MPFR_RNDU);
            } else {
                int syb = mpfr_sgn(y->b);
                if(syb <= 0) {
                    mpfr_mul(b1, x->b, y->a, MPFR_RNDD);
                    mpfr_mul(res->b, x->a, y->a, MPFR_RNDU);
                    mpfr_set(res->a, b1, MPFR_RNDD);
                } else { // here 0 is also in the interior of y
                    mpfr_mul(b1, x->b, y->a, MPFR_RNDD);
                    mpfr_mul(b2, x->a, y->b, MPFR_RNDD);
                    if(mpfr_cmp(b2, b1) < 0) {   // b1 = min(b1, b2)
                        mpfr_set(b1, b2, MPFR_RNDD);
                    }
                    
                    mpfr_mul(b2, x->a, y->a, MPFR_RNDU);
                    mpfr_set(res->a, b1, MPFR_RNDD); // x->a and x->b are not used anymore
                    
                    mpfr_mul(b1, x->b, y->b, MPFR_RNDU);
                    if(mpfr_cmp(b1, b2) > 0) {
                        mpfr_set(res->b, b1, MPFR_RNDU);
                    } else {
                        mpfr_set(res->b, b2, MPFR_RNDU);
                    }
                }
            }
        }
    }
}

void mpi_sqr(mpi res, mpi x) {
    // allocate a buffer on the stack
    def_mpfr(mpi_prec(res), b1);
    
    int sxa = mpfr_sgn(x->a);
    if(sxa >= 0) {
        mpfr_sqr(res->a, x->a, MPFR_RNDD);
        mpfr_sqr(res->b, x->b, MPFR_RNDU);
    } else {
        int sxb = mpfr_sgn(x->b);
        if(sxb <= 0) {
            mpfr_sqr(b1, x->b, MPFR_RNDD);
            mpfr_sqr(res->b, x->a, MPFR_RNDU);
            mpfr_set(res->a, b1, MPFR_RNDD);
        } else { // 0 is in (a, b), the interior of the interval x
            mpfr_neg(b1, x->a, MPFR_RNDU);
            if(mpfr_cmp(b1, x->b) > 0) {
                mpfr_sqr(res->b, b1, MPFR_RNDU);
            } else {
                mpfr_sqr(res->b, x->b, MPFR_RNDU);
            }
            
            mpfr_set_zero(res->a, 1);
        }
    }
}

void mpi_disk(mpfr_t center, mpfr_t radius, mpi x) {
    def_mpfr(mpi_prec(x), b1);
    
    mpfr_add(center, x->a, x->b, MPFR_RNDN);
    mpfr_mul_2si(center, center, -1, MPFR_RNDN);
    
    mpfr_sub(radius, x->b, center, MPFR_RNDU);
    mpfr_sub(b1, center, x->a, MPFR_RNDU);
    if(mpfr_cmp(b1, radius) > 0) {
        mpfr_set(radius, b1, MPFR_RNDU);
    }
}
