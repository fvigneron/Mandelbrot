//
//  mpi.h
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

/**
 \file mpi.h
 \brief Definition of a real interval and basic functions for interval arithmetic computations. Based on  [mpfr] (https://www.mpfr.org).
*/

#ifndef mpi_h
#define mpi_h

#include <mpfr.h>

#include "ntypes.h"

typedef struct {
    mpfr_t a; ///< lowest value in the interval
    mpfr_t b; ///< highest value in the interval
} mpi_struct;

/// Pointer to @c mpi_struct with easy allocation as local variable.
typedef mpi_struct mpi_t[1];

/// Pointer to @c mpi_struct.
typedef mpi_struct *mpi;

mpi mpi_new(long prec);
void mpi_init(mpi x, long prec);
void mpi_free(mpi x);
void mpi_clear(mpi x);

#define mpi_set_zero(x)     mpfr_set_zero(x->a, 1);\
                            mpfr_set_zero(x->b, 1)

#define mpi_prec(x)         mpfr_get_prec(x->a)

#define mpi_set(res, x)     mpfr_set(res->a, x->a, MPFR_RNDD); \
                            mpfr_set(res->b, x->b, MPFR_RNDU)

#define mpi_setr(res, a, b) mpfr_set(res->a, a, MPFR_RNDD);\
                            mpfr_set(res->b, b, MPFR_RNDU)

#define mpi_set_si(res, v)  mpfr_set_si(res->a, v, MPFR_RNDD);\
                            mpfr_set_si(res->b, v, MPFR_RNDU)

#define mpi_add(res, x, y)  mpfr_add(res->a, x->a, y->a, MPFR_RNDD);\
                            mpfr_add(res->b, x->b, y->b, MPFR_RNDU)

#define mpi_sub(res, x, y)  mpfr_sub(res->a, x->a, y->b, MPFR_RNDD);\
                            mpfr_sub(res->b, x->b, y->a, MPFR_RNDU)

#define mpi_mul_2si(x, pow) mpfr_mul_2si(x->a, x->a, pow, MPFR_RNDD);\
                            mpfr_mul_2si(x->b, x->b, pow, MPFR_RNDU)

#define mpi_set_ld(res, v)  mpfr_set_ld(res->a, v, MPFR_RNDD);\
                            mpfr_set_ld(res->b, v, MPFR_RNDU)

void mpi_mul(mpi res, mpi x, mpi y);
void mpi_sqr(mpi res, mpi x);
void mpi_neg(mpi res, mpi x);

void mpi_disk(mpfr_t center, mpfr_t radius, mpi x);

#endif /* mpi_h */
