//
//  mpd.h
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
 \file mpd.h
 \brief A collection of basic functions with multi precision complex disks.
 
 Slower than complex numbers @c mpc_struct, operation with disks produce certified (or proven) results.
 Based on  [mpfr] (https://www.mpfr.org).
*/

#ifndef mpDisk_h
#define mpDisk_h

#include <mpfr.h>

#include "ntypes.h"

#define MPD_EXTRA_PREC 6

typedef struct {
    // center coordinates, modulus upper and lower bounds
    mpfr_t x, y, mu, md;
    
    // radius, initialized to precision precR by mpDisk_init
    mpfr_t r;
} mpd_struct;

// to avoid to constantly use * and &; pointers are passed to methods
// examples of use:
//      mpDisk c;
//      mpDisk_init(c, 120, 120);
//      ...
//      mpDisk_clear(c);
typedef mpd_struct mpd[1];

#include "mpc.h"
#include "fp80d.h"

/// The number of limbs needed to store disks with coords precision @c prec and radius precision @c rad_prec.
#define mpd_limbs(prec, rad_prec)  ((4 * mpfr_limbs(prec)) + mpfr_limbs(rad_prec))

/// @brief Querry for the precision the center of an @c mpd disk.
///
/// This method assumes that the precision on the real and imaginary parts are equal, which is the case in
/// normal use (i.e. when using the public methods of @c mpd.h).
///
/// @param d @c mpd the disk
///
/// @return precision in bits (@c long)
#define mpd_prec(d)     (mpfr_get_prec(d->x))

/// @brief Querry for the precision the radius of an @c mpd disk.
///
/// @param d @c mpd the disk
///
/// @return precision in bits (@c long)
#define mpd_rad_prec(d)     (mpfr_get_prec(d->r))

///  Convenience local definition of a local mpd variable, initialized to (0, 0, 0); no cleanup is needed at the end of the code block
///  @warning to be used only once in the same code block, exclusive of @c defs_mpd()
#define def_mpd(prec, precR, d) \
    ulong dlimbs[mpd_limbs(prec, precR)];\
    mpd d;\
    mpd_iniz(d, prec, precR, climbs)

///  Convenience local definition of local mpd variables, initialized to (0, 0, 0); no cleanup is needed at the end of the code block
///  @warning to be used only once in the same code block, exclusive of @c def_mpd(), with at most 15 arguments
#define defs_mpd(prec, precR, ...)  \
    ulong dlimbs[mpd_limbs(prec, precR) * ARG_COUNT(__VA_ARGS__)];\
    mpd __VA_ARGS__;\
    mpd_initz(prec, precR, dlimbs, __VA_ARGS__, NULL)

///  Convenience local definition of a vector of mpd variables, initialized to (0, 0, 0); no cleanup is needed at the end of the code block
///  @warning to be used only once in the same code block
#define defv_mpd(vect, len, prec, precR)  \
    long dvsize = mpd_limbs(prec, precR); \
    ulong dvlimbs[dvsize * (len)];\
    mpd_struct vect[len]; \
    for(int i = 0; i < len; i++) { \
        mpd_iniz(vect + i, prec, precR, dvlimbs + dvsize * i); \
    }

///  Convenience local definition of a vector of mpd variables, initialized to (0, 0, 0); no cleanup is needed at the end of the code block
///  @warning to be used only once in the same code block
#define defv2_mpd(vect, len, prec, precR)  \
    long dvsize2 = mpd_limbs(prec, precR); \
    ulong dvlimbs2[dvsize2 * (len)];\
    mpd_struct vect[len]; \
    for(int i = 0; i < len; i++) { \
        mpd_iniz(vect + i, prec, precR, dvlimbs2 + dvsize2 * i); \
    }

bool mpd_init(mpd d, long prec, long precR);
bool mpd_set_prec(mpd c, long prec, long precR);

bool mpd_iniz(mpd d, long prec, long precR, ulong *limbs);
int mpd_initz(long prec, long precR, ulong *limbs, mpd d, ...) __MPFR_SENTINEL_ATTR;

bool mpd_clear(mpd d);
int mpd_clears(mpd d, ...) __MPFR_SENTINEL_ATTR;

void mpd_set(mpd d, mpd s);
void mpd_set80(mpd d, fp80d s);
void mpd_setr(mpd d, mpc c, mpfr_t r);
void mpd_setrl(mpd d, mpc c, long double r);
void mpd_set80r(mpd d, fp80 c, long double r);
void mpd_setd(mpd d, double x, double y, double r);
void mpd_setld(mpd d, long double x, long double y, long double r);

void mpd_get80(fp80d d, mpd s);

// ulp, zeroExp is the largest exponent of operands (sub, add), in case x == 0
// if x == 0, its exponent is -(2^63 - 1), thus ulp cannot be computed directly
void mpd_ulp(mpfr_t u, mpfr_t x, long zeroExp);
void mpd_hulp(mpfr_t u, mpfr_t x, long prevExp); // half ulp

void mpd_set_ulp(mpd d, mpfr_t x, mpfr_t y);
void mpd_set_half_ulp(mpd d, mpfr_t x, mpfr_t y);
void mpd_set_exact(mpd d, mpfr_t x, mpfr_t y);
void mpd_set_center(mpd d, mpd a);

void mpd_add(mpd d, mpd a, mpd b);
void mpd_addc(mpd d, mpd a, mpc b);

void mpd_sub(mpd d, mpd a, mpd b);
void mpd_mul(mpd d, mpd a, mpd b);
bool mpd_div(mpd d, mpd a, mpd b);

void mpd_sqr(mpd d, mpd a);

void mpd_add_si(mpd d, mpd a, long b);
void mpd_mull(mpd d, mpd a, long b);
void mpd_scale(mpd d, mpd a, int twoPow);

// in the following functions, if guarantee, the affirmative answer is guaranteed to be correct
// otherwise, the negative answer is guaranteed to be correct
bool mpd_intersect(mpd a, mpd b, bool guarantee);
bool mpd_contains(mpd d, mpd a, bool guarantee);
bool mpd_contains_c(mpd d, mpc a, bool guarantee);
bool mpd_contains_c80(mpd d, fp80 a, bool guarantee);
bool mpd_contains_0(mpd a, bool guarantee);

bool mpd_is_inf(mpd a);
bool mpd_valid(mpd a);

void mpd_min_mod(mpfr_t d, mpd s);
void mpd_max_mod(mpfr_t d, mpd s);

#endif /* mpDisk_h */
