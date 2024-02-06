//
//  mpc.h
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
 \file mpc.h
 \brief A collection of basic functions with arbitrary precision [mpfr] (https://www.mpfr.org) complex numbers.
 
 The emphasis is on execution speed.
 
 \warning Do not use for large lists, prefer @c mpv_struct (defined in @c mpv.h) instead.
 
 \warning Buffers for computation are created on the stack, do not use very high precision, say > 10M bits.
 
 Use @c mpd_struct (defined in @c mpd.h) to perform operations with certified results.
*/

#include <mpfr.h>

#include "fp80.h"

#ifndef mpComplex_h
#define mpComplex_h

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Structures
// //////////////////////////////////////////////////////////////////////////////////////////

/// The min precision of complex numbers
#define MP_MIN_PREC          32

/// Extra bits of precision for internal buffers (especially for quicker multiplications).
#define MPC_EXTRA_PREC 6

/// The number of limbs needed to store real numbers with precision @c prec.
#define mpfr_limbs(prec)  ((((long) prec - 1) / 64) + 1)

/// The number of limbs needed to store complex numbers with precision @c prec.
#define mpc_limbs(prec)  (2 * mpfr_limbs(prec))

/// @struct mpc_struct
/// @brief Complex numbers with arbitrary precision.
///
/// Coordinates are initialized to precision @c prec by @c mpc_init.
///
/// The structure carries 5 buffers, initialized to the precision \f$ prec + MP_COMPLEX_EXTRA_PREC \f$ by @c mpc_init.
///
/// @see [Extended precision] (https://en.wikipedia.org/wiki/Extended_precision)
/// and [mpfr] (https://www.mpfr.org)
typedef struct {
    mpfr_t x, y;
} mpc_struct;

/// Practical wrapper for @c mpc_struct
///
/// To avoid the constant use @c * and @c & the type  @c mpc is a pointer.
///
/// Example of use:
/// @code
/// mpc c;
/// mpc_init(c, 120);
/// mpc_clear(c);
/// @endcode
typedef mpc_struct mpc[1];

/// Convenience type for a pointer to a @c mpc_struct.
typedef mpc_struct *mpc_ptr;

/// @brief Querry for the precision of an @c mpc number.
///
/// This method assumes that the precision on the real and imaginary parts are equal, which is the case in
/// normal use (i.e. when using the public methods of @c mpc.h).
///
/// @param d @c mpc complex number
///
/// @return precision in bits (@c long)
#define mpc_prec(d)     (mpfr_get_prec(d->x))

#define NTH_ARGUMENT(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, ...) (a17)
#define ARG_COUNT(...) NTH_ARGUMENT(, ## __VA_ARGS__, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

///  Convenience local definition of a local mpc variable, initialized to (0, 0); no cleanup is needed at the end of the code block
///  @warning to be used only once in the same code block, exclusive of @c defs_mpc()
#define def_mpc(prec, x)  \
    ulong climbs[mpc_limbs(prec)];\
    mpc x;\
    mpc_iniz(x, prec, climbs)

///  Convenience local definition of local mpc variables, initialized to (0, 0); no cleanup is needed at the end of the code block
///  @warning to be used only once in the same code block, exclusive of @c def_mpc(), with at most 15 arguments
#define defs_mpc(prec, ...)  \
    ulong climbz[mpc_limbs(prec) * ARG_COUNT(__VA_ARGS__)];\
    mpc __VA_ARGS__;\
    mpc_initz(prec, climbz, __VA_ARGS__, NULL)

///  Convenience local definition of local mpc variables, initialized to (0, 0); no cleanup is needed at the end of the code block
///  @warning to be used only once in the same code block, exclusive of @c def_mpc(), with at most 15 arguments
#define defs2_mpc(prec, ...)  \
    ulong climbs2[mpc_limbs(prec) * ARG_COUNT(__VA_ARGS__)];\
    mpc __VA_ARGS__;\
    mpc_initz(prec, climbs2, __VA_ARGS__, NULL)

///  Convenience local definition of a local mpfr variable, initialized to 0; no cleanup is needed at the end of the code block
///  @warning to be used only once in the same code block, exclusive of @c defs_mpfr()
#define def_mpfr(prec, x)  \
    ulong flimbs[mpfr_limbs(prec)];\
    mpfr_t x;\
    mpfr_iniz(x, prec, flimbs)

///  Convenience local definition of local mpfr variables, initialized to 0; no cleanup is needed at the end of the code block
///  @warning to be used only once in the same code block, exclusive of @c def_mpfr(), with at most 15 arguments
#define defs_mpfr(prec, ...)  \
    ulong flimbz[mpfr_limbs(prec) * ARG_COUNT(__VA_ARGS__)];\
    mpfr_t __VA_ARGS__;\
    mpfr_initz(prec, flimbz, __VA_ARGS__, NULL)

///  Convenience local definition of local mpfr variables, initialized to 0; no cleanup is needed at the end of the code block
///  @warning to be used only once in the same code block, with at most 15 arguments
#define defs2_mpfr(prec, ...)  \
    ulong flimbs2[mpfr_limbs(prec) * ARG_COUNT(__VA_ARGS__)];\
    mpfr_t __VA_ARGS__;\
    mpfr_initz(prec, flimbs2, __VA_ARGS__, NULL)

#define mpc_is_zero(a) (mpfr_zero_p(a->x) && mpfr_zero_p(a->y))

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Constructors
// //////////////////////////////////////////////////////////////////////////////////////////

mpc_struct *mpc_new(long prec);
void mpc_free(mpc c);

/// @brief Initialize an @c mpc complex number.
///
/// Malloc and initialize the mpfr numbers and buffers within @c mpcs. The values are set to NaN.
///
/// The precision @c prec is expressed in bits and should be at least 32. It is the precision of the real and imaginary parts of @c c.
/// The buffers are initialized with an additional precision @c MP_COMPLEX_EXTRA_PREC.
///
/// Normally, an @c mpc number should only be initialized once. Use @c mpc_setPrec to change the precision at a later time.
/// @param c  @c mpc complex number
/// @param prec precision of the significand, in bits
/// @return success status
bool mpc_init(mpc c, long prec);

int mpc_inits(long prec, mpc c, ...) __MPFR_SENTINEL_ATTR;

bool mpfr_iniz(mpfr_ptr x, long prec, ulong *limbs);
int mpfr_initz(long prec, ulong *limbs, mpfr_ptr x, ...) __MPFR_SENTINEL_ATTR;

bool mpc_iniz(mpc c, long prec, ulong *limbs);
int mpc_initz(long prec, ulong *limbs, mpc c, ...) __MPFR_SENTINEL_ATTR;

/// @brief Change the precision and reset an @c mpc complex number.
///
/// One should not expect the previous values to survive. They will usually be set to NaN after the call, except if the new precision
/// is invalid (smaller then 32 or identical to the current one) in which case @c c is not changed.
///
/// @param c  @c mpc complex number
/// @param prec new precision of the significand, in bits
///
/// @return success status
bool mpc_set_prec(mpc c, long prec);

/// @brief Clear an @c mpc complex number.
///
/// Free the space used by the mantissa. The number @c mpc can then be re-intialized with @c mpc_init.
///
/// This method is not as costly as freeing the structure and initializing a new number.
/// In most cases, using @c mpc_setPrec should be sufficient.
/// Before calling @c mpc_initAndSet, use @c mpc_clear on the target if it is not a fresh number.
/// @param c @c mpc complex number
///
/// @return success status
bool mpc_clear(mpc c);

int mpc_clears(mpc c, ...) __MPFR_SENTINEL_ATTR;

/// Copy an @c mpc complex number into a fresh one.
///
/// The target @c d should not have been initialized yet.
///
/// If the target is not fresh, use @c mpc_clear if necessary, or consider using @c mpc_set instead.
///
/// @param d destination complex number (uninitialized @c mpc)
/// @param s source complex number (@c mpc)
void mpc_init_set(mpc d, mpc s);

/// Copy an @c mpc complex number into an existing one.
///
/// The target @c d should already have been initialized.
///
/// If the target is fresh, consider using @c mpc_initAndSet instead.
///
/// @param d destination complex number (@c mpc)
/// @param s source complex number (@c mpc)
void mpc_set(mpc d, mpc s);

/// Set an @c mpc complex number to zero.
///
/// The real and imaginary parts of @c d will be set to (plus) zero.
/// @param d target complex number (@c mpc)
void mpc_set0(mpc d);

/// Set an @c mpc complex number from multi-precision components.
///
/// Regardless of the precision of @c re and @c im, the values will be rounded and assigned to
/// the components of @c d  with the precision defined at the initialization time.
///
/// @param d target complex number (@c mpc)
/// @param re real part, given in multi-precision (@c mpfr_t)
/// @param im imaginary part, given in multi-precision (@c mpfr_t)
void mpc_setmp(mpc d, mpfr_t re, mpfr_t im);

/// Set an @c mpc complex number from a multi-precision real number.
///
/// @param d target complex number (@c mpc)
/// @param re real part, given in multi-precision (@c mpfr_t)
void mpc_setr(mpc d, mpfr_t re);

/// Set an @c mpc complex number from an @c fp80 complex number.
///
/// The components will be rounded and assigned to the components of @c d  with the precision defined at the initialization time.
///
/// @param d target number (@c mpc)
/// @param s source number (@c fp80)
void mpc_set80(mpc d, fp80 s);

/// Set an @c mpc complex number from @c long_double components.
///
/// The components will be rounded and assigned to the components of @c d  with the precision defined at the initialization time.
///
/// @param d target complex number (@c mpc)
/// @param x real part (@c long_double)
/// @param y imaginary part (@c long_double)
void mpc_setl(mpc d, long double x, long double y);

/// Set an @c mpc complex number from @c double components.
///
/// The components will be rounded and assigned to the components of @c d  with the precision defined at the initialization time.
///
/// @param d target complex number (@c mpc)
/// @param x real part (@c double)
/// @param y imaginary part (@c double)
void mpc_setd(mpc d, double x, double y);

/// Set an @c mpc complex number from integer @c long components.
///
/// The components will be rounded and assigned to the components of @c d  with the precision defined at the initialization time.
///
/// @param d target complex number (@c mpc)
/// @param x real part, integer (@c long)
/// @param y imaginary part, integer (@c long)
void mpc_seti(mpc d, long x, long y);

/// Set an @c mpc complex number from string components.
///
/// The strings of the components are converted into numbers, rounded and assigned to the components of @c d  with the precision defined at the initialization time.
///
/// @param d target complex number (@c mpc)
/// @param re real part (string)
/// @param im imaginary part (string)
bool mpc_set_str(mpc d, char *re, char *im);

/// Convert an @c mpc complex number into an @c fp80 complex number.
///
/// The components will be rounded and assigned to the components of @c d .
///
/// @param d target number (@c fp80)
/// @param s source number (@c mpcs)
void mpc_get80(fp80 d, mpc s);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Additions
// //////////////////////////////////////////////////////////////////////////////////////////

/// Complex additon among @c mpc numbers.
///
/// Computes \f$ (a1 + i b1) + (a2 + i b2) = (a1 + a2) + i (b1 + b2) \f$ as an @c mpc number
/// and stores the result in @c d. The numbers @c a, @c b and @c ds do not need to be distinct.
/// @param d complex number to store the result in (@c fp80)
/// @param a first complex number (@c fp80) x = a1 + i b1
/// @param b second complex number (@c fp80) y = a2 + i b2
void mpc_add(mpc d, mpc a, mpc b);

/// Complex substraction among @c mpc numbers.
///
/// Computes \f$ (a1 + i b1) - (a2 + i b2) = (a1 - a2) + i (b1 - b2) \f$ as an @c mpc number
/// and stores the result in @c d. The numbers @c a, @c b and @c ds do not need to be distinct.
/// @param d complex number to store the result in (@c fp80)
/// @param a first complex number (@c fp80) x = a1 + i b1
/// @param b second complex number (@c fp80) y = a2 + i b2
void mpc_sub(mpc d, mpc a, mpc b);

/// Add an integer to an @c mpc complex number.
///
/// Computes \f$ (a1 + i b1) + b = (a1 + b) + i b1 \f$ as an @c mpc number
/// and stores the result in @c d. The numbers @c a and @c d do not need to be distinct.
/// @param d complex number to store the result in (@c fp80)
/// @param a first complex number (@c fp80) x = a1 + i b1
/// @param b real integer (@c long)
void mpc_addi(mpc d, mpc a, long b);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Multiplications
// //////////////////////////////////////////////////////////////////////////////////////////

// only the internal buffers of d are used for computations
void mpc_mul(mpc d, mpc a, mpc b);

void mpc_mulr(mpc d, mpc a, mpfr_t b);

void mpc_sqr(mpc d, mpc a);

void mpc_sqrt(mpc d, mpc a);

void mpc_muli(mpc d, mpc a, long b);

void mpc_muld(mpc d, mpc a, double b);

void mpc_scale(mpc d, mpc a, long twoPow);

void mpc_neg(mpc d, mpc a);

void mpc_conj(mpc d, mpc a);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Divisions
// //////////////////////////////////////////////////////////////////////////////////////////

bool mpc_div(mpc d, mpc a, mpc b);

bool mpc_divr(mpc d, mpc a, mpfr_t b);

int mpc_divi(mpc d, mpc a, long b);

bool mpc_inv(mpc d, mpc a);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Distance related - all upper bounds
// //////////////////////////////////////////////////////////////////////////////////////////

void mpc_mod(mpfr_t d, mpc s);

void mpc_mod2(mpfr_t d, mpc s);

ldbl mpc_modl(mpc s);

ldbl mpc_mod2l(mpc s);

void mpc_dist(mpfr_t d, mpc a, mpc b);

void mpc_dist2(mpfr_t d, mpc a, mpc b);

ldbl mpc_distl(mpc a, mpc b);

ldbl mpc_dist2l(mpc a, mpc b);

long mpc_2exp(mpc d);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Trigonometric
// //////////////////////////////////////////////////////////////////////////////////////////

void mpc_clear_all_buffers(void);

void mpc_exp_2Pi_i(mpc expiTheta, mpfr_t theta);

bool mpc_eq(mpc a, mpc b);
bool mpc_close(mpc a, mpc b, mpfr_t error);

bool mpc_is_exact(mpc c, mpfr_t err);

bool mpc_ulp(mpfr_t ulp, mpc c);

bool mpc_is_number(mpc z);

bool mpc_print(mpc c, int digits);

bool mpc_snprint(char *str, int len, mpc c, int digits);

#endif /* mpComplex_h */
