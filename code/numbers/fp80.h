//
//  fp80.h
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
 \file fp80.h
 \brief A collection of basic functions with hardware accelerated @c fp80 (@c long @c double) complex numbers.
 
  The notation @c fp80 is a shortcut for @code floating point on extended 80-bit precision @endcode.
  
  @see [long double] (https://en.wikipedia.org/wiki/Long_double) and
  [Extended precision] (https://en.wikipedia.org/wiki/Extended_precision)
*/

#ifndef fp80_h
#define fp80_h

#include <stdlib.h>
#include "ntypes.h"

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Structures
// //////////////////////////////////////////////////////////////////////////////////////////

#define FP80_TPM64       (5.4211e-20)
#define FP80_TPM64P5     (3.8333e-20)

/// @struct fp80_struct
/// @brief Complex numbers with 80-bit of precision.
///
/// The \f$ ulp(x) = 2^{exp(x)-63} \f$ or about \f$ 10^{-19} \cdot |x| \f$,
/// where \f$ exp(x) \f$ is the 15-bit exponent in the representation of the number \c x and \c |x| is the absolute value of \c x.
///
/// All apps in this project check that the compiler implements long double as an 80-bit extended precision format.
///
/// @see [long double] (https://en.wikipedia.org/wiki/Long_double) and
/// [Extended precision] (https://en.wikipedia.org/wiki/Extended_precision)
typedef struct {
    long double x;  ///< the real part of the complex number
    long double y;  ///< the imaginary part of the complex number
} fp80_struct;

/// \brief Practical wrapper for @c fp80_struct
///
/// To avoid the constant use @c * and @c & the type  @c fp80 is a pointer.
///
/// Example of use:
/// @code
///    fp80 c;
///    c->x = 0.5L;
/// @endcode
typedef fp80_struct fp80[1];

/// Convenience pointer to fp80_struct
typedef fp80_struct *fp80_ptr;

/// @struct fp64_struct
/// @brief Complex numbers with 64-bits of precision.
///
/// The \f$ ulp(x) = 2^{exp(x)-53} \f$ or about \f$ 10^{-17} \cdot |x| \f$,
/// where \f$ exp(x) \f$ is the 11-bit exponent in the representation of the number \c x and \c |x| is the absolute value of \c x.
typedef struct {
    double x;  ///< the real part of the complex number
    double y;  ///< the imaginary part of the complex number
} fp64_struct;

/// \brief Practical wrapper for @c fp64_struct
///
/// To avoid the constant use @c * and @c & the type  @c fp64 is a pointer.
///
/// Example of use:
/// @code
///    fp64 c;
///    c->x = 0.5;
/// @endcode
typedef fp64_struct fp64[1];

/// Convenience pointer to fp64_struct
typedef fp64_struct *fp64_ptr;

/// Upper bound of the @c ulp of the number @c c
#define fp80_uulp(c) ((((c)->x < 0 ? -(c)->x : (c)->x) + ((c)->y < 0 ? -(c)->y : (c)->y)) * FP80_TPM64)

/// Lower bound of the @c ulp of the number @c c
#define fp80_lulp(c) ((((c)->x < 0 ? -(c)->x : (c)->x) + ((c)->y < 0 ? -(c)->y : (c)->y)) * FP80_TPM64P5)


// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Conversions
// //////////////////////////////////////////////////////////////////////////////////////////

void fp80_setl(fp80 d, fp80 s);
void fp80_set(fp80 d, ldbl x, ldbl y);
void fp80_setd(fp80 d, fp64 s);
void fp64_setl(fp64 d, fp80 s);
void fp64_setd(fp64 d, fp64 s);
void fp64_set(fp80 d, double x, double y);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Additions
// //////////////////////////////////////////////////////////////////////////////////////////

/// @brief Complex additon among @c fp80 numbers.
///
/// Computes \f$ (a1 + i b1) + (a2 + i b2) = (a1 + a2) + i (b1 + b2) \f$ as an @c fp80 number
/// and stores the result in @c v. The numbers @c x, @c y and @c v do not need to be distinct.
/// @param v complex number to store the result in (@c fp80)
/// @param x first complex number (@c fp80) x = a1 + i b1
/// @param y second complex number (@c fp80) y = a2 + i b2
void fp80_add(fp80 v, fp80 x, fp80 y);

void fp64_add(fp64 v, fp64 x, fp64 y);

/// @brief Complex substraction among @c fp80 numbers.
///
/// Computes \f$ (a1 + i b1) - (a2 + i b2) = (a1 - a2) + i (b1 - b2) \f$ as an @c fp80 number
/// and stores the result in @c v. The numbers @c x, @c y and @c v do not need to be distinct.
/// @param v complex number to store the result in (@c fp80)
/// @param x first complex number (@c fp80) x = a1 + i b1
/// @param y second complex number (@c fp80) y = a2 + i b2
void fp80_sub(fp80 v, fp80 x, fp80 y);

void fp64_sub(fp64 v, fp64 x, fp64 y);

void fp80_neg(fp80 v, fp80 x);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Multiplications
// //////////////////////////////////////////////////////////////////////////////////////////

/// @brief Complex multiplication among @c fp80 numbers.
///
/// Computes \f$ (a1 + i b1) * (a2 + i b2) = a1.a2 - b1.b2 + i (a1.b2 + a2.b1) \f$ as an @c fp80 number
/// and stores the result in @c v. The numbers @c x, @c y and @c v do not need to be distinct.
/// @param v complex number to store the result in (@c fp80)
/// @param x first complex number (@c fp80)  x = a1 + i b1
/// @param y second complex number (@c fp80)  y = a2 + i b2
void fp80_mul(fp80 v, fp80 x, fp80 y);

void fp64_mul(fp64 v, fp64 x, fp64 y);

/// @brief Multiplies an @c fp80 complex number by a @c double real number.
///
/// Computes \f$ (a + i b) * y = a*y + i b*y \f$ as an @c fp80 number
/// and stores the result in @c v. The numbers @c x and @c v do not need to be distinct.
/// @param v complex number to store the result in (@c fp80)
/// @param x complex number (@c fp80 ) x = a + i b
/// @param y real number (@c double)
void fp80_muld(fp80 v, fp80 x, double y);

void fp64_muld(fp64 v, fp64 x, double y);

/// @brief Multiplies an @c fp80 complex number by a power of two.
///
/// @param v complex number to store the result to
/// @param x complex number to multiply
/// @param tp the power of two to multiply by
void fp80_scale(fp80 v, fp80 x, int tp);

/// @brief Multiplies an @c fp80 complex number by a \ref ldbl real number.
///
/// Computes \f$ (a + i b) * y = a*y + i b*y \f$ as an @c fp80 number
/// and stores the result in @c v. The numbers @c x and @c v do not need to be distinct.
/// @param v complex number to store the result in (@c fp80)
/// @param x complex number (@c fp80 ) x = a + i b
/// @param y real number (\ref ldbl)
void fp80_mull(fp80 v, fp80 x, ldbl y);

/// @brief Multiplies an @c fp80 complex number by a @c long integer.
///
/// Computes \f$ (a + i b) * y = a*y + i b*y \f$ as an @c fp80 number
/// and stores the result in @c v. The numbers @c x and @c v do not need to be distinct.
/// @param v complex number to store the result in (@c fp80)
/// @param x complex number (@c fp80 ) x = a + i b
/// @param y integer (@c long)
void fp80_muli(fp80 v, fp80 x, long y);

void fp64_muli(fp64 v, fp64 x, long y);

/// @brief Square of an @c fp80 complex number
///
/// Computes \f$ (a + i b)^2 = a^2 - b^2 + i 2*a*b \f$  as an @c fp80 number
/// and stores the result in @c v. The numbers @c x and @c v do not need to be distinct.
///
/// Optimized code compared to @c fp80_mul(v,x,x);
/// @param v complex number to store the result in (@c fp80)
/// @param x complex number (@c fp80) x = a + i b
void fp80_sqr(fp80 v, fp80 x);

void fp64_sqr(fp64 v, fp64 x);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Divisions
// //////////////////////////////////////////////////////////////////////////////////////////

bool fp80_inv(fp80 v, fp80 z);

/// @brief Complex division among @c fp80 numbers.
///
/// Computes \f$ (a1 + i b1) / (a2 + i b2) = (a1 + i b1) * (a2 + i b2) / |a1 + i b1|^2 \f$ as an @c fp80 number
/// and stores the result in @c v. The numbers @c x, @c y and @c v do not need to be distinct.
///
/// If a division by zero is attempted, v is left untouched, no error is generated and @c fp80_div returns 0.
/// Otherwise, it returns 1.
/// @param v complex number to store the result in (@c fp80)
/// @param x first complex number (@c fp80) x = a1 + i b1
/// @param y second complex number (@c fp80) y = a2 + i b2
///
/// @return @ref true if successfull, @ref false otherwise
bool fp80_div(fp80 v, fp80 x, fp80 y);

bool fp64_div(fp64 v, fp64 x, fp64 y);

void fp80_quick_div(fp80 v, fp80 x, fp80 y);

void fp64_quick_div(fp64 v, fp64 x, fp64 y);

/// @brief Divides an @c fp80 complex number by a @c double real number.
///
/// Computes \f$ (a + i b) / y = a/y + i b/y \f$ as an @c fp80 number
/// and stores the result in @c v. The numbers @c x and @c v do not need to be distinct.
///
/// If a division by zero is attempted, v is left untouched, no error is generated and @c fp80_div returns 0.
/// Otherwise, it returns 1.
/// @param v complex number to store the result in (@c fp80)
/// @param x complex number (@c fp80) x = a + i b
/// @param y real number (@c double)
///
/// @return @ref true if successfull, @ref false otherwise
bool fp80_divd(fp80 v, fp80 x, double y);

/// @brief Divides an @c fp80 complex number by a \ref ldbl real number.
///
/// Computes \f$ (a + i b) / y = a/y + i b/y \f$ as an @c fp80 number
/// and stores the result in @c v. The numbers @c x and @c v do not need to be distinct.
///
/// If a division by zero is attempted, v is left untouched, no error is generated and @c fp80_div returns 0.
/// Otherwise, it returns 1.
/// @param v complex number to store the result in (@c fp80)
/// @param x complex number (@c fp80) x = a + i b
/// @param y real number (\ref ldbl)
///
/// @return @ref true if successfull, @ref false otherwise
bool fp80_divl(fp80 v, fp80 x, ldbl y);

/// @brief Divides an @c fp80 complex number by a @c long integer.
///
/// Computes \f$ (a + i b) / y = a/y + i b/y \f$ as an @c fp80 number
/// and stores the result in @c v. The numbers @c x and @c v do not need to be distinct.
///
/// If a division by zero is attempted, v is left untouched, no error is generated and @c fp80_div returns 0.
/// Otherwise, it returns 1.
/// @param v complex number to store the result in (@c fp80)
/// @param x complex number (@c fp80) @c x=a+ib
/// @param y integer (@c long)
///
/// @return @ref true if successfull, @ref false otherwise
bool fp80_divi(fp80 v, fp80 x, long y);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Distance related
// //////////////////////////////////////////////////////////////////////////////////////////

/// @brief Complex modulus.
///
/// Computes \f$ |x+iy| = \sqrt{ x^2 + y^2 } \f$ as a \ref ldbl.
///
/// Indicative only, as there is no rounding control for the square root (and costly).
/// Adapt the code to use @c fp80_mod2 instead, whenever possible.
/// @param p  complex number x+iy (@c fp80)
///
/// @return modulus of @c p (\ref ldbl)
ldbl fp80_mod(fp80 p);

double fp64_mod(fp64 p);

/// @brief Square of complex modulus.
///
/// Computes \f$ |x+iy|^2 = x^2+y^2 \f$ as a \ref ldbl.
/// @param p complex number x+iy (@c fp80)
///
/// @return squared modulus of @c p (\ref ldbl)
ldbl fp80_mod2(fp80 p);

double fp64_mod2(fp64 p);

/// @brief Distance between two complex numbers.
///
/// Computes \f$ |a-b| \f$ as a \ref ldbl.
///
/// Indicative only, as there is no rounding control for the square root (and costly).
/// Adapt the code to use @c fp80_dist2 instead, whenever possible.
/// @param a first complex number (@c fp80)
/// @param b second complex number (@c fp80)
///
/// @return distance between @c a and @c b (\ref ldbl)
ldbl fp80_dist(fp80 a, fp80 b);

double fp64_dist(fp64 a, fp64 b);

/// @brief Square of the distance between two complex numbers.
///
/// Computes \f$ |a-b|^2 \f$ as a \ref ldbl.
/// @param a first complex number (@c fp80)
/// @param b second complex number (@c fp80)
///
/// @return squared distance between @c a and @c b (\ref ldbl)
ldbl fp80_dist2(fp80 a, fp80 b);

double fp64_dist2(fp64 a, fp64 b);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Trigonometric
// //////////////////////////////////////////////////////////////////////////////////////////

/// @brief Computes an  @c fp80 complex number on the unit circle.
///
/// Computes \f$ cos(2i*pi*theta) + i sin(2i*pi*theta) \f$ as an @c fp80 complex number
/// and stores the result in @c v.
/// @param v complex number to store the result in (@c fp80)
/// @param theta angle, normalized by 2*pi ( \ref ldbl)
///
/// @return @ref true if successfull, @ref false otherwise
bool fp80_exp_2Pi_i(fp80 v, ldbl theta);

bool fp64_exp_2Pi_i(fp64 v, double theta);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Other
// //////////////////////////////////////////////////////////////////////////////////////////

/// @brief Square root of an @c fp80 complex number
///
/// Computes \f$ a + ib = \sqrt{z} \f$   with @c a>=0 (and @c b>=0 if @c a==0)as an @c fp80 number
/// and stores the result in @c v. The numbers @c x and @c v do not need to be distinct.
///
/// @param v complex number to store the result in (@c fp80)
/// @param z complex number (@c fp80) x = a + i b
void fp80_sqrt(fp80 v, fp80 z);

bool fp80_is_exact(fp80 c, ldbl err);

void fp80_min(fp80 min, fp80 a, fp80 b);
void fp80_max(fp80 max, fp80 a, fp80 b);
void fp80_min_max(fp80 min, fp80 max, fp80 a);
ldbl fp80_det(fp80 a, fp80 b);

void fp80_print(fp80 z);

bool fp80_is_number(fp80 z);

bool fp64_is_number(fp64 z);

#endif /* fp80_h */
