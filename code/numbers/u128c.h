//
//  u128c.h
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
 \file u128c.h
 \brief Comparison and conversion functions for complex numbers in the square \f$ [-2,2] \times [0,4] \f$
 stored as couples of \ref uint128.
 
 Their reperesentation is quite specific for points in the upper half plane of the Mandelbrot set. They are used mainly in the
 class \ref nset to achieve good precision, high speed and low memory and file footprint.
*/

#ifndef u128c_h
#define u128c_h

#include <stdio.h>
#include <limits.h>
#include <mpfr.h>

#include "fp80.h"
#include "mpc.h"

// MARK: u128 definitions for points representation

/// Convenience name of a 128-bit unsigned int.
#define uint128 __uint128_t

/// Max value of a uint128
#define U128_MAX ((uint128) ((uint128) -1))

/// ldexpl(1, -126)
#define TWO_M126 (1.175494350822287507968737e-38L)

/// ldexpl(1, 126)
#define TWO_126  (8.507059173023461586584365e+37L)

/// Conversion of a uint128 to a long double.
#define uint128_to_ldbl(a) (a * TWO_M126)

/// Conversion of a long double to a uint128.
#define ldbl_to_uint128(x) ((uint128)(x * TWO_126))

/// Safe subtraction of uint128.
#define uint128_sub(a, b) (a > b ? a - b : 0)

/// Distance of uint128.
#define uint128_dist(a, b) (a > b ? a - b : b - a)

/// Safe addition of uint128.
#define uint128_add(a, b) (a > U128_MAX - b ? U128_MAX : a + b)

/// @union u128c_struct
/// @brief This struct stores positive integer coordinates for complex values.
///
/// The values are divided by \f$ 2^{126} \f$ and then \f$ 2 \f$ is subtracted to obtain the complex number.
typedef union {
    struct {
        uint128 x;         ///< real part
        uint128 y;         ///< imaginary part
    };
    
    uint128 c[2];          ///< view of the coordiantes as an array
    unsigned long ul[4];   ///< the four components of the coordinates, as an array
    
    struct {
        unsigned long xl;  ///< low part of @c x
        unsigned long xh;  ///< high part of @c x
        unsigned long yl;  ///< low part of @c y
        unsigned long yh;  ///< high part of @c y
    };
} u128c_struct;

/// Convenience type for a pointer to @c u128c_struct.
typedef u128c_struct *u128_ptr;

/// Convenience type for a pointer to @c u128c_struct for easy allocation as a local variable.
typedef u128c_struct u128[1];


// MARK: u128 comparisons

/// @brief Checks if @c a<b.
///
/// @param a the first point
/// @param b the second point
///
/// @return @ref true if @c a<b, @ref false otherwise
bool u128_sless(u128 a, u128 b);

/// @brief Checks if @c a>b.
///
/// @param a the first point
/// @param b the second point
///
/// @return @ref true if @c a>b, @ref false otherwise
bool u128_smore(u128 a, u128 b);

/// @brief Checks if @c a<=b.
///
/// @param a the first point
/// @param b the second point
///
/// @return @ref true if @c a<=b, @ref false otherwise
bool u128_leq(u128 a, u128 b);

/// @brief Checks if @c a>=b.
///
/// @param a the first point
/// @param b the second point
///
/// @return @ref true if @c a>=b, @ref false otherwise
bool u128_geq(u128 a, u128 b);

/// @brief Checks if @c a==b up to error @c err on each coordinate.
///
/// @param a the first point
/// @param b the second point
/// @param err the max error
///
/// @return @ref true if @c a==b up to error @c eps on each coordinate, @ref false otherwise
bool u128_eq(u128 a, u128 b, unsigned long err);

// MARK: u128 conversions

/// @brief Converts a fixed point complex @c u128 value to a standard complex number @c mpc.
///
/// Sets @c c to be equal to @c p.
///
/// @param c the result
/// @param p the operand
///
/// @return @ref true if successfull, @ref false otherwise
bool u128_get(mpc c, u128 p);

/// @brief Converts a standard complex number @c mpc to a fixed point complex @c u128 value.
///
/// @param p the result
/// @param c the operand
///
/// @return @ref true if successfull, @ref false otherwise
bool u128_set(u128 p, mpc c);

/// @brief Converts a fixed point complex @c u128 value to its real @c mpfr coordinates.
///
/// @param x the @c x coordiante of the result
/// @param y x the @c y coordiante of the result
/// @param p the operand
///
/// @return @ref true if successfull, @ref false otherwise
bool u128_getr(mpfr_t x, mpfr_t y, u128 p);

/// @brief Converts two real numbers @c mpfr to a fixed point complex @c u128 value.
///
/// @param p the result
/// @param x the @c x coordiante of the operand
/// @param y the @c y coordiante of the operand
///
/// @return @ref true if successfull, @ref false otherwise
bool u128_setr(u128 p, mpfr_t x, mpfr_t y);

/// @brief Converts a low precision complex number @c fp80 to a fixed point complex @c u128 value.
///
/// @param p the result
/// @param z the operand
///
/// @return @ref true if successfull, @ref false otherwise
bool u128_setl(u128 p, fp80 z);

/// @brief Converts a fixed point complex @c u128 value to a low precision complex number @c fp80.
///
/// @param z the result
/// @param p the operand
///
/// @return @ref true if successfull, @ref false otherwise
bool u128_getl(fp80 z, u128 p);

#endif /* u128c_h */
