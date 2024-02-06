//
//  fp80d.h
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
 \file fp80d.h
 \brief A collection of basic functions with hardware accelerated @c fp80 (@c long @c double) complex
 disks.
 
 While they cannot produce certified (or proven) results as it is the case of @c mpd_struct, they could be
 used for optimizations in low precision settings.
*/

#ifndef fp80Disk_h
#define fp80Disk_h

#include "fp80.h"

/// Complex disks with 80bit of precision : ulp = 2^{-63} or about 10^{-19}.
/// @c x,y : coordinates of center
/// @c r : radius (when not specified, it defaults to barely enclosing the origin)
/// @c m : the modulus of x + i y is updated automatically by the methods below.
typedef struct {
    ldbl x, y, r, m;
} fp80d_struct;

/// To avoid to constantly use @c * and @c & the type  @c fp80Disk is a pointer.
/// Example of use:
/// @code   fp80Disk d; @endcode
/// @code   double dx = d -> x; @endcode
typedef fp80d_struct fp80d[1];

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Constructors of disks
// //////////////////////////////////////////////////////////////////////////////////////////

/// @brief Transfer one @c fp80d disk into another one.
///
/// @param d target @c fd80d disk
/// @param s source @c fd80d disk
void fp80d_set(fp80d d, fp80d s);

/// @brief Initializes an @c fp80d disk by specifying its center and radius.
///
/// The modulus of x + i y is computed automatically.
/// @param d target @c fd80d disk
/// @param x  real part of the center (@c long_double)
/// @param y imaginary part of the center (@c long_double)
/// @param r radius of the disk (@c long_double)
void fp80d_setl(fp80d d, ldbl x, ldbl y, ldbl r);

/// Initialize an @c fp80d disk by specifying its center and radius.
///
/// The modulus of x + i y is computed automatically.
///
/// @param d target @c fd80d disk
/// @param c  the center
/// @param r radius of the disk (@c long_double)
void fp80d_set80(fp80d d, fp80 c, ldbl r);

/// @brief Initializes a (tiny) @c fp80d disk by specifying its center.
///
/// The radius is initialized to 1 ULP of | x + i y |.
///
/// The modulus of x + i y is computed automatically.
/// @param d target @c fd80d disk
/// @param x  real part of the center (@c long_double)
/// @param y imaginary part of the center (@c long_double)
void fp80d_setUlp(fp80d d, ldbl x, ldbl y);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Disk operations
// //////////////////////////////////////////////////////////////////////////////////////////

/// @brief Sum of two @c fp80d disks.
///
/// The geometric sum of two disks is the disk centered at the sum of the centers,
/// and whose radius is the sum of the radii.
///
/// The result is stored in @c d. The disks @c d, @c op1 and @c op2 do not need to be distinct.
///
/// The modulus of the new center is computed automatically.
///
/// @param d disk to store the result in (@c fp80d)
/// @param op1 first disk (@c fp80d)
/// @param op2 second disk (@c fp80d)
void fp80d_add(fp80d d, fp80d op1, fp80d op2);

/// @brief Add a @c double real number to an @c fp80d disk.
///
/// Translate the disk @c op1 along the real axis. Its radius is not changed.
///
/// The result is stored in @c d. The disks @c d and @c op1 do not need to be distinct.
///
/// The modulus of the new center is computed automatically.
///
/// @param d disk to store the result in (@c fp80d)
/// @param op1 first disk (@c fp80d)
/// @param a real number (@c long_double)
void fp80d_addd(fp80d d, fp80d op1, ldbl a);

/// @brief Difference of two @c fp80d disks.
///
/// The geometric difference of two disks is the disk centered at the difference of the centers,
/// and whose radius is the sum of the radii.
///
/// The result is stored in @c d. The disks @c d, @c op1 and @c op2 do not need to be distinct.
///
/// The modulus of the new center is computed automatically.
///
/// @param d disk to store the result in (@c fp80d)
/// @param op1 first disk (@c fp80d)
/// @param op2 second disk (@c fp80d)
void fp80d_sub(fp80d d, fp80d op1, fp80d op2);

/// Multiplication of two @c fp80d disks.
///
/// The geometric product of D(z1,r1) and D(z2,r2) is the disk centered at z1 * z2,
/// and whose radius is r1*|z2| + |z1|*r2 + r1*r2.
///
/// The result is stored in @c d. The disks @c d, @c op1 and @c op2 do not need to be distinct.
/// The modulus of the new center is computed automatically.
///
/// @param d disk to store the result in (@c fp80d)
/// @param op1 first disk (@c fp80d) D(z1,r1)
/// @param op2 second disk (@c fp80d) D(z2,r2)
void fp80d_mul(fp80d d, fp80d op1, fp80d op2);

/// @brief Inversion of a @c fp80d disk.
///
/// The geometric quotient of 1 and D(z, r) is the disk centered at 1 / z,
/// and whose radius is r / (|z|^2 - |z|r).
///
/// The result is stored in @c d. The disks @c d and @c op do not need to be distinct.
/// The modulus of the new center is computed automatically.
///
/// @param d disk to store the result in (@c fp80d)
/// @param op1 first disk (@c fp80d) D(z1,r1)
/// @param op2 second disk (@c fp80d) D(z2,r2)
///
/// @return @ref true if @c op does not contain zero, @ref false otherwise
bool fp80d_inv(fp80d d, fp80d op);

/// @brief Division of two @c fp80d disks.
///
/// The geometric quotient of D(z1,r1) and D(z2,r2) is the disk centered at z1 / z2,
/// and whose radius is ???.
///
/// The result is stored in @c d. The disks @c d, @c op1 and @c op2 do not need to be distinct.
/// The modulus of the new center is computed automatically.
///
/// @param d disk to store the result in (@c fp80d)
/// @param op1 first disk (@c fp80d) D(z1,r1)
/// @param op2 second disk (@c fp80d) D(z2,r2)
///
/// @return @ref true if @c op2 does not contain zero, @ref false otherwise
bool fp80d_div(fp80d d, fp80d op1, fp80d op2);

/// @brief Multiplication of two @c fp80d disks followed by a scaling transform (homothety)
///
/// The geometric product of @c D(z1,r1) and @c D(z2,r2) is the disk centered at @c z1*z2,
/// and whose radius is @c r1*|z2|+|z1|*r2+r1*r2. The resulting disk is then
/// rescaled by a factor @c a, i.e. the center and the radius are both multiplied by @c a and
/// respectively by @c |a|.
///
/// The result is stored in @c d. The disks @c d, @c op1 and @c op2 do not need to be distinct.
///
/// The modulus of the new center is computed automatically.
/// @param d disk to store the result in (@c fp80d)
/// @param op1 first disk (@c fp80d) D(z1,r1)
/// @param op2 second disk (@c fp80d) D(z2,r2)
/// @param a scaling factor (@c long_double)
void fp80d_muld(fp80d d, fp80d op1, fp80d op2, ldbl a);

/// @brief Square of an @c fp80d disk.
///
/// The geometric square of D(z,r) is the disk centered at z^2
/// and whose radius is 2*r*|z| + r^2.
///
/// The result is stored in @c d. The disks @c d and @c op1 do not need to be distinct.
///
/// The modulus of the new center is computed automatically.
///
/// Optimized code compared to @c fp80Disk_mul(d,op1,op1).
/// @param d disk to store the result in (@c fp80d)
/// @param op1 disk (@c fp80d) D(z,r)
void fp80d_sqr(fp80d d, fp80d op1);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Distance related
// //////////////////////////////////////////////////////////////////////////////////////////

/// Minimal modulus of points in the disk
///
/// Returns a lower bound of the minimal modulus of points in the disk,
/// given by the positive part of |z| - r.
///
/// @param op1 disk (@c fp80d) D(z,r)
ldbl fp80d_min_mod(fp80d op1);

/// Maximal modulus of points in the disk
///
/// Returns an upper bound to the maximal modulus of points in the disk,
/// given by |z| + r.
///
/// @param op1 the disk
ldbl fp80d_max_mod(fp80d op1);

/// Test if the @c fp80d disk @c op1 contains the disk @c op2
///
/// The result is valid if both disks are considered open or if both are considered as closed.
///
/// Returns 1 if the second disk is included in the first one. A return value 0 suggests that
/// the inclusion is not true; however, false negatives are possible.
///
/// TODO: check false negative statement
/// @param op1 first disk (@c fp80d)
/// @param op2 second disk (@c fp80d)
///
/// @return @ref true if @c op1 contains @c op2, @ref false otherwise
bool fp80d_contains(fp80d op1, fp80d op2);

bool fp80d_intersect(fp80d op1, fp80d op2);

#endif /* fp80Disk_h */
