//
//  ntypes.h
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
 \dir numbers
 \brief This folder contains type definitions, data structures and basic functions to compute with complex numbers.
 
 There are two types of numbers : @c fp80 numbers are used for quick computations, and arbitary precision numbers that are
 based on [mpfr] (https://www.mpfr.org). For each number type, the computations can either be done with basic numbers
 or with disks, which are slower. The only way to obtain certified results is to use the slowest version of disks, @c mpd_struct,
 which are defined in @c mpd.h.
 
 The type @ref u128 is a fixed point complex number with fixed precision of 128 bits, of which 2 are reserved for the integer part.
 They represent points from the square \f$ [0,4) \times [0,4) \f$ that is mapps by the conversion functions to the square
 \f$ [-2,2) \times [0,4) \f$.
 
 The type @ref mpi are real intervals for performing interval arithmetic.
 
 The type @ref mpv is a vector of numbers with the same precision, which can be used for real numbers @c mpfr_t, complex numbers
 @ref mpc or real intervals @ref mpi. This type of vector trades a bit of speed for lower memory usage. Real or complex numbers can
 @b attach to some position in the vector, that is, the memory used for their digits is common with the vector. This is faster if a few operations
 are performed for each element of the vector, as for example when performing the FFT transform, see fft.h. The price to pay is that
 the exponent and the sign has to be synchronized after each operation.
*/

/**
 \file ntypes.h
 \brief Definition of basic types and a function to check that the operating system and the compiler are compatible with this library.
*/

#ifndef ntypes_h
#define ntypes_h

#include <math.h>

/// @c byte is @c uint8
typedef unsigned char byte;

/// @c word is @c uint16
typedef unsigned short word;

/// @c uint is @c uint32
typedef unsigned int uint;

/// @c ulong is @c uint64
typedef unsigned long ulong;

/// @c ldbl is @c long @c double, that is @c fp80 real
typedef long double ldbl;

/// Logic type @c bool can take values \ref true or \ref false
typedef byte bool;

/// Easy to copy ten bytes.
typedef struct {
    word b[5]; ///< the ten bytes
} tenb;

/// A @ref ldbl squeezed into ten bytes. Slower, only use when memory [bandwidth, chache] is at a premium.
typedef union {
    ldbl val;  ///< the value to compute with
    tenb bulk; ///< the ten bytes actually used to store the value
} sq80;

/// \f$ \pi \f$ with fp80 precision
#define PI (3.14159265358979323846264338327950288L)

/// \f$ log(2) \f$ with fp80 precision
#define LOG2 (0.693147180559945309417232121458176568L)

/// \f$ log(10) \f$ with fp80 precision
#define LOG10 (2.302585092994045684017991454684364208L)

/// \f$ sqrt(2) \f$ with fp80 precision
#define SQRT2 (1.414213562373095048801688724210L)

/// \f$ sqrt(3) \f$ with fp80 precision
#define SQRT3 (1.732050807568877293527446341506L)

/// Boolean value true.
#define true   1

/// Boolean value false.
#define false  0

/// @brief Check if the operating system and the compiler are compatible with this library.
///
/// @return @ref true if successfull, @ref false otherwise
bool ntypes_check(void);

bool is_number(ldbl x);

#define aabs(x) (x < 0 ? -x : x)
#define sgn(x) (x < 0 ? -1 : x > 0 ? 1 : 0)

#endif /* ntypes_h */
