//
//  arithmetic.h
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2022.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2022.
//
//  Copyright © 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the Lesser GNU Public Licence v3.0
//

/**
 \dir algebra
 \brief This folder contains functions related to prime numbers and polynomials.
 */

/**
 \file arithmetic.h
 \brief A collection of basic arithmetic functions.
 
 The functions in this class are needed for counting hyperbolic centers and pre-periodic (or Misiurewicz) parameters. They
 are @b not optimized for general use, in number theoretic context, for example.
*/

#ifndef arithmetic_h
#define arithmetic_h

#include <stdio.h>

#include "ntypes.h"

/// @brief Returns the value of Mobius' function @c mu(n) for @c n<100.
///
/// @param n the argument
/// @return @c mu(n)
int mu(ulong n);

/// @brief Checks if @c p is a prime number.
///
/// @param p the number to check
///
/// @return @ref true if @c p is prime, @ref false otherwise
bool is_prime(ulong p);

/// @brief Returns the smallest prime factor of @c n, or @c n if @c n < 2.
///
/// @param n the number to factor
///
/// @return the smallest prime factor of @c n
ulong factor(ulong n);

/// @brief Returns Euler's totient function of @c n, @c phi(n):=\#k=1,...,n-1
/// that are coprime with @c n, i.e. @c gcd(k,n)=1.
///
/// @param n the number
///
/// @return the number of @c k=1,...,n-1 that are coprime with @c n, i.e. @c gcd(k,n)=1
ulong phi(ulong n);

/// @brief Returns the @c gcd(a, b).
///
/// @param a first argument
/// @param b second argument
///
/// @return the greatest common divisor of @c a and @c b
ulong gcd(ulong a, ulong b);

/// @brief Prints the first @c count values of @c mu().
///
/// @param count the number of values
void print_mu(int count);

/// @brief Prints the first @c count prime numbers.
///
/// @param count the number of prime numbers
void print_primes(int count);

/// @brief Prints the first @c count values of @c phi().
///
/// @param count the number of values
void print_phi(int count);

#endif /* arithmetic_h */
