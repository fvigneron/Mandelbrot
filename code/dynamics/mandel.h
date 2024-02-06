//
//  mandel.h
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
  \dir dynamics
  \brief This folder contains data structures and functions related to the Mandelbrot set and Julia sets.
 */

/**
 \file mandel.h
 \brief A collection of methods related to the Mandelbrot set M.
 
  Two precisions are used: @c long @c double (@c ldbl) and arbitrary (@c mpfr, @c mpc).
 
  \warning This module is not thread-safe, as it uses global buffers for optimization.
*/

#ifndef Mandelbrot_h
#define Mandelbrot_h

#include <math.h>
#include <mpfr.h>

#include "fp80.h"
#include "fp80d.h"
#include "mpc.h"
#include "mpd.h"

// some powers of two
#define T50 1.125899906842624e15 
#define T54 1.801439850948198e16
#define T56 7.205759403792794e16
#define T60 1.152921504606847e18
#define T64 1.844674407370955e19
#define T66 7.378697629483821e19
#define T70 1.180591620717411e21

#define OUT2 1E10
#define MANDEL_MAX_MIS  50

// MARK: stats & housekeeping

/// @brief Returns the number of iterates for the last operation, when applicable.
///
/// @b see mandel_distl(), mandel_dist(), mandel_rootl(), mandel_root(), mandel_rootRefl(), mandel_rootRef(),
/// mandel_solutionl(), mandel_solution(), mandel_solutionRefl(), mandel_solutionRef(),
ulong mandel_last_iter(void);

/// @brief Return the total number of Newton terms computed, including in fp80 precision.
/// Can be used for statistics.
///
/// @return the total number of Newton terms computed
ulong mandel_tot_nt(void);

// MARK: quick versions with low precision (long double, fp80)

/// @brief Computes quickly and with low precision a lower bound for the distance to the Mandelbrot set of the point @c c,
/// with at most @c maxIter iterates.
///
/// @param c the parameter
/// @param maxIter maximal number of iterates to use
///
/// @return the distance to the Mandelbrot set, which is @c 0 if @c c does not escape in at most @c maxIter iterates.
ldbl mandel_distl(fp80 c, int maxIter);

/// @brief Computes quickly and with low precision the value @c v:=p_{per}(c).
///
/// @param v the resulted value at c
/// @param c the argument
/// @param per the period, or index of @c p_n
void mandel_vall(fp80 v, fp80 c, int per);

/// @brief Computes quickly with low precision the value @c v:=p_{pp+per}(c)-p_{pp}(c).
///
/// @param v the resulted value at c
/// @param c the argument
/// @param pp the pre-period
/// @param per the period, or index of @c p_n
void mandel_mis_vall(fp80 v, fp80 c, int pp, int per);

/// @brief Computes quickly with low precision the value
/// @c v:=(p_{pp+per}(c)-p_{pp}(c))/(p_{pp+per-1}(c)-p_{pp-1}(c)).
///
/// @param v the resulted value at c
/// @param c the argument
/// @param pp the pre-period
/// @param per the period
void mandel_miss_vall(fp80 v, fp80 c, int pp, int per);

/// @brief Computes quickly and with low precision the value @c v:=p_{per}(c) and the derivative @c d:=p'_{per}(c).
///
/// @param v the resulted value at c
/// @param d the derivative at c
/// @param c the argument
/// @param per the period, or index of @c p_n
///
/// @return @ref true the computation succeeded, @ref false otherwise
bool mandel_val_derl(fp80 v, fp80 d, fp80 c, int per);

/// @brief Computes quickly and with low precision the value @c v:=p_{per}(c) and the derivatives
/// @c d_k:=p^(k)_{per}(c), for k = 1 ... @c ders.
///
/// @warning @c ders should be at least @c 2 and at most @c 67.
///
/// @param c the argument
/// @param per the period, or index of @c p_n
/// @param ders the number of derivatives to compute
///
/// @return the list of the successive derivatives, @c NULL if some error occurred
fp80_ptr mandel_val_dersl(fp80 c, int per, int ders);

/// @brief Computes quickly with low precision the value @c v:=p_{pp+per}(c)-p_{pp}(c)
/// and the derivative @c d:=p'_{pp+per}(c)-p'_{pp}(c).
///
/// @param v the resulted value at c
/// @param d the derivative at c
/// @param c the argument
/// @param pp the pre-period
/// @param per the period, or index of @c p_n
void mandel_mis_val_derl(fp80 v, fp80 d, fp80 c, int pp, int per);

/// @brief Computes quickly with low precision the Newton term of @c p:=p_{per} at c, that is @c nt:=c-N_p(c)=p(c)/p'(c).
///
/// It computes correctly even if the iterates would grow to numbers outside the range of representable numbers,
/// that is if @c |log_10(|p_{per}(c))| > 33.
///
/// @param nt the result value of the Newton term if the returned value is 1
/// @param c the point in which to compute
/// @param per the period, or index of @c p_n
///
/// @return @ref false if the derivative is @c 0, @ref true otherwise
bool mandel_ntl(fp80 nt, fp80 c, int per);

/// @brief Computes quickly with low precision the Newton term of @c p:=p_{pp+per}(c)-p_{pp}(c)-t at c,
/// that is @c nt:=c-N_p(c)=p(c)/p'(c).
///
/// @param nt the result value of the Newton term if the returned value not @c 0
/// @param c the point in which to compute
/// @param t the target value
/// @param pp the pre-period
/// @param per the period
///
/// @return the square of the modulus of the derivative of @c p at @c c, @c |p'(c)|^2, @c 0 if some error occurred
ldbl mandel_mis_ntl(fp80 nt, fp80 c, fp80 t, int pp, int per);

/// @brief Computes quickly with low precision the Newton term of
/// @c p:=(p_{pp+per}(c)-p_{pp}(c))/(p_{pp+per-1}(c)-p_{pp-1}(c))-t at c,
/// that is @c nt:=c-N_p(c)=p(c)/p'(c).
///
/// @param nt the result value of the Newton term if the returned value is not @c 0
/// @param c the point in which to compute
/// @param t the target value, zero is used if @c NULL
/// @param pp the pre-period
/// @param per the period
///
/// @return the square of the modulus of the denominator, @c 0 if some error occurred
ldbl mandel_miss_ntl(fp80 nt, fp80 c, fp80 t, int pp, int per);

/// @brief Computes quickly with low precision the root to which the Newton method of
/// @c p_n:=p_{per} converges with starting point @c c.
///
/// @param root the root if the @c return value is @ref true
/// @param c the strating point
/// @param per the period, or index of @c p_n
/// @param maxIter the maximum number of iterates
/// @param error the maximum modulus of the last Newton term for the root to be validated
///
/// @return @ref true  if a root was found (up to @c error), @ref false otherwise
bool mandel_rootl(fp80 root, fp80 c, int per, int maxIter, ldbl error);

/// @brief Computes quickly with low precision the root to which the Newton method of
/// @c p_n:=p_{per} converges with starting point @c c.
/// It performs @c refine extra Newton iterates to refine the result
///
/// @param root the root if the @c return value is @ref true
/// @param c the strating point
/// @param per the period, or index of @c p_n
/// @param maxIter the maximum number of iterates
/// @param error the maximum modulus of the last Newton term for the root to be validated
/// @param refine the number of extra iterates after the root was found
/// @param maxDist max distance allowed from @c to @c root
///
/// @return @ref true if a root was found (up to @c error), @ref false otherwise
bool mandel_root_refl(fp80 root, fp80 c, int per, int maxIter, ldbl error, int refine, double maxDist);

/// @brief Computes quickly with low precision the newton term of @c p_n-t at the point @c c and stores it in @c nt.
///
/// @param nt the value of the Newton term
/// @param c the argument of @c p_n
/// @param t the target value
/// @param n the period
///
/// @return the square of the modulus of the derivative of @c p_n at @c c, @c |p_n'(c)|^2
///
/// If we define \f$ p_n(c) = a_n + i b_n \f$ and \f$ p_n'(c) = a_n' + i b_n' \f$ the recurrence relations are
/// \f[ a_{n+1} = a_n^2 - b_n^2 + c_x \f]
/// \f[ b_{n+1} = 2 a_n b_n + c_y \f]
/// \f[ a_{n+1}' = 2(a_n a_n' - b_n b_n') + 1 \f]
/// \f[ b_{n+1}' = 2(a_n b_n' + b_n a_n')\f]
/// and the Newton term is given by the expansion of:
/// \f[ N = (a_n + i b_n - t)(a_n' - i b_n') / ( (a_n')^2 + (b_n')^2) \f]
ldbl mandel_nt_soll(fp80 nt, fp80 c, fp80 t, int n);

/// @brief Computes quickly with low precision the solution of @c p_{per}(v)=t to which the Newton method
/// converges with starting point @c c.
///
/// @param sol the root if the @c return value is @ref true
/// @param err the modulus of the last Newton term, normally larger than the error of the result 
/// @param c the strating point
/// @param t the target value
/// @param per the period, or index of @c p_n
/// @param maxIter the maximum number of iterates
/// @param error the maximum modulus of the last Newton term for the root to be validated
///
/// @return @ref true if a root was found (up to @c error), @ref false otherwise
bool mandel_solutionl(fp80 sol, ldbl *err, fp80 c, fp80 t, int per, int maxIter, ldbl error);

/// @brief Computes quickly with low precision the solution of @c p_{per}(v)=t to which the Newton method
/// converges with starting point @c c. It performs @c refine extra Newton iterates to refine the result
///
/// @param root the root if the @c return value is @ref true
/// @param err the modulus of the last Newton term, normally larger than the error of the result
/// @param c the starting point
/// @param t the target value
/// @param per the period, or index of @c p_n
/// @param maxIter the maximum number of iterates
/// @param error the maximum modulus of the last Newton term for the root to be validated
/// @param refine the number of extra iterates after the root was found
///
/// @return @ref true if a root was found (up to @c error), @ref false otherwise
bool mandel_sol_refl(fp80 root, ldbl *err, fp80 c, fp80 t,
                        int per, int maxIter, ldbl error, int refine);

/// @brief Computes quickly with low precision the root to which the Newton method of
/// @c p_{pp+per}-p_{pp} converges with starting point @c c.
///
/// @param root the root if the @c return value is @ref true
/// @param c the strating point
/// @param pp the pre-period
/// @param per the period
/// @param maxIter the maximum number of iterates
/// @param error the maximum modulus of the last Newton term for the root to be validated
///
/// @return @ref true  if a root was found (up to @c error), @ref false otherwise
bool mandel_mis_rootl(fp80 root, fp80 c, int pp, int per, int maxIter, ldbl error);

/// @brief Computes quickly with low precision the root to which the Newton method of
/// @c p_{pp+per}-p_{pp}-t converges with starting point @c c.
///
/// It performs @c refine extra Newton iterates to refine the result
///
/// @param root the root if the @c return value is @ref true
/// @param c the strating point
/// @param t the target value
/// @param pp the pre-period
/// @param per the period
/// @param maxIter the maximum number of iterates
/// @param error the maximum modulus of the last Newton term for the root to be validated
/// @param refine the number of extra iterates after the root was found
/// @param maxDist max distance allowed from @c to @c root
///
/// @return @ref true if a root was found (up to @c error), @ref false otherwise
bool mandel_mis_root_refl(fp80 root, fp80 c, fp80 t, int pp, int per, int maxIter, ldbl error, int refine, double maxDist);

/// @brief Computes quickly with low precision the root to which the Newton method of
/// @c (p_{pp+per}(c)-p_{pp}(c))/(p_{pp+per-1}(c)-p_{pp-1}(c))-t converges with starting point @c c.
///
/// @param root the root if the @c return value is @ref true
/// @param c the strating point
/// @param pp the pre-period
/// @param per the period
/// @param maxIter the maximum number of iterates
/// @param error the maximum modulus of the last Newton term for the root to be validated
///
/// @return @ref true  if a root was found (up to @c error), @ref false otherwise
bool mandel_miss_rootl(fp80 root, fp80 c, int pp, int per, int maxIter, ldbl error);

/// @brief Computes quickly with low precision the root to which the Newton method of
/// @c (p_{pp+per}(c)-p_{pp}(c))/(p_{pp+per-1}(c)-p_{pp-1}(c))-t converges with starting point @c c.
///
/// It performs @c refine extra Newton iterates to refine the result
///
/// @param root the root if the @c return value is @ref true
/// @param c the strating point
/// @param t the target value
/// @param pp the pre-period
/// @param per the period
/// @param maxIter the maximum number of iterates
/// @param error the maximum modulus of the last Newton term for the root to be validated
/// @param refine the number of extra iterates after the root was found
/// @param maxDist max distance allowed from @c to @c root
///
/// @return @ref true if a root was found (up to @c error), @ref false otherwise
bool mandel_miss_root_refl(fp80 root, fp80 c, fp80 t, int pp, int per, int maxIter, ldbl error, int refine, double maxDist);

/// @brief Quick version of the proof given by mandel_isHyp(). Used only for optimization, not for certification.
///
/// @param root the center of the disk where the root should be
/// @param per the period, or the index n of the polynomial p_n
/// @param error the radious of the disk where the root should be
///
/// @return @ref true if the disk contains a root, @ref false otherwise
bool mandel_is_hypl(fp80 root, int per, ldbl error);

/// @brief Quick version of the proof given by  mandel_isMis(). Used only for optimization, not for certification.
///
/// @param c the center of the disk where the pre-periodic point should be
/// @param pp the pre-period
/// @param per the period
/// @param error the radious of the disk where the root should be
///
/// @return @ref true if the disk contains a pre-periodic point of exact pre-period pp and period per, @ref false otherwise
bool mandel_is_misl(fp80 c, int pp, int per, ldbl error);

/// @brief Quick version of the proof given by mandel_convNp(). Used only for optimization, not for certification.
///
/// This method is based on Theorem 2 of [1].
///
/// @see see also isHypl()
///
/// @param root the best known approximation of the root of p_n
/// @param per the period n of the polynomial p_n
/// @param rootError the error with which the root is known (proven)
/// @param radius the radius of the disk to check
///
/// @return @ref true if  D(root, radius) is included in the attraction bassin of N_{p_per}, @ref false otherwise
bool mandel_conv_npl(fp80 root, int per, ldbl rootError, ldbl radius);

/// @brief Quick version of the proof given by mandel_convNpp(). Used only for optimization, not for certification.
///
/// This method is based on Theorem 2 of [1].
///
/// @see see also isHypl()
///
/// @param root the best known approximation of the root of p_{pp + per} - p_{pp}
/// @param pp the pre-period
/// @param per the period
/// @param rootError the error with which the root is known (proven)
/// @param radius the radius of the disk to check
/// @param checkMultiple true to check convergence for multiple roots (hyperbolic centers)
///
/// @return @ref true if  D(root, radius) is included in the attraction bassin of N_{p_per}, @ref false otherwise
bool mandel_conv_nppl(fp80 c, int pp, int per, ldbl rootError, ldbl radius, bool checkMultiple);

/// @brief Quicker and less reliable version of isHypl(). Used only for optimization, not for certification.
///
/// @param root the center of the disk where the root should be
/// @param per the period, or the index n of the polynomial p_n
/// @param error the radious of the disk where the root should be
/// @param refine the number of Newton iterates to refine the root
///
/// @return @ref true if the disk contains a root, @ref false otherwise
bool mandel_is_prob_hypl(fp80 root, int per, ldbl error, int refine);

/// @brief Returns @c max(|x-nx|,|y-ny|)<=eps.
bool mandel_eql(ldbl x, ldbl y, ldbl nx, ldbl ny, ldbl eps);

/// @brief Computes quickly with low precision the @c deg+1 trailing coefficients of @c p_{per}.
///
/// @param per the period, or the index of @c p_n
/// @param deg the degree of the truncated polynomial of @c p_{per}
///
/// If @c deg exceeds the degree of @c p_{per}, the extra coefficients are zero.s
/// @warning @c per should be at least @c 2 and at most @c 50.
///
/// @return the list of coefficients in a @c deg+1 array of @c ldbl.
ldbl* mandel_coeffs_hypl(int per, int deg);

/// @brief Computes quickly with low precision the @c deg+1 trailing coefficients
/// of @c m_{pp,per} @c = @c p_{pp+per}-p_{pp}.
///
/// @param pp the pre-period, or the first index of @c m_{k,n}
/// @param per the period, or the second index of @c m_{k,n}
/// @param deg the degree of the truncated polynomial of @c m_{pp,per}
///
/// If @c deg exceeds the degree of @c m_{pp,per}, the extra coefficients are zero.
/// @warning @c pp+per should be at least @c 2 and at most @c 50.
///
/// @return the list of coefficients in a @c deg+1 array of @c ldbl.
ldbl* mandel_coeffs_misl(int pp, int per, int deg);

/// @brief For each @c k=1,...,deg, computes quickly with low precision the sum of @c r^{-k},
/// where @c r ranges over all non zero roots of @c p_{per}.
///
/// Returns
/// \f[
/// \sum_{p_{per}(r)=0, r\neq0} r^{-k}
/// \f]
///
/// The sum is computed from the trailing coefficients, using @c deg + 2 trailing coefficients of @c p_{per}
/// and Newton's identities.
/// @warning @c per should be at least @c 1 and at most @c 50.
///
/// @param per the period, or the index of @c p_n
/// @param deg the largest power to use for the sums
///
/// @return the list of @c deg sums
ldbl* mandel_sum_neg_pows_hypl(int per, int deg);

/// @brief For each @c k=1,...,deg, computes quickly with low precision the sum of @c r^{-k},
/// where @c r ranges over all non zero roots of @c m_{pp,per} @c = @c p_{pp+per}-p_{pp}.
///
/// \f[
/// \text{Returns } \sum_{m_{pp,per}(r)=0, r\neq0} r^{-k}
/// \f]
///
/// The sum is computed from the trailing coefficients, using @c deg + 2 trailing coefficients of @c m_{pp,per}
/// and Newton's identities.
/// @warning @c pp+per should be at least @c 1 and at most @c 50.
///
/// @param pp the pre-period, or the first index of @c m_{k,n}
/// @param per the period, or the second index of @c m_{k,n}
/// @param deg the largest power to use for the sums
///
/// @return the list of @c deg sums
ldbl* mandel_sum_neg_pows_misl(int pp, int per, int deg);

// MARK: arbitrary precision versions

/// @brief Computes a lower bound for the distance to the Mandelbrot set of the point @c c, with at most @c maxIter iterates.
/// Stores the value in @c d, which is @c 0 if @c c does not escape in at most @c maxIter iterates.
/// This method does not certify the result.
///
/// @param d the distance to the Mandelbrot set
/// @param c the parameter
/// @param maxIter maximal number of iterates to use
///
/// @return the number of iterates needed to escape (a larger disk than D(0,25))
long mandel_dist(mpfr_t d, mpc c, long maxIter);

/// @brief Computes the value @c v:=p_{per}(c).
///
/// @param v the resulted value at c
/// @param c the argument
/// @param per the period, or index of @c p_n
void mandel_val(mpc v, mpc c, int per);

/// @brief Computes the value @c v:=p_{pp+per}(c)-p_{pp}(c).
///
/// @param v the resulted value at c
/// @param c the argument
/// @param pp the pre-period
/// @param per the period, or index of @c p_n
void mandel_mis_val(mpc v, mpc c, int pp, int per);

/// @brief Computes the value @c v:=p_{per}(c) and the derivative @c d:=p'_{per}(c).
///
/// @param v the resulted value at c
/// @param d the derivative at c
/// @param c the argument
/// @param per the period, or index of @c p_n
void mandel_val_der(mpc v, mpc d, mpc c, long per);

/// @brief Computes the value @c v:=p_{pp+per}(c)-p_{pp}(c) and the derivative @c d:=p'_{pp+per}(c)-p'_{pp}(c).
///
/// @param v the resulted value at c
/// @param d the derivative at c
/// @param c the argument
/// @param pp the pre-period
/// @param per the period, or index of @c p_n
void mandel_mis_val_der(mpc v, mpc d, mpc c, int pp, int per);

/// @brief Computes the value @c v:=p_{pp+per-1}(c)+p_{pp-1}(c) and the derivative
/// @c d:=p'_{pp+per-1}(c)+p'_{pp-1}(c).
///
/// @param v the resulted value at c
/// @param d the derivative at c
/// @param c the argument
/// @param pp the pre-period
/// @param per the period, or index of @c p_n
void mandel_miss_val_der(mpc v, mpc d, mpc c, int pp, int per);

/// @brief Computes the Newton term of @c p:=p_{per} at @c c, that is @c c-N_p(c)=p(c)/p'(c).
///
/// It computes correctly even if the iterates would grow to numbers outside the range of representable numbers,
/// that is if @c |log_2(|p_{per}(c))| > 10^9.
///
/// @param v the result value of the Newton term
/// @param c the point in which to compute
/// @param per the period, or index of @c p_n
///
/// @return @ref false if the derivative is @c 0, @ref true otherwise
bool mandel_nt(mpc v, mpc c, int per);

/// @brief Computes the Newton term of @c p:=p_{pp+per}(c)-p_{pp}(c) at c, that is @c c-N_p(c)=p(c)/p'(c).
///
/// @param v the result value of the Newton term
/// @param c the point in which to compute
/// @param pp the pre-period
/// @param per the period
///
/// @return @ref false if the derivative is @c 0, @ref true otherwise
bool mandel_mis_nt(mpc v, mpc c, int pp, int per);

/// @brief Computes the root to which the Newton method of @c p_n:=p_{per} converges with starting point @c c.
///
/// @param v the root
/// @param c the strating point
/// @param per the period, or index of @c p_n
/// @param maxIter the maximum number of iterates
///
/// @return @ref true if a root was found (up to 3/4 of the precision of @c v), @ref false otherwise
bool mandel_root(mpc v, mpc c, int per, int maxIter);

/// @brief Computes the root to which the Newton method of @c p_n:=p_{per} converges with starting point @c c.
///
/// It performs @c refine extra Newton iterates to refine the result
///
/// @param v the root
/// @param c the strating point
/// @param per the period, or index of @c p_n
/// @param maxIter the maximum number of iterates
/// @param refine the number of extra iterates after the root was found
///
/// @return @ref true if a root was found (up to 3/4 of the precision of @c v), @ref false otherwise
bool mandel_root_ref(mpc v, mpc c, int per, int maxIter, int refine);

/// @brief Computes the Newton term of @c p_{per}-t at c, that is @c c-N_{p_n+t}}(c)=(p_n(c)-t)/p'_n(c).
///
/// @param v the result value of the Newton term
/// @param c the point in which to compute
/// @param t the target value
/// @param per the period, or index of @c p_n
///
/// @return @ref false if the derivative is @c 0, @ref true otherwise
bool mandel_nt_sol(mpc v, mpc c, mpc t, long per);

/// @brief Computes the Newton term of @c p:=p_{pp+per}(c)-p_{pp}(c)-t at c,
/// that is @c nt:=c-N_p(c)=p(c)/p'(c).
///
/// @param nt the result value of the Newton term if the returned value is 1
/// @param c the point in which to compute
/// @param t the target value, @c NULL to use @c 0
/// @param pp the pre-period
/// @param per the period
///
/// @return @ref false if the derivative is @c 0, @ref true otherwise
bool mandel_mis_nt_sol(mpc nt, mpc c, mpc t, int pp, int per);

/// @brief Computes the Newton term of @c p:=(p_{pp+per}(c)-p_{pp}(c))/(p_{pp+per-1}(c)-p_{pp-1}(c))-t at c,
/// that is @c nt:=c-N_p(c)=p(c)/p'(c).
///
/// @param nt the result value of the Newton term if the returned value is 1
/// @param c the point in which to compute
/// @param t the target value
/// @param pp the pre-period
/// @param per the period
///
/// @return @ref false if the derivative is @c 0, @ref true otherwise
bool mandel_miss_nt_sol(mpc nt, mpc c, mpc t, int pp, int per);

/// @brief Computes the solution of @c p_{per}(v)=t to which the Newton method converges with starting point @c c.
///
/// @param v the root
/// @param c the strating point
/// @param t the target value
/// @param per the period, or index of @c p_n
/// @param maxIter the maximum number of iterates
///
/// @return @ref true if a root was found (up to 3/4 of the precision of @c v), @ref false otherwise
bool mandel_solution(mpc v, mpc c, mpc t, int per, int maxIter);

/// @brief Computes the solution of @c p_{per}(v)=t to which the Newton method converges with starting point @c c.
///
/// It performs @c refine extra Newton iterates to refine the result
///
/// @param v the root
/// @param err the modulus of the last Newton term, normally >> the @c error of the result
/// @param c the strating point
/// @param t the target value
/// @param per the period, or index of @c p_n
/// @param maxIter the maximum number of iterates
/// @param refine the number of extra iterates after the root was found
///
/// @return @ref true if a root was found (up to 3/4 of the precision of @c v), @ref false otherwise
bool mandel_sol_ref(mpc v, mpfr_t err, mpc c, mpc t, long per, int maxIter, int refine);

/// @brief Computes the root to which the Newton method of
/// @c p_{pp+per}-p_{pp}-t converges with starting point @c c.
///
/// It performs @c refine extra Newton iterates to refine the result
///
/// @param v the root if the @c return value is @ref true
/// @param c the strating point
/// @param t the target value, @c NULL to use @c 0
/// @param pp the pre-period
/// @param per the period
/// @param maxIter the maximum number of iterates
/// @param refine the number of extra iterates after the root was found
///
/// @return @ref true if a root was found (up to @c error), @ref false otherwise
bool mandel_mis_root_ref(mpc v, mpc c, mpc t, int pp, int per, int maxIter, int refine);

/// @brief Computes the root to which the Newton method of
/// @c (p_{pp+per}(c)-p_{pp}(c))/(p_{pp+per-1}(c)-p_{pp-1}(c))-t converges with starting point @c c.
///
/// It performs @c refine extra Newton iterates to refine the result
///
/// @param v the root if the @c return value is @ref true
/// @param c the strating point
/// @param t the target value
/// @param pp the pre-period
/// @param per the period
/// @param maxIter the maximum number of iterates
/// @param refine the number of extra iterates after the root was found
///
/// @return @ref true if a root was found (up to @c error), @ref false otherwise
bool mandel_miss_root_ref(mpc v, mpc c, mpc t, int pp, int per, int maxIter, int refine);

/// @brief Returns \f$ |a-b|_1<error \f$, where \f$ |x+iy|_1:=max(|x|,|y|). \f$
///
/// @param a a complex number
/// @param b another complex number
/// @param error max distance between @c a and @c b
///
/// @return @ref true if @c equals @c b, up to @c error on each coordinate, @ref false otherwise
bool mandel_eq(mpc a, mpc b, mpfr_t error);

/// @brief Computes the @c deg+1 trailing coefficients of @c p_{per}.
///
/// @param per the period, or the index of @c p_n
/// @param deg the degree of the truncated polynomial of @c p_{per}
/// @param prec the precistion of the results, in bits
///
/// If @c deg exceeds the degree of @c p_{per}, the extra coefficients are zero.
/// @warning @c per should be at least @c 1 and at most @c 1000.
///
/// @return the list of coefficients  in a @c deg+1 array of @c __mpfr_struct.
// FIXME: should it be mpfr_ptr instead ? same below
__mpfr_struct* mandel_coeffs_hyp(int per, int deg, int prec);

/// @brief Computes the @c deg+1 trailing coefficients
/// of @c m_{pp,per} @c = @c p_{pp+per}-p_{pp}.
///
/// @param pp the preperiod, or the first index of @c p_{k,n}
/// @param per the period, or the second index of @c p_{k,n}
/// @param deg the degree of the truncated polynomial of @c p_{pp,per}
/// @param prec the precistion of the results, in bits
///
/// If @c deg exceeds the degree of @c m_{pp,per}, the extra coefficients are zero.
/// @warning @c pp+per should be at least @c 1 and at most @c 1000.
///
/// @return the list of coefficients  in a @c deg+1 array of @c __mpfr_struct.
__mpfr_struct* mandel_coeffs_mis(int pp, int per, int deg, int prec);

/// @brief For each @c k=1,...,deg, computes the sum of @c r^{-k}, where @c r ranges over all non-zero roots of @c p_{per}.
///
/// Returns
/// \f[
/// \sum_{p_{per}(r)=0, r\neq0} r^{-k}
/// \f]
///
/// The sum is computed from the trailing coefficients, using @c deg + 2 trailing coefficients of @c p_{per}
/// and Newton's identities.
///
/// @param per the period, or the index of @c p_n
/// @param deg the largest power to use for the sums
/// @param prec the precistion of the results, in bits
///
/// @warning @c per should be at least @c 1 and at most @c 1000.
///
/// @return the list of @c deg sums
__mpfr_struct* mandel_sum_all_neg_pows_hyp(int per, int deg, int prec);

/// @brief For each @c k=1,...,deg, computes the sum of @c r^{-k}, where @c r ranges over non-zero roots of
/// @c p_{per} that are not roots of some other @c p_k.
///
/// Returns
/// \f[
/// \sum_{p_{per}(r)=0, r\neq0} r^{-k}
/// \f]
///
/// The sum is computed from the trailing coefficients, using @c deg + 2 trailing coefficients of @c p_{per}
/// and Newton's identities.
///
/// @param per the period, or the index of @c p_n
/// @param deg the largest power to use for the sums
/// @param prec the precistion of the results, in bits
///
/// @warning @c per should be at least @c 1 and at most @c 1000.
///
/// @return the list of @c deg sums
__mpfr_struct* mandel_sum_neg_pows_hyp(int per, int deg, int prec);

/// @brief For each @c k=1,...,deg, computes the sum of @c r^{-k}, where @c r ranges over all non-zero roots of @c m_{pp,per} @c = @c p_{pp+per}-p_{pp}.
///
/// Returns
/// \f[
/// \sum_{m_{pp,per}(r)=0, r\neq0} r^{-k}
/// \f]
/// The sum is computed from the trailing coefficients, using @c deg + 2 trailing coefficients of @c m_{pp,per}
/// and Newton's identities.
///
/// @warning @c per should be at least @c 1 and at most @c 1000.
///
/// @param pp the pre-period, or the first index of @c m_{k,n}
/// @param per the period, or the second index of @c m_{k,n}
/// @param deg the largest power to use for the sums
/// @param prec the precistion of the results, in bits
__mpfr_struct*  mandel_sum_all_neg_pows_mis(int pp, int per, int deg, int prec);

/// @brief For each @c k=1,...,deg, computes the sum of @c r^{-k}, where @c r ranges over non-zero roots of
/// @c m_{pp,per} @c = @c p_{pp+per}-p_{pp} that are not roots of some other @c m_{k,n} with @c k<=pp and
/// @c n<=per.
///
/// Returns
/// \f[
/// \sum_{m_{pp,per}(r)=0, r\neq0} r^{-k}
/// \f]
/// The sum is computed from the trailing coefficients, using @c deg + 2 trailing coefficients of @c m_{pp,per}
/// and Newton's identities.
///
/// @warning @c per should be at least @c 1 and at most @c 1000.
///
/// @param pp the pre-period, or the first index of @c m_{k,n}
/// @param per the period, or the second index of @c m_{k,n}
/// @param deg the largest power to use for the sums
/// @param prec the precistion of the results, in bits
__mpfr_struct*  mandel_sum_neg_pows_mis(int pp, int per, int deg, int prec);

// MARK: proofs with disk arithmetics

/// @brief Certifies that a given disk D(root, error) contains an unique root of p_n (a center of a hyperbolic component).
///
/// This is based on Theorem 1 in [1].
///
/// @param root the center of the disk where the root should be
/// @param per the period, or the index n of the polynomial p_n
/// @param error the radious of the disk where the root should be
///
/// @return @ref true if the disk contains a root, @ref false otherwise
bool mandel_is_hyp(mpc root, int per, ldbl error);

/// @brief Certifies that a given disk D(root, error) contains a Misiurewicz (or pre-periodic) point of pre-period pp and period per.
///
/// This is based on Theorem 1 in [1].
///
/// @param c the center of the disk where the pre-periodic point should be
/// @param pp the pre-period
/// @param per the period
/// @param error the radious of the disk where the root should be
///
/// @return @ref true if the disk contains a pre-periodic point of exact pre-period pp and period per, @ref false otherwise
bool mandel_is_mis(mpc c, int pp, int per, ldbl error);

/// @brief Certifies that the Newton method of p_{per} converges in the disk D(root, radius), assuming there is a root in
/// the disk D(root, rootError).
///
/// This method is based on Theorem 2 of [1]
///
/// @see see also isRoot()
///
/// @param root the best known approximation of the root of p_n
/// @param per the period n of the polynomial p_n
/// @param rootError the error with which the root is known (proven)
/// @param radius the radius of the disk to check
///
/// @return @ref true if  D(root, radius) is included in the attraction bassin of N_{p_per}, @ref false otherwise
bool mandel_conv_np(mpc root, int per, ldbl rootError, ldbl radius);

/// @brief Certifies that the Newton method of p_{pp + per} - p_{pp} converges in the disk D(root, radius), assuming that
/// it has a root in the disk D(root, rootError).
///
/// This method is based on Theorem 2 of [1]
///
/// @see see also isRoot()
///
/// @param root the best known approximation of the root of p_{pp + per} - p_{pp}
/// @param pp the pre-period
/// @param per the period
/// @param rootError the error with which the root is known (proven)
/// @param radius the radius of the disk to check
/// @param checkMultiple true to check convergence for multiple roots (hyperbolic centers) 
///
/// @return @ref true if  D(root, radius) is included in the attraction bassin of N_{p_{pp + per} - p_{pp}}, @ref false otherwise
bool mandel_conv_npp(mpc root, int pp, int per, ldbl rootError, ldbl radius, bool checkMultiple);

/// @brief Certifies that p_{per} is injective in the given disk.
///
/// This is based on Lemma 2.5 in [2].
///
/// @param disk the disk
/// @param per the period
///
/// @return @ref true if  p_per is injective on disk, @ref false otherwise
bool mandel_univalent(mpd disk, int per);

// MARK: counting functions

/// @brief Returns the number of hyperbolic components of period per (at least 1, at most 65).
///
/// @param per the period
///
/// @return the number of hyperbolic components of period per in [1, 65], 0 otherwise
ulong mandel_hyp_count(int per);

/// @brief Returns the number of pre-periodic (or Misiurewicz) parameters of period per and pre-period pp, for  pp + per in [1, 65].
///
/// This is based on Corollary 3.3 in
/// B. Hutz & A. Towsey. Misiurewicz points for polynomial maps and transversality. New York J. Math. 21 (2015) 297–319.
///
/// @param pp the pre-period
/// @param per the period
///
/// @return the number of pre-periodic (or Misiurewicz) parameters of period per and pre-period pp, pp + per in [1, 65], 0 otherwise
ulong mandel_mis_count(int pp, int per);

/// @brief Returns the number of primitive hyperbolic components of period per (at least 1, at most 65).
///
/// @param per the period
///
/// @return the number of primitive hyperbolic components of period per in [1, 65], 0 otherwise
ulong mandel_prim_hyp_count(int per);

/// @brief Returns the number of non-primitive hyperbolic components of period per (at least 2, at most 129).
///
/// @param per the period
///
/// @return the number of non-primitive hyperbolic components of period per in [2, 129], 0 otherwise
ulong mandel_non_prim_hyp_count(int per);

/// @brief Returns the number of critical points of the multiplier of period @c per (at least 2, at most 129).
///
/// @param per the period
///
/// @return the number of critical points of the multiplier of period @c per in [2, 129], 0 otherwise
ulong mandel_mult_crit_count(int per);

#endif /* Mandelbrot_h */
