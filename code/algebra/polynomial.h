//
//  polynomial.h
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
 \file polynomial.h
 \brief A data structure and a collection of basic functions to deal with polynomials using arbitrary precision complex numbers.
 
 The complex numbers used are define in @c mpc.h as @c mpc_struct based on [mpfr] (https://www.mpfr.org).
 
 The functions are not fully optimized, they are used only for verification purposes in our project.
*/

#ifndef polynomial_h
#define polynomial_h

#include <stdio.h>

#include "mpc.h"
#include "mpv.h"

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Structures
// //////////////////////////////////////////////////////////////////////////////////////////

/// @struct poly_struct
/// @brief Polynomial structure
///
/// The structure contains the following data:
///
/// @c deg : degree,
///
/// @c prec : precision,
///
/// @c a : coefficient list, initialized to precision @c prec  by @c poly_init.
///
/// @c x : list of powers of the evaluation point, computed to precision @c prec
///
/// @c s,t,u,n : buffers, initialized to precision @c prec
///
/// Coefficients are given in the natural order (trailing degrees first):
/// @c a[j] is the coefficient of @z z^j.
///
/// @warning Optimizations for real coeff or arguments are sorted out by MPFR
typedef struct {
    int deg, prec;
    
    mpc_struct *a, *x;
    mpc s, t, u, n;
} poly_struct;

/// \brief Practical wrapper for @c poly_struct
///
/// To avoid the constant use @c * and @c & the type  @c poly is a pointer.
typedef poly_struct poly[1];

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Constructors
// //////////////////////////////////////////////////////////////////////////////////////////

/// @brief Initializes the polynomial @c p of given degree @c deg and precision @c prec, with the coefficients @c a.
///
/// @param p the polymial
/// @param deg the degree
/// @param a the complex coefficients (array of @c mpc)
/// @param prec the precision of coefficients in bits
///
/// @return 1 if successfull, 0 otherwise
int poly_init(poly p, int deg, mpc_struct *a, int prec);

/// @brief Initializes the polynomial @c p of given degree @c deg and precision @c prec, with the coefficients @c v.
///
/// @param p the polymial
/// @param deg the degree
/// @param v the complex coefficients (@c mpVector)
/// @param prec the precision of coefficients in bits
///
/// @return 1 if successfull, 0 otherwise
int poly_init_vect(poly p, int deg, mpv_t v, int prec);

/// @brief Initializes the polynomial @c p of given degree @c deg and precision @c prec, with the coefficients @c a.
///
/// @param p the polymial
/// @param deg the degree
/// @param a the real coefficients (array of @c mpfr)
/// @param prec the precision of coefficients in bits
///
/// @return 1 if successfull, 0 otherwise
int poly_init_real(poly p, int deg, __mpfr_struct *a, int prec);

/// @brief Construct the complex polynomial @c p defined implicitely by the values @c pk
/// of the sum of the first (posivite) powers of its roots, starting with power 1.
///
/// @param p the polynomial
/// @param deg the degree
/// @param pk the sums of the (positive) powers of the roots, starting with power 1 (array of @c mpc)
/// @param prec the precision of coefficients in bits
///
/// @return 1 if successfull, 0 otherwise
int poly_init_root_powers(poly p, int deg, mpc_struct *pk, int prec);

/// @brief Construct the real polynomial @c p defined implicitely by the values @c pk
/// of the sum of the first (posivite) powers of its roots, starting with power 1.
///
/// @param p the polynomial
/// @param deg the degree
/// @param pk the sums of the (positive) powers of the roots, starting with power 1 (array of @c mpfr)
/// @param prec the precision of coefficients in bits
/// @param eps a bound below which, in absolute value, the coefficients are replaced with @c 0 and the degree recomputed
///
/// @return @ref true if successfull, @ref otherwise otherwise
bool poly_init_root_powers_real(poly p, int deg, __mpfr_struct *pk, int prec, ldbl eps);

/// @brief Clear the polynomial @c p
///
/// @param p the polynomial
///
/// @return 1 if successfull, 0 otherwise
///
/// Clear each @c mpc member of @c p and free the memory.
/// @c p can either be re-initialized with @c poly_init or safely discarded.
int poly_clear(poly p);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Evaluations
// //////////////////////////////////////////////////////////////////////////////////////////

/// @brief Evaluate the polynomial @c p at a point @c z.
///
/// @param p the polynomial
/// @param v an @c mpc number where the value will be written
/// @param z the evaluation point (@c mpc number)
///
/// @return 1 if successfull, 0 otherwise
///
/// The taget @c v should already have been initialized with @c mpc_init.
///
/// The evaluation computes the powers of @c z and stores them in the @c x member
/// of the polynomial @c p.
///
/// Overall cost is O(@c 2*deg) multiplications.
///
/// @warning buffers @c p->s and @c p->t are modified.
int poly_val(poly p, mpc v, mpc z);

/// @brief Evaluate the derivate of the polynomial @c p at a point @c z.
///
/// @param p the polynomial
/// @param d an @c mpc number where the value of the derivative will be written
/// @param z the evaluation point (@c mpc number)
///
/// @return 1 if successfull, 0 otherwise
///
/// The taget @c v should already have been initialized with @c mpc_init.
///
/// The evaluation computes the powers of @c z and stores them in the @c x member
/// of the polynomial @c p.
///
/// Overall cost is O(@c 3*deg) multiplications including a @c long  scaling.
///
/// @warning buffers @c p->u and @c p->t are modified.
int poly_der(poly p, mpc d, mpc z);

/// @brief Evaluate the derivate of the polynomial @c p at a point @c z.
///
/// @param p the polynomial
/// @param v an @c mpc number where the value will be written
/// @param d an @c mpc number where the value of the derivative will be written
/// @param z the evaluation point (@c mpc number)
///
/// @return 1 if successfull, 0 otherwise
///
/// The taget @c v should already have been initialized with @c mpc_init.
///
/// The evaluation computes the powers of @c z and stores them in the @c x member
/// of the polynomial @c p.
///
/// Overall cost is O(@c 4*deg) multiplications  including a @c long  scaling.
///
/// @warning buffers @c p->s, @c p->u and @c p->t are modified.
int poly_val_der(poly p, mpc v, mpc d, mpc z);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Root finding and level sets using Newton's method
// //////////////////////////////////////////////////////////////////////////////////////////

/// @brief Compute the newton term of @c p at the point @c z.
///
/// @param p the polynomial
/// @param nt an @c mpc number where the newton term will be written
/// @param z the evaluation point (@c mpc number)
///
/// @return 1 if successfull, 0 otherwise
///
/// The newton term is
/// \f[
/// nt = \frac{ p(z) }{ p'(z) }
/// \f]
///
/// @warning buffers @c p->s, @c p->u and @c p->t are modified.
int poly_newton_term(poly p, mpc nt, mpc z);

/// @brief Compute a newton term of @c p to solve p(z) = t.
///
/// @param p the polynomial
/// @param nt an @c mpc number where the newton term will be written
/// @param z the evaluation point (@c mpc number)
/// @param t the target value
///
/// @return 1 if successfull, 0 otherwise
///
/// To solve \f$ p(z) = t,\f$ the newton term is
/// \f[
/// nt = \frac{ p(z) - t }{ p'(z) }
/// \f]
///
/// @warning buffers @c p->s, @c p->u and @c p->t are modified.
int poly_newton_term_sol(poly p, mpc nt, mpc z, mpc t);

/// @brief Attempts to find a roots of @c p.
///
/// Attempts to find a roots of @c p using Newton method starting from @c z.
///
/// Convergence is reached when the modulus of the value of @c p
/// at the current point drops below @c err.
///
/// Divergence is claimed if values larger than @c 1E5 (in modulus) are reached.
///
/// Once convergence has been reached, an additional @c refine Newton steps
/// are performed.
/// 
/// @param p the polynomial
/// @param r an @c mpc number where the root will be written
/// @param z the starting point (@c mpc number)
/// @param iter the maximal number of Newton iterations
/// @param err the convergence error
/// @param refine the refinement steps
///
/// @return @ref true if successfull, @ref otherwise otherwise
bool poly_root(poly p, mpc r, mpc z, int iter, long double err, int refine);

mpv poly_roots(poly p, ldbl r, ldbl eps);

/// @brief Attempts to find a roots of @c p-t.
///
/// @param p the polynomial
/// @param r an @c mpc number where the root will be written
/// @param z the starting point (@c mpc number)
/// @param t the target value
/// @param iter the maximal number of Newton iterations
/// @param err the convergence error
/// @param refine the refinement steps
///
/// @return 1 if successfull, 0 otherwise
///
/// @c poly_sol attempts to find a solution of \f$ p(z) = t\f$ using Newton method starting from @c z.
///
/// Convergence is reached when the modulus of the value of @c p
/// at the current point drops below @c err.
///
/// Divergence is claimed if values larger than @c 1E5 (in modulus) are reached.
///
/// Once convergence has been reached, an additional @c refine Newton steps
/// are performed.
int poly_sol(poly p, mpc r, mpc z, mpc t, int iter, long double err, int refine);

#endif /* polynomial_h */
