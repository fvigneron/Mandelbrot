//
//  iterates.h
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2022.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2022.
//
//  Copyright 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the Lesser GNU Public Licence v3.0
//

/**
 \file iterates.h
 \brief Functions to produce bitmaps by iteration.
 
 We are interested in images of the Mandelbrot set \f$ \mathcal M \f$, of the Julia sets \f$ \mathcal {J_c} \f$ for
 parameters  \f$ c \in \mathcal M \f$ and of the Newton methods \f$ N_{p_n} \f$ and \f$ N_{p_{n+k}-p_k} \f$ of the
 polynomials \f$ p_n \f$ and respectively \f$ p_{n+k}-p_k \f$. They correspond to the convergence of the Newton method
 toward centers of hyperbolic components of \f$ \mathcal M \f$ and respectively Misurewicz or pre-preriod points in
 \f$ \partial \mathcal M \f$.
*/

#ifndef iterates_h
#define iterates_h

#include <stdio.h>

#include "bitmap.h"
#include "mpc.h"

// MARK: constants and data types

#define ITER_TYPE_M     1    ///< image representing the Mandelbrot set
#define ITER_TYPE_J     2    ///< image representing the Julia set
#define ITER_TYPE_NH    3    ///< image representing the Julia set of the Newtom method for some @c p_n
#define ITER_TYPE_NM    4    ///< image representing the Julia set of the Newtom method for some @c p_n-p_k
#define ITERT_COUNT     4    ///< the number of iterative bitmap types

#define ITERM_ESCAPE    1    ///< matrix of escape times to radius @c 5
#define ITERM_DIST      2    ///< matrix of distance to \f$ \mathcal M \f$, the Mandelbrot set
#define ITERM_DERZ      4    ///< matrix of derivatives in @c z (phase space, geometric mean)
#define ITERM_DERC      8    ///< matrix of derivatives in the parameter @c c (parameter space, geometric mean)
#define ITERM_ESCPOW   16    ///< matrix of -log(escape times) / log(dist to @c M)
#define ITERM_SMOD     32    ///< matrix of the modulus of @c S=lim(dc/dz)f_c^n
#define ITERM_SARG     64    ///< matrix of the argument of @c S=lim(dc/dz)f_c^n
#define ITERM_COUNT     7    ///< the number of possible image types of \f$ \mathcal M \f$

#define ITERJ_ESCAPE    1    ///< matrix of escape times to radius @c 5
#define ITERJ_DIST      2    ///< matrix of distance to \f$ \mathcal {J_c} \f$, the Julia set of \f$ f_c \f$
#define ITERJ_DER       4    ///< matrix of modulus derivatives, geometric mean
#define ITERJ_COUNT     3    ///< the number of possible image types of \f$ \mathcal {J_c} \f$

#define ITER_DIST_BITS    62   ///< the number of bits after the decimal point in bitmaps that represent a distance
#define ITER_DERZ_BITS    53   ///< the number of bits after the decimal point in bitmaps that represent a @c z derivative
#define ITER_DERC_BITS    53   ///< the number of bits after the decimal point in bitmaps that represent a @c c derivative
#define ITER_POW_BITS     60   ///< the number of bits after the decimal point in bitmaps that represent a power
#define ITER_MOD_BITS     32   ///< the number of bits after the decimal point in bitmaps that represent a modulus
#define ITER_ARG_BITS     53   ///< the number of bits after the decimal point in bitmaps that represent an angle
#define ITER_NEWTON_BITS  32   ///< the number of bits that represent the iterate in Newton iterates images

#define ITER_EPSL     1E-16
#define ITER_EPSL2    (ITER_EPSL * ITER_EPSL)

#define ITER_IOTAL    1E-21
#define ITER_IOTAL2   (ITER_IOTAL * ITER_IOTAL)

#define ITER_EPS      1E-35
#define ITER_EPS2     (ITER_EPS * ITER_EPS)

#define ITER_IOTA     1E-220
#define ITER_IOTA2    (ITER_IOTA * ITER_IOTA)

/// The maximum number of threads that will be used by this module.
#define ITER_MAX_THREADS    64

/// The maximum number of iterates for bitmaps of Newton maps.
#define MAX_ITERN      1000000

/// @struct iterBmap
/// \brief This data structure caracterizes an iterative bitmap request. It is only used by the most
/// general method @ref iter_bitmaps().
///
/// All bitmaps in this module can be obtained by a direct metod, with multi-threading.
///
/// It can represent any type of bitmap produced by this module. Internally, it is used to parallelize the computing tasks.
typedef struct {
    volatile int prec;          ///< the precision, in bits, to be used for computations, <= 64 to use @ref fp80
    volatile int prePer;        ///< the pre-period, used only for the Newton method of @c p_per-p_prePer
    volatile int per;           ///< the period, used only for the Newton method
    volatile ulong maxIter;     ///< the number of iterates to stop at
    
    volatile fp80 c80;          ///< the parameter @c c, used only for Julia sets when @c prec<=64
    volatile mpc c;             ///< the parameter @c c, used only for Julia sets when @c prec>64
    
    volatile dyadic_rect r;      ///< the dyadic rectangle of the bitmaps
    volatile bool usesDelta;    ///< @c true to use delta when @c prec>64
    volatile mpc delta;         ///< a translation of @c r when @c prec>64 and @c usesDelta
    
    volatile int type;          ///< the type of the bitmap, one of @c ITER_TYPE_...
    volatile int subType;       ///< the sub-type[s] of the bitmap, one of @c ITERX_..., where @c X=M,J,NH,NM,S
} iterBmap;

/// A pointer to the iterBmap struct, for convenience.
typedef iterBmap *iter;

// MARK: bitmaps of the Mandelbrot set

/// @brief Computes with fp80 and returns the bitmap of parameters @c c in the dyadic rectangle @c r that escape under iteration.
///
/// Paramaters that do not escape are marked with @c 0, the others with the escape time. The escape time is the
/// smallest @c n such that \f$ | p_n(c) | \geq 5 \f$.
///
/// @param r the dyadic rectangle that defines the range and the resolution of the bitmap
/// @param maxIter the max number of iterates
/// @param threads the number of threads to use for iteration
///
/// @return the bitmap or @c NULL if some error occurred
bmap iterM_bitmap80(drect r, int maxIter, int threads);

/// @brief Computes with mpc and returns the bitmap of parameters @c c in the dyadic rectangle @c r that escape under iteration.
///
/// Paramaters that do not escape are marked with @c 0, the others with the escape time. The escape time is the
/// smallest @c n such that \f$ | p_n(c) | \geq 5 \f$.
///
/// @param r the dyadic rectangle that defines the range and the resolution of the bitmap
/// @param delta a translation of @c r, can be @c NULL
/// @param maxIter the max number of iterates
/// @param prec the precision (in bits) to be used for computation, min value @c 32
/// @param threads the number of threads to use for iteration
///
/// @return the bitmap or @c NULL if some error occurred
bmap iterM_bitmap(drect r, mpc delta, ulong maxIter, int prec, int threads);

/// @brief Computes with fp80 and returns the bitmaps of parameters @c c in the dyadic rectangle @c r that
/// are described by the parameter @c types, a bitwise choice among @c ITERM_...
///
/// Those bitmaps are related to the Mandelbrot set.
///
/// @param r the dyadic rectangle that defines the range and the resolution of the bitmap
/// @param maxIter the max number of iterates
/// @param types bitwise OR of the desired types of bitmaps
/// @param threads the number of threads to use for iteration
///
/// @return the vector of requested bitmaps or @c NULL if some error occurred
bmap *iterM_bitmaps80(drect r, int maxIter, int types, int threads);

/// @brief Computes with mpc and returns the bitmaps of parameters @c c in the dyadic rectangle @c r that
/// are described by the parameter @c types, a bitwise choice among @c ITERM_...
///
/// Those bitmaps are related to the Mandelbrot set.
///
/// @param r the dyadic rectangle that defines the range and the resolution of the bitmap
/// @param delta a translation of @c r, can be @c NULL
/// @param maxIter the max number of iterates
/// @param types bitwise OR of the desired types of bitmaps
/// @param prec the precision (in bits) to be used for computation, min value @c 32
/// @param threads the number of threads to use for iteration
///
/// @return the vector of requested bitmaps or @c NULL if some error occurred
bmap *iterM_bitmaps(drect r, mpc delta, ulong maxIter, int types, int prec, int threads);

// MARK: bitmaps of the Julia set

/// @brief Computes with fp80 and returns the bitmap of @c z in the dyadic rectangle @c r that escape under iteration of
/// \f$ f_c(z)=z^2+c \f$. Paramaters that do not escape are marked with @c 0, the others with the escape time.
///
/// The escape time is the smallest @c n such that \f$ | f_c^n(z) | \geq 5 \f$.
///
/// @param c the parameter fo the map @c f_c
/// @param r the dyadic rectangle that defines the range and the resolution of the bitmap
/// @param maxIter the max number of iterates
/// @param threads the number of threads to use for iteration
///
/// @return the bitmap or @c NULL if some error occurred
bmap iterJ_bitmap80(fp80 c, drect r, int maxIter, int threads);

/// @brief Computes with mpc and returns the bitmap of @c z in the dyadic rectangle @c r that escape under iteration of
/// \f$ f_c(z)=z^2+c \f$. Paramaters that do not escape are marked with @c 0, the others with the escape time.
///
/// The escape time is the smallest @c n such that \f$ | f_c^n(z) | \geq 5 \f$.
///
/// @param c the parameter fo the map @c f_c
/// @param r the dyadic rectangle that defines the range and the resolution of the bitmap
/// @param delta a translation of @c r, can be @c NULL
/// @param maxIter the max number of iterates
/// @param prec the precision (in bits) to be used for computation, min value @c 32
/// @param threads the number of threads to use for iteration
///
/// @return the bitmap or @c NULL if some error occurred
bmap iterJ_bitmap(mpc c, drect r, mpc delta, ulong maxIter, int prec, int threads);

/// @brief Computes with fp80 and returns the bitmaps of the dyadic rectangle @c r that are described by the
/// parameter @c types, a bitwise choice among @c ITERJ_...
///
/// Those bitmaps are related to the Julia set of \f$ f_c(z)=z^2+c \f$.
///
/// @param c the parameter fo the map @c f_c
/// @param r the dyadic rectangle that defines the range and the resolution of the bitmap
/// @param maxIter the max number of iterates
/// @param types bitwise OR of the desired types of bitmaps
/// @param threads the number of threads to use for iteration
///
/// @return the bitmap or @c NULL if some error occurred
bmap *iterJ_bitmaps80(fp80 c, drect r, int maxIter, int types, int threads);

/// @brief Computes with mpc and returns the bitmaps of the dyadic rectangle @c r that are described by the
/// parameter @c types, a bitwise choice among @c ITERJ_...
///
/// Those bitmaps are related to the Julia set of \f$ f_c(z)=z^2+c \f$.
///
/// @param c the parameter fo the map @c f_c
/// @param r the dyadic rectangle that defines the range and the resolution of the bitmap
/// @param delta a translation of @c r, can be @c NULL
/// @param maxIter the max number of iterates
/// @param types bitwise OR of the desired types of bitmaps
/// @param prec the precision (in bits) to be used for computation, min value @c 32
/// @param threads the number of threads to use for iteration
///
/// @return the bitmap or @c NULL if some error occurred
bmap *iterJ_bitmaps(mpc c, drect r, mpc delta, ulong maxIter, int types, int prec, int threads);

// MARK: bitmaps of the Newton methods (Julia set and Fatou set components)

/// @brief Computes with fp80 and returns the bitmap of @c z in the dyadic rectangle @c r that converge under iteration of
/// \f$ N_{p_n}(z)=z-\frac{p_n(z)}{p'_n(z)} \f$ to a fixed point of \f$ N_{p_n} \f$ (a root of \f$ p_n \f$, the center
/// of a hyperbolic component of \f$ \mathcal M \f$ of period \f$ n \f$).
///
/// Paramaters that do not converge are marked with @c 0, the others with the first entry time to a auto-detected disk that is
/// sent into itself by \f$ N_{p_n}(z)=z-\frac{p_n(z)}{p'_n(z)} \f$.
///
/// @param n the index of @c p_n
/// @param r the dyadic rectangle that defines the range and the resolution of the bitmap
/// @param maxIter the max number of iterates, at most @c MAX_ITERN
/// @param quick @ref true for maximal optimization, @ref false for better quality
///
/// @return the bitmap or @c NULL if some error occurred
bmap iterN_hyp_bitmap80(int n, drect r, int maxIter, bool quick);

/// @brief Computes with mpc and returns the bitmap of @c z in the dyadic rectangle @c r that converge under iteration of
/// \f$ N_{p_n}(z)=z-\frac{p_n(z)}{p'_n(z)} \f$ to a fixed point of \f$ N_{p_c} \f$ (a root of \f$ p_n \f$, the center
/// of a hyperbolic component of \f$ \mathcal M \f$ of period \f$ n \f$).
///
/// Paramaters that do not converge are marked with @c 0, the others with the first entry time to a auto-detected disk
/// that is sent into itself by \f$ N_{p_n}(z)=z-\frac{p_n(z)}{p'_n(z)} \f$.
///
/// @param n the index of @c p_n
/// @param r the dyadic rectangle that defines the range and the resolution of the bitmap
/// @param delta a translation of @c r, can be @c NULL
/// @param maxIter the max number of iterates, at most @c MAX_ITERN
/// @param prec the precision (in bits) to be used for computation, min value @c 32
///
/// @return the bitmap or @c NULL if some error occurred
bmap iterN_hyp_bitmap(int n, drect r, mpc delta, int maxIter, int prec);

/// @brief Computes with fp80 and returns the bitmap of @c z in the dyadic rectangle @c r that converge under iteration of
/// \f$ N_{p_{n+k} -p_k}(z)=z-\frac{p_{n+k}(z) - p_k(z)}{p'_{n+k}(z)-p'_k(c)} \f$ to a fixed point of \f$ N_{p_{n+k}-p_k} \f$
/// (a root of \f$ p_{n+k} - p_k \f$, a Misiurewicz point of period \f$ n \f$ and pre-periof \f$ k \f$).
///
/// Paramaters that do not converge are marked with @c 0, the others with the first entry time to a auto-detected
/// disk that is sent into itself by \f$ N_{p_{n+k} -p_k}(z)=z-\frac{p_{n+k}(z) - p_k(z)}{p'_{n+k}(z)-p'_k(c)} \f$.
///
/// @param k the pre-period of the Misiurewicz points
/// @param n the period of the Misiurewicz points
/// @param r the dyadic rectangle that defines the range and the resolution of the bitmap
/// @param maxIter the max number of iterates, at most @c MAX_ITERN
/// @param quick @ref true for maximal optimization, @ref false for better quality 
///
/// @return the bitmap or @c NULL if some error occurred
bmap iterN_mis_bitmap80(int k, int n, drect r, int maxIter, bool quick);

/// @brief Computes with mpc and returns the bitmap of @c z in the dyadic rectangle @c r that converge under iteration of
/// \f$ N_{p_{n+k} -p_k}(z)=z-\frac{p_{n+k}(z) - p_k(z)}{p'_{n+k}(z)-p'_k(c)} \f$ to a fixed point of \f$ N_{p_{n+k}-p_k} \f$
/// (a root of \f$ p_{n+k} - p_k \f$, a Misiurewicz point of period \f$ n \f$ and pre-periof \f$ k \f$).
///
/// Paramaters that do not converge are marked with @c 0, the others with the first entry time to a auto-detected
/// disk that is sent into itself by \f$ N_{p_{n+k} -p_k}(z)=z-\frac{p_{n+k}(z) - p_k(z)}{p'_{n+k}(z)-p'_k(c)} \f$.
///
/// @param k the pre-period of the Misiurewicz points
/// @param n the period of the Misiurewicz points
/// @param r the dyadic rectangle that defines the range and the resolution of the bitmap
/// @param delta a translation of @c r, can be @c NULL
/// @param maxIter the max number of iterates, at most @c MAX_ITERN
/// @param prec the precision (in bits) to be used for computation, min value @c 32
///
/// @return the bitmap or @c NULL if some error occurred
bmap iterN_mis_bitmap(int k, int n, drect r, mpc delta, int maxIter, int prec);

// MARK: highest level of abstraction of this module

/// @brief Computes and returns the bitmaps requested by @c it, using @c threads threads for iterations.
///
/// The main thread will distribute the computing tasks (sub-bitmaps) and will reassemble the bitmaps.
///
/// @warning it should only be called by the main thread, when no other multi-threading is executing.
///
/// @param it the iterative bitmaps request
/// @param threads the number of threads to use for iteration
///
/// @return the vector of requested bitmaps or @c NULL if some error occurred
bmap *iter_bitmaps(iter it, int threads);

#endif /* iterates_h */
