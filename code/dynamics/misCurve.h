//
//  misCurve.h
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
 \file misCurve.h
 \brief A data structure and a collection of basic functions to work with level curves of the polynomials \f$ p_{n,k} := p_{n+k} - p_n \f$.
 
 Every opportunity to optimize speed is used, including
 computations with fp80 until its precision limit LEVM_MAX_FP80_STEP is attained. The switch to precision @c prec
 (see the definition of misCurve) is performed automatically. No more than LEVM_MAX_ITER  Newton steps are performed
 for computing one new point. If the convergence is slower, the operation fails and results in an invalid curve. The errors
 can be controlled by the user, an appropriate @c guard can be chosen (depending on the @c prePeriod, @c period,
 @c prec and even @c radius) to force smaller steps in the search for new points of the curve. The larger @c guard,
 the higher the probability that the curve is computable, with the penalty of execution time.
 
 For explicit formulas of the points of a level curve, see the definition of the data structure @ref misCurve.
 
 The real and complex  numbers  with arbitrary precision are based on [mpfr] (https://www.mpfr.org).
*/

#ifndef misCurve_h
#define misCurve_h

#include <mpfr.h>

#include "fp80.h"
#include "mpc.h"
#include "mpv.h"
#include "memFile.h"
#include "levCurve.h"
#include "mandel.h"

/// The @c levm file ID.
#define LEVM_FILE_ID     "levm v01"

/// The min length of the file header.
#define LEVM_MIN_HEADER_LEN   ((int)(HEADER_MIN_LEN + 40 + mpv_mpfrs(120)))

/// The min precision of numbers. With this exact precision, computations will be performed with \ref fp80.
#define LEVM_MIN_PREC      64

/// The max period when computations are performed with \ref fp80.
#define LEVM_MAX_FP80_PREC 25

/// Max iterates for non-certified curves with low precision.
#define LEVM_MAX_ITER      40

/// Max period for the level curves.
#define LEVM_MAX_PERIOD    50

/// Max power of @c 2 in misCurve.
#define LEVM_MAX_2POW      62

/// The precision of @c eps in misCurve.
#define LEVM_PREC_EPS     120

/// The min @c guard in misCurve.
#define LEVM_MIN_GUARD      6.1

/// The min modulus of a valid fp80 Newton term.
#define LEVM_MIN_FP80_STEP  1E-13

/// @struct misCurve
/// \brief An arc of a level curve of \f$ p_n \f$ (or equivalently
/// of \f$ \Psi_n \f$), computed with arbitrary precision floating point numbers.
///
/// The points on the level curve are such that
/// \f[ p_{\mathrm{period}}(\mathrm{points}[k]) = \mathrm{radius} \cdot \exp\left(
/// 2^{\mathrm{period - ang2pow}}\pi i (k + \mathrm{startAngle}) \right),\f]
/// or
/// \f[ \mathrm{points}[k] = \Psi_{\mathrm{period}}\left(\mathrm{radius}^{2^{1 - \mathrm{period}}} \cdot \exp\left(
/// \frac{2\pi i (k + \mathrm{startAngle})}{2^{\mathrm{ang2pow}}} \right)  \right) ,\f]
/// where \f$ 0 \leq k \leq \mathrm{endAngle - startAngle} \f$.
typedef struct {
    uint prePeriod;    ///< the pre-pretion @c n of the polynomial @c p_{n,k}=p_{n+k}-p_k
    uint period;       ///< the period @c k of the polynomial @c p_{n,k}=p_{n+k}-p_k
    uint ang2pow;      ///< the power of @c 2 by which to divide the integer angles like @c stAng and @c endAng, at least @c 1
    uint prec;         ///< the precision (in bits) of the stored coordinates, min \ref LEVM_MIN_PREC
    ulong startAngle;  ///< the lowest angle (or argument) of a point in the curve
    ulong endAngle;    ///< the highest angle (or argument) of a point in the curve
    mpfr_t radius;     ///< the level of the curve, or \f$ |p_{n,k}(z)| \f$ for every @c z on the curve
    double guard;      ///< all searches from \f$ z_0 \f$ to \f$ z_1 \f$ by the Newton method satisfy \f$ |p_{n,k}(z_0)| - 4 \geq \mathrm{guard} \cdot |p_{n,k}(z_0)-p_{n,k}(z_1)| \f$
    mpv_t points;      ///< the vector of points
} misCurve;

/// Convenince type for a pointer to misCurve.
typedef misCurve *levm;

/// The number of segments in the misCurve.
#define levm_segs(lm) ((lm)->endAngle - (lm)->startAngle)

/// The number of points in the misCurve.
#define levm_count(lm) (levm_segs(lm) + 1)

/// The length of the header of the fully initialized misCurve, given by a pointer.
#define levm_headerLen(lm) ((int)(HEADER_MIN_LEN + 40 + mpv_mpfrSize((lm)->radius)))

/// The length of the file of the fully initialized misCurve, given by a pointer.
#define levm_fileLen(lm) (levm_headerLen(lm) + mpv_pointsLen((lm)->points))

// MARK: Creation, computation and destruction

/// @brief  Computes the @c misCurve of given pre-period @c pp and period @c per, with the given precision
/// @c prec and @c 2^{tpow-1}+1 points, representing the angles from @c 0 to @c PI.
///
/// The number of intermediary points computed is a
/// \f[ 2^p > \frac{2\pi r \cdot \mathrm{guard}}{r-4}.\f]
/// It is recommended to use \f$ r \geq 10 \f$.
///
/// @param per the pre-period of the curve
/// @param per the period of the curve
/// @param tpow the power of two by which to divide the unit circle
/// @param prec the precision of the points in bits
/// @param r the level of the curve, larger than @c 4
/// @param guard the guard, see the definition of misCurve
///
/// @return the new curve if it could be computed, @c NULL if some error occurred
levm levm_new(uint pp, uint per, uint tpow, uint prec, mpfr_t r, double guard);

/// @brief  Computes the @c misCurve of the simplified polynomial of given pre-period @c pp and period @c per,
/// with the given precision @c prec and @c 2^{tpow-1}+1 points, representing the angles from @c 0 to @c PI.
///
/// The number of intermediary points computed is a
/// \f[ 2^p > \frac{2\pi r \cdot \mathrm{guard}}{r-4}.\f]
/// It is recommended to use \f$ r \geq 10 \f$.
///
/// @param per the pre-period of the curve
/// @param per the period of the curve
/// @param tpow the power of two by which to divide the unit circle
/// @param prec the precision of the points in bits
/// @param r the level of the curve, larger than @c 4
/// @param guard the guard, see the definition of misCurve
///
/// @return the new curve if it could be computed, @c NULL if some error occurred
levm levm_simple_new(uint pp, uint per, uint tpow, uint prec, mpfr_t r, double guard);

/// @brief Frees all the memory used by the misCurve @c lm. It does not fail if @c lm==NULL.
///
/// @param lm the level curve to dispose of
void levm_free(levm lm);

// MARK: Pointwise operations

/// @brief Computes the @c angle (divided by @c 2*PI) of the point with index @c ind in the curve @c lm.
///
/// @param angle the result
/// @param lm the curve
/// @param ind the index
///
/// @return @ref true if the @c angle holds the desired value, @ref false otherwise
bool levm_angle(mpfr_t angle, levm lm, ulong ind);

/// @brief Computes with low precision the angle of the point with index @c ind in the curve @c lm.
///
/// @param lm the curve
/// @param ind the index
///
/// @return the angle (divided by @c 2*PI), @c -1 if some error occurred
ldbl levm_anglel(levm lm, ulong ind);

/// @brief Returns the point with index @c ind in the curve @c lm.
///
/// @param p the point
/// @param lm the level curve
/// @param ind the index of the point in @c lm
///
/// @return @ref true if the @c p holds the desired value, @ref false otherwise
bool levm_point(mpc p, levm lm, ulong ind);

/// @brief Returns a low precision value of the point with index @c ind in the curve @c lm.
///
/// @param p the point
/// @param lm the level curve
/// @param ind the index of the point in @c lm
///
/// @return @ref true if the @c p holds the desired value, @ref false otherwise
bool levm_pointl(fp80 p, levm lm, ulong ind);

// MARK: Refinements and derivations

/// @brief Returns a sub-curve of @c lm.
///
/// The first point that is considered is at the position @c start in the curve @c lm, the second at position
/// @c start+step and so on, until @c count values have been read. If @c count is too large,
/// the function fails and returns @c NULL.
///
/// \warning \c step must be a power of @c 2 so that the new misCurve has a valid dyadic structure; @c start
/// must also be a multiple of @c step
///
/// @param lm the level curve to cut from
/// @param start the start position in @c lm
/// @param step the step in @c lm
/// @param count the number of points of the new curve
///
/// @return the sub-curve of @c lm, @c NULL if the parameters are not valid
levm levm_sub_curve(levm lm, ulong start, ulong step, ulong count);

/// @brief Computes a curve with the same @c prePeriod,@c period and @c radius as @c lm, with different structure.
/// If @c tpow>lm->ang2pow then more points are computed. Conversely, if @c tpow<lm->ang2pow
/// a sub-curve will be returned. The same logic applies to @c prec.
/// The computation is performed with the same @c guard as for the original curve.
///
/// @param lm the starting level curve
/// @param tpow the power of @c 2 of the new curve
/// @param prec the precision (in bits) of the new curve
///
/// @return the new curve if it could be computed, @c NULL if some error occurred
levm levm_refine(levm lm, uint tpow, uint prec);

/// @brief Computes a curve with the same @c period and structure as @c lm, with a different @c radius.
/// The computation is performed with the same @c guard as for the original curve.
///
/// @see @c radius and @c guard in misCurve
///
/// @param lm the original curve
/// @param r the new @c radius
///
/// @return the new curve if it could be computed, @c NULL if some error occurred
levm levm_level(levm lm, mpfr_t r);

/// @brief Computes a curve with the same structure and constants as the hyperbolic curve @c lc, by setting
/// the pre-period to @c pp. If the @c period of @c lc is not equal to @c pp+per, the computation fails and
/// @c NULL is returned.
///
/// @param lm the original hyperbolic curve
/// @param pp the pre-period of the returned curve
/// @param the period of the returned curve
/// @param the guard of the returned curve
///
/// @return the new curve if it could be computed, @c NULL if some error occurred
levm levm_from_hyp(levc lc, uint pp, uint per, double guard);

/// @brief Computes a simple curve with the same structure and constants as the hyperbolic curve @c lc, by setting
/// the pre-period to @c pp. If the @c period of @c lc is not equal to @c pp+per, the computation fails and
/// @c NULL is returned.
///
/// Simple curves correspond to the simplified Misiurewicz polynomial @c p_{pp,per}/p_{pp-1,per}.
///
/// @param lm the original hyperbolic curve
/// @param pp the pre-period of the returned curve
/// @param the period of the returned curve
/// @param the guard of the returned curve
///
/// @return the new curve if it could be computed, @c NULL if some error occurred
levm levm_simple_from_hyp(levc lc, uint pp, uint per, double guard);

// MARK: IO functions

/// @brief Writes the level curve @c lm to the file given by its path  @c fileName in CSV format, that is in human readable text.
///
/// The parameters of the curve will appear in the first few lines, one per line of text.
///
/// @param lm the level curve
/// @param fileName the path to the file
/// @param digits the number of digits after the decimal point to write, for the points of the curve
///
/// @return @ref true if the operation completed successfully, @ref false otherwise
bool levm_write_csv(levm lm, char *fileName, uint digits);

/// @brief Reads the content of the CSV file given by the path @c fileName into a new misCurve.
///
/// @param fileName the path of the file
///
/// @return the new level curve with the content of the file, @c NULL if the file is not valid
levm levm_read_csv(char *fileName);

/// @brief Writes the level curve @c lm to the file @c fileName.
///
/// @param lm the level curve to write
/// @param fileName the path of the file
///
/// @return @ref true if the operation completed successfully, @ref false otherwise
bool levm_write_to(levm lm, char *fileName);

/// @brief Reads the contentof the file @c f into a new misCurve.
///
/// @param fileName the path of the file
///
/// @return the new level curve with the content of the file, @c NULL if the file is not valid
levm levm_read_from(char *fileName);

/// @brief Writes the level curve @c lm to the file @c fileName.
///
/// @param lm the level curve to write
/// @param f the file
/// @param pos the postion in the file or @c -1 to use the current position
///
/// @return the number of bytes written, @c 0 if some error occurred
ulong levm_write(levm lm, FILE *f, long pos);

/// @brief Reads the content starting at position @c pos in the file @c f into a new misCurve.
///
/// @param f the file
/// @param pos the absolute postion in the file or @c -1 to use the current position
///
/// @return the new level curve with the content of the file, @c NULL if the file is not valid
levm levm_read(FILE *f, long pos);

/// @brief Partially reads the content starting at position @c pos in the file @c f into a new misCurve.
///
/// The first point that is read is at the position @c start in the curve stored in the file, the second at position
/// @c start+step and so on, until @c count values have been read. If any of those values is not accessible,
/// the function fails and returns @c NULL.
///
/// \warning \c step must be a power of @c 2 so that the new misCurve has a valid dyadic structure; @c start
/// must also be a multiple of @c step
///
/// @param f the file
/// @param pos the absolute postion in the file or @c -1 to use the current position
/// @param start the start position in the misCurve stored in the file
/// @param step the step in stored misCurve
/// @param count the number of values to read from the file
///
/// @return the new vector with the partial content of the file, @c NULL if the file or the parameters are not valid
levm levm_read_partial(FILE *f, long pos, ulong start, ulong step, ulong count);


#endif /* misCurve_h */
