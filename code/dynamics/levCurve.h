//
//  levCurve.h
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
  \file levCurve.h
  \brief A data structure and a collection of basic functions to work with level curves of the polynomials \f$ p_n \f$.
  
  There are two types of level curves represented by this type: @b proven and @b approximated, see the @c proven field of levelCurve.
  
  When working with @b approximated level curves (see levc_new()), every opportunity to optimize speed is used, including
  computations with fp80 until its precision limit LEVC_MAX_FP80_STEP is attained. The switch to precision @c prec
  (see the definition of levelCurve) is performed automatically. No more than LEVC_MAX_ITER   Newton steps are performed
  for computing one new point. If the convergence is slower, the operation fails and results in an invalid curve. The errors
  can be controlled by the user, an appropriate @c guard can be chosen (depending on the @c period, @c prec and
  even @c radius) to force smaller steps in the search for new points of the curve. The larger @c guard, the higher
  the probability that the curve is computable, with the penalty of execution time.
  
  When working with @b proven curves (see levc_prove()), each step of the computation is governed by the class psi.h,
  with respect to the precision @c prec of the levelCurve and the bound @c eps for the error. If @c eps is too small
  compared to @c prec, the curve may not compute.
  
  The advanced operations levc_subCurve(), levc_refine(), levc_level() and levc_nextPeriod() inherit the type (@b proven or
  @b approximated) and other charcteristics of the original curve, so the results are in the same class.
  
  For explicit formulas of the points of a level curve, see the definition of the data structure levelCurve. Note that for
  \f$ \mathrm{period}=1 \f$, both \f$ p_1 \f$ and \f$ \Psi_1 \f$ are the identity and the points of the curve are
  equal to \f$ \mathrm{radius} \cdot u \f$,  where \f$ u \f$ is some root of unity of order \f$ 2^{ \mathrm{ang2pow}} \f$.
  
  The real and complex  numbers  with arbitrary precision are based on [mpfr] (https://www.mpfr.org).
 */

#ifndef levCurve_h
#define levCurve_h

#include <mpfr.h>

#include "fp80.h"
#include "mpc.h"
#include "mpv.h"

// MARK: Constants

/// The @c levc file ID.
#define LEVC_FILE_ID     "levc v01"

/// The min length of the file header.
#define LEVC_MIN_HEADER_LEN  113

/// The min precision of numbers. With this exact precision, computations will be performed with \ref fp80.
#define LEVC_MIN_PREC         64

/// The min precision of numbers for proven curves.
#define LEVC_MIN_PROVEN_PREC 120

/// The max period when computations are performed with \ref fp80.
#define LEVC_MAX_FP80_PREC    30

/// Max iterates for non-certified curves with low precision.
#define LEVC_MAX_ITER         40

/// Max period for the level curves.
#define LEVC_MAX_PERIOD       50

/// Max power of @c 2 in levelCurve.
#define LEVC_MAX_2POW         62

/// The precision of @c eps in levelCurve.
#define LEVC_PREC_EPS        120

/// The precision of numbers in computing intermediary steps.
#define LEVC_PREC_INTERN     120

/// The min @c guard in levelCurve.
#define LEVC_MIN_GUARD         6.1

/// The min modulus of a valid fp80 Newton term.
#define LEVC_MIN_FP80_STEP     1E-15

/// The number of extra bits of precision compared to the scale of @c eps.
#define LEVC_EXTRA_PREC        5

// MARK: Data types

/// @struct levelCurve
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
    uint period;       ///< the index @c n of @c psi_n or equivalently @c p_n used to compute the level curve, at least @c 1
    uint ang2pow;      ///< the power of @c 2 by which to divide the integer angles like @c startAngle and @c endAngle, at least @c 1
    uint prec;         ///< the precision (in bits) of the stored coordinates, min \ref LEVC_MIN_PREC
    ulong startAngle;  ///< the lowest angle (or argument) of a point in the curve
    ulong endAngle;    ///< the highest angle (or argument) of a point in the curve
    mpfr_t radius;     ///< the level of the curve, or \f$ |p_n(z)| \f$ for every @c z on the curve
    mpfr_t eps;        ///< max error of the coordinates, certified only if @c proven
    double guard;      ///< if not @c proven, all searches from \f$ z_0 \f$ to \f$ z_1 \f$ by the Newton method satisfy \f$ |p_n(z_0)| - 2 \geq \mathrm{guard} \cdot |p_n(z_0)-p_n(z_1)| \f$
    mpv_t points;      ///< the vector of points
} levelCurve;

/// Convenince type for a pointer to @ref levelCurve.
typedef levelCurve *levc;

// MARK: Macros

/// The number of segments in the levelCurve.
#define levc_segs(lc) ((lc)->endAngle - (lc)->startAngle)

/// The number of points in the levelCurve.
#define levc_count(lc) (levc_segs(lc) + 1)

/// The length of the header of the fully initialized levelCurve, given by a pointer.
#define levc_headerLen(lc) ((int)(73 + mpv_mpfrSize((lc)->radius) + mpv_mpfrSize((lc)->eps)))

/// The length of the file of the fully initialized levelCurve, given by a pointer.
#define levc_fileLen(lc) (levc_headerLen(lc) + mpv_pointsLen((lc)->points))

/// The length of the file of the fully initialized levelCurve, given by parameters.
#define levc_file_size(r_prec, eps_prec, p_count, p_prec) (73 + mpv_mpfrs(r_prec) + mpv_mpfrs(eps_prec) + \
                                                           (long)(mpv_point_len(p_prec)) * 2 * (p_count))

// MARK: Creation, computation and destruction

/// @brief  Computes the levelCurve of given period @c per, with the given precision @c prec and @c 2^{tpow-1}+1
/// points, representing the angles from @c 0 to @c PI.
///
/// The type of the curve is @b approximative. As an exception, if @c per==1, @c 2^{tpow} points are computed
/// as a @b proven level set, as they are simply the points equidistributed on the circle of radius @c r.
///
/// The number of intermediary points computed is a
/// \f[ 2^p > \frac{2\pi r \cdot \mathrm{guard}}{r-2},\f]
/// so for the minimum \f$ \mathrm{guard} = 2 \f$, it is recommended to use \f$ r \geq 10 \f$ to get the minimum
/// \f$ p = 4 \f$.
/// FIXME: I do not understand the computation here.
///
/// @see levc_prove()
///
/// @param per the period of the curve
/// @param tpow the power of two by which to divide the unit circle
/// @param prec the precision of the points in bits
/// @param r the level of the curve, larger than @c 2
/// @param guard the guard, see the definition of levelCurve
///
/// @return the new curve if it could be computed, @c NULL if some error occurred
levc levc_new(uint per, uint tpow, uint prec, mpfr_t r, double guard);

/// @brief Computes the roots of unity multiplied by @c r as a levelCurve.
///
/// @param tpow the power of two to obtain the smalles dyadic angle
/// @param prec the precision
/// @param r the radius
///
/// @return the circle
levc levc_circle(int tpow, uint prec, mpfr_t r);

/// @brief Computes the levelCurve of given period @c per, with the given precision @c prec and @c 2^{tpow-1}+1
/// points, representing the angles from @c 0 to @c PI.
///
/// It performs a proof of both the continuity (the correct
/// choice of pre-images of @c p_n) and that the maximum error of each point is at most @c eps. If this cannot
/// be guaranteed, the function fails and return @c NULL.
///
/// The type of the curve is @b proven.
///
/// As an exception, if @c per==1, @c 2^{tpow} points are computed, as they are simply the points equidistributed
/// on the circle of radius @c r.
///
/// @see levc_new()
///
/// @param per the period of the curve
/// @param tpow the power of two by which to divide the unit circle
/// @param prec the precision of the points in bits
/// @param r the level of the curve, larger than @c 2
/// @param eps the guard, see the definition of levelCurve
///
/// @return the new curve if it could be computed, @c NULL if some error occurred
levc levc_prove(uint per, uint tpow, long prec, mpfr_t r, mpfr_t eps);

/// @brief Frees all the memory used by the levelCurve @c lc. It does not fail if @c lc==NULL.
///
/// @param lc the level curve to dispose of
void levc_free(levc lc);

// MARK: Pointwise operations

/// @brief Computes the @c angle (divided by @c 2*PI) of the point with index @c ind in the curve @c lc.
///
/// @param angle the result
/// @param lc the curve
/// @param ind the index
///
/// @return @ref true if the @c angle holds the desired value, @ref false otherwise
bool levc_angle(mpfr_t angle, levc lc, ulong ind);

/// @brief Computes with low precision the angle of the point with index @c ind in the curve @c lc.
///
/// @param lc the curve
/// @param ind the index
///
/// @return the angle (divided by @c 2*PI), @c -1 if some error occurred
ldbl levc_anglel(levc lc, ulong ind);

/// @brief Returns the point with index @c ind in the curve @c lc.
///
/// @param p the point
/// @param lc the level curve
/// @param ind the index of the point in @c lc
///
/// @return @ref true if the @c p holds the desired value, @ref false otherwise
bool levc_point(mpc p, levc lc, ulong ind);

/// @brief Returns a low precision value of the point with index @c ind in the curve @c lc.
///
/// @param p the point
/// @param lc the level curve
/// @param ind the index of the point in @c lc
///
/// @return @ref true if the @c p holds the desired value, @ref false otherwise
bool levc_pointl(fp80 p, levc lc, ulong ind);

// MARK: Refinements and derivations

/// @brief Returns a sub-curve of @c lc.
///
/// The first point that is considered is at the position @c start in the curve @c lc, the second at position
/// @c start+step and so on, until @c count values have been read. If @c count is too large,
/// the function fails and returns @c NULL.
///
/// \warning \c step must be a power of @c 2 so that the new levelCurve has a valid dyadic structure; @c start
/// must also be a multiple of @c step
///
/// @param lc the level curve to cut from
/// @param start the start position in @c lc
/// @param step the step in @c lc
/// @param count the number of points of the new curve
///
/// @return the sub-curve of @c lc, @c NULL if the parameters are not valid
levc levc_sub_curve(levc lc, ulong start, ulong step, ulong count);

/// @brief Computes a curve with the same @c period and @c radius as @c lc, with different structure.
/// If @c tpow>lc->ang2pow then more points are computed. Conversely, if @c tpow<lc->ang2pow
/// a sub-curve will be returned. The same logic applies to @c prec.
/// The computation is performed with the same constraints as for the original curve.
///
/// @param lc the starting level curve
/// @param tpow the power of @c 2 of the new curve
/// @param prec the precision (in bits) of the new curve
///
/// @return the new curve if it could be computed, @c NULL if some error occurred
levc levc_refine(levc lc, uint tpow, uint prec);

/// @brief Computes a curve with the same @c period and structure as @c lc, with a different @c radius.
/// The computation is performed with the same constraints as for the original curve.
///
/// @see @c radius, @c eps, @c proven and @c guard in levelCurve
///
/// @param lc the original curve
/// @param r the new @c radius
///
/// @return the new curve if it could be computed, @c NULL if some error occurred
levc levc_level(levc lc, mpfr_t r);

/// @brief Computes a curve with the same structure and constants as @c lc, by increasing the period by @c 1.
/// The computation is performed with the same constraints as for the original curve.
///
/// @see @c period, @c radius, @c eps, @c proven and @c guard in levelCurve
///
/// @param lc the original curve
/// @param squareLevel @ref true to square the level of the curve, @ref false to compute the new curve to the same level as @c lc
///
/// @return the new curve if it could be computed, @c NULL if some error occurred
levc levc_next_period(levc lc, bool squareLevel);

// MARK: IO functions

/// @brief Writes the level curve @c lc to the file given by its path  @c fileName in CSV format, that is in human readable text.
///
/// The parameters of the curve will appear in the first few lines, one per line of text.
///
/// @param lc the level curve
/// @param fileName the path to the file
/// @param digits the number of digits after the decimal point to write, for the points of the curve
///
/// @return @ref true if the operation completed successfully, @ref false otherwise
bool levc_write_csv(levc lc, char *fileName, uint digits);

/// @brief Reads the content of the CSV file given by the path @c fileName into a new levelCurve.
///
/// @param fileName the path of the file
///
/// @return the new level curve with the content of the file, @c NULL if the file is not valid
levc levc_read_csv(char *fileName);

/// @brief Writes the level curve @c lc to the file @c fileName.
///
/// @param lc the level curve to write
/// @param fileName the path of the file
///
/// @return @ref true if the operation completed successfully, @ref false otherwise
bool levc_write_to(levc lc, char *fileName);

/// @brief Reads the contentof the file @c f into a new levelCurve.
///
/// @param fileName the path of the file
///
/// @return the new level curve with the content of the file, @c NULL if the file is not valid
levc levc_read_from(char *fileName);

/// @brief Writes the level curve @c lc to the file @c fileName.
///
/// @param lc the level curve to write
/// @param f the file
/// @param pos the postion in the file or @c -1 to use the current position
///
/// @return the number of bytes written, @c 0 if some error occurred
ulong levc_write(levc lc, FILE *f, long pos);

/// @brief Reads the content starting at position @c pos in the file @c f into a new levelCurve.
///
/// @param f the file
/// @param pos the absolute postion in the file or @c -1 to use the current position
///
/// @return the new level curve with the content of the file, @c NULL if the file is not valid
levc levc_read(FILE *f, long pos);

/// @brief Partially reads the content starting at position @c pos in the file @c f into a new levelCurve.
///
/// The first point that is read is at the position @c start in the curve stored in the file, the second at position
/// @c start+step and so on, until @c count values have been read. If any of those values is not accessible,
/// the function fails and returns @c NULL.
///
/// \warning \c step must be a power of @c 2 so that the new levelCurve has a valid dyadic structure; @c start
/// must also be a multiple of @c step
///
/// @param f the file
/// @param pos the absolute postion in the file or @c -1 to use the current position
/// @param start the start position in the levelCurve stored in the file
/// @param step the step in stored levelCurve
/// @param count the number of values to read from the file
///
/// @return the new vector with the partial content of the file, @c NULL if the file or the parameters are not valid
levc levc_read_partial(FILE *f, long pos, ulong start, ulong step, ulong count);

#endif /* levCurve_h */
