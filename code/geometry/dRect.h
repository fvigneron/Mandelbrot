//
//  dRect.h
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
 \file dRect.h
 \brief A dyadic rectangle data structure and basic operations.
*/

#ifndef dRect_h
#define dRect_h

#include <stdio.h>

#include "fp80.h"
#include "mpc.h"

/// @struct dyadic_rect
/// \brief This struct caraterizes and stores a dyadic rectangles.
///
/// All coordinates are to be divided by @c 2^{tpow} to obtain the absolute coordinates.
typedef struct {
    long x;   ///< left coordinate
    long y;   ///< bottom coordinate
    uint w;   ///< width
    uint h;   ///< height
    int tpow; ///< the power of two by which to divide all the coordinates
} dyadic_rect;

/// Convenience type for @c dyadic_rect.
typedef dyadic_rect drect_t[1];

/// Convenience type for a pointer to @c dyadic_rect.
typedef dyadic_rect *drect;

/// @brief Initialize a new rectangle
///
/// @param r the rectangle to be initialized
/// @param x the @c x coordinate of the bottom left corner
/// @param y the @c y coordinate of the bottom left corner
/// @param w the width of the rectangle
/// @param h the height of the rectangle
/// @param tpow the power of two by which to divide all the coordinates
void drect_init(drect r, long x, long y, uint w, uint h, int tpow);

// //////////////////////////////////////////////////
/// MARK: Coordinates
// //////////////////////////////////////////////////

/// @brief Relative to absolute coordinates
///
/// Computes the complex number that represents the point
/// at relative coordinates @c x and @c y from the left bottom point
/// of the dyadic rectangle.
///
/// Relative coordinates are  multiple of \f$ 2^{tpow} \f$.
///
/// @param c the result
/// @param r the rectangle
/// @param x the @c x coordinate, relative to the rectangle r
/// @param y the @c y coordinate, relative to the rectangle r
///
/// @return @ref true if successfull, @ref false otherwise
bool drect_rel_to_abs_coords80(fp80 c, drect r, long x, long y);

/// @brief Absolute to relative coordinates
///
/// @param x the @c x absolute coordinate
/// @param r the rectangle
///
/// @return the relative @c x coordinate with respect to  @c r
long drect_abs_to_rel_x80(drect r, ldbl x);

/// @brief Absolute to relative coordinates
///
/// @param y the @c y absolute coordinate
/// @param r the rectangle
///
/// @return the relative @c y coordinate with respect to  @c r
long drect_abs_to_rel_y80(drect r, ldbl y);


ldbl drect_rel_to_abs_x80(drect r, int x);;
ldbl drect_rel_to_abs_y80(drect r, int y);

/// @brief Relative to absolute coordinates
///
/// Computes the complex number that represents the point
/// at relative coordinates @c x and @c y from the left bottom point
/// of the dyadic rectangle.
///
/// Relative coordinates are  multiple of \f$ 2^{tpow} \f$.
///
/// @param c the result
/// @param r the rectangle
/// @param x the @c x coordinate, relative to the rectangle r
/// @param y the @c y coordinate, relative to the rectangle r
///
/// @return @ref true if successfull, @ref false otherwise
bool drect_rel_to_abs_coords(mpc c, drect r, long x, long y);

/// @brief Absolute to relative coordinates
///
/// \warning The @c mpfr number @c x is immediately converted to @c long @c double
///
/// @param x the @c x absolute coordinate
/// @param r the rectangle
///
/// @return the relative @c x coordinate with respect to  @c r
long drect_abs_to_rel_x(drect r, mpfr_t x);

/// @brief Absolute to relative coordinates
///
/// \warning The @c mpfr number @c y is immediately converted to @c long @c double
///
/// @param y the @c y absolute coordinate
/// @param r the rectangle
///
/// @return the relative @c y coordinate with respect to  @c r
long drect_abs_to_rel_y(drect r, mpfr_t y);

// //////////////////////////////////////////////////
/// MARK: Geometric operations on rectangles
// //////////////////////////////////////////////////

/// @brief Computes the minimal rectangle that contains both given rectangles and stores the result in
/// @c uni.
///
/// Fails if the rectangles do not have the same @c tpow.
///
/// @param uni the first operand and place to store the result
/// @param r the second operand
///
/// @return @ref true if successfull, @ref false otherwise
bool drect_union(drect uni, drect r);

/// @brief Enlarges the rectangle @c uni to contain the point with coordinates @c x and @c y, divided by @c 2^{uni.tpow}.
///
/// \warning for speed, there are no overflow checks performed, nor the validity of @c r
///
/// @param r the rectangle
/// @param x the @c x coordiante of the point
/// @param y the @c y coordinate of the point
void drect_add(drect r, long x, long y);

/// @brief Recomputes the coordinates of the given rectangle @c uni to scale @c tpow, so that the new rectangle
/// is the smallest that contains the original rectangle.
///
/// \warning for speed, there are no overflow checks performed
///
/// @param r the rectangle to be modified
/// @param tpow the new scale of the rectangle
void drect_rescale(drect r, int tpow);

/// @brief Translates @c r by @c 2^{-tpow}(dx,dy).
///
/// \warning for speed, there are no overflow checks performed, nor the validity of @c r
///
/// @param r the rectangle to be translated
/// @param dx delta x
/// @param dy delta y
/// @param tpow the power of two in which @c dx and @c dy are expressed
void drect_translate(drect r, long dx, long dy, int tpow);

/// @brief Translates @c r so that the center of @c r is as close as possible to @c c.
///
/// @param r the partially computed rectangle
/// @param c the desired center
///
/// @return @c true if successfull, @c false otherwise
bool drect_center80(drect r, fp80 c);

/// @brief Translates @c r so that the center of @c r is as close as possible to @c c.
///
/// Due to the limited precision of the coordiantes of @ref drect_struct, the center cannot be exact.
/// This method optionally computes the delta by which the resulted rectangle should be translated
/// so that its center is exactly @c c. If not @c NULL, the @c prec(delta) should be at least
/// @c prec(c)-62.
///
/// If the coordinates of the bottom left corner of the resulted rectangle are out of range, they are set
/// to @c 0 and @c delta is set accordingly. If, in this case, @c delta==NULL, the method fails.
///
/// @param r the partially computed rectangle
/// @param delta if not @c NULL, it is set to @c c-center(r)
/// @param c the desired center
///
/// @return @c true if successfull, @c false otherwise
bool drect_center(drect r, mpc delta, mpc c);

/// @brief Computes the intersection of two rectangles
///
/// @param inter the result
/// @param r1 the first rectangle
/// @param r2 the second rectangle
///
/// @return @c true if successfull, @c false otherwise
bool drect_intersection(drect inter, drect r1, drect r2);

// //////////////////////////////////////////////////
/// MARK: Membership and overlap querries
// //////////////////////////////////////////////////

/// @brief Test membership
///
/// @param r the rectangle
/// @param c the @c fp80 complex number
///
/// @return @c true if @c c belongs to @c r, @c false otherwise
///
/// Rectangles are closed on the left and bottom side, and open on the top and right side.
bool drect_contains80(drect r, fp80 c);

/// @brief Test membership
///
/// Rectangles are closed on the left and bottom side, and open on the top and right side.
///
/// \warning The @c mpc number @c y is immediately converted to @c long @c double.
/// The result may be inaccurate.
///
/// @param r the rectangle
/// @param c the @c mpc complex number
///
/// @return @c true if @c c belongs to @c r, @c false otherwise
bool drect_contains(drect r, mpc c);

/// @brief Test if @c r2 is included in @c r1
///
/// Rectangles are closed on the left and bottom side, and open on the top and right side.
///
/// @param r1 the first rectangle
/// @param r2 the second rectangle
///
/// @return @c true if @c r2 is included in @c r1, @c false otherwise
bool drect_contains_rect(drect r1, drect r2);

/// @brief Test if two rectangles overlap
///
/// Rectangles are closed on the left and bottom side, and open on the top and right side.
///
/// @param r1 the first rectangle
/// @param r2 the second rectangle
///
/// @return @c true if @c r1 and @c r2 intersect, @c false otherwise
bool drect_intersect(drect r1, drect r2);

bool drect_intersects_disk(drect dr, fp80 c, ldbl r);

bool drect_equals(drect r1, drect r2);

// //////////////////////////////////////////////////
/// MARK: Pretty print
// //////////////////////////////////////////////////

/// @brief Prints the definition of the dyadic rectangle @c r into the string @c str.
///
/// @param r the rectangle
/// @param str the string
/// @param maxLen the max number of characters to output in @c str
///
/// @return the number of characters printed in @c str
int drect_print(drect r, char *str, int maxLen);

#endif /* dRect_h */
