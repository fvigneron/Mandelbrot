//
//  bitmap.h
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
  \dir bitmap
  \brief This folder contains data structures and functions to work with bitmaps.
 
  The bitmaps to represent the Mandelbrot set, Julia sets of quadratic polynomials or of the Newton method.
  The bitmap files can be colloured and collated into images in a post-processing step (with external tools).
 */

/**
 \file bitmap.h
 \brief A data structure and a collection of basic functions to manipulate @c 64-bit bitmaps.
 
 This is an intermediary file format between certified numerical results and actual images. Images can be produced by iteration or
 by matrix representation of [portions of] a \ref treeMap.
 
 @see treeMap.h
 @see iterates.h
 @see dRect.h
*/

#ifndef bitmap64_h
#define bitmap64_h

#include <stdio.h>

#include "fp80.h"
#include "dRect.h"

/// The @c bitmap file ID.
#define BMAP_FILE_ID       "bm64 v01"

/// The length of the header of the file.
#define BMAP_HEADER_MIN_LEN     113

/// Types of bitmaps
#define BMAP_TYPE_GENERIC       0
#define BMAP_TYPE_MANDEL        1
#define BMAP_TYPE_JULIA         2
#define BMAP_TYPE_NEWTON_HYP    3
#define BMAP_TYPE_NEWTON_MIS    4
#define BMAP_TYPE_NEWTON_MISS   5
#define BMAP_TYPE_MEASURE       6
#define BMAP_TYPE_S_FUNCTION    7
#define BMAP_TYPE_VECTORS       8
#define BMAP_TYPE_P_FUNC        9
#define BMAP_TYPE_M_FUNC       10

/// Sub-types of bitmaps
#define BMAP_SUB_TYPE_GENERIC   0
#define BMAP_SUB_TYPE_ESCAPE    1
#define BMAP_SUB_TYPE_DIST      2
#define BMAP_SUB_TYPE_DZ_MOD    3
#define BMAP_SUB_TYPE_DC_MOD    4
#define BMAP_SUB_TYPE_ESC_POW   5
#define BMAP_SUB_TYPE_F_MOD     6    // for functions
#define BMAP_SUB_TYPE_F_ARG     7    // for functions
#define BMAP_SUB_TYPE_CONVERGE  8    // for newton maps
#define BMAP_SUB_TYPE_DENSITY   9    // for measures
#define BMAP_SUB_TYPE_DZ_ARG   10
#define BMAP_SUB_TYPE_DC_ARG   11
#define BMAP_SUB_TYPE_VDENSITY 12    // for vertical measures (each column is a measure)
#define BMAP_SUB_TYPE_HDENSITY 13    // for horizontal measures (each line is a measure)

/// Pixel types
#define BMAP_PIXEL_TYPE_LINEAR  0
#define BMAP_PIXEL_TYPE_ARGB    1
#define BMAP_PIXEL_TYPE_BOOLEAN 2
#define BMAP_PIXEL_TYPE_VECTOR  3

/// Color maps
#define BMAP_COLOR_MAP_LINEAR   0
#define BMAP_COLOR_MAP_POWER    1
#define BMAP_COLOR_MAP_LOG      2
#define BMAP_COLOR_MAP_EXP      3
#define BMAP_COLOR_MAP_ARCTAN   4

/// Pre-defined colors
#define BMAP_COLOR_TRANSPARENT  0x00000000
#define BMAP_COLOR_WHITE        0xFFFFFFFF
#define BMAP_COLOR_BLACK        0xFF000000
#define BMAP_COLOR_RED          0xFFFF0000
#define BMAP_COLOR_GREEN        0xFF00FF00
#define BMAP_COLOR_BLUE         0xFF0000FF
#define BMAP_COLOR_GRAY         0xFF808080
#define BMAP_COLOR_DARK_GRAY    0xFF404040
#define BMAP_COLOR_LIGHT_GRAY   0xFFC0C0C0
#define BMAP_COLOR_YELLOW       0xFFFFFF00
#define BMAP_COLOR_MAGENTA      0xFFFF00FF
#define BMAP_COLOR_CYAN         0xFF00FFFF
#define BMAP_COLOR_RANGE        1024

/// @struct bitMap
/// \brief This data structure caracterizes and stores a bitmap with 64-bit values for each pixel.
///
/// @warning Do not create variables of this type directly, use @ref bmap_new() instead. Otherwise,
/// the size of the member @c pix is unknown or @c 0.
///
/// It represents a rectangle with dyadic coordinates in the plane.
typedef struct {
    byte type;         ///< the type of the image
    byte subType;      ///< the sub-type of the image
    dyadic_rect r;     ///< the bounding rectangle of the bitmap, pixels are squares of size @c 2^{-r.tpow}
    byte pixelType;    ///< the type of the pixel data
    byte sgn;          ///< @c 0 if the values are @c ulong, @c 1 if the values are @c long
    /// @c 2^typeBits types of points are stored in this bitmap; their numerical values is in the highest @c typeBits bits
    byte typeBits;
    /// real values can be obtaind for a pixel by a * pix[i] + b
    byte zeroTransp;   ///< @c 0 if the 0 pixels are normal, @c 1 if they are transparent
    byte colorMap;     ///< the recommended color map
    double a;          ///< multiplicative coefficient
    double b;          ///< additive coefficient
    double power;      ///< the power if the map type is BMAP_COLOR_MAP_POWER
    /// the color is chosen by t = mapA * colorMap(a * pix[i] + b) / BMAP_COLOR_RANGE
    /// if t < 0, the color is @c colorLow
    /// if t > 1, the color is @c colorHigh
    /// otherwise, the color is @c (1-t)*colorFirst+t*colorLast
    double mapA;
    double mapB;
    uint colorLow;     ///< the ARGB color for values below the lowest of the color map
    uint colorFirst;   ///< the ARGB color for the lowest value of the color map
    uint colorLast;    ///< the ARGB color for the highest value of the color map
    uint colorHigh;    ///< the ARGB color for values above the highest of the color map
    bool useHD;        ///< @ref true to use the high definition delta(x, y)
    mpfr_t hdx;        ///< delta x for higher resolution / zoom, w.r.t. r
    mpfr_t hdy;        ///< delta y for higher resolution / zoom, w.r.t. r
    ulong pix[];       ///< the matrix of points / pixels
} bitMap;

/// A pointer to the bitMap struct, for convenience.
typedef bitMap *bmap;

/// A method that maps @c ulong values to @c RGB colors.
typedef uint (* rgb_map)(bmap, ulong);

/// A method that produces an image pixel by pixel. The return value represents and @c RGB color on 24 bits.
typedef uint (* rgb_pix)(int x, int y, void *context);

/// @brief Creates and returns a new bitmap which represents the given dyadic rectangle.
///
/// Use free() to deallocate a \ref bmap, it has no internal pointers.
///
/// @param r the dyadic bounding rectangle
///
/// @return the new bitmap, @c NULL if @c r is invalid
bmap bmap_new(drect r);

bool bmap_fill(bmap m, drect r, ulong pix);

/// @brief Frees the memory used by the @c bitmap.
///
///  If @c bitmap does not use high definition delta(x, y), this can be done by @c free(bitmap).
///
/// @param bitmap the bitmap
void bmap_free(bmap bitmap);

/// @brief Sets the pixel with relative coordinates @c (x,y) to @c val.
///
/// @param bm the bitmap
/// @param val the value
/// @param x the @c x coordinate
/// @param y the @c y coordinate
///
/// @return @ref true if successfull, @ref false otherwise
bool bmap_set_pixel(bmap bm, ulong val, int x, int y);

/// @brief Adds @c 1 to the pixel with relative coordinates @c (x,y).
///
/// @param bm the bitmap
/// @param x the @c x coordinate
/// @param y the @c y coordinate
///
/// @return @ref true if successfull, @ref false otherwise
bool bmap_increment_pixel(bmap bm, int x, int y);

/// @brief Adds @c 1 to the pixel with absolute coordinates @c p to @c val.
///
/// @param bm the bitmap
/// @param p the point
///
/// @return @ref true if successfull, @ref false otherwise
bool bmap_increment_coords(bmap bm, fp80 p);

/// @brief Returns the value of the pixel  with relative coordinates @c (x,y).
///
/// @param bm the bitmap
/// @param x the @c x coordinate
/// @param y the @c y coordinate
///
/// @return the value of the pixel @c (x,y) if successfull, @c 0 otherwise
ulong bmap_get_pixel(bmap bm, int x, int y);

bool bmap_draw(bmap dst, bmap src);
ulong *bmap_histogram(bmap bm, int boxes);
ulong bmap_count(bmap bm, ulong v);

/// @brief If there are no @c typeBits in @c bm, it maps all pixels in the interval @c [min,max] to @c val.
///
/// @warning It fails if @c bm->typBits>0 or fi @c maxVal==0 !
///
/// @param bm the bitmap
/// @param min the min value (included)
/// @param max the max value (included)
/// @param val the end value
///
/// @return @ref true if successfull, @ref false otherwise
bool bmap_filter(bmap bm, ulong min, ulong max, ulong val);

/// @brief If there are no @c typeBits in @c bm, it rescales all the values in @c bm in the range @c 0..maxVal-1 .
///
/// @warning It fails if @c bm->typBits>0 or fi @c maxVal==0 !
///
/// @param bm the bitmap
/// @param maxVal the max value (excluded) after the operation
///
/// @return @ref true if successfull, @ref false otherwise
bool bmap_norm_lin(bmap bm, ulong maxVal);

/// @brief If there are no @c typeBits in @c bm, it rescales all the values in @c bm in the range @c 0..maxVal-1 .
///
/// @warning It fails if @c bm->typBits>0 or fi @c maxVal==0 !
///
/// @param bm the bitmap
/// @param maxVal the max value (excluded) after the operation
/// @param pow the power
/// @param min the minimum (original) value that will be mapped to values > 0
/// @param floor the value to which original values from @c 1 to @c min will be mapped
///
/// @return @ref true if successfull, @ref false otherwise
bool bmap_norm_pow(bmap bm, ulong maxVal, ldbl pow, ulong min, ulong floor);

/// @brief If there are no @c typeBits in @c bm, it rescales all the values in @c bm in the range @c 0..maxVal-1 .
///
/// @warning It fails if @c bm->typBits>0 or fi @c maxVal==0 !
///
/// @param bm the bitmap
/// @param maxVal the max value (excluded) after the operation
///
/// @return @ref true if successfull, @ref false otherwise
bool bmap_norm_log(bmap bm, ulong maxVal);

/// @brief If there are no @c typeBits in @c bm, it rescales the values in the column @c x  of @c bm in the range @c 0..maxVal-1 .
///
/// @warning It fails if @c bm->typBits>0 or fi @c maxVal==0 !
///
/// @param bm the bitmap
/// @param c the column
/// @param maxVal the max value (excluded) after the operation
///
/// @return @ref true if successfull, @ref false otherwise
bool bmap_normalize_column(bmap bm, int x, ulong maxVal);

/// @brief If there are no @c typeBits in @c bm, it rescales the values in the line @c x of @c bm in the range @c 0..maxVal-1 .
///
/// @warning It fails if @c bm->typBits>0 or fi @c maxVal==0 !
///
/// @param bm the bitmap
/// @param c the column
/// @param maxVal the max value (excluded) after the operation
///
/// @return @ref true if successfull, @ref false otherwise
bool bmap_normalize_line(bmap bm, int y, ulong maxVal);

/// @brief Reads and returns a bitmap from the file with given @c fileName.
///
/// @param fileName the file name
/// @param ignoreMD5 @ref true to ignore the @c MD5 checksum stored in the file, @ref false otherwise
///
/// @return the bitmap stored in the file, @c NULL if some error occurred of the @c MD5 checksum is incorrect
bmap bmap_read(char *fileName, bool ignoreMD5);

/// @brief Writes the bitmap to the file with given @c fileName.
///
/// @param bm the bitmap
/// @param fileName the file name
/// @param writeMD5 @ref true to write and @c MD5 checksum to the header of the file, @ref false otherwise
///
/// @return @ref true if successfull, @ref false otherwise
bool bmap_write(bmap bm, char *fileName, bool writeMD5);

/// @brief Writes the given bitmap @c bm to a standard @c BMP file with given @c fileName, using a custom
/// color map from the values stored in @c bm.
///
/// @param bm the bitmap
/// @param fileName the file name
/// @param colors the colors map
///
/// @return @ref true if successfull, @ref false otherwise
bool bmap_write_bmp(bmap bm, char *fileName, rgb_map colors);

/// @brief Writes the given bitmap @c bm to a standard @c BMP file with given @c fileName, using the scaling given by
/// @c bm->a and @c bm->b to the interval @c [0,1].
///
/// The convention is that @c 0 is white, @c 1 is black and intermediary values are lineary mapped to the grayscale.
///
/// @param bm the bitmap
/// @param fileName the file name
/// @param inv @ref true to invert colors, @ref false otherwise
///
/// @return @ref true if successfull, @ref false otherwise
bool bmap_write_bw_bmp(bmap bm, char *fileName, bool inv);

/// @brief Writes a standard @c BMP file with given @c fileName, using a color map from coordinates.
///
/// @param r the bounding rectangle of the image
/// @param fileName the file name
/// @param colors the colors map
/// @param context the @c context to be transmitted to the @c colors map
///
/// @return @ref true if successfull, @ref false otherwise
bool bmp_write(drect r, char *fileName, rgb_pix colors, void *context);

#endif /* bitmap64_h */
