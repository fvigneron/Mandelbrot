//
//  mpv.h
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2022.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2021.
//
//  Copyright 2019 - 2021 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the GNU Public Licence v3.0
//

 /**
  \file mpv.h
  \brief A data structure and a collection of basic functions to use a memory efficient vector of real and complex numbers.
  
  The functions @c mpv_use... and @c mpv_sync... are performance oriented and extra care should be
  exerced when using them.
  
  The numbers with arbitrary precision are based on [mpfr] (https://www.mpfr.org).
  
  When read-only, it is an immuable object and it is @c thread-safe.
 */

#ifndef mpVector_h
#define mpVector_h

#include <mpfr.h>
#include <stdio.h>

#include "ntypes.h"
#include "mpc.h"
#include "mpi.h"
#include "memFile.h"


// MARK: Constants, data structure and types

/// The @c mpv file ID.
#define MPV_FILE_ID     "mpv  v01"

/// The length of the header of the file.
#define MPV_HEADER_LEN    48

/// The min exponent of numbers in @c mpv (power of two).
#define MPV_MIN_EXP   ((int) -1000000000)

/// The max exponent of numbers in @c mpv (power of two).
#define MPV_MAX_EXP   ((int) 1000000000)

/// The exponent of @c 0 in @c mpv.
#define MPV_ZERO_EXP  ((uint) (MPV_MIN_EXP - 1))

/// The exponent of @c infinity in @c mpv.
#define MPV_INF_EXP   ((uint) (MPV_MAX_EXP + 1))

/// The exponent of @c NaN in @c mpv.
#define MPV_NAN_EXP   ((uint) (MPV_MIN_EXP - 2))

/// The min precision, in bits.
#define MPV_MIN_PREC  32

/// @struct mpVector
/// \brief This data structure caracterizes and stores a vector of arbitrary precision floating point numbers.
///
/// All numbers in a vector have the same precision.
typedef struct {
    uint prec;       ///< the precision (in bits) of the numbers, at least @c 32
    uint limbs;      ///< for quicker conversions and access
    ulong count;     ///< the number of numbers in this vector
    ulong *vals;     ///< the limbs (or digits) of the numbers
    int *sexp;       ///< signs and exponents of the numbers
} mpVector;

/// Convenience type for a pointer to @c mpVector for easy allocation as a local variable.
typedef mpVector mpv_t[1];

/// Convenince type for a pointer to @c mpVector.
typedef mpVector *mpv;

/// The number of limbs needed to store numbers with precision @c prec.
#define mpv_limbs(prec)  (((prec - 1) / 64) + 1)

/// The size in bytes of one stored number with given precision @c prec.
#define mpv_point_len(prec) (4 + 8 * mpv_limbs(prec))

/// The number of bytes this mpfr_t number will take to be written to a binary file.
#define mpv_mpfrs(prec) (12 + 8 * mpv_limbs(prec))

/// The number of bytes this mpfr_t number will take to be written to a binary file.
#define mpv_mpfrSize(x) (12 + 8 * mpv_limbs(x->_mpfr_prec))

/// The length in the file of the vector, excluding the header.
#define mpv_pointsLen(v) ((v)->count * (4 + 8 * mpv_limbs((v)->prec)))

/// The length of the file of the fully initialized vector, given by a pointer.
#define mpv_fileLen(v) (MPV_HEADER_LEN + mpv_pointsLen(v))

/// The length of the file of a vector of given precision and length.
#define mpv_file_size(prec, count) (MPV_HEADER_LEN + (count) * (4 + 8 * mpv_limbs(prec)))

// MARK: Creation, destruction, cut and paste

/// @brief Creates a new vector with given precision @c prec and length @c count.
///
/// The precision should be at least @c MPV_MIN_PREC and the length larger than @c 0.
/// All numbers are set to @c NaN.
///
/// @param prec the precision (in bits) of the numbers
/// @param count the length of the vector
///
/// @return the new vector, @c NULL if the parameters are invalid
mpv mpv_new(uint prec, ulong count);

/// @brief Initializes an existing vector with given precision @c prec and length @c count.
///
/// The precision should be at least @c MPV_MIN_PREC and the length larger than @c 0.
/// All numbers are set to @c NaN.
///
/// @param vect the  vector
/// @param prec the precision (in bits) of the numbers
/// @param count the length of the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_init(mpv vect, uint prec, ulong count);

/// @brief Frees the memory used by the vector @c vect, but not the struct @c vect itself (which may be allocated on the stack).
///
/// @param vect the vector
void mpv_clear(mpv vect);

/// @brief Frees all the memory used by the vector @c vect, assuming the struct has been allocated with @c malloc(), for example
/// with @c mpv_new(). It does @b not fail if @c vect==NULL.
///
/// @param vect the vector
void mpv_free(mpv vect);

/// @brief Copies the content of the vector @c src to the vector @c dst, if they are compatible, at the specified position.
///
/// @warning Fails if there is not enough capacity in @c dst, use @c mpv_resize() if needed.
///
/// @param dst the destination vector
/// @param dpos the position in @c dst
/// @param src the source vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_copy(mpv dst, ulong dpos, mpv_t src);

/// @brief Copies some of the content of the vector @c src to the vector @c dst, if they are compatible, at the specified position.
///
/// @warning Fails if there is not enough capacity in @c dst, use @c mpv_resize() if needed.
///
/// @param dst the destination vector
/// @param dpos the position in @c dst
/// @param src the source vector
/// @param spos the position in @c src
/// @param count the number of points
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_copy_partial(mpv dst, ulong dpos, mpv_t src, ulong spos, ulong count);

/// @brief Creates an exact copy of the vector @c vect without shared memort.
///
/// @param vect the vector to clone
///
/// @return the new vector, @c NULL if the @c vect is invalid
mpv mpv_clone(mpv vect);

/// @brief Resizes the vector @c vect to the new length @c count. If The new size is smaller than the current size,
/// the trailing values are lost. If the new size is larger, the new numbers are set to @c NaN.
///
/// @param vect the vector
/// @param count the new length
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_resize(mpv vect, ulong count);

/// @brief Resizes the vector @c dst to allow all values of the vector with the same precision @c src to be appended to it.
///
/// @param dst the destination vector that contains the unoin of the two
/// @param src the vector to append
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_concat(mpv dst, mpv src);

/// @brief Creates a new vector containg the concatenation of the two compatible vectors (same precision).
///
/// @param v1 the leading vector
/// @param v2 the trailing vector
///
/// @return the new vector, @c NULL if the parameters are invalid
mpv mpv_join(mpv v1, mpv_t v2);

/// @brief Returns the sub-vector, with numbers starting at position @c pos, advvancing with step @c step and
/// with @c count elements.
///
/// @param vect the original vector
/// @param start the first element to retain
/// @param step the step in indexes
/// @param count the total count of retained elements
///
/// @return the sub-vector or @c NULL if the parameters are incompatible
mpv mpv_sub_vector(mpv vect, ulong start, ulong step, ulong count);

/// @brief Return the sub-vector, with complex numbers starting at position @c pos, advvancing with
/// step @c step and with @c count elements.
///
/// @param vect the original vector
/// @param start the first element to retain
/// @param step the step in indexes
/// @param count the total count of retained elements
///
/// @return the sub-vector or @c NULL if the parameters are incompatible
mpv mpv_sub_vectorc(mpv vect, ulong start, ulong step, ulong count);

// MARK: Storage and recovery of numbers

/// @brief Sets the element with index @c pos in the vector @c vect to the value of the number @c src.
///
/// @param vect the vector
/// @param pos the position in the vector
/// @param src the number
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_set(mpv vect, ulong pos, mpfr_t src);

/// @brief Changes the sign of the number at position @c pos in the vector @c vect.
///
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_neg(mpv vect, ulong pos);

/// @brief Replaces the complex number at position @c pos in the vector @c vect by its conjugate.
///
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_conj(mpv vect, ulong pos);

/// @brief Sets the element with index @c pos in the vector @c vect to the value of the number @c src.
///
/// Low precision version of mpv_set().
///
/// @param vect the vector
/// @param pos the position in the vector
/// @param src the number
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_setl(mpv vect, ulong pos, ldbl src);

/// @brief Sets the element with index @c pos in the vector @c vect to the value of the integer @c src.
///
/// @param vect the vector
/// @param pos the position in the vector
/// @param src the integer
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_set_si(mpv vect, ulong pos, long src);

long mpv_get_si(mpv vect, ulong pos);

/// @brief Sets the element with index @c pos in the vector @c vect to @c zero.
///
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_set_zero(mpv vect, ulong pos);

/// @brief Sets all the elements in the vector @c vect to @c zero.
///
/// @param vect the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_set_all_zero(mpv vect);

/// @brief Sets the element with index @c pos in the vector @c vect to @c NaN.
///
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_set_nan(mpv vect, ulong pos);

/// @brief Sets the element with index @c pos in the vector @c vect to @c infinity.
///
/// @param vect the vector
/// @param pos the position in the vector
/// @param positive @ref true for @c +inf, @ref false for @c -inf
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_set_inf(mpv vect, ulong pos, bool positive);

/// @brief Sets the elements with indexes @c 2*pos and @c 2*pos+1 in the vector @c vect
/// to the real and respectively imaginary part of the complex number @c src.
///
/// @param vect the vector
/// @param pos the position in the vector (as a complex number)
/// @param src the complex number
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_setc(mpv vect, ulong pos, mpc src);

/// @brief Sets the elements with indexes @c 2*pos and @c 2*pos+1 in the vector @c vect
/// to boundary points of the interval @c src.
///
/// @param vect the vector
/// @param pos the position in the vector (as a complex number)
/// @param src the interval
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_seti(mpv vect, ulong pos, mpi src);

/// @brief Sets the elements with indexes @c 2*pos and @c 2*pos+1 in the vector @c vect
/// to the real and respectively imaginary part of the complex number @c src.
///
/// Low precision version of mpv_setc().
///
/// @param vect the vector
/// @param pos the position in the vector (as a complex number)
/// @param src the complex number
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_setcl(mpv vect, ulong pos, fp80 src);

/// @brief Checks if the value of the element with index @c pos in the vector @c vect is @c NaN.
///
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return @ref true if @c vect[pos] is @c NaN, @ref false otherwise
bool mpv_is_nan(mpv vect, ulong pos);

/// @brief Checks if the value of the element with index @c pos in the vector @c vect is @c zero.
///
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return @ref true if @c vect[pos] is @c zero, @ref false otherwise
bool mpv_is_zero(mpv vect, ulong pos);

/// @brief Checks if the value of the element with index @c pos in the vector @c vect is @c infinite.
///
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return @ref true if @c vect[pos] is @c infinite, @ref false otherwise
bool mpv_is_inf(mpv vect, ulong pos);

/// @brief Checks if the value of the element with index @c pos in the vector @c vect is @c regular.
///
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return @ref true if @c vect[pos] is @c regular, @ref false otherwise
bool mpv_is_reg(mpv vect, ulong pos);

/// @brief Checks if all the values in the vector @c vect are finite regular numbers.
///
/// @param vect the vector
///
/// @return @ref true if all values in @c vect are finite, @ref false otherwise
bool mpv_all_finite(mpv vect);

/// @brief Gets the value of the element with index @c pos in the vector @c vect into the number @c src.
///
/// @param dst the number
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_get(mpfr_t dst, mpv vect, ulong pos);

/// @brief Returns the base @c 2 exponent of the element with index @c pos in the vector @c vect.
///
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return the exponent, or @c INT_MIN if @c NaN or some other error occurred
int mpv_get_2exp(mpv vect, ulong pos);

/// @brief Returns the value of the element with index @c pos in the vector @c vect into the number @c src.
///
/// Low precision version of mpv_get().
///
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return the value @c vect[pos], @c NaN if some error occurred
ldbl mpv_getl(mpv vect, ulong pos);

/// @brief Gets the value of the elements with indexes @c 2*pos and @c 2*pos+1 in the vector @c vect into the
/// real and respectively imaginary parts of the complex number @c src.
///
/// @param dst the complex number
/// @param vect the vector
/// @param pos the position in the vector (as a complex number)
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_getc(mpc dst, mpv vect, ulong pos);

/// @brief Gets the value of the elements with indexes @c 2*pos and @c 2*pos+1 in the vector @c vect into the
/// boundary points of the interval @c src.
///
/// @param dst the interval
/// @param vect the vector
/// @param pos the position in the vector (as a complex number)
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_geti(mpi dst, mpv vect, ulong pos);

/// @brief Gets the value of the elements with indexes @c 2*pos and @c 2*pos+1 in the vector @c vect into the
/// real and respectively imaginary parts of the complex number @c src.
///
/// Low precision version of mpv_getc().
///
/// @param dst the complex number
/// @param vect the vector
/// @param pos the position in the vector (as a complex number)
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_getcl(fp80 dst, mpv vect, ulong pos);

/// @brief Combined get and set from onve vector @c src at position @c spos to another @c dst at position @c dpos.
///
/// @param dst the destination vector
/// @param dpos the destination position
/// @param src the source vector
/// @param spos the source position
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_get_set(mpv dst, ulong dpos, mpv src, ulong spos);

/// @brief Combined get and set of complex values from onve vector @c src at position @c spos
/// to another @c dst at position @c dpos.
///
/// @see mpv_getc() and mpv_setc()
///
/// @param dst the destination vector
/// @param dpos the destination position
/// @param src the source vector
/// @param spos the source position
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_get_setc(mpv dst, ulong dpos, mpv src, ulong spos);

// MARK: Shared limbs with external mpfr_t variables

/// @brief Attaches the number @c dst to the position @c pos of the vector @c vect.
///
/// First, the value @c vect[pos] is stored in @c dst, which should @b not be initialized with mpfr_init2().
/// The memory to store the lims (or digits) is common between @c dst and @c vect[pos].
/// The new sign and exponent of @c dst after any operation is @b not automatically updated in
/// @c vect[pos], @ref mpv_sync() should be used to keep them in synch.
/// This is faster (and uses less memory) than using mpv_set() and mpv_get(), as only the sign and exponent
/// need updates.
///
/// @param dst the uninitialized mpfr number
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_use(mpfr_t dst, mpv vect, ulong pos);

/// @brief Sets @c dst to be the number stored in @c vect at position @c pos.
///
/// @warning For speed, it does not check the compatibility of the vector, use @ref mpv_use() for this.
/// The variable @c dst should not be initialized (nor cleared afterwards), otherwise the memory of the
/// limbs will leak !
///
/// @param dst the number
/// @param vect the vector
/// @param pos the position in the vector
void mpv_fuse(mpfr_t dst, mpv vect, ulong pos);

/// @brief Checks if the number @c src is attached to @c vect[pos] by mpv_use() and if so, updates
/// the sign and the exponent of @c scr into @c vect[pos].
///
/// @param vect the vector
/// @param pos the position in the vector
/// @param src the attached mpfr number
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_sync(mpv vect, ulong pos, mpfr_t src);

/// @brief Copies the sign and exponent of the number @c src to the position @c pos of the vector @c vect.
///
/// @warning For speed, it does not check the compatibility of the parameters, use @ref mpv_sync() for this.
///
/// @param vect the vector
/// @param pos the position in the vector
/// @param src the number
///
/// @return @ref true if @c dst is a regular number, @ref false otherwise
bool mpv_qsync(mpv vect, ulong pos, mpfr_t src);

/// @brief Attaches the real and imaginary parts of the complex number @c dst to @c vect[2*pos] and respectively
/// to @c vect[2*pos+1], in the same way as mpv_use().
///
/// A complex number has other internal buffers and cannot be used if it is not inialized with @c mpc_initBuffers().
/// After external operations with @c dst, @c mpv_syncc() shoul be used.
///
/// @warning Do not use with a complex number that was not initialized or which was initialized with
/// @c mpc_init() or created with @c mpc_new().
///
/// @see The @ref mpc complex numbers
///
/// @param dst the complex number to attach
/// @param vect the vector
/// @param pos the position (multiplied by @c 2) in the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_usec(mpc dst, mpv vect, ulong pos);

/// @brief Attaches the bounds of the interval @c dst to @c vect[2*pos] and respectively
/// to @c vect[2*pos+1], in the same way as mpv_use().
///
/// An interval has other internal buffers and cannot be used if it is not inialized mpi_initBuffers().
/// After external operations with @c dst, mpv_synci() shoul be used.
///
/// @warning Do not use with an interval that was not initialized or which was initialized with @c mpi_init()
/// or created with @c mpi_new().
///
/// @see The \ref mpi interval arithmetic
///
/// @param dst the interval to attach
/// @param vect the vector
/// @param pos the position (multiplied by @c 2) in the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_usei(mpi dst, mpv vect, ulong pos);

/// @brief Checks if the complex number @c src is attached to @c vect[pos] by mpv_usec() and if so, updates
/// the signs and the exponents of the real and imaginary parts of @c src into @c vect[pos].
///
/// @param vect the vector
/// @param pos the position in the vector
/// @param src the attached mpc number
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_syncc(mpv vect, ulong pos, mpc src);

/// @brief Checks if the interval @c src is attached to @c vect[pos] by mpv_usei() and if so, updates
/// the signs and the exponents of the boundary points of @c src into @c vect[pos].
///
/// @param vect the vector
/// @param pos the position in the vector
/// @param src the attached mpc number
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_synci(mpv vect, ulong pos, mpi src);

// MARK: Miscellaneous functions

/// @brief Creates a low precision fp80 copy of this vector.
///
/// @param v the vector
///
/// @return the low precision copy of @c v, @c NULL is @c v is not valid
fp80_ptr mpv_to_fp80(mpv v);

/// @brief Multiplies by 2^tpow the number v[pos]. Fails if the new exponent is out of range.
///
/// @param v the vecror
/// @param pos the index in the vector
/// @param tpow the power of 2
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_scale(mpv v, ulong pos, int tpow);

/// @brief Return the largest @c ulp of the elements in the vector @c vect.
///
/// @param ulp the result
/// @param vect the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_ulp(mpfr_t ulp, mpv vect);

// MARK: IO operations

/// @brief Reads the content of the file given by the path @c fileName into a new vector.
///
/// If @c multi and the file contains several several vectors with the same precision, the concatenation of
/// the vectors is returned. If not @c multi, only the first vector in the file is read and returned. If @c multi, but
/// the file contains incompatible vectors, @c NULL is returned.
///
/// @param fileName the path of the file
/// @param multi @ref true to search if the file contains several vectors of the same precision
///
/// @return the new vector with the content of the file, @c NULL if the file is not valid
mpv mpv_read(char *fileName, bool multi);

/// @brief Writes the vector @c vect to the file given by its path  @c fileName.
///
/// @param vect the vector
/// @param fileName the path to the file
/// @param append @c 1 to append to an existing file, @c 0 to erase the conentent of the existing file
///
/// @return the number of bytes written, @c 0 if some error occurred
ulong mpv_write(mpv vect, char *fileName, int append);

bool mpv_write80(fp80_ptr v, ulong len, uint prec, char *fileName);
fp80_ptr mpv_read80(ulong *len, char *fileName);

/// @brief Writes the vector @c vect as a mini-file intto the file given by its path  @c fileName.
///
/// It considers that all mini-files in the file @c fileName have the same size. It may create gaps in the file
/// (padded with @c 0), as the mini-files may be written in arbitrary order.
///
/// @param vect the vector
/// @param fileName the path to the file
/// @param index the index as a mini-file in the file @c fileName
///
/// @return the number of bytes written, @c 0 if some error occurred
ulong mpv_write_mini(mpv vect, char *fileName, int index);

/// @brief Reads the content starting at position @c pos in the file @c f into a new vector.
///
/// In case of succeess, leaves the pointer of the physical file @c f at the end of the logical file (or mini-file) that
/// has been read.
///
/// @param f the file
/// @param pos the absolute postion in the file or @c -1 to use the current position
///
/// @return the new vector with the content of the file, @c NULL if the file is not valid or other error occurred
mpv mpv_read_from(FILE *f, long pos);

/// @brief Reads the content starting at position @c pos in the file @c f into a existing vector,
/// at a give position.
///
/// It attepts to read the entire logical (or min-) file. In case of succeess, leaves the pointer of the physical file @c f
/// at the end of the logical file (or mini-file) that has been read. It checks the @C MD5 sum of the file.
///
/// @param f the file
/// @param fpos the absolute postion in the file or @c -1 to use the current position
///
/// @return the number of values read @c dst if successfull, @c 0 if the file is not valid or other error occurred
ulong mpv_read_from_to(FILE *f, long fpos, mpv dst, ulong dpos);

/// @brief Partially reads the content starting at position @c pos in the file @c f into a new vector.
///
/// The first element that is read is at the position @c start in the vector stored in the file, the second at
/// @c start+step and so on, until @c count values have been read. If any of those values is not accessible,
/// the function fails and returns @c NULL.
///
/// @param f the file
/// @param pos the absolute postion in the file or @c -1 to use the current position
/// @param start the start position in the vector stored in the file
/// @param step the step in stored vector
/// @param count the number of values to read from the file
/// @param complex @ref true to interpret the vector as complex numbers, @ref false to view it an array of real numbers
///
/// @return the new vector with the partial content of the file, @c NULL if the file or the parameters are not valid
mpv mpv_read_partial(FILE *f, long pos, ulong start, ulong step, ulong count, bool complex);

/// @brief Writes the vector @c vect to the file @c f, starting at absolute position @c pos.
///
/// @param vect the vector
/// @param f the file
/// @param pos the postion in the file or @c -1 to use the current position
///
/// @return the number of bytes written, @c 0 if some error occurred
ulong mpv_write_to(mpv vect, FILE *f, long pos);

/// @brief Reads the points of the vector @c vect from the file @c f starting at the current position.
///
/// @param vect the vector
/// @param vpos the start position in the vector
/// @param count the number of points to read
/// @param f the file
/// @param pos the postion in the file or @c -1 to use the current position
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_read_points(mpv vect, ulong vpos, ulong count, FILE *f, long pos);

/// @brief Reads a subset of the points of the vector @c vect from the file @c f.
///
/// @c start, @c step and @c count are in number of real or complex numbers (depending on @c complex),
/// not in absolute address in the file.
///
/// @param vect the vector
/// @param totCount the number of values present in the file for the vector, interpreted as real values
/// @param f the file
/// @param pos the postion in the file or @c -1 to use the current position
/// @param start the start position in the vector stored in the file
/// @param step the step in stored vector
/// @param count the number of values to read from the file
/// @param complex @ref true to interpret the vector as complex numbers, @ref false to view it an array of real numbers
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_read_points_partial(mpv vect, ulong totCount, FILE *f, long pos,
                           ulong start, ulong step, ulong count, bool complex);

/// @brief Writes the points of vector @c vect to the file @c f, starting at the current position.
///
/// @param vect the vector
/// @param f the file
/// @param pos the postion in the file or @c -1 to use the current position
///
/// @return the number of bytes written, @c 0 if some error occurred
ulong mpv_write_points(mpv vect, FILE *f, long pos);

/// @brief Writes a segment of the vector @c v to a human readable @c CSV file.
///
/// @param v the vector
/// @param fn the output file name
/// @param complex @ref true to write complex numbers, @ref false to write real numbers
/// @param digits the number of digits after the decimal point
/// @param start the start position
/// @param count the number of points to export
/// @param append @ref true to append to the output file, @false otherwise
///
/// @return @ref true if successfull, @ref false otherwise
bool mpv_write_csv(mpv v, char *fn, bool complex, int digits, ulong start, ulong count, bool append);

/// @brief Reads the content of the @c CSV file @c fn into a new vector.
///
/// In case of succeess, leaves the pointer of the physical file @c f at the end of the logical file (or mini-file) that
/// has been read.
///
/// @warning The numbers should NOT have more than 20 000 total digits, otherwise they will not be read correctly.
///
/// @param fn the input file name
/// @param prec the precision in bits
/// @param complex @ref true to read complex numbers, @ref false to read real numbers (one per line)
/// @param max the maximum numbers of lines to read
///
/// @return the new vector with the content of the file, @c NULL if the file is not valid or other error occurred
mpv mpv_read_csv(char *fn, int prec, bool complex, long max);

/// @brief Adds the limbs and exponents of the numbers stored in the vector @c vect to the @c MD5 checksum @c md5.
///
/// @param vect the vector of numbers
/// @param md5 the @c MD5 checksum
///
/// @return @ref true if successful, @ref false otherwise
bool mpv_update_md5(mpv vect, MD5_CTX *md5);

/// @brief Adds the limbs and exponents of a segment of the numbers stored in the vector @c vect to the @c MD5 checksum @c md5.
///
/// @param vect the vestor of numbers
/// @param md5 the @c MD5 checksum
/// @param pos the starting position in the vector of numbers
/// @param count the number of elements to add to the sum
///
/// @return @ref true if successful, @ref false otherwise
bool mpv_update_md5_partial(mpv vect, MD5_CTX *md5, ulong pos, ulong count);

/// @brief Writes the number @c x to the binary file @c f at the current position.
///
/// @param x the number
/// @param f the file
///
/// @return the number of bytes written, @c 0 if some error occurred
ulong fwrite_mpfr(mpfr_t x, FILE *f);

/// @brief Reads the @c mpfr number @c x from the binary file @c f at the current position.
///
/// @param x the number
/// @param f the file
///
/// @return @ref true if successfull, @ref false otherwise
bool fread_mpfr(mpfr_t x, FILE *f);

/// @brief Writes the number @c x to the memory file @c m at the current position.
///
/// @param x the number
/// @param m the memory file
///
/// @return the number of bytes written, @c 0 if some error occurred
ulong mwrite_mpfr(mpfr_t x, mfile m);

/// @brief Reads the @c mpfr number @c x from the memory file @c f at the current position.
///
/// @param x the number
/// @param m the memory file
///
/// @return @ref true if successfull, @ref false otherwise
bool mread_mpfr(mpfr_t x, mfile m);

bool dot_prod(mpfr_t p, mpv x, mpv y, uint prec);

bool dot_prodc(mpc p, mpv x, mpv y, uint prec);

bool mpv_norm(mpfr_t norm, mpv v, uint prec);

bool mpv_norm2(mpfr_t norm, mpv v, uint prec);

bool mpv_dist(mpfr_t dist, mpv v, mpv u, uint prec);

bool mpv_dist2(mpfr_t dist2, mpv v, mpv u, uint prec);

bool mpv_normalize(mpv v);

bool mpv_print(mpv v, ulong pos, int digits);
bool mpv_printc(mpv v, ulong pos, int digits);
bool mpv_snprint(char *str, int len, mpv v, ulong pos, int digits);
bool mpv_snprintc(char *str, int len, mpv v, ulong pos, int digits);

#endif /* mpVector_h */
