//
//  nSet.h
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
  \dir geometry
  @brief This folder contains data structures and functions to work with finte planar sets.
 
  The subfolder \link_to_geometry_bitmap contains the tools to represent the Mandelbrot set and Julia sets of quadratic
  polynomial and of the Newton method as abstract bitmaps. The bitmap can be colloured and collated into images
  in a post-processing step (with external tools).
 */

/**
 \file nSet.h
 @brief A planar set with @ref uint128 coordinates, with a distance @c eps under which two points are considered equal.
 
 For more details about the distance, see the definition of nSet_struct.
 
 It trades [a bit of] speed for a low memory footprint. For the points to be stored in increasing order in memory, the set has to
 be locked by nset_lock(). A locked set cannot add points one by one with nset_add(), but can perform unions
 with nset_union(), provide \a intervals of points (those which project on the given real interval) by nset_interval() and
 be written to binary files with nset_write(). A locked set can unlocked by nset_unlock().
 
 A simpler version of this class is @c planarSet.h, with lower precision (but the same memory footprint).
*/

#ifndef n2Set_h
#define n2Set_h

#include <stdio.h>
#include <mpfr.h>

#include "fp80.h"
#include "mpc.h"
#include "u128c.h"

// MARK: Constants, data structure and types

/// The @c nset file ID.
#define NSET_FILE_ID     "nset v01"

/// The length of the header of the file.
#define NSET_HEADER_LEN    60

/// The default size of the last bar
#define NSET_DEF_SIZE      13

/// The max number of bars
#define NSET_MAX_BARS      32

/// @struct nSet_struct
/// @brief This struct caraterizes and stores a planar set with @c u128c_struct coordinates.
///
/// The first operation should be @c nset_init() and the last operation @c nset_clear(). Two points
/// \f$ p=(x,y) \f$ and \f$ p'=(x',y')\f$ are considered equal if their Chebyshev distance is smaller than @c eps, that is
/// \f[ ||p-p'||_\infty = \max(|x-x'|,|y-y'|)<\mathrm{eps}. \f]
///
/// \warning This is not an equivalence relation, which may cause real issues when using this class in a general setting. We use it
/// for sets which have the minimum distance between points much larger than @c eps. The identification up to @c eps is used
/// to group several approximations of the same point in the abstract set.
typedef struct {
    long count;      ///< the number of points in the set
    long realCount;  ///< the number of real points
    long rejected;   ///< the count of rejected points, as they were already present in the set
    long lastCount;  ///< the number of points present in the last bar
    int barCount;    ///< the number of vectors used to store the points, here named bars
    bool locked;     ///< the state of the set, @ref true means locked, the points are ordered but the set is immutable
    ulong eps;       ///< the distance below which the points are considered equal
    ulong barLen[NSET_MAX_BARS];      ///< the respective lengths of the bars
    u128_ptr pts[NSET_MAX_BARS];      ///< the list of bars that contain the points
    ulong maxMem;
} nSet_struct;

/// Convenience type for a pointer to @c nSet_struct for easy allocation as a local variable.
typedef nSet_struct nset_t[1];

/// Convenience type for a pointer to @c nSet_struct.
typedef nSet_struct *nset; 

// MARK: initialization, validation, lock, unlock and move

/// @brief Initializes the the set with the given minimal distance to distinguish points @c eps (strict inequality).
///
/// A set should be initialized before any other operation. When not in use anymore, call @c nset_clear().
///
/// @param ps the planar set
/// @param eps epsilon @c x2^{126}
void nset_init(nset ps, ulong eps);

/// @brief Allocates and initializes the the set with the given minimal distance to distinguish points @c eps (strict inequality).
///
/// When not in use anymore, call @c nset_clearFree().
///
/// @param eps epsilon @c x2^{126}
/// @param locked the state of the set after initialization
///
/// @return the new set
nset nset_new(ulong eps, bool locked);

/// @brief Frees the memory used by the set @c ps, but not the struct @c ps itself (which may be allocated on the stack).
///
/// Clears all the fields except for @c eps which is unchanged and @c locked which is set to @ref true.
///
/// @param ps the set
///
/// @return @ref true if successfull, @ref false otherwise.
bool nset_clear(nset ps);

/// Frees all the memory occupied by the sets in the list. There may be repetitions in the list.
///
/// @param list the list of sets to clear
/// @param count the number of sets in the list
void nset_clears(nset list[], int count);

/// @brief Frees all the memory used by the set @c ps, assuming the struct has been allocated with @c malloc(), for example
/// with @c nset_new().
///
/// It does @b not fail if @c ps==NULL.
///
/// @param ps the set
void nset_free(nset ps);

/// @brief Checks if the set @c ps is valid.
///
/// @param ps the set
///
/// @return @ref true if @c ps is valid, @ref false otherwise
bool nset_valid(nset ps);

/// @brief Locks the set @c ps. Points cannot be added to a locked set. The operation is reversible.
///
/// @b see nset_unlock()
///
/// @param ps the set
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_lock(nset ps);

/// @brief Unlocks the set @c ps.
///
/// @param ps the set
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_unlock(nset ps);

/// @brief Removes all points from the set @c dst, copies @c src to @c dst and prepares @c src for normal use,
/// that is, it is not locked.
///
/// @param dst the destination set
/// @param src the source set
/// @param lock @ref true to lock both sets, @ref false to leave @c dst in the previous state of @c src and @c src unlocked
void nset_move(nset dst, nset src, bool lock);

/// @brief Searches the point @c p in the given @c locked set @c ps, according to the rule of @c ps:
///
/// if @c |a-b|_1<=ps->eps, then @c a is considered equal to @c b. Returns the position in the @c bar
/// if the point is found, otherwise a negative index @c pos defined as follows: if @c p would be inserted at position
/// @c (-pos-1), then the list would still be increasing.
///
/// The point @c p is checked to be inside the safe rectangle, thus adding @c ps->eps to any
/// of its coordinates does not result in integer overflow.
///
/// @warning Fails if the set is not locked.
///
/// @param ps the set
/// @param p the point to search for
///
/// @return the position of @c p, positive is it is found in the @c bar, otherwise negative, as explained above, @c LONG_MIN if some error occurred
long nset_search(nset ps, u128 p);

// MARK: basic operations

/// @brief Returns the point at position @c pInd in the set @c ps.
///
/// @param ps the set
/// @param pInd the index
/// @param p the point
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_point(u128 p, nset ps, long pInd);

/// @brief Returns the point at position @c pInd in the set @c ps.
///
/// @param ps the set
/// @param pInd the index
/// @param p the point
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_get(mpc p, nset ps, long pInd);

/// @brief Tests if the set @c ps contains the point @c p.
///
/// @param ps the set
/// @param p the point
///
/// @return @ref true if the set @c ps contains the point @c p, @ref false otherwise
bool nset_contains(nset ps, u128 p);

/// @brief If the point @c p is not contained in the set @c ps, then it is added.
///
/// @param ps the set
/// @param p the point
///
/// @return @ref true if the point @c p was added to the set @c ps, @ref false otherwise
bool nset_add(nset ps, u128 p);

/// @brief If the point @c p is not contained in the set @c ps, then it is added.
///
/// @param ps the set
/// @param p the point
///
/// @return @ref true if the point @c p was added to the set @c ps, @ref false otherwise
bool nset_put(nset ps, mpc p);

/// @brief If the point @c p is contained in the set @c ps, then its coordinates are replaced by those of @c p.
///
/// @param ps the set
/// @param p the point
///
/// @return @ref true if the point @c p was replaced in the set @c ps, @ref false otherwise
bool nset_replace(nset ps, u128 p);

/// @brief Reduces the minimal distance between the points of the set @c ps.
///
/// Fails if the set is not valid or the new minimal distance is larger than that of the set.
/// Recounts the number of real points in the set using the new minimal distance @c eps.
///
/// @param ps the set
/// @param eps the new minimal distance
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_set_eps(nset ps, ulong eps);

/// @brief Computes the union of the two sets @c ps and @c ns and stores the result in @c ps. Both sets
/// should be @b locked for this operation to succeed (or @c ps to be empty). Fails if the min dist of the two sets are @b not equal.
///
/// @param ps the fisrt set and the place where the result is stored
/// @param ns the second set
/// @param quick @ref true for a quicker version that may waste some memory
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_union(nset ps, nset ns, bool quick);

/// @brief Computes and returns the intersection of the two sets @c ps and @c ns.
///
/// @warning Both sets should be @b locked for this operation to succeed.
/// Fails if the min dist of the two sets are @b not equal.
///
/// @param ps the fisrt set
/// @param ns the second set
///
/// @return the resulted set, @c NULL if conditions are not satisfied
nset nset_intersection(nset ps, nset ns);

/// @brief Sets @c p to be the rightmost point of the @b locked set @c ps. Fails if @c ps is not locked.
///
/// @param p the point
/// @param ps the set
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_last(u128 p, nset ps);

/// @brief Tests if the point @c p is strictly to the left of the @b locked set @c ps, computed as if @c ps->eps==0.
///
/// @param p the point
/// @param ps the set
///
/// @return @ref true if @c p is smaller than all points of @c ps, @ref false otherwise
bool nset_left(u128 p, nset ps);

/// @brief Tests if the point @c p is strictly to the right of the @b locked set @c ps, computed as if @c ps->eps==0.
///
/// @param p the point
/// @param ps the set
///
/// @return @ref true if @c p is greater than all points of @c ps, @ref false otherwise
bool nset_right(u128 p, nset ps);

/// @brief Removes the points from @c dst and then adds all points from the @b locked set @c src that are in the interval @c [l,r).
/// After this operation @c dst is locked. @b Fails if @c dst==src. The min distance is also copied from @c src to @c dst.
///
/// @param dst the destination set
/// @param src the source set
/// @param l the leftmost point
/// @param r the rightmost point, excluded if it is in @c src
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_interval(nset dst, nset src, u128 l, u128 r);

/// @brief Splits the set @c src into two sets: to the @c left of the point @c m, and to the @c right of @c m.
///
/// After this operation both @c left and @c right are locked. The min distance is also copied from @c src to
/// both @c left and @c right.
///
/// @param left the set on the left of @c m (strict inequality)
/// @param right the set on the right of @c m (large inequality)
/// @param src the source set
/// @param m the split point
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_split(nset left, nset right, nset src, u128 m);

/// @brief Splits the @c locked set @c src into @c parts sets.
///
/// Fails if @c parts==0 or if the number of points in @c src is less than @c parts.
///
/// @param src the source set
/// @param parts the number of parts
///
/// @return a list of @c parts sets, @c NULL if some error has occurred.
nset *nset_divide(nset src, uint parts);

/// @brief Removes the points from @c dst and then adds all points from the @b locked set @c src.
///
/// After this operation @c dst is locked. The min distance is also copied from @c src to @c dst.
///
/// @param dst the destination set
/// @param src the source set
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_copy(nset dst, nset src);

/// @brief Returns the minimal Euclidean distance between two points in the set, @c infinity if the set contains at most one point.
///
/// @param ps the set
///
/// @return the minimal distance between two points in the set, @c infinity if the set contains at most one point,
///         @c -1 if there is an error
ldbl nset_min_dist(nset ps);

/// @brief Returns the minimal Euclidean distance between a point in the set @c ls and a point in the set @c rs,
///
/// with the following conditions: both sets are @b locked, all points in @c ls are stricly smaller than all points in @c rs
/// (first the @c x coordinates are compared, if equal, then the @c y coordiantes are comared) and each set contains
/// at least one point.
///
/// @param ls the set on the left
/// @param rs the set on the right
///
/// @return the minimal distance between one point in @c ls and one point in @c rs, @c infinity if one set is empty,
///  @c -1 if there is an error
ldbl nset_min_dist_lr(nset ls, nset rs);

// MARK: IO operations

/// @brief Reads all points stored in a CSV file and adds them to the @b unlocked and @b initialized set @c ps.
///
/// @param ps the set
/// @param fileName the name of the CSV file
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_read_csv(nset ps, char *fileName);

/// @brief Writes the points of the set @c ps to a CSV file with the given precision in @c digits.
///
/// @param ps the set
/// @param fileName the name of the CSV file
/// @param digits the number of digits after the decinal point
/// @param append @ref true to append points to the CSV file, @ref false to erase existing points
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_write_csv(nset ps, char *fileName, int digits, bool append);

/// @brief Writes a subset of the points of the set @c ps to a CSV file with the given precision in @c digits.
///
/// @param ps the set
/// @param fileName the name of the CSV file
/// @param digits the number of digits after the decinal point
/// @param start the index of the first point to write
/// @param end the index @c +1 of the first last to write
/// @param append @ref true to append points to the CSV file, @ref false to erase existing points
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_write_partial_csv(nset ps, char *fileName, int digits, long start, long end, bool append);

/// @brief Removes the points of the set @c ps and then loads the content of the binary file with the given name @c fileName.
///
/// After this operation, the set @c ps is @b locked.
///
/// @param ps the set of points
/// @param fileName the name of the file
/// @param checkMD5 @ref true to check the MD5 sum of the file
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_read(nset ps, char *fileName, bool checkMD5);

/// @brief Loads the content of the binary file with the given name @c fileName.
///
/// The returned set is @b locked.
///
/// @param fileName the name of the file
/// @param checkMD5 @ref true to check the MD5 sum of the file
///
/// @return the set of points if successfull, @ref NULL otherwise
nset nset_load(char *fileName, bool checkMD5);

/// @brief Partially loads the content of the binary file with the given name @c fileName into a new nSet_struct.
///
/// After this operation, the set @c ps is @b locked. If @c start<0, then it is relative to end of the nSet in the file.
/// The numbers @c start and @c count represent the number of points, not absolute addresses in the binary file.
///
/// @param ps the set
/// @param fileName the name of the file
/// @param pos the absolute postion in the file or @c -1 to use the current position
/// @param start the start position in the vector stored in the file
/// @param count the number of values to read from the file
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_read_partial(nset ps, char *fileName, long pos, long start, ulong count);

/// @brief Removes the points of the set @c ps and then loads the content of the raw / old format
/// binary file with the given name @c fileName.
///
/// After this operation, the set @c ps is @b locked.
///
/// @deprecated
/// @warning This function will not be part of the final release.
///
/// @param ps the set of points
/// @param fileName the name of the file
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_read_old(nset ps, char *fileName);

/// @brief Erases the content and writes the points of the set @c ps to a binary file. Fails if the set @c ps is not @b locked
/// or if it is empty.
///
/// @param ps the set
/// @param fileName the name of the CSV file
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_write(nset ps, char *fileName);

/// @brief Removes the points from @c ps and reads the content starting at position @c pos in the file @c f into @c ps.
///
/// @param ps the set
/// @param f the file
/// @param pos the absolute postion in the file or @c -1 to use the current position
/// @param checkMD5 @ref true to check the MD5 sum of the file 
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_read_from(nset ps, FILE *f, long pos, bool checkMD5);

/// @brief Writes the set @c ps to the file @c f, starting at absolute position @c pos.
/// Fails if the set @c ps is not locked or if it is empty.
///
/// @param ps the set
/// @param f the file
/// @param pos the postion in the file or @c -1 to use the current position
///
/// @return the number of bytes written, @c 0 if some error occurred
ulong nset_write_to(nset ps, FILE *f, long pos);

/// @brief Partially loads the content of the binary file @c f into a new nSet_struct.
///
/// After this operation, the set @c ps is @b locked. If @c start<0, then it is relative to end of the nSet in the file.
/// The numbers @c start and @c count represent the number of points, not absolute addresses in the binary file.
///
/// @param ps the set
/// @param f the file
/// @param pos the absolute postion in the file or @c -1 to use the current position
/// @param start the start position in the vector stored in the file
/// @param count the number of values to read from the file
///
/// @return @ref true if successfull, @ref false otherwise
bool nset_read_segment(nset ps, FILE *f, long pos, long start, ulong count);

#endif /* n2Set_h */
