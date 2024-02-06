//
//  planarSet.h
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
 \file planarSet.h
 \brief A planar set with @c fp80 coordinates, with a distance @c eps under which two points are considered
 equal.
 
 This is only used by @c hypQuickMain(), as a low resolution version of @c nSet. It is simpler version of @c nSet_struct,
 with lower precision (but the same memory footprint).
*/

#ifndef planarSet_h
#define planarSet_h

#include "fp80.h"

/// the defeult size of the last bar
#define PSET_DEF_SIZE  13
/// the max number of bars
#define PSET_MAX_BARS  30

// MARK: tha main structure definition

/// @struct planarSet_struct
/// @brief This struct caraterizes and stores a planar set with @c fp80 coordinates.
///
/// The first operation should be @c pset_init() and the last operation @c pset_clear().
///
/// @see nSet_struct, nSet.h
typedef struct {
    long count;       ///< the number of points in the set
    long realCount;   ///< the number of real points
    long rejected;    ///< the count of rejected points, as they were already present in the set
    long lastCount;   ///< the number of points present in the last bar
    int barCount;     ///< the number of vectors used to store the points, here named bars
    bool locked;      ///< the state of the set, @c true means locked, the points are ordered but the set is immutable
    double eps;       ///< the distance below which the points are considered equal
    bool absoluteEps; ///< @ref false to use relative eps size
    long barLen[PSET_MAX_BARS];      ///< the respective lengths of the bars
    fp80_ptr pts[PSET_MAX_BARS];    ///< the list of bars that contain the points
} planar_set;

/// A convenience type for easy alocation as local variable and passing as function parameters.
typedef planar_set pset[1];

/// A convenience pointer type.
typedef planar_set *pset_ptr;

// MARK: initialization, validation, lock

/// @brief Initializes the the set with the given minimal distance to distinguish points @c eps (strict inequality).
///
/// A set should be initialized before any other operation. When not in use anymore, call @c pset_clear().
///
/// @param ps the planar set
/// @param eps epsilon
///
/// @return @ref true if all necessary memory could be allocated, @ref false otherwise
bool pset_init(pset ps, double eps);

/// @brief Frees the memory used by the set @c ps, but not the struct @c ps itself (which may be allocated on the stack).
///
/// @param ps the set
void pset_clear(pset ps);

/// @brief Checks if the set @c ps is valid.
///
/// @param ps the set
///
/// @return @ref true if @c ps is valid, @ref false otherwise
bool pset_valid(pset ps);

/// Locks the set @c ps. Points cannot be added to a locked set. The operation is \b irreversible.
///
/// @param ps the set
void pset_lock(pset ps);

// MARK: basic operations

/// @brief Returns the point at position @c pInd in the set @c ps.
///
/// @param p the point
/// @param ps the set
/// @param pInd the index
///
/// @return @ref true if successfull, @ref false otherwise
bool pset_point(fp80 p, pset ps, long pInd);

/// @brief Tests if the set @c ps contains the point @c p.
///
/// @param ps the set
/// @param p the point
///
/// @return @ref true if the set @c ps contains the point @c p, @ref false otherwise
bool pset_contains(pset ps, fp80 p);

/// @brief If the point @c p is not contained in the set @c ps, then it is added.
///
/// @param ps the set
/// @param p the point
///
/// @return @ref true if the point @c p was added to the set @c ps, @ref false otherwise
bool pset_add(pset ps, fp80 p);

/// @brief Retrun the position of the point @c p in the set @c ps, if @c ps is \b locked and @c p is a point in the set.
///
/// @param ps the set
/// @param p the point
///
/// @return the index of @c p in the \b locked set @c ps, @c -1 if @c ps does not contain @c p or some error occurred
long pset_index(pset ps, fp80 p);

// MARK: IO operations

/// @brief Reads all points stored in a CSV file and adds them to the @b unlocked and @b initialized set @c ps.
///
/// @param ps the set
/// @param fileName the name of the CSV file
///
/// @return @ref true if successfull, @ref false otherwise
bool pset_read(pset ps, char *fileName);

/// @brief Erases the content of the file and writes the points of the set @c ps to a CSV file with the given @c 19 digits
/// after the decimal point.
///
/// @param ps the set
/// @param fileName the name of the CSV file
///
/// @return @ref true if successfull, @ref false otherwise
bool pset_write(pset ps, char *fileName, bool append);

#endif /* planarSet_h */
