//
//  misMinDist.h
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2021.
//
//  Copyright 2019 - 2021 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the Lesser GNU Public Licence v3.0
//

/**
  \file misMinDist.h
  \brief Computation of the min distance between Misiurewicz points,
 with and without the divisor (pre-)periods.
 
 Gives extra details on vertical strips that accumulate to @c -2.
*/

#ifndef misMinDist_h
#define misMinDist_h

#include "nSet.h"

/// @brief Interprets the command line arguments and computes the minimum distance
/// between Misiurewicz points of requested periods.
///
/// Produces a CSV file that contains the results from all
/// periods with several columns, grouped by their proximity to @c -2. For each location,
/// inlcuding the globality of roots of a given (pre-)period, computes two version: with or without
/// the roots of divisor (pre-)period.
///
/// @b Important: it needs that the roots are checked for unicity and reordered by
/// @c misCountMain() and stored in files that can be accessed using @c mis_fileName()
/// and @c mis_loadResults().
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int mis_min_dist_main(int argc, const char * argv[]);

nset mis_new_set(bool locked);

#endif /* misMinDist_h */
