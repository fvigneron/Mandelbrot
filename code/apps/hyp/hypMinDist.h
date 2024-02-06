//
//  hypMinDist.h
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
  \file hypMinDist.h
  \brief Computation of the min distance between hyperbolic centers of each period,
 with and without the divisor periods.
 
 Gives extra details on vertical strips that accumulate to @c -2.
*/

#ifndef hypMinDist_h
#define hypMinDist_h

#include "nSet.h"

/// Returns the set of hyperbolic centers of period smaller than and dividing @c per.
///
/// @param per the set of roots of period dividing @c per
nset hyp_load_div(int per);

/// @brief Interprets the command line arguments and computes the minimum distance
/// between hyperbolic centers of requested periods. 
///
/// Produces a CSV file that contains the results from all
/// periods with several columns, grouped by their proximity to @c -2. For each location,
/// inlcuding the globality of roots of a given period, computes two version: with or without
/// the roots of divisor period.
///
/// @b Important: it needs that the roots are checked for unicity and reordered by
/// @c hypCountMain() and stored in files that can be accessed using @c hyp_fileName()
/// and @c hyp_loadResults().
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int hyp_min_dist_main(int argc, const char * argv[]);

#endif /* hypMinDist_h */
