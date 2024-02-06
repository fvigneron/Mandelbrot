//
//  hypSeparation.h
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

#ifndef hypSeparation_h
#define hypSeparation_h

/// @brief Interprets the command line arguments and computes a histogram of distances
/// and separation (distance * derivative) between hyperbolic centers of requested periods.
///
/// Produces two CSV files that contain the results from all periods with several columns, each
/// containing the count of values in a dyadic interval. All roots of each polynomial @c p_n are
/// taken into account.
///
/// @b Important: it needs that the roots are checked for unicity and reordered by
/// @c hypCountMain() and stored in files that can be accessed using @c hyp_fileName()
/// and @c hyp_loadResults().
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int hyp_separation_main(int argc, const char * argv[]);

#endif /* hypSeparation_h */
