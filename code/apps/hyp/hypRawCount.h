//
//  hypRawCount.h
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
  \file hypRawCount.h
  \brief This is an app to count hyperbolic centers of given period.
 
 It rearranges the centers in increasing order, without duplicates, and saves them to
 final binary @c .nset files that can be accessed by @c hyp_fileName() and
 @c hyp_loadResults().
*/

#ifndef hypRawCount_h
#define hypRawCount_h

#include <stdio.h>

#include "nSet.h"


/// @brief Provides the name of the file of hyperbolic centers with index @c ind and period @c per.
///
/// @param fileName the file name
/// @param maxChars the max length of the file name
/// @param per the period
/// @param ind the index of the file (@c 0 for periods samller then @c 28)
///
/// @return the number of characters put into @c fileName, negative if there was an error 
int hyp_fileName(char fileName[], int maxChars, int per, int ind);

/// @brief Returns the number of final nset files that contain the hyperbolic centers of period @c per.
///
/// @see @c hyp_loadResults()
///
/// @param per the period
///
/// @return the number of nset files containing the hyperbolic centers of period @c per
int hyp_resultsCount(int per);

/// @brief Loads and returns the set of hyperbolic centers of period @c per stored in the file  with
/// index @c ind by @c hypRawCountMain().
///
/// The returned set is @b locked. The total number of files containing hyperbolic centers of each period can
/// be obtained by @c hyp_resultsCount().
///
/// @param per the period
/// @param ind the  index of the file
///
/// @return the set of hyperbolic centers, @c NULL if an error occurred
nset hyp_loadResults(int per, int ind);

/// @brief  Interprets the command line arguments and launches the count of the hyperbolic centers.
///
/// Optionally, it recomputes the tree map and stores the unique results in increasing order.
/// Those results can be accessed using @c hyp_fileName() and @c hyp_loadResults().
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
/// 
/// @return @c 0 if successfull, @c 1 otherwise
int hyp_raw_count_main(int argc, const char * argv[]);

#endif /* hypRawCount_h */
