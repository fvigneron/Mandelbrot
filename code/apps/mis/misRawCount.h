//
//  misRawCount.h
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
  \file misRawCount.h
  \brief This is an app to count Misiurewicz points of given type.
 
 It rearranges the pre-periodic parameters in increasing order, without duplicates, and saves them to
 final binary @c .nset files that can be accessed by @c mis_fileName() and
 @c mis_loadResults().
*/

#ifndef misRawCount_h
#define misRawCount_h

#include <stdio.h>

#include "nSet.h"


/// @brief Provides the name of the file of pre-periodic parameters  with index @c ind and type @c (pp,per).
///
/// @param fileName the file name
/// @param maxChars the max length of the file name
/// @param pp the pre-period
/// @param per the period
/// @param ind the index of the file (@c 0 for periods samller then @c 28)
///
/// @return the number of characters put into @c fileName, negative if there was an error
int mis_file_name(char fileName[], int maxChars, int pp, int per, int ind);

/// @brief Returns the number of final nset files that contain the pre-periodic parameters of type @c (pp,per).
///
/// @see @c hyp_loadResults()
///
/// @param pp the pre-period
/// @param per the period
///
/// @return the number of nset files containing the pre-periodic parameters of period @c per
int mis_results_count(int pp, int per);

/// @brief Loads and returns the set of pre-periodic points of type @c (pp,per) stored in the file  with
/// index @c ind by @c misRawCountMain().
///
/// The returned set is @b locked. The total number of files containing pre-periodic parameters of each period can
/// be obtained by @c mis_resultsCount().
///
/// @param pp the pre-period
/// @param per the period
/// @param ind the  index of the file
///
/// @return the set of pre-periodic parameters, @c NULL if an error occurred
nset mis_load_results(int pp, int per, int ind);

/// @brief PAtially loads and returns the set of hyperbolic centers of type @c (pp,per) stored in the file  with
/// index @c ind by @c misRawCountMain().
///
/// The returned set is @b locked. The total number of files containing pre-periodic parameters of each period can
/// be obtained by @c mis_resultsCount(). If @c start<0, then it is relative to end of the nSet in the file.
/// The numbers @c start and @c count represent the number of points, not absolute addresses in the binary file.
///
/// @param pp the pre-period
/// @param per the period
/// @param ind the  index of the file
/// @param start the start position in the vector stored in the file
/// @param count the number of values to read from the file
///
/// @return the set of pre-periodic parameters, @c NULL if an error occurred
nset mis_load_partial_results(int pp, int per, int ind, long start, ulong count);

/// @brief Interprets the command line arguments and launches the count of the pre-periodic parameters.
///
/// Optionally, it recomputes the tree map and stores the unique results in increasing order.
/// Those results can be accessed using @c mis_fileName() and @c mis_loadResults().
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int mis_raw_count_main(int argc, const char * argv[]);

#endif /* misRawCount_h */
