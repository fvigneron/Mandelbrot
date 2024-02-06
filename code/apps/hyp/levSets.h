//
//  levSets.h
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the reference below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
//
//  Copyright 2019 - 2021 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the GNU Public Licence v3.0
//

/**
  \file levSets.h
  \brief Computation of the level curves files needed by @c hypRawMain() for high periods.
 
 These level curves are needed  to render the search for hyperbolic centers highly parallelizable.
*/

#ifndef levSets_h
#define levSets_h

#include <stdio.h>

#include "levCurve.h"

/// The default folder of the level sets' files.
#define LEVS_FOLDER      "lev"

/// The power of @c 2 of dyadic angles of the level sets.
#define LEVS_2POW           15

/// The precision of computation of the level sets.
#define LEVS_PREC         120

/// The level of the level sets.
#define LEVS_R              5

/// The default guard for distortion control.
#define LEVS_GUARD          6.1

/// The level of the level sets.
#define LEVS_MAX_PER       41


/// @brief Prints the name of the file containing the level set of period @c per.
///
/// @param fileName the name of the file
/// @param len the max length of the file name
/// @param per the period
///
/// @return @ref true if successfull, @ref false otherwise
bool levs_file_name(char *fileName, int len, int per);

/// @brief Loads the level set of period @c per.
///
/// @param per the period
///
/// @return the level set of period @c per or @c NULL if some error occurred
levc levs_load(int per);

/// @brief Ignores the command line and computes all level sets needed for computing the
/// hyperbolic centers. This should take a few minutes.
///
/// @param argc ignored
/// @param argv ignored
///
/// @return @c 0 if successfull, @c 1 otherwise
int lev_sets_main(int argc, const char * argv[]);

#endif /* levSets_h */
