//
//  misSets.h
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
  \brief Computation of the level curves files needed by @c misRawMain() for high degrees.
 
 These level curves are needed  to render the search for Misiurewicz points highly parallelizable.
*/

#ifndef misSets_h
#define misSets_h

#include <stdio.h>

#include "misCurve.h"

/// The default folder of the level sets' files.
#define MISS_FOLDER      "lev"

/// The power of @c 2 of dyadic angles of the level sets.
#define MISS_2POW          15

/// The precision of computation of the level sets.
#define MISS_PREC         120

/// The level of the level sets.
#define MISS_R            100

/// The default guard for distortion control.
#define MISS_GUARD          6.1

/// The first level of the level sets.
#define MISS_MIN_PER       28

/// The max level of the level sets.
#define MISS_MAX_PER       35


/// @brief Writes to @c fileName the name of the file containing the level set of pre-period @c pp and period @c per.
///
/// @param fileName the name of the file
/// @param len the max length of the file name
/// @param pp the pre-period
/// @param per the period
///
/// @return @ref true if successfull, @ref false otherwise
bool miss_file_name(char *fileName, int len, int pp, int per);

/// @brief Loads the level set of pre-period @c pp and period @c per.
///
/// @param pp the pre-period
/// @param per the period
///
/// @return the level set or @c NULL if some error occurred
levm miss_load(int pp, int per);

/// @brief Ignores the command line and computes all level sets needed for computing the
/// Misiurewicz points. Total running time is a quarter of an hour or less, depending on the system.
///
/// @param argc ignored
/// @param argv ignored
///
/// @return @c 0 if successfull, @c 1 otherwise
int mis_sets_main(int argc, const char * argv[]);

#endif /* misSets_h */
