//
//  misSimpleSets.h
//  Mandel
//
//  Created by MIHALACHE Nicolae on 12/26/22.
//  Copyright © 2022 MIHALACHE Nicolae. All rights reserved.
//

#ifndef misSimpleSets_h
#define misSimpleSets_h

#include <stdio.h>

#include "misCurve.h"
#include "misSets.h"

/// @brief Writes to @c fileName the name of the file containing the level set of the simplified polynomial of
/// pre-period @c pp and period @c per.
///
/// @param fileName the name of the file
/// @param len the max length of the file name
/// @param pp the pre-period
/// @param per the period
///
/// @return @ref true if successfull, @ref false otherwise
bool missi_file_name(char *fileName, int len, int pp, int per);

/// @brief Loads the level set of the simplified polynomial of pre-period @c pp and period @c per.
///
/// @param pp the pre-period
/// @param per the period
///
/// @return the level set or @c NULL if some error occurred
levm missi_load(int pp, int per);

/// @brief Ignores the command line and computes all level sets needed for computing the
/// Misiurewicz points with simplified polynomials. Total running time is a quarter of an hour or less, depending on the system.
///
/// @param argc ignored
/// @param argv ignored
///
/// @return @c 0 if successfull, @c 1 otherwise
int mis_simple_sets_main(int argc, const char * argv[]);

#endif /* misSimpleSets_h */
