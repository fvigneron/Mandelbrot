//
//  mpv2mpv.h
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2022.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2022.
//
//  Copyright 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the Lesser GNU Public Licence v3.0
//

/**
  \file mpv2csv.h
  \brief This app [partially] exports real or complex numbers [with reduced precision] from binary a @c .mpv file
 [containing mini-files] to another @c mpv binary files [containing mini-files].
*/

#ifndef mpv2mpv_h
#define mpv2mpv_h

#include "mpv.h"

/// @brief Interprets the command line arguments and launches the export from a .mpv file to another.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int mpv2mpv_main(int argc, const char * argv[]);

#endif /* mpv2mpv_h */
