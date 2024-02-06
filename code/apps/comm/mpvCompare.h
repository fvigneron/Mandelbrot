//
//  mpvCompare.h
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2022.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2021.
//
//  Copyright 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the GNU Public Licence v3.0
//

/**
  \file csvCompare.h
 
  \brief This app compares two binary @c .mpv files, viewed as vectors of real or complex numbers.
*/

#ifndef mpvCompare_h
#define mpvCompare_h

#include "mpv.h"

/// @brief Interprets the command line arguments and launches the comparison of two .csv files.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int mpv_compare_main(int argc, const char * argv[]);

#endif /* mpvCompare_h */
