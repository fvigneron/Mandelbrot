//
//  mpv2csv.h
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
  \brief This app exports real or complex numbers from binary @c .mpv files to human readable CSV files.
*/

#ifndef mpv2csv_h
#define mpv2csv_h

#include "mpv.h"

/// @brief Interprets the command line arguments and launches the conversion of the .mpv file to a .csv file.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int mpv2csv_main(int argc, const char * argv[]);

#endif /* mpv2csv_h */
