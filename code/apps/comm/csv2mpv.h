//
//  csv2mpv.h
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2022.
//
//  Copyright 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the GNU Public Licence v3.0
//

/**
  \file mpv2csv.h
 
  \brief This app packs real or complex numbers from text @c CSV files to binary @c .mpv files.
*/

#ifndef csv2mpv_h
#define csv2mpv_h

#include "mpv.h"

bool csv_read_real(mpfr_t x, FILE *f, int maxDigits);
bool csv_read_complex(mpfr_t x, mpfr_t y, FILE *f, int maxDigits);

/// @brief Interprets the command line arguments and launches the conversion of the .csv file to a .mpv file.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int csv2mpv_main(int argc, const char * argv[]);

#endif /* csv2mpv_h */
