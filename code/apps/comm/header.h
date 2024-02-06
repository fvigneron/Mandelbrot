//
//  header.h
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
  \file nset2csv.h
  \brief This app explains the content of the header of a binary type like @c nset, @c tmap, @c levc or @c bm64.
*/

#ifndef header_h
#define header_h

/// @brief Reads the binary file name from the first parameter and prints the content of its header.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int header_main(int argc, const char * argv[]);

#endif /* header_h */
