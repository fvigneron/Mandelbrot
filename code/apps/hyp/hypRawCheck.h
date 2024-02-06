//
//  hypRawCheck.h
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
  \file hypRawCheck.h
  \brief   This is an app to check the hyperbolic centers files produced by @c hypRawMain().
 
 Optionally, it refines the results with better precision. As a limitation of the .nset files, all results are stored
 with 2 bits precision for the integer part and 126 bits for the fractionary (dyadic) part.
*/

#ifndef hypCheck_h
#define hypCheck_h

/// @brief Interprets the command line arguments, checks the files produced by @c hypRawMain().
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
/// @return @c 0 if successfull, @c 1 otherwise
int hyp_check_main(int argc, const char * argv[]);

#endif /* hypCheck_h */
