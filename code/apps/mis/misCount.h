//
//  misCount.h
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2021.
//
//  Copyright 2019 - 2021 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the Lesser GNU Public Licence v3.0
//

/**
  \file misCount.h
  \brief This app checks the unicity and the total number of the Misiurewicz points of a given type.
 
  It uses the files generated by misCountMain() that normally contain all Misiurewicz points without repetitions. It checks the MD5
  checksum of @c nset files and cetifies that the correct number of Misiurewicz points are present. It does not attempt to prove
  the corectness of each parameter, use @c misProveMain() for that purpose.
 
  The code is writen for transparency, that is, ease of verification. It is neither optimal, nor the shortest possible: that would be using
  the functions in nSet.h. Uses @c loadAndCheckNSet() that is written for this purpose only.
*/

#ifndef misCount_h
#define misCount_h


/// @brief Interprets the command line arguments and launches the final count of the Misiurewicz points.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int mis_count_main(int argc, const char * argv[]);

#endif /* misCount_h */