//
//  misTreeMap.h
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
  \file misTreeMap.h
  \brief This app computes the tree map of the Misiurewicz points of a given type.
 
  It uses the files generated by misCountMain() that contain all hyperbolic centers without repetitions. Several maps
  with different resolutions (power of two or the max level of the tree) can be produced in one pass. As the @ref nset files
  produced by misCountMain() have an MD5 checksum included, this is also a test that those files are correctly written to
  the storage medium.
 
  Each map generated by this app contains an MD5 checksum and is reloaded after it has been saved. If no error is reported,
  this is a guarantee that the file is not corrupted and the measure (or planar histogram) of the Misiurewicz points can be published.
*/

#ifndef misTreeMap_h
#define misTreeMap_h

/// @brief Interprets the command line arguments, and produces the requested tree maps of the results.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int mis_tree_main(int argc, const char * argv[]);

#endif /* misTreeMap_h */
