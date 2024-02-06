//
//  misQuick.h
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
  \dir mis
  \brief Apps that search, check, refine and export the Misiurewicz (or pre-periodic) parameters of the Mandelbrot set.
 
  For quick results with precision of about \f$ 10^{-18} \f$ and periods \f$ n < 30 \f$, use @c misQuickMain().
  To invoke this app, use @c -misQuick as the first parameter in the command line.
 
 For more precise results (proven to have error smaller than \f$ 10^{-30} \f$, but refined to about @c 36 correct digits
 after the decinal point), use @c misRawMain() (command line parameter @c -misRaw). For periods \f$ n < 28 \f$, it
 produces a unique set of points, stored both in CSV format and in a binary file that is easier to manipulate. Starting with
 period \f$ n = 28 \f$, the search is divided into \f$ 2^{n-27} \f$ @c jobs, each producing a binary file of unique
 pre-periodic parameters. However, global unicity cannot be guaranteed anymore. Several instances of @c misRawMain() can run
 simultaneously, or even on distinct machines, the results should be gathered together before be following step.
 
 The app @c misCountMain() (command line parameter @c -misRawCount) will analyse all the raw results for a given period
 and produce a list with unique pre-periodic parameters of period \f$ n \f$, ordered increasingly (by the @c x coordiante first,
 and by the @c y coordinate in case of equality).
 
 The app @c nset2csvMain() can be used to export resulte from the binary files to CSV files (both human readable, but also standard
 text / spreadsheet data files). Optionally, it can refine those results by the Newton method (proven to correctly converge for each
 point by @c misRawMain()) to a precision up to @c 1500 digits.
 
 The app @c levSetsMain() has to be executed before using @c misRawMain() for periods @c 28 or higher, to produce the level curves,
 needed for parallelization.
 
 The app @c misCheckMain() checks the integrity of the files (but not the corectness of the results) produced by @c misRawMain().
 */

 /**
  \file misQuick.h
  \brief This app quickly finds pre-periodic parameters without proofs.
  
  As a hard limit of the @c fp80 precision, the sum of the @c period and the @c pre-period is limited to @c 30.
  
  We have checked the results to coincide with the proven results, the maximal error is bounded by \f$ x.y \cdot 10^{-t} \f$,
  which is close the representation limit of the @c fp80 numbers (@c long @c double on @c x86 machines).
  If one wants to understand the algorithm behind @c hypRawMain() and of @c misRawMain(), this is
  the perfect place to start, as a complement to the reference [1] below.
  
  [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
 */

#ifndef misQuick_h
#define misQuick_h

#include <stdio.h>

/// @brief Interprets the command line arguments and computes the Misiurewicz (or pre-periodic)
/// parameters of requested (pre-)periods, with low precision.
///
/// To compute a the points of pre-period @c pp and period @c per, the points for which @c pp+per is strictly
/// smaller have to be available.
///
/// Produces a CSV file for each type @c (pp,per), with the approximative parameter in increasing order. It works for
/// types (here the sum @c pp+per) from @c 3 to @c 30, higer types need better precision than the @c fp80.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int mis_quick_main(int argc, const char * argv[]);

#endif /* misQuick_h */
