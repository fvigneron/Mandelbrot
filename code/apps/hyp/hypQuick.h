//
//  hypQuick.h
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
  \dir hyp
  \brief This folder contains the apps that search, check, refine and export the hyperbolic centers of the Mandelbrot set.
 
  For quick results with precision of about \f$ 10^{-18} \f$ and periods \f$ n < 30 \f$, use @c hypQuickMain().
  To invoke this app, use @c -hypQuick as the first parameter in the command line.
 
 For more precise results (proven to have error smaller than \f$ 10^{-30} \f$, but refined to about @c 36 correct digits
 after the decinal point), use @c hypRawMain() (command line parameter @c -hypRaw). For periods \f$ n < 28 \f$, it
 produces a unique set of points, stored both in CSV format and in a binary file that is easier to manipulate. Starting with
 period \f$ n = 28 \f$, the search is divided into \f$ 2^{n-27} \f$ @c jobs, each producing a binary file of unique
 hyperbolic centers. However, global unicity cannot be guaranteed anymore. Several instances of @c hypRawMain() can run
 simultaneously, or even on distinct machines, the results should be gathered together before be following step.
 
 The app @c hypCountMain() (command line parameter @c -hypRawCount) will analyse all the raw results for a given period
 and produce a list with unique hyperbolic centers of period \f$ n \f$, ordered increasingly (by the @c x coordiante first,
 and by the @c y coordinate in case of equality).
 
 The app @c nset2csvMain() can be used to export resulte from the binary files to CSV files (both human readable, but also standard
 text / spreadsheet data files). Optionally, it can refine those results by the Newton method (proven to correctly converge for each
 point by @c hypRawMain()) to a precision up to @c 1500 digits.
 
 The app @c levSetsMain() has to be executed before using @c hypRawMain() for periods @c 28 or higher, to produce the level curves,
 needed for parallelization.
 
 The app @c hypCheckMain() checks the integrity of the files (but not the corectness of the results) produced by @c hypRawMain().
 */

/**
 \file hypQuick.h
 \brief This app quickly finds hyperbolic centers without proofs.
 
 As a hard limit of the @c fp80 precision, only periods up to @c 31 are supported.
 
 We have checked the results to coincide with the proven results, the maximal error is bounded by \f$ 5.4 \cdot 10^{-19} \f$,
 which is close the representation limit of the @c fp80 numbers (@c long @c double on @c x86 machines).
 If one wants to understand the algorithm behind @c hypRawMain() and of @c misRawMain(), this is
 the perfect place to start, as a complement to the reference [1] below.
 
 [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
*/

#ifndef hypQuick_h
#define hypQuick_h

#include <stdio.h>

/// @brief Interprets the command line arguments and computes the hyperbolic centers of requested periods,
/// with low precision.
///
/// To compute a period @c p, the periods diving @c p should alredy be computed
/// and the points list in the respective @c .txt file available in the current folder.
///
/// Produces a CSV file for each period, with the approximative roots in increasing order. It works for periods
/// from @c 3 to @c 31, higer periods need better precision than the @c fp80.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int hypQuickMain(int argc, const char * argv[]);

#endif /* hypQuick_h */ 
