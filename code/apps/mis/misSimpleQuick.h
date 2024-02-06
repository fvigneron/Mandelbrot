//
//  misSimpleQuick.h
//  Mandel
//
//  Created by MIHALACHE Nicolae on 12/24/22.
//  Copyright © 2022 MIHALACHE Nicolae. All rights reserved.
//

#ifndef misSimpleQuick_h
#define misSimpleQuick_h

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
/// It uses reduced polynomial, in contrast to @ref misQuickMain().
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int mis_simple_quick_main(int argc, const char * argv[]);

#endif /* misSimpleQuick_h */
