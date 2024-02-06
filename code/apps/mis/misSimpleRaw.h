//
//  misSimpleRaw.h
//  Mandel
//
//  Created by MIHALACHE Nicolae on 12/26/22.
//  Copyright © 2022 MIHALACHE Nicolae. All rights reserved.
//

#ifndef misSimpleRaw_h
#define misSimpleRaw_h

#include "misRaw.h"

/// @brief Runs the search for pre-periodic points, as indicated by the command line parameters.
///
/// The static function @c help() provides information about these parameters.
///
/// @param argc the number of command line parameters
/// @param argv the command line parmeters
///
/// @return @c 0 if successfull, @c 1 otherwise
int mis_simple_raw_main(int argc, const char * argv[]);

#endif /* misSimpleRaw_h */
