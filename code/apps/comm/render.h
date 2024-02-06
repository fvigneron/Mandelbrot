//
//  render.h
//  Mandel
//
//  Created by MIHALACHE Nicolae on 10/19/22.
//  Copyright © 2022 MIHALACHE Nicolae. All rights reserved.
//

#ifndef render_h
#define render_h

#include <stdio.h>

/// @brief Interprets the command line arguments and launches the bitmap renderings described in a CSV file.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int render_main(int argc, const char * argv[]);

#endif /* render_h */
