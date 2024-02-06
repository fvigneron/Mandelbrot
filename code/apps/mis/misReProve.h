//
//  misReProve.h
//  Mandel_v0.6
//
//  Created by MIHALACHE Nicolae on 5/11/21.
//  Copyright © 2021 MIHALACHE Nicolae. All rights reserved.
//

#ifndef misReProve_h
#define misReProve_h

int mis_is_hyp(mpc c, int per);
bool mis_add_to_results(int pp, int per, nset all);

/// @brief Interprets the command line arguments and launches the final proofs of the Misiurewicz points.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int mis_re_prove_main(int argc, const char * argv[]);

#endif /* misReProve_h */
