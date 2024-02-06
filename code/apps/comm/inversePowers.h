//
//  inversePowers.h
//  Mandel
//
//  Created by MIHALACHE Nicolae on 1/9/23.
//  Copyright © 2023 MIHALACHE Nicolae. All rights reserved.
//

#ifndef inversePowers_h
#define inversePowers_h

/// @brief Stores the names of the @c nset input file and the @c mpv output file for invers powers computation.
///
/// @param inFile the input file name
/// @param outFile the output file name
/// @param len the max number of characters to write in each file name
/// @param pp the pre-perid
/// @param per the period
/// @param index the index of the final @c nset file
/// @param pows the number of powers stored in the output files
///
/// @return @ref true if successful, @ref false otherwise
bool pows_file_names(char *inFile, char *outFile, int len, int pp, int per, int index, int pows);

/// @brief Computes the inverse powers of points in a collection of final @c nset files.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int powers_main(int argc, const char * argv[]);

#endif /* inversePowers_h */
