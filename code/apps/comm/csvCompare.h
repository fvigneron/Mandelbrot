//
//  csvCompare.h
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2022.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2021.
//
//  Copyright 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the GNU Public Licence v3.0
//

/**
  \file csvCompare.h
 
  \brief This app compares two @c CSV files storing only real values (not limited to one or two per line).
*/

// FIXME: if error is smaller than the precision of the number used, add a warning

#ifndef csvCompare_h
#define csvCompare_h

#include <stdio.h>
#include <mpfr.h>

#include "ntypes.h"

/// The precision of @ref ldbl numbers.
#define LDBL_PREC    64

/// @brief Convenience type for a number of unknown precision.
typedef union {
    ldbl    ld;    ///< the @ref ldbl value
    mpfr_t  mp;    ///< the @c mpfr value
} number_uni;

/// @brief A vector of numbers of the same precision.
typedef struct {
    uint prec;     ///< the precision of the numbers, @ref ldbl if @c prec<=LDBL_PREC, @c mpfr otherwise
    uint count;    ///< the number of numbers in this vector
    ulong line;    ///< the line number in the file, starting with @c 0
    number_uni vals[]; ///< the values
} csvLine;

/// Convenience type for @c CSV lines.
typedef csvLine *csvl;


/// @brief Reads a line of numerical values from a @c CSV file.
///
/// @param f the file
/// @param prec the precision of the results, @c LDBL_PREC or less for @ref ldbl, more for @c mpfr numbers
///
/// @return the line, @c NULL if some error occurred
csvl csvl_read(FILE *f, uint prec);

/// @brief Frees the meory used by the @c CSV line.
///
/// @param line the @c CSV line.
void csvl_free(csvl line);

/// @brief Reads a line of numerical values from a @c CSV file into a @ref fp80 number.
///
/// @param c the complex number
/// @param f the file
///
/// @return @ref true if successfull, @ref false  if some error occurred
bool csv_read_fp80(fp80 c, FILE *f);

/// @brief Compares if two vectors are equal up to given @c error on each coordinate.
///
/// @param l the first vector
/// @param r the second vector
/// @param maxError the max error for values if the two lines are equal, unchanged otherwise
/// @param error the mac error
///
/// @return @ref true is the vectors are valid and equal, @ref false otherwise
bool csvl_eq(csvl l, csvl r, ldbl error, ldbl *maxError);

/// @brief Compare two @c CSV files containing numerical values.
///
/// @param fn1 the path to the first file
/// @param fn2 the path to the second file
/// @param err the max error for each numbers
/// @param prec the precision in bits to be used
/// @param lines the max distance of out of order lines
/// @param maxError the max error for lines that were considered equal
/// @param verbose @ref true to print info to @c stdout and check all lines, @ref false to compare quietly and
///               stop when the first difference is found
///
/// @return @ref true is the files are readable and equal, @ref false otherwise
bool csv_compare(char *fn1, char *fn2, ldbl err, int prec, int lines, ldbl *maxError,
                 bool verbose, int diffCount);

/// @brief Interprets the command line arguments and launches the comparison of two .csv files.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int csv_compare_main(int argc, const char * argv[]);

#endif /* csvCompare_h */
