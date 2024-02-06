//
//  nset2csv.h
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2022.
//
//  Copyright 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the GNU Public Licence v3.0
//

/**
  \file nset2csv.h
  \brief This app exports and refines hyperbolic centers and Misiurewicz points from  binary @c .nset files to
 human readable CSV files.
*/

#ifndef nset2csv_h
#define nset2csv_h 

#include "nSet.h"

#define U128_DIGITS     36

/// @brief Refines the root @c sp of @c p_{per} by an adaptive precision Newton method to have at least
/// @c digits of precision.
///
/// This methos does not perform a proof of the result, the error is bounded by the modulus of the last Newton term,
/// which normaly is comparable to the square root of the error.
///
/// @param res the result, that will be modified if needed to store the precise result
/// @param sp the starting point, a good aproximation of the root
/// @param per the period
/// @param digits the number of exact digits after the decimal point
///
/// @return @ref true if the result was stored, @ref false if some error occurred
bool hyp_refine(mpc res, mpc sp, int per, int digits);

/// @brief Refines the roots of @c p_{per} present in the set @c ps to the max precision of @ref nset.
///
/// This methos does not perform proofs of the results, the error is estimated by the modulus of the last Newton term.
/// The results are stored back in @c ps.
///
/// @param ps the set of points, that are refined in place
/// @param per the period
/// @param maxError the max error found, in @c ulp, or simply units of the coordinates of @ref u128.
///
/// @return @ref true if all points of @c ps could be refined and stored back in @c ps, @ref false if some error occurred
bool hyp_refine_nset(nset_t ps, int per, ulong *maxError);

/// @brief Refines the root @c sp of @c p_{pp+per}-p_{pp} by an adaptive precision Newton method to have at least
/// @c digits of precision.
///
/// This methos does not perform a proof of the result, the error is bounded by the modulus of the last Newton term,
/// which normaly is comparable to the square root of the error.
///
/// @param res the result, that will be modified if needed to store the precise result
/// @param sp the starting point, a good aproximation of the root
/// @param pp the pre-period
/// @param per the period
/// @param digits the number of exact digits after the decimal point
///
/// @return @ref true if the result was stored, @ref false if some error occurred
bool mis_refine(mpc res, mpc sp, int pp, int per, int digits);

/// @brief Refines the roots of @c p_{pp+per}-p_{pp} present in the set @c ps to the max precision of @ref nset.
///
/// This methos does not perform proofs of the results, the error is estimated by the modulus of the last Newton term.
/// The results are stored back in @c ps.
///
/// @param ps the set of points, that are refined in place
/// @param pp the pre-period
/// @param per the period
/// @param maxError the max error found, in @c ulp, or simply units of the coordinates of @ref u128.
///
/// @return @ref true if all points of @c ps could be refined and stored back in @c ps, @ref false if some error occurred
bool mis_refine_nset(nset_t ps, int pp, int per, ulong *maxError);

/// @brief Reads a portion of hyperbolic centers file @c input and writes it to a CSV file @c output.
///
/// Outputs the two coordinates of the points with @c dig digits after the decimal point.
/// If @c deg>DIGITS, the hyperbolic centers are refined by the Newton method.
///
/// @param input the name of the input file
/// @param output the name of the output file
/// @param digits the number of digits after the decimal point
/// @param st the index of the first point
/// @param en the index of the last point @c +1
/// @param per the period
/// @param append @ref true to append points to the CSV file, @ref false to erase existing points
///
/// @return @ref true if successfull, @ref false otherwise
bool hyp_export2csv(char *input, char *output, int digits, long st, long en, int per, bool append);

/// @brief Reads a portion of a text file @c input and writes it to another text file @c output.
///
/// @param input the name of the input file
/// @param output the name of the output file
/// @param st the index of the first line, starting with @c 1
/// @param len the number of lines
/// @param lineLen the max length of a line in the input file
/// @param append @ref true to append the lines to the output file, @ref false to erase its existing content
/// @param verbose @ref true to print info to @c stdout and check all lines, @ref false to compare quietly and
///               stop when the first difference is found
///
/// @return @ref true if successfull, @ref false otherwise
bool txt_export_lines(char *input, char *output, long st, long len, int lineLen, bool append, bool verbose);

/// @brief Counts the lines of a text file @c input.
///
/// @param input the name of the input file
/// @param lineLen the max length of a line in the input file 
/// @param verbose @ref true to print info to @c stdout and check all lines, @ref false to compare quietly and
///               stop when the first difference is found 
///
/// @return the number of lines that can be read from the text file @c input
ulong txt_line_count(char *input, int lineLen, bool verbose);

/// @brief Interprets the command line arguments and launches the conversion of the .nset file to a .csv file.
///
/// Optionally, it refines the hyperbolic points if the desired number of digits is larger than the stored precision.
///
/// @param argc the number of command line arguments
/// @param argv the command line arguments
///
/// @return @c 0 if successfull, @c 1 otherwise
int nset2csv_main(int argc, const char * argv[]);

#endif /* nset2csv_h */
