//
//  misRaw.h
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
  \file misRaw.h
  \brief   This app computes and saves the pre-periodic points.
 
 It needs level curves produced by @c levSetsMain() to be available and decomposes
 the task for types = pre-period + period larger than @c 27 into jobs.
 
 Each job produces a set of unique pre-periodic points and saves it into a binary @c .nset file. However, there
 is always some overlap between different files. The unicity and the presence of all hyperbolic centers
 has to be checked by @c misCountMain().
*/

#ifndef misRaw_h
#define misRaw_h

#include "nSet.h"

#define MIS_RAW_MAX_TYPE    35
#define MIS_RAW_SET_2POW    16
#define MIS_RAW_SET_EPS    (1L << MIS_RAW_SET_2POW) // nset eps = MIS_RAW_SET_EPS >> 126
#define MIS_FOLDER         "mis"
#define MIS_RAW_DIGITS      28                      // default export to csv
#define MIS_RAW_PROOF       1E-35                   // should be < nset eps
#define MIS_RAW_CONV        1E-32                   // should be > 5 * MIS_RAW_PROOF
#define MIS_RAW_ANG_2POW    3 // at least 3, at most 5
#define MIS_RAW_ANG         (1 << MIS_RAW_ANG_2POW)
#define MIS_RAW_STEP_2POW   1 // >= 0, at most (MIS_RAW_ANG_2POW - 2)
#define MIS_RAW_STEP        (1 << MIS_RAW_STEP_2POW)
#define MIS_RAW_R         100

#define MIS_RAW_LC_PREC   120
#define MIS_RAW_RT_PREC   120
#define MIS_RAW_REF_PREC  170

#define MIS_RAW_QCONV_2POW  (128 - MIS_RAW_SET_2POW)
#define MIS_RAW_RCONV_2POW   128

#define MIS_RAW_TMAP_MAX_LEVEL  18
#define MIS_RAW_TMAP_LEVEL_STEP  2

/// @brief Provides the name of the file of pre-periodic points found by the job with index
///  @c job and of type @c (pp,per).
///
/// @param fileName the file name
/// @param maxChars the max length of the file name
/// @param pp the pre-period
/// @param per the period
/// @param job the index of the file (@c 0 for periods samller then @c 28)
///
/// @return the number of characters put into @c fileName, negative if there was an error
int mis_raw_job_file_name(char fileName[], int maxChars, int pp, int per, int job);

/// @brief Provides the name of the tree map file of pre-periodic points found by the job with index
///  @c job and of type @c (pp,per).
///
/// @param fileName the file name
/// @param maxChars the max length of the file name
/// @param pp the pre-period
/// @param per the period
/// @param job the index of the file (@c 0 for periods samller then @c 28)
///
/// @return the number of characters put into @c fileName, negative if there was an error
int mis_raw_map_file_name(char fileName[], int maxChars, int pp, int per, int job);

/// @brief Loads and returns the set of pre-periodic points found for period @c per by the job with index
/// @c job by @c hypRawMain().
///
/// @param pp the pre-period
/// @param per the period
/// @param job the job index
///
/// @return the set of hyperbolic centers, @c NULL if an error occurred
nSet_struct *mis_raw_load_job(int pp, int per, int job);

/// @brief Returns the number of jobs needed to compute the pre-periodic points of type @c (pp,per).
///
/// @param pp the pre-period
/// @param per the period
///
/// @return the number of jobs
int mis_raw_jobs_count(int pp, int per);

/// @brief Allocates memory and prepares an empty set with the distance / epsilon
/// used for pre-periodic points sets.
///
/// @param locked the state of the set after initialization
///
/// @return an empty set.
nset mis_new_set(bool locked);

/// @brief Computes the pre-periodic points of type @c (pp,per) using as starting points only the segment @c seg
/// of the level curve.
///
/// @param pp the pre-period
/// @param per the period
/// @param seg the segment of the level curve
///
/// @return the set of points, @c NULL if some error occurred.
nset mis_raw_mini(int pp, int per, int seg);

/// @brief Runs the search for pre-periodic points, as indicated by the command line parameters.
///
/// The static function @c help() provides information about these parameters.
///
/// @param argc the number of command line parameters
/// @param argv the command line parmeters
///
/// @return @c 0 if successfull, @c 1 otherwise 
int mis_raw_main(int argc, const char * argv[]);

#endif /* misRaw_h */
