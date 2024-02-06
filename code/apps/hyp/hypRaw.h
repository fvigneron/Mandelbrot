//
//  hypRaw.h
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
  \file hypRaw.h
  \brief This app computes and saves the hyperbolic centers of periods @c 3 to @c 41.
 
 It needs level curves produced by @c levSetsMain() to be available and decomposes
 the task for periods larger than @c 27 into jobs.
 
 Each job produces a set of unique centers and saves it into a binary @c .nset file. However, there
 is always some overlap between different files. The unicity and the presence of all hyperbolic centers
 has to be checked by @c hypCountMain().
*/

#ifndef hypRaw_h 
#define hypRaw_h

#include "nSet.h"

#define HYP_RAW_MAX_PER     41
#define HYP_RAW_SET_2POW    38
#define HYP_RAW_SET_EPS    (1L << HYP_RAW_SET_2POW) // nset eps = MIS_RAW_SET_EPS >> 126

/// The default folder of the hyperbolic centers files.
#define HYP_FOLDER         "hyp"

#define HYP_RAW_DIGITS      26                      // default export to csv
#define HYP_RAW_PROOF       1E-30                   // should be < nset eps
#define HYP_RAW_CONV        1E-25                   // should be > 5 * MIS_RAW_PROOF
#define HYP_RAW_ANG         8
#define HYP_RAW_ANG_2POW    3
#define HYP_RAW_R           5

#define HYP_RAW_LC_PREC   120
#define HYP_RAW_RT_PREC   120
#define HYP_RAW_REF_PREC  170

#define HYP_RAW_QCONV_2POW  (128 - HYP_RAW_SET_2POW)
#define HYP_RAW_RCONV_2POW   128

#define HYP_RAW_TMAP_MAX_LEVEL  18
#define HYP_RAW_TMAP_LEVEL_STEP  2

/// @brief Provides the name of the file of hyperbolic centers found by the job with index
///  @c job and of period @c per.
///
/// @param fileName the file name
/// @param maxChars the max length of the file name
/// @param per the period
/// @param job the index of the file (@c 0 for periods samller then @c 28)
///
/// @return the number of characters put into @c fileName, negative if there was an error
int hyp_raw_job_file_name(char fileName[], int maxChars, int per, int job);

/// @brief Provides the name of the tree map file of hyperbolic centers found by the job with index
///  @c job and of period @c per.
///
/// @param fileName the file name
/// @param maxChars the max length of the file name
/// @param per the period
/// @param job the index of the file (@c 0 for periods samller then @c 28)
///
/// @return the number of characters put into @c fileName, negative if there was an error 
int hyp_raw_map_file_name(char fileName[], int maxChars, int per, int job);

/// @brief Loads and returns the set of hyperbolic centers found for period @c per by the job with index
/// @c job by @c hypRawMain().
///
/// @param per the period
/// @param job the job index
///
/// @return the set of hyperbolic centers, @c NULL if an error occurred
nSet_struct *hyp_raw_load_job(int per, int job);

/// @brief Returns the number of jobs needed to compute the hyperbolic centers for the period @c per.
///
/// @param per the period
///
/// @return the number of jobs
int hyp_raw_jobs_count(int per);

/// @brief Allocates memory and prepares an empty and @b locked set with the distance / epsilon used for
/// hyperbolic centers sets.
///
/// @param locked the state of the set after initialization
///
/// @return an empty and locked set.
nset hyp_new_set(bool locked);

/// @brief Runs the search for hyperbolic centers, as indicated by the command line parameters.
///
/// The static function @c help() provides information about these parameters.
/// 
/// The final results for period 1 and 2 are hardcoded.
///
/// @param argc the number of command line parameters
/// @param argv the command line parmeters
///
/// @return @c 0 if successfull, @c 1 otherwise 
int hyp_raw_main(int argc, const char * argv[]);

#endif /* hypRaw_h */
