//
//  main.h
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
    \mainpage The Mandelbrot project
 
    This project contains a library of functions and data types to treat numerical problems in complex
    dynamics and more generally, in complex analysis. We have tried to use opimal implementation, both
    for memory footprint and speed of execution. The main focus is the @b corectness of the algorithms.
    
    An important example in this direction is the class (or collection of functions) @c mpd, for computation
    with disks instead of floating point complex numbers. This is the optimal planar variant of interval arithmetic
    for complex analysis, as discussed in [1]. With the theoretical support from [1] and [2], all our numerical
    result are \a certified (or numerically proven).
 
    The project contains a collection of apps that can be used to obtain numerical results described in the
    papers [1] and [2] below. They can be used on personal computers. They are also highly parallelizable
    and have been extensively used on the super-computer
    [ROMEO] (https://romeo.univ-reims.fr/pages/aboutUs) of Univ. de Reims Champagne-Ardenne.
 
 \section sectMaths A bit of maths
 
 Despite its ultimate complexity, the Mandelbrot set admits a very simple definition. For some complex number
 \f$ c \in \mathbb C \f$, let \f$ f_c \f$ be the complex map \f$ f_c:\mathbb C \rightarrow \mathbb C \f$ defined
 by \f$ f_c(z):=z^2+c \f$ for all complex numbers \f$ z \in \mathbb C \f$. Then we can define the Mandelbrot Set
 \f[  \mathcal M := \{c \in \mathbb C\ :\ |f_c^n(0)| \leq 2 \}, \f]
 where \f$ |z| \f$ is the modulus of \f$ z \f$ and \f$ f_c^n := f_c \circ f_c \circ \ldots \circ f_c, \ n\f$ times. That is
 the \f$n\f$-th iterate (or composition) of \f$ f_c \f$ (with itself).
 
 This definition renders \f$ \mathcal M \f$ quite accessible numerically. The first images of \f$ \mathcal M \f$ obtained by
 B. Mandelbrot in 1979 have since fueled the imagination of the general public and of many scientists. Important mathematical
 questions about \f$ \mathcal M \f$ remain open, e.g. the Fatou conjecture which states that the set of hyperbolic parameters
 is open and dense is one century old.
 
 Our tools are to be seen as an improved microscope to observe finer \a phenomenae than what we believe was possible
 with previous algorithms and numerical approaches.
 \section sectDetails Some more maths
 
 \image html mandelRays.png
 \image latex mandelRays.png
 
 \section sectCode The code
 
 The code is structured like a collection of inter-dependent libraries, grouped in @c packages (or folders) by the type
 of @c classes (pairs of header @c .h and implementation @c .c files) they contain.
 The folder @link_to_apps contains main.h and main.c that launch the application. The first command line argument selects
 which of the apps will run. If the first argument is missing or is not in the predefiend list, a help screen is printed and the
 main application quits. 
 
 In the subfolders of @link_to_apps, there are several applications, grouped by the problem to solve:
  - @link_to_apps_hyp search for the centers of hyperbolic compenents of \f$ \mathcal M \f$
  - the apps in @link_to_apps_mis search the pre-periodic or Misiurewicz points that are in the boundary \f$ \partial \mathcal M \f$
  - @link_to_apps_psi are used to compute the coefficients of holomorphic maps related to \f$ \mathcal M \f$
  - the folder @link_to_apps_tools contains tools that are related to the implementation and exploitation of the code
 
 The library part of this project cosists of the following main packages:
 
 \warning The clasees have not been developped to be thread safe, some minor modifications would be needed for that.
 The problems we solve can be divided convenably and the apps can be run in parallel processes. The results are easily
 checked, refined and cleared afterwards.
 
 @see [<b> Files List </b>] (files.html)
  
 \section sectCopy Copyright and references
 
  @authors Nicolae Mihalache\n François Vigneron
 
  Please cite the references below if you use or distribute this software.
 
  - [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.\n
  - [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2021.
 
  \copyright Copyright 2019 - 2021 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.\n
             This software is released under the GNU Lesser Public Licence (LGPL) v3.0\n
 */

/**
 \dir apps
 \brief This folder contains the apps available to users and a few general tools.
  
 All apps have a @c main function that satisfies the signature mainFunc() which will be called by main app,
 according to the first parameter in the command line. If the main app is used without any parameter, a short
 help page will be printed.
 */

/**
 \file main.h
 \brief This is the entry point in the main app.
 
 Depending on the first command line parameter,
 it launches an internal app or presents a help screen to guide the user.
*/

#ifndef main_h
#define main_h

#include <stdio.h>
#include "ntypes.h"
#include "io.h"
#include "stopWatch.h"
#include "nSet.h"

// this converts to string
#define STR_(X) #X

// this makes sure the argument is expanded before converting to string
#define STR(X) STR_(X)

#define AREA_FOLDER           "area"

#define NON_PRIM_FOLDER       "satellite"
#define HYP_AREA_FOLDER       "hyperbolic"
#define PRIM_FOLDER           "primitive"
#define PRIM_NSET_FOLDER      "nset"
#define AREA_BIN_FOLDER       "bin"
#define ERRORS_FOLDER         "errors"
#define LOGS_FOLDER           "logs"
#define MAP_FOLDER            "mmap"
#define TEMP_FOLDER           "temp"
#define BITMAP_FOLDER         "bmp"
#define DEPTH_FOLDER          "depth"
#define RUNS_FOLDER           "runs"
#define STATS_FOLDER          "stats"
#define PARTIAL_FOLDER        "partial"
#define SPECTRUM_FOLDER       "spectrum"
#define SPEC_MAP_FOLDER       "map"
#define SPEC_PSI_FOLDER       "psi"
#define SPEC_POLY_FOLDER      "poly"
#define SPEC_COEFF_FOLDER     "coeffs"

/// Comment for production version.
//#define MANDEL_DEBUG

/// The signature of a main function of an app in this project.
typedef int (* mainFunc)(int, const char **);

bool prep_folders(int depth);

bool sat_file_name(char *file_name, int max, int per, int log2size);
bool area_folder_name(char *file_name, int max);
bool area_file_name(char *file_name, int max, int prec, int per, int log2size, int depth, bool folder);
bool area_stats_file_name(char *file_name, int max, int prec, int per, int depth);
bool areas_stats_file_name(char *file_name, int max, int prec, int per, int max_per, int depth);
bool areas_limits_file_name(char *file_name, int max, int prec, int per, int max_per, int depth, int max_depth);

bool prim_file_name(char *file_name, int max, int per, int deep);
bool prim_partial_file_name(char *file_name, int max, int per, int index);
bool prim_distotion_histogram_file_name(char *file_name, int max, int per, int deep);
bool prim_partial_distotion_histogram_file_name(char *file_name, int max, int per, int index);
bool prim_distoted_list_file_name(char *file_name, int max, int per, int deep);
bool prim_partial_distoted_list_file_name(char *file_name, int max, int per, int index);
bool prim_nset_folder_name(char *file_name, int max, int deep);
bool prim_nset_file_name(char *file_name, int max, int per, int deep);

bool spec_map_hist_file_name(char *file_name, int max, int deep, int lev);
bool means_abs_graph_file_name(char *file_name, int max, int deep, int grid);
bool means_rel_graph_file_name(char *file_name, int max, int deep, int grid);
bool means_abs_bmp_file_name(char *file_name, int max, int st, int en);
bool means_rel_bmp_file_name(char *file_name, int max, int st, int en);

bool hyp_area_stats_file_name(char *file_name, int max, int prec, int log2size, int depth);

bool map_file_name(char *file_name, int max, int x, int y, int tp, int deep, bool sub_maps);
bool bitmap_file_name(char *file_name, int max, int tp);
bool map_escape_hist_file_name(char *file_name, int max, int depth);
bool map_escape_list_file_name(char *file_name, int max, int depth);

bool hyp_stat_sat_size_file_name(char *file_name, int max, int log2size, int depth);
bool hyp_stat_prim_size_file_name(char *file_name, int max, int log2size, int depth);
bool hyp_stat_prim_area_file_name(char *file_name, int max, int log2size, int depth);

bool set_sat_text_files(void);
bool set_prim_text_files(int deep);
bool set_area_text_files(int deep);
bool set_map_text_files(int deep);
bool set_spec_text_files(int deep);

void line_separator(int len, bool out_file);
bool out_progress(char *sign, int count);

int hyp_file_name(char fn[], int max, int per, int ind);
int hyp_files_count(int per);
nset hyp_load(int per, int ind);


#endif /* main_h */
