//
//  main.c
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2022.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2022.
//
//  Copyright 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the Lesser GNU Public Licence v3.0
//

#include <stdio.h>
#include <string.h>

#include "hypRaw.h"
#include "levSets.h"
#include "nset2csv.h"
#include "hypRawCheck.h"
#include "hypRawCount.h"
#include "hypMinDist.h"
#include "hypQuick.h"
#include "fp80.h"
#include "misQuick.h"
#include "misRaw.h"
#include "hypTreeMap.h"
#include "csvCompare.h"
#include "hypCount.h"
#include "hypProve.h"
#include "header.h"
#include "misRawCount.h"
#include "misReProve.h"
#include "misCount.h"
#include "misProve.h"
#include "misSets.h"
#include "mpv2csv.h"
#include "csv2mpv.h"
#include "mpvCompare.h"
#include "mpv2mpv.h"
#include "render.h"
#include "misSimpleQuick.h"
#include "playNicu.h"
#include "misSimpleSets.h"
#include "misSimpleRaw.h"
#include "inversePowers.h"
#include "missingRoots.h"
#include "hypSeparation.h"
#include "misMinDist.h"
#include "compile_time.h"

#include "main.h"

// MARK: the list of available utilities and the help system


static const int commandsCount = 32;
static const char *commands[] = {
    "-hypLevelSets",
    "-nset2csv",
    "-hypRawCheck",
    "-hypRawCount",
    "-hypRaw",
    "-hypMinDist",
    "-hypQuick",
    "-misQuick",
    "-misRaw",
    "-misMinDist",
    "-hypTree",
    "-csvCompare",
    "-hypCount",
    "-hypProve",
    "-header",
    "-misRawCount",
    "-misReProve",
    "-misCount",
    "-misProve",
    "-misLevelSets",
    "-mpv2csv",
    "-csv2mpv",
    "-mpvCompare",
    "-mpv2mpv",
    "-bitmap",
    "-misSimpleQuick",
    "-play",
    "-misSimpleLevelSets",
    "-misSimpleRaw",
    "-inversePowers",
    "-missingRoots",
    "-hypSeparation"
};






static mainFunc programs[] = {
    &lev_sets_main,
    &nset2csv_main,
    &hyp_check_main,
    &hyp_raw_count_main,
    &hyp_raw_main,
    &hyp_min_dist_main,
    &hypQuickMain,
    &mis_quick_main,
    &mis_raw_main,
    &mis_min_dist_main,
    &hyp_tree_main,
    &csv_compare_main,
    &hyp_count_main,
    &hyp_prove_main,
    &header_main,
    &mis_raw_count_main,
    &mis_re_prove_main,
    &mis_count_main,
    &mis_prove_main,
    &mis_sets_main,
    &mpv2csv_main,
    &csv2mpv_main,
    &mpv_compare_main,
    &mpv2mpv_main,
    &render_main,
    &mis_simple_quick_main,
    &playNicuMain,
    &mis_simple_sets_main,
    &mis_simple_raw_main,
    &powers_main,
    &missing_main,
    &hyp_separation_main
};






static const char *parameters[] = {
    "",
    "inFile outFile [digits] [start end]",
    "period [start end]",
    "period cuts [buffLen]",
    "period [start end]",
    "start [end] [intervals]",
    "start [end]",
    "start [end]",
    "pre-period period [start] [end]",
    "start [end] [intervals]",
    "period [minLevel] [maxLevel]",
    "file1 file2 [error] [prec] [lines]",
    "period",
    "period [start end]",
    "inFile",
    "pre-period period cuts [buffLen]",
    "pre-period period [radius] [conv] [prec]",
    "pre-period period",
    "pre-period period [start] [end]",
    "",
    "src dst comp [dig] [start end]",
    "src dst comp [dig] [pack] [append]",
    "src dst comp",
    "src dst comp [dig] [start end] [tpack]",
    "descFile [threads start end]",
    "start [end refine]",
    "depends",
    "",
    "pre-period period [start] [end]",
    "pre-period period [pows] [err] [start] [end]",
    "pre-period period [pows] [err]",
    "start [end interSep interDist]"
};





static const char *descriptions[] = {
    "prepares the level sets for roots search",
    "exports data from binary .nset files to .csv files",
    "checks the results of hypRaw",
    "counts the results of hypRaw",
    "computes the hyperbolic centers, saves results in ~1GB binary files",
    "computes minimum distance between hyperbolic centers",
    "computes hyperbolic centers quickly, with low precision",
    "computes pre-periodic points quickly, with low precision",
    "computes the pre-periodic points, saves results in ~1GB binary files",
    "computes minimum distance between pre-periodic points",
    "computes the maps of the hyperbolic centers with different resolutions",
    "compare CSV files with numerical values",
    "counts and checks the unicity of hyperbolic centers",
    "re-proves the hyperbolic centers and the convergence of the Newton map",
    "explains the header of the input file",
    "counts the results of misRaw",
    "re-checks the failed proofs of misRaw",
    "counts and checks the unicity of Misiurewicz parameters",
    "re-proves the Misiurewicz points and the convergence of the Newton map",
    "prepares the level sets for Misiurewicz points search",
    "exports data from binary .mpv files to .csv files",
    "packs numbers from .csv files to binary .mpv files",
    "compare the values stored in two binary .mpv files",
    "exports [and converts] [partial] data from binary .mpv files to .mpv files",
    "renders a range of the bitmaps described in the input CSV file",
    "computes pre-periodic points quickly, with low precision, using simplified polynomials",
    "quick tests",
    "prepares the simple level sets for Misiurewicz points search",
    "computes pre-periodic points using simplified polynomials",
    "computes the inverse powers of final .nset files",
    "searches for missing roots computed by -inversePowers",
    "computes the separation and min distance histograms"
};

static const char *headers[] = {
    "Task",
    "Parameters",
    "Short description"
};

static const char* before = "\nThe following tasks can be performed with the corresponding command line arguments below:\n\n";
static const char* after = "\n\nFor each of the tasks above, use -task -help for more detailed information.\n\n";
static const int columnWidths[] = {22, 50};
static const int perm[] = { 6,  7, 25,                                      // quick
                            0,  4,  2,  3, 12, 13,                          // hyp
                           19,  8, 15, 17, 18, 16,                          // mis
                           27, 28,                                          // simple mis
                            1, 11, 14, 20, 21, 22, 23,                      // general
                            5,  9, 10, 31,                                  // stats
                           29, 30,                                          // specialized
                           24, 26};                                         // misc

/// Prints instructions for usage and some details about the command line arguments.
static void help(void) {
    printf("%s", before);
    char format[50];
    snprintf(format, 45, "    %%-%ds %%-%ds %%s\n", columnWidths[0], columnWidths[1]);
    
    printf(format, headers[0], headers[1], headers[2]);
    printf("\n");
    for(int i = 0; i < commandsCount; i++) {
        int p = perm[i];
        
        printf(format, commands[p], parameters[p], descriptions[p]);
    }
    
    printf("%s", after);
}

/// The entry point in the app. Depending on the first agument, it launches one of the incorporated
/// utilisites or prints a @c help() screen.
///
/// @param argc the number of command line argument
/// @param argv the command line arguments
/// @return @c 0 if successfull, @c 1 otherwise
int main(int argc, const char * argv[]) {
    printf("Mandelbrot app, compile time : %d-%02d-%02d @ %02d:%02d:%02d\n\n",
           __TIME_YEARS__, __TIME_MONTH__, __TIME_DAYS__, __TIME_HOURS__, __TIME_MINUTES__, __TIME_SECONDS__);
    fflush(stdout);
    
    if(! ntypes_check()) {
        printf("Sorry, your system is not 64-bit or long double compiled as fp64 instead of fp80.\nPlease check and unpdate.\n");
        
        return 1;
    }
    
    if(argc < 2) {
        help();
        
        return 1;
    }
    
    const char *command = argv[1];
    for(int i = 0; i < commandsCount; i++) {
        if(! strncmp(command, commands[i], 30)) {
            return programs[i](argc - 2, argv + 2);
        }
    }
    
    help();
}
