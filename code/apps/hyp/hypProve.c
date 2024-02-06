//
//  hypProve.c
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

#include <stdlib.h>
#include <sys/timeb.h>

#include "stopWatch.h"
#include "hypRaw.h"
#include "hypRawCount.h"
#include "io.h"
#include "hypProve.h"
#include "hypCount.h"
#include "mandel.h"


// MARK: static functions

static bool prove(int per, int start, int end, double maxError, double convRadius) {
    int all = hyp_resultsCount(per);
    
    if(start == 0 && end == all) {
        printf("Re-proving all hyperbolic centers of period %d:\n", per);
    } else {
        if(end - start == 1) {
            printf("Re-proving hyperbolic centers of period %d, file %d:\n",
                   per, start);
        } else {
            printf("Re-proving hyperbolic centers of period %d, files %d to %d:\n",
                   per, start, end - 1);
        }
    }
    
    printf(" - each point is a hyperbolic center with max error %lg on each coordiante,\n",
           maxError);
    printf(" - the Newton method of the corresponding polynomial converges in the disk of radius %lg.\n",
           convRadius);
    
    char time[80];
    struct timeb ats, ts;
    ftime(&ats);
    
    time_stamp(time, 80, true, true);
    printf("Started at %s\n\n", time);
    fflush(stdout);
    
    u128_ptr list = NULL;
    ulong count = 0, real = 0, eps = 0, tot = 0;
    char fn[100];
    
    mpc c;
    mpc_init(c, 128);
    
    for (int i = start; i < end; i++) {
        if(hyp_fileName(fn, 99, per, i) > 10) {
            list = count_load_and_check_nset(fn, &count, &real, &eps, true);
        }
        
        if(list == NULL) {
            printf("\nCould not load results file with index %d, will stop here.\n", i);
            fflush(stdout);
            
            return false;
        }

        ftime(&ts);
        
        for (ulong i = 0; i < count; i++) {
            u128_get(c, list + i);
            
            bool proven = mandel_is_hyp(c, per, maxError);
            proven = proven && mandel_conv_np(c, per, maxError, convRadius);
            
            if(! proven) {
                char *pt = stu(list + i);
                
                printf("The point (%s) is not a hyperbolic center of period %d with max error %lg\n\
                       or its Newton method does not converge in the disk of radius %lg\nwill stop here.\n",
                       pt, per, maxError, convRadius);
                
                free(pt);                
                free(list);
                mpc_clear(c);
                
                return false;
            }
        }
        
        free(list);
        tot += count;
        
        lapse(&ts, time);
        printf("  - %ld points of file %s checked in %s\n", count, fn, time);
        fflush(stdout);
    }
    
    mpc_clear(c);
    
    lapse(&ats, time);
    printf("\nAll %ld points checked are hyperbolic centers of period %d. Total time %s\n",
           tot, per, time);
    fflush(stdout);
    
    return true;
}

// MARK: the help system and the main function

static const char* before = "This task re-proves the hyperbolic centers of given period.\n\n";
static const char* after = "\nCan be parralelized by specifying a range of files to check.\n\n";

static const char *parameters[] = {
    "period",
    "start",
    "end"
};

static const char *types[] = {
    "required",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "0",
    "# of files"
};

static const char *descriptions[] = {
    "the period, integer, at least 1, at most 41",
    "the first file index to check",
    "the first file index that will not be checked"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 3;
static const int columnWidths[] = {18, 18, 18};

/// Prints instructions for usage and some details about the command line arguments.
static void help(void) {
    printf("%s", before);
    char format[50];
    snprintf(format, 45, "    %%-%ds %%-%ds %%-%ds %%s\n",
             columnWidths[0], columnWidths[1], columnWidths[2]);
    
    printf(format, headers[0], headers[1], headers[2], headers[3]);
    printf("\n");
    for(int i = 0; i < paramCount; i++) {
        printf(format, parameters[i], types[i], defaults[i], descriptions[i]);
    }
    
    printf("%s", after);
    fflush(stdout);
}

int hyp_prove_main(int argc, const char * argv[]) {
    int per, start = 0, end;
    
    if(argc < 1 || sscanf(argv[0], "%d", &per) < 1 || per < 1 || per > HYP_RAW_MAX_PER) {
        help();
        
        return 1;
    }
    
    end = hyp_resultsCount(per);
    int mend = end;
    if(argc >= 3) {
        if(sscanf(argv[1], "%d", &start) < 1 || start < 0 || start >= end) {
            help();
            
            return 1;
        }
        
        if(sscanf(argv[2], "%d", &end) < 1 || end <= start || end > mend) {
            help();
            
            return 1;
        }
    }
    
    return ! prove(per, start, end, HYP_RAW_PROOF, HYP_RAW_CONV);
}
