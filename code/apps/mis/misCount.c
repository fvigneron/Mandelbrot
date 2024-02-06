//
//  misCount.c
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2021.
//
//  Copyright 2019 - 2021 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the Lesser GNU Public Licence v3.0
//

#include <stdlib.h>
#include <sys/timeb.h>

#include "stopWatch.h"
#include "misRaw.h"
#include "misRawCount.h"
#include "misCount.h"
#include "hypCount.h"
#include "mandel.h"

// MARK: static functions

static bool count(int pp, int per) {
    printf("Counting Misiurewicz points of type (%d, %d).", pp, per);
    
    char time[80];
    struct timeb ats, ts;
    ftime(&ats);
    
    time_stamp(time, 80, true, true);
    printf(" Started at %s\n\n", time);
    fflush(stdout);
    
    int res = mis_results_count(pp, per);
    u128_ptr prev = NULL, list = NULL;
    ulong prevc = 0, count = 0, real = 0, eps = 0, totc = 0, totr = 0, prevEps = 0;
    char fn[100];
    
    int n = pp + per;
    int files = n < 28 ? 1 : 1 << (n - 27);
    int digits = 0;
    while(files > 0) {
        files /= 10;
        digits ++;
    }
    
    char pf[20];
    snprintf(pf, 20, " %%0%dd", digits);
    
    for (int i = 0; i < res; i++) {
        list = NULL;
        
        if(mis_file_name(fn, 99, pp, per, i) > 10) {
            list = count_load_and_check_nset(fn, &count, &real, &eps, true);
        }
        
        ftime(&ts);
        
        if(list == NULL) {
            printf("\nCould not load results file with index %d, will stop here.\n", i);
            fflush(stdout);
            
            if(prev != NULL) {
                free(prev);
            }
            
            return false;
        }
        
        // announce next files for external buffering
        printf("        prepare:");
        for (int j = i + 1; j < i + 21 && j < res; j++) {
            printf(pf, j);
        }
        printf("\n");
                
        if(! count_check_increasing(list, count)) {
            printf("The points in %s are not in strictly increasing order, count failed.\n", fn);
            fflush(stdout);
            
            free(list);
            if(prev != NULL) {
                free(prev);
            }
            
            return false;
        }
        
        if(count > 1 && list[count - 1].x - list[0].x <= eps) { // no overflow, as list has been checked to be increasing
            printf("The points in %s are contained in very thin vertical strip, count failed.\n", fn);
            fflush(stdout);
            
            free(list);
            if(prev != NULL) {
                free(prev);
            }
            
            return false;
        }
        
        ulong meps = eps < prevEps ? prevEps : eps;
        if(i > 0 && ! count_check_boundary(prev, prevc, list, count, meps)) {
            printf("Boundary intersection of files of %d and %d, count failed.\n",
                   i -1, i);
            fflush(stdout);
            
            free(list);
            if(prev != NULL) {
                free(prev);
            }
            
            return false;
        }
        
        if(! count_check_min_dist(list, count, eps)) {
            printf("The points in %s are not separated by %lg on at least one coordinate, count failed.\n",
                   fn, ldexp(eps, -126));
            fflush(stdout);
            
            free(list);
            if(prev != NULL) {
                free(prev);
            }
            
            return false;
        }
        
        if(real != count_real_count(list, count, eps)) {
            printf("The number of real parameters stored in %s is incorrect, count failed.\n", fn);
            fflush(stdout);
            
            free(list);
            if(prev != NULL) {
                free(prev);
            }
            
            return false;
        }
        
        if(i > 0) {
            free(prev);
        }
        
        totc += count;
        totr += real;
        
        prev = list;
        prevc = count;
        prevEps = eps;
        
        lapse(&ts, time);
        printf("  - %ld points of file %s checked in %s\n", count, fn, time);
        fflush(stdout);
    }
    
    free(list);
    ulong tr = 2 * totc - totr;
    ulong hc = mandel_mis_count(pp, per);
    
    if(hc != tr) {
        printf("\nThe files contain a list of unique points, but the total count is %ld instead of %ld\n",
               tr, hc);
        fflush(stdout);
        
        return false;
    }
    
    printf("\nThe files contain all %ld Misiurewicz points (%ld stored) of type (%d, %d),\nin strictly increasing order. Any two points are at distance at least %lg from each other.\n\n",
          hc, totc, pp, per, ldexp(eps, -126));
    
    lapse(&ats, time);
    printf("Total loading an computing time %s\n", time);
    fflush(stdout);
    
    return true;
}

// MARK: the help system and the main function

static const char* before = "This task counts the Misiurewicz points of given type, checking for unicity.\n\n";
static const char* after = "\nIt also checks that the nset files are correct with an MD5 checksum.\n\n";

static const char *parameters[] = {
    "pre-period",
    "period"
};

static const char *types[] = {
    "required",
    "required"
};

static const char *defaults[] = {
    "",
    ""
};

static const char *descriptions[] = {
    "the pre-period, integer, at least 2, at most 34",
    "the period, integer, at least 1, at most 35 - pre-period"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 2;
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

int mis_count_main(int argc, const char * argv[]) {
    int pp, per;
    
    if(argc < 2 || sscanf(argv[0], "%d", &pp) < 1  || sscanf(argv[1], "%d", &per) < 1 ||
       pp < 2 || per < 1 || pp + per > MIS_RAW_MAX_TYPE) {
        help();
        
        return 1;
    }
    
    return ! count(pp, per);
}
