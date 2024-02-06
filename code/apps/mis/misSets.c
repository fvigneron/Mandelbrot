//
//  misSets.c
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the reference below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
//
//  Copyright 2019 - 2021 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the GNU Public Licence v3.0
//

/**
 \file misSets.c
 \brief Implements the functions defined in misSets.h
 */

#include <stdlib.h>
#include <sys/timeb.h>

#include "misSets.h"
#include "levSets.h"
#include "stopWatch.h"
#include "io.h"

// MARK: static functions

/// @brief Effectively computes the level sets for Misiurewicz points.
///
/// @return @ref true if successfull, @ref false otherwise
static bool compute_mis_sets(void) {
    struct timeb te, ta, tt;
    char time[80], fn[100];
    
    ftime(&ta);
    ftime(&tt);

    for (int tp = MISS_MIN_PER; tp <= MISS_MAX_PER; tp ++ ) {
        printf("Loading the level set of period %d ... ", tp);
        ftime(&te);
        
        levc lc = levs_load(tp);
        if(lc == NULL) {
            printf("failed !");
            
            return false;
        }
        
        lapse(&te, time);
        printf("done in %s.\n", time);
        
        printf("Lifting it to level %d ... ", MISS_R);
        ftime(&te);

        mpfr_t r;
        mpfr_init2(r, MISS_PREC);
        mpfr_set_ld(r, MISS_R, MPFR_RNDN);

        levc llc = levc_level(lc, r);
        levc_free(lc);
        
        mpfr_clear(r);
        
        if(llc == NULL) {
            printf("failed !");

            return false;
        }
        
        lapse(&te, time);
        printf("done in %s.\n", time);
        
        for (int pp = 2; pp < tp; pp ++) {
            int per = tp - pp;
            printf("Deriving the level curve for Mis(%d, %d) ... ", pp, per);
            ftime(&te);
            
            levm lm = levm_from_hyp(llc, pp, per, MISS_GUARD);
            
            if(lm == NULL) {
                printf("failed !\n");
                levc_free(llc);
                
                return false;
            }
            
            lapse(&te, time);
            printf("done in %s. ", time);
            
            ftime(&te);
            if(miss_file_name(fn, 99, pp, per) && levm_write_to(lm, fn)) {
                lapse(&te, time);
                printf("Saved in %s.\n", time);
            } else {
                printf("\nCould NOT write the level set !!!\n");
            }
        }
        
        levc_free(llc);
        
        lapse(&ta, time);
        printf("All level sets of type %d computed and saved in %s.\n\n", tp, time);
        
    }
    
    lapse(&tt, time);
    printf("Total computing time: %s.\n", time);
    
    return true;
}

// MARK: IO operations

bool miss_file_name(char *fileName, int len, int pp, int per) {
    int n = pp + per;
    if(len < 30 || n < MISS_MIN_PER || n > MISS_MAX_PER) {
        return false;
    }
    
    return snprintf(fileName, len, "%s/misSet%d_%d.levm", MISS_FOLDER, pp, per) < len;
}

levm miss_load(int pp, int per) {
    if(pp + per < 3) {
        return NULL;
    }
    
    if(pp + per < MISS_MIN_PER) {
        mpfr_t r;
        mpfr_init2(r, MISS_PREC);
        mpfr_set_ld(r, MISS_R, MPFR_RNDN);
        
        levm lm = levm_new(pp, per, 1, 120, r, MISS_GUARD);
        
        mpfr_clear(r);
        
        return lm;
    }
    
    char fn[100];
    
    return miss_file_name(fn, 99, pp, per) ? levm_read_from(fn) : NULL;
}


// MARK: the main function

int mis_sets_main(int argc, const char * argv[]) {
    // prepare the output folder
    dir(MISS_FOLDER);
    
    return ! compute_mis_sets();
}

