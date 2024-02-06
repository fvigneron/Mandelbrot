//
//  misSimpleSets.c
//  Mandel
//
//  Created by MIHALACHE Nicolae on 12/26/22.
//  Copyright © 2022 MIHALACHE Nicolae. All rights reserved.
//

#include <stdlib.h>
#include <sys/timeb.h>

#include "misSimpleSets.h"
#include "levSets.h"
#include "stopWatch.h"
#include "io.h"

// MARK: static functions

/// @brief Effectively computes the level sets for simple Misiurewicz points.
///
/// @return @ref true if successfull, @ref false otherwise
static bool compute_misi_sets(void) {
    printf("Computing the simple level sets for the search of Misiurewicz points.\n\n");
    
    struct timeb te, ta, tt;
    char time[80], fn[100];
    
    ftime(&ta);
    ftime(&tt);

    for (int tp = MISS_MIN_PER; tp <= MISS_MAX_PER; tp ++ ) {
        bool first = tp == MISS_MIN_PER;
        int n = tp - 1;
        printf("Loading the level set of period %d ... ", first ? tp : n);
        ftime(&te);
        
        levc lc = levs_load(first ? tp : n);
        if(lc == NULL) {
            printf("failed !");
            
            return false;
        }
        
        lapse(&te, time);
        printf("done in %s.\n", time);
        
        printf("Lifting it to level %d ... ", first ? MISS_R * MISS_R : MISS_R);
        ftime(&te);

        mpfr_t r;
        mpfr_init2(r, MISS_PREC);
        mpfr_set_ld(r, MISS_R, MPFR_RNDN);
        if(first) {
            mpfr_sqr(r, r, MPFR_RNDN);
        }

        levc llc = levc_level(lc, r);
        levc_free(lc);
        
        mpfr_clear(r);
        
        if(llc == NULL) {
            printf("failed !");

            return false;
        }
        
        mpfr_set_ld(llc->radius, MISS_R, MPFR_RNDN);
        
        lapse(&te, time);
        printf("done in %s.\n", time);
        
        for (int pp = 2; pp <= n; pp ++) {
            int per = n + 1 - pp;
            printf("Deriving the level curve for simple Mis(%d, %d) ... ", pp, per);
            ftime(&te);
            
            levm lm = levm_simple_from_hyp(llc, pp, per, MISS_GUARD);
            
            if(lm == NULL) {
                printf("failed !\n");
                levc_free(llc);
                
                return false;
            }
            
            lapse(&te, time);
            printf("done in %s. ", time);
            
            ftime(&te);
            if(missi_file_name(fn, 99, pp, per) && levm_write_to(lm, fn)) {
                lapse(&te, time);
                printf("Saved in %s.\n", time);
            } else {
                printf("\nCould NOT write the level set !!!\n");
            }
        }
        
        levc_free(llc);
        
        lapse(&ta, time);
        printf("All simple level sets of type %d computed and saved in %s.\n\n", tp, time);
    }
    
    lapse(&tt, time);
    printf("Total computing time: %s.\n", time);
    
    return true;
}

// MARK: IO operations

bool missi_file_name(char *fileName, int len, int pp, int per) {
    int n = pp + per;
    if(len < 30 || n < MISS_MIN_PER || n > MISS_MAX_PER) {
        return false;
    }
    
    return snprintf(fileName, len, "%s/misSimpleSet%d_%d.levm", MISS_FOLDER, pp, per) < len;
}

levm missi_load(int pp, int per) {
    if(pp + per < 3) {
        return NULL;
    }
    
    if(pp + per < MISS_MIN_PER) {
        mpfr_t r;
        mpfr_init2(r, MISS_PREC);
        mpfr_set_ld(r, MISS_R, MPFR_RNDN);
        
        levm lm = levm_simple_new(pp, per, 1, 120, r, MISS_GUARD);
        
        mpfr_clear(r);
        
        return lm;
    }
    
    char fn[100];
    
    return missi_file_name(fn, 99, pp, per) ? levm_read_from(fn) : NULL;
}


// MARK: the main function

int mis_simple_sets_main(int argc, const char * argv[]) {
    // prepare the output folder
    dir(MISS_FOLDER);
    
    return ! compute_misi_sets();
}
