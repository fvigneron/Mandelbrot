//
//  levSets.c
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
 \file levSets.c
 \brief Implements the functions defined in levSets.h
 */

#include <stdlib.h>
#include <sys/timeb.h>

#include "levSets.h"
#include "stopWatch.h"
#include "io.h"

// MARK: static functions

/// @brief Effectively computes the level sets.
///
/// @return @ref true if successfull, @ref false otherwise
static bool computeLevelSets(void) {
    struct timeb te, ta;
    char time[80], fn[100];
    
    ftime(&ta);

    int stp = LEVS_2POW + 2;
    printf("Computing level set of period %d with %d bits ... ", stp, LEVC_MIN_PREC);
    ftime(&te);

    mpfr_t r;
    mpfr_init2(r, LEVS_PREC);
    mpfr_set_ld(r, LEVS_R, MPFR_RNDN);
    levc lc = levc_new(stp, LEVS_2POW, LEVC_MIN_PREC, r, LEVS_GUARD);
    
    mpfr_clear(r);
    
    if(lc == NULL) {
        printf("failed !");

        return false;
    }
    
    printf("refining to %d bits ... ", LEVS_PREC);
    levc mplc = levc_refine(lc, LEVS_2POW, LEVS_PREC);
    levc_free(lc);
    
    if(mplc == NULL) {
        printf("failed !\n");

        return false;
    }

    lapse(&te, time);
    printf("done in %s.\n", time);

    for (int per = stp + 1; per <= LEVS_MAX_PER; per++) {
        printf("Down the rays toward period %d ... ", per);
        ftime(&te);

        lc = levc_next_period(mplc, false);
        
        if(lc == NULL) {
            printf("failed !\n");

            return false;
        }

        levc_free(mplc);
        mplc = lc;
        
        lapse(&te, time);
        printf("done in %s.\n", time);

        if(per > 27) {
            ftime(&te);

            if(levs_file_name(fn, 99, per) && levc_write_to(mplc, fn)) {
                lapse(&te, time);
                printf("Saved in %s.\n", time);
            } else {
                printf("Could NOT write the level set !!!\n");
            }
        }
    }
    
    lapse(&ta, time);
    printf("All level sets computed and saved in %s.\n", time);

    return true;
}

// MARK: IO operations

bool levs_file_name(char *fileName, int len, int per) {
    if(len < 30 || per < 28 || per > LEVS_MAX_PER) {
        return false;
    }
    
    return snprintf(fileName, len - 1, "%s/levSet%d.levc", LEVS_FOLDER, per) < len - 1;
}

levc levs_load(int per) {
    if(per < 3) {
        return NULL;
    }
    
    if(per < 28) {
        mpfr_t r;
        mpfr_init2(r, 120);
        mpfr_set_ld(r, 5, MPFR_RNDN);
        
        levc lc = levc_new(per, 1, 120, r, LEVS_GUARD);
        
        mpfr_clear(r);
        
        return lc;
    }
    
    char fn[100];
    
    return levs_file_name(fn, 99, per) ? levc_read_from(fn) : NULL;
}

// MARK: the main function

int lev_sets_main(int argc, const char * argv[]) {
    // prepare the output folder
    dir(LEVS_FOLDER);
    
    return ! computeLevelSets();
}
