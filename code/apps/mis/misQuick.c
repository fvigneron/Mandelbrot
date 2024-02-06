//
//  misQuick.c
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

#include "misQuick.h"

#include <sys/timeb.h>
#include <stdlib.h>

#include "stopWatch.h"
#include "hypQuick.h"
#include "planarSet.h"
#include "mandel.h"
#include "fp80.h"

#define SET_EPS         1E-18
#define QUICK_LC_CONV   1E-16
#define QUICK_RT_CONV   1E-18
#define QUICK_ANG       8
#define QUICK_R       250
#define LC_ITER        20
#define MIN_ITER       25
#define MULT_ITER       2.6
#define MAX_PER        30

static pset actSet, divSet;

/// Initializes buffers.
///
/// @param per the period
static bool initSets(int pp, int per) {
    pset_init(actSet, SET_EPS);
    pset_init(divSet, SET_EPS);
    
    fp80 z = {0, 0}, m1 = {-1, 0};
    pset_add(divSet, z);
    if(per % 2 == 0) {
        pset_add(divSet, m1);
    }
    
    char fn[100];
    bool ok = true;
    for (int dper = 1; dper <= per && ok; dper++) {
        if(per % dper != 0) {
            continue;
        }

        if(dper > 2) {
            snprintf(fn, 99, "hyp%02d_19.csv", dper);
            ok = ok && pset_read(divSet, fn);
        }
        
        for (int spp = dper == per ? pp - 1 : pp; spp > 1 && ok; spp --) {
            snprintf(fn, 99, "mis%02d-%02d_19.csv", spp, dper);
            ok = ok && pset_read(divSet, fn);
        }
    }
    
    pset_lock(divSet); // quicker searches
    
    return ok;
}

static bool compute(int pp, int per) {
    int n = pp + per;
    
    fp80 sp = {-2, 0}, t = {QUICK_R, 0}, left, right;
    bool ok = mandel_sol_refl(left, NULL, sp, t, n, 1000, QUICK_LC_CONV, 1);
    
    sp->x = 0.25 + 8 / pow(n, 1.9);
    ok = ok && mandel_sol_refl(right, NULL, sp, t, n, 1000, QUICK_LC_CONV, 1);
    
    int qa = n >= 30 ? 2 * QUICK_ANG : QUICK_ANG;
    int sa = n >= 30 ? 4 : 2;
    fp80_struct scl[qa];
    
    // initialization of circle targets
    for (int i = 0; i < qa; i++) {
        scl[i].x = QUICK_R * cosl(2 * PI * i / qa);
        scl[i].y = QUICK_R * sinl(2 * PI * i / qa);
    }
    
    if(! ok) {
        return false;
    }
    
    long spCount = (1L << (n - 2)) * qa / sa;
    fp80 p8, c, op;
    *p8 = *left;
    int iter = n * MULT_ITER;
    iter = iter < MIN_ITER ? MIN_ITER : iter;
    for (long i = spCount; i > 0 && ok; i--) {
        // compute a root from the starting point p8 from the level curve
        if(ok && mandel_mis_rootl(c, p8, pp, per, iter, QUICK_RT_CONV)) {
            // if not of lower period
            c->y = c->y < 0 ? -c->y : c->y;
            c->y = c->y < 1E-20 ? 0 : c->y;
            if(! pset_contains(divSet, c)) {
                pset_add(actSet, c);
            }
        }
        
        // compute the next starting point on the level curve with @c sa intermediary points
        for (int j = sa - 1; j >= 0 && ok; j--) {
            int ang = (int) ((sa * i + j - sa + qa) % qa);
            
            *op = *p8;
            ok = ok && mandel_sol_refl(p8, NULL, p8, scl + ang, n, LC_ITER, QUICK_LC_CONV, 1);
            
            if(! ok) {
                printf("Could not compute point %ld starting at %.19Lf %.19Lf\n",
                       2 * i + j, op->x, op->y);
            }
        }
    }
    
    ok = fp80_dist(p8, right) < 2 * QUICK_LC_CONV;
    if(! ok) {
        printf("At the right, it did not end well, \nstart: %.19Lf %.19Lf\n  end: %.19Lf %.19Lf\n",
               op->x, op->y, p8->x, p8->y);
    }
    
    return ok;
}

static bool findQuickMis(int st, int en) {
    if(st < en) {
        printf("\nComputing pre-periodic points of types %d to %d (R = %.1lf, Lc Err = %.1le, \
                Rt Err = %.1le, Set Err = %.1le):\n\n",
               st, en, (double) QUICK_R, QUICK_LC_CONV, QUICK_RT_CONV, SET_EPS);
        fflush(stdout);
    } else {
        printf("\nComputing pre-periodic points of type %d (R = %.1lf, Lc Err = %.1le, Rt Err = %.1le):\n\n", st, (double) QUICK_R, QUICK_LC_CONV, QUICK_RT_CONV);
        fflush(stdout);
    }

    char time[80], fn[100];
    struct timeb ats, ts;
    ftime(&ats);
    
    bool ok = true, all = true;
    for (int n = st; n <= en && ok; n ++) {
        // TODO: restore from 1
        for (int per = n - 2; per < n - 1; per ++) {
            int pp = n - per;

            printf("Computing pre-periodic points of type (%d, %d) ... ", pp, per);
            fflush(stdout);
            ftime(&ts);
            
            // prepare sets
            ok = initSets(pp, per);
            if(! ok) {
                printf("could not read the pre-periodic points of types dividing\
                       (%d, %d), will stop here.\n", pp, per);
                       fflush(stdout);
                
                return false;
            }
            
            ok = compute(pp, per);
            if(! ok) {
                printf("could not compute the level curve, will stop here.\n");
                fflush(stdout);
            }
            
            if(ok) {
                pset_lock(actSet);
                
                snprintf(fn, 99, "mis%02d-%02d_19.csv", pp, per);
                ok = pset_write(actSet, fn, false);
                
                if(! ok) {
                    printf("could not write the results to %s, will stop here.\n", fn);
                    fflush(stdout);
                }
            }
            
            long rr = actSet->realCount;
            long tr = 2 * actSet->count - rr;
            long roots = mandel_mis_count(pp, per);
            
            pset_clear(actSet);
            pset_clear(divSet);
            
            if(ok) {
                lapse(&ts, time);
                
                if(tr == roots) {
                    printf("done in %s, found all %ld pre-periodic points (%ld of them real).\n",
                           time, tr, rr);
                    fflush(stdout);
                } else {
                    printf("done in %s, found only %ld pre-periodic points (%ld of them real) - out of %ld.\n",
                           time, tr, rr, roots);
                           fflush(stdout);
                    
                    all = false;
                }
            } else {
                printf("failed, will stop here.\n");
                fflush(stdout);
            }
        }
    }
    
    if(ok && en > st) {
        lapse(&ats, time);
        if(all) {
            printf("\nAll searches succesfully completed in %s\n", time);
            fflush(stdout);
        } else {
            printf("\nSome searches have not been successfull. Total time: %s\n", time);
            fflush(stdout);
        }
    }
    
    return ok;
}

// MARK: the help system and the main function

static const char* before = "This task computes and saves to ./mis#PP_#PER_19.txt the pre-periodic points (#PP, #PER), as follows:\n\n";
static const char* after = "\nThe precision is reduced and no proofs are performed.\nThe type, that is #PP + #PER, is limited by the lack of precision to 30.\n\n";

static const char *parameters[] = {
    "start",
    "end"
};

static const char *types[] = {
    "required",
    "optional"
};

static const char *defaults[] = {
    "",
    "start"
};

static const char *descriptions[] = {
    "the lowest type, integer, at least 3, at most 30",
    "the highest type, optional; at least start, at most 30"
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
    snprintf(format, 45, "    %%-%ds %%-%ds %%-%ds %%s\n", columnWidths[0], columnWidths[1], columnWidths[2]);
    
    printf(format, headers[0], headers[1], headers[2], headers[3]);
    printf("\n");
    for(int i = 0; i < paramCount; i++) {
        printf(format, parameters[i], types[i], defaults[i], descriptions[i]);
    }
    
    printf("%s", after);
}

int mis_quick_main(int argc, const char * argv[]) {
    int start, end;
    
    if(argc < 1 || sscanf(argv[0], "%d", &start) < 1 || start < 3 || start > MAX_PER) {
        help();
        
        return 1;
    }
    
    end = start;
    if(argc >= 2 && (sscanf(argv[1], "%d", &end) < 1 || end < start || end > MAX_PER)) {
        help();
        
        return 1;
    }
    
    return ! findQuickMis(start, end);
}
