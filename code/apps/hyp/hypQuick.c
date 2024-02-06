//
//  hypQuick.c
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

#include <sys/timeb.h>
#include <stdlib.h>

#include "stopWatch.h"
#include "hypQuick.h"
#include "planarSet.h"
#include "mandel.h"
#include "fp80.h"

#define SET_EPS        1E-18 // should be large than the distance between two roots, >= 1e-19
#define QUICK_LC_CONV  1E-16
#define QUICK_RT_CONV  1E-18
#define QUICK_ANG      8 // min val 8, only powers of 2
#define QUICK_R       50 // min val 5, may need to be >= 50 for periods >= 30
#define LC_ITER       20
#define MIN_ITER      25
#define MULT_ITER      1.6
#define MAX_PER       33

static pset actSet, divSet;

/// Initializes buffers.
///
/// @param per the period
static bool initSets(int per) {
    pset_init(actSet, SET_EPS);
    pset_init(divSet, SET_EPS);
    
    fp80 z = {0, 0}, m1 = {-1, 0};
    pset_add(divSet, z);
    if(per % 2 == 0) {
        pset_add(divSet, m1);
    }
    
    char fn[100];
    bool ok = true;
    for (int i = 3; i <= per / 2 && ok; i++) {
        if(per % i != 0) {
            continue;
        }

        snprintf(fn, 99, "hyp%02d_19.csv", i);
        ok = pset_read(divSet, fn);
    }
    
    pset_lock(divSet); // quicker searches
    
    return ok;
}

static ldbl dst = 0, area = 0;
static bool compute(int per) {
    fp80 sp = {-2, 0}, t = {QUICK_R, 0}, left, right;
    bool ok = mandel_sol_refl(left, NULL, sp, t, per, 1000, QUICK_LC_CONV, 1);
    
    sp->x = 0.25 + 8 / pow(per, 1.9);
    ok = ok && mandel_sol_refl(right, NULL, sp, t, per, 1000, QUICK_LC_CONV, 1);
    
    int qa = per >= 30 ? 2 * QUICK_ANG : QUICK_ANG;
    long spCount = 1L << per;
    int sa = (int) ((1L << (per - 2)) * qa / spCount);
    
    fp80_struct scl[qa];
    
    // initialization of circle targets
    for (int i = 0; i < qa; i++) {
        scl[i].x = QUICK_R * cosl(2 * PI * i / qa);
        scl[i].y = QUICK_R * sinl(2 * PI * i / qa);
    }
    
    if(! ok) {
        return false;
    }
    
    fp80 p8, c, op;
    *p8 = *left;
    int iter = per * MULT_ITER;
    iter = iter < MIN_ITER ? MIN_ITER : iter;
    ldbl testX;
    
    // main loop of roots and level curve computation
    for (long i = spCount; i > 0 && ok; i--) {
        // compute a root from the starting point p8 from the level curve
        if(ok && mandel_root_refl(c, p8, per, iter, QUICK_RT_CONV, 1, 100)) {
            c->y = c->y < 0 ? -c->y : c->y;
            
            // check if not of lower period
            if(! pset_contains(divSet, c)) {
                // minor corrections at period 33, the last frontier for ldbl precision
                if(per == 33) {
                    testX = c->x + 2;
                    
                    if(testX < 3E-18 && c->y < 3E-18) {
                        c->y = 0;
                    }
                    
                    if(testX > 0) {
                        pset_add(actSet, c);
                    }
                } else {
                    // add it to the set of roots, with unicity verification up to
                    // error @c SET_EPS on each coordinate
                    c->y = c->y < 1E-20 ? 0 : c->y;
                    
                    pset_add(actSet, c);
                }
            }
        }
        
        // compute the next starting point on the level curve with @c sa-1 intermediary points
        for (int j = sa - 1; j >= 0 && ok; j--) {
            int ang = (int) ((sa * i + j - sa + qa) % qa);
            
            *op = *p8;
            ok = ok && mandel_sol_refl(p8, NULL, p8, scl + ang, per, LC_ITER, QUICK_LC_CONV, 1);
            
            if(! ok) {
                printf("Could not compute point %ld starting at %.19Lf %.19Lf\n",
                       2 * i + j, op->x, op->y);
                       fflush(stdout);
            }
            
            dst += fp80_dist(op, p8);
            area += (p8->x - op->x) * (op->y + p8->y);
        }
    }
    
    // check if the level curve arrives at the correct point on the real line
    ok = fp80_dist(p8, right) < 2 * QUICK_LC_CONV;
    if(! ok) {
        printf("At the right, it did not end well, \nstart: %.19Lf %.19Lf\n  end: %.19Lf %.19Lf\n",
               op->x, op->y, p8->x, p8->y);
               fflush(stdout);
    }
    
    return ok;
}

static bool findQuickHyp(int st, int en) {
    if(st < en) {
        printf("\nComputing hyperbolic centers of periods %d to %d (R = %.1lf, Lc Err = %.1le, Rt Err = %.1le):\n\n", st, en, (double) QUICK_R, QUICK_LC_CONV, QUICK_RT_CONV);
        fflush(stdout);
    } else {
        printf("\nComputing hyperbolic centers of period %d (R = %.1lf, Lc Err = %.1le, Rt Err = %.1le):\n\n", st, (double) QUICK_R, QUICK_LC_CONV, QUICK_RT_CONV);
        fflush(stdout);
    }

    char time[80], fn[100];
    struct timeb ats, ts;
    ftime(&ats);
    
    bool ok = true;
    for (int per = st; per <= en && ok; per ++) {
        printf("Computing hyperbolic centers of period %2d ... ", per);
        fflush(stdout);
        ftime(&ts);
        
        // prepare sets
        ok = initSets(per);
        if(! ok) {
            printf("could not read the hyperbolic centers of periods dividing %d, will stop here.\n", per);
            fflush(stdout);
            
            return false;
        }
        
        ok = compute(per);
        if(! ok) {
            printf("could not compute the level curve, will stop here.\n");
            fflush(stdout);
        }
        
        if(ok) {
            pset_lock(actSet);
            
            snprintf(fn, 99, "hyp%02d_19.csv", per);
            ok = pset_write(actSet, fn, false);
            
            if(! ok) {
                printf("could not write the results to %s, will stop here.\n", fn);
                fflush(stdout);
            }
        }
        
        long rr = actSet->realCount;
        long tr = 2 * actSet->count - rr;
        long roots = mandel_hyp_count(per);
        pset_clear(actSet);
        pset_clear(divSet);
        
        if(ok) {
            lapse(&ts, time);
            char *all[] = {"all", "only"};
            
            printf("done in %s, found %s %ld hyperbolic centers (%ld of them real).\n",
                   time, all[tr == roots ? 0 : 1], tr, rr);
                   fflush(stdout);
        } else {
            printf("failed, will stop here.\n");
            fflush(stdout);
        }
        
        printf("The length of the level %.2Lf curve of period %d is >= %.19Lf, enclosed area ~ %.19Lf\n\n",
               (ldbl) QUICK_R, per, dst, area);
        
        dst = 0;
        area = 0;
    }
    
    if(ok && en > st) {
        lapse(&ats, time);
        printf("\nAll searches completed in %s\n", time);
        fflush(stdout);
    }
    
    return ok;
}

// MARK: the help system and the main function

static const char* before = "This task computes and saves to ./hyp#PER_19.txt the hyperbolic centers of given period #PER, as follows:\n\n";
static const char* after = "\nThe precision is reduced and no proofs are performed.\nThe period is limited by the lack of precision to 33.\n\n";

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
    "the lowest period, integer, at least 3, at most 33",
    "the highest period, optional; at least start, at most 33"
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

int hypQuickMain(int argc, const char * argv[]) {
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
    
    return ! findQuickHyp(start, end);
}
