//
//  misSimpleQuick.c
//  Mandel
//
//  Created by MIHALACHE Nicolae on 12/24/22.
//  Copyright © 2022 MIHALACHE Nicolae. All rights reserved.
//

#include "misSimpleQuick.h"

#include <sys/timeb.h>
#include <stdlib.h>

#include "stopWatch.h"
#include "hypQuick.h"
#include "planarSet.h"
#include "nSet.h"
#include "mandel.h"
#include "fp80.h"
#include "misRaw.h"

#define SET_EPS         6E-19
#define DIV_SET_EPS     1.3E-18
#define REFINED_SET_EPS 1E-25
#define REFINE_PREC     140
#define REFINE_CONV     1E-27
#define REFINE_MAX      8
#define QUICK_LC_CONV   1E-16
#define QUICK_RT_CONV   1E-18
#define QUICK_ANG       8
#define QUICK_R       250
#define LC_ITER        20
#define MIN_ITER       25
#define MULT_ITER       2.8
#define MAX_PER        30
#define LEFT_INTERVAL   1e-10

static pset actSet, divSet;
static nset_t leftSet;

/// Initializes buffers.
///
/// @param per the period
static bool initSets(int pp, int per, bool refine) {
    if(refine) {
        nset_init(leftSet, MIS_RAW_SET_EPS);
    }
    pset_init(actSet, SET_EPS);
    pset_init(divSet, DIV_SET_EPS);
    
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
        
        if(dper > 2) {
            snprintf(fn, 99, "hyp%02d_19.csv", dper);
            ok = ok && pset_read(divSet, fn);
        }
    }
    
    pset_lock(divSet); // quicker searches
    
    return ok;
}

static bool compute(int pp, int per, bool refine) {
    int n = pp + per - 1;
    
    fp80 sp = {-2, 0}, t = {QUICK_R, 0}, left, right;
    bool ok = mandel_sol_refl(sp, NULL, sp, t, n, 1000, QUICK_LC_CONV, 1);
    ok = ok && mandel_miss_root_refl(left, sp, t, pp, per, 1000, QUICK_LC_CONV, 1, 1);
    
    sp->x = 0.25 + 8 / pow(n, 1.9);
    ok = ok && mandel_sol_refl(sp, NULL, sp, t, n, 1000, QUICK_LC_CONV, 1);
    ok = ok && mandel_miss_root_refl(right, sp, t, pp, per, 1000, QUICK_LC_CONV, 1, 1);
    
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
    
    def_mpc(REFINE_PREC, mc);
    
    long spCount = (1L << (n - 2)) * qa / sa;
    fp80 p8, c, op, z = {0, 0};
    *p8 = *left;
    int iter = n * MULT_ITER;
    iter = iter < MIN_ITER ? MIN_ITER : iter;
    ldbl md = n < 10 ? 100 : n < 20 ? 0.5 : 0.1;
    for (long i = spCount; i > 0 && ok; i--) {
        // compute a root from the starting point p8 from the level curve
        if(refine && p8->x <= -2 + LEFT_INTERVAL) {
            mpc_set80(mc, p8);
            if(mandel_miss_root_ref(mc, mc, NULL, pp, per, iter + REFINE_MAX, 0)) {
                mpfr_abs(mc->y, mc->y, MPFR_RNDN);
                mpc_get80(c, mc);
                
                // if not of lower period
                if(! pset_contains(divSet, c)) {
                    if(mpfr_cmp_ld(mc->x, -2 + LEFT_INTERVAL) <= 0) {
                        nset_put(leftSet, mc);
                    } else {
                        pset_add(actSet, c);
                    }
                }
            }
        } else if(mandel_miss_root_refl(c, p8, z, pp, per, iter, QUICK_RT_CONV, 1, md)) {
            // if not of lower period
            c->y = c->y < 0 ? -c->y : c->y;
            if(! pset_contains(divSet, c)) {
                if(refine && c->x <= -2 + LEFT_INTERVAL) {
                    mpc_set80(mc, c);
                    if(mandel_miss_root_ref(mc, mc, NULL, pp, per, REFINE_MAX, 0)) {
                        fp80 nc;
                        mpc_get80(nc, mc);
                        if(fp80_dist(c, nc) < SET_EPS) {
                            nset_put(leftSet, mc);
                        }
                    }
                } else {
                    pset_add(actSet, c);
                }
            }
        }
        
        // compute the next starting point on the level curve with @c sa intermediary points
        for (int j = sa - 1; j >= 0 && ok; j--) {
            int ang = (int) ((sa * i + j - sa + qa) % qa);
            
            *op = *p8;
            ok = ok && mandel_miss_root_refl(p8, p8, scl + ang, pp, per, LC_ITER, QUICK_LC_CONV, 1, md);
            
            if(! ok) {
                printf("Could not compute point %ld starting at %.19Lf %.19Lf\n",
                       2 * i + j, op->x, op->y);
            }
        }
    }
    
    if(! ok) {
        return false;
    }
    
    ok = fp80_dist(p8, right) < 2 * QUICK_LC_CONV;
    if(! ok) {
        printf("At the right, it did not end well, \nstart: %.19Lf %.19Lf\n  end: %.19Lf %.19Lf\n",
               op->x, op->y, p8->x, p8->y);
    }
    
    return ok;
}

static void msq_refine(pset ps, int pp, int per) {
    pset rs;
    pset_init(rs, REFINED_SET_EPS);
    fp80 p;
    ldbl eps = 2 * ps->eps;
    
    defs_mpc(REFINE_PREC, nt, c);
    for (int i = 0; i < ps->count; i++) {
        pset_point(p, ps, i);
        if(p->y >= -eps && p->y <= eps) {
            mpc_set80(c, p);
            mandel_mis_root_ref(c, c, NULL, pp, per, REFINE_MAX, 0);
            mpc_get80(p, c);
            
            p->y = p->y < 0 ? -p->y : p->y;
            if(p->y < REFINED_SET_EPS) {
                p->y = 0;
            }
        }
        
        pset_add(rs, p);
    }
    
    pset_clear(ps);
    *ps = *rs;
}

static bool findSimpleQuickMis(int st, int en, bool refine) {
    if(st < en) {
        printf("\nComputing pre-periodic points of types %d to %d (R = %.1lf, Lc Err = %.1le, Rt Err = %.1le, Set Err = %.1le):\n\n", st, en, (double) QUICK_R, QUICK_LC_CONV, QUICK_RT_CONV, SET_EPS);
        fflush(stdout);
    } else {
        printf("\nComputing pre-periodic points of type %d (R = %.1lf, Lc Err = %.1le, Rt Err = %.1le):\n\n", st, (double) QUICK_R, QUICK_LC_CONV, QUICK_RT_CONV);
        fflush(stdout);
    }
    if(refine) {
        printf("Will refine with MPFR results that are closer to the real line than twice Set Err.\n\n");
    }

    char time[80], fn[100];
    struct timeb ats, ts;
    ftime(&ats);
    
    bool ok = true, all = true;
    for (int n = st; n <= en && ok; n ++) {
        for (int per = 1; per < n - 1; per ++) {
            int pp = n - per;

            printf("Computing pre-periodic points of type (%d, %d) ... ", pp, per);
            fflush(stdout);
            ftime(&ts);
            
            // prepare sets
            ok = initSets(pp, per, refine);
            if(! ok) {
                printf("could not read the pre-periodic points of types dividing\
                       (%d, %d), will stop here.\n", pp, per);
                       fflush(stdout);
                
                return false;
            }
            
            ok = compute(pp, per, refine);
            if(! ok) {
                printf("could not compute the level curve, will stop here.\n");
                fflush(stdout);
            }
            
            if(ok) {
                nset_lock(leftSet);
                pset_lock(actSet);
                
                if(refine) {
                    msq_refine(actSet, pp, per);
                }
                
                snprintf(fn, 99, "mis%02d-%02d_19.csv", pp, per);
                ok = ! refine || leftSet->count == 0 || nset_write_csv(leftSet, fn, 26, false);
                ok = ok && pset_write(actSet, fn, refine && leftSet->count > 0);
                
                if(! ok) {
                    printf("could not write the results to %s, will stop here.\n", fn);
                    fflush(stdout);
                }
            }
            
            long rr = actSet->realCount + (refine ? leftSet->realCount : 0);
            long tr = 2 * (actSet->count + (refine ? leftSet->count : 0)) - rr;
            long roots = mandel_mis_count(pp, per);
            
            if(refine) {
                nset_clear(leftSet);
            }
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
        printf("\n");
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
static const char* after = "\nThe precision is reduced and no proofs are performed.\nThe type, that is #PP + #PER, is limited by the lack of precision to 30.\nUses reduced polynomials which have only simple roots.\n\n";

static const char *parameters[] = {
    "start",
    "end",
    "refine"
};

static const char *types[] = {
    "required",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "start",
    "0"
};

static const char *descriptions[] = {
    "the lowest type, integer, at least 3, at most 30",
    "the highest type, optional; at least start, at most 30",
    "1 to refine almost real parameters with MPFR"
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
    snprintf(format, 45, "    %%-%ds %%-%ds %%-%ds %%s\n", columnWidths[0], columnWidths[1], columnWidths[2]);
    
    printf(format, headers[0], headers[1], headers[2], headers[3]);
    printf("\n");
    for(int i = 0; i < paramCount; i++) {
        printf(format, parameters[i], types[i], defaults[i], descriptions[i]);
    }
    
    printf("%s", after);
}

int mis_simple_quick_main(int argc, const char * argv[]) {
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
    
    int ref = 0;
    if(argc >= 3 && (sscanf(argv[2], "%d", &ref) < 1)) {
        help();
        
        return 1;
    }
    
    return ! findSimpleQuickMis(start, end, ref != 0);
}
