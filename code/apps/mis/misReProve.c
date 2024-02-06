//
//  misReProve.c
//  Mandel_v0.6
//
//  Created by MIHALACHE Nicolae on 5/11/21.
//  Copyright © 2021 MIHALACHE Nicolae. All rights reserved.
//

#include <stdio.h>

#include "ntypes.h"
#include "mpv.h"
#include "nSet.h"
#include "misRaw.h"
#include "mandel.h"
#include "misReProve.h"
#include "misRawCount.h"
#include "nset2csv.h"
#include "io.h"

nset mis_loadProofs(int pp, int per, int job) {
    char fn[200];
    snprintf(fn, 99, "%s/%02d-%02d/misRaw%02d-%02d_%d_proof.csv", MIS_FOLDER, pp, per, pp, per, job);
    
    nset ps = nset_new(MIS_RAW_SET_EPS, false);
    nset_read_csv(ps, fn);
    nset_lock(ps);
    
    return ps;
}

int mis_is_hyp(mpc c, int per) {
    fp80 c8, m1 = {-1, 0};
    mpc_get80(c8, c);
    
    if(fp80_mod(c8) < 0.125 || fp80_dist(c8, m1) < 0.125) {
        return true;
    }
    
    mpc hc;
    mpc_init(hc, mpc_prec(c));
    
    int pf = 0;
    for (int d = 3; d <= per && ! pf; d ++) {
        if(per % d != 0) {
            continue;
        }
        
        if(! mandel_root(hc, c, d, 10)) {
            continue;
        }
        
        ldbl dst = mpc_distl(hc, c);
        pf = dst < 1E-28 ? d : 0;
    }
    
    mpc_clear(hc);
    
    return pf;
}

nset mis_prove(nset pr, int pp, int per, ldbl rad, ldbl conv, int prec) {
    nset ms = nset_new(MIS_RAW_SET_EPS, false);
    mpc c;
    mpc_init(c, prec);
    int hyp = 0;
    
    for (int i = 0; i < pr->count; i++) {
        u128 p;
        nset_point(p, pr, i);
        u128_get(c, p);
        
        if(mis_is_hyp(c, per) > 0) {
            hyp ++;
            
            continue;
        }
        
        if(! mis_refine(c, c, pp, per, (int) (prec / 3.2))) {
            char *nb = stu(p);
            printf("Could not refine %s\n", nb);
            free(nb);
            
            continue;
        }
        
        if(mandel_is_mis(c, pp, per, rad) && mandel_conv_npp(c, pp, per, rad, conv, false)) {
            nset_add(ms, p);
        }
    }
    
    mpc_clear(c);
    nset_lock(ms);
    
    if(hyp > 0) {
        printf("Found and ignored %d hyperbolic centers.\n", hyp);
    }
    
    return ms;
}

bool mis_add_to_results(int pp, int per, nset all) {
    nset_lock(all);
    u128 l, r;
    fp80 m2 = {-2, 0};
    u128_setl(l, m2);
    
    nset_t uni;
    ulong tot = 0, real = 0;
    
    int nsetCount = mis_results_count(pp, per);
    for (int i = 0; i < nsetCount; i++) {
        if(i < nsetCount - 1) {
            nset right = mis_load_partial_results(pp, per, i, -1, 1);
            if(right == NULL) {
                printf("Could not load the results file Mis(%d, %d) with index %d !\n",
                       pp, per, i);
                
                return false;
            }
            
            nset_point(r, right, 0);
        } else {
            m2->x = 1;
            u128_setl(r, m2);
        }
        
        nset_init(uni, all->eps);
        bool ok = nset_interval(uni, all, l, r);
        *l = *r;
        
        if(! ok) {
            printf("Could not intersect a nset with an interval !\n");
            
            continue;
        }
        
        if(uni->count > 0) {
            nset res = mis_load_results(pp, per, i);
            if(res == NULL) {
                printf("Could not load the results file Mis(%d, %d) with index %d !\n",
                       pp, per, i);
                
                nset_free(uni);
                
                return false;
            }
            
            ulong t = res->count, rl = res->realCount;
            if(uni->eps < res->eps) {
                res->eps = uni->eps;
            }
            
            if(! nset_union(res, uni, true)) {
                printf("Could not perform the union with the results file Mis(%d, %d) with index %d !\n",
                       pp, per, i);
                
                nset_free(uni);
                
                return false;
            }
            
            tot += res->count - t;
            real += res->realCount - rl;
                        
            char fn[200];
            mis_file_name(fn, 199, pp, per, i);
            
            char cmd[1000];
            snprintf(cmd, 1000, "mv -n %s %s_old", fn, fn); // FIXME: unix only
            if(system(cmd) == -1) {
                printf("Could not execute system command: %s\n\n", cmd);
            }
            
            nset_write(res, fn);
            nset_free(res);
        }
        
        nset_clear(uni);
    }
    
    if(tot == 0) {
        printf("No new parameters in Mis(%d, %d) were founnd.\n", pp, per);
    } else {
        printf(tot == 1 ? "Added %lu parameter (%lu real) to Mis(%d, %d) results files.\n" :
                          "Added %lu parameters (%lu of them real) to Mis(%d, %d) results files.\n",
           tot, real, pp, per);
    }
    
    return true;
}

bool mis_re_prove(int pp, int per, ldbl rad, ldbl conv, int prec) {
    printf("Searching for failed proofs of points in Mis(%d, %d), reproving with (%Lg, %Lg) and precision %d :\n\n",
           pp, per, rad, conv, prec);
    
    nset all = nset_new(MIS_RAW_SET_EPS, true);
    
    int n = pp + per;
    int jobs = n < 28 ? 1 : 1 << (n - 27);
    
    bool ok = true;
    for (int i = 0; i < jobs; i++) {
        nset pr = mis_loadProofs(pp, per, i);
        
        if(pr->count == 0) {
            nset_free(pr);
            
            continue;
        }
        
        printf(pr->count == 1 ? "Loaded %lu point with failed proofs from job %d ... " :
                                "Loaded %lu points with failed proofs from job %d ... ", pr->count, i);
        nset ms = mis_prove(pr, pp, per, rad, conv, prec);
        
        nset_free(pr);
        
        if(ms == NULL) {
            nset_free(all);
            
            return false;
        }
        
        if(ms->count == 0) {
            printf("None is in Mis(%d, %d) !\n\n", pp, per);
            nset_free(ms);
            
            continue;
        }
        
        printf("%lu of them are in Mis(%d, %d).\n", ms->count, pp, per);
        ok = ok && nset_union(all, ms, true);
        
        nset_free(ms);
    }
    
    printf("\nFound %lu distinct new parameters from Mis(%d, %d).\n", all->count, pp, per);
    
    ok = ok && mis_add_to_results(pp, per, all);
    nset_free(all);
    
    return ok;
}

// MARK: the help system and the main function

static const char* before = "This task re-proves the Misiurewicz points of given period, as follows:\n\n";
static const char* after = "\nAuto-searches for failed proofs CSV files; uses the given radius for the precision\nof the parameter and the radius of convergence of the Newton method. Refines the results if\nthe precision is larger than 128. Adds the new results to the corresponding final results files.\n\n";

static const char *parameters[] = {
    "pre-period",
    "period",
    "radius",
    "conv",
    "prec"
};

static const char *types[] = {
    "required",
    "required",
    "optional",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "",
    "1E-35",  // Should match MIS_RAW_PROOF
    "1E-32",  // Should match MIS_RAW_CONV
    "128"
};

static const char *descriptions[] = {
    "the pre-period, integer, at least 2",
    "the period, integer, at least 1",
    "the error of the parameter value in [1E-25, 1E-50]",
    "the radius of convergence of the Newton method, in (radius * 10, 1E-20]",
    "the precision to use, in bits, at least 128, at most 1024"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 5;
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

int mis_re_prove_main(int argc, const char * argv[]) {
    int pper, per, prec = 128;
    ldbl rad = MIS_RAW_PROOF, conv = MIS_RAW_CONV;
    
    if(argc < 2 || sscanf(argv[0], "%d", &pper) < 1 || sscanf(argv[1], "%d", &per) < 1 ||
       pper < 2 || per < 1 || pper + per < 28 || pper + per > 35) {
        help();
        
        return 1;
    }
    
    if(argc >= 3) {
        if(sscanf(argv[2], "%Lg", &rad) < 1 || rad < 1E-50 || rad > 1E-25) {
            help();
            
            return 1;
        }
    }
    
    if(argc >= 4) {
        if(sscanf(argv[3], "%Lg", &conv) < 1 || conv < 1E-50 || conv > 1E-20 || conv <= 10 * rad) {
            help();
            
            return 1;
        }
    }
        
    if(argc >= 5) {
        if(sscanf(argv[4], "%d", &prec) < 1 || prec < 128 || prec > 1024) {
            help();
            
            return 1;
        }
    }
    
    return ! mis_re_prove(pper, per, rad, conv, prec);
}
