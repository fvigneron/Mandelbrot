//
//  inversePowers.c
//  Mandel
//
//  Created by MIHALACHE Nicolae on 1/9/23.
//  Copyright © 2023 MIHALACHE Nicolae. All rights reserved.
//

#include <stdio.h>

#include "nSet.h"
#include "mpv.h"
#include "hypRawCount.h"
#include "misRawCount.h"
#include "mpv.h"
#include "mandel.h"

#include "inversePowers.h"

bool invPows(char *fn, char *out, int pp, int per, int pows, ldbl err) {
    int l2e = - log2l(err) + 1;
    int prec = 16 + 2 * pows + l2e;
    prec = prec < 170 ? 170 : prec;
    int gprec = 3 * prec / 4;
    
    mpv cp = mpv_read(out, false);
    if(cp != NULL) {
        bool found = cp->prec >= prec;
        mpv_free(cp);
        
        if(found) {
            printf("%s already exists, containing results of the requested precision\n", out);
            
            return true;
        }
    }
    
    nset ns = nset_load(fn, true);
    if(ns == NULL) {
        printf("Could not load %s\n", fn);
        
        return false;
    }
    
    printf("Read %ld points from %s\n", ns->count, fn);

    mpc c, nt;
    mpc_inits(prec, c, nt, NULL);
    
    mpfr_t s, buf;
    mpfr_init2(s, prec);
    mpfr_init2(buf, prec);
    mpfr_set_zero(s, 1);
    
    mpv pc = mpv_new(prec, pows);
    for (int i = 0; i < pows; i++) {
        mpv_set(pc, i, s);
    }
    
    u128 uc;
    for (int i = 0; i < ns->count; i++) {
        nset_point(uc, ns, i);
        u128_get(c, uc);
    
        if(125 + pows * log2l(mpc_modl(c)) < l2e || 125 < l2e) { // refine
            do {
                if(pp == 0) {
                    mandel_nt(nt, c, per);
                } else {
                    mandel_mis_nt(nt, c, pp, per);
                }
                
                mpc_sub(c, c, nt);
            } while(mpc_2exp(nt) > -gprec);
        }
        
        bool real = mpfr_zero_p(c->y) || c->y->_mpfr_exp < -gprec;
        bool nSetReal = uc->y <= (ns->eps >> 1);
        if(real != nSetReal) {
            printf("**** The point ");
            mpc_print(c, 40);
            printf("\n\tis considered %sreal in the nset file %s, but it is in fact %sreal !\n",
                   nSetReal ? "" : "NOT ", fn, real ? "" : "NOT ");
        }
        
        mpc_inv(c, c);
        mpc_set(nt, c);
        for (int i = 0; i < pows; i++) {
            if(i > 0) {
                mpc_mul(nt, nt, c);
            }
            
            mpv_get(s, pc, i);
            mpfr_set(buf, nt->x, MPFR_RNDN);
            mpfr_mul_2si(buf, buf, real ? 0 : 1, MPFR_RNDN);
            mpfr_add(s, s, buf, MPFR_RNDN);
            mpv_set(pc, i, s);
        }
    }
    
    mpc_clears(c, nt, NULL);
    mpfr_clear(s);
    mpfr_clear(buf);
    
    nset_free(ns);
    
    if(out != NULL) {
        if(mpv_write(pc, out, false)) {
            printf("Powers writtten to %s\n", out);
        } else {
            printf("Could not write powers to %s\n", out);
        }
    }
    
    mpv_free(pc);
    
    return true;
}

bool pows_file_names(char *in, char *out, int len, int pp, int per, int index, int pows) {
    if(len < 50) {
        return false;
    }
    
    char tmp[len];
    
    if(pp == 0) {
        hyp_fileName(in, len, per, index);
        hyp_fileName(tmp, len, per, index);
    } else {
        mis_file_name(in, len, pp, per, index);
        mis_file_name(tmp, len, pp, per, index);
    }
    
    char *pt = strchr(tmp, '.');
    if(pt == NULL) {
        return false;
    }
    
    pt[0] = 0;
    snprintf(out, len, "%sInvPows%d.mpv", tmp, pows);
    
    return true;
}

// MARK: the help system and the main function

static const char* before = "This task computes the sums of inverse powers in one or several final .nset files.\nThe starting points are refind to the desired precision, unsing their type\nprovided in the command line arguments.\n\n";
static const char* after = "\nThe results are stored as .mpv files in the same folder.\nTo be used with the -missingRoots task.\n\n";

static const char *parameters[] = {
    "pre-period",
    "period",
    "powers",
    "max error",
    "start",
    "count"
};

static const char *types[] = {
    "required",
    "required",
    "optional",
    "optional",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "",
    "30",
    "1e-60",
    "0",
    "0"
};

static const char *descriptions[] = {
    "the pre-period, 0 for hyperbolic polynomials",
    "the period, at most 41 in the hyperbolic case, 35 - pre-period otherwise",
    "the number of negative powers to compute, at most 1000",
    "the max error for the any power of any point, from 1e-10 to 1e-200",
    "the first file to compute the powers from",
    "the number of files to analyse, 0 for all remaining files"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 6;
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

int powers_main(int argc, const char * argv[]) {
    int pp = 0, per = 3, pows = 30, st = 0, count = 0;
    ldbl err = 1e-60;
    
    bool bp = argc < 2 || sscanf(argv[0], "%d", &pp) < 1 || sscanf(argv[1], "%d", &per) < 1;
    bp = bp || pp < 0 || per > 41 || pp == 1 || (pp > 1 && pp + per > 35);
    bp = bp || (argc > 2 && (sscanf(argv[2], "%d", &pows) < 1 || pows < 1 || pows > 1000));
    bp = bp || (argc > 3 && (sscanf(argv[3], "%Lg", &err) < 1 || err < 1e-200 || err > 1e-10));
    bp = bp || (argc > 4 && (sscanf(argv[4], "%d", &st) < 1 || st < 0));
    bp = bp || (argc > 5 && (sscanf(argv[5], "%d", &count) < 1 || count < 0));
    int totCount = pp == 0 ? hyp_resultsCount(per) : mis_results_count(pp, per);
    bp = bp || st + count > totCount;
    
    if(bp) {
        help();
        
        return 1;
    }
    
    if(count == 0) {
        count = totCount - st;
    }
    
    char mfn[120], pwfn[120];
    int done = 0;
    for (int i = st; i < st + count; i++) {
        if(pows_file_names(mfn, pwfn, 120, pp, per, i, pows)) {
            done += invPows(mfn, pwfn, pp, per, pows, err) ? 1 : 0;
        }
    }
    
    if(done == 0) {
        printf("\nCould not find any final .nset file of the requsted type.\n\n");
    } else if(done == totCount) {
        printf("\nAll %d files have been found and processed.\n\n", done);
    } else {
        printf("\nOnly %d out of %d files have been found and processed.\n\n", done, totCount);
    }
    
    return 0;
}
