//
//  missingRoots.c
//  Mandel
//
//  Created by MIHALACHE Nicolae on 1/9/23.
//  Copyright © 2023 MIHALACHE Nicolae. All rights reserved.
//

#include <stdio.h>

#include "mpv.h"
#include "hypRawCount.h"
#include "misRawCount.h"
#include "mpv.h"
#include "mandel.h"
#include "polynomial.h"
#include "hypRaw.h"
#include "misRaw.h"
#include "misReProve.h"
#include "nset2csv.h"
#include "misCurve.h"
#include "levCurve.h"
#include "misSets.h"
#include "levSets.h"

#include "inversePowers.h"
#include "missingRoots.h"

static void free_mpfr_ptr(mpfr_ptr v, int len) {
    for (int i = 0; i < len; i++) {
        mpfr_clear(v + i);
    }
    
    free(v);
}

static mpv write_mpv_csv(char *out, char *csv, mpfr_ptr tpows, int pows, int prec) {
    mpfr_t es;
    mpfr_init2(es, prec);
    
    mpfr_set_zero(es, 1);
    
    mpv v = mpv_new(prec, pows);
    for (int i = 0; i < pows; i++) {
        mpv_set(v, i, tpows + i);
        mpfr_abs(tpows + i, tpows + i, MPFR_RNDN);
        mpfr_add(es, es, tpows + i, MPFR_RNDU);
    }
    
    if(es->_mpfr_exp > -2 * prec / 3) {
        mpfr_printf("**** Errors sum %Rg\n", es);
    } else {
        mpfr_printf("\t Errors sum %Rg\n", es);
    }
    
    mpfr_clear(es);
    free_mpfr_ptr(tpows, pows);
    
    bool ok = mpv_write(v, out, false) > 0;
    if(ok) {
        printf("Written results to %s\n", out);
    }
    ok = ok & mpv_write_csv(v, csv, false, 1 + prec / 3.2, 0, pows, false);
    if(ok) {
        printf("Written text results to %s\n", csv);
    } else {
        printf("Some error occurred while writing results !\n");
    }
    
    printf("\n");
    
    return v;
}

mpv get_pows(char *pfn, int pp, int per, int pows) {
    mpv cp = mpv_read(pfn, false);
    if(cp != NULL) {
        printf("Read powers from %s\n", pfn);
        
        return cp;
    } else {
        printf("Could not read %s\n", pfn);
    }
    
    return NULL;
}

static void pows_files_names(char *out, char *csv, int len, int pp, int per, int pows) {
    if(pp == 0) {
        snprintf(out, len, "hyp/%02d/hyp%02d_missingPows%d.mpv", per, per, pows);
        snprintf(csv, len, "hyp/%02d/hyp%02d_missingPows%d.csv", per, per, pows);
    } else {
        snprintf(out, len, "mis/%02d-%02d/mis%02d-%02d_missingPows%d.mpv", pp, per, pp, per, pows);
        snprintf(csv, len, "mis/%02d-%02d/mis%02d-%02d_missingPows%d.csv", pp, per, pp, per, pows);
    }
}

static void roots_files_names(char *out, char *csv, int len, int pp, int per) {
    if(pp == 0) {
        snprintf(out, len, "hyp/%02d/hyp%02d_missingRoots.mpv", per, per);
        snprintf(csv, len, "hyp/%02d/hyp%02d_missingRoots.csv", per, per);
    } else {
        snprintf(out, len, "mis/%02d-%02d/mis%02d-%02d_missingRoots.mpv", pp, per, pp, per);
        snprintf(csv, len, "mis/%02d-%02d/mis%02d-%02d_missingRoots.csv", pp, per, pp, per);
    }
}

mpv compute_and_save_missing_powers(int pp, int per, int pows, ldbl err) {
    if(pp < 0 || pp == 1 || per < 1 || pp + per > 35 || pows < 1 || pows > 1000) {
        return NULL;
    }
    
    int l2e = - log2l(err) + 65;
    int prec = 16 + 2 * pows + l2e;
    prec = prec < 170 ? 170 : prec;
    
    char fn[120], pfn[120], out[120], csv[120];
    int count = 0;
    mpfr_ptr tpows;
    if(pp == 0) {
        count = hyp_resultsCount(per);
        tpows = mandel_sum_neg_pows_hyp(per, pows, prec);
    } else {
        count = mis_results_count(pp, per);
        tpows = mandel_sum_neg_pows_mis(pp, per, pows, prec);
    }
    
    mpfr_t x;
    mpfr_init2(x, prec);
    for (int i = 0; i < count; i++) {
        if(! pows_file_names(fn, pfn, 120, pp, per, i, pows)) {
            free_mpfr_ptr(tpows, pows);
            mpfr_clear(x);
            
            return NULL;
        }
        
        mpv fpows = get_pows(pfn, pp, per, pows);
        if(fpows == NULL) {
            free_mpfr_ptr(tpows, pows);
            mpfr_clear(x);
            
            return NULL;
        }
        
        for (int p = 0; p < pows; p++) {
            mpv_get(x, fpows, p);
            mpfr_sub(tpows + p, tpows + p, x, MPFR_RNDN);
        }
    }
    
    mpfr_clear(x);
    
    pows_files_names(out, csv, 120, pp, per, pows);
    mpv v = write_mpv_csv(out, csv, tpows, pows, prec); // tpows is cleared here
    
    return v;
}

static bool refine_and_prove(nset nr, int pp, int per, mpc r, int prec) {
    mpc c;
    ulong limbs[mpc_limbs(prec)];
    mpc_iniz(c, prec, limbs);
    
    if(! mis_refine(c, r, pp, per, (int) (prec / 3.2))) {
        mpfr_printf("Could not refine (%.36Rg, %.36Rg)\n", r->x, r->y);
        
        return false;
    }
    
    int hp = mis_is_hyp(c, per);
    if(hp > 0) {
        mpfr_printf("(%.40Rg, %.40Rg) is hyperbolic of period %d\n", c->x, c->y, hp);
        
        return false;
    }
    
    ldbl proof = MIS_RAW_PROOF;
    if(! mandel_is_mis(c, pp, per, proof)) {
        mpfr_printf("(%.40Rg, %.40Rg) cannot be certified with radius %Lg\n", c->x, c->y, proof);
        
        bool proven = false;
        for (int i = 0; i < 30 && ! proven; i++) {
            proof /= 2;
            proven = mandel_is_mis(c, pp, per, proof);
        }
        
        if(! proven) {
            mpfr_printf("(%.40Rg, %.40Rg) cannot be certified with radius %Lg, will give up on this one\n", c->x, c->y, proof);
            
            return false;
        }
    }
    mpfr_printf("(%.40Rg, %.40Rg) is certified with radius %Lg\n", c->x, c->y, proof);
    
    ldbl conv = MIS_RAW_CONV;
    if(! mandel_conv_npp(c, pp, per, proof, conv, false)) {
        mpfr_printf("(%.40Rg, %.40Rg) cannot be separated at scale %Lg\n", c->x, c->y, conv);
        
        bool proven = false;
        for (int i = 0; i < 30 && ! proven && conv > 3 * proof; i++) {
            conv /= 2;
            proven = mandel_conv_npp(c, pp, per, proof, conv, false);
        }
        
        if(! proven) {
            mpfr_printf("(%.40Rg, %.40Rg) cannot be cannot be separated at scale %Lg, will give up on this one\n", c->x, c->y, proof);
            
            return false;
        }
    }
    mpfr_printf("The parameter is separated at scale %Lg\n", c->x, c->y, conv);
    
    if(! nset_put(nr, c)) {
        mpfr_printf("(%.40Rg, %.40Rg) is a double root ?\n", c->x, c->y, proof);
        
        return false;
    }
    printf("It has been added to the list of new Mis(%d, %d) parameters\n", pp, per);
    
    return true;
}

static void search_levSet(levc hc, mpc r, int jobsCount) {
    long prec = mpc_prec(r);
    mpc p;
    ulong limbs[mpc_limbs(prec)];
    mpc_iniz(p, prec, limbs);
    
    ldbl md = 1e100, d;
    ulong pos = 0;
    ulong points = hc->endAngle - hc->startAngle + 1;
    for (ulong i = 0; i < points; i ++) {
        levc_point(p, hc, i);
        d = mpc_distl(p, r);
        
        if(d < md) {
            md = d;
            pos = i;
        }
    }
    
    int jobSize = (int) (levc_segs(hc) / jobsCount);
    printf("The closest point on the level curve is at position %ld, at distance %Lg and belongs to job %ld\n",
           pos, md, pos / jobSize);
}

static void search_misSet(levm mc, mpc r, int jobsCount) {
    long prec = mpc_prec(r);
    mpc p;
    ulong limbs[mpc_limbs(prec)];
    mpc_iniz(p, prec, limbs);
    
    ldbl md = 1e100, d;
    ulong pos = 0;
    ulong points = mc->endAngle - mc->startAngle + 1;
    for (ulong i = 0; i < points; i ++) {
        levm_point(p, mc, i);
        d = mpc_distl(p, r);
        
        if(d < md) {
            md = d;
            pos = i;
        }
    }
    
    int jobSize = (int) (levm_segs(mc) / jobsCount);
    printf("The closest point on the level curve is at position %ld, at distance %Lg and belongs to job %ld\n",
           pos, md, pos / jobSize);
}

bool find_inv_roots(int pp, int per, mpv ir, ldbl mod0, int prove, int add) {
    printf("Computing missing parameters from inverse powers files for type ");
    if(pp == 0) {
        printf("Hyp(%d)", per);
    } else {
        printf("Mis(%d, %d)", pp, per);
    }
    printf(":\n- coefficients of modulus smaller than %Lg are set to 0\n", mod0);
    if(prove) {
        printf("- will attempt to prove the new roots\n");
        
        if(add) {
            printf("- will attempt to add the proven parameters to the final nset files\n");
        }
    }
    printf("\n");
    
    int n = (int) ir->count;
    int prec = ir->prec;
        
    mpfr_ptr irl = malloc(sizeof(__mpfr_struct) * n);
    for (int i = 0; i < ir->count; i++) {
        mpfr_init2(irl + i, prec);
        mpv_get(irl + i, ir, i);
    }
    
    poly P;
    poly_init_root_powers_real(P, n, irl, prec, mod0);
    
    printf("Polynomial coeffs:\n");
    for (int i = 0; i <= P->deg; i++) {
        mpfr_printf("%.40Rg\n", P->a[i].x);
    }
    
    if(P->deg == 0) {
        printf("\n**** The polynomial is of degree 0 !!!\n\n");
        poly_clear(P);
        
        return false;
    }
    
    mpv roots = poly_roots(P, 4, 1e-32);
    
    mpc r;
    mpc_init(r, prec);
    printf("\nRoots of the polynomial:\n");
    for (int i = 0; i < roots->count / 2; i++) {
        mpv_getc(r, roots, i);
        mpfr_printf("%.34Rg, %.34Rg\n", r->x, r->y);
        
        mpc_inv(r, r);
        mpfr_abs(r->y, r->y, MPFR_RNDN);
        mpv_setc(roots, i, r);
    }
    
    nset_t nr;
    nset_init(nr, pp == 0 ? HYP_RAW_SET_EPS : MIS_RAW_SET_EPS);
    
    levm mc = NULL;
    levc hc = NULL;
    
    if(pp) {
        mc = miss_load(pp, per);
        
        if(mc == NULL) {
            printf("Could not load the level set !\n");
        }
    } else {
        hc = levs_load(per);
        
        if(hc == NULL) {
            printf("Could not load the level set !\n");
        }
    }
    
    int jobs = pp == 0 ? hyp_raw_jobs_count(per) : mis_raw_jobs_count(pp, per);
    printf("\nPotential missing parameters:\n");
    for (int i = 0; i < roots->count / 2; i++) {
        mpv_getc(r, roots, i);
        mpfr_printf("%.34Rg, %.34Rg\n", r->x, r->y);
        
        if(prove) {
            if(refine_and_prove(nr, pp, per, r, 180)) {
                if(pp == 0 && hc != 0) {
                    search_levSet(hc, r, jobs);
                } else {
                    search_misSet(mc, r, jobs);
                }
            }
        }
        
        printf("\n");
    }
    
    if(prove && add) {
        if(pp == 0) {
            printf("Adding new hyperbolic parameters is not implemented !\n");
        } else if(nr->count > 0) {
            printf("Found %ld new parameters, adding to the final nset files:\n", nr->count);
            
            if(! mis_add_to_results(pp, per, nr)) {
                printf("Some error occurred while adding new Misiurewicz parameters to the final results !\n");
            }
        } else {
            printf("None of the roots were of the desired type !\n");
        }
    }
    
    mpc_clear(r);
    printf("\n");
    
    char out[120], csv[120];
    roots_files_names(out, csv, 120, pp, per);
    
    bool ok = mpv_write(roots, out, false) > 0;
    if(ok) {
        printf("Written the missing roots to %s\n", out);
    }
    ok = ok & mpv_write_csv(roots, csv, true, 1 + prec / 3.2, 0, roots->count / 2, false);
    if(ok) {
        printf("Written the missing roots as text to %s\n", csv);
    } else {
        printf("Some error occurred while writing the missing roots !\n");
    }
    printf("\n");
    
    mpv_clear(roots);
    
    return ok;
}

// MARK: the help system and the main function

static const char* before = "This task computes the missing roots from the sums of their inverse powers obtained by the -inversePowers task.\n\n";
static const char* after = "\nThe results are stored as a .mpv file and exported to a .csv file in the same folder.\n\n";

static const char *parameters[] = {
    "pre-period",
    "period",
    "powers",
    "max error",
    "prove",
    "add"
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
    "1e-45",
    "0",
    "0"
};

static const char *descriptions[] = {
    "the pre-period, 0 for hyperbolic polynomials",
    "the period, at most 41 in the hyperbolic case, 35 - pre-period otherwise",
    "the number of negative powers already computed, at most 1000",
    "the max modulus of the coefficients of the resulted polynomial that should be replaced by 0",
    "1 to prove that the new roots are indeed of the desired type",
    "1 to add the new parameters to the final .nset files"
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

int missing_main(int argc, const char * argv[]) {
    int pp = 0, per = 3, pows = 30, prove = 0, add = 0;
    ldbl err = 1e-45;
    
    bool bp = argc < 2 || sscanf(argv[0], "%d", &pp) < 1 || sscanf(argv[1], "%d", &per) < 1;
    bp = bp || pp < 0 || per > 41 || pp == 1 || (pp > 1 && pp + per > 35);
    bp = bp || (argc > 2 && (sscanf(argv[2], "%d", &pows) < 1 || pows < 1 || pows > 1000));
    bp = bp || (argc > 3 && (sscanf(argv[3], "%Lg", &err) < 1));
    bp = bp || (argc > 4 && (sscanf(argv[4], "%d", &prove) < 1));
    bp = bp || (argc > 5 && (sscanf(argv[5], "%d", &add) < 1));
    
    if(bp) {
        help();
        
        return 1;
    }
    
    printf("Loading inverse powers ...\n");
    mpv mp = compute_and_save_missing_powers(pp, per, pows, err);
    
    if (mp != NULL ) {
        return find_inv_roots(pp, per, mp, err, prove, add) ? 0 : 1;
    } else {
        printf("Some files could not be read !\n");
    }
    
    return 2;
}
