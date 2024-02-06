//
//  hypSeparation.c
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

#include <stdlib.h>
#include <sys/timeb.h>

#include "stopWatch.h"
#include "nSet.h"
#include "hypRaw.h"
#include "hypRawCount.h"
#include "hypMinDist.h"
#include "mandel.h"

#include "hypSeparation.h"

ldbl min_dist(nset ps, long pos) {
    long il = pos;
    long ir = pos;
    ldbl md = 10, dx, dy, d;
    u128 p, q;
    nset_point(p, ps, pos);
    
    bool inl = true, inr = true;
    while(inl || inr) {
        il --;
        inl = inl && il >= 0;
        if(inl) {
            nset_point(q, ps, il);
            dx = uint128_to_ldbl(uint128_sub(p->x, q->x));
            if(dx > md) {
                inl = false;
            } else {
                dy = uint128_to_ldbl(uint128_dist(p->y, q->y));
                d = sqrtl(dx * dx + dy * dy);
                
                md = d < md ? d : md;
            }
        }
        
        ir ++;
        inr = inr && ir < ps->count;
        if(inr) {
            nset_point(q, ps, ir);
            dx = uint128_to_ldbl(uint128_sub(q->x, p->x));
            if(dx > md) {
                inr = false;
            } else {
                dy = uint128_to_ldbl(uint128_dist(p->y, q->y));
                d = sqrtl(dx * dx + dy * dy);
                
                md = d < md ? d : md;
            }
        }
    }
    
    return md;
}

ldbl min_ext_dist(nset ps, u128 p, ldbl maxd) {
    if(ps == NULL || ps->count == 0) {
        return maxd;
    }
    
    u128 q;
    nset_point(q, ps, 0);
    if(q->x > p->x && uint128_to_ldbl(uint128_sub(q->x, p->x)) >= maxd) {
        return maxd;
    }
    
    nset_point(q, ps, ps->count - 1);
    if(q->x < p->x && uint128_to_ldbl(uint128_sub(p->x, q->x)) >= maxd) {
        return maxd;
    }
    
    long pos = nset_search(ps, p);
    pos = pos < 0 ? -pos - 1 : pos;
    long il = pos + 1;
    long ir = pos;
    ldbl md = maxd, dx, dy, d;
    
    bool inl = true, inr = true;
    while(inl || inr) {
        il --;
        inl = inl && il >= 0;
        if(inl) {
            nset_point(q, ps, il);
            dx = uint128_to_ldbl(uint128_dist(p->x, q->x));
            if(dx > md) {
                inl = false;
            } else {
                dy = uint128_to_ldbl(uint128_dist(p->y, q->y));
                d = sqrtl(dx * dx + dy * dy);
                
                md = d < md ? d : md;
            }
        }
        
        ir ++;
        inr = inr && ir < ps->count;
        if(inr) {
            nset_point(q, ps, ir);
            dx = uint128_to_ldbl(uint128_sub(q->x, p->x));
            if(dx > md) {
                inr = false;
            } else {
                dy = uint128_to_ldbl(uint128_dist(p->y, q->y));
                d = sqrtl(dx * dx + dy * dy);
                
                md = d < md ? d : md;
            }
        }
    }
    
    return md;
}

bool addHist(ldbl v, ulong h[], int hl, bool printFirst, ldbl *minSep, mpc c) {
    int k = -log2l(fabsl(v));
    k = k < 0 ? 0 : k >= hl ? hl - 1 : k;
    
    bool printed = false;
    if(minSep != NULL && c != NULL && v < *minSep) {
        *minSep = v;
        
        mpfr_printf("\n    -min separation : %.20Lg at (%.40Rf, %.40Rf)", v, c->x, c->y);
        printed = true;
    }
    
    if(printFirst && h[k] == 0 && ! printed && c != NULL) {
        mpfr_printf("\n    -new sep inter  : %.20Lg at (%.40Rf, %.40Rf)", v, c->x, c->y);
        printed = true;
    }
    
    h[k] ++;
    
    return printed;
}

bool hist(int per, nset pps, nset ps, nset nps, nset div, ulong hs[], ulong hp[], int is, ulong hd[], int id, ldbl *minSep) {
    if(ps == NULL || ps->count <= 1 || ps->barCount != 1 || div == NULL || div->barCount != 1) {
        return false;
    }
    
    defs_mpc(128, c, v, d);
    
    bool ok = true;
    for (long i = 0; i < ps->count; i++) {
        ldbl minDist = min_dist(ps, i);
        u128 p;
        nset_point(p, ps, i);
        
        minDist = min_ext_dist(pps, p, minDist);
        minDist = min_ext_dist(nps, p, minDist);
        ldbl mdnd = minDist;
        minDist = min_ext_dist(div, p, minDist);
        
        addHist(minDist, hd, id, false, NULL, NULL);
        
        u128_get(c, p);
        mandel_val_der(v, d, c, per);
        ldbl der = mpc_modl(d);
        ldbl sep = minDist * der;
        
        bool print = addHist(sep, hs, is, true, minSep, c);
        if(mdnd == minDist) { // the closest center is not from a divisor, probably primitive
            addHist(sep, hp, is, ! print, minSep, c);
        }
    }
    
    return ok;
}

/// Computes the min dist between hyperbolic centers of period @c per and writes the results in the CSV file @c f.
///
/// @param per the period
/// @param inter the number of small intervals near @c -2
/// @param mds min distances for all intervals, of length @c 2*inter+2
///
/// @return @ref true if successful, @ref false otherwise
static bool histPer(int per, int is, ulong hs[], ulong hp[], int id, ulong hd[]) {
    for (int i = 0; i < is; i++) {
        hs[i] = 0;
    }
    for (int i = 0; i < is; i++) {
        hp[i] = 0;
    }
    for (int i = 0; i < id; i++) {
        hd[i] = 0;
    }
    
    // prepare the set of hyperbolic centers of period dividing per
    nset div = hyp_load_div(per);
    if(div == NULL) {
        return false;
    }
    
    // planarSet, prevPlanarSet, nextPlanarSet
    nset ps = NULL, pps = NULL, nps = NULL;
     
    int res = hyp_resultsCount(per);
    ldbl minSep = 10;
    bool ok = true;
    if(res == 1) {
        ps = hyp_loadResults(per, 0);
        
        ok = hist(per, pps, ps, nps, div, hs, hp, is, hd, id, &minSep);
    } else if(res == 2) {
        ps = hyp_loadResults(per, 0);
        nps = hyp_loadResults(per, 1);
        
        ok = hist(per, pps, ps, nps, div, hs, hp, is, hd, id, &minSep);
        ok = ok && hist(per, ps, nps, NULL, div, hs, hp, is, hd, id, &minSep);
    } else {
        nps = hyp_loadResults(per, 0);
        
        for (int i = 0; ok && i < res; i++) {
            nset_free(pps);
            pps = ps;
            ps = nps;
            
            nps = hyp_loadResults(per, i);
            
            ok = ok && hist(per, pps, ps, nps, div, hs, hp, is, hd, id, &minSep);
        }
    }
    
    nset clr[] = {div, ps, pps, nps};
    nset_clears(clr, 4);
        
    return ok;
}

/// Writes the fist line of the CVS output file, dynamically built for @c 2+2*inter results.
///
/// @param f the file
/// @param inter the number of small intervals
///
/// @return @ref true if successful, @ref false otherwise
static int writeCsvHeader(FILE *f, int inter) {
    char line[1500];
    int pos = snprintf(line, 50, "Period");
    for (int i = inter - 1; i >= 0; i--) {
        pos += snprintf(line + pos, 50, ", I(%d)", i);
        
        if(pos > 900) {
            return false;
        }
    }
    
    pos = fprintf(f, "%s\n", line);
    
    return pos > 15;
}

/// Writes a line of the CVS output file, corresponding to the given @c period.
///
/// @param f the file
/// @param period the period
/// @param inter the number of small intervals
/// @param mds the list of results
///
/// @return @ref true if successful, @ref false otherwise
static bool writeCsvLine(FILE *f, int period, int inter, ulong h[]) {
    char line[2000];
    int pos = snprintf(line, 50, "%d", period);
    
    for (int i = inter - 1; i >= 0; i--) {
        pos += snprintf(line + pos, 50, ", %lu", h[i]);
        
        if(pos > 1900) {
            return false;
        }
    }
    
    pos = fprintf(f, "%s\n", line);
    fflush(f);
    
    return pos > 15;
}

/// Conputes the min separation between hyperbolic centers of periods @c st to @c en.
///
/// @param st the first period
/// @param en the last period
/// @param interS the number of small intervals for separation
/// @param interD the number of small intervals for separation
///
/// @return @ref true if successful, @ref false otherwise
static bool minHypSep(int st, int en, int interS, int interD) {
    printf("Computing min separation between hyperbolic centers of periods %d to %d.\n\n", st, en);

    char time[80];
    struct timeb ats, ts;
    ftime(&ats);
    
    time_stamp(time, 80, true, true);
    printf("Started at %s\n\n", time);
    
    char fns[100], fnd[100], fnp[100];
    if(en > st) {
        snprintf(fns, 99, "%s/hypSepHist%d-%d.csv", HYP_FOLDER, st, en);
        snprintf(fnp, 99, "%s/hypSepPrimHist%d-%d.csv", HYP_FOLDER, st, en);
        snprintf(fnd, 99, "%s/hypDistHist%d-%d.csv", HYP_FOLDER, st, en);
    } else {
        snprintf(fns, 99, "%s/hypSepHist%d.csv", HYP_FOLDER, st);
        snprintf(fnp, 99, "%s/hypSepPrimHist%d.csv", HYP_FOLDER, st);
        snprintf(fnd, 99, "%s/hypDistHist%d.csv", HYP_FOLDER, st);
    }
    
    FILE *fs = fopen(fns, "w");
    bool ok = fs != NULL && writeCsvHeader(fs, interS);
    if(! ok) {
        fclose(fs);
        printf("Could not write results to %s, will stop here.\n", fns);
        
        return false;
    }
    
    FILE *fd = fopen(fnd, "w");
    ok = fd != NULL && writeCsvHeader(fd, interD);
    if(! ok) {
        fclose(fs);
        fclose(fd);
        printf("Could not write results to %s, will stop here.\n", fnd);
        
        return false;
    }
    
    FILE *fp = fopen(fnp, "w");
    ok = fp != NULL && writeCsvHeader(fp, interS);
    if(! ok) {
        fclose(fs);
        fclose(fp);
        fclose(fd);
        printf("Could not write results to %s, will stop here.\n", fnp);
        
        return false;
    }
    
    int pok;
    ulong hs[interS], hd[interD], hp[interS];
    for (int per = st; per <= en; per++) {
        printf("Analysing centers of period %d : ", per);
        
        ftime(&ts);
        
        pok = histPer(per, interS, hs, hp, interD, hd);
        
        if(pok) {
            writeCsvLine(fs, per, interS, hs);
            writeCsvLine(fp, per, interS, hp);
            writeCsvLine(fd, per, interD, hd);
            
            lapse(&ts, time);
            printf("\ndone in %s\n\n", time);
        } else {
            printf("\n\n********** failed ! **********\n\n\n");
        }
        
        ok &= pok;
    }

    fclose(fs);
    fclose(fp);
    fclose(fd);
    
    if(ok) {
        lapse(&ats, time);
        printf("\nComputed min separation for all periods in %s\n", time);
    } else {
        printf("Analysis failed for some of the periods, please check the availability of the files.\n");
    }
    
    return ok;
}

// MARK: the help system and the main function

static const char* before = "This task computes the minimal separation between hyperbolic centers, as follows:\n\n";
static const char* after = "\nFor each considered period, one CSV line is produced. On each line, several intervals are considered, with\nthe centers of period which divide the current period. The interval I_k is [2^{-k}, 2^{-k-1}).\n\n";

static const char *parameters[] = {
    "start",
    "end",
    "interSep",
    "interDist"
};

static const char *types[] = {
    "required",
    "optional",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "start",
    "30",
    "50"
};

static const char *descriptions[] = {
    "the first period period, integer, at least 3, at most 41",
    "the last period, default value: the start period",
    "the number of intervals for the separation = distance * derivative",
    "the number of intervals for the distance"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 4;
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
}

int hyp_separation_main(int argc, const char * argv[]) {
    int st, en, interS = 30, interD = 50;
    
    
    if(argc < 1 || sscanf(argv[0], "%d", &st) < 1 || st < 3 || st > 41) {
        help();
        
        return 1;
    }
    
    en = st;
    if(argc >= 2) {
        if(sscanf(argv[1], "%d", &en) < 1 || en < st || en > 41) {
            help();
            
            return 1;
        }
    }
    
    if(argc >= 3) {
        if(sscanf(argv[2], "%d", &interS) < 1 || interS < 2 || interS > 50) {
            help();
            
            return 1;
        }
    }
    
    if(argc >= 4) {
        if(sscanf(argv[3], "%d", &interD) < 1 || interD < 2 || interD > 100) {
            help();
            
            return 1;
        }
    }
    
    return ! minHypSep(st, en, interS + 1, interD + 1);
}
