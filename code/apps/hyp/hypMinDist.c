//
//  hypMinDist.c
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
#include "hypMinDist.h"
#include "nSet.h"
#include "hypRaw.h"
#include "hypRawCount.h"

// MARK: static functions

nset hyp_load_div(int per) {
    nset div = hyp_new_set(false), tmp = NULL;
    
    u128 c;
    fp80 c80 = {0, 0}; 
    u128_setl(c, c80);
    bool ok = nset_add(div, c);
    
    if(per % 2 == 0) {
        c80->x = -1;
        u128_setl(c, c80);
        ok = ok && nset_add(div, c);
    }
    
    ok = ok && nset_lock(div);
        
    for (int i = 3; i <= per / 2 && ok; i++) {
        if(per % i != 0) {
            continue;
        }
        
        tmp = hyp_loadResults(i, 0);
        ok = ok && tmp != NULL;
        if(! ok) {
            continue;
        }
        
        ok = ok && nset_union(div, tmp, 0);
        
        nset_clear(tmp);
        free(tmp);
    }
    
    if(! ok) {
        nset_clear(div);
        free(div);
        div = NULL;
    }
    
    return div;
}

/// Returns the minimal distance of points in the set @c ps if smaller than @c md, otherwise @c md.
///
/// @param ps the set
/// @param md the previous min distance
static ldbl minDist(nset_t ps, ldbl md) {
    if(ps == NULL) {
        return md;
    }
    
    ldbl d = nset_min_dist(ps);
    
    return d < md ? d : md;
}

/// Returns the minimal distance of pairs of points in the sets @c ls and @c rs if smaller than @c md, otherwise @c md.
///
/// @param ls the set at the left
/// @param rs the set at the right
/// @param md the previous min distance
static ldbl minDistLR(nset_t ls, nset_t rs, ldbl md) {
    if(ls == NULL || rs == NULL) {
        return md;
    }
    
    ldbl d = nset_min_dist_lr(ls, rs);
    
    return d < md ? d : md;
}

static ldbl minDistLL(nset_t ls, nset_t rs, ldbl md) {
    u128 p, q;
    bool ok = nset_point(p, ls, 0);
    ok = ok && nset_point(q, rs, 0);
    
    return ok && q->x >= p->x ? uint128_to_ldbl(q->x - p->x) : md;
}

/// Returns the left side of the interval @c I_k out of @c inter intervals.
///
/// @param k the index of the interval
/// @param inter the total number of intervals
///
/// @return the @c x coordiante of the left of the interval @c I_k
static uint128 left(int k, int inter) {
    return k >= inter - 1 ? 0 : ((uint128) 1) << (126 - 2 * (k + 1));
}

/// Returns the right side of the interval @c I_k.
///
/// @param k the index of the interval
/// @param inter the total number of intervals
///
/// @return the @c x coordiante of the right of the interval @c I_k
static uint128 right(int k, int inter) {
    return k >= inter ? 0 : ((uint128) 1) << (126 - 2 * k);
}

/// Updates the min distances in in the vector @c mds inside the main interval and inside the intervals I_k.
///
/// @param ls the set of points on the left
/// @param rs the set of points on the right
/// @param ik a buffer for the sets of points in ik
/// @param buf another buffer st
/// @param mds the vector of min distances
/// @param inter the number of small intervals near @c -2
///
/// @return @ref true if successful, @ref false otherwise
static int withoutDivisors(nset_t ls, nset_t rs, nset_t ik, nset_t buf, ldbl mds[], int inter) {
    mds[0] = minDist(rs, mds[0]);
    mds[0] = minDistLR(ls, rs, mds[0]);
    mds[2 * inter + 2] = minDistLL(ls, rs, mds[2 * inter + 2]);

    nset_clear(ik);
    
    // interval I_k, no divisors
    bool ok = true;
    for (int k = inter - 1; k >= 0 && ok; k--) {
        u128c_struct li = {.x = left(k, inter), .y = 0}, ri = {.x = right(k, inter), .y = 0};
        int ind = 2 + 2 * (inter - 1 - k);
        
        nset_move(buf, ik, 1);
        ok = ok && nset_interval(ik, rs, &li, &ri);
        
        mds[ind] = minDist(ik, mds[ind]);
        mds[ind] = minDistLR(buf, ik, mds[ind]);
    }
    
    nset_clear(buf);
    
    return ok;
}

/// Updates the min distances in in the vector @c mds inside the main interval and inside the intervals I_k,
/// taking into account the divisor periods.
///
/// @param ls the set of points on the left
/// @param rs the set of points on the right
/// @param ik a buffer for the sets of points in ik
/// @param buf another buffer st
/// @param mds the vecor of min distances
/// @param inter the number of small intervals near @c -2
///
/// @return @ref true if successful, @ref false otherwise
static bool withDivisors(nset_t ls, nset_t rs, nset_t div, nset_t ik, nset_t buf, ldbl mds[], int inter) {
    nset ld = hyp_new_set(true), rd = hyp_new_set(true);
    u128 p;
    bool ok = nset_point(p, rs, 0);
    ok = ok && nset_split(ld, rd, div, p);
    
    ok = ok && (ls == NULL || nset_union(ld, ls, 1));
    ok = ok && nset_union(rd, rs, 1);
    
    ok = ok && withoutDivisors(ld, rd, ik, buf, mds + 1, inter);
    
    nset_free(ld);
    nset_free(rd);
    
    return ok;
}

/// Computes the min dist between hyperbolic centers of period @c per and writes the results in the CSV file @c f.
///
/// @param per the period
/// @param inter the number of small intervals near @c -2
/// @param mds min distances for all intervals, of length @c 2*inter+2
///
/// @return @ref true if successful, @ref false otherwise
static bool minDistPer(int per, int inter, ldbl mds[]) {
    for (int i = 0; i <= 2 * inter + 2; i++) {
        mds[i] = HUGE_VALL;
    }
    
    // prepare the set of hyperbolic centers of period dividing per
    nset div = hyp_load_div(per);
    if(div == NULL) {
        return false;
    }
    
    // planarSet, prevPlanarSet
    nset ps = NULL, pps = NULL;
    // I_k, buffer
    nset ik = hyp_new_set(true), buf = hyp_new_set(true);
     
    int res = hyp_resultsCount(per);
    bool ok = true;
    for (int i = 0; i < res && ok; i++) {
        ps = hyp_loadResults(per, i);
        ok = ps != NULL && ps->count > 0;
        if(! ok) {
            continue;
        }
        
        ok = ok && withoutDivisors(pps, ps, ik, buf, mds, inter);
        ok = ok && withDivisors(pps, ps, div, ik, buf, mds, inter);
        
        nset_free(pps);
        pps = ps;
    }
    
    nset clr[] = {div, buf, pps, ik};
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
    char line[1000];
    int pos = snprintf(line, 50, "Period, All centers, All + div");
    for (int i = inter - 1; i >= 0; i--) {
        pos += snprintf(line + pos, 50, ", I(%d), I(%d) + div", i, i);
        
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
static bool writeCsvLine(FILE *f, int period, int inter, ldbl mds[]) {
    char line[2000];
    int pos = snprintf(line, 50, "%d", period);
    
    for (int i = 0; i <= inter; i++) {
        pos += snprintf(line + pos, 50, ", %.8Le, %.8Le", mds[2 * i], mds[1 + 2 * i]);
        
        if(pos > 1900) {
            return false;
        }
    }
    
    pos = fprintf(f, "%s\n", line);
    fflush(f);
    
    return pos > 15;
}

/// Conputes the min dist between hyperbolic centers of periods @c st to @c en.
///
/// @param st the first period
/// @param en the last period
/// @param inter the number of small intervals near @c -2
/// 
/// @return @ref true if successful, @ref false otherwise
static bool minHypDist(int st, int en, int inter) {
    printf("Computing min distance between hyperbolic centers of periods %d to %d.\n\n", st, en);
    
    char time[80];
    ptime ats, ts;
    lap(CLOCK_MONOTONIC, &ats, NULL);
    
    time_stamp(time, 80, true, true);
    printf("Started at %s\n\n", time);
    
    char fn[100];
    if(en > st) {
        snprintf(fn, 99, "%s/hypMinDist%d-%d.csv", HYP_FOLDER, st, en);
    } else {
        snprintf(fn, 99, "%s/hypMinDist%d.csv", HYP_FOLDER, st);
    }
    
    FILE *f = fopen(fn, "w");
    if(f == NULL) {
        return false;
    }
    
    bool ok = writeCsvHeader(f, inter);
    if(! ok) {
        fclose(f);
        printf("Could not write results to %s, will stop here.\n", fn);
        
        return false;
    }
    
    for (int per = st; per <= en; per++) {
        printf("Analysing centers of period %d ... ", per);
        lap(CLOCK_MONOTONIC, &ts, NULL);
        
        ldbl mds[2 * inter + 3];
        bool pok = minDistPer(per, inter, mds);
        pok = pok && writeCsvLine(f, per, inter, mds);
        
        if(pok) {
            lap(CLOCK_MONOTONIC, &ts, time);
            printf("done in %s\n", time);
        } else {
            printf("failed !\n");
        }
        
        ok = ok && pok;
        
        if(mds[0] > mds[2 * inter + 2]) {
            printf("The minimal distance %Lg is larger than the minimal width %Lg of an nset file !\n",
                   mds[0], mds[2 * inter + 2]);
        }
    }

    fclose(f);
    
    if(ok) {
        lap(CLOCK_MONOTONIC, &ats, time);
        printf("\nComputed min dist for all periods in %s\n", time);
    } else {
        printf("Analysis failed for some of the periods, please check the availability of the files.\n");
    }
    
    return ok;
}

// MARK: the help system and the main function

static const char* before = "This task computes the minimal distance between hyperbolic centers, as follows:\n\n";
static const char* after = "\nFor each considered period, one CSV line is produced. On each line, several intervals are considered, with or without\nthe centers of period which divide the current period. The interval I_k is [-2, -2 + 2^{-4k}].\n\n";

static const char *parameters[] = {
    "start",
    "end",
    "intervals"
};

static const char *types[] = {
    "required",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "start",
    "30"
};

static const char *descriptions[] = {
    "the first period period, integer, at least 1, at most 41",
    "the last period, default value: the start period",
    "the number of 2^{-4k} size intervals to consider"
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
    snprintf(format, 45, "    %%-%ds %%-%ds %%-%ds %%s\n",
             columnWidths[0], columnWidths[1], columnWidths[2]);
    
    printf(format, headers[0], headers[1], headers[2], headers[3]);
    printf("\n");
    for(int i = 0; i < paramCount; i++) {
        printf(format, parameters[i], types[i], defaults[i], descriptions[i]);
    }
    
    printf("%s", after);
}

int hyp_min_dist_main(int argc, const char *argv[]) {
    int st, en, inter = 30;
    
    if(argc < 1 || sscanf(argv[0], "%d", &st) < 1 || st < 1 || st > 41) {
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
        if(sscanf(argv[2], "%d", &inter) < 0 || inter > 40) {
            help();
            
            return 1;
        }
    }
    
    return ! minHypDist(st, en, inter);
}
