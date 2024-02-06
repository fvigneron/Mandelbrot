//
//  misMinDist.c
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2021.
//
//  Copyright 2019 - 2021 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the Lesser GNU Public Licence v3.0
//

#include <stdlib.h>
#include <sys/timeb.h>

#include "stopWatch.h"
#include "misMinDist.h"
#include "nSet.h"
#include "misRaw.h"
#include "misRawCount.h"

// MARK: static functions


/// Writes the fist line of the CVS output file, dynamically built for @c 2+2*inter results.
///
/// @param f the file
/// @param inter the number of small intervals
///
/// @return @ref true if successful, @ref false otherwise
static int writeCsvHeader(FILE *f, int inter) {
    char line[1000];
    int pos = snprintf(line, 100, "Pre-period, Period, All centers");
    for (int i = inter - 1; i >= 0; i--) {
        pos += snprintf(line + pos, 50, ", I(%d)", i);
        
        if(pos > 900) {
            return false;
        }
    }
    
    pos = fprintf(f, "%s\n", line);
    
    return pos > 15;
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

static ldbl minDistLL(nset_t ls, nset_t rs, ldbl md) {
    u128 p, q;
    bool ok = nset_point(p, ls, 0);
    ok = ok && nset_point(q, rs, 0);
    
    return ok && q->x >= p->x ? uint128_to_ldbl(q->x - p->x) : md;
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
    mds[inter + 1] = minDistLL(ls, rs, mds[inter + 1]);

    nset_clear(ik);
    
    // interval I_k, no divisors
    bool ok = true;
    for (int k = inter - 1; k >= 0 && ok; k--) {
        u128c_struct li = {.x = left(k, inter), .y = 0}, ri = {.x = right(k, inter), .y = 0};
        int ind = inter - k;
        
        nset_move(buf, ik, 1);
        ok = ok && nset_interval(ik, rs, &li, &ri);
        
        mds[ind] = minDist(ik, mds[ind]);
        mds[ind] = minDistLR(buf, ik, mds[ind]);
    }
    
    nset_clear(buf);
    
    return ok;
}

/// Writes a line of the CVS output file, corresponding to the given @c period.
///
/// @param f the file
/// @param period the period
/// @param inter the number of small intervals
/// @param mds the list of results
///
/// @return @ref true if successful, @ref false otherwise
static bool writeCsvLine(FILE *f, int pp, int period, int inter, ldbl mds[]) {
    char line[2000];
    int pos = snprintf(line, 100, "%d, %d", pp, period);
    
    for (int i = 0; i <= inter; i++) {
        pos += snprintf(line + pos, 50, ", %.8Le", mds[i]);
        
        if(pos > 1900) {
            return false;
        }
    }
    
    pos = fprintf(f, "%s\n", line);
    fflush(f);
    
    return pos > 15;
}

/// Computes the min dist between hyperbolic centers of period @c per and writes the results in the CSV file @c f.
///
/// @param pp the pre-period
/// @param per the period
/// @param inter the number of small intervals near @c -2
/// @param mds min distances for all intervals, of length @c 2*inter+2
///
/// @return @ref true if successful, @ref false otherwise
static bool minDistPer(int pp, int per, int inter, ldbl mds[]) {
    for (int i = 0; i <= inter; i++) {
        mds[i] = HUGE_VALL;
    }
        
    // planarSet, prevPlanarSet
    nset ps = NULL, pps = NULL;
    // I_k, buffer
    nset ik = mis_new_set(true), buf = mis_new_set(true);
     
    int res = mis_results_count(pp, per);
    bool ok = true;
    for (int i = 0; i < res && ok; i++) {
        ps = mis_load_results(pp, per, i);
        ok = ps != NULL && ps->count > 0;
        if(! ok) {
            continue;
        }
        
        ok = ok && withoutDivisors(pps, ps, ik, buf, mds, inter);
        
        nset_free(pps);
        pps = ps;
    }
    
    nset clr[] = {buf, pps, ik};
    nset_clears(clr, 3);
        
    return ok;
}

/// Computes the min dist between hyperbolic centers of periods @c st to @c en.
///
/// @param st the first period
/// @param en the last period
///
/// @return @ref true if successful, @ref false otherwise
static bool minMisDist(int st, int en, int inter) {
    printf("Computing min distance between Misiurewicz points of periods %d to %d.\n\n", st, en);

    char time[80];
    ptime ats, ts;
    lap(CLOCK_MONOTONIC, &ats, NULL);

    time_stamp(time, 80, true, true);
    printf("Started at %s\n\n", time);

    char fn[100];
    if(en > st) {
        snprintf(fn, 99, "%s/misMinDist%d-%d.csv", MIS_FOLDER, st, en);
    } else {
        snprintf(fn, 99, "%s/misMinDist%d.csv", MIS_FOLDER, st);
    }

    FILE *f = fopen(fn, "w");
    if(f == NULL) {
        return false;
    }

    bool ok =  writeCsvHeader(f, inter);
    if(! ok) {
        fclose(f);
        printf("Could not write results to %s, will stop here.\n", fn);

        return false;
    }

    for (int tot = st; tot <= en; tot++) {
        for (int pp = 2; pp < tot; pp++) {
            int per = tot - pp;
            printf("Analysing Mis(%d, %d) ... ", pp, per);
            lap(CLOCK_MONOTONIC, &ts, NULL);
            
            ldbl mds[inter + 2];
            bool pok = minDistPer(pp, per, inter, mds);
            pok = pok && writeCsvLine(f, pp, per, inter, mds);
            
            if(pok) {
                lap(CLOCK_MONOTONIC, &ts, time);
                printf("done in %s\n", time);
            } else {
                printf("failed !\n");
            }
            
            ok = ok && pok;
            
            if(mds[0] > mds[inter + 1]) {
                printf("The minimal distance %Lg is larger than the minimal width %Lg of an nset file !\n",
                       mds[0], mds[inter + 1]);
            }
        }
        
        printf("\n");
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

static const char* before = "This task computes the minimal distance between Misiurewicz parameters, as follows:\n\n";
static const char* after = "\nFor each considered type, one CSV line is produced. On each line, several intervals are considered, with or without\nthe points of types which divide the current type. The interval I_k is [-2, -2 + 2^{-4k}].\n\n";

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
    "the first period period, integer, at least 3, at most 41",
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

int mis_min_dist_main(int argc, const char *argv[]) {
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
    
    return ! minMisDist(st, en, inter);
}
