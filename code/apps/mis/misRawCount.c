//
//  misRawCount.c
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
#include "mandel.h"
#include "nSet.h"
#include "io.h"
#include "misRaw.h"
#include "treeMap.h"
#include "misRawCount.h"

// MARK: constants and buffers definitions

#define MINI_BUF      1000
#define MAX_UNI_SIZE  35000000

#define NEXT_FILES_COUNT 20

static int jobs, pp, per, bufLen, loaded = 0, nuniCount = 0, actCut = 0;
static nSet_struct *n2b, *n2l, *n2r;
static unsigned long *used, totCount = 0, totRealCount = 0, timer = 0;
static nset_t n2cut, n2buf, right, nuni;
static mpc bl, br;
static tmap tr;
static ulong baseMem = 0, maxMem = 0;

static bool *jleft;
static int lastJob;

// MARK: static functions

/// @brief Checks if the given @c job is loaded in the memory buffer.
///
/// @param job the index of the raw Misiurewicz job nset file
///
/// @return @ref true if the mis raw nset file is loaded into the memory, @ref false otherwise
static bool isLoaded(int job) {
    if(job < 0 || job >= jobs) {
        return false;
    }
    
    return n2b[job].count > 0;
}

/// @brief Removes the specified job from the memory buffer.
///
/// @param job the index of the job file
static bool unload(int job) {
    if(! isLoaded(job)) {
        return false;
    }
    
    nset_clear(n2b + job);
    loaded --;
    
    return true;
}

// computes left and right boundaries of the cut
static void interval(u128 li, u128 ri, int cut, int cuts) {
    mpc_setl(bl, cut, 0);
    mpc_divi(bl, bl, cuts);
    mpc_addi(bl, bl, -2);
    u128_set(li, bl);
    
    mpc_setl(br, cut + 1, 0);
    mpc_divi(br, br, cuts);
    mpc_addi(br, br, -2);
    u128_set(ri, br);
}

static int prevList[NEXT_FILES_COUNT];
static int prevLLen = 0;

/// @brief Predicts which misRaw nset files would be needed next.
///
/// Used for optimizing networked disk access in a large, distributed computing machine.
/// Checks if the previous prediction is identical, only prints new forecasts.
///
/// @param cut the current interval or cut
/// @param cuts the total number of intervals
/// @param count the length of the files list
static void printNext(int cut, int cuts, int count) {
    if(count > NEXT_FILES_COUNT) {
        count = NEXT_FILES_COUNT;
    }
    
    int len = 0;
    int list[count];
    u128 li, ri;
    
    for (int c = cut; c < 4 * cuts && len < count; c++) {
        interval(li, ri, c, cuts);
        
        for (int i = lastJob; i >= 0 && len < count; i--) {
            // check if the file will not be loaded for the interval c
            if(jleft[i] || nset_left(ri, n2l + i) || nset_right(li, n2r + i) ||
               ! nset_right(ri, n2l + i) || ! nset_left(li, n2r + i) || isLoaded(i)) {
                continue;
            }
            
            bool found = false;
            for (int j = 0; j < len && ! found; j ++) {
                found = list[j] == i;
            }
            
            if(! found) {
                list[len ++] = i;
            }
        }
    }

    // check for reapeated list
    bool equal = prevLLen == len;
    for (int i = 0; i < len && equal; i++) {
        equal = prevList[i] == list[i];
    }
    
    if(equal) {
        return;
    }
    
    // store the current list
    prevLLen = len;
    for (int i = 0; i < len; i++) {
        prevList[i] = list[i];
    }
    
    // announce next files for external buffering
    printf("        prepare:");
    for (int i = 0; i < len; i++) {
        printf(" %d", list[i]);
    }
    printf("\n");
}

/// @brief Provides a set of points corrsponding to a .nset file. It frees the memory buffers by inverse chronological
/// order of usage, when needed.
///
/// @param job the index of the computing job (or .nset file)
///
/// @return the set of points produced by @c job
static nset getJob(int job, int cut, int cuts) {
    if(! isLoaded(job)) {    // load
        if(loaded == bufLen) {   // unload first accessed job
            int pos = -1;
            unsigned long min = -1L; // max value
            
            for (int i = 0; i < jobs; i++) {
                if(n2b[i].count > 0 && used[i] < min) {
                    pos = i;
                    min = used[i];
                }
            }
            
            nset_clear(n2b + pos);
            loaded --;
        }

        struct timeb ts;
        ftime(&ts);
        
        char fn[100], time[80];
        mis_raw_job_file_name(fn, 99, pp, per, job);
        
        printf("Loading %s ... ", fn);
        fflush(stdout);
        if(! nset_read(n2b + job, fn, false)) {
            printf("FAILED !\n");
            fflush(stdout);
            
            return NULL;
        } else {
            lapse(&ts, time);
            printf("done in %s (cut: %d)\n", time, actCut);
            fflush(stdout);

            printNext(cut, cuts, NEXT_FILES_COUNT);
        }
        
        loaded ++;
    }
    
    used[job] = ++ timer; // keeps track of usage
    
    return n2b + job;
}

/// @brief Initializes buffers and partially reads the pre-periodic points files.
///
/// @param pper the pre-period
/// @param period the period
/// @param maxBuf the max number of .nset files stored in memory at any time
///
/// @return @ref true if successfull, @ref false otherwise
static bool init(int pper, int period, int maxBuf) {
    pp = pper;
    per = period;
    bufLen = maxBuf;
    
    jobs = mis_raw_jobs_count(pp, per);
    jleft = malloc(sizeof(bool) * jobs);
    for (int i = 0; i < jobs; i++) {
        jleft[i] = false;
    }
    lastJob = jobs - 1;
    
    n2b = malloc(jobs * sizeof(nSet_struct));
    n2l = malloc(jobs * sizeof(nSet_struct));
    n2r = malloc(jobs * sizeof(nSet_struct));
    used = malloc(jobs * sizeof(long));

    mpc_init(bl, 150);
    mpc_init(br, 150);
    
    char fn[100];
    for (int i = 0; i < jobs; i++) {
        nset_init(n2b + i, MIS_RAW_SET_EPS);
        nset_init(n2l + i, MIS_RAW_SET_EPS);
        nset_init(n2r + i, MIS_RAW_SET_EPS);
        used[i] = 0;
    }
    
    for (int i = 0; i < jobs; i++) {
        mis_raw_job_file_name(fn, 99, pp, per, i);
        bool ok = nset_read_partial(n2l + i, fn, 0, 0, MINI_BUF);
        ok = ok && nset_read_partial(n2r + i, fn, 0, -MINI_BUF, MINI_BUF);
        
        if(! ok) {
            return false;
        }
    }
    
    baseMem = jobs * (8 + 2 * n2l->maxMem) + 2 * (sizeof(mpc_struct) + 7 * 24);
    
    nset_init(n2cut, MIS_RAW_SET_EPS);
    nset_init(n2buf, MIS_RAW_SET_EPS);
    nset_init(right, MIS_RAW_SET_EPS);
    nset_init(nuni, MIS_RAW_SET_EPS);
    
    tr = tmap_new(0, 0, -2, MIS_RAW_TMAP_MAX_LEVEL, MIS_RAW_TMAP_LEVEL_STEP, -2);
    
    return true;
}

static void updateMaxMem(ulong extra) {
    ulong mem = baseMem + extra;
    
    for (int i = 0; i < jobs; i++) {
        mem += n2b[i].maxMem;
    }
    
    mem += n2cut->maxMem;
    mem += n2buf->maxMem;
    mem += right->maxMem;
    mem += nuni->maxMem;
    
    maxMem = mem > maxMem ? mem : maxMem;
}

/// Frees the buffers.
static void clearNSets(void) {
    updateMaxMem(0);
    
    for (int i = 0; i < jobs; i++) {
        nset_clear(n2b + i);
        nset_clear(n2l + i);
        nset_clear(n2r + i);
    }
    
    free(n2b);
    free(n2l);
    free(n2r);
    free(used);
    free(jleft);
    
    mpc_clear(bl);
    mpc_clear(br);
    
    nset_clear(n2cut);
    nset_clear(n2buf);
    nset_clear(right);
    nset_clear(nuni);
}

/// @brief This function construct one interval of points.
///
/// @param cut the index of the interval.
/// @param cuts the number of intervals per unit of length @c 1
/// @param map @c 1 to add point to the global map, @c 0 to skip this step
/// @param unique @c 1 to reorganize all the pre-periodic points, 0 to skip this step
///
/// @return @ref true if successful, @ref false otherwise
static bool addCut(int cut, int cuts, bool map, bool unique) {
    bool ok = true;
    u128 li, ri;
    
    // clear sets
    nset_clear(n2cut);
    nset_clear(n2buf);
    
    interval(li, ri, cut, cuts);
    
    nSet_struct *tps;
    for (int i = lastJob; i >= 0 && ok; i--) {
        // if no intersection, do nothing
        if(jleft[i]) { // already checked
            continue;
        }
        
        if(nset_right(li, n2r + i)) { // mark it, will not come back
            jleft[i] = true;
            unload(i);
            
            while(jleft[lastJob] && lastJob > 0) {
                lastJob --;
            }
            
            continue;
        }
        
        if(nset_left(ri, n2l + i)) {
            continue;
        }
        
        // there is an intersection with the left part of the job i only
        if(! nset_right(ri, n2l + i)) {
            ok = ok && nset_interval(n2buf, n2l + i, li, ri);
            ok = ok && nset_union(n2cut, n2buf, 1);
            updateMaxMem(sizeof(treeMap_struct));
            
            nset_clear(n2buf);
            
            continue;
        }
        
        // there is an intersection with the right part of the job i only
        if(! nset_left(li, n2r + i)) {
            ok = ok && nset_interval(n2buf, n2r + i, li, ri);
            ok = ok && nset_union(n2cut, n2buf, 1);
            updateMaxMem(sizeof(treeMap_struct));
            
            nset_clear(n2buf);
            
            continue;
        }
        
        tps = getJob(i, cut, cuts);
        if(tps == NULL || tps->count == 0) {
            return false;
        }
        
        ok = ok && nset_interval(n2buf, tps, li, ri);
        ok = ok && nset_union(n2cut, n2buf, 1);
        updateMaxMem(sizeof(treeMap_struct));
        
        nset_clear(n2buf);
    }
    
    totCount += n2cut->count;
    totRealCount += n2cut->realCount;
    
    if(ok && right->count > 0) {
        u128 p;
        for (long i = 0; i < right->count; i++) {
            nset_point(p, right, i);
            
            if(nset_contains(n2cut, p)) {
                totCount --;
                totRealCount -= (p->y <= MIS_RAW_SET_EPS >> 1) ? 1 : 0;
            }
        }
    }
    
    ulong memTree = 0;
    if(ok && map) {
        u128 p;
        bool added = false;
        for (long i = 0; i < n2cut->count; i++) {
            nset_point(p, n2cut, i);
            
            if(! nset_contains(right, p)) {
                added = true;
                
                tmap_add(tr, p);
            }
        }
        
        if(added) {
            memTree = tmap_memory(tr);
        }
    }
    
    if(ok && unique) {
        char fn[100];
        u128 p;
        
        if(n2cut->eps < nuni->eps) {
            ok = ok && nset_set_eps(nuni, n2cut->eps);
        }
        
        for (long i = 0; i < n2cut->count; i++) {
            nset_point(p, n2cut, i);
            if(nset_contains(right, p)) {
                continue;
            }
            
            nset_add(nuni, p);
            if(nuni->count >= MAX_UNI_SIZE) {
                mis_file_name(fn, 99, pp, per, nuniCount ++);
                
                nset_lock(nuni);
                ok &= nset_write(nuni, fn);

                updateMaxMem(memTree);
                
                nset_clear(nuni);
                nset_unlock(nuni);
            }
        }
        
        if(cut == 4 * cuts - 1 && nuni->count > 0) { // write whatever is in nuni
            mis_file_name(fn, 99, pp, per, nuniCount ++);

            nset_lock(nuni);
            ok &= nset_write(nuni, fn);

            updateMaxMem(memTree);
            
            nset_clear(nuni);
            nset_unlock(nuni);
        }
    }
    
    li->x = ri->x - MIS_RAW_SET_EPS - 2;
    ok = ok && nset_interval(right, n2cut, li, ri);

    updateMaxMem(memTree);
        
    return ok;
}

/// @brief Checks if a partial count has been attempted and recovers as much data as possible from it.
///
/// Prepares the set @c right with the righmost points before the cut that is returned.
///
/// @param cuts the number of cuts per unit
/// @param period the period
///
/// @return the fist cut that has to be checked.
static int searchPartial(int cuts, int pper, int period, bool verbose) {
    int first = 0, last = period < 28 ? 1 : 1 << (period - 27), mid = (first + last) >> 1;
    nset ps = mis_new_set(false);
    char fn[100];
    mis_file_name(fn, 100, pper, period, first);
    
    // quckly test if there is any final nset file
    if(! nset_read_partial(ps, fn, 0, -1, 1)) {
        nset_free(ps);
        
        if(verbose) {
            printf("Could not read %s, starting from cut %d.\n", fn, 0);
        }
        
        return 0;
    }
    
    while(last - first > 1) {
        mis_file_name(fn, 100, pper, period, mid);
        if(nset_read_partial(ps, fn, 0, -1, 1)) {
            if(verbose) {
                printf("%s is readable.\n", fn);
            }
            
            first = mid;
        } else {
            if(verbose) {
                printf("%s is not readable.\n", fn);
            }
            
            last = mid;
        }
        
        mid = (first + last) >> 1;
    }

    if(verbose) {
        mis_file_name(fn, 100, pper, period, mid);
        
        printf("%s is the last readable file.\n", fn);
    }
    
    // here @c first is the last complete final nset file
    // count all points found in previous files
    for (int i = 0; i < first; i++) {
        mis_file_name(fn, 100, pper, period, i);
        mfile h = mfile_load_header(fn, 0, NSET_FILE_ID, NSET_HEADER_LEN);
        
        if(h == NULL) {
            mfile_free(h);
            
            if(verbose) {
                printf("Could not load header of %s ?!\n", fn);
            }
            
            return 0;
        }
        
        h->pos = 20;
        totCount += mfile_getl(h);
        totRealCount += mfile_getl(h);
        
        mfile_free(h);
    }
    
    // load the last complete nset file and analyse it, prepare the search
    mis_file_name(fn, 100, pper, period, first);
    if(! nset_read(ps, fn, true)) {
        nset_free(ps);
        totCount = 0;
        totRealCount = 0;

        if(verbose) {
            printf("Could not load %s ?!\n", fn);
        }
        
        return 0;
    }
    
    // the file is indeed correct, search the righmost complete interval in ps
    int fc = 0, lc = 4 * cuts, mc = 2 * cuts;
    u128 li, ri;
    while(lc - fc > 1) {
        interval(li, ri, mc, cuts);
        if(nset_right(ri, ps)) {
            lc = mc;
        } else {
            fc = mc;
        }
        
        mc = (fc + lc) / 2;
    }
    
    // here @c fc is included in @c ps, but fc + 1 == lc is not
    // initialize the sets @c right, @c nuni, and return the next cut
    interval(li, ri, fc, cuts);
    
    li->x = ri->x - MIS_RAW_SET_EPS - 2;
    nset_interval(right, ps, li, ri);
    
    li->x = 0;
    nset_interval(nuni, ps, li, ri);
    nset_unlock(nuni);
    nuniCount = first;
    
    nset_free(ps);
    
    totCount += nuni->count;
    totRealCount += nuni->realCount;

    if(verbose) {
        printf("All OK, starting from cut %d.\n", lc);
    }

    return lc;
}

/// @brief Counts the points (Misiurewicz parameters) of given type @c (pper,period).
///
/// The procedure is to divide the interval @c [-2,2] in @c 4*cuts equal length intervals, starting from the left, at @c -2.
/// For each interval, it is quickly determined what are the intersections with @c misRaw*.nset files, as the buffers
/// @c n2l and @c n2r store the beginning and the end of each such set stored in a ~ @c 1GB file. The common points to the
/// previous interval are excluded. Then the cardinal of the sets of points in each such interval a summed up and compared to the
/// total number of points, for which there is a formula, see countinf functions in @c mandel.h. The more intervals used, the lower
/// the maximal number of points in any interval, thus the peak memory use.
///
/// The needed files are loaded, or retrieved in memory, if they are among the last @c maxBuf used sets. The more buffers,
/// the fewer reads, but larger memory footprint.
///
/// The points map (a geometrical tree stored as a binary .tmap file) can
/// be recomputed and compared to the union of map files of each collection of points.
///
/// The points can be reorganized to avoid duplicated and in increasing order (increasing in @c x coordinate, if equal,
/// increasing in @ y coordinate), thus further verification and exploatation are much quicker. The total size of files will
/// temporarily double.
///
/// @param pper the pre-period
/// @param period the period
/// @param cuts used intervals per unit of length @c 1
/// @param maxBuf the maximal number of .nset files stored in memory at any time
/// @param map @c 1 to compute a global map and compare it to the existing local maps, @c 0 to skip this step
/// @param unique @c 1 to reorganize all the points, 0 to skip this step
///
/// @return @ref true1 if successful, @ref false otherwise
static bool count(int pper, int period, int cuts, int maxBuf, bool map, bool unique, bool testOnly) {
    printf("Counting Misiurewicz points of type (%d, %d), with ~ %d cuts (%d cuts per unit) and %d buffers.\n",
           pper, period, 5 * cuts / 2, cuts, maxBuf);
    fflush(stdout);

    char time[80];
    struct timeb ats, ts;
    ftime(&ats);
    ftime(&ts);
    
    time_stamp(time, 80, true, true);
    printf("Started at %s\n\n", time);
    fflush(stdout);
    
    if(! init(pper, period, maxBuf)) {
        printf("Could not find the misRaw.nset files !\n");
        clearNSets();
        tmap_free(tr);
        
        return false;
    }

    lapse(&ts, time);
    printf("Initialized in %s\n", time);
    fflush(stdout);
    
    int firstCut = 0;
    if(! map && unique) {
        firstCut = searchPartial(cuts, pper, period, testOnly);
        actCut = firstCut;
        
        if(firstCut == 0) {
            printf("Could not find a partially completed count, starting from scratch.\n");
        } else {
            printf("The first %d intervals have already been checked; they contain %ld points (%ld real); continuing to the right.\n",
                   firstCut, totCount, totRealCount);
        }
    }
    
    if(testOnly) {
        printf("Only a test, will stop here.\n");
        
        return firstCut > 0;
    }
    
    bool ok = true;
    for (int i = firstCut; ok && i < 4 * cuts; i++) {
        ok = ok && addCut(i, cuts, map, unique);
        if(n2cut->count > 0) {
            actCut = i;
        }
    }
    
    // upper bound, otherwise too slow
    maxMem += tmap_memory(tr);
    clearNSets();
    
    if(! ok) {
        printf("\n\nCould not complete the count !\n");
        fflush(stdout);
    } else {
        lapse(&ats, time);
                
        long tot = 2 * totCount - totRealCount;
        
        ulong hc = mandel_mis_count(pp, per);
        if(tot == hc) {
            printf("\nFound all %ld points (%ld real) in Mis(%d, %d) in %s\n",
                   tot, totRealCount, pp, per, time);
        } else {
            printf("\nThere are %ld points of Mis(%d, %d) missing (T: %ld, R: %ld), total time %s\n",
                   hc - tot, pp, per, totCount, totRealCount, time);
        }
        
        char mem[30];
        mem_size(maxMem, mem);
        printf("Peak memory usage: %s\n\n", mem);
        fflush(stdout);
        
        if(map) {
            char fileName[100];
            snprintf(fileName, 99, "%s/%02d-%02d/mis%02d-%02d.tmap", MIS_FOLDER, pp, per, pp, per);
            
            if(tmap_save(tr, fileName)) {
                printf("Points map written to %s.\n", fileName);
            } else {
                printf("Could not write the treeMap to %s.\n", fileName);
                fflush(stdout);
                
                return false;
            }
            
            printf("Loading jobs points maps ... ");
            fflush(stdout);
            
            bool lok = true;
            tmap utr, ctr;
            utr = tmap_new(0, 0, -2, MIS_RAW_TMAP_MAX_LEVEL, MIS_RAW_TMAP_LEVEL_STEP, -2);
            for (int i = 0; lok && i < jobs; i++) {
                mis_raw_map_file_name(fileName, 99, pper, period, i);
                ctr = tmap_load(fileName, 0, false, NULL);
                
                if(ctr == NULL) {
                    printf("failed to load %s\n", fileName);
                    fflush(stdout);
                    
                    lok = false;
                } else {
                    tmap_union(utr, ctr);
                    
                    tmap_free(ctr);
                }
            }
            
            if(lok) {
                printf("done.\n");
                fflush(stdout);
                
                if(tmap_leq(tr, utr)) {
                    printf("The global map is included in the union of local maps, as it should be (%ld points are repeated).\n", utr->count - tr->count);
                    fflush(stdout);
                } else {
                    printf("The global points maps do not match !!!.\n");
                    fflush(stdout);
                }
            }
            
            tmap_free(utr);
            fflush(stdout);
        }
    }
    
    tmap_free(tr);
    
    return ok;
}

// MARK: public functions

int mis_file_name(char fn[], int max, int pper, int per, int ind) {
    if(pper + per < 28) {
        return snprintf(fn, max - 1, "%s/%02d-%02d/mis%02d-%02d.nset",
                 MIS_FOLDER, pper, per, pper, per);
    }
    
    int jobs = 1 << (pper + per - 27);
    int digits = 0;
    while(jobs > 0) {
        jobs /= 10;
        digits ++;
    }
    
    char format[100];
    snprintf(format, 99, "%%s/%%02d-%%02d/mis%%02d-%%02d_%%0%dd.nset", digits);
    
    return snprintf(fn, max - 1, format, MIS_FOLDER, pper, per, pper, per, ind);
}

int mis_results_count(int pper, int period) {
    int n = pper + period;
    
    if(pper < 2 || period < 1 || n > 35) {
        return 0;
    }
    
    if(n < 28) {
        return 1;
    }
    
    // Number of nset files actually generated for types n = 28 ... 35
    
    int res28[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    int res29[] = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2};
    
    int res30[] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2, 4};
    
    int res31[] = {8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 6, 6, 4, 8};
    
    int res32[] = {16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 15, 16, 14, 15, 12, 12, 8, 16};
    
    int res33[] = {31, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 31, 32, 31, 32, 30, 32, 27, 30, 24, 24, 16, 32};
    
    int res34[] = {62, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 63, 64, 63, 63, 63, 64, 62, 64, 61, 63, 60, 63, 54, 60, 48, 48, 32, 63};
    
    int res35[] = {123, 127, 127, 127, 127, 127, 127, 127, 127, 126, 127, 126, 127, 126, 127, 126, 127, 126, 126, 126, 127, 124, 126, 122, 125, 119, 125, 107, 119, 95, 95, 64, 125};
    
    int *res[] = {res28, res29, res30, res31, res32, res33, res34, res35};
    
    return res[n - 28][pper - 2];
}

nset mis_load_results(int pper, int period, int ind) {
    nset ps = nset_new(MIS_RAW_SET_EPS, false);
    
    char fn[100];
    mis_file_name(fn, 99, pper, period, ind);
    
    if(nset_read(ps, fn, true)) {
        return ps;
    } else {
        nset_clear(ps);
        free(ps);
        
        return NULL;
    }
}

nset mis_load_partial_results(int pper, int period, int ind, long start, ulong count) {
    char fn[100];
    mis_file_name(fn, 99, pper, period, ind);
    
    FILE *f = fopen(fn, "r");
    if(f == NULL) {
        return NULL;
    }
    
    nset ps = nset_new(MIS_RAW_SET_EPS, false);
    if(nset_read_segment(ps, f, 0, start, count)) {
        fclose(f);
        
        return ps;
    } else {
        fclose(f);
        
        nset_clear(ps);
        free(ps);
        
        return NULL;
    }
}

// MARK: the help system and the main function

static const char* before = "This task counts the Misiurewicz points of given period, as follows:\n\n";
static const char* after = "\nThe more cuts used and fewer buffers, the lower the memory footprint. More buffers speed things up!\nThere needs to be enough RAM for buffers * jobs plus twice the largest cut, normally near -2.\n\n";

static const char *parameters[] = {
    "pre-period",
    "period",
    "cuts",
    "buffSize",
    "map",
    "unique"
};

static const char *types[] = {
    "required",
    "required",
    "required",
    "optional",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "",
    "",
    "3",
    "0",
    "0"
};

static const char *descriptions[] = {
    "the pre-period, integer, at least 2",
    "the period, integer, at least 1",
    "2^cuts intervals per unit, at least 1, at most 24",
    "the max number of jobs stored in the RAM at any given time, at least 2, max 100",
    "1 to create the map file of all points of this type",
    "1 to rewrite all points without repetition, in a collection of .nset files"
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

int mis_raw_count_main(int argc, const char *argv[]) {
    int pper, per, cuts, buf = 3, map = 0, unique = 0, test = 0;
    
    if(argc < 3 || sscanf(argv[0], "%d", &pper) < 1 || sscanf(argv[1], "%d", &per) < 1 ||
       sscanf(argv[2], "%d", &cuts) < 1 || pper + per < 28 || pper + per > 35 || cuts < 1 || cuts > 24) {
        help();
        
        return 1;
    }
    
    if(argc >= 4) {
        if(sscanf(argv[3], "%d", &buf) < 1 || buf < 2 || buf > 100) {
            help();
            
            return 1;
        }
    }
        
    if(argc >= 5) {
        if(sscanf(argv[4], "%d", &map) < 1) {
            help();
            
            return 1;
        }
    }
    
    if(argc >= 6) {
        if(sscanf(argv[5], "%d", &unique) < 1) {
            help();
            
            return 1;
        }
    }
    
    if(argc >= 7) {
        if(sscanf(argv[6], "%d", &test) < 1) {
            help();
            
            return 1;
        }
    }
    
    return ! count(pper, per, 1 << cuts, buf, map, unique, test);
}
