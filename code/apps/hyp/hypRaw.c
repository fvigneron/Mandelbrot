//
//  hypRaw.c
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

#include "hypRaw.h"

#include <sys/timeb.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>

#include "levCurve.h"
#include "stopWatch.h"
#include "mandel.h"
#include "nSet.h"
#include "io.h"
#include "treeMap.h"
#include "computeNodeId.h"
#include "levSets.h"

// MARK: constants and buffers definitions

#define LC_ITER      15
#define LC_REF        1
#define LCL_EPS       1E-18
#define RT_MIN_ITER  20
#define MP_MIN_ITER   5
#define LC_ERR        1E-15
#define NT_ERR        1E-15
#define LC_DER2       1E30
#define LC_ERR2      (LC_ERR * LC_ERR)
#define LC_CONV2      1E-36
#define LC_PREC      120

// max distace of root from the starting point, by period, starting with period = 3
static const ldbl JUMP[] = {0.8, 0.6, 0.4, 0.3, 0.2, 0.15, 0.1, 0.09,
    0.07, 0.06, 0.05, 0.05, 0.04,  0.04, 0.03, 0.03, 0.025, 0.025,
    0.02, 0.02, 0.018, 0.017, 0.015,  0.015, 0.015, 0.014, 0.013, 0.012,
    0.011, 0.011, 0.01, 0.01, 0.01,  0.01, 0.01, 0.01, 0.01, 0.01,  0.01};

// statistics
static const char  *States[] = {"Running", "Complete", "Failed"};
static const char  *COLUMNS[] = {"Job index", "Job state",
    "Start time", "End time", "Duration", "Roots", "Real roots",
    "Start points", "Lev Newton 80", "Lev Newton", "Failed proofs", "Divisor roots",
    "Divergent", "Tot Newton 80", "Conv Newton 80", "Tot Newton", "Conv Newton",
    "Ref Newton", "Prev job", "Repeated", "Jumps"};

#define columnsLen   21
#define statsLen    (columnsLen - 2)

static ulong stats[statsLen];

// sets of roots
static nset_t prevSet, actSet, divSet, rejSet;

// mp version of circle targets
static mpc_struct sc[HYP_RAW_ANG];

// fp80 version of circle targets
static fp80_struct scl[HYP_RAW_ANG];

static mpc lcp, c, d, nt, c170, nt170, d170, c130;
static fp80 m1 = {-1, 0};
static ldbl maxDist = 0;

// folder name
static char folder[100];

// MARK: initialization functions

/// @brief Initializes buffers and creates the folders on the disk, if necessary.
///
/// @param per the period
///
/// @return @ref true if successfull, @ref false otherwise
static bool initSets(int per) {
    nset_init(prevSet, HYP_RAW_SET_EPS);
    nset_init(actSet, HYP_RAW_SET_EPS);
    nset_init(divSet, HYP_RAW_SET_EPS);
    nset_init(rejSet, HYP_RAW_SET_EPS);
    
    // prepare the output folder
    snprintf(folder, 99, "%s/%02d", HYP_FOLDER, per);
    if(! dir(HYP_FOLDER) || ! dir(folder)) {
        printf("Could not create dir %s\n", folder);
        
        return false;
    }
    
    return true;
}

/// @brief Initializes the the target points used to construct the detailed level curve.
static void initCircles(void) {
    mpfr_t th;
    mpfr_init2(th, 128);
    for (int i = 0; i < HYP_RAW_ANG; i++) {
        mpc_init(sc + i, HYP_RAW_LC_PREC);
        mpfr_set_ui(th, i, MPFR_RNDN);
        mpfr_mul_2si(th, th, -HYP_RAW_ANG_2POW, MPFR_RNDN);
        
        mpc_exp_2Pi_i(sc + i, th);
        mpc_muld(sc + i, sc + i, HYP_RAW_R);
    }
    mpfr_clear(th);
    
    // fp80 version of circle targets
    for (int i = 0; i < HYP_RAW_ANG; i++) {
        scl[i].x = HYP_RAW_R * cosl(2 * PI * i / HYP_RAW_ANG);
        scl[i].y = HYP_RAW_R * sinl(2 * PI * i / HYP_RAW_ANG);
    }
    
    mpc_init(lcp, HYP_RAW_LC_PREC);
    mpc_init(c, HYP_RAW_RT_PREC);
    mpc_init(d, HYP_RAW_RT_PREC);
    mpc_init(nt, HYP_RAW_RT_PREC);
    mpc_init(c170, HYP_RAW_REF_PREC);
    mpc_init(d170, HYP_RAW_REF_PREC);
    mpc_init(nt170, HYP_RAW_REF_PREC);
    mpc_init(c130, 130);
}

// MARK: level sets and roots computing functions

/// @brief Computes the next point to @c sp on the level curve with given angle @c ang and high precision.
///
/// @param res the result
/// @param sp the starting point
/// @param per the period
/// @param ang the angle
/// @param iter the max number of iterates of the Newton method
///
/// @return @ref true if the point has been found and stored in @c res, @c 0 otherwise
static inline bool findLc(mpc res, mpc sp, int per, int ang, int iter) {
    ulong nt = mandel_tot_nt();
    
    def_mpfr(mpc_prec(res) + MPC_EXTRA_PREC, b1);
    
    int ok = mandel_sol_ref(res, b1, sp, &sc[(ang + HYP_RAW_ANG) % HYP_RAW_ANG], per, iter, LC_REF);
    
    stats[7] += mandel_tot_nt() - nt;
    
    return ok;
}

/// Quickly checks the unicity of the root @c c and that it is not a root of a period that divides @c per.
///
/// @param c the root to check
///
/// @return @ref true if @c c may be a new root, @ref false otherwise
static inline bool quickCheck(mpc c) {
    u128 u;
    if(! u128_set(u, c) || nset_contains(divSet, u)) {
        stats[9] ++;
        
        return false;
    }
    
    if(nset_contains(actSet, u)) {
        stats[17] ++;
        
        return false;
    }
    
    if(nset_contains(prevSet, u)) {
        stats[16] ++;
        
        return false;
    }
    
    return true;
}

/// @brief Verifies the unicity of a root and attempts to prove that it is a root with error at most @c HYP_RAW_PROOF
/// and that the Newton method converges in a disk centered at @c cr with radius @c HYP_RAW_CONV to @c cr.
///
/// @param per the period
/// @param cr the root
///
/// @return @ref true if @c c is a new root, @ref false otherwise
static bool checkAndAddRoot(int per, mpc cr) {
    u128 u;
    if(! u128_set(u, cr)) {
        return false;
    }
    
    mpc_set(c, cr);
    
    if(mandel_is_hyp(c, per, HYP_RAW_PROOF) && mandel_conv_np(c, per, HYP_RAW_PROOF, HYP_RAW_CONV)) {
        nset_add(actSet, u);
        
        return true;
    }
    
    for (int i = per / 2; i >= 3; i--) {
        if(per % i != 0) {
            continue;
        }
        
        if(mandel_is_hyp(c, i, HYP_RAW_PROOF)) {
            nset_add(divSet, u);
            
            stats[9] ++;
            
            return false;
        }
    }
    
    nset_add(rejSet, u);
    stats[8] ++;
    
    return false;
}

/// @brief Iterates the Newton method with low precision, as long as the steps are at least of length @c err.
///
/// @param v the last value reached with steps larger then @c err
/// @param sp the starting point
/// @param per the period
/// @param iter the maximum number of iterates
/// @param err the longest valid step, or modulus of a Newton term
///
/// @return the number of Newton steps performed
static int accelRoot(fp80 v, fp80 sp, int per, int iter, double err) {
    fp80 nt, op = {sp->x, sp->y};
    double e2 = err * err;
    
    for (int i = 0; i < iter; i++) {
        bool ntOk = mandel_ntl(nt, op, per);
        
        stats[11] ++;
        
        if(nt->x * nt->x + nt->y * nt->y < e2 || ! ntOk) {
            *v = *op;
            
            return i;
        }
        
        op->x -= nt->x;
        op->y -= nt->y;
    }

    *v = *op;
    
    return iter;
}

/// @brief Computes the Newton term of @c p:=p_{per} at @c c, that is @c c-N_p(c)=p(c)/p'(c).
///
/// @param nt the result
/// @param c the starting point
/// @param d a buffer for the derivative
/// @param per the period
static inline void newton(mpc nt, mpc c, mpc d, int per) {
    mpc_set(nt, c);
    mpc_seti(d, 1, 0);
    
    for (long k = 1; k < per; k++) {
        mpc_mul(d, d, nt);
        mpc_scale(d, d, 1);
        mpc_addi(d, d, 1);
        
        mpc_sqr(nt, nt);
        mpc_add(nt, nt, c);
    }
    
    mpc_div(nt, nt, d);
}

/// @brief Refines a root @c sp of @c p_{per} and checks as soon as possible for unicity.
///
/// Returns @ref true if the Newton method converged to a new root, @c 0 otherwise.
///
/// @param rc the result
/// @param sp the starting point
/// @param per the period
/// @param qIter the max number of iterates with normal precision
/// @param qPrec the minimal precision of the root to validate the convergence is @c 2^{-qPrec}
/// @param rIter the maximum number of iterates with extended precision
/// @param rPrec the minimal precision of the root to validate the convergence is @c 2^{-rPrec}
///
/// @return @ref true if the Newton method converged to a new root, @c 0 otherwise
static inline bool refine(mpc rc, mpc sp, int per, int qIter, int qPrec, int rIter, int rPrec) {
    mpc_set(c, sp);
    
    int conv = 0;
    long exp = 0, prexp = 0, checked = 0;
    for (int i = 0; ! conv && i < qIter; i++) {
        prexp = exp;
        
        newton(nt, c, d, per);
        mpc_sub(c, c, nt);
        
        stats[13] ++;
        
        exp = -mpc_2exp(nt);
        conv = exp > qPrec;
        checked = checked && conv; // valid only one step, if the prediction below was correct
        
        if(i == 0) {
            continue;
        }
        
        // normally already has the precision, only not confirmed yet by the scale of nt
        if(2 * exp - prexp > qPrec && ! checked) {
            if(quickCheck(c)) {
                checked = 1;
            } else {
                return false;
            }
        }
    }
    
    if(! conv) {
        stats[10] ++;
        
        return false;
    }
    
    mpc_set(c170, c);
    
    conv = 0;
    long rexp;
    for (int i = 0; ! conv && i < rIter; i++) {
        newton(nt170, c170, d170, per);
        mpc_sub(c170, c170, nt170);
        
        stats[15] ++;
        
        rexp = -mpc_2exp(nt170);
        
        // we assume, that with enough representation precision,
        // the consecutive Newton terms ratio gets smaller in log scale
        conv = rexp > rPrec || (prexp > 0 && rexp + exp - prexp > rPrec);
    }
    
    if(! conv) {
        stats[10] ++;
        
        return false;
    }
    
    mpc_set(rc, c170);
    
    return true;
}

/// Searches for a root of @c p_{per} which, if found, is checked for unicity, refined, proven to be correct
/// and eventually added to the set of points.
///
/// @param per the period
/// @param sp the starting point
/// @param sp8 the starting point with low precision
/// @param prevSp8 the previous starting point with low precision
/// @param iter the ma number of iterates
/// @param jump2 the square of the max allowed distance between the starting point and the root
///
/// @return @ref true if a new root has been found, @ref false otherwise
static int searchRoot(int per, mpc sp, fp80 sp8, fp80 prevSp8, int iter, ldbl jump2) {
    int it = iter + MP_MIN_ITER;
    int rf = 1;
    
    fp80 c8;
    long nt80 = stats[11];
    if(fp80_dist(sp8, prevSp8) >= LC_ERR) { // try fp80 acceleration
        it -= accelRoot(c8, sp8, per, iter, NT_ERR);

        if(it <= MP_MIN_ITER) { // did not converge
            rf = 0;
            
            stats[10] ++;
        } else if(fp80_dist2(c8, sp8) >= jump2) { // jumped too far
            rf = 0;
            
            stats[18] ++;
        } else if(fp80_mod(c8) <= 0.2 && fp80_dist(c8, m1) <= 0.15) { // close to 0 or -1
            rf = 0;
            
            stats[9] ++;
        } else {
            mpc_set80(c130, c8);
        }
    } else {
        mpc_set(c130, sp);
    }
    
    // mp search and validation
    if(rf) {
        ulong nt = stats[13];
        
        if(refine(c170, c130, per, it, HYP_RAW_QCONV_2POW, 3, HYP_RAW_RCONV_2POW)) {
            if(checkAndAddRoot(per, c170)) {
                mpc_get80(c8, c170);
                ldbl d = fp80_dist(c8, sp8);

                maxDist = d > maxDist ? d : maxDist;
                
                stats[12] += stats[11] - nt80;
                stats[14] += stats[13] - nt;
                
                return true;
            }
        }
    }
    
    return false;
}

/// @brief Computes the next point on the level curve using low prcision for as long as possible, then switching
/// to high precision.
///
/// @param per the period
/// @param ep the end point
/// @param ep8 the end point in low precision
/// @param sp the starting point
/// @param sp8 the starting point with low precision
/// @param prevSp8 the previous starting point with low precision
/// @param ang the angle
///
/// @return @ref true ifsuccessfull, @ref false otherwise
static bool nextLc(int per, mpc ep, fp80 ep8, mpc sp, fp80 sp8, fp80 prevSp8, int ang) {
    int iter = 3 * per / 4;
    iter = iter < LC_ITER ? LC_ITER : iter;
    
    if(fp80_dist(sp8, prevSp8) < LC_ERR) { // consecutive points on the LC too close, use mp
        if(! findLc(ep, sp, per, ang, iter)) {
            return false;
        }
        
        mpc_get80(ep8, ep);
        
        return true;
    }
    
    // perform at least a few steps with fp80, then hadover to mp
    fp80 nt, p8 = {sp8->x, sp8->y}, op8;
    fp80_ptr t = &scl[(ang + HYP_RAW_ANG) % HYP_RAW_ANG];
    
    ldbl d2, m2;
    int it = 0;
    do {
        *op8 = *p8;
        d2 = mandel_nt_soll(nt, p8, t, per);
        fp80_sub(p8, p8, nt);
        m2 = fp80_mod2(nt);
        
        it ++;
        
        stats[6] ++;
    } while(it < iter && d2 < LC_DER2 && m2 >= LC_ERR2);
    
    if(it > 1 && d2 < LC_DER2 && m2 < LC_ERR2) { // it may converge in fp80, let's see
        d2 = mandel_nt_soll(nt, p8, t, per);
        m2 = fp80_mod2(nt);
        
        stats[6] ++;
        
        if(d2 < LC_DER2 && m2 < LC_CONV2) { // yes, keep it
            fp80_sub(ep8, p8, nt);
            mpc_set80(ep, ep8);
            
            return true;
        }
    }
    
    if(it == 1) { // if cannot start in fp80, do everything in mp
        if(! findLc(ep, sp, per, ang, iter)) {
            return false;
        }
        
        mpc_get80(ep8, ep);
        
        return true;
    }
    
    // otherwise, continue from op8
    mpc_set80(ep, op8);
    if(! findLc(ep, ep, per, ang, iter - it + 1)) {
        return 0;
    }
    
    mpc_get80(ep8, ep);
    
    return true;
}

/// @brief Searches the roots using starting points from an arc of a level curve that it generates from the
/// starting point towards the enpoint.
///
/// If per <= 27, sp and ep are in R+ and respectively R-. Starting with period 28, there are 2^14 such points and
/// @c sp and @c ep are consecutive in this list.
///
/// @param per the period
/// @param sp the starting point
/// @param ep the end point
/// @param spCount the number of starting points between @c sp and @c ep
///
/// @return @ref true if successfull, @ref false otherwise
static bool roots(int per, mpc sp, mpc ep, int spCount) {
    mpc_set(lcp, sp);
    fp80 p8, op8 = {0, 0}, tp8 = {0, 0};
    mpc_get80(p8, lcp);
    
    ldbl jump2 = JUMP[per - 3] * JUMP[per - 3];
    int iter = per < RT_MIN_ITER ? RT_MIN_ITER : per;
    bool ok = true;
    
    for (int i = 1; i <= spCount && ok; i++) {
        for (int j = 0; j < 2 && ok; j++) {
            *tp8 = *p8;
            ok = ok && nextLc(per, lcp, p8, lcp, p8, op8, 2 * i + j - 1);
            *op8 = *tp8;
        }
        
        if(ok) {
            searchRoot(per, lcp, p8, op8, iter, jump2);
        }
    }
    
    // refine lcp, it may come from fp80 computation
    
    // allocate a buffer on the stack
    def_mpfr(mpc_prec(lcp) + MPC_EXTRA_PREC, b1);
    
    ok = ok && mandel_sol_ref(lcp, b1, lcp, sc, per, LC_ITER, LC_REF);
    ok = ok && mpc_distl(ep, lcp) < 1E-30;
    
    return ok;
}

/// @brief Frees the buffers.
static void clear(void) {
    nset_clear(prevSet);
    nset_clear(actSet);
    nset_clear(divSet);
    nset_clear(rejSet);
    
    for (int i = 0; i < HYP_RAW_ANG; i++) {
        mpc_clear(sc + i);
    }
    
    mpc_clear(lcp);
    mpc_clear(c);
    mpc_clear(d);
    mpc_clear(nt);
    mpc_clear(c170);
    mpc_clear(d170);
    mpc_clear(nt170);
    mpc_clear(c130);
}

/// @brief Computes and writes the tree map of the set of points @c ps.
///
/// @param fileName the file name
/// @param ps the set of point
/// @param check @ref true to check that the file reads correctly, @c 0 to ignore this step
///
/// @return @ref true if successfull, @ref false otherwise
static bool writeETree(char *fileName, nset_t ps, int check) {
    tmap tr = tmap_map(ps, HYP_RAW_TMAP_MAX_LEVEL, HYP_RAW_TMAP_LEVEL_STEP);
    
    if(tr == NULL) {
        printf("Could not construct the eTree.\n");
    } else {
        if(tmap_save(tr, fileName)) {
            printf("Points map written to %s.\n", fileName);
        } else {
            printf("Could not write the eTree to %s.\n", fileName);
            
            return false;
        }
    }
    
    bool ok = true;
    if(check) {
        tmap tc = tmap_load(fileName, 0, true, NULL);
        
        if(tc == NULL) {
            printf("Could not read %s.\n", fileName);
            
            ok = false;
        } else {
            if(tmap_eq(tr, tc)) {
//                printf("The map file is correct.\n");
            } else {
                printf("The map file is DAMAGED !!!\n");
                
                ok = false;
            }
        }
        
        tmap_free(tc);
    }
    
    tmap_free(tr);
    
    return true;
}

int hyp_raw_job_file_name(char fn[], int max, int per, int job) {
    if(per < 28) {
        return snprintf(fn, max - 1, "%s/%02d/hyp%02d.nset", HYP_FOLDER, per, per);
    } else {
        return snprintf(fn, max - 1, "%s/%02d/hypRaw%02d_%d.nset", HYP_FOLDER, per, per, job);
    }
}

int hyp_raw_map_file_name(char fn[], int max, int per, int job) {
    if(per < 28) {
        return snprintf(fn, max - 1, "%s/%02d/hyp%02d.tmap", HYP_FOLDER, per, per);
    } else {
        if(job >= (1 << (per - 27))) {
            return snprintf(fn, max - 1, "%s/%02d/hypRaw%02d.tmap", HYP_FOLDER, per, per);
        }
        
        return snprintf(fn, max - 1, "%s/%02d/hypRaw%02d_%d.tmap", HYP_FOLDER, per, per, job);
    }
}

/// @brief Provides the name of the stats file for the job with index  @c job and of period @c per.
///
/// @param fn the file name
/// @param max the max length of the file name
/// @param per the period
/// @param job the index of the file (@c 0 for periods samller then @c 28)
///
/// @return the number of characters put into @c fn
static int raw_statsFileName(char fn[], int max, int per, int job) {
    if(per < 28) {
        return snprintf(fn, max - 1, "%s/%02d/hyp%02d_stats.csv", HYP_FOLDER, per, per);
    } else {
        return snprintf(fn, max - 1, "%s/%02d/hypRaw%02d_%d_stats.csv", HYP_FOLDER, per, per, job);
    }
}

/// @brief Called to mark this job with state @c state in the stats file.
///
/// @param per the period
/// @param job the job
static void writeStats(int per, int job, int state) {
    char fn[100];
    raw_statsFileName(fn, 99, per, job);
    
    FILE *f = fopen(fn, "w");
    if(f == NULL) {
        return;
    }
    
    for (int i = 0; i < columnsLen; i++) {
        fprintf(f, "%s", COLUMNS[i]);
        fprintf(f, i < columnsLen - 1 ? ", " : "\n");
    }
    
    fprintf(f, "%d, %s, ", job, States[state]);
    
    for (int i = 0; i < statsLen; i++) {
        fprintf(f, "%lu", stats[i]);
        fprintf(f, i < statsLen - 1 ? ", " : "\n");
    }
    
    fclose(f);
}

/// @brief Called to mark this job as @c Running in the stats file.
///
/// @param per the period
/// @param job the job
static void jobStarted(int per, int job, struct timeb *ts) {
    ftime(ts);
    char date[65];
    time_stamp(date, 65, true, false);
    
    printf("\nJob %d started at %s\n", job, date);
    fflush(stdout);
    
    stats[0] = ts->time;
    stats[0] = stats[0] * 1000 + ts->millitm;      // start time, see COLUMS[2 ... 17]
    for (int i = 1; i < columnsLen - 2; i++) {
        stats[i] = 0;
    }
    
    writeStats(per, job, 0);
}

/// @brief Called to mark this job as @c Completed in the stats file.
///
/// @param per the period
/// @param job the job
static void jobDone(int per, int job, struct timeb *ts) {
    char time[80];
    
    lapse(ts, time);
    printf("Job %d completed, %ld roots found in %s.\n", job, actSet->count, time);
    fflush(stdout);
    
    struct timeb tts;
    ftime(&tts);
    
    stats[1] = tts.time;
    stats[1] = stats[1] * 1000 + tts.millitm;      // end time, see COLUMS[2 ... 17]
    stats[2] = stats[1] - stats[0];  // duration
    stats[3] = actSet->count;
    stats[4] = actSet->realCount;
    stats[5] = per < 27 ? 1L << per : 1L << 27;
    
    writeStats(per, job, 1);
}

/// @brief Called to mark this job as @c Failed in the stats file.
///
/// @param per the period
/// @param job the job
static void jobFailed(int per, int job) {
    writeStats(per, job, 2);
}

int hyp_raw_jobs_count(int per) {
    return per < 28 || per > 50 ? 1 : 1 << (per - 27);
}

nSet_struct *hyp_raw_load_job(int per, int job) {
    nSet_struct *ps = malloc(sizeof(nSet_struct));
    nset_init(ps, HYP_RAW_SET_EPS);
    
    char fn[100];
    hyp_raw_job_file_name(fn, 99, per, job);
    
    if(nset_read(ps, fn, false)) {
        return ps;
    } else {
        nset_clear(ps);
        free(ps);
        
        return NULL;
    }
}

/// @brief Searches for the roots of the polynomial @c p_{per} and saves the results.
///
/// Constructs a tree
/// map of the points found. Starting with period @c 28, the serach is divided in @c 2^{per-27} jobs,
/// for parallelization.
///
/// The names of the corresponding files are provided by @c raw_jobFileName() and the results
/// can be loaded with @c raw_loadJob(). The map file name is given by @c raw_mapFileName().
///
/// @param per the period
/// @param start the first job
/// @param end last job @c +1
///
/// @return @ref true if successfull, @ref false otherwise
static bool findHypRaw(int per, int start, int end) {
    char hostname[_SC_HOST_NAME_MAX + 1];
    gethostname(hostname, _SC_HOST_NAME_MAX + 1);
    
    int cpu = get_cpu_id();
    
    // prepare sets
    initSets(per);
    
    if(per == 1) {
        u128 p;
        fp80 c = {0, 0};
        u128_setl(p, c);
        
        nset_add(actSet, p);
    }
    
    if(per == 2) {
        u128 p;
        fp80 c = {-1, 0};
        u128_setl(p, c);
        
        nset_add(actSet, p);
    }
    
    if(per < 3) {
        printf("Creating the files for hyperbolic centers of period %d:\n\n", per);
        
        char fn[100];
        
        bool ok = true;
        
        hyp_raw_job_file_name(fn, 100, per, 0);
        nset_lock(actSet);
        if(! nset_write(actSet, fn)) {
            printf("Could not write the points to %s.\n", fn);
            ok = false;
        } else {
            printf("Written the point to %s\n", fn);
        }
        
        hyp_raw_map_file_name(fn, 100, per, 0);
        if(! writeETree(fn, actSet, 1)) {
            printf("Could not write the tree to %s.\n", fn);
            ok = false;
        }
        
        snprintf(fn, 99, "%s/hyp%02d_%d.csv", folder, per, HYP_RAW_DIGITS);
        if(nset_write_csv(actSet, fn, HYP_RAW_DIGITS, 0)) {
            printf("Exported to %s.\n\n", fn);
        } else {
            printf("Could not export the list to %s.\n\n", fn);
            ok = false;
        }
        
        return ok;
    }
    
    printf("\nComputing centers of period %d, jobs %d to %d (max jump %.4Lf) on %s with cpu %d:\n\n",
           per, start, end, JUMP[per - 3], hostname, cpu);
    
    fflush(stdout); // Flushing the buffer from time to time ensures that the logs survive an app kill
    
    struct timeb ats, ts;
    ftime(&ats);
    
    levc cc = levs_load(per);
    int jobs = per < 28 ? 1 : 1 << (per - 27);
    
    if(cc == NULL || levc_segs(cc) % jobs != 0) {
        printf("Could not prepare the level curve, will stop here.\n\n");
        fflush(stdout);
        
        return 0;
    }
    
    int jobSize = (int) (levc_segs(cc) / jobs);
    levc nc = levc_sub_curve(cc, jobSize * start, 1, jobSize * (end - start) + 1);
    levc_free(cc);
    cc = nc;
    
    char time[80], fn[100];
    lapse(&ats, time);
    
    printf("Level curve prepared in %s.\n", time);
    fflush(stdout);
    
    // init level curve targets
    initCircles(); 
    
    // try to load previous job
    if(start > 0) {
        ftime(&ts);
        
        hyp_raw_job_file_name(fn, 100, per, start - 1);
        if(nset_read(prevSet, fn, false)) {
            lapse(&ts, time);
            
            printf("Results of job %d loaded in %s\n", start - 1, time);
            fflush(stdout);
        } // otherwise, no problem, just a few extra points !
    }
    
    bool ok = true;
    long troots = 0;
    mpc sp, ep;
    mpc_init(sp, cc->prec);
    mpc_init(ep, cc->prec);
    
    for (int job = start; ok && job < end; job++) {
        jobStarted(per, job, &ts);
        
        int spCount = per < 28 ? 1 << per : 1 << (per - 14); // ~4 starting points per root
        for (int seg = 0; ok && seg < jobSize; seg++) {
            int pos = jobSize * (job - start) + seg;
            levc_point(sp, cc, pos);
            levc_point(ep, cc, pos + 1);
            
            ok = roots(per, sp, ep, spCount);
            
            if(! ok) {
                printf("Could not compute the roots from the segment %d of the level curve, will stop here.\n\n", pos);
            }
        }
        
        if(ok) {
            if(rejSet->count > 0) {
                snprintf(fn, 99, "%s/hypRaw%02d_%d_proof.csv", folder, per, job);
                nset_write_csv(rejSet, fn, 30, 0);
                
                nset_clear(rejSet);
                nset_init(rejSet, HYP_RAW_SET_EPS);
            }
            
            hyp_raw_job_file_name(fn, 100, per, job);
            nset_lock(actSet);
            ok = nset_write(actSet, fn);
            
            if(! ok) {
                printf("Could not write results to \"%s\", will stop here.\n\n", fn);
            } else {
                troots += actSet->count;
                
                jobDone(per, job, &ts);

                hyp_raw_map_file_name(fn, 100, per, job);
                writeETree(fn, actSet, 1);
            }
            
            nset_move(prevSet, actSet, 0);
        }
        
        if(! ok) {
            jobFailed(per, job);
        }
    }
    
    if(ok) {
        lapse(&ats, time);
        int j = end - start;
        
        if(per < 28) {
            long rt = prevSet->count;
            long rrt = prevSet->realCount;
            
            printf("\nMax dist from the level curve to a root: %.4Lf\n", maxDist);
            int all = mandel_hyp_count(per) == (2 * rt - rrt);
            if(all) {
                printf("\nFound all %ld roots (%ld real) of period %d in %s.\n", rt, rrt, per, time);
                fflush(stdout); // before nset_writeCsv
                
                char fn[100];
                snprintf(fn, 99, "%s/hyp%02d_%d.csv", folder, per, HYP_RAW_DIGITS);
                if(nset_write_csv(prevSet, fn, HYP_RAW_DIGITS, 0)) {
                    printf("Exported to %s.\n\n", fn);
                } else {
                    printf("Could not export the list to %s.\n\n", fn);
                }
            } else {
                printf("\nFound only %ld roots (%ld real) of period %d in %s.\n\n", rt, rrt, per, time);
            }
        } else {
            printf(j == 1 ? "\nCompleted %d job.\nFound and saved %ld centers in %s.\n\n" :
                   "\nCompleted %d jobs.\nFound and saved %ld centers in %s.\n\n", j, troots, time);
        }
    }
    
    fflush(stdout); // likely superfluous (at least flush error messages if not(ok) induced a quick jump to here)
    free(cc);
    
    mpc_clear(sp);
    mpc_clear(ep);
    
    clear();
    
    return ok;
}

// MARK: convenience function

nset hyp_new_set(bool locked) {
    return nset_new(HYP_RAW_SET_EPS, locked);
}

// MARK: the help system and the main function

static const char* before = "This task computes and saves to ./hyp/#PER/hyp#PER_job.nset the hyperbolic centers of period #per, as follows:\n\n";
static const char* after = "\nFor each increment in the period, the total computing time and the disk space required double.\nEach job is saved to a binary .nset file of size ~1.1 GB, up to 2.15 GB.\nThe RAM used should be less than 7 GB.\n\n";

static const char *parameters[] = {
    "period",
    "start",
    "end"
};

static const char *types[] = {
    "required",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "",
    ""
};

static const char *descriptions[] = {
    "the period, integer, at least 1, at most 41",
    "start job >= 0, inclusive; default value = 0",
    "end job <= 1 to pow(2,period - 27), exclusive; default value = jobs count"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 3;
static const int columnWidths[] = {18, 18, 18};

/// @brief Prints instructions for usage and some details about the command line arguments.
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

int hyp_raw_main(int argc, const char * argv[]) {
    int per, start, end;
    
    if(argc < 1 || sscanf(argv[0], "%d", &per) < 1 || per < 1 || per > HYP_RAW_MAX_PER) {
        help();
        
        return 1;
    }
    
    int jobs = per < 28 ? 1 : 1 << (per - 27);
    if(argc >= 3) {
        if(sscanf(argv[1], "%d", &start) < 1 || sscanf(argv[2], "%d", &end) < 1 ||
           start < 0 || end < start || end > jobs) {
               help();
               
               return 1;
        }
    } else {
        start = 0;
        end = jobs;
    }
    
    return ! findHypRaw(per, start, end);
}
