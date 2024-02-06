//
//  misRaw.c
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

#include "misRaw.h"
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
#include "misSets.h"

// MARK: constants and buffers definitions

#define LC_ITER      15
#define LC_REF        1
#define LCL_EPS       1E-18
#define RT_MIN_ITER  20
#define MP_MIN_ITER   5
#define MULT_ITER     2.6
#define LC_ERR        1E-13
#define NT_ERR        1E-13
#define MAX_DER2      1E30
#define MIN_DER2      1E-25
#define LC_ERR2      (LC_ERR * LC_ERR)
#define LC_CONV2      1E-36
#define LC_PREC      120
#define MIN_DER2      1E-25

// max distace of root from the starting point, by period, starting with period = 3
static const ldbl JUMP[] = {0.8, 0.6, 0.4, 0.3, 0.2, 0.15, 0.1, 0.09, 0.07, 0.06,
    0.05, 0.05, 0.04, 0.04, 0.03, 0.03, 0.025, 0.025, 0.02, 0.02,
    0.018, 0.017, 0.017, 0.017, 0.017, 0.017, 0.016, 0.016, 0.015, 0.015,
    0.015, 0.014, 0.014, 0.013, 0.013, 0.013, 0.013, 0.012, 0.011};

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
static mpc_struct sc[MIS_RAW_ANG];

// fp80 version of circle targets
static fp80_struct scl[MIS_RAW_ANG];

static mpc lcp, c, d, nt, vp, dp, c170, nt170, d170, c130, vp170, dp170;
static mpc_struct tc[MIS_RAW_MAX_TYPE];
static fp80 m1 = {-1, 0};
static ldbl maxDist = 0;

// folder name
static char folder[100];

// MARK: initialization functions

/// @brief Initializes buffers and creates the folders on the disk, if necessary.
///
/// @param pp the pre-period
/// @param per the period
///
/// @return @ref true if successfull, @ref false otherwise
static bool initSets(int pp, int per) {
    nset_init(prevSet, MIS_RAW_SET_EPS);
    nset_init(actSet, MIS_RAW_SET_EPS);
    nset_init(divSet, MIS_RAW_SET_EPS);
    nset_init(rejSet, MIS_RAW_SET_EPS);
    
    // prepare the output folder
    snprintf(folder, 99, "%s/%02d-%02d", MIS_FOLDER, pp, per);
    if(! dir(MIS_FOLDER) || ! dir(folder)) {
        printf("Could not create dir %s\n", folder);
        
        return false;
    }
    
    return true;
}

/// @brief Initializes the the target points used to construct the detailed level curve.
static void initCircles(void) {
    mpfr_t th;
    mpfr_init2(th, 128);
    for (int i = 0; i < MIS_RAW_ANG; i++) {
        mpc_init(sc + i, MIS_RAW_LC_PREC);
        mpfr_set_ui(th, i, MPFR_RNDN);
        mpfr_mul_2si(th, th, -MIS_RAW_ANG_2POW, MPFR_RNDN);
        
        mpc_exp_2Pi_i(sc + i, th);
        mpc_muld(sc + i, sc + i, MIS_RAW_R);
    }
    mpfr_clear(th);
    
    // fp80 version of circle targets
    for (int i = 0; i < MIS_RAW_ANG; i++) {
        scl[i].x = MIS_RAW_R * cosl(2 * PI * i / MIS_RAW_ANG);
        scl[i].y = MIS_RAW_R * sinl(2 * PI * i / MIS_RAW_ANG);
    }
    
    mpc_init(lcp, MIS_RAW_LC_PREC);
    mpc_init(c, MIS_RAW_RT_PREC);
    mpc_init(d, MIS_RAW_RT_PREC);
    mpc_init(nt, MIS_RAW_RT_PREC);
    mpc_init(vp, MIS_RAW_RT_PREC);
    mpc_init(dp, MIS_RAW_RT_PREC);
    mpc_init(c170, MIS_RAW_REF_PREC);
    mpc_init(d170, MIS_RAW_REF_PREC);
    mpc_init(nt170, MIS_RAW_REF_PREC);
    mpc_init(vp170, MIS_RAW_REF_PREC);
    mpc_init(dp170, MIS_RAW_REF_PREC);
    mpc_init(c130, 130);
    
    for (int i = 0; i < MIS_RAW_MAX_TYPE; i++) {
        mpc_init(tc + i, MIS_RAW_RT_PREC);
    }
}

// MARK: level sets and roots computing functions

/// @brief Computes the next point to @c sp on the level curve with given angle @c ang and high precision.
///
/// @param res the result
/// @param sp the starting point
/// @param pp the pre-period
/// @param per the period
/// @param ang the angle
/// @param iter the max number of iterates of the Newton method
///
/// @return @ref true if the point has been found and stored in @c res, @ref false otherwise
static inline bool mis_raw_find_misc(mpc res, mpc sp, int pp, int per, int ang, int iter) {
    ulong nt = mandel_tot_nt();
    
    bool ok = mandel_mis_root_ref(res, sp, &sc[(ang + MIS_RAW_ANG) % MIS_RAW_ANG], pp, per, iter, LC_REF);
    
    stats[7] += mandel_tot_nt() - nt;
    
    return ok;
}

/// @brief Quickly checks the unicity of the root @c c and that it is not a root of a period that divides @c per.
///
/// @param c the root to check
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

/// @brief Computes the most probable pre-periodic or periodic type of the parameter @c cr and then attempts a proof for this guess.
///
/// If proven to be o a lower type than @c pp+per, it is added to the list of divisors and return @ref true. Otherwise it returns @ref false.
///
/// @param cr the parameter
/// @param pp the pre-period
/// @param per the period
///
/// @return @ref true if proven to be of lower type than @c maxType, @ref false.
static bool checkLowerType(mpc cr, int pp, int per) {
    int maxType = pp + per;
    
    // cr is already checked to not be in {0, -1}
    mpc_set(tc, cr); // tc[i] = p_{i+1}(cr), here i == 0
    mpc_sqr(tc + 1, tc);
    mpc_add(tc + 1, tc + 1, tc); // p_2(cr)
    
    mpc_struct *buf = tc + (MIS_RAW_MAX_TYPE - 1);
    fp80 lb;
    
    int t, p;
    ldbl d;
    u128 u;
    for (int type = 3; type < maxType; type++) {
        t = type - 1;
        p = t - 1;
        
        mpc_sqr(tc + t, tc + p);
        mpc_add(tc + t, tc + t, tc);
        
        mpc_get80(lb, tc + t);
        d = fp80_mod(lb);
        if(d < HYP_RAW_PROOF) { // probably hyperbolic
            if(mandel_is_hyp(cr, type, HYP_RAW_PROOF)) {
                u128_set(u, cr);
                nset_add(divSet, u);
                
                stats[9] ++;
                
                return true;
            }
            
            return false;
        }
        
        for (int fp = 1; fp < t; fp ++) {
            int p = t - fp;
            if(per % p != 0) {
                continue;
            }
            
            mpc_sub(buf, tc + t, tc + fp);
            mpc_get80(lb, buf);
            
            d = fp80_mod(lb);
            // not enough precision to pass a proof with MIS_RAW_PROOF
            if(d < HYP_RAW_PROOF) { // probably pre-periodic of type (fp + 1, type - (fp + 1))
                if(mandel_is_mis(cr, fp + 1, p, HYP_RAW_PROOF)) {
                    u128_set(u, cr);
                    nset_add(divSet, u);
                    
                    stats[9] ++;
                    
                    return true;
                }
                
                return false;
            }
        }
    }
    
    return false;
}

/// @brief Verifies the unicity of a root and attempts to prove that it is a root with error at most @c MIS_RAW_PROOF
/// and that the Newton method converges in a disk centered at @c cr with radius @c MIS_RAW_CONV to @c cr.
///
/// @param pp the pre-period
/// @param per the period
/// @param cr the root
///
/// @return @ref true if @c c is a new root, @ref false otherwise
static bool checkAndAddRoot(int pp, int per, mpc cr) {
    u128 u;
    if(! u128_set(u, cr)) {
        return false;
    }
    
    mpc_set(c, cr);
    
    // for speed, check the probable type of c, lower than pp + per
    if(checkLowerType(c, pp, per)) {
        return false;
    }
    
    // if c is Mis(pp, per) and the Newton method converges in the disc of radius MIS_RAW_CONV
    // add it to the set of results, while checking for unicity
    if(mandel_is_mis(cr, pp, per, MIS_RAW_PROOF) &&
       mandel_conv_npp(cr, pp, per, MIS_RAW_PROOF, MIS_RAW_CONV, false)) {
        nset_add(actSet, u);
        
        return true;
    }
    
    // for speed, check the most probable cases first
    if(mandel_is_mis(c, pp - 1, per, MIS_RAW_PROOF)) {
        nset_add(divSet, u);
        
        stats[9] ++;
        
        return false;
    }
    
    if(per == 2 && mandel_is_mis(c, pp, 1, MIS_RAW_PROOF)) {
        nset_add(divSet, u);
        
        stats[9] ++;
        
        return false;
    }
    
    // check all the other cases
    for (int dper = per; dper >= 1; dper--) {
        if(per % dper != 0) {
            continue;
        }
        
        // special cases already checked
        for (int spp = dper == per ? pp - 2 : dper == per - 1 ? pp - 1 : pp; spp > 1; spp --) {
            if(mandel_is_mis(c, spp, dper, MIS_RAW_PROOF)) {
                nset_add(divSet, u);
                
                stats[9] ++;
                
                return false;
            }
        }

        if(dper > 2 && mandel_is_hyp(c, dper, MIS_RAW_PROOF)) {
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
/// @param pp the pre-period
/// @param per the period
/// @param iter the maximum number of iterates
/// @param err the longest valid step, or modulus of a Newton term
///
/// @return the number of Newton steps performed
static int accelRoot(fp80 v, fp80 sp, int pp, int per, int iter, double err) {
    fp80 nt, op = {sp->x, sp->y}, z = {0, 0};
    double e2 = err * err;
    ldbl den;
    
    for (int i = 0; i < iter; i++) {
        den = mandel_mis_ntl(nt, op, z, pp, per) >= MIN_DER2;
        
        stats[11] ++;
        
        if(nt->x * nt->x + nt->y * nt->y < e2 || den < 1e-25) {
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
/// @param d a buffer
/// @param vp a buffer
/// @param dp a buffer
/// @param pp the pre-period
/// @param per the period
static inline void newton(mpc nt, mpc c, mpc d, mpc vp, mpc dp, int pp, int per) {
    mpc_set(nt, c);
    mpc_seti(d, 1, 0);
    
    for (long k = 1; k < pp; k++) {
        mpc_mul(d, d, nt);
        mpc_scale(d, d, 1);
        mpc_addi(d, d, 1);
        
        mpc_sqr(nt, nt);
        mpc_add(nt, nt, c);
    }
    
    mpc_set(vp, nt);
    mpc_set(dp, d);
    
    for (int i = 0; i < per; i++) {
        mpc_mul(d, d, nt);
        mpc_scale(d, d, 1);
        mpc_addi(d, d, 1);
        
        mpc_sqr(nt, nt);
        mpc_add(nt, nt, c);
    }
    
    mpc_sub(nt, nt, vp);
    mpc_sub(d, d, dp);
    
    mpc_div(nt, nt, d);
}

/// @brief Refines a root @c sp of @c p_{per} and checks as soon as possible for unicity.
///
/// @param rc the result
/// @param sp the starting point
/// @param pp the pre-period
/// @param per the period
/// @param qIter the max number of iterates with normal precision
/// @param qPrec the minimal precision of the root to validate the convergence is @c 2^{-qPrec}
/// @param rIter the maximum number of iterates with extended precision
/// @param rPrec the minimal extended precision of the root to validate the convergence is @c 2^{-rPrec}
///
/// @return @ref true if the Newton method converged to a new root, @ref false otherwise
static inline bool refine(mpc rc, mpc sp, int pp, int per, int qIter, int qPrec, int rIter, int rPrec) {
    mpc_set(c, sp);
    
    bool conv = false;
    long exp = 0, prexp = 0, checked = 0;
    for (int i = 0; ! conv && i < qIter; i++) {
        prexp = exp;
        
        newton(nt, c, d, vp, dp, pp, per);
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
    
    conv = false;
    long rexp;
    for (int i = 0; ! conv && i < rIter; i++) {
        newton(nt170, c170, d170, vp170, dp170, pp, per);
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

/// @brief Searches for a root of @c p_{per} which, if found, is checked for unicity, refined, proven to be correct
/// and eventually added to the set of points.
///
/// @param pp the pre-period
/// @param per the period
/// @param sp the starting point
/// @param sp8 the starting point with low precision
/// @param prevSp8 the previous starting point with low precision
/// @param iter the ma number of iterates
/// @param jump2 the square of the max allowed distance between the starting point and the root
///
/// @return @ref true ifsuccessfull, @ref false otherwise
static bool searchRoot(int pp, int per, mpc sp, fp80 sp8, fp80 prevSp8, int iter, ldbl jump2) {
    int it = iter + MP_MIN_ITER;
    bool rf = true; // root found
    
    fp80 c8;
    long nt80 = stats[11];
    if(fp80_dist(sp8, prevSp8) >= LC_ERR) { // try fp80 acceleration
        int acIt = accelRoot(c8, sp8, pp, per, iter, NT_ERR);
        it -= acIt;

        if(acIt == 0) { // no useful Newton term
            mpc_set(c130, sp);
        } else if(it <= MP_MIN_ITER) { // did not converge
            rf = false;
            
            stats[10] ++;
        } else if(fp80_dist2(c8, sp8) >= jump2) { // jumped too far
            rf = false;
            
            stats[18] ++;
        } else if(fp80_mod(c8) <= 0.2 && fp80_dist(c8, m1) <= 0.15) { // close to 0 or -1
            rf = false;
            
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
        
        if(refine(c170, c130, pp, per, it, MIS_RAW_QCONV_2POW, 3, MIS_RAW_RCONV_2POW)) {
            if(checkAndAddRoot(pp, per, c170)) {
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

/// @brief Computes the next point on the level curve using low precision for as long as possible, then switching
/// to high precision.
///
/// @param pp the pre-period
/// @param per the period
/// @param ep the end point
/// @param ep8 the end point in low precision
/// @param sp the starting point
/// @param sp8 the starting point with low precision
/// @param prevSp8 the previous starting point with low precision
/// @param ang the angle
///
/// @return @ref true if sucessfull, @ref false otherwise
static bool nextLc(int pp, int period, mpc ep, fp80 ep8, mpc sp, fp80 sp8, fp80 prevSp8, int ang) {
    int n = pp + period;
    int iter = 3 * n / 4;
    iter = iter < LC_ITER ? LC_ITER : iter;
    
    if(fp80_dist(sp8, prevSp8) < LC_ERR) { // consecutive points on the LC too close, use mp
        if(! mis_raw_find_misc(ep, sp, pp, period, ang, iter)) {
            return false;
        }
        
        mpc_get80(ep8, ep);
        
        return true;
    }
    
    // perform at least a few steps with fp80, then hadover to mp
    fp80 nt, p8 = {sp8->x, sp8->y}, op8;
    fp80_ptr t = &scl[(ang + MIS_RAW_ANG) % MIS_RAW_ANG];
    
    ldbl d2, m2;
    int it = 0;
    do {
        *op8 = *p8;
        d2 = mandel_mis_ntl(nt, p8, t, pp, period);
        fp80_sub(p8, p8, nt);
        m2 = fp80_mod2(nt);
        
        it ++;
        
        stats[6] ++;
    } while(it < iter && d2 <= MAX_DER2 && d2 >= MIN_DER2 && m2 >= LC_ERR2);
    
    if(it > 1 && d2 <= MAX_DER2 && d2 >= MIN_DER2 && m2 < LC_ERR2) { // it may converge in fp80, let's see
        d2 = mandel_mis_ntl(nt, p8, t, pp, period);
        m2 = fp80_mod2(nt);
        
        stats[6] ++;
        
        if(d2 <= MAX_DER2 && d2 >= MIN_DER2 && m2 < LC_CONV2) { // yes, keep it
            fp80_sub(ep8, p8, nt);
            mpc_set80(ep, ep8);
            
            return true;
        }
    }
    
    if(it == 1) { // if cannot start in fp80, do everything in mp
        if(! mis_raw_find_misc(ep, sp, pp, period, ang, iter)) {
            return false;
        }
        
        mpc_get80(ep8, ep);
        
        return true;
    }
    
    // otherwise, continue from op8
    mpc_set80(ep, op8);
    if(! mis_raw_find_misc(ep, ep, pp, period, ang, iter - it + 1)) {
        return false;
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
static bool roots(int pp, int per, mpc sp, mpc ep, int spCount, int step) {
    int n = pp + per;
    
    mpc_set(lcp, sp);
    fp80 p8, op8, tp8 = {0, 0};
    mpc_get80(p8, lcp);
    *op8 = *p8;
    
//    ldbl jump2 = JUMP[n - 3] * JUMP[n - 3];
    
    // test if there are large jumps
    ldbl jump2 = JUMP[n - 3] * JUMP[n - 3] * 4;
    int iter = n * MULT_ITER;
    iter = iter < RT_MIN_ITER ? RT_MIN_ITER : iter;
    bool ok = true;
    
    for (int i = spCount - 1; i >= 0 && ok; i--) {
        if(ok) {
            searchRoot(pp, per, lcp, p8, op8, iter, jump2);
        }
        
        for (int j = step - 1; j >= 0 && ok; j--) {
            *tp8 = *p8;
            ok = ok && nextLc(pp, per, lcp, p8, lcp, p8, op8, step * i + j);
            *op8 = *tp8;
        }
    }
    
    // refine lcp, it may come from fp80 computation
    ok = ok && mandel_mis_root_ref(lcp, lcp, sc, pp, per, LC_ITER, LC_REF);
    ok = ok && mpc_distl(ep, lcp) < 1E-27;
    
    return ok;
}

/// @brief Frees the buffers.
static void clear(void) {
    nset_clear(prevSet);
    nset_clear(actSet);
    nset_clear(divSet);
    nset_clear(rejSet);
    
    for (int i = 0; i < MIS_RAW_ANG; i++) {
        mpc_clear(sc + i);
    }
    
    mpc_clear(lcp);
    mpc_clear(c);
    mpc_clear(d);
    mpc_clear(nt);
    mpc_clear(vp);
    mpc_clear(dp);
    mpc_clear(c170);
    mpc_clear(d170);
    mpc_clear(nt170);
    mpc_clear(c130);
    mpc_clear(vp170);
    mpc_clear(dp170);
    
    for (int i = 0; i < MIS_RAW_MAX_TYPE; i++) {
        mpc_clear(tc + i);
    }
}

/// @brief Computes and writes the tree map of the set of points @c ps.
///
/// @param fileName the file name
/// @param ps the set of point
/// @param check @ref true to check that the file reads correctly, @ref false to ignore this step
static int writeETree(char *fileName, nset_t ps, int check) {
    tmap tr = tmap_map(ps, MIS_RAW_TMAP_MAX_LEVEL, MIS_RAW_TMAP_LEVEL_STEP);
    
    if(tr == NULL) {
        printf("Could not construct the eTree.\n");
    } else {
        if(tmap_save(tr, fileName)) {
            printf("Points map written to %s.\n", fileName);
        } else {
            printf("Could not write the eTree to %s.\n", fileName);
            
            return 0;
        }
    }
    
    int ok = 1;
    if(check) {
        tmap tc = tmap_load(fileName, 0, true, NULL);
        
        if(tc == NULL) {
            printf("Could not read %s.\n", fileName);
            
            ok = 0;
        } else {
            if(tmap_eq(tr, tc)) {
//                printf("The map file is correct.\n");
            } else {
                printf("The map file is DAMAGED !!!\n");
                
                ok = 0;
            }
        }
        
        tmap_free(tc);
    }
    
    tmap_free(tr);
    
    return 1;
}

int mis_raw_job_file_name(char fn[], int max, int pp, int per, int job) {
    if(pp + per < 28) {
        return snprintf(fn, max - 1, "%s/%02d-%02d/mis%02d-%02d.nset",
                        MIS_FOLDER, pp, per, pp, per);
    } else {
        return snprintf(fn, max - 1, "%s/%02d-%02d/misRaw%02d-%02d_%d.nset",
                        MIS_FOLDER, pp, per, pp, per, job);
    }
}

int mis_raw_map_file_name(char fn[], int max, int pp, int per, int job) {
    if(pp + per < 28) {
        return snprintf(fn, max - 1, "%s/%02d-%02d/mis%02d-%02d.tmap",
                        MIS_FOLDER, pp, per, pp, per);
    } else {
        if(job >= (1 << (pp + per - 27))) {
            return snprintf(fn, max - 1, "%s/%02d-%02d/misRaw%02d-%02d.tmap",
                            MIS_FOLDER, pp, per, pp, per);
        }
        
        return snprintf(fn, max - 1, "%s/%02d-%02d/misRaw%02d-%02d_%d.tmap",
                        MIS_FOLDER, pp, per, pp, per, job);
    }
}

/// @brief Provides the name of the stats file for the job with index  @c job and of period @c per.
///
/// @param fn the file name
/// @param max the max length of the file name
/// @param per the period
/// @param job the index of the file (@c 0 for periods samller then @c 28)
static int mis_statsFileName(char fn[], int max, int pp, int per, int job) {
    if(pp + per < 28) {
        return snprintf(fn, max - 1, "%s/%02d-%02d/mis%02d-%02d_stats.csv",
                        MIS_FOLDER, pp, per, pp, per);
    } else {
        return snprintf(fn, max - 1, "%s/%02d-%02d/misRaw%02d-%02d_%d_stats.csv",
                        MIS_FOLDER, pp, per, pp, per, job);
    }
}

/// @brief Called to mark this job with state @c state in the stats file.
///
/// @param per the period
/// @param job the job
static void writeStats(int pp, int per, int job, int state) {
    char fn[100];
    mis_statsFileName(fn, 99, pp, per, job);
    
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
/// @param pp the pre-period
/// @param per the period
/// @param job the job
static void jobStarted(int pp, int per, int job, struct timeb *ts) {
    ftime(ts);
    char date[80];
    time_stamp(date, 80, true, true);
    
    printf("\nJob %d started at %s\n", job, date);
    fflush(stdout);
    
    stats[0] = ts->time;
    stats[0] = stats[0] * 1000 + ts->millitm;      // start time, see COLUMS[2 ... 17]
    for (int i = 1; i < columnsLen - 2; i++) {
        stats[i] = 0;
    }
    
    writeStats(pp, per, job, 0);
}

/// @brief Called to mark this job as @c Completed in the stats file.
///
/// @param pp the pre-period
/// @param per the period
/// @param job the job
static void jobDone(int pp, int per, int job, struct timeb *ts) {
    char time[80];
    
    lapse(ts, time);
    printf("Job %d completed, %ld points found in %s.\n", job, actSet->count, time);
    fflush(stdout);
    
    struct timeb tts;
    ftime(&tts);
    
    stats[1] = tts.time;
    stats[1] = stats[1] * 1000 + tts.millitm;      // end time, see COLUMS[2 ... 17]
    stats[2] = stats[1] - stats[0];  // duration
    stats[3] = actSet->count;
    stats[4] = actSet->realCount;
    stats[5] = per < 27 ? 1L << per : 1L << 27;
    
    writeStats(pp, per, job, 1);
}

/// @brief Called to mark this job as @c Failed in the stats file.
///
/// @param pp the pre-period
/// @param per the period
/// @param job the job
static void jobFailed(int pp, int per, int job) {
    writeStats(pp, per, job, 2);
}

int mis_raw_jobs_count(int pp, int per) {
    int type = pp + per;
    
    return type < 28 || type > 50 ? 1 : 1 << (type - 27);
}

nSet_struct *mis_raw_load_job(int pp, int per, int job) {
    nSet_struct *ps = malloc(sizeof(nSet_struct));
    nset_init(ps, MIS_RAW_SET_EPS);
    
    char fn[100];
    mis_raw_job_file_name(fn, 99, pp, per, job);
    
    if(nset_read(ps, fn, false)) {
        return ps;
    } else {
        nset_clear(ps);
        free(ps);
        
        return NULL;
    }
}

/// @brief Searches for the roots of the polynomial @c p_{pp+per}-p_{pp} and saves the results.
///
/// Constructs a tree
/// map of the points found. Starting with period @c 28, the serach is divided in @c 2^{per-27} jobs,
/// for parallelization.
///
/// The names of the corresponding files are provided by @c mis_raw_jobFileName() and the results
/// can be loaded with @c mis_raw_loadJob(). The map file name is given by @c mis_raw_mapFileName().
///
/// @param pp the pre-period
/// @param per the period
/// @param start the first job
/// @param end last job @c +1
///
/// @return @ref true if the operation completed successfully, @ref false otherwise
static bool find_mis_raw(int pp, int per, int start, int end) {
    char hostname[_SC_HOST_NAME_MAX + 1];
    gethostname(hostname, _SC_HOST_NAME_MAX + 1);
    
    int cpu = get_cpu_id();
    
    printf("\nComputing pre-periodic points of type (%d, %d), jobs %d to %d (max jump %.4Lf) on %s with cpu %d:\n\n",
           pp, per, start, end, JUMP[pp + per - 3], hostname, cpu);
    
    fflush(stdout); // Flushing the buffer from time to time ensures that the logs survive an app kill
    
    struct timeb ats, ts;
    ftime(&ats);
    
    levm cc = miss_load(pp, per);
    
    // TODO: remove
    mpc a, b;
    mpc_init(a, cc->prec);
    mpc_init(b, cc->prec);
    levm_point(a, cc, 0);
    levm_point(b, cc, 1);
    
    int jobs = mis_raw_jobs_count(pp, per);
    
    if(cc == NULL || levm_segs(cc) % jobs != 0) {
        printf("Could not prepare the level curve, will stop here.\n\n");
        fflush(stdout);
        
        return false;
    }
    
    int jobSize = (int) (levm_segs(cc) / jobs);
    cc = levm_sub_curve(cc, jobSize * start, 1, jobSize * (end - start) + 1);
    if(mpfr_get_ld(cc->radius, MPFR_RNDN) != MIS_RAW_R) {
        mpfr_t r;
        mpfr_init2(r, 120);
        mpfr_set_ld(r, MIS_RAW_R, MPFR_RNDN);
        
        cc = levm_level(cc, r);
        
        mpfr_clear(r);
        
        if(cc == NULL) {
            printf("Could not prepare the level curve, will stop here.\n\n");
            fflush(stdout);
            
            return false;
        }
    }
    
    char time[80], fn[100];
    lapse(&ats, time);
    
    printf("Level curve prepared in %s.\n", time);
    fflush(stdout);
    
    // prepare sets
    initSets(pp, per);
    
    // init level curve targets
    initCircles();
    
    // try to load previous job
    if(start > 0) {
        ftime(&ts);
        
        mis_raw_job_file_name(fn, 100, pp, per, start - 1);
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
        jobStarted(pp, per, job, &ts);
        
        int n = pp + per;
        // ~ MIS_RAW_ANG / MIS_RAW_STEP starting points per root
        int spCount = n < MISS_MIN_PER ? 1 << n :
                      1 << (n - 1 - MISS_2POW + MIS_RAW_ANG_2POW - MIS_RAW_STEP_2POW);
        for (int seg = 0; ok && seg < jobSize; seg++) {
            int pos = jobSize * (job - start) + seg;
            levm_point(ep, cc, pos);
            levm_point(sp, cc, pos + 1);
            
            ok = roots(pp, per, sp, ep, spCount, MIS_RAW_STEP);
            
            if(! ok) {
                printf("Could not compute the pre-periodic points from the segment %d of the level curve, will stop here.\n\n", pos);
            }
        }
        
        if(ok) {
            if(rejSet->count > 0) {
                snprintf(fn, 99, "%s/misRaw%02d-%02d_%d_proof.csv", folder, pp, per, job);
                nset_write_csv(rejSet, fn, 38, 0);
                
                nset_clear(rejSet);
                nset_init(rejSet, MIS_RAW_SET_EPS);
            }
            
            mis_raw_job_file_name(fn, 100, pp, per, job);
            nset_lock(actSet);
            ok = nset_write(actSet, fn);
            
            if(! ok) {
                printf("Could not write results to \"%s\", will stop here.\n\n", fn);
            } else {
                troots += actSet->count;
                
                jobDone(pp, per, job, &ts);

                mis_raw_map_file_name(fn, 100, pp, per, job);
                writeETree(fn, actSet, 1);
            }
            
            nset_move(prevSet, actSet, 0);
        }
        
        if(! ok) {
            jobFailed(pp, per, job);
        }
    }
    
    if(ok) {
        lapse(&ats, time);
        int j = end - start;
        
        if(per < 28) {
            long rt = prevSet->count;
            long rrt = prevSet->realCount;
            
            printf("\nMax dist from the level curve to a pre-periodic point: %.4Lf\n", maxDist);
            int all = mandel_mis_count(pp, per) == (2 * rt - rrt);
            if(all) {
                printf("\nFound all %ld pre-periodic points (%ld real) of type (%d, %d) in %s.\n",
                       rt, rrt, pp, per, time);
                fflush(stdout); // before nset_writeCsv
                
                char fn[100];
                snprintf(fn, 99, "%s/mis%02d-%02d_%d.csv", folder, pp, per, MIS_RAW_DIGITS);
                if(nset_write_csv(prevSet, fn, MIS_RAW_DIGITS, 0)) {
                    printf("Exported to %s.\n\n", fn);
                } else {
                    printf("Could not export the list to %s.\n\n", fn);
                }
            } else {
                printf("\nFound only %ld pre-periodic points (%ld real) of type (%d, %d) in %s.\n\n",
                       rt, rrt, pp, per, time);
            }
        } else {
            printf(j == 1 ? "\nCompleted %d job.\nFound and saved %ld pre-periodic points in %s.\n\n" :
                   "\nCompleted %d jobs.\nFound and saved %ld pre-periodic points in %s.\n\n", j, troots, time);
        }
    }
    
    fflush(stdout); // likely superfluous (at least flush error messages if not(ok) induced a quick jump to here)
    free(cc);
    
    mpc_clear(sp);
    mpc_clear(ep);
    
    clear();
    
    return ok;
}

nset mis_raw_mini(int pp, int per, int seg) {
    printf("\nComputing pre-periodic points of type (%d, %d), from the segment %d of the levelcurve:\n",
           pp, per, seg);
    printf("Using the folowing parameters:\n\tLC_ERR %le\n\tNT_ERR %le\n\tMIS_RAW_CONV %le\n",
           LC_ERR, NT_ERR, MIS_RAW_CONV);
    printf("\tMIS_RAW_PROOF %le\n\tMIS_RAW_ANG_2POW %d\n\tMIS_RAW_STEP_2POW %d\n\tMULT_ITER %lg\n\n",
           MIS_RAW_PROOF, MIS_RAW_ANG_2POW, MIS_RAW_STEP_2POW, MULT_ITER);
    fflush(stdout); // Flushing the buffer from time to time ensures that the logs survive an app kill
    
    struct timeb ats;
    ftime(&ats);
    
    levm cc = miss_load(pp, per);
    
    if(cc == NULL || levm_segs(cc) < seg) {
        printf("Could not prepare the level curve, will stop here.\n\n");
        fflush(stdout);
        
        return false;
    }
    
    cc = levm_sub_curve(cc, seg, 1, 2);
    if(mpfr_get_ld(cc->radius, MPFR_RNDN) != MIS_RAW_R) {
        mpfr_t r;
        mpfr_init2(r, 120);
        mpfr_set_ld(r, MIS_RAW_R, MPFR_RNDN);
        
        cc = levm_level(cc, r);
        
        mpfr_clear(r);
        
        if(cc == NULL) {
            printf("Could not prepare the level curve, will stop here.\n\n");
            fflush(stdout);
            
            return false;
        }
    }
    
    char time[80];
    lapse(&ats, time);
    
    printf("Level curve prepared in %s.\n", time);
    fflush(stdout);
    
    // prepare sets
    initSets(pp, per);
    
    // init level curve targets
    initCircles();
    
    bool ok = true;
    mpc sp, ep;
    mpc_init(sp, cc->prec);
    mpc_init(ep, cc->prec);
    
    int n = pp + per;
    // ~ MIS_RAW_ANG / MIS_RAW_STEP starting points per root
    int spCount = n < MISS_MIN_PER ? 1 << n :
    1 << (n - 1 - MISS_2POW + MIS_RAW_ANG_2POW - MIS_RAW_STEP_2POW);
    
    levm_point(ep, cc, 0);
    levm_point(sp, cc, 1);
    
    mpfr_printf("Start (%.22Rf, %.22Rf), end (%.22Rf, %.22Rf).\n", sp->x, sp->y, ep->x, ep->y);
    
    ok = roots(pp, per, sp, ep, spCount, MIS_RAW_STEP);
    
    if(! ok) {
        printf("Could not compute the pre-periodic points from the segment %d of the level curve, will stop here.\n\n", seg);
        
        return NULL;
    }
    
    nset_lock(actSet);
    
    lapse(&ats, time);
    printf("Found %ld pre-periodic points in %s\n", actSet->count, time);
    fflush(stdout); // likely superfluous (at least flush error messages if not(ok) induced a quick jump to here)
    free(cc);
    
    mpc_clear(sp);
    mpc_clear(ep);
    
    nset res = nset_new(1, false);
    nset_copy(res, actSet);
    
    clear();
    
    return res;
}

// MARK: convenience function

nset mis_new_set(bool locked) {
    return nset_new(MIS_RAW_SET_EPS, locked);
}

// MARK: the help system and the main function

static const char* before = "This task computes and saves to ./mis/#PP-#PER/mis#PP-#PER_job.nset the pre-periodic points of given type, as follows:\n\n";
static const char* after = "\nFor each increment in the period, the total computing time and the disk space required double.\nEach job is saved to a binary .nset file of size ~1.1 GB, up to 2.15 GB.\nThe RAM used should be less than 7 GB.\n\n";

static const char *parameters[] = {
    "pre-period",
    "period",
    "start",
    "end"
};

static const char *types[] = {
    "required",
    "required",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "",
    "",
    ""
};

static const char *descriptions[] = {
    "the pre-period, integer, at least 2, at most 37",
    "the period, integer, at least 1, at most 38 - pre-period",
    "start job >= 0, inclusive; default value = 0",
    "end job <= 1 to pow(2, period - 27), exclusive; default value = jobs count"
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
    snprintf(format, 45, "    %%-%ds %%-%ds %%-%ds %%s\n", columnWidths[0], columnWidths[1], columnWidths[2]);
    
    printf(format, headers[0], headers[1], headers[2], headers[3]);
    printf("\n");
    for(int i = 0; i < paramCount; i++) {
        printf(format, parameters[i], types[i], defaults[i], descriptions[i]);
    }
    
    printf("%s", after);
}

int mis_raw_main(int argc, const char * argv[]) {
    int pp, per, start, end;
    
    if(argc < 2 || sscanf(argv[0], "%d", &pp) < 1 || sscanf(argv[1], "%d", &per) < 1 ||
       pp + per < 3 || pp + per > MIS_RAW_MAX_TYPE || pp < 2 || per < 1) {
        help();
        
        return 1;
    }
    
    int jobs = mis_raw_jobs_count(pp, per);
    if(argc >= 4) {
        if(sscanf(argv[2], "%d", &start) < 1 || sscanf(argv[3], "%d", &end) < 1 ||
           start < 0 || end < start || end > jobs) {
               help();
               
               return 1;
        }
    } else {
        start = 0;
        end = jobs;
    }
    
    return ! find_mis_raw(pp, per, start, end);
}
