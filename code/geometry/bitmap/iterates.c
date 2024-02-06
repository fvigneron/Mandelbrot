//
//  iterates.c
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2022.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2022.
//
//  Copyright 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the Lesser GNU Public Licence v3.0
//

#include <stdlib.h>
#include <limits.h>

#include "planarSet.h"
#include "mandel.h"
#include "iterates.h"
#include "mthread.h"

// MARK: bitmaps of the Mandelbrot set

bmap iterM_bitmap80(drect r, int maxIter, int threads) {
    if(r == NULL || r->w == 0 || r->h == 0 || maxIter <= 1 ||
       threads < 1 || threads > ITER_MAX_THREADS) {
        return NULL;
    }
    
    if(threads > 1) {
        iterBmap it = {.prec = 64, .type = ITER_TYPE_M, .subType = ITERM_ESCAPE
            ,.maxIter = maxIter, .r = *r};
        
        bmap *bms = iter_bitmaps(&it, threads);
        if(bms == NULL) {
            return NULL;
        }
        
        bmap bm = bms[0];
        free(bms);
        
        return bm;
    }
    
    bmap bm = bmap_new(r);
    bm->type = BMAP_TYPE_MANDEL;
    bm->subType = BMAP_SUB_TYPE_ESCAPE;
    
    bm->colorMap = BMAP_COLOR_MAP_LOG;
    bm->power = 0.5;
    bm->mapA /= log(maxIter + 1);
    bm->colorLow = 0xFFC0C0C0;
    
    fp80 c, v, dz, dc, b;
    
    int w = r->w, h = r->h;
    ldbl ip2 = ldexpl(100, 2 * r->tpow);
    ldbl p2 = 1 / ip2;
    ip2 *= ip2;
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            drect_rel_to_abs_coords80(c, r, x, y);
            *v = *c;
            dz->x = 1;
            dz->y = 0;
            *dc = *dz;
            
            ulong iter = 1;
            ldbl m2 = fp80_mod2(v);
            ldbl dz2 = 1;
            ldbl dc2 = 1;
            for (; iter < maxIter && m2 < 25 && dc2 < ip2 && dz2 > p2; iter++) {
                b->x = 2 * v->x;
                b->y = 2 * v->y;
                
                fp80_mul(dc, dc, b);
                dc->x += 1;
                
                fp80_mul(dz, dz, b);
                
                fp80_sqr(v, v);
                fp80_add(v, v, c);
                
                m2 = fp80_mod2(v);
                dz2 = fp80_mod2(dz);
                dc2 = fp80_mod2(dc);
            }
            
            ulong pix = m2 >= 25 ? iter : dz2 < 1 ? 0 : maxIter + 1;
            bmap_set_pixel(bm, pix, x, y);
        }
    }
        
    return bm;
}

bmap iterM_bitmap(drect r, mpc delta, ulong maxIter, int prec, int threads) {
    if(r == NULL || r->w == 0 || r->h == 0 || prec < 32 || maxIter <= 1 ||
       threads < 1 || threads > ITER_MAX_THREADS) {
        return NULL;
    }
    
    if(threads > 1) {
        iterBmap it = {.prec = prec, .type = ITER_TYPE_M, .subType = ITERM_ESCAPE
            ,.maxIter = maxIter, .r = *r, .usesDelta = delta != NULL};
        
        if(delta != NULL) {
            mpc_init((mpc_ptr) it.delta, mpc_prec(delta));
            mpc_set((mpc_ptr) it.delta, delta);
        }
        
        bmap *bms = iter_bitmaps(&it, threads);
        if(delta != NULL) {
            mpc_clear((mpc_ptr) it.delta);
        }
        
        if(bms == NULL) {
            return NULL;
        }
        
        bmap bm = bms[0];
        free(bms);
        
        return bm;
    }
    
    bmap bm = bmap_new(r);
    bm->type = BMAP_TYPE_MANDEL;
    bm->subType = BMAP_SUB_TYPE_ESCAPE;
    bm->colorLow = 0xFFC0C0C0;
    
    if(delta != NULL) {
        uint prec = (uint) mpc_prec(delta);
        bm->useHD = true;
        mpfr_inits2(prec, bm->hdx, bm->hdy, NULL);
        mpfr_set(bm->hdx, delta->x, MPFR_RNDN);
        mpfr_set(bm->hdy, delta->y, MPFR_RNDN);
    }
    
    mpc c, v;
    mpc_inits(prec, c, v, NULL);
    
    int w = r->w, h = r->h;
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            drect_rel_to_abs_coords(c, r, x, y);
            if(delta != NULL) {
                mpc_add(c, c, delta);
            }
            
            mpc_set(v, c);
            
            ulong iter = 1;
            for (; iter < maxIter && mpc_modl(v) < 5; iter++) {
                mpc_sqr(v, v);
                mpc_add(v, v, c);
            }
            
            ulong pix = mpc_modl(v) >= 5 ? iter : maxIter + 1;
            bmap_set_pixel(bm, pix, x, y);
        }
    }
    
    mpc_clears(c, v, NULL);
    
    bm->a /= maxIter + 1;
    
    return bm;
}

bmap *iterM_bitmaps80(drect r, int maxIter, int types, int threads) {
    if(r == NULL || r->w == 0 || r->h == 0 || types == 0 || maxIter <= 1 ||
       threads < 1 || threads > ITER_MAX_THREADS) {
        return NULL;
    }
    
    if (types == ITERM_ESCAPE) {
        bmap *bm = malloc(sizeof(bmap));
        
        bm[0] = iterM_bitmap80(r, maxIter, threads);
        
        return bm;
    }
    
    if(threads > 1) {
        iterBmap it = {.prec = 64, .type = ITER_TYPE_M, .subType = types
            ,.maxIter = maxIter, .r = *r};
        
        return iter_bitmaps(&it, threads);
    }
    
    int mps[] = {(types & ITERM_ESCAPE) != 0, (types & ITERM_DIST) != 0,
        (types & ITERM_DERZ) != 0, (types & ITERM_DERC) != 0, (types & ITERM_ESCPOW) != 0,
        (types & ITERM_SMOD) != 0, (types & ITERM_SARG) != 0
    };
    
    bmap bms[ITERM_COUNT];
    int bits[] = {ITER_DIST_BITS, ITER_DERZ_BITS, ITER_DERC_BITS,
        ITER_POW_BITS, ITER_MOD_BITS, ITER_ARG_BITS};
    ldbl ia[ITERM_COUNT - 1];
    int imgCount = 0;
    for (int i = 0; i < ITERM_COUNT; i++) {
        if(mps[i]) {
            bms[i] = bmap_new(r);
            bms[i]->type = i >= 5 ? BMAP_TYPE_S_FUNCTION : BMAP_TYPE_MANDEL;
            bms[i]->subType = i + 1;
            
            imgCount ++;
            
            if(i > 0) {
                ia[i - 1] = ldexp(1, bits[i - 1]);
                bms[i]->a = 1 / ia[i - 1];
                
                switch(i) {
                    case 5:
                        bms[i]->sgn = true;
                        bms[i]->mapA = 0.1;
                        bms[i]->mapB = 0.5;
                        break;
                    case 6:
                        bms[i]->sgn = true;
                        bms[i]->mapA = 1 / (2 * PI);
                        bms[i]->mapB = 0.5;
                        break;
                    default:
                        break;
                }
            }
        } else {
            bms[i] = NULL;
        }
    }
    
    if(mps[0]) {
        bms[0]->colorMap = BMAP_COLOR_MAP_LOG;
        bms[0]->power = 0.5;
        bms[0]->mapA /= log(maxIter + 1);
    }
    
    fp80 c, v, dz, dc, b;
    int w = r->w, h = r->h;
    ldbl ip2 = ldexpl(100, 2 * r->tpow);
    ldbl p2 = 1 / ip2;
    ip2 *= ip2;
    
    ldbl minArg = 4;
    ldbl maxArg = -4;
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            drect_rel_to_abs_coords80(c, r, x, y);
            *v = *c;
            
            dc->x = 1;
            dc->y = 0;
            *dz = *dc;
            
            ulong iter = 1;
            ldbl m2 = fp80_mod2(v);
            ldbl dz2 = 1;
            ldbl dc2 = 1;
            for (; iter < maxIter && m2 < 25 && dc2 < ip2 && dz2 > p2; iter++) {
                b->x = 2 * v->x;
                b->y = 2 * v->y;
                
                fp80_mul(dc, dc, b);
                dc->x += 1;
                
                fp80_mul(dz, dz, b);
                
                fp80_sqr(v, v);
                fp80_add(v, v, c);
                
                m2 = fp80_mod2(v);
                dz2 = fp80_mod2(dz);
                dc2 = fp80_mod2(dc);
            }
            
            if(mps[0]) { // escape time
                ulong pix = m2 >= 25 ? iter : dz2 < 1 ? 0 : maxIter + 1;
                
                bmap_set_pixel(bms[0], pix, x, y);
            }
            
            if(mps[2]) { // z derivative
                if(iter <= 1) {
                    bmap_set_pixel(bms[2], 0, x, y);
                } else {
                    ldbl mdz = fp80_mod(dz);
                    ldbl ad = ia[1] * expl(logl(mdz) / (iter - 1));
                    ulong pix = ad > ULONG_MAX ? ULONG_MAX : (ulong) ad;
                    
                    bmap_set_pixel(bms[2], pix, x, y);
                }
            }
            
            if(mps[3]) { // c derivative
                if(iter <= 1) {
                    bmap_set_pixel(bms[3], 0, x, y);
                } else {
                    ldbl mdc = fp80_mod(dc);
                    ldbl ad = ia[2] * expl(logl(mdc) / (iter - 1));
                    ulong pix = ad > ULONG_MAX ? ULONG_MAX : (ulong) ad;
                    
                    bmap_set_pixel(bms[3], pix, x, y);
                }
            }
            
            if(mps[1] || mps[4]) { // distance to M, need to iterate some more
                if(m2 < 25) { // did not escape, distance is 0, mark as max if in the interior
                    if(mps[1]) {
                        bmap_set_pixel(bms[1], dz2 < 1 ? ULONG_MAX : 0, x, y);
                    }
                    
                    if(mps[4]) {
                        bmap_set_pixel(bms[4], 0, x, y);
                    }
                } else {
                    ulong esc = iter;
                    
                    for (; m2 < 1E43; iter++) {
                        fp80_mul(dc, dc, v);
                        fp80_muli(dc, dc, 2);
                        dc->x += 1;
                        
                        fp80_sqr(v, v);
                        fp80_add(v, v, c);
                        
                        m2 = fp80_mod2(v);
                    }
                    
                    int k = (int) (iter - 1);
                    
                    ldbl dm2 = dc->x * dc->x + dc->y * dc->y;
                    ldbl d = sqrtl(m2 / dm2) * 0.25;
                    
                    ldbl lm2p = logl(m2 * 0.25);
                    
                    if(iter < 10000) { // only so that the exponent of ldbl does not overflow
                        d *= ldexpl(expm1l(ldexpl(lm2p, -k - 1)), k);
                    } else {
                        d *= 0.5 * lm2p;
                    }
                    
                    if(mps[1]) {
                        ldbl m = d * ia[0];
                        ulong pix = m > ULONG_MAX ? ULONG_MAX : (ulong) m;
                        
                        bmap_set_pixel(bms[1], pix, x, y);
                    }
                    
                    if(mps[4]) {
                        ldbl m = esc < 10 || d > 0.5 ? 0 : -logl(esc) / logl(d);
                        m *= ia[3];
                        ulong pix = m > ULONG_MAX ? ULONG_MAX : (ulong) m;
                        
                        bmap_set_pixel(bms[4], pix, x, y);
                    }
                }
            }
            
            if(mps[5]) {
                ldbl v = m2 < 25 ? 0 : ia[4] * logl(dc2 / dz2);
                long pix = v > LONG_MAX ? LONG_MAX : v < LONG_MIN ? LONG_MIN : v;
                
                bmap_set_pixel(bms[5], pix, x, y);
            }
            
            if(mps[6]) {
                if(m2 < 25) {
                    bmap_set_pixel(bms[6], 0, x, y);
                } else {
                    fp80_div(v, dc, dz);
                    ldbl m = atan2l(v->y, v->x);
                    minArg = m < minArg ? m : minArg;
                    maxArg = m > maxArg ? m : maxArg;
                    
                    m *= ia[5];
                    long pix = m > LONG_MAX ? LONG_MAX : m < LONG_MIN ? LONG_MIN : m;
                    
                    bmap_set_pixel(bms[6], pix, x, y);
                }
            }
        }
    }
    
    if(mps[6]) {
        printf("Arg(s) in [%.6Lf, %.6Lf] Pi\n", minArg / PI, maxArg / PI);
    }
    
    bmap *bm = malloc(imgCount * sizeof(bmap));
    int pos = 0;
    for (int i = 0; i < ITERM_COUNT; i++) {
        if(! mps[i]) {
            continue;
        }
        
        bm[pos ++] = bms[i];
    }
    
    return bm;
}

bmap *iterM_bitmaps(drect r, mpc delta, ulong maxIter, int types, int prec, int threads) {
    if(r == NULL || r->w == 0 || r->h == 0 || types == 0 || maxIter <= 1 || prec < 32 ||
       threads < 1 || threads > ITER_MAX_THREADS) {
        return NULL;
    }
    
    if (types == ITERM_ESCAPE) {
        bmap *bm = malloc(sizeof(bmap));
        
        bm[0] = iterM_bitmap(r, delta, maxIter, prec, threads);
        
        return bm;
    }
    
    if(threads > 1) {
        iterBmap it = {.prec = prec, .type = ITER_TYPE_M, .subType = types
            ,.maxIter = maxIter, .r = *r};
        
        if(delta != NULL) {
            mpc_init((mpc_ptr) it.delta, mpc_prec(delta));
            mpc_set((mpc_ptr) it.delta, delta);
        }
        
        bmap *bms = iter_bitmaps(&it, threads);
        if(delta != NULL) {
            mpc_clear((mpc_ptr) it.delta);
        }
        
        return bms;
    }
    
    int mps[] = {(types & ITERM_ESCAPE) != 0, (types & ITERM_DIST) != 0,
        (types & ITERM_DERZ) != 0, (types & ITERM_DERC) != 0, (types & ITERM_ESCPOW) != 0,
        (types & ITERM_SMOD) != 0, (types & ITERM_SARG) != 0
    };
    
    bmap bms[ITERM_COUNT];
    
    int bits[] = {ITER_DIST_BITS, ITER_DERZ_BITS, ITER_DERC_BITS,
        ITER_POW_BITS, ITER_MOD_BITS, ITER_ARG_BITS};
    
    ldbl ia[ITERM_COUNT - 1];
    int imgCount = 0;
    for (int i = 0; i < ITERM_COUNT; i++) {
        if(mps[i]) {
            bms[i] = bmap_new(r);
            bms[i]->type = i >= 5 ? BMAP_TYPE_S_FUNCTION : BMAP_TYPE_MANDEL;
            bms[i]->subType = i + 1;
            
            imgCount ++;
            
            if(i > 0) {
                ia[i - 1] = ldexp(1, bits[i - 1]);
                bms[i]->a = 1 / ia[i - 1];
            }
            
            if(delta != NULL) {
                uint prec = (uint) mpc_prec(delta);
                bms[i]->useHD = true;
                mpfr_inits2(prec, bms[i]->hdx, bms[i]->hdy, NULL);
                mpfr_set(bms[i]->hdx, delta->x, MPFR_RNDN);
                mpfr_set(bms[i]->hdy, delta->y, MPFR_RNDN);
            }
        } else {
            bms[i] = NULL;
        }
    }
    
    if(mps[0]) {
        bms[0]->a /= maxIter + 1;
    }
    
    int cdc = mps[1] || mps[3] || mps[4];    // compute the derivative w.r.t. c
    int cdz = mps[2] || mps[5] || mps[6];    // compute the derivative w.r.t. z
    
    mpc c, v, dz, dc;
    mpc_init(c, prec);
    mpc_init(v, prec);
    
    if(cdz) {
        mpc_init(dz, prec);
    }
    
    if(cdc) {
        mpc_init(dc, prec);
    }
    
    mpfr_t d;
    if(mps[1]) {
        mpfr_init2(d, prec);
    }
    
    defs_mpfr(prec + MPC_EXTRA_PREC, b1, b2, b3, b4, b5);
    
    int w = r->w, h = r->h;
    long lsc = 2 * r->tpow + 10;
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            drect_rel_to_abs_coords(c, r, x, y);
            if(delta != NULL) {
                mpc_add(c, c, delta);
            }
            
            mpc_set(v, c);
            
            if(cdz) {
                mpc_seti(dz, 1, 0);
            }
            
            if(cdc) {
                mpc_seti(dc, 1, 0);
            }
            
            ulong iter = 1;
            for (; iter < maxIter && mpc_modl(v) < 5 && (! cdc || mpc_2exp(dc) < lsc) &&
                    (! cdz || mpc_2exp(dz) > -lsc); iter++) {
                if(cdc) {
                    mpc_mul(dc, dc, v);
                    mpc_scale(dc, dc, 1);
                    mpc_addi(dc, dc, 1);
                }
                
                if(cdz) {
                    mpc_mul(dz, dz, v);
                    mpc_scale(dz, dz, 1);
                }
                                
                mpc_sqr(v, v);
                mpc_add(v, v, c);
            }
            
            bool esc = mpc_modl(v) >= 5;
            bool smallDz = cdz && mpc_modl(dz) < 1;
            if(mps[0]) { // escape time
                ulong pix = esc ? iter : smallDz ? 0 : maxIter + 1;
                
                bmap_set_pixel(bms[0], pix, x, y);
            }
            
            if(mps[2]) { // z derivative
                if(iter <= 1) {
                    bmap_set_pixel(bms[2], 0, x, y);
                } else {
                    ldbl mdz = mpc_modl(dz);
                    ldbl ad = ia[1] * expl(logl(mdz) / (iter - 1));
                    ulong pix = ad > ULONG_MAX ? ULONG_MAX : (ulong) ad;
                    
                    bmap_set_pixel(bms[2], pix, x, y);
                }
            }
            
            if(mps[3]) { // c derivative
                if(iter <= 1) {
                    bmap_set_pixel(bms[3], 0, x, y);
                } else {
                    ldbl mdc = mpc_modl(dc);
                    ldbl ad = ia[2] * expl(logl(mdc) / (iter - 1));
                    ulong pix = ad > ULONG_MAX ? ULONG_MAX : (ulong) ad;
                    
                    bmap_set_pixel(bms[3], pix, x, y);
                }
            }
            
            if(mps[1] || mps[4]) { // distance to M, need to iterate some more
                if(! esc) { // did not escape, distance is 0, mark as max if in the interior
                    if(mps[1]) {
                        bmap_set_pixel(bms[1], smallDz ? 0 : ULONG_MAX / 2 - 1, x, y);
                    }
                    
                    if(mps[4]) {
                        bmap_set_pixel(bms[4], 0, x, y);
                    }
                } else {
                    ulong esc = iter;
                    
                    ldbl m2 = mpc_mod2l(v);
                    for (; m2 < 1E430L; iter++) {
                        mpc_mul(dc, dc, v);
                        mpc_scale(dc, dc, 1);
                        mpc_addi(dc, dc, 1);
                        
                        mpc_sqr(v, v);
                        mpc_add(v, v, c);
                        
                        m2 = mpc_mod2l(v);
                    }
                    
                    int k = (int) (iter - 1);
                    
                    mpc_mod(b5, dc);                      // b5 stores |der|
                    
                    mpc_mod(b1, v);                       // b1 stores |v|
                    mpfr_div(b2, b1, b5, MPFR_RNDN);      // b2 stores |v| / |der|
                    mpfr_mul_2si(b2, b2, -2, MPFR_RNDN);  // b2 stores |v| / 4 |der|
                    
                    mpfr_mul_2si(b1, b1, -1, MPFR_RNDN);  // b1 stores |v| / 2
                    mpfr_log(b3, b1, MPFR_RNDN);          // b3 stores log(|v| / 2)
                    
                    mpfr_mul_2si(b3, b3, -k, MPFR_RNDN);  // b3 stores log(|v| / 2) / 2^k
                    mpfr_expm1(b3, b3, MPFR_RNDN);        // b3 stores exp(log(|v| / 2) / 2^k) - 1
                    mpfr_mul_2si(b3, b3, k, MPFR_RNDN);   // b3 stores 2^k (exp(log(|v| / 2) / 2^k) - 1)
                    
                    mpfr_mul(d, b2, b3, MPFR_RNDN);       // finally, the distance to M
                    
                    if(mps[1]) {
                        mpfr_mul_2si(b2, d, r->tpow, MPFR_RNDN);
                        ldbl sdst = ia[0] / mpfr_get_ld(b2, MPFR_RNDN);
                        ulong pix = sdst >= ULONG_MAX / 2 ? ULONG_MAX / 2 - 1 : (ulong) sdst;
                        
                        bmap_set_pixel(bms[1], pix, x, y);
                    }
                    
                    if(mps[4]) {
                        mpfr_log(d, d, MPFR_RNDN);
                        ldbl m = esc < 10 ? 0 : -logl(esc) / mpfr_get_ld(d, MPFR_RNDN);
                        m *= ia[3];
                        ulong pix = m > ULONG_MAX ? ULONG_MAX : (ulong) m;
                        
                        bmap_set_pixel(bms[4], pix, x, y);
                    }
                }
                
                if(mps[5]) {
                    ldbl mdc = mpc_modl(dc);
                    ldbl mdz = mpc_modl(dz);
                    ldbl m = ia[4] * mdc / mdz;
                    ulong pix = m > ULONG_MAX ? ULONG_MAX : m;
                    
                    bmap_set_pixel(bms[5], pix, x, y);
                }
                
                if(mps[6]) {
                    mpc_div(v, dc, dz);
                    ldbl y = mpfr_get_ld(v->y, MPFR_RNDN);
                    ldbl x = mpfr_get_ld(v->x, MPFR_RNDN);
                    ldbl m = ia[5] * atan2l(y, x);
                    ulong pix = m > ULONG_MAX ? ULONG_MAX : m;
                    
                    bmap_set_pixel(bms[6], pix, x, y);
                }
            }
        }
    }
    
    bmap *bm = malloc(imgCount * sizeof(bmap));
    int pos = 0;
    for (int i = 0; i < ITERM_COUNT; i++) {
        if(! mps[i]) {
            continue;
        }
        
        bm[pos ++] = bms[i];
    }
    
    mpc_clear(c);
    mpc_clear(v);
    
    if(cdz) {
        mpc_clear(dz);
    }
    
    if(cdc) {
        mpc_clear(dc);
    }
    
    if(mps[1]) {
        mpfr_clear(d);
    }
    
    return bm;
}

// MARK: bitmaps of the Julia set

bmap iterJ_bitmap80(fp80 c, drect r, int maxIter, int threads) {
    if(c == NULL || r == NULL || r->w == 0 || r->h == 0 || maxIter <= 1 ||
       threads < 1 || threads > ITER_MAX_THREADS) {
        return NULL;
    }
    
    if(threads > 1) {
        iterBmap it = {.prec = 64, .type = ITER_TYPE_J, .subType = ITERM_ESCAPE
            ,.maxIter = maxIter, .r = *r, .c80 = *c};
        
        bmap *bms = iter_bitmaps(&it, threads);
        if(bms == NULL) {
            return NULL;
        }
        
        bmap bm = bms[0];
        free(bms);
        
        return bm;
    }
    
    bmap bm = bmap_new(r);
    bm->type = BMAP_TYPE_JULIA;
    bm->subType = BMAP_SUB_TYPE_ESCAPE;
    
    fp80 z, dz, b;
    
    int w = r->w, h = r->h;
    ldbl ip2 = ldexpl(100, 2 * r->tpow);
    ldbl p2 = 1 / ip2;
    ip2 *= ip2;
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            drect_rel_to_abs_coords80(z, r, x, y);
            dz->x = 1;
            dz->y = 0;
            
            ulong iter = 0;
            ldbl m2 = fp80_mod2(z);
            ldbl dz2 = 1;
            for (; iter < maxIter && m2 < 25 && dz2 < ip2 && dz2 > p2; iter++) {
                b->x = 2 * z->x;
                b->y = 2 * z->y;
                
                fp80_mul(dz, dz, b);
                
                fp80_sqr(z, z);
                fp80_add(z, z, c);
                
                m2 = fp80_mod2(z);
                dz2 = fp80_mod2(dz);
            }
            
            ulong pix = m2 >= 25 ? iter : dz2 < 1 ? 0 : maxIter + 1;
            bmap_set_pixel(bm, pix, x, y);
        }
    }
    
    return bm;
}

bmap iterJ_bitmap(mpc c, drect r, mpc delta, ulong maxIter, int prec, int threads) {
    if(c == NULL || r == NULL || r->w == 0 || r->h == 0 || prec < 32 || maxIter <= 1 ||
       threads < 1 || threads > ITER_MAX_THREADS) {
        return NULL;
    }
    
    if(threads > 1) {
        iterBmap it = {.prec = prec, .type = ITER_TYPE_J, .subType = ITERM_ESCAPE
            ,.maxIter = maxIter, .r = *r};
        
        mpc_init((mpc_ptr) it.c, mpc_prec(c));
        mpc_set((mpc_ptr) it.c, c);
        
        if(delta != NULL) {
            mpc_init((mpc_ptr) it.delta, mpc_prec(delta));
            mpc_set((mpc_ptr) it.delta, delta);
        }
        
        bmap *bms = iter_bitmaps(&it, threads);
        mpc_clear((mpc_ptr) it.c);
        if(delta != NULL) {
            mpc_clear((mpc_ptr) it.delta);
        }
        
        if(bms == NULL) {
            return NULL;
        }
        
        bmap bm = bms[0];
        free(bms);
        
        return bm;
    }
    
    bmap bm = bmap_new(r);
    bm->type = BMAP_TYPE_JULIA;
    bm->subType = BMAP_SUB_TYPE_ESCAPE;
    
    if(delta != NULL) {
        uint prec = (uint) mpc_prec(delta);
        bm->useHD = true;
        mpfr_inits2(prec, bm->hdx, bm->hdy, NULL);
        mpfr_set(bm->hdx, delta->x, MPFR_RNDN);
        mpfr_set(bm->hdy, delta->y, MPFR_RNDN);
    }
    
    mpc z;
    mpc_init(z, prec);
    
    int w = r->w, h = r->h;
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            drect_rel_to_abs_coords(z, r, x, y);
            if(delta != NULL) {
                mpc_add(z, z, delta);
            }
            
            ulong iter = 0;
            for (; iter < maxIter && mpc_modl(z) < 5; iter++) {
                mpc_sqr(z, z);
                mpc_add(z, z, c);
            }
            
            ulong pix = mpc_modl(z) >= 5 ? iter : maxIter + 1;
            bmap_set_pixel(bm, pix, x, y);
        }
    }
    
    mpc_clear(z);
    
    return bm;
}

bmap *iterJ_bitmaps80(fp80 c, drect r, int maxIter, int types, int threads) {
    if(c == NULL || r == NULL || r->w == 0 || r->h == 0 || types == 0 || maxIter <= 1 ||
       threads < 1 || threads > ITER_MAX_THREADS) {
        return NULL;
    }
    
    if (types == ITERJ_ESCAPE) {
        bmap *bm = malloc(sizeof(bmap));
        
        bm[0] = iterJ_bitmap80(c, r, maxIter, threads);
        
        return bm;
    }
    
    if(threads > 1) {
        iterBmap it = {.prec = 64, .type = ITER_TYPE_J, .subType = types
            ,.maxIter = maxIter, .r = *r, .c80 = *c};
        
        return iter_bitmaps(&it, threads);
    }
    
    int mps[] = {(types & ITERJ_ESCAPE) != 0, (types & ITERJ_DIST) != 0,
        (types & ITERJ_DER) != 0};
    
    bmap bms[ITERJ_COUNT];
    int imgCount = 0;
    for (int i = 0; i < ITERJ_COUNT; i++) {
        if(mps[1]) {
            bms[i] = bmap_new(r);
            bms[i]->type = BMAP_TYPE_JULIA;
            bms[i]->subType = i + 1;
            
            imgCount ++;
        } else {
            bms[i] = NULL;
        }
    }
    
    ldbl iad = ldexp(1, ITER_DIST_BITS);
    if(mps[1]) {
        bms[1]->a = 1 / iad;
    }
    
    ldbl iaz = ldexp(1, ITER_DERZ_BITS);
    if(mps[2]) {
        bms[2]->a = 1 / iaz;
        bms[2]->sgn = 1;
    }
        
    fp80 v, dz, b;
    
    int w = r->w, h = r->h;
    ldbl ip2 = ldexpl(100, 2 * r->tpow);
    ldbl p2 = 1 / ip2;
    ip2 *= ip2;
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            drect_rel_to_abs_coords80(v, r, x, y);
            
            dz->x = 1;
            dz->y = 0;
            
            ulong iter = 0;
            ldbl m2 = fp80_mod2(v);
            ldbl dz2 = 1;
            for (; iter < maxIter && m2 < 25 && dz2 < ip2 && dz2 > p2; iter++) {
                b->x = 2 * v->x;
                b->y = 2 * v->y;
                
                fp80_mul(dz, dz, b);
                
                fp80_sqr(v, v);
                fp80_add(v, v, c);
                
                m2 = fp80_mod2(v);
                dz2 = fp80_mod2(dz);
            }
            
            if(mps[0]) { // escape time
                ulong pix = m2 >= 25 ? iter : dz2 < 1 ? 0 : maxIter + 1;
                
                bmap_set_pixel(bms[0], pix, x, y);
            }
            
            if(mps[2]) { // z derivative
                if(iter <= 1) {
                    bmap_set_pixel(bms[2], 0, x, y);
                } else {
                    ldbl mdz = fp80_mod(dz);
                    ldbl ad = iaz * logl(mdz) / (iter - 1);
                    // works fine if mdz <= 0
                    long pix = ad > LONG_MAX ? LONG_MAX : ad < LONG_MIN ? LONG_MIN : (long) ad;
                    
                    bmap_set_pixel(bms[2], pix, x, y);
                }
            }
            
            if(mps[1]) { // distance to M, need to iterate some more
                if(m2 < 25) { // did not escape, distance is 0
                    bmap_set_pixel(bms[1], 0, x, y);
                } else {
                    for (; m2 < 1E43; iter++) {
                        fp80_mul(dz, dz, v);
                        fp80_muli(dz, dz, 2);
                        
                        fp80_sqr(v, v);
                        fp80_add(v, v, c);
                        
                        m2 = fp80_mod2(v);
                    }
                    
                    int k = (int) (iter - 1);
                    
                    ldbl dm2 = dz->x * dz->x + dz->y * dz->y;
                    ldbl d = sqrtl(m2 / dm2) * 0.25;
                    
                    ldbl lm2p = logl(m2 * 0.25);
                    
                    if(iter < 10000) { // only so that the exponent of ldbl does not overflow
                        d *= ldexpl(expm1l(ldexpl(lm2p, -k - 1)), k);
                    } else {
                        d *= 0.5 * lm2p;
                    }
                    
                    d *= iad;
                    long pix = d > LONG_MAX ? LONG_MAX : (long) d;
                    
                    bmap_set_pixel(bms[1], pix, x, y);
                }
            }
        }
    }
    
    bmap *bm = malloc(imgCount * sizeof(bmap));
    int pos = 0;
    for (int i = 0; i < ITERJ_COUNT; i++) {
        if(! mps[i]) {
            continue;
        }
        
        bm[pos ++] = bms[i];
    }
    
    return bm;
}

bmap *iterJ_bitmaps(mpc c, drect r, mpc delta, ulong maxIter, int types, int prec, int threads) {
    if(c == NULL || r == NULL || r->w == 0 || r->h == 0 || types == 0 ||
       maxIter <= 1 || prec < 32 || threads < 1 || threads > ITER_MAX_THREADS) {
        return NULL;
    }
    
    
    if (types == ITERJ_ESCAPE) {
        bmap *bm = malloc(sizeof(bmap));
        
        bm[0] = iterJ_bitmap(c, r, delta, maxIter, prec, threads);
        
        return bm;
    }
    
    if(threads > 1) {
        iterBmap it = {.prec = prec, .type = ITER_TYPE_J, .subType = types
            ,.maxIter = maxIter, .r = *r};
        
        mpc_init((mpc_ptr) it.c, mpc_prec(c));
        mpc_set((mpc_ptr) it.c, c);
        
        if(delta != NULL) {
            mpc_init((mpc_ptr) it.delta, mpc_prec(delta));
            mpc_set((mpc_ptr) it.delta, delta);
        }
        
        bmap *bms = iter_bitmaps(&it, threads);
        mpc_clear((mpc_ptr) it.c);
        if(delta != NULL) {
            mpc_clear((mpc_ptr) it.delta);
        }
        
        return bms;
    }
    
    int mps[] = {(types & ITERJ_ESCAPE) != 0, (types & ITERJ_DIST) != 0,
        (types & ITERJ_DER) != 0};
    
    bmap bms[ITERJ_COUNT];
    int imgCount = 0;
    for (int i = 0; i < ITERJ_COUNT; i++) {
        if(mps[1]) {
            bms[i] = bmap_new(r);
            bms[i]->type = BMAP_TYPE_JULIA;
            bms[i]->subType = i + 1;
            
            if(delta != NULL) {
                uint prec = (uint) mpc_prec(delta);
                bms[i]->useHD = true;
                mpfr_inits2(prec, bms[i]->hdx, bms[i]->hdy, NULL);
                mpfr_set(bms[i]->hdx, delta->x, MPFR_RNDN);
                mpfr_set(bms[i]->hdy, delta->y, MPFR_RNDN);
            }
            
            imgCount ++;
        } else {
            bms[i] = NULL;
        }
    }
    
    ldbl iad = ldexp(1, ITER_DIST_BITS);
    if(mps[1]) {
        bms[1]->a = 1 / iad;
    }
    
    ldbl iaz = ldexp(1, ITER_DERZ_BITS);
    if(mps[2]) {
        bms[2]->a = 1 / iaz;
    }
    
    mpc z, dz;
    mpc_init(z, prec);
    
    if(mps[2]) {
        mpc_init(dz, prec);
        bms[2]->sgn = 1;
    }
    
    mpfr_t d;
    if(mps[1]) {
        mpfr_init2(d, prec);
    }
    
    defs_mpfr(prec + MPC_EXTRA_PREC, b1, b2, b3, b4, b5);
    
    int w = r->w, h = r->h;
    long lsc = 2 * r->tpow + 10;
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            drect_rel_to_abs_coords(z, r, x, y);
            if(delta != NULL) {
                mpc_add(z, z, delta);
            }
            
            if(mps[2]) {
                mpc_seti(dz, 1, 0);
            }
            
            ulong iter = 0;
            bool esc = false, small = false, large = false;
            for (; iter < maxIter && ! esc && ! small && ! large; iter++) {
                if(mps[2]) {
                    mpc_mul(dz, dz, z);
                    mpc_scale(dz, dz, 1);
                    
                    long sc = mpc_2exp(dz);
                    small = sc < -lsc;
                    large = sc > lsc;
                }
                
                mpc_sqr(z, z);
                mpc_add(z, z, c);
                
                esc = mpc_modl(z) >= 5;
            }
            
            if(mps[0]) { // escape time
                ulong pix = esc ? iter : small ? 0 : maxIter + 1;
                
                bmap_set_pixel(bms[0], pix, x, y);
            }
            
            if(mps[2]) { // z derivative
                if(iter <= 1) {
                    bmap_set_pixel(bms[2], 0, x, y);
                } else {
                    ldbl mdz = mpc_modl(dz);
                    ldbl ad = iaz * logl(mdz) / (iter - 1);
                    // works fine if mdz <= 0
                    long pix = ad > LONG_MAX ? LONG_MAX : ad < LONG_MIN ? LONG_MIN : (long) ad;
                    
                    bmap_set_pixel(bms[2], pix, x, y);
                }
            }
            
            if(mps[1]) { // distance to J, need to iterate some more
                if(! esc) { // did not escape, distance is 0
                    bmap_set_pixel(bms[1], 0, x, y);
                } else {
                    ldbl m2 = 0;
                    for (; m2 < 1E430L; iter++) {
                        mpc_mul(dz, dz, z);
                        mpc_scale(dz, dz, 1);
                        mpc_addi(dz, dz, 1);
                        
                        mpc_sqr(z, z);
                        mpc_add(z, z, c);
                        m2 = mpc_mod2l(z);
                    }
                    
                    int k = (int) (iter - 1);
                    
                    mpc_mod(b5, dz);                      // b5 stores |der|
                    
                    mpc_mod(b1, z);                       // b1 stores |v|
                    mpfr_div(b2, b1, b5, MPFR_RNDN);      // b2 stores |v| / |der|
                    mpfr_mul_2si(b2, b2, -2, MPFR_RNDN);  // b2 stores |v| / 4 |der|
                    
                    mpfr_mul_2si(b1, b1, -1, MPFR_RNDN);  // b1 stores |v| / 2
                    mpfr_log(b3, b1, MPFR_RNDN);          // b3 stores log(|v| / 2)
                    
                    mpfr_mul_2si(b3, b3, -k, MPFR_RNDN);  // b3 stores log(|v| / 2) / 2^k
                    mpfr_expm1(b3, b3, MPFR_RNDN);        // b3 stores exp(log(|v| / 2) / 2^k) - 1
                    mpfr_mul_2si(b3, b3, k, MPFR_RNDN);   // b3 stores 2^k (exp(log(|v| / 2) / 2^k) - 1)
                    
                    mpfr_mul(d, b2, b3, MPFR_RNDN);       // finally, the distance to J
                    
                    ldbl dst = iad * mpfr_get_ld(d, MPFR_RNDN);
                    long pix = dst > LONG_MAX ? LONG_MAX : (long) dst;
                    
                    bmap_set_pixel(bms[1], pix, x, y);
                }
            }
        }
    }
    
    bmap *bm = malloc(imgCount * sizeof(bmap));
    int pos = 0;
    for (int i = 0; i < ITERJ_COUNT; i++) {
        if(! mps[i]) {
            continue;
        }
        
        bm[pos ++] = bms[i];
    }
    
    mpc_clear(c);
    mpc_clear(z);
    
    if(mps[2]) {
        mpc_clear(dz);
    }
    
    if(mps[1]) {
        mpfr_clear(d);
    }
    
    return bm;
}

// MARK: bitmaps of the Newton methods (Julia set and Fatou set components)

// iterates the Newton method of p_n starting at c
// root is used to compute, it may be equal to c
long itern_iter_hyp(fp80 root, fp80 c, int n, int maxIter, ldbl r) {
    if(root != c) {
        *root = *c;
    }
    
    long iter = 1;
    
    fp80 nt;
    mandel_ntl(nt, root, n);
    fp80_sub(root, root, nt);
    
    ldbl r2 = r * r;
    ldbl m2 = fp80_mod2(nt);
    for (; iter < maxIter && m2 >= r2 && isfinite(m2); iter++) {
        if(mandel_ntl(nt, root, n)) {
            fp80_sub(root, root, nt);
        }
        
        m2 = fp80_mod2(nt);
    }
    
    if(m2 < r2) {
        if(mandel_ntl(nt, root, n)) {
            fp80_sub(root, root, nt);
        }
        
        return iter;
    }
    
    return -iter;
}

// find the maximal radius of convergence of the N_{p_n} around c, up to resolution res
ldbl itern_rconvH(fp80 c, int n, ldbl res) {
    ldbl r = 3.01 * ITER_EPSL;                  // convergent radius
    if(mandel_conv_npl(c, n, ITER_EPSL, r)) {
        ldbl dr = 2 * r;                        // divergent radius
        while(mandel_conv_npl(c, n, ITER_EPSL, dr)) {
            r = dr;
            dr *= 2;
        }
        
        while(dr - r > res) {
            ldbl mr = (r + dr) / 2;
            if(mandel_conv_npl(c, n, ITER_EPSL, mr)) {
                r = mr;
            } else {
                dr = mr;
            }
        }
        
        return r;
    }
    
    return 0;
}

static void mark_bmap(bmap bm, int type, fp80 c, ldbl r) {
    ulong pix = type + 1;
    pix <<= ITER_NEWTON_BITS;
    pix ++;
    
    // set the center pixel
    bmap_set_pixel(bm, pix, (int) drect_abs_to_rel_x80(&bm->r, c->x),
                   (int) drect_abs_to_rel_y80(&bm->r, c->y));
    
    long lx = ldexpl(c->x - r, bm->r.tpow) - bm->r.x;
    lx = lx < 0 ? 0 : lx;
    
    long rx = ldexpl(c->x + r, bm->r.tpow) - bm->r.x;
    rx = rx >= bm->r.w ? bm->r.w - 1 : rx;
    
    long ly = ldexpl(c->y - r, bm->r.tpow) - bm->r.y;
    ly = ly < 0 ? 0 : ly;
    
    long hy = ldexpl(c->y + r, bm->r.tpow) - bm->r.y;
    hy = hy >= bm->r.h ? bm->r.h - 1 : hy;
    
    fp80 z;
    
    // the pixel has to be included in the disk
    ldbl rm = r - sqrtl(2) * ldexpl(1, -bm->r.tpow);
    if(rm <= 0) {
        return;
    }
    
    for (int i = (int) lx; i <= rx; i++) {
        for (int j = (int) ly; j <= hy; j++) {
            drect_rel_to_abs_coords80(z, &bm->r, i, j);
            if(fp80_dist(z, c) >= rm) {
                continue;
            }
            
            bmap_set_pixel(bm, pix, i, j);
        }
    }
}

bmap iterN_hyp_bitmap80(int n, drect r, int maxIter, bool quick) {
    if(n <= 1 || n > 100 || r == NULL || r->w <= 0 || r->h <= 0 || maxIter <= 1 ||
       maxIter > MAX_ITERN) {
        return NULL;
    }
    
    long tn = mandel_tot_nt();
    
    bmap bm = bmap_new(r);
    bm->type = BMAP_TYPE_NEWTON_HYP;
    bm->subType = BMAP_SUB_TYPE_CONVERGE;
    bm->typeBits = 64 - ITER_NEWTON_BITS;
    
    fp80 c, nt;
    pset ps;
    pset_init(ps, 2 * ITER_EPSL);
    
    int w = r->w, h = r->h;
    
    // bottom and top edges of the bitmap
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y += h - 1) {
            drect_rel_to_abs_coords80(c, r, x, y);
            long iter = itern_iter_hyp(c, c, n, maxIter, ITER_EPSL);
                        
            if(iter > 0 && drect_contains80(r, c)) {
                pset_add(ps, c);
            }
        }
    }
    
    // left and right edges of the bitmap
    for (int x = 0; x < w; x += w - 1) {
        for (int y = 1; y < h - 1; y++) {
            drect_rel_to_abs_coords80(c, r, x, y);
            long iter = itern_iter_hyp(c, c, n, maxIter, ITER_EPSL);
                        
            if(iter > 0 && drect_contains80(r, c)) {
                pset_add(ps, c);
            }
        }
    }
    
    if(ps->count == 0 || ps->count > 1E6) { // only Julia set points or too many Fatou components
        pset_clear(ps);
        free(bm);
        
        return NULL;
    }
    
    // for each root c of p_n, search for the largest disk around c on which N_{p_n} converges to c
    pset_lock(ps);
        
    double *rs = malloc(ps->count * sizeof(double));
    double res = ldexp(1, -r->tpow - 2);
    for (int i = 0; i < ps->count; i++) {
        pset_point(c, ps, i);
        
        rs[i] = itern_rconvH(c, n, res);
        
        mark_bmap(bm, i, c, rs[i]);
    }
    
    double *tr = malloc(4 * sizeof(double) * maxIter);
    int *it = malloc(2 * sizeof(int) * maxIter);
    int trp = 0;
        
    long skip = 0, partial = 0;
    
    // parse the points in the Fatou set to compute the first entry time (+ 1) in the disk
    // of convergence and to associate the index of the basin of attraction
    fp80 v, root, cent;
    ldbl r2 = ITER_EPSL * ITER_EPSL;
    ulong mask = (1L << ITER_NEWTON_BITS) - 1;
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            ulong pix = bmap_get_pixel(bm, x, y);
            if(pix != 0) {   // already computed
                skip ++;
                
                continue;
            }
            
            drect_rel_to_abs_coords80(v, r, x, y);
            ldbl m2 = 1;
            
            trp = 0;
            int iter;
            ulong p = 0;
            for (iter = 0; iter < 2 * maxIter && m2 >= r2 && isfinite(m2) && p == 0; iter++) {
                if(drect_contains80(r, v)) { // store the trace, relative to the bitmap rectangle
                    fp80_scale(root, v, r->tpow);
                    int ind = trp << 1;
                    
                    tr[ind++] = root->x - r->x;
                    tr[ind] = root->y - r->y;
                    it[trp++] = iter;
                    
                    int xloc = (int) drect_abs_to_rel_x80(r, v->x);
                    int yloc = (int) drect_abs_to_rel_y80(r, v->y);
                    p = bm->pix[xloc + yloc * r->w];

                    if((p & mask) == 1) { // has just entered the marked basin of attraction
                        int pos = (p >> ITER_NEWTON_BITS) - 1;
                        if(pset_point(cent, ps, pos) && fp80_dist(v, cent) < rs[pos]) {
                            partial ++;

                            continue;
                        }
                    } else if(quick && (p != 0)) { // has found a pre-computed point
                        partial ++;

                        continue;
                    }

                    p = 0;
                }
                
                if(mandel_ntl(nt, v, n)) {
                    fp80_sub(v, v, nt);
                    
                    m2 = fp80_mod2(nt);
                } else {
                    m2 = INFINITY;
                }
            }
            
            ulong pi = p & mask;
            if(quick && pi > 1) { // if partially computed, go back and mark all the trace
                ulong ty = p & (~mask);
                if(pi == 2) { // v may be in the marked contraction basin, but marked as time 1  (+ 1)
                    int pos = (p >> ITER_NEWTON_BITS) - 1;
                    if(pset_point(cent, ps, pos) && fp80_dist(v, cent) < rs[pos]) {
                        pi = 1;
                    }
                }
                
                for (int i = 0; i < trp; i++) {
                    ulong pit = pi + (iter - it[i]) - 1;
                    pit = pit > maxIter ? maxIter + 1 : pit;
                    int cx = floor(tr[i << 1]);
                    int cy = floor(tr[(i << 1) + 1]);

                    if(bmap_get_pixel(bm, cx, cy) == 0) {
                        bmap_set_pixel(bm, pit | ty, cx, cy);
                    }
                }

                continue;
            }
            
            int pos = -1;
            if(drect_contains80(r, v)) { // check if it entered a marked basin of attraction
                int xloc = (int) drect_abs_to_rel_x80(r, v->x);
                int yloc = (int) drect_abs_to_rel_y80(r, v->y);
                ulong lp = bm->pix[xloc + yloc * r->w];
                
                if((lp & mask) == 1) {
                    int posg = (lp >> ITER_NEWTON_BITS) - 1;
                    if(pset_point(cent, ps, posg) && fp80_dist(v, cent) < rs[posg]) {
                        pos = posg;
                    }
                }
            }
                
            bool crit = ! isfinite(m2);
            if(pos < 0 && (m2 >= r2 || crit)) { // first maxIter iterates from tr are in the Julia set, mark them
                for (int i = 0; i < trp && (crit || iter - it[i] > maxIter); i++) {
                    int xt = tr[i << 1];
                    int yt = tr[(i << 1) + 1];
                    
                    if(bmap_get_pixel(bm, xt, yt) == 0) {
                        bmap_set_pixel(bm, maxIter + 1, xt, yt);
                    }
                }
                
                continue;
            }
            
            // here c in in the Fatou set; check if it is in the basin of a known root of p_n
            if(pos < 0) {
                pos = (int) pset_index(ps, v);
            }
            
            if(pos < 0) { // unknown basin
                ulong pix = ps->count;
                pix <<= ITER_NEWTON_BITS;
                
                if(bmap_get_pixel(bm, x, y) == 0) {
                    bmap_set_pixel(bm, pix | iter, x, y);
                }
                
                int ti, et;
                for (int i = 0; i < trp; i++) {
                    ti = i << 1;
                    int xt = floor(tr[ti ++]);
                    if(xt < 0 || xt >= r->w) {
                        continue;
                    }
                    
                    int yt = floor(tr[ti]);
                    if(yt < 0 || yt >= r->h) {
                        continue;
                    }
                    
                    et = iter - it[i];
                    
                    if(bmap_get_pixel(bm, xt, yt) == 0) {
                        bmap_set_pixel(bm, pix | et, xt, yt);  // entry time + 1
                    }
                }
                
                continue;
            }
            
            // when the basin is known, go back to find the first entry
            pset_point(root, ps, pos);
            ldbl rs2 = rs[pos] * rs[pos];
            
            // count the first entry time + 1
            ldbl d2 = fp80_dist2(v, root);
            int et = iter > 0 ? iter - 1 : 0, let = et, ti;
            int i = trp - 1;
            for (; i >= 0 && d2 < rs2; i--) {
                ti = i << 1;
                
                c->x = tr[ti ++];
                c->x += r->x;
                c->y = tr[ti];
                c->y += r->y;
                fp80_scale(c, c, -r->tpow);
                
                et = let;
                let = it[i];
                d2 = fp80_dist2(c, root);
            }
            
            if(i < 0 && d2 < rs2) {
                et = let;
            }
            
            et = et == 0 ? 1 : et;
            pix = pos + 1;
            pix <<= ITER_NEWTON_BITS;          // the type is the position in the list
            if(bmap_get_pixel(bm, x, y) == 0) {
                bmap_set_pixel(bm, pix | (et + 1), x, y);  // entry time + 1
            }
            
            for (int i = 1; i < trp && it[i] < et; i++) {
                ti = i << 1;
                int xt = floor(tr[ti ++]);
                if(xt < 0 || xt >= r->w) {
                    continue;
                }
                
                int yt = floor(tr[ti]);
                if(yt < 0 || yt >= r->h) {
                    continue;
                }
                
                let = et - it[i];
                
                if(bmap_get_pixel(bm, xt, yt) == 0) {
                    bmap_set_pixel(bm, pix | (let >= 0 ? let + 1 : 1), xt, yt);  // entry time + 1
                }
            }
        }
    }
    
    free(rs);
    free(tr);
    free(it);
    pset_clear(ps);
    
    tn = mandel_tot_nt() - tn;
    long pix = r->w;
    pix *= r->h;
    printf("Computed %ld Newton terms for the bitmap with %ld pixels, of which %ld skipped and %ld partially skipped.\n",
           tn, pix, skip, partial);
    
    return bm;
}

long itern_iter_mis(fp80 root, fp80 c, int k, int n, int maxIter, ldbl r) {
    if(root != c) {
        *root = *c;
    }
    
    long iter = 1;
    
    fp80 nt;
    mandel_mis_ntl(nt, root, NULL, k, n);
    fp80_sub(root, root, nt);
    
    ldbl r2 = r * r;
    ldbl m2 = fp80_mod2(nt), d2 = 1;
    for (; iter < maxIter && m2 >= r2 && isfinite(m2) && d2 > 1e-10; iter++) {
        d2 = mandel_mis_ntl(nt, root, NULL, k, n);
        if(d2 > 1e-10) {
            fp80_sub(root, root, nt);
            
            m2 = fp80_mod2(nt);
        }
    }
    
    // if it is hyperbolic, it is a multiple root of the mis(k, n), slow convergence
    fp80 th;
    long it = itern_iter_hyp(th, root, n, maxIter, r);
    if(it > 0) {
        ldbl d = fp80_dist(th, root);
        ldbl rc = itern_rconvH(th, n, d / 100);
        if(rc >= d) {
            *root = *th;
            
            return it;
        }
    }
        
    mandel_mis_ntl(nt, root, NULL, k, n);
    fp80_sub(root, root, nt);
    
    return m2 < r2 && d2 > 1e-10 ? iter : -iter;
}

ldbl itern_rconvM(fp80 c, int k, int n, ldbl res) {
    ldbl r = 3.01 * ITER_EPSL;                  // convergent radius
    
    if(mandel_conv_nppl(c, k, n, ITER_EPSL, r, true)) {
        ldbl dr = 2 * r;                        // divergent radius
        while(mandel_conv_nppl(c, k, n, ITER_EPSL, dr, true)) {
            r = dr;
            dr *= 2;
        }
        
        while(dr - r > res) {
            ldbl mr = (r + dr) / 2;
            if(mandel_conv_nppl(c, k, n, ITER_EPSL, mr, true)) {
                r = mr;
            } else {
                dr = mr;
            }
        }
        
        return r;
    }
    
    return 0;
}

bmap iterN_mis_bitmap80(int k, int n, drect r, int maxIter, bool quick) {
    if(k <= 1 || n < 1 || n > 100 || r == NULL || r->w == 0 || r->h == 0 ||
       maxIter <= 1 || maxIter > MAX_ITERN) {
        return NULL;
    }
    
    long tn = mandel_tot_nt();
    
    bmap bm = bmap_new(r);
    bm->type = BMAP_TYPE_NEWTON_MIS;
    bm->subType = BMAP_SUB_TYPE_CONVERGE;
    bm->typeBits = 64 - ITER_NEWTON_BITS;
    
    fp80 c, nt;
    pset ps;
    pset_init(ps, 2 * ITER_EPSL);
    
    int w = r->w, h = r->h;
    
    // bottom and top edges of the bitmap
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y += h - 1) {
            drect_rel_to_abs_coords80(c, r, x, y);
            long iter = itern_iter_mis(c, c, k, n, maxIter, ITER_EPSL);
                        
            if(iter > 0 && drect_contains80(r, c)) {
                pset_add(ps, c);
            }
        }
    }
    
    // left and right edges of the bitmap
    for (int x = 0; x < w; x += w - 1) {
        for (int y = 1; y < h - 1; y++) {
            drect_rel_to_abs_coords80(c, r, x, y);
            long iter = itern_iter_mis(c, c, k, n, maxIter, ITER_EPSL);
                        
            if(iter > 0 && drect_contains80(r, c)) {
                pset_add(ps, c);
            }
        }
    }
    
    if(ps->count == 0 || ps->count > 1E6) { // only Julia set points or too many Fatou components
        pset_clear(ps);
        free(bm);
        
        return NULL;
    }
    
    // for each root c of mis(k, n), search for the largest disk around c on which N_{mis(k, n)} converges to c
    pset_lock(ps);
    
//    printf("Found convergence on the following disks:\n");
    double *rs = malloc(ps->count * sizeof(double));
    double res = ldexp(1, -r->tpow - 2);
    
    for (int i = 0; i < ps->count; i++) {
        pset_point(c, ps, i);
        
        rs[i] = itern_rconvM(c, k, n, res);
//        printf("B((%Lf, %Lf), %g)\n", c->x, c->y, rs[i]);
        
        mark_bmap(bm, i, c, rs[i]);
    }
    
    double *tr = malloc(4 * sizeof(double) * maxIter);
    int *it = malloc(2 * sizeof(int) * maxIter);
    int trp = 0;
        
    long skip = 0, partial = 0;
    
    // parse the points in the Fatou set to compute the first entry time (+ 1) in the disk
    // of convergence and to associate the index of the basin of attraction
    fp80 v, root, cent;
    ldbl r2 = ITER_EPSL * ITER_EPSL;
    ulong mask = (1L << ITER_NEWTON_BITS) - 1;
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            ulong pix = bmap_get_pixel(bm, x, y);
            if(pix != 0) {   // already computed
                skip ++;
                
                continue;
            }
            
            drect_rel_to_abs_coords80(v, r, x, y);
            ldbl m2 = 1;
            
            trp = 0;
            int iter;
            ulong p = 0;
            for (iter = 0; iter < 2 * maxIter && m2 >= r2 && isfinite(m2) && p == 0; iter++) {
                if(drect_contains80(r, v)) { // store the trace, relative to the bitmap rectangle
                    fp80_scale(root, v, r->tpow);
                    int ind = trp << 1;
                    
                    tr[ind++] = root->x - r->x;
                    tr[ind] = root->y - r->y;
                    it[trp++] = iter;
                    
                    int xloc = (int) drect_abs_to_rel_x80(r, v->x);
                    int yloc = (int) drect_abs_to_rel_y80(r, v->y);
                    p = bm->pix[xloc + yloc * r->w];

                    if((p & mask) == 1) { // has just entered the marked basin of attraction
                        int pos = (p >> ITER_NEWTON_BITS) - 1;
                        if(pset_point(cent, ps, pos) && fp80_dist(v, cent) < rs[pos]) {
                            partial ++;

                            continue;
                        }
                    } else if(quick && p != 0) { // has found a pre-computed point
                        partial ++;

                        continue;
                    }

                    p = 0;
                }
                
                if(mandel_mis_ntl(nt, v, NULL, k, n) > 1e-25) {
                    fp80_sub(v, v, nt);
                    
                    m2 = fp80_mod2(nt);
                } else {
                    m2 = INFINITY;
                }
            }
            
            ulong pi = p & mask;
            if(quick && pi > 1) { // if strictly partially computed, go back and mark all the trace
                ulong ty = p & (~mask);
                if(pi == 2) { // v may be in the marked contraction basin, but marked as time 1  (+ 1)
                    int pos = (p >> ITER_NEWTON_BITS) - 1;
                    if(pset_point(cent, ps, pos) && fp80_dist(v, cent) < rs[pos]) {
                        pi = 1;
                    }
                }
                
                for (int i = 0; i < trp - 1; i++) {
                    ulong pit = pi + (iter - it[i]) - 1;
                    pit = pit > maxIter ? maxIter + 1 : pit;
                    int cx = floor(tr[i << 1]);
                    int cy = floor(tr[(i << 1) + 1]);

                    if(bmap_get_pixel(bm, cx, cy) == 0) {
                        bmap_set_pixel(bm, pit | ty, cx, cy);
                    }
                }

                continue;
            }
            
            int pos = -1;
            if(drect_contains80(r, v)) { // check if it entered a marked basin of attraction
                int xloc = (int) drect_abs_to_rel_x80(r, v->x);
                int yloc = (int) drect_abs_to_rel_y80(r, v->y);
                ulong lp = bm->pix[xloc + yloc * r->w];
                
                if((lp & mask) == 1) {
                    int posg = (lp >> ITER_NEWTON_BITS) - 1;
                    if(pset_point(cent, ps, posg) && fp80_dist(v, cent) < rs[posg]) {
                        pos = posg;
                    }
                }
            }
                
            bool crit = ! isfinite(m2);
            if(pos < 0 && (m2 >= r2 || crit)) { // first maxIter iterates from tr are in the Julia set, mark them
                for (int i = 0; i < trp && (crit || iter - it[i] > maxIter); i++) {
                    int xt = tr[i << 1];
                    int yt = tr[(i << 1) + 1];
                    
                    if(bmap_get_pixel(bm, xt, yt) == 0) {
                        bmap_set_pixel(bm, maxIter + 1, xt, yt);
                    }
                }
                
                continue;
            }
            
            // here c in in the Fatou set; check if it is in the basin of a known root of p_n
            if(pos < 0) {
                pos = (int) pset_index(ps, v);
            }
            
            if(pos < 0) { // unknown basin
                ulong pix = ps->count;
                pix <<= ITER_NEWTON_BITS;
                
                if(bmap_get_pixel(bm, x, y) == 0) {
                    bmap_set_pixel(bm, pix | iter, x, y);
                }
                
                int ti, et;
                for (int i = 0; i < trp; i++) {
                    ti = i << 1;
                    int xt = floor(tr[ti ++]);
                    if(xt < 0 || xt >= r->w) {
                        continue;
                    }
                    
                    int yt = floor(tr[ti]);
                    if(yt < 0 || yt >= r->h) {
                        continue;
                    }
                    
                    et = iter - it[i];
                    
                    if(bmap_get_pixel(bm, xt, yt) == 0) {
                        bmap_set_pixel(bm, pix | et, xt, yt);  // entry time + 1
                    }
                }
                
                continue;
            }
            
            // when the basin is known, go back to find the first entry
            pset_point(root, ps, pos);
            ldbl rs2 = rs[pos] * rs[pos];
            
            // count the first entry time + 1
            ldbl d2 = fp80_dist2(v, root);
            int et = iter > 0 ? iter - 1 : 0, let = et, ti;
            int i = trp - 1;
            for (; i >= 0 && d2 < rs2; i--) {
                ti = i << 1;
                
                c->x = tr[ti ++];
                c->x += r->x;
                c->y = tr[ti];
                c->y += r->y;
                fp80_scale(c, c, -r->tpow);
                
                et = let;
                let = it[i];
                d2 = fp80_dist2(c, root);
            }
            
            if(i < 0 && d2 < rs2) {
                et = let;
            }
            
            et = et == 0 ? 1 : et;
            pix = pos + 1;
            pix <<= ITER_NEWTON_BITS;          // the type is the position in the list
            if(bmap_get_pixel(bm, x, y) == 0) {
                bmap_set_pixel(bm, pix | (et + 1), x, y);  // entry time + 1
            }
            
            for (int i = 1; i < trp && it[i] < et; i++) {
                ti = i << 1;
                int xt = floor(tr[ti ++]);
                if(xt < 0 || xt >= r->w) {
                    continue;
                }
                
                int yt = floor(tr[ti]);
                if(yt < 0 || yt >= r->h) {
                    continue;
                }
                
                let = et - it[i];
                
                if(bmap_get_pixel(bm, xt, yt) == 0) {
                    bmap_set_pixel(bm, pix | (let >= 0 ? let + 1 : 1), xt, yt);  // entry time + 1
                }
            }
        }
    }
    
    free(rs);
    free(tr);
    free(it);
    pset_clear(ps);
    
    tn = mandel_tot_nt() - tn;
    long pix = r->w;
    pix *= r->h;
    printf("Computed %ld Newton terms for the bitmap with %ld pixels, of which %ld skipped and %ld partially skipped.\n",
           tn, pix, skip, partial);
    
    return bm;
}

// MARK: highest level of abstraction of this module

static int bit_count(long mask) {
    ulong m = mask;
    int c = 0;
    for (int i = 0; i < 64 && m; i++) {
        c += m & 1;
        m >>= 1;
    }
    
    return c;
}

static void *iter_bmap(void *task, long task_ID, int thread_index) {
    if(task == NULL) {
        return NULL;
    }
    
    iter it = (iter) task;
    switch (it->type) {
        case ITER_TYPE_M:
            if(it->prec <= 64) {
                return iterM_bitmaps80((drect) &it->r, (uint) it->maxIter, it->subType, 1);
            } else {
                return iterM_bitmaps((drect) &it->r,
                                     it->usesDelta ? (mpc_ptr) it->delta : NULL,
                                     it->maxIter, it->subType, it->prec, 1);
            }
            break;
            
        case ITER_TYPE_J:
            if(it->prec <= 64) {
                return iterJ_bitmaps80((fp80_struct *) it->c80, (drect) &it->r,
                                       (uint) it->maxIter, it->subType, 1);
            } else {
                return iterJ_bitmaps((mpc_struct *) it->c, (drect) &it->r,
                                     it->usesDelta ? (mpc_ptr) it->delta : NULL,
                                     it->maxIter, it->subType, it->prec, 1);
            }
            break;
            
        case ITER_TYPE_NH:
            // TODO: reinstate when methods are implemented
//            if(it->prec <= 64) {
//                return iterN_hypBitmaps80(it->per, (drect) &it->r,
//                                          (uint) it->maxIter, it->subType, 1);
//            }
//            else {
//                return iterN_hypBitmaps(it->per, (drect) &it->r, it->maxIter, it->subType, it->prec);
//            }
            break;
            
        case ITER_TYPE_NM:
            // TODO: reinstate when methods are implemented
//            if(it->prec <= 64) {
//                return iterN_misBitmaps80(it->prePer, it->per, (drect) &it->r,
//                                          (uint) it->maxIter, it->subType);
//            } else {
//                return iterN_misBitmaps(it->prePer, it->per, (drect) &it->r,
//                                        it->maxIter, it->subType, it->prec);
//            }
            break;
    }
    
    return NULL;
}

bmap *iter_bitmaps(iter it, int threads) {
    bool pok = it != NULL && threads > 0 && threads <= ITER_MAX_THREADS;
    pok = pok && it->type > 0 && it->type <= ITERT_COUNT && it->r.w > 0 && it->r.h > 0;
    
    if(pok && threads == 1) {
        return (bmap *) iter_bmap(it, 0, 0);
    }
    
    int bmc = bit_count(it->subType);
    bmap *bm = malloc(sizeof(bmap) * bmc);
    if(bm == NULL) {
        return NULL;
    }
        
    bool ok = true;
    for (int i = 0; i < bmc && ok; i++) {
        bm[i] = bmap_new((drect) &it->r);
        ok = ok && bm[i] != NULL;
        if(! ok) {
            continue;
        }
                
        if(it->usesDelta) {
            uint prec = (uint) mpc_prec(it->delta);
            bm[i]->useHD = true;
            mpfr_inits2(prec, bm[i]->hdx, bm[i]->hdy, NULL);
            mpfr_set(bm[i]->hdx, (mpfr_ptr) it->delta->x, MPFR_RNDN);
            mpfr_set(bm[i]->hdy, (mpfr_ptr) it->delta->y, MPFR_RNDN);
        }
    }
    
    if(! ok) { // not enough memory?
        for (int i = 0; i < bmc && ok; i++) {
            if(bm[i] != NULL) {
                bmap_free(bm[i]);
            }
        }
        
        free(bm);
        
        return NULL;
    }
    
    pok = pok && mth_init(threads, &iter_bmap, 10, 10);
    if(! pok) {
        return NULL;
    }
    
    iterBmap tasks[threads];
    for (int i = 0; i < threads; i++) {
        tasks[i] = *it;
    }
    
    int tx = 7 * sqrt(threads);
    int ty = tx;
    tx = tx > it->r.w ? it->r.w : tx;
    ty = ty > it->r.h ? it->r.h : ty;
    int t = tx * ty;
    
    // submit initial jobs for all threads
    drect r = (drect) &it->r;
    for (int i = 0; i < threads && ok; i++) {
        long l = i / tx;
        long c = i % tx;
        
        long px = c * r->w / tx;
        tasks[i].r.x = r->x + px;
        tasks[i].r.w = (uint) ((c + 1) * r->w / tx - px);
        
        long py = l * r->h / ty;
        tasks[i].r.y = r->y + py;
        tasks[i].r.h = (uint) ((l + 1) * r->h / ty - py);
        
        ok = ok && mth_add_task(i, tasks + i);
    }
    
    ulong ID;
    bmap *res = NULL;
    iter task = NULL;
    int done = 0, sent = threads;
    while(ok && done < t) {
        while(mth_errors() == 0 && mth_results() == 0) {
            mth_wait();
        }
        
        ok = ok && mth_errors() == 0;
        ok = ok && mth_get_result(&ID, (void **) &res, (void **) &task);
        
        for (int i = 0; i < bmc && ok; i++) {
            ok = ok && bmap_draw(bm[i], res[i]);
            
            bm[i]->a = res[i]->a;
            bm[i]->b = res[i]->b;
            
            bm[i]->type = res[i]->type;
            bm[i]->subType = res[i]->subType;
            
            bm[i]->pixelType = res[i]->pixelType;
            bm[i]->sgn = res[i]->sgn;
            bm[i]->zeroTransp = res[i]->zeroTransp;
            bm[i]->colorMap = res[i]->colorMap;
            bm[i]->power = res[i]->power;
            bm[i]->mapA = res[i]->mapA;
            bm[i]->mapB = res[i]->mapB;
            bm[i]->colorLow = res[i]->colorLow;
            bm[i]->colorLast = res[i]->colorLast;
            bm[i]->colorHigh = res[i]->colorHigh;
            bm[i]->colorFirst = res[i]->colorFirst;
            
            bmap_free(res[i]);
        }
        free(res);
        
        if(! ok) {
            continue;
        }
        
        done ++;
        
        if(sent < t) {
            long l = sent / tx;
            long c = sent % tx;
            
            long px = c * r->w / tx;
            task->r.x = r->x + px;
            task->r.w = (uint) ((c + 1) * r->w / tx - px);
            
            long py = l * r->h / ty;
            task->r.y = r->y + py;
            task->r.h = (uint) ((l + 1) * r->h / ty - py);
            
            ok = ok && mth_add_task(sent ++, task);
        }
    }
    
    if(! ok || ! mth_end_threads(false, false)) {
        mth_kill_threads();
        for (int i = 0; i < bmc && ok; i++) {
            bmap_free(bm[i]);
        }
        
        free(bm);
        
        return NULL;
    }
    
    return bm;
}
