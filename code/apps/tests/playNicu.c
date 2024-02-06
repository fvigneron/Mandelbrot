//
//  playNicu.c
//  Mandel
//
//  Created by MIHALACHE Nicolae on 12/24/22.
//  Copyright © 2022 MIHALACHE Nicolae. All rights reserved.
//

#include "playNicu.h"

#include "misSets.h"
#include "misRaw.h"
#include "nSet.h"
#include "u128c.h"
#include "nset2csv.h"
#include "misSimpleSets.h"
#include "misRawCount.h"
#include "hypRawCount.h"
#include "polynomial.h"
#include "arithmetic.h"
#include "treeMap.h"
#include "iterates.h"
#include "stopWatch.h"
#include "levSets.h"
#include "io.h"
#include "planarSet.h"
#include "mthread.h"

int searchLevm(int pp, int per, char *x, char *y) {
    levm lc = miss_load(pp, per);
    
    mpc p, z, d;
    int prec = lc->prec;
    mpc_inits(prec, p, z, d, NULL);
    mpfr_set_str(z->x, x, 10, MPFR_RNDN);
    mpfr_set_str(z->y, y, 10, MPFR_RNDN);
    
    fp80 d8;
    int pos = -1, npos = -1, nnpos = -1;
    ldbl minDist = 5, nDist = 5, nnDist = 5;
    for (int i = 0; i < levm_count(lc); i++) {
        levm_point(p, lc, i);
        mpc_sub(d, p, z);
        mpc_get80(d8, d);
        
        ldbl m = fp80_mod(d8);
        if(m < minDist) {
            nnpos = pos;
            nnDist = nDist;
            npos = pos;
            nDist = minDist;
            
            pos = i;
            minDist = m;
        } else if(m < nDist) {
            nnpos = npos;
            nnDist = nDist;
            
            npos = i;
            nDist = m;
        } else if(m < nnDist) {
            nnpos = i;
            nnDist = m;
        }
    }
    
    levm_point(p, lc, pos);
    mpfr_printf("The closest point to (%s, %s) is (%.36Rf, %.36Rf), on position %d, at distance %.3Lg\n",
           x, y, p->x, p->y, pos, minDist);
    
    levm_point(p, lc, npos);
    mpfr_printf("The next closest point is (%.36Rf, %.36Rf), on position %d, at distance %.3Lg\n",
           p->x, p->y, npos, nDist);
    
    levm_point(p, lc, nnpos);
    mpfr_printf("The next closest point is (%.36Rf, %.36Rf), on position %d, at distance %.3Lg\n\n",
           p->x, p->y, nnpos, nnDist);
    
    levm_free(lc);
    mpc_clears(p, z, d, NULL);
    
    return pos;
}

ldbl levmDist(int pp, int per, int pos, char *x, char *y) {
    levm lc = miss_load(pp, per);
    
    mpc p, z, d;
    int prec = lc->prec;
    mpc_inits(prec, p, z, d, NULL);
    mpfr_set_str(z->x, x, 10, MPFR_RNDN);
    mpfr_set_str(z->y, y, 10, MPFR_RNDN);
    
    fp80 d8;
    
    levm_point(p, lc, pos);
    mpc_sub(d, p, z);
    mpc_get80(d8, d);
    
    ldbl m = fp80_mod(d8);
    
    mpfr_printf("The distance from (%s, %s) to (%.26Rf, %.26Rf), on position %d, is %.5Lg\n",
           x, y, p->x, p->y, pos, m);
    
    levm_free(lc);
    mpc_clears(p, z, d, NULL);
    
    return m;
}

void testMisRaw(int pp, int per, int seg, char *x, char *y, char *x1, char *y1) {
    mpc c;
    u128 u;
    int prec = 170;
    mpc_init(c, prec);
    
    nset mis = mis_raw_mini(pp, per, seg);
    
    mpfr_set_str(c->x, x, 10, MPFR_RNDN);
    mpfr_set_str(c->y, y, 10, MPFR_RNDN);
    
    mis_refine(c, c, pp, per, prec / 3.2);
    
    u128_set(u, c);
    if(nset_contains(mis, u)) {
        mpfr_printf("The parameter (%.40Rf, %.40Rf) was found.\n", c->x, c->y);
    } else {
        mpfr_printf("The parameter (%.40Rf, %.40Rf) was NOT found.\n", c->x, c->y);
    }
    
    mpfr_set_str(c->x, x1, 10, MPFR_RNDN);
    mpfr_set_str(c->y, y1, 10, MPFR_RNDN);
    
    mis_refine(c, c, pp, per, prec / 3.2);
    
    u128_set(u, c);
    if(nset_contains(mis, u)) {
        mpfr_printf("The parameter (%.40Rf, %.40Rf) was found.\n", c->x, c->y);
    } else {
        mpfr_printf("The parameter (%.40Rf, %.40Rf) was NOT found.\n", c->x, c->y);
    }
    
    nset_free(mis);
    mpc_clear(c);
}

void mis_ref(int pp, int per, int prec, char *x, char *y) {
    mpc c;
    mpc_init(c, prec);
    
    mpfr_set_str(c->x, x, 10, MPFR_RNDN);
    mpfr_set_str(c->y, y, 10, MPFR_RNDN);
    
    int dig = prec / 3.2;
    if(mis_refine(c, c, pp, per, dig)) {
        char fmt[130];
        snprintf(fmt, 130, "Mis(%d, %d) : (%%.%dRf, %%.%dRf)\n", pp, per, dig, dig);
        mpfr_printf(fmt, c->x, c->y);
    } else {
        printf("Could not refine (%s, %s) as a Mis(%d, %d) !\n", x, y, pp, per);
    }
}
                 
char *x[] = {"-1.94158865210904662057070991168542614663279061760593011",
             "-1.94158865210904660328362032237132607619084190244713771",
             "-1.941588652109046620302975019289018045"};

char *y[] = {"0.00767776557696507453523590425645090679750852806464792",
             "0.00767776557696507757432207556589471308992779018356888",
             "0.007677765576965153978095324849604831"};

mpc_ptr mpc_str(int prec, char *x, char *y) {
    mpc_ptr c = malloc(sizeof(mpc_struct));
    mpc_init(c, prec);
    
    mpfr_set_str(c->x, x, 10, MPFR_RNDN);
    mpfr_set_str(c->y, y, 10, MPFR_RNDN);
    
    return c;
}

mpc_ptr mpc_zero(int prec) {
    mpc_ptr c = malloc(sizeof(mpc_struct));
    mpc_init(c, prec);
    
    mpfr_set_zero(c->x, 1);
    mpfr_set_zero(c->y, 1);
    
    return c;
}

void mis_root(int pp, int per, int prec, char *x, char *y) {
    mpc_ptr c = mpc_str(prec, x, y);
    mpc_ptr nt = mpc_zero(prec);
    
    mpfr_printf("%.30Rf, %.30Rf\n", c->x, c->y);
    do {
        mandel_mis_nt(nt, c, pp, per);
        mpc_sub(c, c, nt);
        mpfr_printf("%.30Rf, %.30Rf\n", c->x, c->y);
    } while(-mpc_2exp(nt) < 3 * prec / 4);
    
    mandel_mis_nt(nt, c, pp, per);
    mpc_sub(c, c, nt);
    mpfr_printf("%.30Rf, %.30Rf\n", c->x, c->y);
    
    printf("\n");
    mpc_set(nt, c);
    for (int i = 0; i < pp + per; i++) {
        mpc_sqr(nt, nt);
        mpc_add(nt, nt, c);
        mpfr_printf("%.30Rf, %.30Rf\n", nt->x, nt->y);
    }
    
    mpc_clears(c, nt, NULL);
    free(c);
    free(nt);
}

void mis_lc(int pp, int per, int prec, char *xc, char *yc) {
    mpc_ptr c = mpc_str(prec, xc, yc);
    mpc_ptr v = mpc_zero(prec);
    mpc_ptr t = mpc_zero(prec);
    mpc_ptr r = mpc_zero(prec);
    
    mpc_set(v, c);
    mpfr_printf("%.30Rf, %.30Rf\n", v->x, v->y);
    
    mpc_set_str(c, x[0], y[0]);
    
    mpfr_t th;
    mpfr_init2(th, prec);
    
    for (int i = 0; i < 30; i++) {
        mpfr_set_d(th, 0.25 * (i % 8), MPFR_RNDN);
        mpc_exp_2Pi_i(t, th);
        mpc_muli(t, t, 100);
        
        mandel_mis_root_ref(v, v, t, pp, per, 50, 1);
        
        mpc_sub(t, v, c);
        ldbl d = mpc_modl(t);
        
        mpc_set0(t);
        mandel_mis_root_ref(r, v, t, pp, per, 150, 1);
        mpc_sub(t, r, c);
        ldbl dr = mpc_modl(t);
        
        mpfr_printf("%.30Rf, %.30Rf; d = %.4Lg -> %.30Rf, %.30Rf; d = %.4Lg\n",
                    v->x, v->y, d, r->x, r->y, dr);
    }
        
    mpfr_clear(th);
    mpc_clears(c, v, t, r, NULL);
    free(c);
    free(v);
    free(t);
    free(r);
}

void testMisSimpleSet(int pp, int per) {
    levm lc = missi_load(pp, per);
    int prec = lc->prec;
    
    mpc c, v, r;
    mpc_inits(prec, c, v, r, NULL);
    mpc_setr(r, lc->radius);
    
    ldbl md = 0;
    
    for (int i = 0; i <= levm_segs(lc); i++) {
        levm_point(c, lc, i);
        mandel_mis_val(v, c, pp, per);
        mandel_mis_val(c, c, pp - 1, per);
        
        mpc_div(v, v, c);
        ldbl d = mpc_distl(v, r);
        
        if(d > md) {
            md = d;
            
            printf("i = %d, dist = %.4Lg\n", i, d);
        }
    }
    
    mpc_clears(c, v, r, NULL);
}

void testHypLc(void) {
    char *spx = "0.288078790747107486313150279327109954737";
    mpc c, v;
    mpc_inits(120, c, v, NULL);
    
    mpc_set_str(c, spx, "0");
    mandel_val(v, c, 15);
    
    mpfr_printf("%10Rg %10Rg\n", v->x, v->y);
    
    mpc_clears(c, v, NULL);
}

void countMultCrit(void) {
    for (int per = 2; per < 42; per++) {
        ulong cc = mandel_mult_crit_count(per);
        printf("%d, %lu\n", per, cc);
    }
}

void helpInvPows(void) {
    printf("Use params: pre-period period pows first count error\n");
}

void testPolyInvPows(void) {
    int prec = 200;
    __mpfr_struct irl[2];
    mpfr_inits2(prec, irl, irl + 1, NULL);
    mpfr_set_si(irl, 0, MPFR_RNDN);
    mpfr_set_si(irl + 1, -8, MPFR_RNDN);
    
    poly P;
    poly_init_root_powers_real(P, 2, irl, prec, 1e-35);
    
    for (int i = 0; i <= P->deg; i++) {
        mpfr_printf("%.20Rg\n", P->a[i].x);
    }
}

void mpvDiff(char *fn1, char *fn2, char *ofn) {
    mpv a = mpv_read(fn1, false);
    mpv b = mpv_read(fn2, false);
    
    int n = (int) (a->count < b->count ? a->count : b->count);
    int prec = a->prec;
    
    mpfr_t d, buf;
    mpfr_init2(d, prec);
    mpfr_init2(buf, prec);
    
    mpv diff = mpv_new(prec, n);
    
    for (int i = 0; i < n; i++) {
        mpv_get(d, a, i);
        mpv_get(buf, b, i);
        
        mpfr_sub(d, d, buf, MPFR_RNDN);
        mpv_set(diff, i, d);
    }
    
    mpfr_clear(d);
    mpfr_clear(buf);
    
    if(mpv_write(diff, ofn, false)) {
        printf("Difference written to %s\n", ofn);
    } else {
        printf("Could not write the difference to %s\n", ofn);
    }
    
    mpv_free(a);
    mpv_free(b);
    mpv_free(diff);
}

void nsetDiff(char *fn, char *fd) {
    printf("Searching for points in %s not contained in %s:\n", fn, fd);
    
    nset lns = nset_load(fn, true);
    nset_lock(lns);
    printf("%ld points read from %s\n", lns->count, fn);
    
    nset sns = nset_load(fd, true);
    nset_lock(sns);
    printf("%ld points read from %s\n", sns->count, fd);
    
    u128 c;
    mpc mc;
    mpc_init(mc, 128);
    for (int i = 0; i < lns->count; i++) {
        nset_point(c, lns, i);
        
        if((i + 1) % 1000000 == 0) {
            printf((i + 1) % 10000000 == 0 ? ". " : ".");
        }
        
        if(nset_contains(sns, c)) {
            continue;
        }
        
        u128_get(mc, c);
        mpfr_printf("%.30Rg, %.30Rg\n", mc->x, mc->y);
    }
    
    mpc_clear(mc);
}

void nsetSearchNset(int pp, int per, char *fn) {
    printf("Searching for points in %s contained in some raw set of type Mis(%d, %d):\n",
           fn, pp, per);
    
    nset ns = nset_load(fn, true);
    nset_lock(ns);
    long n = ns->count;
    printf("%ld points read from %s\n", n, fn);
    
    int fc = mis_raw_jobs_count(pp, per);
    char rawfn[120];
    
    mpc p;
    mpc_init(p, 128);
    int files = 0;
    
    for (int i = 0; i < fc; i++) {
        mis_raw_job_file_name(rawfn, 120, pp, per, i);
        nset ps = nset_load(rawfn, true);
        if(ps == NULL) {
            continue;
        }
        
        printf("\nLoaded %ld points from %s\n", ps->count, rawfn);
        files ++;
        
        for (int j = 0; j < n; j++) {
            u128 c;
            nset_point(c, ns, j);
            
            if(! nset_contains(ps, c)) {
                continue;
            }
            
            u128_get(p, c);
            mpfr_printf("%.36Rg, %.36Rg was found in this set\n", p->x, p->y);
        }
        
        nset_free(ps);
    }
    
    mpc_clear(p);
    nset_free(ns);
    
    if(files == 0) {
        printf("\nCould not find any misRaw file !\n\n");
    } else {
        printf("\nSearched for the %ld points in %d file.\n\n", n, files);
    }
}

void nsetSearchMpvInNSet(char *mpvFn, char *nSetFn) {
    printf("Searching for points in %s contained in %s:\n", mpvFn, nSetFn);
    
    nset ns = nset_load(nSetFn, true);
    nset_lock(ns);
    long m = ns->count;
    printf("%ld points read from %s\n", m, nSetFn);
    
    mpv v = mpv_read(mpvFn, false);
    long n = v->count / 2;
    printf("%ld points read from %s\n", n, mpvFn);
    
    defs_mpc(v->prec, c, a, b);
    for (int i = 0; i < n; i++) {
        mpv_getc(c, v, i);
        u128 u;
        u128_set(u, c);
        if(nset_contains(ns, u)) {
            int pos = -1;
            ldbl minDist = 5;
            
            for (int j = 0; j < m; j++) {
                nset_get(a, ns, j);
                ldbl d = mpc_distl(c, a);
                
                if(d < minDist) {
                    pos = j;
                    minDist = d;
                }
            }
            
            nset_get(a, ns, pos);
            
            printf("The %d-th point int the mpv ", i);
            mpc_print(c, 40);
            printf("\nis closest to the %d-th of the nset ", pos);
            mpc_print(a, 36);
            printf("\ndistance = %Lg\n\n", minDist);
        }
    }
    
    printf("All the other points in the mpv are not containd in the nset.\n");
    
    nset_free(ns);
    mpv_free(v);
}

bool refineAndSave(int pp, int per, mpv roots, char *fn) {
    int prec = roots->prec;
    int n = (int) roots->count / 2;
    
    mpc c;
    mpc_init(c, prec);
    u128 u;
    
    nset_t ns;
    nset_init(ns, MIS_RAW_SET_EPS);
    
    for (int i = 0; i < n; i++) {
        mpv_getc(c, roots, i);
        mpfr_printf("\n%.40Rg, %.40Rg refined to\n", c->x, c->y);
        mandel_mis_root_ref(c, c, NULL, pp, per, 10, 1);
        mpfr_printf("%.40Rg, %.40Rg\n", c->x, c->y);
        
        u128_set(u, c);
        nset_add(ns, u);
    }
    
    nset_lock(ns);
    bool ok = nset_write(ns, fn);
    if(ok) {
        printf("\n%ld refined parameters written to %s, of which %ld are real\n", ns->count, fn, ns->realCount);
    } else {
        printf("\nCould not write the refined parameters to %s\n", fn);
    }
    
    nset_clear(ns);
    
    return ok;
}

void tmapTest(int pp, int per) {
    char fn[100];
    snprintf(fn, 100, "mis/%02d-%02d/mis%02d-%02d.tmap", pp, per, pp, per);
    
    tmap t = tmap_load(fn, 0, true, NULL);
    if(t == NULL) {
        printf("Could not load the theeMap from %s\n", fn);
        
        return;
    }
    
    printf("Loaded the treeMap with %ld points\n", t->count);
    drect_t d = {{.x = -2, .y = 0, .w = 4, .h = 2, .tpow = 0}};
    
    long c = tmap_count_abs(t, d);
    drect_print(d, fn, 100);
    printf("There are %ld points in the rectangle %s\n", c, fn);
}

void resetEpsilon(int pp, int per, bool raw) {
    int count = raw ? mis_raw_jobs_count(pp, per) : mis_results_count(pp, per);
    printf("Setting epsilon to %Lg for %s nset files of type (%d, %d):\n",
           ldexpl(1, -126 + MIS_RAW_SET_2POW), raw ? "raw" : "final", pp, per);
    
    char fn[200];
    for (int i = 0; i < count; i++) {
        if(raw) {
            mis_raw_job_file_name(fn, 200, pp, per, i);
        } else {
            mis_file_name(fn, 200, pp, per, i);
        }
        
        nset ns = nset_load(fn, true);
        if(ns == NULL || MIS_RAW_SET_EPS >= ns->eps) {
            if(ns != NULL) {
                printf("Loaded %s, no need for changes.\n", fn);
                
                nset_free(ns);
            } else {
                printf("Could not load %s !\n", fn);
            }
                
            continue;
        }
        
        nset_set_eps(ns, MIS_RAW_SET_EPS);
        nset_write(ns, fn);
        printf("Changed epsilon in %s.\n", fn);
        
        nset_free(ns);
    }
}

static double colBeta = 20;
uint col_n(bmap bm, ulong iter) {
    if(iter == 0) {
        return 0;
    }
    
    int tb = bm->typeBits;
//    int type = (int) (iter >> (64 - tb));
    long it = (iter << tb) >> tb;
    
    if(it == 1) {
        return 0xFF0000;
    }
    
    long cl = 255 - colBeta * (it - 2);
    uint c = cl < 20 ? 0 : (uint) cl;
    
    return c | (c << 8) | (c << 16);
}

void compImg(int per, int iter, int size, int twop) {
    drect_t r = {-size / 2, -size / 2, size, size, twop};
    
    bmap bmq = iterN_hyp_bitmap80(per, r, iter, true);
    bmap bm = iterN_hyp_bitmap80(per, r, iter, false);
    
    for (int y = 0; y < size; y++) {
        for (int x = 0; x < size; x++) {
            ulong pq = bmap_get_pixel(bmq, x, y);
            ulong p = bmap_get_pixel(bm, x, y);
            
            if(p != pq) {
                printf("at (%d, %d), bmq = %lX and bm = %lX\n", x, y, pq, p);
                
                return;
            }
        }
    }
}

void testImg(int pp, int per, int iter, int size, int dx, int dy, int twop, bool quick) {
    drect_t r = {-size / 2 + dx, -size / 2 + dy, size, size, twop};
    char pr[100];
    drect_print(r, pr, 100);
    printf("Computing the image of the Newton method in the rectangle %s\n", pr);

    if(pp == 0) {
        printf("of p_%d, with at most %d iterates\n", per, iter);
    } else {
        printf("of mis(%d, %d), with at most %d iterates\n", pp, per, iter);
    }
    if(quick) {
        printf("Using a quick and dirty method !\n\n");
    }
    
    bmap bmn = NULL;
    char fn[100];
    if(pp <= 0) {
        bmn = iterN_hyp_bitmap80(per, r, iter, quick);
        if(bmn == NULL) {
            printf("Could not compute the image, too few iterates?\n\n");
            
            return;
        }
        
        printf("Image rendered.\n");
        
        if(dx == 0 && dy == 0) {
            sprintf(fn, quick ? "newton%d_%d_%dq.bmap"
                    : "newton%d_%d_%d.bmap", per, twop, size);
        } else {
            sprintf(fn, quick ? "newton%d_%d_%d_%d_%dd.bmap"
                    : "newton%d_%d_%d_%d_%d.bmap", per, dx, dy, twop, size);
        }
        bmap_write(bmn, fn, false);
        
        if(dx == 0 && dy == 0) {
            sprintf(fn, quick ? "newton%d_%d_%dq.bmp"
                    : "newton%d_%d_%d.bmp", per, twop, size);
        } else {
            sprintf(fn, quick ? "newton%d_%d_%d_%d_%dq.bmp"
                    : "newton%d_%d_%d_%d_%d.bmp", per, dx, dy, twop, size);
        }
    } else {
        bmn = iterN_mis_bitmap80(pp, per, r, iter, quick);
        if(bmn == NULL) {
            printf("Could not compute the image, too few iterates?\n\n");
            
            return;
        }
        
        printf("Image rendered.\n");
        
        if(dx == 0 && dy == 0) {
            sprintf(fn, quick ? "newtonMis%d_%d_%d_%dq.bmap"
                    : "newtonMis%d_%d_%d_%d.bmap", pp, per, twop, size);
        } else {
            sprintf(fn, quick ? "newtonMis%d_%d_%d_%d_%d_%dq.bmap"
                    : "newtonMis%d_%d_%d_%d_%d_%d.bmap", pp, per, dx, dy, twop, size);
        }
        bmap_write(bmn, fn, false);
        
        if(dx == 0 && dy == 0) {
            sprintf(fn, quick ? "newtonMis%d_%d_%d_%dq.bmp"
                    : "newtonMis%d_%d_%d_%d.bmp", pp, per, twop, size);
        } else {
            sprintf(fn, quick ? "newtonMis%d_%d_%d_%d_%d_%dq.bmp"
                    : "newtonMis%d_%d_%d_%d_%d_%d.bmp", pp, per, dx, dy, twop, size);
        }
    }
    
    if(bmap_write_bmp(bmn, fn, &col_n)) {
        printf("Written the image to %s.\n", fn);
    } else {
        printf("Could not write the image to %s.\n", fn);
    }
    
    free(bmn);
}

bool scurve(levc lc) {
    int per = lc->period - 4;
    
    struct timeb te;
    char time[80], fn[100];
    
    printf("    Refining the level curve of period %u: ", lc->period);
    fflush(stdout);
    ftime(&te);
    
    lc = levc_next_period(lc, true);
    levc tlc = lc;
    printf(".");
    
    while(lc != NULL && lc->radius->_mpfr_exp < 125) {
        tlc = levc_next_period(lc, true);
        levc_free(lc);
        lc = tlc;
        
        printf(".");
        fflush(stdout);
    }
    
    lapse(&te, time);
    printf(" done in %s\n", time);
    fflush(stdout);
    
//    defs_mpc(LEVS_PREC, c, v, dz, idz, s, ls);
    mpc c, v, dz, idz, s, ls;
    mpc_inits(LEVS_PREC, c, v, dz, idz, s, ls, NULL);
    
    ulong n = levc_count(lc);
    mpv sc = mpv_new(LEVS_PREC, 2 * n);
    mpv lsc = mpv_new(LEVS_PREC, 2 * n);
    
    printf("    Computing the S-curve %u: ", lc->period);
    fflush(stdout);
    ftime(&te);
    
    for (int i = 0; i < n; i++) {
        levc_point(c, lc, i);
        
        // compute s = S(c)
        mpc_set0(v);
        mpc_seti(s, 1, 0);
        mpc_seti(dz, 1, 0);
        
        do {
            mpc_sqr(v, v);
            mpc_add(v, v, c);
            
            mpc_mul(dz, dz, v);
            mpc_scale(dz, dz, 1);
            
            mpc_inv(idz, dz);
            mpc_add(s, s, idz);
        } while(mpc_2exp(dz) - mpc_2exp(s) < 125);
        
        // store the value and its Log
        mpv_setc(sc, i, s);
        
        mpc_mod2(ls->x, s);
        mpfr_log(ls->x, ls->x, MPFR_RNDN);
        mpfr_mul_2si(ls->x, ls->x, -1, MPFR_RNDN);
        
        mpfr_atan2(ls->y, s->y, s->x, MPFR_RNDN);
        mpv_setc(lsc, i, ls);
        
        if(i > 0 && i % 1024 == 0) {
            printf(i / 1024 == 8 ? ". " : ".");
            fflush(stdout);
        }
    }
    
    lapse(&te, time);
    printf(" done in %s\n", time);
    fflush(stdout);
    
    int digits = 20;
    snprintf(fn, 100, "levc_%d_%d.csv", per, digits);
    if(mpv_write_csv(lc->points, fn, true, digits, 0, n, false)) {
        printf("    Written the level curve to %s\n", fn);
    } else {
        printf("****Could not write the level curve to %s !!!\n", fn);
    }
    
    snprintf(fn, 100, "S_%d_%d.csv", per, digits);
    if(mpv_write_csv(sc, fn, true, digits, 0, n, false)) {
        printf("    Written the S-curve to %s\n", fn);
    } else {
        printf("****Could not write the S-curve to %s !!!\n", fn);
    }
    
    snprintf(fn, 100, "LogS_%d_%d.csv", per, digits);
    if(mpv_write_csv(lsc, fn, true, digits, 0, n, false)) {
        printf("    Written the Log(S)-curve to %s\n", fn);
    } else {
        printf("****Could not write the Log(S)-curve to %s !!!\n", fn);
    }
    
    return true;
}

bool slc(int maxPer) {
    struct timeb te, ta;
    char time[80];
    
    ftime(&ta);

    int stp = LEVS_2POW + 2;
    printf("Computing level set of period %d ... ", stp);
    fflush(stdout);
    ftime(&te);

    def_mpfr(LEVS_PREC, r);
    mpfr_set_ld(r, 8, MPFR_RNDN);
    levc lc = levc_new(stp, LEVS_2POW, LEVC_MIN_PREC, r, LEVS_GUARD);
        
    if(lc == NULL) {
        printf("failed !");

        return false;
    }
    
    printf("refining ... ");
    fflush(stdout);
    levc mplc = levc_refine(lc, LEVS_2POW, LEVS_PREC);
    levc_free(lc);
    
    if(mplc == NULL) {
        printf("failed !\n");

        return false;
    }

    lapse(&te, time);
    printf("done in %s.\n", time);
    fflush(stdout);

    for (int per = stp + 1; per <= maxPer; per++) {
        scurve(mplc);
        
        printf("\nDown the rays toward period %d ... ", per);
        fflush(stdout);
        ftime(&te);

        lc = levc_next_period(mplc, false);
        
        if(lc == NULL) {
            printf("failed !\n");

            return false;
        }

        levc_free(mplc);
        mplc = lc;
        
        lapse(&te, time);
        printf("done in %s.\n", time);
        fflush(stdout);
    }
    
    scurve(mplc);
    
    lapse(&ta, time);
    printf("\nAll computed and saved in %s.\n", time);
    fflush(stdout);
    
    return true;
}

bool hyp_count(char *fn) {
    io_write(fn, "Period, Primitive, Non-primitive\n");
    for (int i = 1; i < 66; i++) {
        ulong pr = mandel_prim_hyp_count(i);
        ulong np = mandel_non_prim_hyp_count(i);
        
        io_append(fn, "%d, %lu, %lu\n", i, pr, np);
    }
    
    return true;
}

void ifs_dim(void) {
    ulong max = 1000000;
    ldbl s = 0, la = -2.48;
    
    for (ulong n = 2; n <= max; n++) {
        ldbl p = phi(n) * powl(n, la);
        s += p;
    }
    
    printf("S_%ld(%.3Lf) = %.10Lf\n", max, -la / 2, s);
}

bool nt_npn(fp80 nt, fp80 nt2, fp80 c, int n) {
    if(fabsl(c->x) <= 0.25 && fabsl(c->y)) {
        *nt = *c;
        
        return true;
    }
    
    fp80 v = {c->x, c->y}, b, d = {1, 0}, d2 = {0, 0};
    bool small = true;
    for (int k = 1; k < n && small; k++) {
        fp80_sqr(b, d);
        fp80_mul(d2, d2, v);
        fp80_add(d2, d2, b);
        fp80_muli(d2, d2, 2);
        
        fp80_mul(d, d, v);
        fp80_muli(d, d, 2);
        d->x += 1;
        
        fp80_sqr(v, v);
        fp80_add(v, v, c);
        
        ldbl mv2 = fp80_mod2(v), md2 = fp80_mod2(d);
        small = mv2 < 1e21 || mv2 * md2 < 1e41;
    }
    
    fp80 id;
    if(! fp80_inv(id, d) || ! fp80_div(b, d, v)) {
        return false;
    }
    
    fp80_mul(nt, v, id);
    fp80_mul(id, d2, id);
    fp80_sub(b, b, id);
    
    return fp80_inv(nt2, b);
}

bool root_pn(fp80 s, fp80 c, int n, int max_iter, ldbl eps2, double speed) {
    fp80 nt, nt2, v;
    if(! nt_npn(nt, nt2, c, n)) {
        return false;
    }
    
    fp80_muld(nt2, nt2, speed);
    
    ldbl m2 = fp80_mod2(nt), l2 = m2;
    fp80_sub(v, c, nt2);
    
    bool conv = false, div = false;
    for (int i = 1; i < max_iter && ! conv && ! div; i++) {
        if(! nt_npn(nt, nt2, v, n)) {
            return false;
        }
        
        fp80_muld(nt2, nt2, speed);
        fp80_sub(v, v, nt2);
        
        conv = fp80_mod2(nt2) < eps2;
        m2 = fp80_mod2(nt);
        conv = conv || m2 < eps2;
        
        div = m2 >= l2;
    }
    
    if(conv) {
        *s = *v;
        
        return true;
    }
    
    return false;
}

bool roots(int n) {
    ldbl eps = 1e-18, eps2 = eps * eps;
    pset ps;
    pset_init(ps, eps);
    
    ulong sc = 1L << n;
    for (ulong i = 0; i <= sc; i++) {
        fp80 st = {5 * cosl(i * PI / sc), 5 * sinl(i * PI / sc)}, c;
        if(root_pn(c, st, n, 300, eps2, 0.25)) {
            pset_add(ps, c);
        }
    }
    
    ulong rc = 2 * ps->count - ps->realCount;
    pset_clear(ps);
    
    if(rc == sc >> 1) {
        printf("Found all %lu roots of order %d.\n", rc, n);
    } else {
        printf("Found only %lu roots of order %d !\n", rc, n);
    }
    
    return true;
}

#define hsq2 (SQRT2 / 2)
static fp80_struct tg[8] = {{1, 0}, {hsq2, hsq2}, {0, 1}, {-hsq2, hsq2},
                           {-1, 0}, {-hsq2, -hsq2}, {0, -1}, {hsq2, -hsq2}};

typedef struct {
    int per;
    ldbl lev;
    uint circles;
    u128 st;
    u128 en;
    bool ok;
    ldbl len;
    ldbl area;
} struct_task_ll;

typedef struct_task_ll *tll;

void *tll_run(void *task, long ID, int index) {
    tll t = (tll) task;
    
    t->ok = true;
    fp80 en, pc, nc;
    u128_getl(pc, t->st);
    u128_getl(en, t->en);
    
    if(pc->x <= -2 + 1e-4) { // may exceed the precision and not gain much length, nor area
        *nc = *en;
        t->len = fp80_dist(pc, nc);
        t->area = (pc->x - nc->x) * (nc->y + pc->y);
        
        return task;
    }
    
    fp80_struct tgs[8];
    for (int i = 0; i < 8; i++) {
        fp80_mull(tgs + i, tg + i, t->lev);
    }
    
    t->len = 0;
    t->area = 0;
    for (int i = 0; i < t->circles && t->ok; i++) {
        for (int j = 1; j <= 8 && t->ok; j++) {
            if(! mandel_sol_refl(nc, NULL, pc, tgs + (j % 8), t->per, 20, 1e-16, 1)) {
                t->ok = false;
                
                continue;
            }
            
            t->len += fp80_dist(pc, nc);
            t->area += (pc->x - nc->x) * (nc->y + pc->y);
            
            *pc = *nc;
        }
    }
    
    t->ok = t->ok && fp80_dist(nc, en) < 1e-16;
    
    return t->ok ? task : NULL;
}

bool comp_len(levc lc, int threads, ldbl *pll, ldbl *plr) {
    struct timeb ts;
    ftime(&ts);
    
    int per = lc->period;
    uint seg = (uint) levc_segs(lc);
    tll tasks[seg];
    
    u128 a;
    def_mpc(128, c);
    levc_point(c, lc, 0);
    u128_set(a, c);
    ldbl lev = mpfr_get_ld(lc->radius, MPFR_RNDN);
    
    uint circ = 1 << (per - lc->ang2pow - 1);
    for (int i = 0; i < seg; i++) {
        tasks[i] = malloc(sizeof(struct_task_ll));
        
        *tasks[i]->st = *a;
        levc_point(c, lc, i + 1);
        u128_set(a, c);
        *tasks[i]->en = *a;
        
        tasks[i]->per = lc->period;
        tasks[i]->lev = lev;
        tasks[i]->circles = circ;
        tasks[i]->ok = false;
        tasks[i]->len = 0;
        tasks[i]->area = 0;
    }
    
    bool ok = true;
    if(threads == 1) {
        for (int i = 0; i < seg && ok; i++) {
            ok = ok && tll_run(tasks[i], i, 0) != NULL;
        }
    } else {
        ok = mth_batch((void **)tasks, seg);
    }
    
    ldbl len = 0, area = 0;
    for (int i = 0; i < seg; i++) {
        ok = ok && tasks[i]->ok;
        if(ok) {
            len += tasks[i]->len;
            area += tasks[i]->area;
        }
        
        free(tasks[i]);
    }
    
    if(! ok) {
        printf("Could not compute the curve of period %d !\n", per);
        
        return false;
    }
    
    ldbl rm1 = expm1l(logl(lev) / (1L << (per - 1)));
    ldbl lrm1 = -logl(rm1);
    ldbl ll = logl(len / PI);
    
    char time[100];
    lapse(&ts, time);
    
    if(*pll < 1) {
        ldbl b1 = ll / lrm1;
        printf("Period %2d, len = %9.6Lf, area = %.6Lf, B(1) = %.6Lf, r - 1 = %Lg, computed in %s\n",
               per, len, area, b1, rm1, time);
    } else {
        ldbl b1 = (ll - *pll) / (lrm1 - *plr);
        printf("Period %2d, log(len ratio) = %9.6Lf, area = %.6Lf, B(1) = %.6Lf, log(r ratio) = %Lg, computed in %s\n",
               per, ll - *pll, area, b1, lrm1 - *plr, time);
    }
    fflush(stdout);
    
    *pll = ll;
    *plr = lrm1;
    
    return true;
}

void levCurves(ldbl lev, int max_per, int threads) {
    printf("Computing level curves of level %Lg, up to period %d, with %d threads:\n",
           lev, max_per, threads);
    fflush(stdout);
    
    def_mpfr(128, r);
    mpfr_set_ld(r, lev, MPFR_RNDN);
    
    if(threads > 1) {
        mth_init(threads, &tll_run, 30, 30);
    }
    
    int n = 14;
    printf("Computing the level curve, please wait ... ");
    fflush(stdout);
    levc lc = levc_new(n, n - 1, 128, r, LEVC_MIN_GUARD);
    printf("done\n\n");
    fflush(stdout);
    
    ldbl pll = 0, plr = 0;
    comp_len(lc, threads, &pll, &plr);
    
    for (int i = n + 1; i <= max_per; i++) {
        levc l = levc_next_period(lc, false);
        levc_free(lc);
        lc = l;
        
        comp_len(lc, threads, &pll, &plr);
    }
    
    if(threads > 1) {
        mth_end_threads(true, true);
    }
}

void proveMis(int pp, int per, ldbl r, ldbl cr, char *x, char *y, bool refine) {
    defs_mpc(256, c);
    mpc_set_str(c, x, y);
    mpfr_printf("(%.37Rf, %.37Rf) ->\n", c->x, c->y);
    if(refine) {
        mandel_mis_root_ref(c, c, NULL, pp, per, 50, 3);
        mpfr_printf("(%.37Rf, %.37Rf)\n\n", c->x, c->y);
    }
    
    if(mandel_is_mis(c, pp, per, r)) {
        mpfr_printf("(%.37Rf, %.37Rf) is in Mis(%d, %d)\n", c->x, c->y, pp, per);
    } else {
        mpfr_printf("(%.37Rf, %.37Rf) is NOT in Mis(%d, %d) !\n", c->x, c->y, pp, per);
    }
    
    if(mandel_conv_npp(c, pp, per, r, cr, false)) {
        mpfr_printf("(%.37Rf, %.37Rf, %Lg) is in the basin of the center\n", c->x, c->y, cr);
    } else {
        mpfr_printf("(%.37Rf, %.r37Rf, %Lg) is in NOT the basin of the center !\n", c->x, c->y, cr);
    }
}

void testMinDist(char *fn) {
    char time[100];
    struct timeb ts;
    ftime(&ts);
    nset ps = nset_load(fn, false);

    lapse(&ts, time);
    printf("Loaded in %s.\n", time);

    ldbl omd = nset_min_dist(ps);

    lapse(&ts, time);
    printf("Min dist in %s.\n", time);
}

void refine(int pp, int per, int ind) {
    nset ns = mis_load_results(pp, per, ind);
    if(ns == NULL) {
        printf("Could not load mis(%d, %d, %d)\n", pp, per, ind);
        
        return;
    } else {
        printf("Loaded mis(%d, %d, %d)\n", pp, per, ind);
    }
    
    u128 p;
    nset rs = mis_new_set(false);
    def_mpc(180, c);
    
    bool ok = true;
    for (int i = 0; i <  ns->count && ok; i++) {
        ok = ok && nset_point(p, ns, i);
        ok = ok && u128_get(c, p);
        ok = ok && mandel_mis_root_ref(c, c, NULL, pp, per, 50, 2);
        ok = ok && u128_set(p, c);
        ok = ok && nset_add(rs, p);
    }
    
    nset_free(ns);
    char fn[150];
    mis_file_name(fn, 150, pp, per, ind);
    ok = ok && nset_lock(rs);
    ok = ok && nset_write(rs, fn);
    
    if(ok) {
        printf("Refined all %lu points from %s\n", rs->count, fn);
    } else {
        printf("Some error occurred !\n");
    }
    
    nset_free(rs);
}

void doublons(int pp, int per, ldbl eps, ldbl del, ldbl dst) {
    pset ps, ns;
    pset_init(ps, eps);
    pset_init(ns, del);
    
    char fn[100];
    mis_file_name(fn, 100, pp, per, 0);
    nset ms = mis_load_results(pp, per, 0);
    if(ms == NULL) {
        printf("Cannot read %s\n", fn);
        
        return;
    }
    printf("Read %lu points from %s (total %lu, real %lu)\n",
           ms->count, fn, 2 * ms->count - ms->realCount, ms->realCount);
    
    snprintf(fn, 100, "quick/mis%02d-%02d_19.csv", pp, per);
    if(! pset_read(ps, fn)) {
        printf("Could not read %s !\n", fn);
        
        return;
    }
    pset_lock(ps);
    printf("Read %lu points from %s (total %lu, real %lu)\n\n",
           ps->count, fn, 2 * ps->count - ps->realCount, ps->realCount);
    
    defs_mpc(128, u, v);
    fp80 c, p;
    u128 s, t;
    for (int i = 0; i < ps->count; i++) {
        pset_point(p, ps, i);
        
        // check for distance
        if(! pset_add(ns, p)) {
            long ind = pset_index(ns, p);
            if(ind < 0 || ! pset_point(c, ns, ind)) {
                printf("Cannot add and cannot find (%.19Lf, %.19Lg)\n", p->x, p->y);
                
                continue;
            }
            
            mpc_set80(u, c);
            mpc_set80(v, p);
            if(! mandel_mis_root_ref(u, u, NULL, pp, per, 30, 2)) {
                mpfr_printf("Cannot refine (%.22Rf, %.22Rf)\n", u->x, u->y);
                
                continue;
            }
            
            if(! mandel_mis_root_ref(v, v, NULL, pp, per, 30, 2)) {
                mpfr_printf("Cannot refine (%.22Rf, %.22Rf)\n", u->x, u->y);
                
                continue;
            }
            
            if(mpc_distl(u, v) < 1e-22) {
                mpfr_printf("(%.24Rf, %.24Rg) has two copies: (%.19Lf, %.19Lg) and (%.19Lf, %.19Lg)\n",
                            u->x, u->y, c->x, c->y, p->x, p->y);
            } else {
                mpfr_printf("(%.19Lf, %.19Lg) -> (%.24Rf, %.24Rf)\n(%.19Lf, %.19Lg) -> (%.24Rf, %.24Rf)\n\n",
                            c->x, c->y, u->x, u->y, p->x, p->y, v->x, v->y);
            }
        }
        
        // check if not extra
        u128_setl(s, p);
        long pos = nset_search(ms, s);
        if(pos < 0) {
            if(nset_point(t, ms, -pos - 1)) {
                u128_getl(c, t);
                if(fp80_dist(c, p) < dst) {
                    continue;
                }
            }
            
            if(nset_point(t, ms, -pos - 2)) {
                u128_getl(c, t);
                if(fp80_dist(c, p) < dst) {
                    continue;
                }
            }
                   
            printf("The point (%.19Lf, %.19Lg) is not in Mis(%d, %d)\n",
                   p->x, p->y, pp, per);
            mpc_set80(u, p);
            if(! mandel_mis_root_ref(u, u, NULL, pp, per, 30, 2)) {
                mpfr_printf("Cannot refine (%.22Rf, %.22Rf)\n", u->x, u->y);
                
                continue;
            }
            
            mpfr_printf("       -> (%.26Rf, %.26Rf)\n", u->x, u->y);
        }
    }
    
    printf("\n");
    for (int i = 0; i < ms->count; i++) {
        nset_point(s, ms, i);
        u128_getl(c, s);
        
        if(! pset_contains(ps, c)) {
            u128_get(u, s);
            mpfr_printf("The point (%.26Rf, %.26Rf) is not in %s\n", u->x, u->y, fn);
        }
    }
    
    printf("\nThe new set has %lu points (total %lu, real %lu)\n",
           ns->count, 2 * ns->count - ns->realCount, ns->realCount);
    
    pset_clear(ps);
    pset_clear(ns);
}

void oldMain(void) {
    
//    searchLevm(33, 1, x[0], y[0]);
//    searchLevm(33, 1, x[1], y[1]);
//
//    levmDist(33, 1, 15358, "-1.9415886521090466205707099100003", "0.0076777655769650745352359000000");
//    levmDist(33, 1, 15359, "-1.9415886521090466205707099100003", "0.0076777655769650745352359000000");
//    levmDist(33, 1, 15360, "-1.9415886521090466205707099100003", "0.0076777655769650745352359000000");
//    levmDist(33, 1, 15361, "-1.9415886521090466205707099100003", "0.0076777655769650745352359000000");
//    levmDist(33, 1, 15362, "-1.9415886521090466205707099100003", "0.0076777655769650745352359000000");
    
//    testMisRaw(33, 1, 15359, x[0], y[0], x[1], y[1]);
//    testMisRaw(33, 1, 15360, x[0], y[0], x[1], y[1]);
    
//    mis_ref(33, 1, 170, "-1.9415886521090466205707099100003", "0.0076777655769650745352359000000");
//    mis_ref(33, 1, 170, "-1.9415886521090466032836203199993", "0.0076777655769650775743220800000");
    
//    mis_root(33, 1, 120, x[2], y[2]);
//    mis_lc(33, 1, 120, x[2], y[2]);
    
//    testMisSimpleSet(28, 1);
    
//    testHypLc();
//    countMultCrit();
//    compInvPowers(argc, argv);
//    computeMissingPowers(argc, argv);
    //    computeallHypMissing(22, 30, 1e-60);
//        computeallMisMissing(22, 30, 1e-60);
//
//    mpv r = solveInvRoots("mis/33-01/mis33-01_missingPows30o.mpv");
//    refineAndSave(33, 1, r, "mis/33-01/mis33-01_missing.nset");
//    mpv_free(r);
//
//    nsetSearchNset(33, 1, "mis/33-01/mis33-01_missing.nset");
    
//    mpvDiff("mis/33-01/mis33-01_004InvPows30.mpv", "mis/33-01/mis33-01_004InvPows30o.mpv",
//            "mis/33-01/mis33-01_004InvPows30diff.mpv");
//    testPolyInvPows();
//    nsetDiff("mis/33-01/mis33-01_004n.nset", "mis/33-01/mis33-01_004.nset");
    
//    tmapTest(4, 1);
    
//    printf("We have %d arguments\n", ARG_COUNT(1,5,67,4));
    
//    nsetSearchMpvInNSet(argv[0], argv[1]);
    
//    int pp, per, raw;
//    sscanf(argv[0], "%d", &pp);
//    sscanf(argv[1], "%d", &per);
//    sscanf(argv[2], "%d", &raw);
//
//    resetEpsilon(pp, per, raw == 1 ? true : false);
    
//    testImg(3, 3, 1024, 1024, 0, 0, 8, true);
//    compImg(6, 50, 180, 8);
    
//    slc(44);
//    testRay(240, 1, 6);
    
//    hyp_count("hyp_count.csv");
    
//    ifs_dim();
    
//    roots(3);
    
//    if(argc < 3) {
//        printf("Two arguments: max_per threads level\n");
//    } else {
//        ldbl lev;
//        int max_per, threads;
//        sscanf(argv[0], "%d", &max_per);
//        sscanf(argv[1], "%d", &threads);
//        sscanf(argv[2], "%Lf", &lev);
//        levCurves(lev, max_per, threads);
//    }
    
//    testMinDist("hyp30_7.nset");
//    proveMis(31, 1, 1e-35, 1e-32, "-1.771257023356245744857734568401999999999",
//             "0.06616150907984357256001719991799999999661", true);
    
}

void testMiss(void) {
    fp80 p8 = {-16, 0}, c8, nt8, snt8;
    defs_mpc(128, c, p, nt, snt, v2, v3, d2, d3);
    mpc_set80(p, p8);
    
//    mandel_mis_ntl(nt8, p8, NULL, 2, 1);
    mandel_miss_ntl(snt8, p8, NULL, 2, 1);
//    
//    mandel_mis_nt(nt, p, 2, 1);
    
//    mpfr_printf("(%.20Lf, %.20Lf)\n(%.20Rf, %.20Rf)\n(%.20Lf, %.20Lf)\n(%.20Rf, %.20Rf)",
//                nt8->x, nt8->y, nt->x, nt->y, snt8->x, snt8->y, snt->x, snt->y);
    
    mandel_val_der(v2, d2, p, 1);
    mandel_val_der(v3, d3, p, 2);
    mpc_add(v2, v2, v3);
    mpc_add(d2, d2, d3);
    mpc_div(nt, v2, d2);
    
    mandel_miss_nt_sol(snt, p, NULL, 2, 1);
    
    mpfr_printf("(%.20Rf, %.20Rf)\n(%.20Lf, %.20Lf)\n(%.20Rf, %.20Rf)",
                   nt->x, nt->y, snt8->x, snt8->y, snt->x, snt->y);
}

int playNicuMain(int argc, const char * argv[]) {
    char time[100];
    time_stamp(time, 100, true, true);
    printf("Playground started at %s :\n\n", time);
    
    struct timeb ts;
    ftime(&ts);
    
//    refine(31, 1, 1);
    
//    doublons(3, 23, 6e-19, 1e-18, 7e-19);
    testMiss();
    
    lapse(&ts, time);
    printf("\nTotal time %s.\n", time);
    
    return 0;
}
