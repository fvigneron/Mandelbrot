//
//  Mandelbrot.c
//
//  A collection of methods related to the Mandelbrot set M.
//  Two precisions are used: long double (ldbl) and arbitrary (mpfr, mpc).
//
//  TODO: recheck all the code and maths behind this
// - that would be nice, yeah !
//
//  TODO: lower bound of the distance to M with proof
// - actually, we should better remove this for the first two papers
//
//  TODO: auto search for the error / radius that provides a proof for the existence of a root
// - not sure about that one: it would be slower, and we want uniform epsilon as we do not store it
//
//  The following two references explain the maths behind some of the
//  functions below.
//
//  [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
//  [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2021.
//
//  Created by MIHALACHE Nicolae on 12/19/19.
//  Copyright © 2019 MIHALACHE Nicolae. All rights reserved.
//
// TODO: look up which GNU licence to use; we want people to be able to use the
// software for non-comercial purposes, to be able to modify it etc., but with
// one CONDITION: to cite the authors, both the code repository and the papers
// [1] and [2] on which this is based

#include <stdlib.h>
#include <string.h>

#include "mandel.h"
#include "arithmetic.h"

// MARK: local variables

// needed for stats
static ulong lastIter = 0;
static ulong totNewton = 0;

// extra details from mandel_isMisl()
static int ppFail, perFail;

// MARK: stats

ulong mandel_last_iter(void) {
    return lastIter;
}

ulong mandel_tot_nt(void){
    return totNewton;
}

// MARK: quick versions with low precision (long double, fp80)

ldbl mandel_distl(fp80 c, int maxIter) {
    lastIter = 0;
    
    ldbl cx = c->x, cy = c->y;
    ldbl vx = cx, vy = cy, dx = 1, dy = 0;
    ldbl x2 = vx * vx, y2 = vy * vy, b1 = 0;
    
    int small = 1;
    
    int k = 1;
    for (; k < maxIter + 7 && small; k++) {
        b1 = dx;
        dx = 2 * (dx * vx - dy * vy) + 1;
        dy = 2 * (b1 * vy + dy * vx);
        
        b1 = vx;
        vx = x2 - y2 + cx;
        vy = 2 * b1 * vy + cy;
        
        x2 = vx * vx;
        y2 = vy * vy;
        
        small = x2 + y2 <= 1E44;
        lastIter = x2 + y2 <= 25 ? k + 1 : lastIter; // for Carleson's conjecture, escape from D(0, 5)
    }
        
    ldbl m2dp = dx * dx + dy * dy;
    
    if(small || m2dp < 1) { // bounded
        return -1;
    }
    
    ldbl m2pk = x2 + y2;
    ldbl d = sqrtl(m2pk / m2dp);
    
    ldbl lm2p = logl(m2pk * 0.25);
    
    if(k < 10000) { // only so that the exponent of ldbl does not overflow
        d *= ldexpl(-expm1l(-ldexpl(lm2p, -k)), k - 2);
    } else {
        d *= 0.25 * lm2p;
    }
    
    return d;
}

void mandel_vall(fp80 v, fp80 c, int per) {
    ldbl cx = c->x, cy = c->y;
    ldbl vx = cx, vy = cy;
    ldbl x2, y2, b1;
    
    for (int k = 1; k < per; k++) {
        x2 = vx * vx;
        y2 = vy * vy;
        
        b1 = vx;
        vx = x2 - y2 + cx;
        vy = 2 * b1 * vy + cy;
    }
    
    v->x = vx;
    v->y = vy;
}

void mandel_mis_vall(fp80 v, fp80 c, int pp, int per) {
    ldbl cx = c->x, cy = c->y, vx = cx, vy = cy, vpx = 0, vpy = 0, bx;
    
    int n = pp + per;
    for (int k = 1; k < n; k ++) {
        if(k == pp) {
            vpx = vx;
            vpy = vy;
        }
        
        bx = vx * vx - vy * vy;
        vy *= 2 * vx;
        vx = bx + cx;
        vy += cy;
    }
    
    v->x = vx - vpx;
    v->y = vy - vpy;
}

void mandel_miss_vall(fp80 v, fp80 c, int pp, int per) {
    ldbl cx = c->x, cy = c->y, vx = cx, vy = cy, vpx = 0, vpy = 0, bx;
    
    int n = pp + per - 1;
    for (int k = 1; k < n; k ++) {
        if(k == pp - 1) {
            vpx = vx;
            vpy = vy;
        }
        
        bx = vx * vx - vy * vy;
        vy *= 2 * vx;
        vx = bx + cx;
        vy += cy;
    }
    
    v->x = vx + vpx;
    v->y = vy + vpy;
}

bool mandel_val_derl(fp80 v, fp80 d, fp80 c, int per) {
    if(v == NULL || d == NULL || c == NULL || per < 1) {
        return false;
    }

    ldbl cx = c->x, cy = c->y;
    ldbl vx = cx, vy = cy, dx = 1, dy = 0, b1;
    
    for (int k = 1; k < per; k++) {
        b1 = dx;
        dx = 2 * (dx * vx - dy * vy) + 1;
        dy = 2 * (b1 * vy + dy * vx);
                
        b1 = vx;
        vx = vx * vx - vy * vy + cx;
        vy = 2 * b1 * vy + cy;
    }
    
    v->x = vx;
    v->y = vy;
    
    d->x = dx;
    d->y = dy;
    
    return true;
}

fp80_ptr mandel_val_dersl(fp80 c, int per, int ders) {
    if(c == NULL || per < 1 || ders < 2 || ders > 67) {
        return NULL;
    }
    
    int h = ders >> 1;
    ulong ms[ders - 1][h + 1];
    
    // M(2)
    ulong *m = ms[0], *pm;
    for (int i = 0; i <= h; i++) {
        m[i] = i < 2 ? 1 : 0;
    }
    
    // M(ders)
    for (int n = 3; n <= ders; n++) {
        pm = m;
        m = ms[n - 2];
        
        ulong ov = 1, v;
        m[0] = 1;
        
        int k = n >> 1;
        for (int i = 1; i < k; i++) {
            v = pm[i];
            m[i] = ov + v;
            ov = v;
        }
        
        m[k] = ov + (n & 1 ? pm[k] << 1 : 0);
        
        for (int i = k + 1; i <= h; i++) {
            m[i] = 0;
        }
    }

    fp80_struct v1[ders + 1], v2[ders + 1];
    fp80_ptr p = v1, a = v2, b;
    a[0] = *c;
    a[1].x = 1;
    a[1].y = 0;
    for (int i = 2; i <= ders; i++) {
        a[i].x = 0;
        a[i].y = 0;
    }
    
    for (int k = 1; k < per; k++) {
        b = p;
        p = a;
        a = b;
        
        // p_k(c)
        fp80_sqr(a, p);
        fp80_add(a, a, c);
        
        // p_k'(c)
        fp80_mul(a + 1, p, p + 1);
        a[1].x = 2 * a[1].x + 1;
        a[1].y *= 2;
        
        // all other derivatives
        for (int n = 2; n <= ders; n++) {
            m = ms[n - 2];
            fp80 s = {0, 0}, d;
            ulong cf;
            int hn = n >> 1;
            
            for (int i = 0; i <= hn; i++) {
                cf = m[i];
                
                fp80_mul(d, p + i, p + (n - i));
                d->x *= cf;
                d->y *= cf;
                
                fp80_add(s, s, d);
            }
            
            s->x *= 2;
            s->y *= 2;
            a[n] = s[0];
        }
    }
    
    fp80_ptr v = malloc(sizeof(fp80_struct) * (ders + 1));
    for (int i = 0; i <= ders; i++) {
        v[i] = a[i];
    }
    
    return v;
}

void mandel_mis_val_derl(fp80 v, fp80 d, fp80 c, int pp, int per) {
    ldbl cx = c->x, cy = c->y, vx = cx, vy = cy, vpx = 0, vpy = 0, bx;
    ldbl dx = 1, dy = 0, dpx = 0, dpy = 0;
    
    int n = pp + per;
    for (int k = 1; k < n; k ++) {
        if(k == pp) {
            vpx = vx;
            vpy = vy;
            
            dpx = dx;
            dpy = dy;
        }
        
        bx = dx;
        dx = 2 * (dx * vx - dy * vy) + 1;
        dy = 2 * (bx * vy + dy * vx);
        
        bx = vx * vx - vy * vy;
        vy *= 2 * vx;
        vx = bx + cx;
        vy += cy;
    }
    
    v->x = vx - vpx;
    v->y = vy - vpy;
    
    d->x = dx - dpx;
    d->y = dy - dpy;
}

bool mandel_ntl(fp80 v, fp80 c, int per) {
    ldbl cx = c->x, cy = c->y;
    ldbl vx = cx, vy = cy, dx = 1, dy = 0;
    ldbl x2 = vx * vx, y2 = vy * vy, b1 = 0;
    
    int small = 1;
    
    int k = 1;
    for (; k < per && small; k++) {
        b1 = dx;
        dx = 2 * (dx * vx - dy * vy) + 1;
        dy = 2 * (b1 * vy + dy * vx);
        
        b1 = vx;
        vx = x2 - y2 + cx;
        vy = 2 * b1 * vy + cy;
        
        x2 = vx * vx;
        y2 = vy * vy;
        
        small = x2 + y2 <= T66 || dx * dx + dy * dy <= T60;
    }
    
    totNewton ++;
    
    x2 = dx * dx + dy * dy;
    
    if(x2 == 0) {
        v->x = INFINITY;
        v->y = INFINITY;
        
        return false;
    }
    
    if(per > k) {
        x2 = ldexpl(x2, per - k);
    }
    
    b1 = vx;
    vx = vx * dx + vy * dy;
    vy = dx * vy - dy * b1;
    
    v->x = vx / x2;
    v->y = vy / x2;
    
    return true;
}

ldbl mandel_mis_ntl(fp80 v, fp80 c, fp80 t, int pp, int per) {
    ldbl cx = c->x, cy = c->y;
    ldbl vx = cx, vy = cy, dx = 1, dy = 0, b1, x2 = vx * vx, y2 = vy * vy;
       
    int k = 1;
    for (; k < pp && x2 + y2 <= T66; k++) {
        b1 = dx;
        dx = 2 * (dx * vx - dy * vy) + 1;
        dy = 2 * (b1 * vy + dy * vx);
        
        b1 = vx;
        vx = x2 - y2 + cx;
        vy = 2 * b1 * vy + cy;
        
        x2 = vx * vx;
        y2 = vy * vy;
    }
    
    int n = pp + per;
    if(x2 + y2 > T66) {
        x2 = ldexpl(dx * dx + dy * dy, n - k);
        
        b1 = vx;
        vx = vx * dx + vy * dy;
        vy = dx * vy - dy * b1;
        
        v->x = vx / x2;
        v->y = vy / x2;
        
        totNewton ++;
        
        return x2;
    }
    
    ldbl dpx = dx, dpy = dy, vpx = vx, vpy = vy;
    for (; k < n && x2 + y2 <= T66; k++) {
        b1 = dx;
        dx = 2 * (dx * vx - dy * vy) + 1;
        dy = 2 * (b1 * vy + dy * vx);
        
        b1 = vx;
        vx = x2 - y2 + cx;
        vy = 2 * b1 * vy + cy;
        
        x2 = vx * vx;
        y2 = vy * vy;
    }
        
    totNewton ++;
    
    dx -= dpx;
    dy -= dpy;
    
    ldbl d2 = dx * dx + dy * dy;
    
    if(d2 == 0) {
        v->x = INFINITY;
        v->y = INFINITY;
        
        return 0;
    }
    
    if(k < n) {
        x2 = ldexpl(d2, n - k);
        
        b1 = vx;
        vx = vx * dx + vy * dy;
        vy = dx * vy - dy * b1;
        
        v->x = vx / x2;
        v->y = vy / x2;
        
        totNewton ++;
        
        return x2;
    }
    
    if(t != NULL) {
        vx -= vpx + t->x;
        vy -= vpy + t->y;
    } else {
        vx -= vpx;
        vy -= vpy;
    }
    
    b1 = vx;
    vx = vx * dx + vy * dy;
    vy = dx * vy - dy * b1;
    
    v->x = vx / d2;
    v->y = vy / d2;
    
    return d2;
}

ldbl mandel_miss_ntl(fp80 nt, fp80 c, fp80 t, int pp, int per) {
    ldbl cx = c->x, cy = c->y;
    ldbl vx = cx, vy = cy, dx = 1, dy = 0, b1, x2 = vx * vx, y2 = vy * vy;
       
    for (int k = 1; k < pp - 1 && x2 + y2 < OUT2; k++) {
        b1 = dx;
        dx = 2 * (dx * vx - dy * vy) + 1;
        dy = 2 * (b1 * vy + dy * vx);
        
        b1 = vx;
        vx = x2 - y2 + cx;
        vy = 2 * b1 * vy + cy;
        
        x2 = vx * vx;
        y2 = vy * vy;
    }
    
    ldbl dpx = dx, dpy = dy, vpx = vx, vpy = vy;
    for (int i = 0; i < per && x2 + y2 < OUT2; i++) {
        b1 = dx;
        dx = 2 * (dx * vx - dy * vy) + 1;
        dy = 2 * (b1 * vy + dy * vx);
        
        b1 = vx;
        vx = x2 - y2 + cx;
        vy = 2 * b1 * vy + cy;
        
        x2 = vx * vx;
        y2 = vy * vy;
    }
    
    if(x2 + y2 >= OUT2) {
        return 0;
    }
    
    totNewton ++;
    
    dx += dpx;
    dy += dpy;
    
    ldbl d2 = dx * dx + dy * dy;
    
    if(t != NULL) {
        vx += vpx - t->x;
        vy += vpy - t->y;
    } else {
        vx += vpx;
        vy += vpy;
    }
    
    b1 = vx;
    vx = vx * dx + vy * dy;
    vy = dx * vy - dy * b1;
    
    nt->x = vx / d2;
    nt->y = vy / d2;
    
    return d2;
}

bool mandel_rootl(fp80 v, fp80 c, int per, int maxIter, ldbl error) {
    return mandel_root_refl(v, c, per, maxIter, error, 3, 3);
}

bool mandel_root_refl(fp80 v, fp80 c, int per, int maxIter, ldbl error, int refine, double maxDist) {
    fp80 lc = {c->x, c->y};
    fp80 nc;
    ldbl e2 = error * error, m2 = maxDist * maxDist, dx, dy;
    ldbl sx = c->x, sy = c->y;
    
    int div = 0;
    int conv = 0;
    
    ulong tn = totNewton;
    for (int i = 0; i < maxIter && ! div && ! conv; i++) {
        if(! mandel_ntl(nc, lc, per)) {
            return false;
        }
        
        lc->x -= nc->x;
        lc->y -= nc->y;
        
        dx = lc->x - sx;
        dy = lc->y - sy;
        
        div = dx * dx + dy * dy > m2 || lc->x * lc->x + lc->y * lc->y >= OUT2;
        conv = nc->x * nc->x + nc->y * nc->y < e2; 
    }
    
    if(div || ! conv) {
        lastIter = totNewton - tn;
        
        return false;
    }
    
    for (int i = 0; i < refine; i++) {
        if(! mandel_ntl(nc, lc, per)) {
            return false;
        }
        
        lc->x -= nc->x;
        lc->y -= nc->y;
    }
    
    v->x = lc->x;
    v->y = lc->y;

    lastIter = totNewton - tn;
    
    return true;
}

ldbl mandel_nt_soll(fp80 v, fp80 c, fp80 t, int per) {
    fp80 lc = {c->x, c->y};
    fp80 vc = {lc->x, lc->y};
    fp80 lt = {t->x, t->y};
    ldbl dx = 1, dy = 0;
    ldbl x2 = vc->x * vc->x, y2 = vc->y * vc->y, b1 = 0, b2 = 0;
    
    ldbl large = lt->x * lt->x + lt->y * lt->y;
    large = large < 1 ? 1 : large;
    large *= T70;
    
    int small = 1;
    
    int k = 1;
    for (; k < per && small; k++) {
        b1 = vc->x;
        b2 = dx;
        
        dx = 2 * (dx * b1 - dy * vc->y) + 1;
        dy = 2 * (b2 * vc->y + dy * b1);
        
        vc->x = x2 - y2 + lc->x;
        vc->y = 2 * b1 * vc->y + lc->y;
        
        x2 = vc->x * vc->x;
        y2 = vc->y * vc->y;
        
        small = x2 + y2 <= large || dx * dx + dy * dy <= T60;
    }
    
    if(per == k) {
        x2 = dx * dx + dy * dy;
        vc->x -= lt->x;
        vc->y -= lt->y;
    } else {
        x2 = ldexpl(dx * dx + dy * dy, per - k);
    }
    
    b1 = vc->x;
    vc->x = vc->x * dx + vc->y * dy;
    vc->y = dx * vc->y - dy * b1;
    
    v->x = vc->x / x2;
    v->y = vc->y / x2;
    
    totNewton ++;
    
    return x2;
}

bool mandel_solutionl(fp80 v, ldbl *err, fp80 c, fp80 t, int per, int maxIter, ldbl error) {
    return mandel_sol_refl(v, err, c, t, per, maxIter, error, 3);
}

bool mandel_sol_refl(fp80 v, ldbl *err, fp80 c, fp80 t, int per, int maxIter, ldbl error, int refine) {
    fp80 lc = {c->x, c->y};
    fp80 nc;
    ldbl e2 = error * error;
    
    int div = 0;
    int conv = 0;
    
    ulong tn = totNewton;
    for (int i = 0; i < maxIter && ! div && ! conv; i++) {
        mandel_nt_soll(nc, lc, t, per);
        
        lc->x -= nc->x;
        lc->y -= nc->y;
        
        div = lc->x * lc->x + lc->y * lc->y >= OUT2;
        conv = nc->x * nc->x + nc->y * nc->y < e2;
    }
    
    if(div || ! conv) {
        v->x = lc->x;
        v->y = lc->y;
        
        lastIter = totNewton - tn;
        
        return false;
    }
    
    for (int i = 0; i < refine; i++) {
        mandel_nt_soll(nc, lc, t, per);
        
        lc->x -= nc->x;
        lc->y -= nc->y;
    }
    
    v->x = lc->x;
    v->y = lc->y;

    if(err != NULL) {
        *err = fp80_mod(nc);
    }
    
    lastIter = totNewton - tn;
    
    return true;
}

bool mandel_mis_rootl(fp80 root, fp80 c, int pp, int per, int maxIter, ldbl error) {
    fp80 z = {0, 0};
    
    return mandel_mis_root_refl(root, c, z, pp, per, maxIter, error, 1, 3);
}

bool mandel_miss_rootl(fp80 root, fp80 c, int pp, int per, int maxIter, ldbl error) {
    fp80 z = {0, 0};
    
    return mandel_miss_root_refl(root, c, z, pp, per, maxIter, error, 2, 100);
}

bool mandel_mis_root_refl(fp80 root, fp80 c, fp80 t, int pp, int per, int maxIter, ldbl error, int refine, double maxDist) {
    fp80 lc = {c->x, c->y};
    fp80 nc;
    ldbl e2 = error * error, m2 = maxDist * maxDist, dx, dy;
    ldbl sx = c->x, sy = c->y;
    
    int div = 0;
    int conv = 0;
    
    ulong tn = totNewton;
    for (int i = 0; i < maxIter && ! conv; i++) {
        if(mandel_mis_ntl(nc, lc, t, pp, per) < 1E-25) {
            return false;
        }
        
        lc->x -= nc->x;
        lc->y -= nc->y;
        conv = nc->x * nc->x + nc->y * nc->y < e2;
    }

    // divergence is very rare, it is faster to check it here
    dx = lc->x - sx;
    dy = lc->y - sy;
    
    div = dx * dx + dy * dy > m2 || lc->x * lc->x + lc->y * lc->y >= OUT2;
    
    if(div || ! conv) {
        lastIter = totNewton - tn;
        
        return false;
    }
    
    // check if not close to hyperbolic and do one refinement
    fp80 vp, dp;
    mandel_val_derl(vp, dp, lc, pp);
    
    ldbl cx = lc->x, cy = lc->y;
    ldbl vx = vp->x, vy = vp->y, b1, dzx = 1, dzy = 0, x2, y2;
    dx = dp->x;
    dy = dp->y;
    
    for (int k = 0; k < per; k++) {
        b1 = dx;
        dx = 2 * (dx * vx - dy * vy) + 1;
        dy = 2 * (b1 * vy + dy * vx);
        
        b1 = dzx;
        dzx = 2 * (dzx * vx - dzy * vy);
        dzy = 2 * (b1 * vy + dzy * vx);
                
        b1 = vx;
        x2 = vx * vx;
        y2 = vy * vy;
        vx = x2 - y2 + cx;
        vy = 2 * b1 * vy + cy;
        
        m2 = x2 + y2;
        if(x2 + y2 < error || m2 > OUT2) { // hyperbolic or divergent
            return false;
        }
    }
    
    ldbl dz2 = dzx * dzx + dzy * dzy;
    if(dz2 < 1) { // in hyperbolic component
        return false;
    }
    
    dx -= dp->x;
    dy -= dp->y;
    
    vx -= vp->x + t->x;
    vy -= vp->y + t->y;
    
    ldbl d2 = dx * dx + dy * dy;
    if(d2 == 0 || vx * vx + vy * vy > d2 * e2) { // cannot compute or not pre-periodic
        return false;
    }
    
    b1 = vx;
    vx = vx * dx + vy * dy;
    vy = dx * vy - dy * b1;
    
    lc->x -= vx / d2;
    lc->y -= vy / d2;
    
    for (int i = 0; i < refine - 1; i++) {
        if(! (mandel_mis_ntl(nc, lc, t, pp, per) >= 1E-25)) {
            return false;
        }
        
        lc->x -= nc->x;
        lc->y -= nc->y;
    }

    lastIter = totNewton - tn;
    
    *root = *lc;
    
    return true;
}

bool mandel_miss_root_refl(fp80 v, fp80 c, fp80 t, int pp, int per, int maxIter, ldbl error,
                          int refine, double maxDist) {
    fp80 lc = {c->x, c->y};
    fp80 nc;
    ldbl e2 = error * error, m2 = maxDist * maxDist, dx, dy;
    ldbl sx = c->x, sy = c->y;
    
    int div = 0;
    int conv = 0;
    
    ulong tn = totNewton;
    for (int i = 0; i < maxIter && ! div && ! conv; i++) {
        ldbl d2 = mandel_miss_ntl(nc, lc, t, pp, per);
        if(d2 < 1E-25) {
            return false;
        }
        
        lc->x -= nc->x;
        lc->y -= nc->y;
        
        div = lc->x * lc->x + lc->y * lc->y >= OUT2;
        conv = nc->x * nc->x + nc->y * nc->y < e2;
    }
    
    // divergence is very rare, it is faster to check it here
    dx = lc->x - sx;
    dy = lc->y - sy;
    
    div = dx * dx + dy * dy > m2 || lc->x * lc->x + lc->y * lc->y >= OUT2;
    
    if(div || ! conv) {
        v->x = lc->x;
        v->y = lc->y;
        
        lastIter = totNewton - tn;
        
        return false;
    }
    
    for (int i = 0; i < refine; i++) {
        ldbl d2 = mandel_miss_ntl(nc, lc, t, pp, per);
        if(d2 < 1E-25) {
            return false;
        }
        
        lc->x -= nc->x;
        lc->y -= nc->y;
    }
    
    v->x = lc->x;
    v->y = lc->y;
    
    lastIter = totNewton - tn;
    
    return true;
}

bool mandel_is_hypl(fp80 root, int per, ldbl error) {
    fp80 lc = {root->x, root->y};
    ldbl m2 = lc->x * lc->x + lc->y * lc->y;
    if(per < 1 || m2 >= 4 || error < 0) { // incorrect parameters or outside of D(0, 2)
        return false;
    }
    
    // check period 1
    ldbl e2 = error * error;
    if(per == 1 || m2 <= e2) {
        return per == 1 && m2 <= e2;
    }
    
    // check period 2
    ldbl dm1 = hypot(lc->x + 1, lc->y);
    if(per == 2 || dm1 <= error) {
        return per == 2 && dm1 <= error;
    }
    
    // attempt a formal proof, within the limits of fp80Disk
    fp80d c0;
    fp80d_setUlp(c0, lc->x, lc->y); // minimal error
    fp80d v0;
    fp80d_set(v0, c0);
    
    fp80d c;
    fp80d_setl(c, lc->x, lc->y, error);
    fp80d v;
    fp80d_set(v, c);
    fp80d d;
    fp80d_setUlp(d, 1, 0);
    for (int i = 1; i < per - 1; i++) {
        fp80d_muld(d, d, v, 2);
        fp80d_addd(d, d, 1);
        
        fp80d_sqr(v, v);
        fp80d_add(v, v, c);
        
        if(fp80d_min_mod(v) <= 0) {
            return false;
        }
        
        fp80d_sqr(v0, v0);
        fp80d_add(v0, v0, c0);
    }
    
    fp80d_muld(d, d, v, 2);
    fp80d_addd(d, d, 1);
    
    fp80d_sqr(v, v);
    fp80d_add(v, v, c);
    
    if(fp80d_min_mod(v) > 0) {
        return false;
    }
    
    fp80d_sqr(v0, v0);
    fp80d_add(v0, v0, c0);
    
    // by Theorem 1 in [1]
    return error * fp80d_min_mod(d) > fp80d_max_mod(v0);
}

bool mandel_is_misl(fp80 z, int pp, int per, ldbl error) {
    ppFail = -1;
    perFail = -1;
    
    if(pp < 2 || per < 1 || error < 0) {
        return false;
    }
    
    ldbl x = z->x;
    ldbl y = z->y;
    ldbl m2 = x * x + y * y;
    
    // (two plus error)^2
    ldbl tpe2 = 2 + 1E-15 + error;
    tpe2 = tpe2 * tpe2;
    
    // (quarter minus error)^2
    ldbl qme2 = 0.25 - 1E-15 - error;
    qme2 = qme2 * qme2;
    
    // quick check for outside M or in main cardiodid or main disk
    if(m2 > tpe2 || m2 < qme2 || (x + 1) * (x + 1) + y * y < qme2) {
        return false;
    }
    
    // attempt a formal proof, within the limits of fp80d
    fp80d c0;
    fp80d_setUlp(c0, x, y); // minimal error
    fp80d v0, v0p;
    fp80d_set(v0, c0);

    int n = pp + per;
    fp80d c;
    fp80d_struct v[n];
    
    fp80d_setl(c, x, y, error);
    fp80d_set(v, c);
    fp80d d, dpp;
    fp80d_setUlp(d, 1, 0);
    
    for (int i = 1; i < n; i++) {
        fp80d_muld(d, d, v + (i - 1), 2);
        fp80d_addd(d, d, 1);
        
        fp80d_sqr(v + i, v + (i - 1));
        fp80d_add(v + i, v + i, c);
        
        if(fp80d_min_mod(v) <= 0) {
            ppFail = 0;
            perFail = i + 1;
            
            return false;
        }
        
        for (int j = 0; j < i; j++) {
            // only those disks should intersect
            if(j == pp - 1 && i == n - 1) {
                continue;
            }
            
            if(fp80d_intersect(v + j, v + i)) {
                ppFail = j + 1;
                perFail = i  - j;
                
                return false;
            }
        }
        
        fp80d_sqr(v0, v0);
        fp80d_add(v0, v0, c0);
        
        if(i == pp - 1) {
            fp80d_set(dpp, d);
            fp80d_set(v0p, v0);
        }
    }
    
    fp80d_sub(v0, v0, v0p);
    fp80d_sub(d, d, dpp);
    
    // by Theorem 1 in [1]
    int ok = error * fp80d_min_mod(d) > fp80d_max_mod(v0);
    
    if(! ok) {
        ppFail = pp;
        perFail = per;
    }
    
    return ok;
}

bool mandel_conv_npl(fp80 c, int per, ldbl rootError, ldbl radius) {
    if (radius <= 3 * rootError) {
         return false;
    }
    
    fp80d dc, dv, dd;
    fp80d_setl(dc, c->x, c->y, radius);
    fp80d_setl(dv, c->x, c->y, radius);
    fp80d_setl(dd, 1, 0, 0);
    
    for (int i = 1; i < per ; i++) {
        // derivative d_n = 2 * v * d_{n-1} + 1
        fp80d_muld(dd, dd, dv, 2);
        fp80d_addd(dd, dd, 1);
        
        if(i < per - 1) {
            // new value v_n = v_{n-1}^2 + c
            fp80d_sqr(dv, dv);
            fp80d_add(dv, dv, dc);
        }
    }
    
    return fp80d_min_mod(dd) > 2 * dd->r;
}

bool mandel_conv_nppl(fp80 c, int pp, int per, ldbl rootError, ldbl radius, bool checkMultiple) {
    if (radius <= 3 * rootError || pp <= 1 || per < 1) {
         return false;
    }
    
    int n = pp + per;
    fp80d dc, d;
    fp80d_struct dv[n - 1], dd[n - 1]; // at position s, p_{s+1} and p'_{s+1}
    
    fp80d_setl(dc, c->x, c->y, radius);
    fp80d_setl(dv, c->x, c->y, radius);
    fp80d_setl(dd, 1, 0, 0);
    
    for (int i = 1; i < n - 1 ; i++) {
        // derivative d_i = 2 * v * d_{i-1} + 1
        fp80d_muld(dd + i, dd + (i - 1), dv + (i - 1), 2);
        fp80d_addd(dd + i, dd + i, 1);
        
        // new value v_i = v_{i-1}^2 + c
        fp80d_sqr(dv + i, dv + (i - 1));
        fp80d_add(dv + i, dv + i, dc);
    }
    
    // compute d = p'_n
    fp80d_muld(d, dd + (n - 2), dv + (n - 2), 2);
    fp80d_addd(d, d, 1);
    
    fp80d_sub(d, d, dd + (pp - 1));
    ldbl mm = fp80d_min_mod(d);
    
    if(! checkMultiple || mm > 0) {
        return mm > 4 * d->r;
    }
    
    // quite probably multiple root, check requested
    
    // find k, the order of the hyperbolic root
    int k = n;
    for (int i = 0; i < k; i++) {
        if(fp80d_min_mod(dv + i) == 0) {
            k = i + 1;
        }
    }
    
    // if not a divisor of per, radius is simply too large
    if(per % k != 0) {
        return false;
    }
    
    ldbl e = 8 * dd[per - 1].r / fp80d_min_mod(dd + (per - 1));
    fp80d sd = {{0, 0, 0, 0}}, buf;
    for (int i = 0; i < pp - 1; i++) {
        fp80d_add(d, dd + i, dd + (i + per));  // d = (p'_{per+i+1} + p'_{i+1})(B(c, radius))
        
        if((i + 1) % k == 0) {
            e += 4 * d->r / fp80d_min_mod(d);
        } else {
            fp80d_add(buf, dv + i, dv + (i + per));  // buf = (p_{per+i+1} + p_{i+1})(B(c, radius))
            if(! fp80d_div(buf, d, buf)) { // should not contain 0, radius too large?
                return false;
            }
            fp80d_add(sd, sd, buf);
        }
    }
    
    ldbl t = 2 * radius * fp80d_max_mod(sd);
    
    return t + e < (pp - 1) / k + 2;
}

bool mandel_is_prob_hypl(fp80 root, int per, ldbl error, int refine) {
    fp80 lc = {root->x, root->y};
    ldbl m2 = lc->x * lc->x + lc->y * lc->y;
    if(per < 1 || m2 >= 4 || error < 0) { // incorrect parameters or outside of D(0, 2)
        return false;
    }
    
    // check period 1
    ldbl e2 = error * error;
    if(per == 1 || m2 <= e2) {
        return per == 1 && m2 <= e2;
    }
    
    // check period 2
    ldbl dm1 = hypotl(lc->x + 1, lc->y);
    if(per == 2 || dm1 <= error) {
        return per == 2 && dm1 <= error;
    }
    
    // use approximations
    fp80 v, d;
    mandel_val_derl(v, d, lc, per);
    ldbl v2 = v->x * v->x + v->y * v->y;
    ldbl d2 = d->x * d->x + d->y * d->y;
    if(v2 > 2 * e2 * d2) {
        return false;
    }
    
    if(refine && v2 >= 1E-26) {
        fp80 c = {lc->x, lc->y};
        mandel_rootl(c, c, per, 40, error);
        
        if(hypotl(c->x - lc->x, c->y - lc->y) > error) {
            return false;
        }
        
        *lc = *c;
    }
    
    for (int i = 2; i <= per / 2; i++) {
        if(per % i != 0) {
            continue;
        }
        
        mandel_val_derl(v, d, lc, i);
        ldbl v2 = v->x * v->x + v->y * v->y;
        ldbl d2 = d->x * d->x + d->y * d->y;
        
        if(v2 <= 2 * e2 * d2) {
            return false;
        }
    }
    
    return true;
}

bool mandel_eql(ldbl x, ldbl y, ldbl nx, ldbl ny, ldbl eps) {
    return x >= nx - eps && x <= nx + eps && y >= ny - eps && y <= ny + eps;
}

ldbl* mandel_coeffs_hypl(int per, int deg) {
    if(deg < 0 || per < 1 || per > 50) {
        return NULL;
    }
    
    ldbl *cf = malloc((deg + 1) * sizeof(ldbl));
    
    for (int i = 0; i <= deg; i++) {
        cf[i] = 0;
    }
    cf[1] = 1;
    
    if(per == 1) {
        return cf;
    }
    
    cf[2] = 1;
    
    if(per == 2) {
        return cf;
    }
    
    ldbl *ncf = malloc((deg + 1) * sizeof(ldbl));
    for(int i = 0; i <= deg; i++) {
        ncf[i] = cf[i];
    }
    
    for (int p = 3; p <= per; p++) {
        long dp = 1L << (p - 1);
        
        for (int k = 3; k <= deg && k <= dp; k++) {
            ldbl c = cf[k - 1];
            
            int m = (k - 1) / 2;
            for (int j = 2; j <= m; j++) {
                c += cf[j] * cf[k - j];
            }
            
            c *= 2;
            
            if(k % 2 == 0) {
                ldbl cm = cf[m + 1];
                c += cm * cm;
            }
            
            ncf[k] = c;
        }
        
        ldbl *b = cf;
        cf = ncf;
        ncf = b;
    }
    
    free(ncf);
    
    return cf;
}

ldbl* mandel_coeffs_misl(int pp, int per, int deg) {
    if(deg < 0 || per < 1 || pp + per > 50 || pp < 0) {
        return NULL;
    }
    
    if( pp == 0 ) {
        return mandel_coeffs_hypl(per, deg);
    }
        
    ldbl* cf = mandel_coeffs_hypl(pp + per, deg);
    ldbl* ncf = mandel_coeffs_hypl(pp, deg);
    
    for (int i = 0; i <= deg; i++) {
        cf[i] -= ncf[i];
    }
    
    free(ncf);
    
    return cf;
}

ldbl* mandel_sum_neg_pows_hypl(int per, int deg) {
    if(deg < 1 || per < 1 || per > 50) {
        return NULL;
    }
    
    ldbl *pk = malloc(deg * sizeof(ldbl));
    ldbl *cf = mandel_coeffs_hypl(per, deg + 1);
    
    ldbl v;
    for (int k = 1; k <= deg; k++) {
        v = -k * cf[k + 1];
        for (int i = 1; i < k; i++) {
            v -= pk[i - 1] * cf[k - i + 1];
        }
        
        pk[k - 1] = v;
    }
    
    free(cf);
    
    return pk;
}

ldbl* mandel_sum_neg_pows_misl(int pp, int per, int deg) {
    if(deg < 1 || per < 1 || pp + per > 50 || pp < 0) {
        return NULL;
    }
    
    if( pp == 0 ) {
        return mandel_sum_neg_pows_hypl(per, deg);
    }
    
    int mult = pp + 1; // zero has multiplicity pp + 1
    ldbl *pk = malloc(deg * sizeof(ldbl));
    ldbl *cf = mandel_coeffs_misl(pp, per, deg + mult);
    ldbl prod = cf[mult]; // product of non-zero roots
    
    ldbl v;
    for (int k = 1; k <= deg; k++) {
        v = -k * cf[k + mult];
        for (int i = 1; i < k; i++) {
            v -= pk[i - 1] * cf[k - i + mult];
        }
        
        pk[k - 1] = v / prod;
    }
    
    free(cf);
    
    return pk;
}

// MARK: arbitrary precision versions

long mandel_dist(mpfr_t dist, mpc c, long maxIter) {
    lastIter = 0;
    
    long prec = mpfr_get_prec(dist);
    defs_mpc(prec, bcV, bcD);
    mpc_set(bcV, c);
    mpc_setd(bcD, 1, 0);

    int small = 1;
    long k = 1, expV;
    for (; k < maxIter + 7 && small; k++) {
        mpc_mul(bcD, bcD, bcV);
        mpc_scale(bcD, bcD, 1);
        mpc_addi(bcD, bcD, 1);

        mpc_sqr(bcV, bcV);
        mpc_add(bcV, bcV, c);

        expV = mpc_2exp(bcV);
        small = expV < 1000;
        lastIter = expV < 5 ? k + 1 : lastIter;        // for Carleson's conjecture, escape from D(0, 5)
    }
    
    // allocate 4 buffers on the stack
    defs_mpfr(prec + MPC_EXTRA_PREC, b1, b2, b3, b4);
    
    mpc_mod2(b4, bcD);                                 // b4 stores |d|^2
    if(small || mpfr_cmp_d(b4, 1) <= 0) {              // bounded
        mpfr_set_zero(dist, 0);

        return -1;
    }
    
    mpc_mod2(b1, bcV);                       // b1 stores |v|^2
    mpfr_div(b2, b1, b4, MPFR_RNDD);         // b2 stores |v|^2 / |d|^2
    mpfr_sqrt(b2, b2, MPFR_RNDD);            // b2 stores |v| / |d|
    
    mpfr_mul_2si(b1, b1, -2, MPFR_RNDD);     // b1 stores |v|^2 / 4
    mpfr_log(b3, b1, MPFR_RNDD);             // b3 stores log(|v|^2 / 4)
    
    mpfr_mul_2si(b3, b3, -k, MPFR_RNDD);     // b3 stores log(|v|^2 / 4) / 2^k
    mpfr_neg(b3, b3, MPFR_RNDU);             // b3 stores -log(|v|^2 / 4) / 2^k
    mpfr_expm1(b3, b3, MPFR_RNDU);           // b3 stores exp(-log(|v|^2 / 4) / 2^k) - 1
    mpfr_mul_2si(b3, b3, k - 2, MPFR_RNDU);  // b3 stores 2^{k-2} (exp(-log(|v|^2 / 4) / 2^k) - 1)
    mpfr_neg(b3, b3, MPFR_RNDD);             // b3 stores 2^{k-2} (1 - exp(-log(|v|^2 / 4) / 2^k))
    
    mpfr_mul(dist, b2, b3, MPFR_RNDD);       // finally, the distance to M
    
    return k;
}

void mandel_val(mpc v, mpc c, int per) {
    def_mpc(mpc_prec(v), bcV);
    mpc_set(bcV, c);
    
    for (int k = 1; k < per; k ++) {
        mpc_sqr(bcV, bcV);
        mpc_add(bcV, bcV, c);
    }
    
    mpc_set(v, bcV);
}

void mandel_mis_val(mpc v, mpc c, int pp, int per) {
    defs_mpc(mpc_prec(v), bcV, bcV0);
    mpc_set(bcV, c);
    mpc_set0(bcV0);
    
    int n = pp + per;
    for (int k = 1; k < n; k ++) {
        if(k == pp) {
            mpc_set(bcV0, bcV);
        }
        
        mpc_sqr(bcV, bcV);
        mpc_add(bcV, bcV, c);
    }
    
    mpc_sub(v, bcV, bcV0);
}

void mandel_val_der(mpc v, mpc d, mpc c, long per) {
    defs_mpc(mpc_prec(v), bcV, bcD);
    mpc_set(bcV, c);
    mpc_seti(bcD, 1, 0);
    
    for (long k = 1; k < per; k ++) {
        mpc_mul(bcD, bcD, bcV);
        mpc_scale(bcD, bcD, 1);
        mpc_addi(bcD, bcD, 1);
        
        mpc_sqr(bcV, bcV);
        mpc_add(bcV, bcV, c);
    }
    
    mpc_set(v, bcV);
    mpc_set(d, bcD);
}

void mandel_mis_val_der(mpc v, mpc d, mpc c, int pp, int per) {
    defs_mpc(mpc_prec(v), bcV, bcD, bcV0, bcD0);
    mpc_set(bcV0, c);
    mpc_setd(bcD0, 1, 0);
    
    int n = pp + per;
    for (int k = 1; k < n; k ++) {
        mpc_mul(bcD0, bcD0, bcV0);
        mpc_scale(bcD0, bcD0, 1);
        mpc_addi(bcD0, bcD0, 1);
        
        mpc_sqr(bcV0, bcV0);
        mpc_add(bcV0, bcV0, c);
        
        if(k == pp - 1) {
            mpc_set(bcV, bcV0);
            mpc_set(bcD, bcD0);
        }
    }
    
    mpc_sub(v, bcV0, bcV);
    mpc_sub(d, bcD0, bcD);
}

void mandel_miss_val_der(mpc v, mpc d, mpc c, int pp, int per) {
    defs_mpc(mpc_prec(v), bcV, bcD, bcV0, bcD0);
    mpc_set(bcV0, c);
    mpc_setd(bcD0, 1, 0);
    
    if(pp == 2) {
        mpc_set(bcV, bcV0);
        mpc_set(bcD, bcD0);
    }
    
    int n = pp + per - 1;
    for (int k = 1; k < n; k ++) {
        mpc_mul(bcD0, bcD0, bcV0);
        mpc_scale(bcD0, bcD0, 1);
        mpc_addi(bcD0, bcD0, 1);
        
        mpc_sqr(bcV0, bcV0);
        mpc_add(bcV0, bcV0, c);
        
        if(k == pp - 2) {
            mpc_set(bcV, bcV0);
            mpc_set(bcD, bcD0);
        }
    }
    
    mpc_add(v, bcV0, bcV);
    mpc_add(d, bcD0, bcD);
}

bool mandel_nt(mpc v, mpc c, int per) {
    defs_mpc(mpc_prec(v), bcV, bcD);
    mpc_set(bcV, c);
    mpc_seti(bcD, 1, 0);
    
    int small = 1;
    long k = 1, ce = mpc_2exp(c);
    
    long prec = mpc_prec(v);
    for (; k < per && small; k++) {
        mpc_mul(bcD, bcD, bcV);
        mpc_scale(bcD, bcD, 1);
        mpc_addi(bcD, bcD, 1);
        
        mpc_sqr(bcV, bcV);
        mpc_add(bcV, bcV, c);
        
        small = mpc_2exp(bcV) <= prec + ce || mpc_2exp(bcD) <= prec + 1;
    }
    
    mpc_scale(bcD, bcD, per - k);
    bool ok = mpc_div(v, bcV, bcD);
    
    totNewton ++;
    
    return ok;
}

// does not compute correctly outside the scale of MPFR numbers
// uses bcV0 and bcD0
bool mandel_mis_nt(mpc v, mpc c, int pp, int per) {
    defs_mpc(mpc_prec(v), bcV, bcD);
    
    mandel_mis_val_der(bcV, bcD, c, pp, per);
    bool ok = mpc_div(v, bcV, bcD);
    
    totNewton ++;
    
    return ok;
}

// does not compute correctly outside the scale of MPFR numbers
// uses bcV0 and bcD0
bool mandel_mis_nt_sol(mpc v, mpc c, mpc t, int pp, int per) {
    defs_mpc(mpc_prec(v), bcV, bcD);
        
    mandel_mis_val_der(bcV, bcD, c, pp, per);
    if(t != NULL) {
        mpc_sub(bcV, bcV, t);
    }
        
    bool ok = mpc_div(v, bcV, bcD);
    
    totNewton ++;
    
    return ok;
}

// does not compute correctly outside the scale of MPFR numbers
// uses bcV, bcD, bcV0, bcD0, bcV1, bcD1
bool mandel_miss_nt_sol(mpc v, mpc c, mpc t, int pp, int per) {
    defs_mpc(mpc_prec(v), bcV, bcD);
        
    mandel_miss_val_der(bcV, bcD, c, pp, per);
    if(t != NULL) {
        mpc_sub(bcV, bcV, t);
    }
        
    bool ok = mpc_div(v, bcV, bcD);
    
    totNewton ++;
    
    return ok;
}

bool mandel_root(mpc v, mpc c, int per, int maxIter) {
    return mandel_root_ref(v, c, per, maxIter < 2 ? 2 : maxIter, 1);
}

bool mandel_root_ref(mpc v, mpc c, int per, int maxIter, int refine) {
    defs_mpc(mpc_prec(v), bcV, bcNt);
    long exp = 0;
    mpc_set(bcV, c);
    
    ulong tn = totNewton;
    int div = 0, conv = 0;
    
    long sc = mpc_prec(v) * 3 / 4;
    for (int i = 0; i < maxIter && ! div && ! conv; i++) {
        mandel_nt(bcNt, bcV, per);
        mpc_sub(bcV, bcV, bcNt);
        
        div = mpc_2exp(bcV) >= 8;
        exp = mpc_2exp(bcNt);
        conv = -exp > sc;
    }
    
    if(div || ! conv) {
        lastIter = totNewton - tn;
        
        return false;
    }
    
    if(-exp > mpc_prec(v) * 9 / 10) { // optimization: good enough, cannot hope for much more
        mpc_set(v, bcV);
        
        return true;
    }
    
    for (int i = 0; i < refine; i++) {
        mandel_nt(bcNt, bcV, per);
        mpc_sub(bcV, bcV, bcNt);
    }
    
    mpc_set(v, bcV);
    
    lastIter = totNewton - tn;
    
    return true;
}

// uses bcV and bcD
bool mandel_nt_sol(mpc v, mpc c, mpc t, long per) {
    defs_mpc(mpc_prec(v), bcV, bcD);
    mpc_set(bcV, c);
    mpc_seti(bcD, 1, 0);
    
    long ce = mpc_2exp(c);
    long te = mpc_2exp(t);
    long exp = ce >= te ? ce : te;
    
    bool small = true;
    long k = 1;
    
    long prec = mpc_prec(v);
    for (; k < per && small; k++) {
        mpc_mul(bcD, bcD, bcV);
        mpc_scale(bcD, bcD, 1);
        mpc_addi(bcD, bcD, 1);
        
        mpc_sqr(bcV, bcV);
        mpc_add(bcV, bcV, c);
        
        small = mpc_2exp(bcV) <= prec + exp || mpc_2exp(bcD) <= prec + 1;
    }
    
    if(small) {
        mpc_sub(bcV, bcV, t);
    }
    
    if(per > k) {
        mpc_scale(bcD, bcD, per - k);
    }
    bool ok = mpc_div(v, bcV, bcD);
    
    totNewton ++;
    
    return ok;
}

bool mandel_solution(mpc v, mpc c, mpc t, int per, int maxIter) {
    return mandel_sol_ref(v, NULL, c, t, per, maxIter, 1);
}

bool mandel_sol_ref(mpc v, mpfr_t err, mpc c, mpc t, long per, int maxIter, int refine) {
    defs_mpc(mpc_prec(v), bcV, bcNt);
    mpc_set(bcV, c);
    
    ulong tn = totNewton;
    int div = 0;
    int conv = 0;
    
    long sc = mpc_prec(v) * 3 / 4;
    for (int i = 0; i < maxIter && ! div && ! conv; i++) {
        mandel_nt_sol(bcNt, bcV, t, per);
        mpc_sub(bcV, bcV, bcNt);
        
        div = mpc_2exp(bcV) >= 8;
        conv = -mpc_2exp(bcNt) > sc;
    }
    
    if(div || ! conv) {
        lastIter = totNewton - tn;
        
        return false;
    }
    
    for (int i = 0; i < refine; i++) {
        mandel_nt_sol(bcNt, bcV, t, per);
        mpc_sub(bcV, bcV, bcNt);
    }
    
    if(err != NULL) {
        mpc_mod(err, bcNt);
    }
    
    mpc_set(v, bcV);
    
    lastIter = totNewton - tn;
    
    return true;
}

bool mandel_mis_root_ref(mpc v, mpc c, mpc t, int pp, int per, int maxIter, int refine) {
    defs_mpc(mpc_prec(v), bcV, bcNt);
    mpc_set(bcV, c);
    
    ulong tn = totNewton;
    int div = 0;
    int conv = 0;
    
    long sc = mpc_prec(v) * 3 / 4;
    for (int i = 0; i < maxIter && ! div && ! conv; i++) {
        if(! mandel_mis_nt_sol(bcNt, bcV, t, pp, per)) {
            return false;
        }
        
        mpc_sub(bcV, bcV, bcNt);
        
        div = mpc_2exp(bcV) >= 8;
        conv = -mpc_2exp(bcNt) > sc;
    }
    
    if(div || ! conv) {
        lastIter = totNewton - tn;
        
        return false;
    }
    
    for (int i = 0; i < refine; i++) {
        if(! mandel_mis_nt_sol(bcNt, bcV, t, pp, per)) {
            return false;
        }
        
        mpc_sub(bcV, bcV, bcNt);
    }
    
    mpc_set(v, bcV);
    
    lastIter = totNewton - tn;
    
    return true;
}

bool mandel_miss_root_ref(mpc v, mpc c, mpc t, int pp, int per, int maxIter, int refine) {
    defs_mpc(mpc_prec(v), bcV, bcNt);
    mpc_set(bcV, c);
    
    ulong tn = totNewton;
    int div = 0;
    int conv = 0;
    
    long sc = mpc_prec(v) * 3 / 4;
    for (int i = 0; i < maxIter && ! div && ! conv; i++) {
        if(! mandel_miss_nt_sol(bcNt, bcV, t, pp, per)) {
            return false;
        }
        
        mpc_sub(bcV, bcV, bcNt);
        
        div = mpc_2exp(bcV) >= 8;
        conv = -mpc_2exp(bcNt) > sc;
    }
    
    if(div || ! conv) {
        lastIter = totNewton - tn;
        
        return false;
    }
    
    for (int i = 0; i < refine; i++) {
        if(! mandel_miss_nt_sol(bcNt, bcV, t, pp, per)) {
            return false;
        }
        
        mpc_sub(bcV, bcV, bcNt);
    }
    
    mpc_set(v, bcV);
    
    lastIter = totNewton - tn;
    
    return true;
}

mpfr_ptr mandel_coeffs_hyp(int per, int deg, int prec) {
    if(deg < 0 || per < 1 || per > 1000) {
        return NULL;
    }
    
    mpfr_ptr cf = malloc((deg + 1) * sizeof(__mpfr_struct));
    
    for (int i = 0; i <= deg; i++) {
        mpfr_init2(cf + i, prec);
        mpfr_set_si(cf + i, 0, MPFR_RNDN);
    }
    mpfr_set_si(cf + 1, 1, MPFR_RNDN);
    
    if(per == 1) {
        return cf;
    }
    
    mpfr_set_si(cf + 2, 1, MPFR_RNDN);
    
    if(per == 2) {
        return cf;
    }
    
    mpfr_ptr ncf = malloc((deg + 1) * sizeof(__mpfr_struct));
    for (int i = 0; i <= deg; i++) {
        mpfr_init2(ncf + i, prec);
        mpfr_set(ncf + i, cf + i, MPFR_RNDN);
    }
    
    mpfr_t c, cm;
    mpfr_init2(c, prec);
    mpfr_init2(cm, prec);
    
    for (int p = 3; p <= per; p++) {
        long dp = 1L << (p - 1);
        
        for (int k = 3; k <= deg && k <= dp; k++) {
            mpfr_set(c, &cf[k - 1], MPFR_RNDN);
            
            int m = (k - 1) / 2;
            for (int j = 2; j <= m; j++) {
                mpfr_mul(cm, &cf[j], &cf[k - j], MPFR_RNDN);
                
                mpfr_add(c, c, cm, MPFR_RNDN);
            }
            
            mpfr_mul_2si(c, c, 1, MPFR_RNDN);
            
            if(k % 2 == 0) {
                mpfr_set(cm, &cf[m + 1], MPFR_RNDN);
                mpfr_sqr(cm, cm, MPFR_RNDN);
                
                mpfr_add(c, c, cm, MPFR_RNDN);
            }

            mpfr_set(&ncf[k], c, MPFR_RNDN);
        }
        
        mpfr_ptr b = cf;
        cf = ncf;
        ncf = b;
    }
    
    for (int i = 0; i <= deg; i ++) {
        mpfr_clear(ncf + i);
    }
    free(ncf);
    
    return cf;
}

mpfr_ptr mandel_coeffs_mis(int pp, int per, int deg, int prec) {
    if(deg < 0 || per < 1 || pp + per > 1000 || pp < 0) {
        return NULL;
    }
    
    if( pp == 0 ) {
        return mandel_coeffs_hyp(per, deg, prec);
    }
    
    mpfr_ptr cf = mandel_coeffs_hyp(pp + per, deg, prec);
    mpfr_ptr ncf = mandel_coeffs_hyp(pp, deg, prec);
    
    for (int i = 0; i <= deg; i++) {
        mpfr_sub(cf + i, cf + i, ncf + i, MPFR_RNDN);
    }
    
    for (int i = 0; i <= deg; i ++) {
        mpfr_clear(ncf + i);
    }
    free(ncf);
    
    return cf;
}

static inline void mandel_free_mpfr_list(mpfr_ptr v, int len) {
    for (int i = 0; i < len; i ++) {
        mpfr_clear(v + i);
    }
    free(v);
}

static inline void mandel_mpfr_list_sub(mpfr_ptr v, mpfr_ptr d, int len) {
    for (int i = 0; i < len; i++) {
        mpfr_sub(v + i, v + i, d + i, MPFR_RNDN);
    }
}

static inline void mandel_mpfr_list_sub_mult(mpfr_ptr v, mpfr_ptr d, int len, int m) {
    if(m == 0) {
        return;
    }
    
    if(m == 1) {
        mandel_mpfr_list_sub(v, d, len);
        
        return;
    }
    
    for (int i = 0; i < len; i++) {
        mpfr_mul_si(d + i, d + i, m, MPFR_RNDN);
        mpfr_sub(v + i, v + i, d + i, MPFR_RNDN);
    }
}

mpfr_ptr mandel_sum_all_neg_pows_hyp(int per, int deg, int prec) {
    if(deg < 1 || per < 1 || per > 1000) {
        return NULL;
    }
    
    mpfr_ptr cf = mandel_coeffs_hyp(per, deg + 1, prec);
    mpfr_ptr pk = malloc(deg * sizeof(__mpfr_struct));
    for (int i = 0; i < deg; i++) {
        mpfr_init2(pk + i, prec);
    }
    
    mpfr_t v;
    mpfr_init2(v, prec);
    for (int k = 1; k <= deg; k++) {
        mpfr_set_si(&pk[k - 1], -k, MPFR_RNDN);
        mpfr_mul(&pk[k - 1], &pk[k - 1], &cf[k + 1], MPFR_RNDN);
        
        for (int i = 1; i < k; i++) {
            mpfr_mul(v, &pk[i - 1], &cf[k - i + 1], MPFR_RNDN);
            mpfr_sub(&pk[k - 1], &pk[k - 1], v, MPFR_RNDN);
        }
    }
    
    mpfr_clear(v);
    mandel_free_mpfr_list(cf, deg + 1);
    
    return pk;
}

mpfr_ptr mandel_sum_neg_pows_hyp(int per, int deg, int prec) {
    mpfr_ptr pk = mandel_sum_all_neg_pows_hyp(per, deg, prec);
    if(pk == NULL) {
        return NULL;
    }
    
    for (int d = 2; d <= per / 2; d++) {
        if(per % d != 0) {
            continue;
        }
        
        mpfr_ptr pd = mandel_sum_neg_pows_hyp(d, deg, prec);
        mandel_mpfr_list_sub(pk, pd, deg);
        mandel_free_mpfr_list(pd, deg);
    }
    
    return pk;
}

mpfr_ptr mandel_sum_all_neg_pows_mis(int pp, int per, int deg, int prec) {
    if(deg < 1 || per < 1 || pp + per > 1000 || pp < 0) {
        return NULL;
    }
    
    if( pp == 0 ) {
        return mandel_sum_all_neg_pows_hyp(per, deg, prec);
    }
    
    int mult = pp + 1; // zero has multiplicity pp + 1
    mpfr_ptr cf = mandel_coeffs_mis(pp, per, deg + mult, prec);
    mpfr_ptr pk = malloc(deg * sizeof(__mpfr_struct));
    for (int i = 0; i < deg; i++) {
        mpfr_init2(pk + i, prec);
    }
    
    mpfr_t prod; // product of non-zero roots
    mpfr_init2(prod, prec);
    mpfr_set(prod, cf + mult, MPFR_RNDN);
    
    mpfr_t v, t;
    mpfr_inits2(prec, v, t, NULL);
    for (int k = 1; k <= deg; k++) {
        mpfr_mul_si(v, cf + (k + mult), -k, MPFR_RNDN);
        for (int i = 1; i < k; i++) {
            mpfr_mul(t, pk + (i - 1), cf + (k - i + mult), MPFR_RNDN);
            mpfr_sub(v, v, t, MPFR_RNDN);
        }

        mpfr_div(pk + (k - 1), v, prod, MPFR_RNDN); //FIXME: bad index
    }
    
    mpfr_clears(v, t, NULL);
    mandel_free_mpfr_list(cf, deg + 1);
    
    return pk;

}

/// @brief Multiplicity of a hyperbolic factor in Misiurewicz polynomials
///
/// @param pp the pre-period
/// @param per the period
/// @param hyp the hyperbolic factor
///
/// @return the multiplicity
///
/// Returns the multiplicity of a hyperbolic factor @c Hyp(hyp) in the polynomial \f$ m_[pp, per}.\f$
/// All hyperbolic centers of period @c hyp share the same multiplicity.
///
/// There is no upper limit on @c per, aside from @c int type limitation to 32767.
static int mandel_mis_factor_hyp_mult(int pp, int per, int hyp) { // private
    if ( pp < 0 || per <= 0 || hyp <= 0) // undefined cases
        return 0;
    
    if ( per % hyp != 0 ) // hyp is not a divisor of per
        return 0;
    
    if ( pp < 2 || hyp == 1 ) { // hyperbolic cases pp=0, pp=1 or divisibility by z (Hyp(1) == {0})
        return( pp + 1 );
    } else { // see our result on the factorization of Misiurewicz polynomials
        return 2 + (pp - 1) / hyp; // integer division implements floor
    }
}

mpfr_ptr mandel_sum_neg_pows_mis(int pp, int per, int deg, int prec) {
    mpfr_ptr pk = mandel_sum_all_neg_pows_mis(pp, per, deg, prec);
    if(pk == NULL) {
        return NULL;
    }
    
    if(pp > 2) {
        mpfr_ptr pd = mandel_sum_all_neg_pows_mis(pp - 1, per, deg, prec);
        mandel_mpfr_list_sub(pk, pd, deg);
        mandel_free_mpfr_list(pd, deg);
    }
    
    for (int d = 1; d <= per; d++) {
        if(per % d != 0) {
            continue;
        }
        
        if(d < per) {
            mpfr_ptr pdiv = mandel_sum_neg_pows_mis(pp, d, deg, prec);
            mandel_mpfr_list_sub(pk, pdiv, deg);
            mandel_free_mpfr_list(pdiv, deg);
        }
        
        if(d == 1) {
            continue;
        }
        
        int mult = mandel_mis_factor_hyp_mult(pp, per, d);
        if(pp > 2) {
            mult -= mandel_mis_factor_hyp_mult(pp - 1, per, d);
        }
        
        if(mult == 0) {
            continue;
        }
        
        mpfr_ptr pd = mandel_sum_neg_pows_hyp(d, deg, prec);
        mandel_mpfr_list_sub_mult(pk, pd, deg, mult);
        mandel_free_mpfr_list(pd, deg);
    }
    
    return pk;
}

// MARK: proofs with disk arithmetics

bool mandel_is_hyp(mpc root, int per, ldbl error) {
    long prec = mpc_prec(root);
    long precR = prec <= 128 ? 61 - MPD_EXTRA_PREC : 113 - MPD_EXTRA_PREC;
    defs_mpd(prec, precR, bdC, bdC0, bdV, bdV0, bdD);
    defs_mpfr(prec + MPD_EXTRA_PREC, bfE, bfV, bfD);
    
    mpfr_set_ld(bfE, error, MPFR_RNDU);
    mpd_setr(bdC, root, bfE);
    
    // check for low periods
    int cont0 = mpd_contains_0(bdC, 0);
    if(per == 1 || cont0) {
        return per == 1 && cont0;
    }
    
    fp80 m1 = {-1, 0};
    int cm1 = mpd_contains_c80(bdC, m1, 0);
    if(per == 2 || cm1) {
        return per == 2 && cm1;
    }
    
    // compute iterates and check if it is a root of lower period
    mpd_set_center(bdC0, bdC);          // precise parameter
    mpd_set(bdV0, bdC0);               // precise iterates
    
    mpd_set(bdV, bdC);                // iterates in the disk c
    mpd_setd(bdD, 1, 0, 0);           // derivatives in the disk c
    
    for (int i = 1; i < per - 1; i++) {
        // derivative d_n = 2 * v * d_{n-1} + 1
        mpd_mul(bdD, bdD, bdV);
        mpd_scale(bdD, bdD, 1);
        mpd_add_si(bdD, bdD, 1);
        
        // new value v_n = v_{n-1}^2 + c
        mpd_sqr(bdV, bdV);
        mpd_add(bdV, bdV, bdC);
        
        if(mpd_contains_0(bdV, 0) || ! mpd_valid(bdV)) {
            return false;
        }
        
        // precise value
        mpd_sqr(bdV0, bdV0);
        mpd_add(bdV0, bdV0, bdC0);
    }
    
    // compute the last iterate; if it does not contain 0, root cannot be valid (a root)
    mpd_mul(bdD, bdD, bdV);
    mpd_scale(bdD, bdD, 1);
    mpd_add_si(bdD, bdD, 1);
    
    mpd_sqr(bdV, bdV);
    mpd_add(bdV, bdV, bdC);
    
    if(! mpd_contains_0(bdV, 1) || ! mpd_valid(bdV) || ! mpd_valid(bdD)) {
        return false;
    }
    
    mpd_sqr(bdV0, bdV0);
    mpd_add(bdV0, bdV0, bdC0);
    
    if(! mpd_valid(bdV0)) {
        return false;
    }
    
    // use the argument principle to prove that a root is in the disk D(root, error)
    // see Theorem 1 in [1]
    mpd_min_mod(bfD, bdD);
    mpfr_mul(bfD, bfD, bfE, MPFR_RNDZ);
    mpd_max_mod(bfV, bdV0);
    
    if(! mpfr_number_p(bfV) || ! mpfr_number_p(bfE)) {
        return false;
    }
        
    return mpfr_cmp(bfD, bfV) > 0;
}

bool mandel_is_mis(mpc c, int pp, int per, ldbl error) {
    if(pp < 2 || per < 1 || error < 0) {
        return false;
    }
    
    // modulus of c
    ldbl mc = mpc_modl(c);
    // c in fp80, minus 1
    fp80 c80, m1 = {-1, 0};
    mpc_get80(c80, c);
    // distance to minus 1
    ldbl dm1 = fp80_dist(c80, m1);
    // qaurter minus error
    ldbl qme = 0.0625L - 1E-15L - error;
    
    // quick check for outside M or in main cardiodid or main disk
    if(mc > 2 + 1E-15L + error || mc < qme || dm1 < qme) {
        return false;
    }
    
    long prec = mpc_prec(c);
    long precR = prec <= 128 ? 61 - MPD_EXTRA_PREC : 113 - MPD_EXTRA_PREC;
    defs_mpd(prec, precR, bdC, bdC0, bdV, bdV0, bdD, bdD0);
    defs_mpfr(prec + MPD_EXTRA_PREC, bfE, bfV, bfD);
    
    mpfr_set_ld(bfE, error, MPFR_RNDU);
    mpd_setr(bdC, c, bfE);
        
    // compute iterates and check if it is a root of lower period
    mpd_set_center(bdC0, bdC);          // precise parameter
    mpd_set(bdV0, bdC0);               // precise iterates
    
    mpd_setd(bdD, 1, 0, 0);           // derivatives in the disk c
        
    int n = pp + per;
    ulong size = mpd_limbs(prec, precR);
    ulong vlimbs[size * n];
    mpd_struct mv[n];     // the list of disks that should be disjoint
    for (int i = 0; i < n; i++) {
        mpd_iniz(mv + i, prec, precR, vlimbs + size * i);
    }
    
    mpd_set(mv, bdC);               // iterates in the disk c
    
    for (int i = 1; i < n; i++) {
        // derivative d_n = 2 * v * d_{n-1} + 1
        mpd_mul(bdD, bdD, mv + (i - 1));
        mpd_scale(bdD, bdD, 1);
        mpd_add_si(bdD, bdD, 1);
        
        // new value v_n = v_{n-1}^2 + c
        mpd_sqr(mv + i, mv + (i - 1));
        mpd_add(mv + i, mv + i, bdC);
        
        if(mpd_contains_0(mv + i, 0) || ! mpd_valid(mv + i)) {
            return false;
        }
        
        for (int j = 0; j < i; j++) {
            // only those disks should intersect
            if(j == pp - 1 && i == n - 1) {
                continue;
            }
            
            if(mpd_intersect(mv + j, mv + i, 0)) {
                return false;
            }
        }
        
        // precise value
        mpd_sqr(bdV0, bdV0);
        mpd_add(bdV0, bdV0, bdC0);
        
        if(i == pp - 1) {
            mpd_set(bdD0, bdD);
            mpd_set(bdV, bdV0);
        }
    }
    
    mpd_sub(bdV0, bdV0, bdV);
    mpd_sub(bdD, bdD, bdD0);
    
    // use the argument principle to prove that a root is in the disk D(root, error)
    // see Theorem 1 in [1]
    mpd_min_mod(bfD, bdD);
    mpfr_mul(bfD, bfD, bfE, MPFR_RNDZ);
    mpd_max_mod(bfV, bdV0);
    
    if(! mpfr_number_p(bfV) || ! mpfr_number_p(bfE)) {
        return false;
    }
        
    return mpfr_cmp(bfD, bfV) > 0;
}

bool mandel_conv_np(mpc root, int per, ldbl rootError, ldbl radius) {
    if (radius <= 3 * rootError) {
         return false;
    }
    
    long prec = mpc_prec(root);
    long precR = prec <= 128 ? 61 - MPD_EXTRA_PREC : 113 - MPD_EXTRA_PREC;
    defs_mpd(prec, precR, bdC, bdV, bdD);
    defs_mpfr(prec + MPD_EXTRA_PREC, bfE, bfV);
    
    mpfr_set_ld(bfE, radius, MPFR_RNDU);
    mpd_setr(bdC, root, bfE);
    
    // compute iterates and derivative
    mpd_set(bdV, bdC);                    // iterates in the disk c
    mpd_setd(bdD, 1, 0, 0);               // derivatives in the disk c
    
    for (int i = 1; i < per ; i++) {
        // derivative d_n = 2 * v * d_{n-1} + 1
        mpd_mul(bdD, bdD, bdV);
        mpd_scale(bdD, bdD, 1);
        mpd_add_si(bdD, bdD, 1);
        
        if(i < per - 1) {
            // new value v_n = v_{n-1}^2 + c
            mpd_sqr(bdV, bdV);
            mpd_add(bdV, bdV, bdC);
        }
    }
    
    // Check if p_per'(D(root,radius)) subset D(d,eps') such that |d| > 5 * eps'
    mpd_min_mod(bfV, bdD);                 // |d|-eps'
    mpfr_mul_2si(bfE, bdD->r, 2, MPFR_RNDU); // 4 * eps'
    
    if(! mpfr_number_p(bfV) || ! mpfr_number_p(bfE)) {
        return false;
    }
    
    // TODO: improvememnt as in the ldbl version?
    return mpfr_cmp(bfV, bfE) > 0;
}

bool mandel_conv_npp(mpc root, int pp, int per, ldbl rootError, ldbl radius, bool checkMultiple) {
    if (radius <= 3 * rootError || pp <= 1 || per < 1) {
         return false;
    }
    
    int n = pp + per;
    
    long prec = mpc_prec(root);
    long precR = prec <= 128 ? 61 - MPD_EXTRA_PREC : 113 - MPD_EXTRA_PREC;
    defv_mpd(dv, n - 1, prec, precR);
    defv2_mpd(dd, n - 1, prec, precR);
    defs_mpd(prec, precR, c, d, sd, buf);
    defs_mpfr(prec + MPD_EXTRA_PREC, mm, r, e);
    
    mpfr_set_ld(r, radius, MPFR_RNDU);
    mpd_setr(c, root, r);
    
    // compute iterates and derivative of p_{pp+per}
    mpd_set(dv, c);                      // iterates in the disk c
    mpd_setd(dd, 1, 0, 0);               // derivatives in the disk c
    
    for (int i = 1; i < n - 1 ; i++) {
        // derivative d_i = 2 * v * d_{i-1} + 1
        mpd_mul(dd + i, dd + (i - 1), dv + (i - 1));
        mpd_scale(dd + i, dd + i, 1);
        mpd_add_si(dd + i, dd + i, 1);
        
        // new value v_i = v_{i-1}^2 + c
        mpd_sqr(dv + i, dv + (i - 1));
        mpd_add(dv + i, dv + i, c);
    }
    
    // compute d = p'_n
    mpd_mul(d, dd + (n - 2), dv + (n - 2));
    mpd_scale(d, d, 1);
    mpd_add_si(d, d, 1);
    
    mpd_sub(d, d, dd + (pp - 1));
    
    // Check if mis'_{pp, per}(D(root,radius)) subset D(d,eps') such that |d| > 5 * eps'
    mpd_min_mod(mm, d);                    // |d|-eps'
    mpfr_mul_2si(r, d->r, 2, MPFR_RNDU);   // 4 * eps'
    
    if(! mpfr_number_p(mm) || ! mpfr_number_p(r)) {
        return false;
    }
    
    if(! checkMultiple || mpfr_cmp_si(mm, 0) > 0) {
        return mpfr_cmp(mm, r) > 0;
    }
    
    // quite probably multiple root, check requested
    
    // find k, the order of the hyperbolic root
    int k = n;
    for (int i = 0; i < k; i++) {
        mpd_min_mod(mm, dv + i);
        if(mpfr_zero_p(mm)) {
            k = i + 1;
        }
    }
    
    // if not a divisor of per, radius is simply too large
    if(per % k != 0) {
        return false;
    }
    
    mpd_min_mod(mm, dd + (per - 1));
    mpfr_mul_2si(r, dd[per - 1].r, 3, MPFR_RNDU);
    mpfr_div(e, r, mm, MPFR_RNDU);
    
    for (int i = 0; i < pp - 1; i++) {
        mpd_add(d, dd + i, dd + (i + per));  // d = (p'_{per+i+1} + p'_{i+1})(B(c, radius))
        
        if((i + 1) % k == 0) {
            mpd_min_mod(mm, d);
            mpfr_mul_2si(r, d->r, 2, MPFR_RNDU);
            mpfr_div(r, r, mm, MPFR_RNDU);
            mpfr_add(e, e, r, MPFR_RNDU);
        } else {
            mpd_add(buf, dv + i, dv + (i + per));  // buf = (p_{per+i+1} + p_{i+1})(B(c, radius))
            if(! mpd_div(buf, d, buf)) { // should not contain 0, radius too large?
                return false;
            }
            mpd_add(sd, sd, buf);
        }
    }
    
    mpd_max_mod(mm, sd);
    mpfr_mul_d(mm, mm, 2 * radius, MPFR_RNDU);
    mpfr_add(e, e, mm, MPFR_RNDU);
    
    return mpfr_cmp_si(e, (pp - 1) / k + 2) < 0;
}

bool mandel_univalent(mpd disk, int per) {
    if(disk == NULL || per < 1) {
        return false;
    }
    
    int prec = (int) mpd_prec(disk);
    long precR = prec <= 128 ? 61 - MPD_EXTRA_PREC : 113 - MPD_EXTRA_PREC;
    defs_mpd(prec, precR, bdC, bdV, bdD);
    
    mpd_set(bdV, disk);               // value of p_n
    mpd_set(bdC, disk);               // value of p_n
    mpd_setd(bdD, 1, 0, 0);           // derivative of p_n
    
    for (int i = 1; i < per ; i++) {
        // derivative d_n = 2 * v * d_{n-1} + 1
        mpd_mul(bdD, bdD, bdV);
        mpd_scale(bdD, bdD, 1);
        mpd_add_si(bdD, bdD, 1);
        
        if(i < per - 1) {
            // new value v_n = v_{n-1}^2 + c
            mpd_sqr(bdV, bdV);
            mpd_add(bdV, bdV, bdC);
        }
    }
    
    return ! mpd_contains_0(bdD, 0);
}

// MARK: counting functions

ulong mandel_hyp_count(int n) {
    if(n < 1 || n > 65) {
        return 0;
    }
    
    if(n < 3) {
        return 1;
    }
    
    ulong c = (1L << (n - 1)) + mu(n);
    for (int i = 2; i <= n / 2; i++) {
        if(n % i != 0)
            continue;
        
        c += mu(n / i) * (1L << (i - 1));
    }
    
    return c;
}

ulong mandel_mis_count(int pp, int per) {
    int n = pp + per;
    if(n < 1 || pp < 0 || pp >= n || n > 65) {
        return 0;
    }
    
    ulong h = mandel_hyp_count(per);
    if(pp == 0) {
        return h;
    }
    
    ulong d = (1L << pp) - (1L << (pp - 1)) - ((pp - 1) % per == 0 ? 1 : 0);
    
    return d * h;
}

ulong mandel_prim_hyp_count(int per) {
    if(per < 1 || per > 65) {
        return 0;
    }
    
    return mandel_hyp_count(per) - mandel_non_prim_hyp_count(per);
}

ulong mandel_non_prim_hyp_count(int per) {
    if(per < 2 || per > 129) {
        return 0;
    }
    
    ulong n = phi(per); // on the main cardioid
    int d, k;
    for (d = 2; d * d < per; d++) {
        if(per % d != 0) {
            continue;
        }
        
        k = per / d;
        n += mandel_hyp_count(d) * phi(k); // attached to components of period d
        n += mandel_hyp_count(k) * phi(d); // attached to components of period k
    }
    
    if(d * d == per) {
        n += mandel_hyp_count(d) * phi(d); // attached to components of period d
    }
    
    return n;
}

ulong mandel_mult_crit_count(int per) {
    ulong nu = 2 * mandel_hyp_count(per);
    
    ulong s = phi(per);
    int d, k;
    for (d = 2; d * d < per; d++) {
        if(per % d != 0) {
            continue;
        }
        
        k = per / d;
        s += mandel_hyp_count(d) * phi(k);
        s += mandel_hyp_count(k) * phi(d);
    }
    
    if(d * d == per) {
        s += mandel_hyp_count(d) * phi(d);
    }
    
    return nu - nu / per - s;
}
