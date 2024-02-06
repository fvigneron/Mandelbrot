//
//  dRect.c
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

#include "dRect.h"

void drect_init(drect r, long x, long y, uint w, uint h, int tpow) {
    r->x = x;
    r->y = y;
    r->w = w;
    r->h = h;
    r->tpow = tpow;
}

// //////////////////////////////////////////////////
/// MARK: Coordinates
// //////////////////////////////////////////////////

inline long drect_abs_to_rel_x80(drect r, ldbl x) {
    ldbl lx = ldexpl(x, r->tpow);
    
    return lx - r->x;
}

inline ldbl drect_rel_to_abs_x80(drect r, int x) {
    return ldexpl(r->x + x, -r->tpow);
}

inline long drect_abs_to_rel_y80(drect r, ldbl y) {
    ldbl ly = ldexpl(y, r->tpow);
    
    return ly - r->y;
}

inline ldbl drect_rel_to_abs_y80(drect r, int y) {
    return ldexpl(r->y + y, -r->tpow);
}

bool drect_rel_to_abs_coords80(fp80 c, drect r, long x, long y) {
    if(c == NULL || r == NULL) {
        return false;
    }
    
    c->x = r->x;
    c->x += x;
    c->x = ldexpl(c->x, -r->tpow);
    
    c->y = r->y;
    c->y += y;
    c->y = ldexpl(c->y, -r->tpow);
    
    return true;
}

bool drect_rel_to_abs_coords(mpc c, drect r, long x, long y) {
    if(c == NULL || r == NULL) {
        return false;
    }
    
    mpfr_set_si(c->x, r->x, MPFR_RNDN);
    mpfr_add_si(c->x, c->x, x, MPFR_RNDN);
    mpfr_mul_2si(c->x, c->x, -r->tpow, MPFR_RNDN);
    
    mpfr_set_si(c->y, r->y, MPFR_RNDN);
    mpfr_add_si(c->y, c->y, y, MPFR_RNDN);
    mpfr_mul_2si(c->y, c->y, -r->tpow, MPFR_RNDN);
    
    return true;
}

long drect_abs_to_rel_x(drect r, mpfr_t x) {//TODO: check implementation
    ldbl xx = mpfr_get_ld(x, MPFR_RNDN);
    return drect_abs_to_rel_x80(r, xx);
}

long drect_abs_to_rel_y(drect r, mpfr_t y) {//TODO: check implementation
    ldbl yy = mpfr_get_ld(y, MPFR_RNDN);
    return drect_abs_to_rel_y80(r, yy);
}

// //////////////////////////////////////////////////
/// MARK: Geometric operations on rectangles
// //////////////////////////////////////////////////

inline bool drect_union(drect uni, drect r) {
    if(uni == NULL || r == NULL || uni->tpow != r->tpow) {
        return false;
    }
    
    long l = uni->x <= r->x ? uni->x : r->x;
    long b = uni->y <= r->y ? uni->y : r->y;
    
    long ur = uni->x + uni->w;
    long rr = r->x + r->w;
    long ri = ur >= rr ? ur : rr;
    
    long ut = uni->y + uni->h;
    long rt = r->y + r->h;
    long t = ut >= rt ? ut : rt;
    
    uni->x = l;
    uni->y = b;
    uni->w = (uint) (ri - l);  // no checks for overflow, for speed; to be used reasonably
    uni->h = (uint) (t - b);   // no checks for overflow, for speed; to be used reasonably
    
    return true;
}

inline void drect_add(drect uni, long x, long y) {
    long l = uni->x <= x ? uni->x : x;
    long b = uni->y <= y ? uni->y : y;
    
    long ur = uni->x + uni->w;
    long r = ur >= x ? ur : x;
    
    long ut = uni->y + uni->h;
    long t = ut >= y ? ut : y;
    
    uni->x = l;
    uni->y = b;
    uni->w = (uint) (r - l);   // no checks for overflow, for speed; to be used reasonably
    uni->h = (uint) (t - b);   // no checks for overflow, for speed; to be used reasonably
}

void drect_rescale(drect r, int tpow) {
    int dp = tpow - r->tpow;
    
    if(dp == 0) {
        return;
    }
    
    if(dp > 0) {
        r->tpow = tpow;
        r->x <<= dp;
        r->y <<= dp;
        r->w <<= dp;
        r->h <<= dp;
        
        return;
    }
    
    dp = -dp; // does not mean dp == 0 :)
    
    r->tpow = tpow;
    long mx = r->x + r->w;
    r->x >>= dp;
    long my = r->y + r->h;
    r->y >>= dp;
    
    uint mask = (1 << dp) - 1;
    mx = (mx >> dp) + (mx & mask ? 1 : 0);
    my = (my >> dp) + (my & mask ? 1 : 0);
    
    r->w = (uint) (mx - r->x);
    r->h = (uint) (my - r->y);
}

inline void drect_translate(drect r, long dx, long dy, int tpow) {
    int tpd = r->tpow - tpow;
    long x = tpd >= 0 ? dx << tpd : dx >> -tpd;
    long y = tpd >= 0 ? dy << tpd : dy >> -tpd;
    
    r->x += x;
    r->y += y;
}

bool drect_center80(drect r, fp80 c) {
    if(r == NULL || c == NULL) {
        return false;
    }
    
    ldbl x = ldexpl(c->x, r->tpow) - r->w / 2;
    ldbl y = ldexpl(c->y, r->tpow) - r->h / 2;
    
    if(x > LONG_MAX || x < LONG_MIN || y > LONG_MAX || y < LONG_MIN) {
        return false;
    }
    
    r->x = roundl(x);
    r->y = roundl(y);
    
    return true;
}

bool drect_center(drect r, mpc delta, mpc c) {
    if(r == NULL || c == NULL || r->w < 0 || r->h < 0) {
        return false;
    }
    
    long ex = c->x->_mpfr_exp;
    long ey = c->y->_mpfr_exp;
    long ec = 1 + ex < ey ? ey : ex;
    
    if(delta == NULL) {
        if(r->tpow + ec > 62) {
            return false;
        }
        
        mpc_scale(c, c, r->tpow);
        r->x = mpfr_get_si(c->x, MPFR_RNDN) - r->w / 2;
        r->y = mpfr_get_si(c->y, MPFR_RNDN) - r->h / 2;
        mpc_scale(c, c, -r->tpow);
        
        return true;
    }
    
    if(r->tpow + ec > 62) {
        
        // translate the rectangle to (0, 0)
        r->x = 0;
        r->y = 0;
    } else {
        // here we can compute r->x and r->y, they are in range
        mpc_scale(delta, c, r->tpow + 1);
        mpfr_sub_ui(delta->x, delta->x, r->w, MPFR_RNDN);
        mpfr_sub_ui(delta->y, delta->y, r->h, MPFR_RNDN);
        mpc_scale(delta, delta, -1);
        
        r->x = mpfr_get_si(delta->x, MPFR_RNDN);
        r->y = mpfr_get_si(delta->y, MPFR_RNDN);
    }
    
    // compute the center of this rectangle
    mpfr_set_ui(delta->x, r->w, MPFR_RNDN);
    mpfr_set_ui(delta->y, r->h, MPFR_RNDN);
    mpc_scale(delta, delta, -1);
    mpfr_add_si(delta->x, delta->x, r->x, MPFR_RNDN);
    mpfr_add_si(delta->y, delta->y, r->y, MPFR_RNDN);
    mpc_scale(delta, delta, -r->tpow);
    
    // compute the translation delta
    mpc_sub(delta, c, delta);
    
    return true;
}

bool drect_intersection(drect inter, drect r1, drect r2) {
    if(inter == NULL || ! drect_intersect(r1, r2)) {
        return false;
    }
    
    inter->tpow = r1->tpow;
    
    inter->x = r1->x < r2->x ? r2->x : r1->x;
    long dx1 = r1->x + r1->w;
    long dx2 = r2->x + r2->w;
    long dx = dx1 < dx2 ? dx1 : dx2;
    inter->w = (uint) (dx - inter->x);
    
    inter->y = r1->y < r2->y ? r2->y : r1->y;
    long dy1 = r1->y + r1->h;
    long dy2 = r2->y + r2->h;
    long dy = dy1 < dy2 ? dy1 : dy2;
    inter->h = (uint) (dy - inter->y);
    
    return true;
}

// //////////////////////////////////////////////////
/// MARK: Inclusion, intersection
// //////////////////////////////////////////////////

inline bool drect_contains80(drect r, fp80 c) {
    ldbl lx = ldexpl(c->x, r->tpow) - r->x;
    ldbl ly = ldexpl(c->y, r->tpow) - r->y;
    
    return lx >= 0 && ly >= 0 && lx < r->w && ly < r->h;
}

bool drect_contains(drect r, mpc c) {
    fp80 cc;
    cc->x = mpfr_get_ld(c->x, MPFR_RNDD);
    cc->y = mpfr_get_ld(c->y, MPFR_RNDD);
    
    return drect_contains80(r, cc);
}

bool drect_contains_rect(drect r1, drect r2) {
    if(r1 == NULL || r2 == NULL) {
        return false;
    }
    
    if(r1->tpow == r2->tpow) {
        return r2->x >= r1->x && r2->y >= r1->y && r2->x + r2->w <= r1->x + r1->w &&
            r2->y + r2->h <= r1->y + r1->h;
    }
    
    if(r1->tpow > r2->tpow) {
        uint dt = r1->tpow - r2->tpow;
        
        return (r2->x << dt) >= r1->x && (r2->y << dt) >= r1->y &&
            ((r2->x + r2->w) << dt) <= r1->x + r1->w &&
            ((r2->y + r2->h) << dt) <= r1->y + r1->h;
    }
    
    uint dt = r2->tpow - r1->tpow;
    
    return r2->x >= (r1->x << dt) && r2->y >= (r1->y << dt) &&
        r2->x + r2->w <= ((r1->x + r1->w) << dt) &&
        r2->y + r2->h <= ((r1->y + r1->h) << dt);
}

inline bool drect_intersect(drect r1, drect r2) { 
    if(r1 == NULL || r2 == NULL) {
        return false;
    }
    
    if(r1->tpow == r2->tpow) {
        bool outx = r1->x >= r2->x + r2->w || r2->x >= r1->x + r1->w;
        bool outy = r1->y >= r2->y + r2->h || r2->y >= r1->y + r1->h;
        
        return ! outx && ! outy;
    }
    
    if(r1->tpow > r2->tpow) {
        uint dt = r1->tpow - r2->tpow;
        
        bool outx = r1->x >= ((r2->x + r2->w) << dt) || (r2->x << dt) >= r1->x + r1->w;
        bool outy = r1->y >= ((r2->y + r2->h) << dt) || (r2->y << dt) >= r1->y + r1->h;
        
        return ! outx && ! outy;
    }
    
    uint dt = r2->tpow - r1->tpow;
    
    bool outx = (r1->x << dt) >= r2->x + r2->w || r2->x >= ((r1->x + r1->w) << dt);
    bool outy = (r1->y << dt) >= r2->y + r2->h || r2->y >= ((r1->y + r1->h) << dt);
    
    return ! outx && ! outy;
}

bool drect_intersects_disk(drect dr, fp80 c, ldbl r) {
    long sc = 1L << dr->tpow;
    ldbl cx = sc * c->x, cy = sc * c->y, lr = sc * r;

    long xr = dr->x + dr->w, yt = dr->y + dr->h;
    if(dr->x > cx + lr || dr->y > cy + lr || xr < cx - lr || yt < cy - lr) {
        return false;
    }
    
    bool left = cx < dr->x, right = ! left && cx > xr;
    bool bott = cy < dr->y, top = ! bott && cy > yt;
    
    if(! left && ! right) {
        ldbl d = top ? cy - yt : bott ? dr->y - cy : 0;
        
        return d <= lr;
    }
    
    if(! bott && ! top) {
        ldbl d = right ? cx - xr : dr->x - cx;
        
        return d <= lr;
    }
    
    ldbl dx = cx - (left ? dr->x : xr);
    ldbl dy = cy - (bott ? dr->y : yt);
        
    return dx * dx + dy * dy <= lr * lr;
}

bool inline drect_equals(drect r1, drect r2) {
    return r1 != NULL && r2 != NULL && r1->x == r2->x && r1->y == r2->y &&
           r1->w == r2->w && r1->h == r2->h && r1->tpow == r2->tpow;
}

// //////////////////////////////////////////////////
/// MARK: Pretty print
// //////////////////////////////////////////////////

int drect_print(drect r, char *str, int maxLen) {
    return snprintf(str, maxLen, "(%ld, %ld, %d, %d) >> %d", r->x, r->y, r->w, r->h, r->tpow);
}
