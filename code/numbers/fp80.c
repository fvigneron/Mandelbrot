//
//  fp80.c
//
//  Created by MIHALACHE Nicolae on 11/18/19.
//  Revised by VIGNERON François on 01/15/21.
//

#include <stdio.h>

#include "fp80.h"

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Conversion
// //////////////////////////////////////////////////////////////////////////////////////////

inline void fp80_setl(fp80 d, fp80 s) {
    d->x = s->x;
    d->y = s->y;
}

inline void fp80_set(fp80 d, ldbl x, ldbl y) {
    d->x = x;
    d->y = y;
}

inline void fp80_setd(fp80 d, fp64 s) {
    d->x = s->x;
    d->y = s->y;
}

inline void fp64_setl(fp64 d, fp80 s) {
    d->x = s->x;
    d->y = s->y;
}

inline void fp64_setd(fp64 d, fp64 s) {
    d->x = s->x;
    d->y = s->y;
}

inline void fp64_set(fp80 d, double x, double y) {
    d->x = x;
    d->y = y;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Additions
// //////////////////////////////////////////////////////////////////////////////////////////

inline void fp80_add(fp80 v, fp80 x, fp80 y) {
    v->x = x->x + y->x;
    v->y = x->y + y->y;
}

inline void fp80_sub(fp80 v, fp80 x, fp80 y) {
    v->x = x->x - y->x;
    v->y = x->y - y->y;
}

inline void fp64_add(fp64 v, fp64 x, fp64 y) {
    v->x = x->x + y->x;
    v->y = x->y + y->y;
}

inline void fp64_sub(fp64 v, fp64 x, fp64 y) {
    v->x = x->x - y->x;
    v->y = x->y - y->y;
}

inline void fp80_neg(fp80 v, fp80 x) {    
    v->x = -x->x;
    v->y = -x->y;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Multiplications
// //////////////////////////////////////////////////////////////////////////////////////////

inline void fp80_mul(fp80 v, fp80 x, fp80 y) {
    ldbl px = x->x * y->x - x->y * y->y;
    ldbl py = x->x * y->y + x->y * y->x;
    
    v->x = px;
    v->y = py;
}

inline void fp64_mul(fp64 v, fp64 x, fp64 y) {
    double px = x->x * y->x - x->y * y->y;
    double py = x->x * y->y + x->y * y->x;
    
    v->x = px;
    v->y = py;
}

inline void fp80_muld(fp80 v, fp80 x, double y) {
    v->x = x->x * y;
    v->y = x->y * y;
}

inline void fp64_muld(fp64 v, fp64 x, double y) {
    v->x = x->x * y;
    v->y = x->y * y;
}

inline void fp80_scale(fp80 v, fp80 x, int tp) {
    if(tp >= 0 && tp < 64) {
        ulong m = 1L << tp;
        v->x = x->x * m;
        v->y = x->y * m;
        
        return;
    }
    
    if(tp < 0 && tp > -64) {
        ulong m = 1L << (-tp);
        v->x = x->x / m;
        v->y = x->y / m;
        
        return;
    }
    
    v->x = ldexpl(x->x, tp);
    v->y = ldexpl(x->y, tp);
}

inline void fp80_mull(fp80 v, fp80 x, ldbl y) {
    v->x = x->x * y;
    v->y = x->y * y;
}

inline void fp80_muli(fp80 v, fp80 x, long y) {
    v->x = x->x * y;
    v->y = x->y * y;
}

inline void fp80_sqr(fp80 v, fp80 x) {
    ldbl px = x->x * x->x - x->y * x->y;
    ldbl py = 2 * x->x * x->y;
    
    v->x = px;
    v->y = py;
}

inline void fp64_sqr(fp64 v, fp64 x) {
    double px = x->x * x->x - x->y * x->y;
    double py = 2 * x->x * x->y;
    
    v->x = px;
    v->y = py;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Divisions
// //////////////////////////////////////////////////////////////////////////////////////////

inline bool fp80_inv(fp80 v, fp80 z) {
    ldbl m2 = z->x * z->x + z->y * z->y;
    if(m2 == 0 || ! isnormal(m2)) {
        return false;
    }
    
    v->x = z->x / m2;
    v->y = -z->y / m2;
    
    return true;
}

inline bool fp80_div(fp80 v, fp80 x, fp80 y) {
    ldbl my2 = y->x * y->x + y->y * y->y;
    if(my2 == 0 || ! isnormal(my2)) {
        return false;
    }
    
    ldbl px = x->x * y->x + x->y * y->y;
    ldbl py = x->y * y->x - x->x * y->y;
    
    v->x = px / my2;
    v->y = py / my2;
    
    return true;
}

inline bool fp64_div(fp64 v, fp64 x, fp64 y) {
    double my2 = y->x * y->x + y->y * y->y;
    if(my2 == 0 || ! isnormal(my2)) {
        return false;
    }
    
    ldbl px = x->x * y->x + x->y * y->y;
    ldbl py = x->y * y->x - x->x * y->y;
    
    v->x = px / my2;
    v->y = py / my2;
    
    return true;
}

inline void fp80_quick_div(fp80 v, fp80 x, fp80 y) {
    ldbl my2 = y->x * y->x + y->y * y->y;
        
    v->x = (x->x * y->x + x->y * y->y) / my2;
    v->y = (x->y * y->x - x->x * y->y) / my2;
}

inline void fp64_quick_div(fp64 v, fp64 x, fp64 y) {
    double my2 = y->x * y->x + y->y * y->y;
        
    v->x = (x->x * y->x + x->y * y->y) / my2;
    v->y = (x->y * y->x - x->x * y->y) / my2;
}

inline bool fp80_divd(fp80 v, fp80 x, double y) {
    if(y == 0) {
        return false;
    }
    
    v->x = x->x / y;
    v->y = x->y / y;
    
    return true;
}

inline bool fp80_divl(fp80 v, fp80 x, ldbl y) {
    if(y == 0) {
        return false;
    }
    
    v->x = x->x / y;
    v->y = x->y / y;
    
    return true;
}

inline bool fp80_divi(fp80 v, fp80 x, long y) {
    if(y == 0) {
        return false;
    }
    
    v->x = x->x / y;
    v->y = x->y / y;
    
    return true;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Distance related
// //////////////////////////////////////////////////////////////////////////////////////////

inline ldbl fp80_mod(fp80 p) {
    return sqrtl(p->x * p->x + p->y * p->y);
}

inline ldbl fp80_dist(fp80 a, fp80 b) {
    ldbl dx = a->x - b->x;
    ldbl dy = a->y - b->y;
    
    return sqrtl(dx * dx + dy * dy);
}

inline ldbl fp80_mod2(fp80 p) {
    return p->x * p->x + p->y * p->y;
}

inline double fp64_mod2(fp64 p) {
    return p->x * p->x + p->y * p->y;
}

inline ldbl fp80_dist2(fp80 a, fp80 b) {
    ldbl dx = a->x - b->x;
    ldbl dy = a->y - b->y;
    
    return dx * dx + dy * dy;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Trigonometric
// //////////////////////////////////////////////////////////////////////////////////////////

bool fp80_exp_2Pi_i(fp80 v, ldbl theta) {
    ldbl ip;
    ldbl th = 2 * PI * modfl(theta, &ip);
    
    v->x = cosl(th);
    v->y = sinl(th);
    
    return true;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Other
// //////////////////////////////////////////////////////////////////////////////////////////

inline void fp80_sqrt(fp80 v, fp80 z) {
    ldbl m = fp80_mod(z);
    ldbl x = z->x;
    ldbl mx = x < 0 ? -x : x;
    
    // rarely, there is a large loss of precision, when |y| << |x|
    if(m - mx < 1E-6 * mx) {
        ldbl ha = z->y / x;
        ha *= ha;
        ha /= 2;
        ldbl qa2 = ha * ha;
        
        if(x > 0) {
            v->x = sqrtl((m + x) / 2);
            
            ldbl mmx = x * (ha + qa2 * (ha - 1) / 2);
            ldbl sqmmx = sqrtl(mmx / 2);
            v->y = z->y >= 0 ? sqmmx : -sqmmx;
        } else {
            ldbl mpx = mx * (ha + qa2 * (ha - 1) / 2);
            v->x = sqrtl(mpx / 2);            
            
            ldbl sqmmx = sqrtl((m - x) / 2);
            v->y = z->y >= 0 ? sqmmx : -sqmmx;
        }
        
        return;
    }
    
    v->x = sqrtl((m + x) / 2);
    v->y = z->y >= 0 ? sqrtl((m - x) / 2) : -sqrtl((m - x) / 2);
}


inline void fp80_min(fp80 min, fp80 a, fp80 b) {
    min->x = a->x <= b->x ? a->x : b->x;
    min->y = a->y <= b->y ? a->y : b->y;
}

inline void fp80_max(fp80 max, fp80 a, fp80 b) {
    max->x = a->x >= b->x ? a->x : b->x;
    max->y = a->y >= b->y ? a->y : b->y;
}

inline void fp80_min_max(fp80 min, fp80 max, fp80 a) {
    fp80_min(min, min, a);
    fp80_max(max, max, a);
}

inline bool fp80_is_number(fp80 z) {
    return is_number(z->x) && is_number(z->y);
}

inline bool fp64_is_number(fp64 z) {
    return is_number((ldbl) z->x) && is_number((ldbl) z->y);
}

inline ldbl fp80_det(fp80 a, fp80 b) {
    return a->x * b->y - a->y * b->x;
}

inline bool fp80_is_exact(fp80 c, ldbl err) {
    if(c == NULL || err <= 0 || isnan(err)) {
        return false;
    }
    
    return err <= fp80_lulp(c);
}

void fp80_print(fp80 z) {
    printf("(%.20Lg, %.20Lg)", z->x, z->y);
}
