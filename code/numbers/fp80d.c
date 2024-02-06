//
//  fp80d.c
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

#include <math.h>

#include "fp80d.h"

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Constructors of disks
// //////////////////////////////////////////////////////////////////////////////////////////

void fp80d_set(fp80d d, fp80d s) {
    d->x = s->x;
    d->y = s->y;
    d->r = s->r;
    d->m = s->m;
}

inline void fp80d_setl(fp80d d, ldbl x, ldbl y, ldbl r) {
    d->x = x;
    d->y = y;
    d->r = r;
    
    d->m = hypotl(x, y);
}

inline void fp80d_set80(fp80d d, fp80 c, ldbl r) {
    d->x = c->x;
    d->y = c->y;
    d->r = r;
    
    d->m = hypotl(d->x, d->y);
}

void fp80d_setUlp(fp80d d, ldbl x, ldbl y) { // ULP = unit on last position
    d->x = x;
    d->y = y;
    
    d->m = hypotl(x, y);
    
    d->r = nextafterl(d->m, d->m + 1) - d->m;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Disk operations
// //////////////////////////////////////////////////////////////////////////////////////////

void fp80d_add(fp80d d, fp80d op1, fp80d op2) { // sum of disks (add centers, add radii)
    d->x = op1->x + op2->x;
    d->y = op1->y + op2->y;
    d->r = op1->r + op2->r;
    
    d->m = hypotl(d->x, d->y);
}

void fp80d_addd(fp80d d, fp80d op1, ldbl a) { // translate along x axis
    d->x = op1->x + a;
    d->y = op1->y;
    d->r = op1->r;
    
    d->m = hypotl(d->x, d->y);
}

void fp80d_sub(fp80d d, fp80d op1, fp80d op2) { // difference of disks (sub centers, add radii)
    d->x = op1->x - op2->x;
    d->y = op1->y - op2->y;
    d->r = op1->r + op2->r;
    
    d->m = hypotl(d->x, d->y);
}

void fp80d_mul(fp80d d, fp80d op1, fp80d op2) { // multiplication of disks
    ldbl x1 = op1->x;
    ldbl x2 = op2->x;
    ldbl y1 = op1->y;
    ldbl y2 = op2->y;
    ldbl m1 = op1->m;
    ldbl m2 = op2->m;
    ldbl r2 = op2->r;
    
    d->x = x1 * x2 - y1 * y2;   // complex multiplication of centers z1 * z2
    d->y = x2 * y1 + x1 * y2;
    d->r = op1->r * (m2 + r2) + m1 * r2; // radius because (|z1|+r1)*(|z2|+r2) - |z1*z2| = r1*(|z2|+r2) + |z1|*r2
    
    d->m = m1 * m2;     // |z1*z2| = |z1| * |z2|
}

bool fp80d_inv(fp80d d, fp80d op) {
    if(op->m <= op->r) {
        d->x = 0;
        d->y = 0;
        d->r = INFINITY;
        d->m = 0;
        
        return false;
    }
    
    ldbl m2 = op->m * op->m;
    d->x = op->x / m2;
    d->y = -op->y / m2;
    d->r = op->r / (m2 - op->r * op->m);
    d->m = 1 / op->m;
    
    return true;
}

bool fp80d_div(fp80d d, fp80d op1, fp80d op2) {
    if(op2->m <= op2->r) {
        d->x = 0;
        d->y = 0;
        d->r = INFINITY;
        d->m = 0;
        
        return false;
    }
    
    ldbl m2 = op2->m * op2->m;
    ldbl dx = op1->x * op2->x + op1->y * op2->y;
    ldbl dy = op1->y * op2->x - op1->x * op2->y;
    
    d->x = dx / m2;
    d->y = dy / m2;
    d->r = (op1->m * op2->r + op2->m * op1->r) / (m2 - op2->r * op2->m);
    d->m = op1->m / op2->m;
    
    return true;
}

void fp80d_muld(fp80d d, fp80d op1, fp80d op2, ldbl a) { // multiplication then scaling by a factor
    ldbl x1 = op1->x;
    ldbl x2 = op2->x;
    ldbl y1 = op1->y;
    ldbl y2 = op2->y;
    ldbl m1 = op1->m;
    ldbl m2 = op2->m;
    ldbl r2 = op2->r;
    ldbl aa = a < 0 ? -a : a;    // aa = |a|
    
    d->x = a * (x1 * x2 - y1 * y2);
    d->y = a * (x2 * y1 + x1 * y2);
    d->r = aa * (op1->r * (m2 + r2) + m1 * r2);
    
    d->m = aa * (m1 * m2);
}

void fp80d_sqr(fp80d d, fp80d op1) { // identical but faster than fp80Disk_mul(d, op1, op1)
    ldbl x = op1->x;
    d->x = x * x - op1->y * op1->y;
    d->y = 2 * x * op1->y;
    d->r = op1->r * (op1->r + 2 * op1->m);
    
    d->m = op1->m * op1->m;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Distance related
// //////////////////////////////////////////////////////////////////////////////////////////

ldbl fp80d_min_mod(fp80d op1) { // minimal modulus of points in the disk
    ldbl m = op1->m - op1->r;
    
    return m < 0 ? 0 : m;
}

ldbl fp80d_max_mod(fp80d op1) { // maximal modulus of points in the disk
     return op1->m + op1->r;
}

bool fp80d_contains(fp80d op1, fp80d op2) { // check if op1 contains the disk op2 (both open or both closed)
    ldbl dr = op1->r - op2->r;
    
    return dr >= 0 && dr >= hypotl(op2->x - op1->x, op2->y - op1->y); // dist between centers + small radius <= big radius
}

// open disks
bool fp80d_intersect(fp80d op1, fp80d op2) {
    ldbl dx = op1->x - op2->x;
    ldbl dy = op1->y - op2->y;
    
    ldbl d2 = dx * dx + dy * dy;
    ldbl sr = op1->r + op2->r;
    
    return d2 < sr * sr;
}
