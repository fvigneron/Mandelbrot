//
//  planarSet.c
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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "fp80.h"
#include "planarSet.h"


bool pset_init(pset ps, double eps) {
    ps->eps = eps;
    ps->absoluteEps = true;
    ps->count = 0;
    ps->barCount = 1;
    ps->lastCount = 0;
    ps->realCount = 0;
    ps->locked = 0;
    ps->rejected = 0;
        
    ps->pts[0] = malloc(PSET_DEF_SIZE * sizeof(fp80_struct));
    ps->barLen[0] = PSET_DEF_SIZE;
    
    for(int i = 1; i < PSET_MAX_BARS; i ++) {
        ps->pts[i] = NULL;
        ps->barLen[i] = 0;
    }
    
    return ps->pts[0] != NULL;
}

void pset_clear(pset ps) {
    for(int i = 0; i < ps->barCount; i ++) {
        if(ps->pts[i] != NULL && ps->barLen[i] > 0) {
            free(ps->pts[i]);
        }
        
        ps->pts[i] = NULL;
        ps->barLen[i] = 0;
    }
    
    ps->count = 0;
    ps->barCount = 0;
    ps->lastCount = 0;
    ps->realCount = 0;
    ps->locked = 0;
    ps->rejected = 0;
    ps->eps = 0;
}

bool pset_valid(pset ps) {
    if(ps == NULL || ps->eps < 0 || ps->eps >= 1E-3 || ps->count < 0 || ps->barCount < 1
       || ps->barCount > PSET_MAX_BARS) {
        return false;
    }
    
    if(ps->lastCount < 0 || ps->realCount < 0 || ps->realCount > ps->count || ps->rejected < 0) {
        return false;
    }
    
    if((ps->locked && ps->barCount != 1) || (! ps->locked && ps->lastCount > PSET_DEF_SIZE)) {
        return false;
    }
    
    for(long i = 0; i < ps->barCount; i ++) {
        if(ps->barLen[i] <= 0 || ps->pts[i] == NULL) {
            return false;
        }
    }
    
    return true;
}

bool pset_point(fp80 p, pset ps, long pInd) {
    if(pInd < 0 || pInd >= ps->count) {
        return false;
    }
    
    long pos = pInd;
    int bar = 0;
    
    while(pos >= ps->barLen[bar] && bar < ps->count) {
        pos -= ps->barLen[bar ++];
    }
    
    if(pos < 0 || bar >= ps->barCount || (bar == ps->barCount - 1 && pos >= ps->lastCount)) {
        return false;
    }
    
    *p = ps->pts[bar][pos];
    
    return true;
}

static inline bool pset_sless(fp80 a, fp80 b) {
    return a->x < b->x || (a->x == b->x && a->y < b->y);
}

static inline bool pset_smore(fp80 a, fp80 b) {
    return a->x > b->x || (a->x == b->x && a->y > b->y);
}

static inline bool pset_eq(fp80 a, fp80 b, double eps) {
    long double dx = b->x - a->x;
    dx = dx < 0 ? -dx : dx;
    
    long double dy = b->y - a->y;
    dy = dy < 0 ? -dy : dy;
    
    return dx <= eps && dy <= eps;
}

static fp80_struct* pset_merge(fp80_struct u[], long uLen, fp80_struct v[], long vLen) {
    long size = (uLen + vLen) * sizeof(fp80_struct);
    fp80_struct* d = malloc(size);
    
    long i = 0, j = 0, dp = 0;
    while(i < uLen || j < vLen) {
        if(i < uLen && (j == vLen || pset_sless(&u[i], &v[j]))) {
            d[dp ++] = u[i ++];
        } else {
            d[dp ++] = v[j ++];
        }
    }
    
    return d;
}

static long pset_pos(fp80 p, fp80_struct* list, long len) {
    long fr = 0, ls = len - 1;
    if(pset_sless(p, list + fr))
        return 0;
    
    if(pset_smore(p, list + ls))
        return len;
    
    long mid;
    do {
        mid = (fr + ls) / 2;
        
        if(pset_sless(p, list + mid)) {
            ls = mid;
            
            continue;
        }
        
        if(pset_smore(p, list + mid)) {
            fr = mid;
            
            continue;
        }
        
        return mid;
    } while(ls - fr > 1);
    
    return pset_sless(p, list + ls) ? fr + 1: ls + 1;
}

static long pset_search(pset ps, long bar, fp80 p) {
    fp80_struct* list = ps->pts[bar];
    long len = bar == ps->barCount - 1 ? ps->lastCount : ps->barLen[bar];
    if(len <= 0)
        return -1;
    
    ldbl eps = ps->eps * (ps->absoluteEps ? 1 : fp80_mod(p));
    fp80 c = {p->x - eps, p->y - eps};
    long l = pset_pos(c, list, len);
    
    c->x += 2 * eps;
    c->y += 2 * eps;
    long r = pset_pos(c, list, len);
    
    c->x = p->x;
    c->y = p->y;
    
    long right = 1;
    long pos = l;
    for (long i = l; i < r; i ++) {
        if(pset_eq(c, list + i, eps))
            return i;
        
        right = right && pset_smore(c, list + i);
        pos += right ? 1 : 0;
    }
    
    return -pos - 1;
}

long pset_index(pset ps, fp80 p) {
    if(ps == NULL || p == NULL || ! ps->locked) {
        return -1;
    }
    
    long pos = pset_search(ps, 0, p);
    
    return pos >= 0 ? pos : -1;
}

bool pset_contains(pset ps, fp80 p) {
    for (int i = ps->barCount - 1; i >= 0; i --) {
        if(pset_search(ps, i, p) >= 0) {
            return true;
        }
    }
    
    return false;
}

static void pset_mergeLast(pset ps) {
    int l = ps->barCount - 1;
    long len1 = ps->barLen[l - 1];
    long len2 = ps->lastCount;
    
    fp80_struct *list = ps->pts[l - 1];
    ps->pts[l - 1] = pset_merge(list, len1, ps->pts[l], len2);
    free(list);
    free(ps->pts[l]);
    ps->pts[l] = NULL;
    
    ps->barLen[l - 1] = len1 + len2;
    ps->barLen[l] = 0;
    
    ps->lastCount = len1 + len2;
    ps->barCount --;
}

void pset_lock(pset ps) {
    while(ps->barCount > 1) {
        pset_mergeLast(ps);
    }
    
    ps->locked = true;
}

static bool pset_pack(pset ps) {
    while(ps->barCount > 1 && ps->lastCount == ps->barLen[ps->barCount - 2]) {
        pset_mergeLast(ps);
    }
    
    if(ps->barCount == PSET_MAX_BARS) {
        return false;
    }
    
    // grow the points list
    ps->barLen[ps->barCount] = PSET_DEF_SIZE;
    ps->pts[ps->barCount ++] = malloc(PSET_DEF_SIZE * sizeof(fp80_struct));
    ps->lastCount = 0;
    
    return true;
}

static bool pset_insert(pset ps, fp80 p) {    
    long pos = pset_search(ps, ps->barCount - 1, p);
    if(pos >= 0) {
        return false;
    }
    
    pos = -pos - 1;
    long len = ps->lastCount;
    fp80_struct* v = ps->pts[ps->barCount - 1];
    for (long i = len; i > pos; i--) {
        v[i] = v[i - 1];
    }
    v[pos] = *p;
    
    ps->lastCount ++;
    ps->count ++;
    
    ldbl py = p->y;
    double eps = ps->eps * (ps->absoluteEps ? 1 : fp80_mod(p)) / 2;
    ps->realCount += py > -eps && py < eps ? 1 : 0;
    
    return true;
}

bool pset_add(pset ps, fp80 p) {
    if(ps->locked || pset_contains(ps, p)) {
        ps->rejected ++;

        return false;
    }
    
    if(ps->lastCount == PSET_DEF_SIZE) {
        if(! pset_pack(ps)) {
            ps->rejected ++;

            return false;
        }
    }
    
    return pset_insert(ps, p);
}

bool pset_read(pset ps, char *fileName) {
    if(! pset_valid(ps)) {
        return false;
    }
    
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return false;
    }
    
    fp80 p;
    char sx[50], sy[50];
    while(! feof(f)) {
        if(fscanf(f, "%s%s", sx, sy) == 2) {
            long lx = strlen(sx);
            if(lx > 1 && sx[lx - 1] == ',') {
                sx[lx - 1] = 0;
            }
            
            int ok = sscanf(sx, "%Lf", &p->x) == 1;
            ok = ok && sscanf(sy, "%Lf", &p->y) == 1;
            
            if(ok) {
                pset_add(ps, p);
            }
        }
    }
    
    fclose(f);
    
    return true;
}

bool pset_write(pset ps, char *fileName, bool append) {
    if(! pset_valid(ps)) {
        printf("The pSet is not valid !\n");
        
        return false;
    }
    
    FILE *f = fopen(fileName, append ? "a" : "w");
    if(f == NULL) {
        printf("Could not open %s for writing!\n", fileName);
        
        return false;
    }
    
    bool ok = true;
    fp80 p;
    for (long i = 0; i < ps->count && ok; i ++) {
        ok = pset_point(p, ps, i);
        ok = ok && fprintf(f, "%.21Lg, %.21Lg\n", p->x, p->y) >= 4;
    }
    
    fclose(f);
    
    return ok;
}
