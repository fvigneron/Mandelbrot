//
//  nSet.c
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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include <math.h>

#include "stopWatch.h"
#include "io.h"
#include "nSet.h"
#include "memFile.h"


// MARK: initialization

void nset_init(nset ps, unsigned long eps) {
    ps->eps = eps;
    ps->count = 0;
    ps->barCount = 1;
    ps->lastCount = 0;
    ps->realCount = 0;
    ps->locked = 0;
    ps->rejected = 0;
        
    ps->pts[0] = malloc(NSET_DEF_SIZE * sizeof(u128c_struct));
    ps->barLen[0] = NSET_DEF_SIZE;
    
    for(int i = 1; i < NSET_MAX_BARS; i ++) {
        ps->pts[i] = NULL;
        ps->barLen[i] = 0;
    }
    
    ps->maxMem = sizeof(nSet_struct) + NSET_DEF_SIZE * sizeof(u128c_struct);
}

nset nset_new(unsigned long eps, bool locked) {
    nset set;
    set = malloc(sizeof(nSet_struct));
    nset_init(set, eps);
    
    if(locked) {
        nset_lock(set);
    }
    
    return set;
}

bool nset_clear(nset ps) {
    if(ps == NULL) {
        return false;
    }
    
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
    ps->locked = true;
    ps->rejected = 0;
    
    ps->maxMem = sizeof(nSet_struct);
    
    return true;
}

void nset_clears(nset list[], int count) {
    for (int i = 0; i < count; i++) {
        int skip = 0;
        for (int j = 0; j < i && ! skip; j++) {
            skip = list[j] == list[i];
        }
        
        if(skip) {
            continue;
        }
        
        nset_free(list[i]);
    }
}

void nset_free(nset ps) {
    if(ps == NULL) {
        return;
    }
    
    nset_clear(ps);
    free(ps);
}

bool nset_valid(nset ps) {
    if(ps == NULL || ps->count < 0 ||  ((! ps->locked || ps->count > 0) && ps->barCount < 1)
       || ps->barCount > NSET_MAX_BARS) {
        return false;
    }
    
    if(ps->lastCount < 0 || ps->realCount < 0 || ps->realCount > ps->count || ps->rejected < 0) {
        return false;
    }
    
    if((ps->locked && ps->barCount > 1) || (! ps->locked && ps->lastCount > NSET_DEF_SIZE)) {
        return false;
    }
    
    for(int i = 0; i < ps->barCount; i ++) {
        if(ps->barLen[i] <= 0 || ps->pts[i] == NULL) {
            return false;
        }
    }
    
    return true;
}

static void nset_mergeLast(nset ps);

bool nset_lock(nset ps) {
    if(ps == NULL) {
        return false;
    }
    
    while(ps->barCount > 1) {
        nset_mergeLast(ps);
    }
    
    ps->locked = true;
    
    return true;
}

static bool nset_pack(nset ps);

bool nset_unlock(nset ps) {
    if(ps == NULL) {
        return false;
    }
    
    if(ps->barCount > 1 || ! ps->locked) {
        return false;
    }
    
    ps->locked = false;
    
    return nset_pack(ps);
}

void nset_move(nset dst, nset src, bool lock) {
    if(dst == src || dst == NULL || src == NULL) {
        return;
    }
    
    nset_clear(dst);
    
    *dst = *src;
    
    nset_init(src, src->eps);
    
    if(lock) {
        nset_lock(dst);
        nset_lock(src);
    }
}

// MARK: static functions

/// @brief Merges two striclly increasing and disjoint internal lists into an increasing list, without checking for points
/// that are equal up to some distance.
///
/// @param u the first list
/// @param uLen the lenght of the first list
/// @param v the second list
/// @param vLen the second of the first list
///
/// @return the new increasing list, of lenght @c uLen+vLen
static u128c_struct* nset_merge(u128c_struct u[], long uLen, u128c_struct v[], long vLen) {
    long size = (uLen + vLen) * sizeof(u128c_struct);
    u128c_struct* d = malloc(size);
    
    long i = 0, j = 0, dp = 0;
    while(i < uLen || j < vLen) {
        if(i < uLen && (j == vLen || u128_sless(&u[i], &v[j]))) {
            d[dp ++] = u[i ++];
        } else {
            d[dp ++] = v[j ++];
        }
    }
    
    return d;
}

/// @brief Returns the position @c pos of @c p in @c list, defined as follows: if @c p>list[len-1] then
/// @c pos:=len, otherwise @c 0<=pos<len is the smallest index such that @c p<=list[pos].
///
/// In other words, @c 0<=pos<=len and it is minial such that @c p<=list[pos].
///
/// @param p the point
/// @param list the list of points
/// @param len the length of the list
///
/// @return the position of @c p in @c list
static ulong nset_pos(u128 p, u128c_struct* list, long len) {
    long fr = 0, ls = len - 1;
    if(u128_leq(p, list + fr)) {
        return 0;
    }
    
    if(u128_smore(p, list + ls)) {
        return len;
    }
    
    if(ls == 1) {
        return 1;
    }
    
    long mid;
    do {
        mid = (fr + ls) / 2;
        
        if(u128_sless(p, list + mid)) {
            ls = mid;
            
            continue;
        }
        
        if(u128_smore(p, list + mid)) {
            fr = mid;
            
            continue;
        }
        
        return mid;
    } while(ls - fr > 1);
    
    return ls;
}

static long nset_bar_search(nset ps, int bar, u128 p) {
    u128c_struct* list = ps->pts[bar];
    long len = bar == ps->barCount - 1 ? ps->lastCount : ps->barLen[bar];
    if(len <= 0) {
        return -1;
    }
    
    ulong eps = ps->eps;
    u128 c = {uint128_sub(p->x, eps), uint128_sub(p->y, eps)};
    long l = nset_pos(c, list, len);
    
    c->x = p->x + eps;
    c->y = p->y + eps;
    long r = nset_pos(c, list, len);
    
    long right = 1;
    long pos = l;
    for (long i = l; i <= r && i < len; i ++) {
        if(u128_eq(p, list + i, eps)) {
            return i;
        }
        
        right = right && u128_smore(p, list + i);
        pos += right ? 1 : 0;
    }
    
    return -pos - 1;
}

long nset_search(nset ps, u128 p) {
    if(ps == NULL || p == NULL || ! ps->locked) {
        return LONG_MIN;
    }
    
    return nset_bar_search(ps, 0, p);
}

/// @brief Same as nset_serach() with the difference that if found, the coordinates of @c p are updated in @ ps.
///
/// @param ps the set
/// @param bar the index of the bar of @c pss
/// @param p the point to search for
///
/// @return the position of @c p, positive is it is found in the @c bar, otherwise negative, as explained above
static long nset_searchAndReplace(nset ps, int bar, u128 p) {
    u128c_struct* list = ps->pts[bar];
    long len = bar == ps->barCount - 1 ? ps->lastCount : ps->barLen[bar];
    if(len <= 0) {
        return -1;
    }
    
    ulong eps = ps->eps;
    u128 c = {uint128_sub(p->x, eps), uint128_sub(p->y, eps)};
    long l = nset_pos(c, list, len);
    
    c->x = p->x + eps;
    c->y = p->y + eps;
    long r = nset_pos(c, list, len);
    
    long right = 1;
    long pos = l;
    for (long i = l; i <= r && i < len; i ++) {
        if(u128_eq(p, list + i, eps)) {
            list[i] = *p;
            
            return i;
        }
        
        right = right && u128_smore(p, list + i);
        pos += right ? 1 : 0;
    }
    
    return -pos - 1;
}

/// @brief Merges the last two bars of the set @c ps.
///
/// @param ps the set
static void nset_mergeLast(nset ps) {
    int l = ps->barCount - 1;
    long len1 = ps->barLen[l - 1];
    long len2 = ps->lastCount;
    
    ulong mem = sizeof(u128c_struct) * (ps->count + ps->barLen[l] + len1) + sizeof(nSet_struct);
    ps->maxMem = mem > ps->maxMem ? mem : ps->maxMem;
    
    u128c_struct *list = ps->pts[l - 1];
    ps->pts[l - 1] = nset_merge(list, len1, ps->pts[l], len2);
    free(list);
    free(ps->pts[l]);
    ps->pts[l] = NULL;
    
    ps->barLen[l - 1] = len1 + len2;
    ps->barLen[l] = 0;
    
    ps->lastCount = len1 + len2;
    ps->barCount --;
}

/// @brief Returns the length of the union of the two ordered lists, up to distance @c eps.
///
/// @param u the first list
/// @param uLen the length of the first list
/// @param v the second list
/// @param vLen the length of the second list
/// @param eps the distance up to which the points are considered equal
///
/// @return the length of the union of the two lists
static long count(u128 u, long uLen, u128 v, long vLen, unsigned long eps) {
    long i = 0, j = 0, dp = 0;
    while(i < uLen || j < vLen) {
        if(i < uLen && j < vLen && u128_eq(u + i, v + j, eps)) {
            j ++;
            
            continue;
        }
        
        if(i < uLen && (j == vLen || u128_sless(u + i, v + j))) {
            i ++;
        } else {
            j ++;
        }
        
        dp ++;
    }
    
    return dp;
}

/// @brief Repeatedly checks if the last bar is at least as long as the previous last bar and merges them in this case.
///
/// In all cases, it adds a new bar with the default capacity @c NSET_DEF_SIZE at the end of the list.
///
/// @param ps the set
///
/// @return @ref true if successfull, @ref false otherwise
static bool nset_pack(nset ps) {
    while(ps->barCount > 1 && ps->lastCount >= ps->barLen[ps->barCount - 2]) {
        nset_mergeLast(ps);
    }
    
    if(ps->barCount == NSET_MAX_BARS) {
        return false;
    }
    
    // grow the points list
    ps->barLen[ps->barCount] = NSET_DEF_SIZE;
    ps->pts[ps->barCount ++] = malloc(NSET_DEF_SIZE * sizeof(u128c_struct));
    ps->lastCount = 0;
    
    ulong mem = sizeof(u128c_struct) * (ps->count + NSET_DEF_SIZE) + sizeof(nSet_struct);
    ps->maxMem = mem > ps->maxMem ? mem : ps->maxMem;
    
    return true;
}

/// @brief Inserts the point @c p into the set @c ps.
///
/// @param ps the set
/// @param p the point
///
/// @return @ref true if successfull, @ref false otherwise
static bool nset_insert(nset ps, u128 p) {
    long pos = nset_bar_search(ps, ps->barCount - 1, p);
    if(pos >= 0) {
        return false;
    }
    
    pos = -pos - 1;
    long len = ps->lastCount;
    u128c_struct* v = ps->pts[ps->barCount - 1];
    for (long i = len; i > pos; i--) {
        v[i] = v[i - 1];
    }
    v[pos] = *p;
    
    ps->lastCount ++;
    ps->count ++;
    
    ulong eps2 = ps->eps >> 1;
    ps->realCount += p->y <= eps2 ? 1 : 0;
        
    return true;
}

/// @brief Returns the distance between @c a and @c b.
///
/// @param a a point
/// @param b another point
///
/// @return the distance between @c a and @c b
static inline uint128 u128_dist(u128 a, u128 b) {
    uint128 dx = a->x <= b->x ? b->x - a->x : a->x - b->x;
    uint128 dy = a->y <= b->y ? b->y - a->y : a->y - b->y;
    
    ldbl x = uint128_to_ldbl(dx);
    ldbl y = uint128_to_ldbl(dy);

    ldbl d = sqrtl(x * x + y * y);
    
    return d >= 4 ? U128_MAX : ldbl_to_uint128(d);
}

// MARK: basic operations

bool nset_point(u128 p, nset ps, long pInd) {
    if(p == NULL || ps == NULL || pInd < 0 || pInd >= ps->count) {
        return false;
    }
    
    long pos = pInd;
    int bar = 0;
    
    while(pos >= ps->barLen[bar] && bar < ps->count)
        pos -= ps->barLen[bar ++];
    
    if(pos < 0 || bar >= ps->barCount || (bar == ps->barCount - 1 && pos >= ps->lastCount)) {
        return false;
    }
    
    *p = ps->pts[bar][pos];
    
    return true;
}

bool nset_get(mpc p, nset ps, long pInd) {
    if(p == NULL) {
        return false;
    }
    
    u128 u;
    if(! nset_point(u, ps, pInd)) {
        return false;
    }
    
    return u128_get(p, u);
}

bool nset_contains(nset ps, u128 p) {
    if(ps == NULL || p == NULL || ps->count == 0) {
        return false;
    }
    
    __uint128_t max = U128_MAX - ps->eps;
    
    if(p->x >= max || p->y >= max) {
        return false;
    }
            
    for (int i = ps->barCount - 1; i >= 0; i --) {
        if(nset_bar_search(ps, i, p) >= 0) {
            return true;
        }
    }
    
    return false;
}

bool nset_add(nset ps, u128 p) {
    __uint128_t max = U128_MAX - ps->eps;
    
    if(ps->locked || nset_contains(ps, p) || p->x >= max || p->y >= max) {
        ps->rejected ++;

        return false;
    }
    
    if(ps->lastCount == NSET_DEF_SIZE) {
        if(! nset_pack(ps)) {
            ps->rejected ++;
            
            return false;
        }
    }
    
    return nset_insert(ps, p);
}

bool nset_put(nset ps, mpc p) {
    if(p == NULL) {
        return false;
    }
    
    u128 c;
    u128_set(c, p);
    
    return nset_add(ps, c);
}

bool nset_replace(nset ps, u128 p) {
    if(ps == NULL || p == NULL || ps->count == 0) {
        return false;
    }
    
    __uint128_t max = U128_MAX - ps->eps;
    
    if(p->x >= max || p->y >= max) {
        return false;
    }
            
    for (int i = ps->barCount - 1; i >= 0; i --) {
        if(nset_searchAndReplace(ps, i, p) >= 0) {
            return true;
        }
    }
    
    return false;
}

bool nset_set_eps(nset ps, ulong eps) {
    if(! nset_valid(ps) || eps == 0 || eps > ps->eps) {
        return false;
    }
    
    if(eps == ps->eps) {
        return true;
    }
    
    ulong hoe = ps->eps >> 1;
    ps->eps = eps;
    if(ps->realCount == 0) {
        return true;
    }
    
    long rc = 0, found = 0;
    u128 p;
    ulong heps = eps >> 1;
    for (long i = 0; i < ps->count && found < ps->realCount; i++) {
        nset_point(p, ps, i);
        
        found += p->y <= hoe ? 1 : 0;
        rc += p->y <= heps ? 1 : 0;
    }
    
    ps->realCount = rc;
    
    return true;
}

bool nset_union(nset ps, nset ns, bool quick) {
    if(! nset_valid(ps) || ! nset_valid(ns) || (! ps->locked && ps->count > 0)
       || ! ns->locked) {
        return false;
    }

    if(ns->count == 0 || ps == ns) {
        return true;
    }

    if(ps->eps > ns->eps && ! nset_set_eps(ps, ns->eps)) {
        return false;
    }
    
    long uLen = ps->count;
    u128c_struct* u = ps->pts[0];

    long vLen = ns->count;
    u128c_struct* v = ns->pts[0];

    long size = quick ? uLen + vLen : count(u, uLen, v, vLen, ps->eps);
    u128c_struct* d = malloc(size * sizeof(u128c_struct));
    
    ulong mem = sizeof(u128c_struct) * (ps->count + size) + sizeof(nSet_struct);
    ps->maxMem = mem > ps->maxMem ? mem : ps->maxMem;
    
    long i = 0, j = 0, dp = 0;
    ulong eps2 = ps->eps >> 1;
    while(i < uLen || j < vLen) {
        while(i < uLen && j < vLen && u128_eq(u + i, v + j, ps->eps)) {
            j ++;
        }

        if(i < uLen && (j == vLen || u128_sless(u + i, v + j))) {
            d[dp] = u[i ++];
        } else {
            d[dp] = v[j ++];
            ps->realCount += d[dp].y <= eps2 ? 1 : 0;
        }

        dp ++;
    }
    
    ps->count = dp;
    ps->barLen[0] = dp;
    ps->lastCount = dp;
    ps->barCount = 1;
    ps->locked = true;
    
    free(ps->pts[0]);
    ps->pts[0] = d;
    
    return true;
}

nset nset_intersection(nset ps, nset ns) {
    if(ps == NULL || ns == NULL || ! ps->locked || ! ns->locked || ps->eps != ns->eps) {
        return NULL;
    }
    // TODO: implement
    
    return NULL;
}

bool nset_last(u128 p, nset ps) {
    if(ps == NULL || p == NULL || ! ps->locked || ps->count == 0) {
        return false;
    }
    
    *p = ps->pts[0][ps->count - 1];
    
    return true;
}

bool nset_left(u128 p, nset ps) {
    if(ps == NULL || p == NULL || ! ps->locked) {
        return false;
    }
    
    if(ps->count == 0) {
        return true;
    }
            
    return u128_sless(p, ps->pts[0]);
}

bool nset_right(u128 p, nset ps) {
    if(ps == NULL || p == NULL || ! ps->locked) {
        return false;
    }
    
    if(ps->count == 0) {
        return true;
    }
            
    return u128_smore(p, ps->pts[0] + (ps->count - 1));
}

bool nset_init_subset(nset dst, nset src, ulong st, ulong en) {
    if(dst == NULL || src == NULL || en < st || en > src->count || ! src->locked) {
        return false;
    }
        
    ulong size = en - st;
    dst->count = size;
    dst->barLen[0] = size;
    dst->lastCount = size;
    dst->barCount = 1;
    dst->eps = src->eps;
    dst->locked = true;
    dst->rejected = 0;
    
    if(size == 0) {
        dst->realCount = 0;
        dst->pts[0] = NULL;
        dst->barCount = 0;
        
        return true;
    }
    
    ulong block = sizeof(u128c_struct);
    u128c_struct *d = malloc(size * block);
    
    memcpy(d, src->pts[0] + st, size * block);
    dst->pts[0] = d;
    
    ulong mem = size * block + sizeof(nSet_struct);
    dst->maxMem = mem;
    
    ulong eps2 = src->eps >> 1;
    dst->realCount = 0;
    for (ulong i = 0; i < size; i++) {
        dst->realCount += d[i].y <= eps2 ? 1 : 0;
    }
    
    return true;
}

bool nset_interval(nset dst, nset src, u128 l, u128 r) {
    if(! nset_valid(src) || ! nset_valid(dst) || l == NULL || r == NULL || ! src->locked ||
       dst == src) {
        return false;
    }
    
    nset_clear(dst);
    
    if(src->count == 0 || nset_left(r, src) || nset_right(l, src)) {
        return true;
    }
    
    ulong st = nset_pos(l, src->pts[0], src->lastCount); // pl < src->lastCount, by the test above
    ulong en = nset_pos(r, src->pts[0], src->lastCount);
    
    return nset_init_subset(dst, src, st, en);
}

bool nset_split(nset left, nset right, nset src, u128 m) {
    u128c_struct z = {.x = 0, .y = 0}, o = {.x = U128_MAX, .y = U128_MAX};
    
    bool ok = nset_interval(left, src, &z, m);
    ok = ok && nset_interval(right, src, m, &o);
    
    if(! ok) {
        nset_clear(left);
        nset_clear(right);
    }
    
    return ok;
}

nset *nset_divide(nset src, uint parts) {
    if(src == NULL || parts == 0 || src->count < parts || ! src->locked) {
        return NULL;
    }
    
    nset *list = malloc(parts * sizeof(nset));
    ulong st = 0, en = src->count / parts;
    for (uint i = 0; i < parts; i++) {
        list[i] = malloc(sizeof(nSet_struct));
        
        if(! nset_init_subset(list[i], src, st, en)) {
            for (int j = 0; j < i; j++) {
                nset_free(list[j]);
            }
            
            free(list);
            
            return NULL;
        }
        
        st = en;
        uint128 e = i + 2;
        e *= src->count;
        e /= parts;
        en = (ulong) e;
    }
    
    return list;
}

bool nset_copy(nset dst, nset src) {
    if(dst == src) {
        return true;
    }
    
    if(! nset_valid(src) || ! nset_valid(dst) || ! src->locked) {
        return false;
    }
    
    nset_clear(dst);
    
    long size = src->count;
    if(size == 0) {
        return true;
    }
    
    dst->count = size;
    dst->barLen[0] = size;
    dst->lastCount = size;
    dst->barCount = 1;
    dst->eps = src->eps;
    dst->realCount = src->realCount;
    
    long block = sizeof(u128c_struct);
    u128c_struct *d = malloc(size * block);
    
    memcpy(d, src->pts[0], size * block);
    dst->pts[0] = d;
    
    ulong mem = size * block + sizeof(nSet_struct);
    dst->maxMem = mem > dst->maxMem ? mem : dst->maxMem;
    
    return true;
}

// b > a
static inline uint128 u128_d(u128 a, u128 b, uint128 md) {
    if(b->x - a->x >= md || (b->y > a->y ? b->y - a->y >= md : a->y - b->y >= md)) {
        return md;
    }
    
    return u128_dist(a, b);
}

ldbl nset_min_dist(nset ps) {
    if(ps == NULL || ! ps->locked) {
        return -1;
    }
    
    if(ps->count < 2) {
        return HUGE_VALL;
    }
    
    u128_ptr list = ps->pts[0];
    ulong n = ps->lastCount;
    ulong k = 0, l = 1;
        
    uint128 md = u128_dist(list + k, list + l), d, xk = list[k].x, xl = list[l].x;
    while(l < n - 1) {
        if(k == l - 1 || xl - xk <= md) {
            l++;
            xl = list[l].x;
            d = u128_d(list + k, list + l, md);
            md = d < md ? d : md;
        } else {
            k++;
            xk = list[k].x;
            for (ulong i = k + 1; i <= l; i++) {
                d = u128_d(list + k, list + i, md);
                md = d < md ? d : md;
            }
        }
    }
    
    for (ulong j = k + 2; j <= l; j++) {
        for (ulong i = k + 1; i < j; i++) {
            d = u128_d(list + i, list + j, md);
            md = d < md ? d : md;
        }
    }
    
    return uint128_to_ldbl(md);
}

ldbl nset_min_dist_lr(nset ls, nset rs) {
    if(ls == NULL || rs == NULL || ! ls->locked || ! rs->locked) {
        return -1;
    }
    
    if(ls->count == 0 || rs->count == 0) {
        return HUGE_VALL;
    }
    
    u128 r, l;
    nset_last(l, ls); 
    nset_point(r, rs, 0);
    if(u128_leq(r, l)) {
        return -1;
    }
    
    uint128 lx = l->x;
    uint128 rx = r->x;
    uint128 md = uint128_add(u128_dist(l, r), 1), d;
    
    bool outl = false, outr = false;
    long il = ls->count - 1, ir = 0;
    u128c_struct *pl, *pr;
    while(! outl || ! outr) {
        if(! outl) {
            if(il > 0) {
                il --;
                pl = ls->pts[0] + il;
                outl = rx - pl->x > md;
                
                for (long i = 0; i <= ir; i++) {
                    pr = rs->pts[0] + ir;
                    d = u128_dist(pl, pr);
                    
                    md = d < md ? d : md;
                }
            } else {
                outl = true;
            }
        }
        
        if(! outr) {
            if(ir < rs->count - 1) {
                ir ++;
                pr = rs->pts[0] + ir;
                outr = pr->x - lx > md;
                
                long last = ls->count - 1;
                for (long i = last; i >= il; i--) {
                    pl = ls->pts[0] + il;
                    d = u128_dist(pl, pr);
                    
                    md = d < md ? d : md;
                }
            } else {
                outr = true;
            }
        }
    }
    
    return uint128_to_ldbl(md);
}

// MARK: IO operations

bool nset_read_csv(nset ps, char *fileName) {
    if(! nset_valid(ps)) {
        return false;
    }
    
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return false;
    }
    
    mpfr_t x, y;
    mpfr_init2(x, 128);
    mpfr_init2(y, 128);
    
    u128 p;
    char sx[1550], sy[1550], line[3200], *ll = line;
    while(! feof(f)) {
        ll = fgets(line, 3200, f);
        if((ll != NULL) && sscanf(line, "%s%s", sx, sy) == 2) {
            long lx = strlen(sx);
            if(lx > 1 && sx[lx - 1] == ',') {
                sx[lx - 1] = 0;
            }
            
            bool ok = 0 == mpfr_set_str(x, sx, 10, MPFR_RNDN);
            ok = ok && 0 == mpfr_set_str(y, sy, 10, MPFR_RNDN);
            
            if(ok && u128_setr(p, x, y)) {
                nset_add(ps, p);
            }
        }
    }
    
    fclose(f);
    
    mpfr_clear(x);
    mpfr_clear(y);
    
    return true;
}

bool nset_write_csv(nset ps, char *fileName, int digits, bool append) {
    if(fileName == NULL || ps == NULL) {
        return false;
    }
    
    return nset_write_partial_csv(ps, fileName, digits, 0, ps->count, append);
}

bool nset_write_partial_csv(nset ps, char *fileName, int digits, long start, long end, bool append) {
    if(fileName == NULL || ! nset_valid(ps) || digits < 15 || digits > 40 ||
       start < 0 || end <= start || end > ps->count) {
        return false;
    }
    
    FILE *f = fopen(fileName, append ? "a" : "w");
    if(f == NULL) {
        return false;
    }
    
    mpfr_t x, y;
    mpfr_init2(x, 128);
    mpfr_init2(y, 128);
    
    u128 p;
    char rf[100], cf[100];
    sprintf(rf, "%%.%dRf, 0\n", digits);
    sprintf(cf, "%%.%dRf, %%.%dRf\n", digits, digits);

    bool err = false;
    ulong eps = ps->eps / 2;
    for (long i = start; i < end && ! err; i ++) {
        nset_point(p, ps, i);

        u128_getr(x, y, p);

        if(p->y < eps) {
            err = mpfr_fprintf(f, rf, x) < digits + 5;
        } else {
            err = mpfr_fprintf(f, cf, x, y) < 2 * digits + 6;
        }
    }

    fclose(f);
    
    mpfr_clear(x);
    mpfr_clear(y);
    
    return ! err;
}

/// @brief Removes the points of the set @c ps and then loads the content of the binary file @c f starting at the current position.
///
/// After this operation, the set @c ps is @b locked.
///
/// @param ps the set of points
/// @param f the file
/// @param count the number of values to read from the file
///
/// @return @ref true if successfull, @ref false otherwise
static bool nset_read_points(nset ps, FILE *f, ulong count) {
    if(f == NULL || ! nset_valid(ps) || count <= 0 || count > (1L << 26)) {
        return false;
    }

    long block = sizeof(u128c_struct);
    nset_clear(ps);
    
    ps->count = count;
    ps->locked = 1;
    ps->barCount = 1;
    ps->barLen[0] = count;
    ps->lastCount = count;
    
    long msize = count * block;
    ps->pts[0] = malloc(msize);
    
    ulong mem = msize + sizeof(nSet_struct);
    ps->maxMem = mem > ps->maxMem ? mem : ps->maxMem;
    
    return fread_block(ps->pts[0], block, count, FILE_BLOCK, f) == count;
}

bool nset_read(nset ps, char *fileName, bool checkMD5) {
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return false;
    }
    
    bool ok = nset_read_from(ps, f, 0, checkMD5);
    
    fclose(f);
    
    return ok;
}

nset nset_load(char *fileName, bool checkMD5) {
    nset ps = nset_new(0, true);
    
    bool ok = nset_read(ps, fileName, checkMD5);
    
    if(ok) {
        return ps;
    }
    
    nset_free(ps);
    
    return NULL;
}

bool nset_read_partial(nset ps, char *fileName, long pos, long start, ulong count) {
    if(fileName == NULL || ! nset_valid(ps)) {
        return false;
    }
    
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return false;
    }
    
    bool ok = nset_read_segment(ps, f, pos, start, count);
    
    fclose(f);
    
    return ok;
}

bool nset_read_old(nset ps, char *fileName) {
    if(ps == NULL || fileName == NULL) {
        return false;
    }
    
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return false;
    }
    
    long block = sizeof(u128c_struct);
    
    // read the points, no header is present in the old format
    long fs = io_file_size(f);
    bool ok = fs > 0 && (fs % block == 0);
    ok = ok && nset_read_points(ps, f, fs / block);
    
    fclose(f);
    
    if(ok) {
        ulong eps2 = ps->eps >> 1;
        ps->realCount = 0;
        
        for (long i = 0; i < ps->count; i++) {
            ps->realCount += ps->pts[0][i].y <= eps2 ? 1 : 0;
        }
    }
    
    return ok;
}

bool nset_write(nset ps, char *fileName) {
    if(! nset_valid(ps) || ! ps->locked || ps->count > ps->barLen[0] || fileName == NULL || ps->count == 0) {
        printf("nset valid = %d, locked = %d, count = %lu, barLen[0] = %lu, fn = \"%s\"\n",
               nset_valid(ps), ps->locked, ps->count, ps->barLen[0], fileName);
        
        return false;
    }

    FILE *f = fopen(fileName, "w");
    ulong wr = nset_write_to(ps, f, 0);
    fclose(f);
        
    ulong exp = NSET_HEADER_LEN + sizeof(u128c_struct) * ps->count;
    if(wr != exp)  {
        printf("Written %lu bytes to \"%s\'instead of expected %lu !\n", wr, fileName, exp);
        
        return false;
    }
    
    return true;
}

bool nset_read_from(nset ps, FILE *f, long pos, bool checkMD5) {
    mfile h = mfile_read_header(f, pos, NSET_FILE_ID, checkMD5 ? NSET_HEADER_LEN : NSET_HEADER_LEN - 16);
    if(h == NULL) {
        return false;
    }
    
    h->pos = 8; // skip the file ID
    ulong flen = mfile_getl(h);
    uint hlen = mfile_geti(h);
    
    ulong count = mfile_getl(h);
    ulong realCount = mfile_getl(h);
    ulong eps = mfile_getl(h);
            
    ulong pslen = hlen + count * sizeof(u128c_struct);
    if(flen != pslen || ! nset_read_points(ps, f, count)) {
        nset_clear(ps);
        mfile_free(h);
        
        return false;
    }
    
    ps->count = count;
    ps->realCount = realCount;
    ps->eps = eps;
    
    ulong dlen = ((long) ps->count) * sizeof(u128c_struct);
    bool ok = ! checkMD5 || mfile_header_check_md5(h, ps->pts[0], dlen);
    
    mfile_free(h);
    
    return ok;
}

/// @brief Writes the file header of the set @c ps to the file @c f at absolute position @c pos, or the current position if
/// @c pos<0.
///
/// @param ps the set
/// @param f the file
/// @param pos the absolute position, the current position if @c pos<0
///
/// @return @ref true if successfull, @ref false otherwise
static bool nset_write_header(nset ps, FILE *f, long pos) {
    if(f == NULL) { // ps already checked
        return false;
    }
        
    uint hlen = NSET_HEADER_LEN;
    ulong flen = hlen + sizeof(u128c_struct) * ps->count;
    mfile h = mfile_header(NSET_FILE_ID, flen, hlen);
    if(h == NULL) {
        return false;
    }
    
    mfile_putl(h, ps->count);
    mfile_putl(h, ps->realCount);
    mfile_putl(h, ps->eps);
    mfile_header_add_md5(h, ps->pts[0], ((long) ps->count) * sizeof(u128c_struct));
        
    bool ok = mfile_write_to(f, pos, h, 0, h->len);
    mfile_free(h);
    
    return ok;
}

ulong nset_write_to(nset ps, FILE *f, long pos) {
    if(! nset_valid(ps) || ! ps->locked || ps->count > ps->barLen[0] || f == NULL || ps->count == 0) {
        return 0;
    }
    
    if(! nset_write_header(ps, f, pos)) {
        return 0;
    }
    
    ulong block = sizeof(u128c_struct);
    ulong wr = fwrite(ps->pts[0], block, ps->count, f);
    
    return NSET_HEADER_LEN + block * wr; // 16 for md5 sum
}

bool nset_read_segment(nset ps, FILE *f, long pos, long start, ulong count) {
    // no MD5 is needed for a partial read, so the header may be 16 bytes shorter
    mfile h = mfile_read_header(f, pos, NSET_FILE_ID, NSET_HEADER_LEN - 16);
    if(h == NULL) {
        return false;
    }
    
    h->pos = 8; // skip de file ID, go to the length of the logical nset file
    ulong flen = mfile_getl(h);
    int hlen = mfile_geti(h);
    ulong plen = flen - hlen; // length of the points data in the file
    
    ps->count = mfile_getl(h);
    ps->realCount = mfile_getl(h);
    ps->eps = mfile_getl(h);
    mfile_free(h);
    
    // check corectness of the data and the size of the file
    long block = sizeof(u128c_struct);
    if(! nset_valid(ps) || count <= 0 || io_file_remainder(f) < plen ||
       plen <= 0 || plen != block * ps->count) {
        return false;
    }
    
    // check the compatibility of the parameters and seek the correct position in f
    long len = plen / block;
    long fp = ftell(f);
    bool ok = fp >= NSET_HEADER_LEN - 16; // no MD5 is needed for a partial read
    ok = ok && (start < 0 || (start + count <= len &&
                              fseek(f, fp + block * start, SEEK_SET) == 0));
    ok = ok && (start >= 0 || (-start <= len && -start <= count &&
                               fseek(f, fp + plen - block * count, SEEK_SET) == 0));
    
    if(! ok) {
        return false;
    }
    
    ok = nset_read_points(ps, f, count);
    
    if(! ok) {
        nset_clear(ps);
    } else if(count < len) { // count the real roots
        ulong eps2 = ps->eps >> 1;
        ps->realCount = 0;
        
        for (long i = 0; i < count; i++) {
            ps->realCount += ps->pts[0][i].y <= eps2 ? 1 : 0;
        }
    }
    
    return ok;
}
