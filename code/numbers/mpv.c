//
//  mpv.c
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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "mpv.h"
#include "memFile.h"
#include "io.h"

static volatile uint mpvInited = 0;

// MARK: Creation, destruction, cut and paste vectors

mpv mpv_new(uint prec, ulong count) {
    if(prec < MPV_MIN_PREC) {
        return NULL;
    }
    
    mpv v = malloc(sizeof(mpVector));
    mpv_init(v, prec, count);
    
    return v;
}

bool mpv_init(mpv v, uint prec, ulong count) {
    if(v == NULL || prec < MPV_MIN_PREC) {
        return false;
    }
    
    if(! mpvInited) {
        mpvInited = 1;
        
        mpfr_set_emax(MPV_MAX_EXP);
        mpfr_set_emin(MPV_MIN_EXP);
    }
    
    v->prec = prec;
    v->count = count;
    
    uint limbs = mpv_limbs(prec);
    v->limbs = limbs;
    
    if(count > 0) {
        v->vals = malloc(8 * limbs * count);
        v->sexp = malloc(4 * count);
    } else {
        v->vals = NULL;
        v->sexp = NULL;
    }
    
    for (ulong i = 0; i < count; i++) {
        v->sexp[i] = MPV_NAN_EXP << 1;
    }
    
    return true;
}

void mpv_clear(mpv vect) {
    if(vect == NULL || vect->prec < MPV_MIN_PREC) {
        return;
    }
    
    vect->prec = 0;
    vect->count = 0;
    vect->limbs = 0;
    
    if(vect->vals != NULL) {
        free(vect->vals);
        
        vect->vals = NULL;
    }
    
    if(vect->sexp != NULL) {
        free(vect->sexp);
        
        vect->sexp = NULL;
    }
}

void mpv_free(mpv vect) {
    if(vect == NULL) {
        return;
    }
    
    mpv_clear(vect);
    free(vect);
}

/// @brief Copies the sign and exponent from  the position @c pos of the vector @c vect to the number @c src.
///
/// @param dst the number
/// @param vect the vector
/// @param pos the position in the vector
///
/// @return @ref true if @c dst is a regular number, @ref false otherwise
static inline bool mpv_quse(mpfr_t dst, mpv vect, ulong pos) {
    int exp = vect->sexp[pos];
    
    dst->_mpfr_sign = exp & 1 ? -1 : 1;
    exp >>= 1;
    
    if(exp >= MPV_MIN_EXP && exp <= MPV_MAX_EXP) {
        dst->_mpfr_exp = exp;
        
        return true;
    }
    
    dst->_mpfr_exp = exp == MPV_ZERO_EXP ? __MPFR_EXP_ZERO :
        exp == MPV_INF_EXP ? __MPFR_EXP_INF : __MPFR_EXP_NAN;
    
    return false;
}

inline void mpv_fuse(mpfr_t dst, mpv vect, ulong pos) {
    dst->_mpfr_prec = vect->prec;
    mpv_quse(dst, vect, pos);
        
    dst->_mpfr_d = (mp_limb_t *) (vect->vals + pos * vect->limbs);
}

bool mpv_copy(mpv dst, ulong dpos, mpv src) {
    if(dst == NULL || src == NULL  || dpos >= dst->count ||
       dst->count - dpos < src->count) {
        return false;
    }
    
    if(dst->prec != src->prec) {
        mpfr_t buf;
        
        for (ulong i = 0; i < src->count; i++) {
            mpv_fuse(buf, src, i);
            mpv_set(dst, i + dpos, buf);
        }
    } else if(src != dst) {
        uint limbs = dst->limbs;
        memcpy(dst->vals + (limbs * dpos), src->vals, 8 * limbs * src->count);
        memcpy(dst->sexp + dpos, src->sexp, 4 * src->count);
    } // if src == dst we are here only if dpos == 0, there is nothing to do !
    
    return true;
}

int mpv_get_2exp(mpv vect, ulong pos) {
    if(vect == NULL || vect->count <= pos) {
        return INT_MIN;
    }
    
    int exp = vect->sexp[pos] >> 1;
    
    return exp == MPV_NAN_EXP ? INT_MIN : exp;
}

bool mpv_copy_partial(mpv dst, ulong dpos, mpv_t src, ulong spos, ulong count) {
    if(dst == NULL || src == NULL  || dpos + count > dst->count ||
       spos + count > src->count) {
        return false;
    }
    
    if(count == 0) {
        return true;
    }
    
    if(dst->prec != src->prec) {
        mpfr_t buf;
        
        for (ulong i = 0; i < count; i++) {
            mpv_fuse(buf, src, i + spos);
            mpv_set(dst, i + dpos, buf);
        }
    } else if(src != dst) {
        uint limbs = dst->limbs;
        memcpy(dst->vals + (limbs * dpos), src->vals + (limbs * spos), 8 * limbs * count);
        memcpy(dst->sexp + dpos, src->sexp + spos, 4 * count);
    } else if(spos > dpos) {
        uint limbs = dst->limbs;
        
        for (ulong i = 0; i < count; i++) {
            memcpy(dst->vals + (limbs * (dpos + i)), src->vals + (limbs * (spos + i)), 8 * limbs);
            dst->sexp[dpos + i] = src->sexp[spos + i];
        }
    } else if(spos < dpos) {
        uint limbs = dst->limbs;
        
        for (ulong j = 0; j < count; j++) {
            ulong i = count - 1 - j;
            memcpy(dst->vals + (limbs * (dpos + i)), src->vals + (limbs * (spos + i)), 8 * limbs);
            dst->sexp[dpos + i] = src->sexp[spos + i];
        }
    } // if no condition above is satisfied, then src == dst and dpos == spos, nothing to do !
    
    return true;
}

mpv mpv_clone(mpv vect) {
    if(vect == NULL) {
        return NULL;
    }
    
    mpv nv = mpv_new(vect->prec, vect->count);
    if(nv == NULL) {
        return NULL;
    }
    
    mpv_copy(nv, 0, vect);
    
    return nv;
}

bool mpv_resize(mpv vect, ulong count) {
    if(vect == NULL || count == 0) {
        return false;
    }
    
    if(count == vect->count) {
        return true;
    }
    
    ulong limbs = vect->limbs;
    
    ulong *vals = vect->vals;
    vect->vals = realloc(vals, 8 * limbs * count);
    if(vect->vals == NULL) {
        vect->vals = vals;
        
        return false;
    }
    
    int *sexp = vect->sexp;
    vect->sexp = realloc(sexp, 4 * count);
    if(vect->sexp == NULL) {
        vect->sexp = sexp;
        
        return false;
    }
    
    for (ulong i = vect->count; i < count; i++) {
        vect->sexp[i] = MPV_NAN_EXP << 1;
    }
    
    vect->count = count;
    
    return true;
}

bool mpv_concat(mpv dst, mpv src) {
    if(dst == NULL || src == NULL || dst->prec != src->prec) {
        return false;
    }
    
    ulong dcount = dst->count;
    if(! mpv_resize(dst, dst->count + src->count)) {
        return false;
    }
    
    ulong limbs = dst->limbs;
    memcpy(dst->vals + limbs * dcount, src->vals, 8 * limbs * src->count);
    memcpy(dst->sexp + dcount, src->sexp, 4 * src->count);
    
    return true;
}

mpv mpv_join(mpv v1, mpv_t v2) {
    if(v1 == NULL || v2 == NULL || v1->prec != v2->prec) {
        return NULL;
    }
    
    mpv nv = mpv_new(v1->prec, v1->count + v2->count);
    if(nv == NULL) {
        return NULL;
    }
    
    ulong limbs = v1->limbs;
    memcpy(nv->vals, v1->vals, 8 * limbs * v1->count);
    memcpy(nv->sexp, v1->sexp, 4 * v1->count);
    memcpy(nv->vals + 8 * limbs * v1->count, v2->vals, 8 * limbs * v2->count);
    memcpy(nv->sexp + 4 * v1->count, v2->sexp, 4 * v2->count);
    
    return nv;
}

mpv mpv_sub_vector(mpv vect, ulong start, ulong step, ulong count) {
    if(vect == NULL || step == 0 || count == 0 || start >= vect->count
       || (vect->count - start - 1) / step < count - 1) {
        return NULL;
    }
    
    mpv v = mpv_new(vect->prec, count);
    ulong limbs = v->limbs;
    
    for (ulong i = 0; i < count; i++) {
        ulong is = i * step;
        ulong ps = start + is;
        v->sexp[i] = vect->sexp[ps];
        
        ulong ls = limbs * ps;
        ulong li = limbs * i;
        for (uint j = 0; j < limbs; j++) {
            v->vals[li + j] = vect->vals[ls + j];
        }
    }
    
    return v;
}

mpv mpv_sub_vectorc(mpv vect, ulong start, ulong step, ulong count) {
    ulong tstart = start << 1;
    
    if(vect == NULL || step == 0 || count == 0 || tstart >= vect->count
       || (vect->count & 1) != 0 || (vect->count - tstart - 2) / step < count - 1) {
        return NULL;
    }
    
    mpv v = mpv_new(vect->prec, count << 1);
    ulong limbs = v->limbs;
    
    ulong tl = limbs << 1;
    ulong ts = start << 1;
    for (ulong i = 0; i < count; i++) {
        ulong ti = i << 1;
        ulong is = ti * step;
        ulong ps = ts + is;
        v->sexp[ti] = vect->sexp[ps];
        v->sexp[ti + 1] = vect->sexp[ps + 1];
        
        ulong ls = limbs * ps;
        ulong li = limbs * ti;
        for (uint j = 0; j < tl; j++) {
            v->vals[li + j] = vect->vals[ls + j];
        }
    }
    
    return v;
}

// MARK: Storage and recovery of numbers

inline bool mpv_qsync(mpv vect, ulong pos, mpfr_t src) {
    int exp;
    if(mpfr_regular_p(src)) {
        exp = (int) src->_mpfr_exp;
        vect->sexp[pos] = (exp << 1) | (src->_mpfr_sign < 0 ? 1 : 0);
        
        return true;
    }
    
    exp = mpfr_zero_p(src) ? MPV_ZERO_EXP :
        mpfr_nan_p(src) ? MPV_NAN_EXP : MPV_INF_EXP;
    
    vect->sexp[pos] = (exp << 1) | (src->_mpfr_sign < 0 ? 1 : 0);
    
    return false;
}

bool mpv_set(mpv vect, ulong pos, mpfr_t src) {
    if(vect == NULL || pos >= vect->count || src == NULL) {
        return false;
    }
        
    int limbs = vect->limbs;
    
    if(src->_mpfr_prec == vect->prec) {
        if(mpv_qsync(vect, pos, src)) {
            memcpy(vect->vals + pos * limbs, src->_mpfr_d, limbs * 8);
        }
    } else {
        mpfr_t buf;
        mpv_fuse(buf, vect, pos);
        
        mpfr_set(buf, src, MPFR_RNDN);
        mpv_qsync(vect, pos, buf);
    }
    
    return true;
}

bool mpv_neg(mpv vect, ulong pos) {
    if(vect == NULL || pos >= vect->count) {
        return false;
    }
    
    vect->sexp[pos] ^= 1;
    
    return true;
}

bool mpv_conj(mpv vect, ulong pos) {
    ulong cpos = (pos << 1) + 1;
    if(vect == NULL || cpos >= vect->count) {
        return false;
    }
    
    vect->sexp[cpos] ^= 1;
    
    return true;
}

bool mpv_setl(mpv vect, ulong pos, ldbl src) {
    if(vect == NULL || pos >= vect->count) {
        return false;
    }
    
    mpfr_t buf;
    mpv_fuse(buf, vect, pos);
    
    mpfr_set_ld(buf, src, MPFR_RNDN);
    mpv_qsync(vect, pos, buf);
    
    return true;
}

bool mpv_set_si(mpv vect, ulong pos, long src) {
    if(vect == NULL || pos >= vect->count) {
        return false;
    }
    
    mpfr_t buf;
    mpv_fuse(buf, vect, pos);
    
    mpfr_set_si(buf, src, MPFR_RNDN);
    mpv_qsync(vect, pos, buf);
    
    return true;
}

long mpv_get_si(mpv vect, ulong pos) {
    if(vect == NULL || pos >= vect->count) {
        return LONG_MIN;
    }
    
    mpfr_t buf;
    mpv_fuse(buf, vect, pos);
    
    return mpfr_get_si(buf, MPFR_RNDN);
}

bool mpv_set_zero(mpv vect, ulong pos) {
    if(vect == NULL || pos >= vect->count) {
        return false;
    }
        
    vect->sexp[pos] = MPV_ZERO_EXP << 1;
    
    return true;
}

bool mpv_set_all_zero(mpv vect) {
    if(vect == NULL) {
        return false;
    }
        
    for (ulong pos = 0; pos < vect->count; pos++) {
        vect->sexp[pos] = MPV_ZERO_EXP << 1;
    }
    
    return true;
}

bool mpv_set_nan(mpv vect, ulong pos) {
    if(vect == NULL || pos >= vect->count) {
        return false;
    }
        
    vect->sexp[pos] = MPV_NAN_EXP << 1;
    
    return true;
}

bool mpv_set_inf(mpv vect, ulong pos, bool positive) {
    if(vect == NULL || pos >= vect->count) {
        return false;
    }
        
    vect->sexp[pos] = (MPV_INF_EXP << 1) | (positive ? 0 : 1);
    
    return true;
}

bool mpv_setc(mpv vect, ulong pos, mpc src) {
    if(src == NULL) {
        return false;
    }
    
    ulong tpos = pos << 1;
    
    bool ok = mpv_set(vect, tpos, src->x);
    ok = ok && mpv_set(vect, tpos + 1, src->y);
    
    return ok;
}

bool mpv_seti(mpv vect, ulong pos, mpi src) {
    if(src == NULL) {
        return false;
    }
    
    ulong tpos = pos << 1;
    
    bool ok = mpv_set(vect, tpos, src->a);
    ok = ok && mpv_set(vect, tpos + 1, src->b);
    
    return ok;
}

bool mpv_setcl(mpv vect, ulong pos, fp80 src) {
    if(src == NULL) {
        return false;
    }
    
    ulong tpos = pos << 1;
    
    bool ok = mpv_setl(vect, tpos, src->x);
    ok = ok && mpv_setl(vect, tpos + 1, src->y);
    
    return ok;
}

bool mpv_is_nan(mpv vect, ulong pos) {
    if(vect == NULL || pos >= vect->count) {
        return false;
    }
    
    return (vect->sexp[pos] & -2) == MPV_NAN_EXP << 1;
}

bool mpv_is_zero(mpv vect, ulong pos) {
    if(vect == NULL || pos >= vect->count) {
        return false;
    }
    
    return (vect->sexp[pos] & -2) == MPV_ZERO_EXP << 1;
}

bool mpv_is_inf(mpv vect, ulong pos) {
    if(vect == NULL || pos >= vect->count) {
        return false;
    }
    
    return (vect->sexp[pos] & -2) == MPV_INF_EXP << 1;
}

bool mpv_is_reg(mpv vect, ulong pos) {
    if(vect == NULL || pos >= vect->count) {
        return false;
    }
    
    int exp = vect->sexp[pos] >> 1;
    
    return exp >= MPV_MIN_EXP && exp <= MPV_MAX_EXP;
}

bool mpv_all_finite(mpv vect) {
    if(vect == NULL) {
        return false;
    }
    
    bool fin = true;
    int exp = 0;
    for (ulong i = 0; i < vect->count && fin; i++) {
        exp = vect->sexp[i] >> 1;
        fin = (exp >= MPV_MIN_EXP && exp <= MPV_MAX_EXP) || exp == MPV_ZERO_EXP;
    }
    
    return fin;
}

bool mpv_get(mpfr_t dst, mpv vect, ulong pos) {
    if(vect == NULL || pos >= vect->count || dst == NULL) {
        return false;
    }
        
    int limbs = vect->limbs;
    
    if(dst->_mpfr_prec == vect->prec) {
        if(mpv_quse(dst, vect, pos)) { // no need to copy if NAN, +/- INF or +/- 0
            memcpy(dst->_mpfr_d, vect->vals + pos * limbs, limbs * 8);
        }
    } else {
        mpfr_t buf;
        mpv_fuse(buf, vect, pos);
        
        mpfr_set(dst, buf, MPFR_RNDN);
    }
    
    return true;
}

ldbl mpv_getl(mpv vect, ulong pos) {
    mpfr_t buf;
    if(! mpv_use(buf, vect, pos)) {
        return NAN;
    }
    
    return mpfr_get_ld(buf, MPFR_RNDN);
}

bool mpv_getc(mpc dst, mpv vect, ulong pos) {
    if(dst == NULL) {
        return false;
    }
    
    ulong tpos = pos << 1;
    
    bool ok = mpv_get(dst->x, vect, tpos);
    ok = ok && mpv_get(dst->y, vect, tpos + 1);
    
    return ok;
}

bool mpv_geti(mpi dst, mpv vect, ulong pos) {
    if(dst == NULL) {
        return false;
    }
    
    ulong tpos = pos << 1;
    
    bool ok = mpv_get(dst->a, vect, tpos);
    ok = ok && mpv_get(dst->b, vect, tpos + 1);
    
    return ok;
}

bool mpv_getcl(fp80 dst, mpv vect, ulong pos) {
    if(dst == NULL) {
        return false;
    }
    
    ulong tpos = pos << 1;
    
    ldbl x = mpv_getl(vect, tpos);
    ldbl y = mpv_getl(vect, tpos + 1);
    
    if(x != NAN && y != NAN) {
        dst->x = x;
        dst->y = y;
        
        return true;
    }
    
    return false;
}

bool mpv_get_set(mpv dst, ulong dpos, mpv src, ulong spos) {
    if(dst == NULL || src == NULL || dpos >= dst->count || spos >= src->count) {
        return false;
    }
    
    uint dprec = dst->prec;
    if(dprec == src->prec) {
        uint limbs = dst->limbs;
        
        if(src != dst || dpos != spos) {
            dst->sexp[dpos] = src->sexp[spos];
            memcpy(dst->vals + limbs * dpos, src->vals + limbs * spos, limbs << 3);
        }
            
        return true;
    }
    
    mpfr_t sbuf, dbuf;
    mpv_fuse(sbuf, src, spos);
    mpv_fuse(dbuf, dst, dpos);
    
    mpfr_set(dbuf, sbuf, MPFR_RNDN);
    mpv_qsync(dst, dpos, dbuf);
    
    return true;
}

bool mpv_get_setc(mpv dst, ulong dpos, mpv src, ulong spos) {
    ulong tspos = spos << 1;
    ulong tdpos = dpos << 1;
    
    if(dst == NULL || src == NULL || tdpos + 1 >= dst->count || tspos + 1 >= src->count) {
        return false;
    }
    
    if(dst->prec == src->prec) {
        uint limbs = dst->limbs;
        
        if(src != dst || dpos != spos) {
            dst->sexp[tdpos] = src->sexp[tspos];
            dst->sexp[tdpos + 1] = src->sexp[tspos + 1];
            memcpy(dst->vals + limbs * tdpos, src->vals + limbs * tspos, limbs << 4);
        }
        
        return true;
    }
    
    mpfr_t sbuf, dbuf;
    mpv_fuse(sbuf, src, tspos);
    mpv_fuse(dbuf, dst, tdpos);
    
    mpfr_set(dbuf, sbuf, MPFR_RNDN);
    mpv_qsync(dst, tdpos, dbuf);
    
    tspos ++;
    tdpos ++;
    mpv_fuse(sbuf, src, tspos);
    mpv_fuse(dbuf, dst, tdpos);
    
    mpfr_set(dbuf, sbuf, MPFR_RNDN);
    mpv_qsync(dst, tdpos, dbuf);
    
    return true;
}

// MARK: Shared limbs with external mpfr_t variables

bool mpv_use(mpfr_t dst, mpv vect, ulong pos) {
    if(vect == NULL || pos >= vect->count || dst == NULL) {
        return false;
    }
    
    mpv_fuse(dst, vect, pos);
    
    return true;
}

bool mpv_sync(mpv vect, ulong pos, mpfr_t src) {
    if(vect == NULL || pos >= vect->count || src == NULL) {
        return false;
    }
    
    int step = vect->limbs;
    if(src->_mpfr_d != (mp_limb_t *) (vect->vals + pos * step)) {
        return false;
    }
    
    mpv_qsync(vect, pos, src);
    
    return true;
}

bool mpv_usec(mpc dst, mpv vect, ulong pos) {
    ulong tpos = pos << 1;
    
    if(dst == NULL || vect == NULL || tpos + 1 >= vect->count) {
        return false;
    }
    
    mpv_fuse(dst->x, vect, tpos ++);
    mpv_fuse(dst->y, vect, tpos);
    
    return true;
}

bool mpv_usei(mpi dst, mpv vect, ulong pos) {
    ulong tpos = pos << 1;
    
    if(dst == NULL || vect == NULL || tpos + 1 >= vect->count ||
       mpi_prec(dst) != vect->prec) {
        return false;
    }
    
    mpv_fuse(dst->a, vect, tpos ++);
    mpv_fuse(dst->b, vect, tpos);
    
    return true;
}

bool mpv_syncc(mpv vect, ulong pos, mpc src) {
    if(src == NULL) {
        return false;
    }
    
    ulong tpos = pos << 1;
    
    bool ok = mpv_sync(vect, tpos, src->x);
    ok = ok && mpv_sync(vect, tpos + 1, src->y);
    
    return ok;
}

bool mpv_synci(mpv vect, ulong pos, mpi src) {
    if(src == NULL) {
        return false;
    }
    
    ulong tpos = pos << 1;
    
    bool ok = mpv_sync(vect, tpos, src->a);
    ok = ok && mpv_sync(vect, tpos + 1, src->b);
    
    return ok;
}

// MARK: Miscellaneous functions

fp80_ptr mpv_to_fp80(mpv v) {
    if(v == NULL || v->count <= 0 || v->count % 2 != 0) {
        return NULL;
    }
    
    fp80_ptr v80 = malloc(v->count * sizeof(fp80_struct));
    for (ulong i = 0; i < v->count / 2; i++) {
        ulong ti = i << 1;
        
        mpfr_t buf;
        
        mpv_use(buf, v, ti);
        v80[i].x = mpfr_get_ld(buf, MPFR_RNDN);
        
        mpv_use(buf, v, ti + 1);
        v80[i].y = mpfr_get_ld(buf, MPFR_RNDN);
    }
    
    return v80;
}

bool mpv_scale(mpv v, ulong pos, int tpow) {
    if(v == NULL || pos >= v->count) {
        return false;
    }
    
    long exp = v->sexp[pos];
    long re = exp >> 1;
    
    if(re == MPV_NAN_EXP) {
        return false;
    }
    
    if(re < MPV_MIN_EXP || re + tpow < MPV_MIN_EXP) {
        v->sexp[pos] = (int) (MPV_ZERO_EXP << 1) | (exp & 1);
        
        return true;
    }
    
    if(re > MPV_MAX_EXP || re + tpow > MPV_MAX_EXP) {
        v->sexp[pos] = (int) (MPV_INF_EXP << 1) | (exp & 1);
        
        return true;
    }
    
    v->sexp[pos] = (int) ((re + tpow) << 1) | (exp & 1);
    
    return true;
}

bool mpv_ulp(mpfr_t ulp, mpv vect) {
    if(ulp == NULL || vect == NULL) {
        return false;
    }
    
    long maxExp = LONG_MIN;
    for (ulong i = 0; i < vect->count; i++) {
        int exp = vect->sexp[i] >> 1;
        maxExp = exp > maxExp ? exp : maxExp;
    }
    
    mpfr_set_ui_2exp(ulp, 1, maxExp - vect->prec, MPFR_RNDU);
    
    return true;
}

// MARK: IO operations

static bool mpv_check_headers(FILE *f, uint *prec, ulong *count) {
    uint pp = 0;
    ulong len = 0;
    
    long pos = 0;
    do {
        mfile h = mfile_read_header(f, pos, MPV_FILE_ID, MPV_HEADER_LEN);
        if(h == NULL) { // EOF
            if(pp == 0 || len == 0) {
                return false;
            }
            
            *prec = pp;
            *count = len;
            
            return true;
        }
        
        h->pos = 8; // skip the file ID
        ulong flen = mfile_getl(h);
        h->pos = 20; // skip header size
        uint fprec = mfile_geti(h);
        ulong pc = mfile_getl(h);
        
        mfile_free(h);
        
        if(pos == 0) {
            pp = fprec;
        } else if(pp != fprec) {
            return false;
        }
        
        if(pc == 0 || flen == 0) {
            return false;
        }
        
        len += pc;
        pos += flen;
    } while(true);
    
    return false;
}

ulong mpv_read_from_to(FILE *f, long fpos, mpv dst, ulong dpos) {
    if(f == NULL || dst == NULL || dpos >= dst->count) {
        return 0;
    }
    
    mfile h = mfile_read_header(f, fpos, MPV_FILE_ID, MPV_HEADER_LEN);
    if(h == NULL) {
        return 0;
    }
    
    h->pos = 20; // skip the file ID, length of the file and the header length
    uint prec = mfile_geti(h);
    ulong count = mfile_getl(h);
    if(dpos  + count > dst->count || dst->prec != prec) {
        mfile_free(h);
        
        return 0;
    }
    
    // compute MD5 checksum
    MD5_CTX md5;
    byte md5sum[16];
    
    bool ok = MD5_Init(&md5) == 1;
    ok = ok && MD5_Update(&md5, h->data, (uint) h->pos) == 1;
    
    byte fmd5sum[16];
    mfile_getbs(h, fmd5sum, 16);
    
    mfile_free(h);
        
    if(! ok || ! mpv_read_points(dst, dpos, count, f, -1)) {
        MD5_Final(md5sum, &md5);
        
        return 0;
    }
    
    ok = ok && mpv_update_md5_partial(dst, &md5, dpos, count);
    ok = ok && MD5_Final(md5sum, &md5);
    
    for (int i = 0; i < 16 && ok; i++) {
        ok = ok && fmd5sum[i] == md5sum[i];
    }
    
    return ok ? count : 0;
}

mpv mpv_read(char *fileName, bool multi) {
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return NULL;
    }
        
    if(! multi) {
        mpv v = mpv_read_from(f, 0);
        fclose(f);
        
        return v;
    }
    
    uint prec;
    ulong count, loaded = 0;
    if(! mpv_check_headers(f, &prec, &count)) {
        fclose(f);
            
        return NULL;
    }
    
    mpv v = mpv_new(prec, count);
    if(v == NULL) {
        fclose(f);
        
        ulong mem = mpv_file_size(prec, count);
        char mems[100];
        mem_size(mem, mems);
        
        io_error_echo("Could not allocate %s memory.\n", mems);
        fflush(stdout);
        
        return NULL;
    }
    
    bool ok = true;
    long pos = 0;
    while(ok && loaded < count) {
        ulong new = mpv_read_from_to(f, pos, v, loaded);
        
        ok = ok && new > 0;
        loaded += new;
        
        pos = -1;
    }
    
    fclose(f);
    
    if(ok) {
        return v;
    }
    
    mpv_free(v);
    
    return NULL;
}

ulong mpv_write(mpv vect, char *fileName, int append) {
    if(vect == NULL || fileName == NULL) {
        return 0;
    }
    
    FILE *f = fopen(fileName, append ? "a" : "w");
    if(f == NULL) {
        io_error_echo("Could not %s \"%s\"\n", append ? "append to" : "create", fileName);
        
        return 0;
    }
    
    ulong bytes = mpv_write_to(vect, f, -1);
    fclose(f);
    
    return bytes;
}

bool mpv_write80(fp80_ptr v, ulong len, uint prec, char *fileName) {
    if(v == NULL || len == 0 || prec < 32 || fileName == NULL) {
        return false;
    }
    
    mpv mv = mpv_new(prec, len << 1);
    bool ok = mv != NULL;
    for (ulong i = 0; ok && i < len; i++) {
        ok = ok && mpv_setcl(mv, i, v + i);
    }
    
    ok = ok && mpv_write(mv, fileName, false);
    mpv_free(mv);
    
    return ok;
}

fp80_ptr mpv_read80(ulong *len, char *fileName) {
    mpv v = mpv_read(fileName, false);
    if(v == NULL || v->count < 2) {
        return NULL;
    }
    
    ulong n = v->count >> 1;
    fp80_ptr v8 = malloc(sizeof(fp80) * n);
    bool ok = v8 != NULL;
    for (ulong i = 0; ok && i < n; i++) {
        ok = ok && mpv_getcl(v8 + i, v, i);
    }
    
    mpv_free(v);
    
    if (! ok) {
        if(v8 != NULL) {
            free(v8);
        }
        
        *len = 0;
        
        return NULL;
    }
    
    *len = n;
    
    return v8;
}

ulong mpv_write_mini(mpv vect, char *fileName, int index) {
    if(fileName == NULL || index < 0 || vect == NULL) {
        return 0;
    }
    
    long miniFileLen = mpv_fileLen(vect);
    
    FILE *f = fopen(fileName, "r+");
    if(f == NULL) {
        f = fopen(fileName, "w");
    }
    if(f == NULL) {
        return 0;
    }
    
    ulong bytes = mpv_write_to(vect, f, miniFileLen * index);
    fclose(f);
    
    return bytes;
}

mpv mpv_read_from(FILE *f, long pos) {
    mfile h = mfile_read_header(f, pos, MPV_FILE_ID, MPV_HEADER_LEN);
    if(h == NULL) {
        return NULL;
    }
    
    h->pos = 20; // skip the file ID, length of the file and the header length
    uint prec = mfile_geti(h);
    ulong count = mfile_getl(h);
    
    // compute MD5 checksum
    MD5_CTX md5;
    byte md5sum[16];
    
    MD5_Init(&md5);
    bool ok = MD5_Update(&md5, h->data, (uint) h->pos);
    
    byte fmd5sum[16];
    mfile_getbs(h, fmd5sum, 16);
    
    mfile_free(h);
    
    mpv v = mpv_new(prec, count);
        
    if(! ok || ! mpv_read_points(v, 0, v->count, f, -1)) {
        mpv_free(v);
        MD5_Final(md5sum, &md5);
        
        return NULL;
    }
    
    ok = ok && mpv_update_md5(v, &md5);
    ok = ok && MD5_Final(md5sum, &md5);
    
    for (int i = 0; i < 16 && ok; i++) {
        ok = ok && fmd5sum[i] == md5sum[i];
    }
    
    if(! ok) {
        mpv_free(v);
        
        return NULL;
    }
    
    return v;
}

mpv mpv_read_partial(FILE *f, long pos, ulong start, ulong step, ulong count, bool complex) {
    mfile h = mfile_read_header(f, pos, MPV_FILE_ID, MPV_HEADER_LEN);
    if(h == NULL) {
        return NULL;
    }
    
    h->pos = 20; // skip the file ID, length of the file and the header length
    uint prec = mfile_geti(h);
    ulong totCount = mfile_getl(h);
    
    // no MD5 check, as this is a quick partial read
    mfile_free(h);
    
    mpv v = mpv_new(prec, count);
    
    if(! mpv_read_points_partial(v, totCount, f, -1, start, step, count, complex)) {
        mpv_free(v);
        
        return NULL;
    }
    
    return v;
}

static bool mpv_write_header(mpv v, FILE *f, long pos) {
    if(f == NULL) { // v already checked
        return false;
    }
        
    ulong flen = mpv_fileLen(v);
    mfile h = mfile_header(MPV_FILE_ID, flen, MPV_HEADER_LEN);
    if(h == NULL) {
        return false;
    }
    
    mfile_puti(h, v->prec);
    mfile_putl(h, v->count);
    
    // compute MD5 checksum
    MD5_CTX md5;
    byte md5sum[16];
    
    MD5_Init(&md5);
    bool ok = MD5_Update(&md5, h->data, (uint) h->len);
    
    ok = ok && mpv_update_md5(v, &md5);
    
    ok = ok && MD5_Final(md5sum, &md5);
    mfile_putbs(h, md5sum, 16);
    
    ok = ok && h->len == MPV_HEADER_LEN && mfile_write_to(f, pos, h, 0, h->len);
    
    mfile_free(h);
    
    return ok;
}

ulong mpv_write_to(mpv v, FILE *f, long pos) {
    if(v == NULL || v->count == 0 || v->vals == NULL || v->sexp == NULL) {
        return 0;
    }
    
    if(! mpv_write_header(v, f, pos)) {
        return 0;
    }
    
    ulong bytes = mpv_write_points(v, f, -1);
    
    return bytes == 0 ? 0 : MPV_HEADER_LEN + bytes;
}

bool mpv_read_points(mpv v, ulong vpos, ulong count, FILE *f, long pos) {
    if(v == NULL || f == NULL || v->count < vpos + count || v->vals == NULL || v->sexp == NULL ||
          (pos >= 0 && fseek(f, pos, SEEK_SET) != 0)) {
        return false;
    }
        
    int step = v->limbs;
    ulong cs = count * step;
    
    bool ok = fread_block(v->vals + vpos * step, 8, cs, FILE_BLOCK, f) == cs;
    ok = ok && fread_block(v->sexp + vpos, 4, count, FILE_BLOCK, f) == count;
    
    return ok;
}

bool mpv_read_points_partial(mpv v, ulong totCount, FILE *f, long pos,
                           ulong start, ulong step, ulong count, bool complex) {
    if(v == NULL || f == NULL || v->count == 0 || v->vals == NULL || v->sexp == NULL ||
          (pos >= 0 && fseek(f, pos, SEEK_SET) != 0)) {
        return false;
    }
    
    if(complex) {
        if(v->count != 2 * count || totCount < 2 * (start + (count - 1) * step + 1)) {
            return false;
        }
    } else if(v->count != count || totCount < start + (count - 1) * step + 1) {
        return false;
    }
        
    ulong fp = ftell(f);
    int limbs = v->limbs;
    int numSize = complex ? limbs << 4 : limbs << 3;
    
    ulong first = fp + numSize * start;
    ulong fstep = numSize * step;
    
    bool ok = true;
    for (ulong i = 0; ok && i < count; i++) {
        ulong ipos = first + i * fstep;
        ok = ok && fseek(f, ipos, SEEK_SET) == 0;
        ok = ok && fread(v->vals + i * numSize, 1, numSize, f) == numSize;
    }
    
    first = fp + numSize * totCount;
    int expc = complex ? 2 : 1;
    fstep = 4 * expc * step;
    for (ulong i = 0; ok && i < count; i++) {
        ulong ipos = first + i * fstep;
        ok = ok && fseek(f, ipos, SEEK_SET) == 0;
        ok = ok && fread(v->sexp + i * expc, 4, expc, f) == expc;
    }
    
    return ok;
}

ulong mpv_write_points(mpv v, FILE *f, long pos) {
    if(v == NULL || f == NULL || v->count == 0 || v->vals == NULL || v->sexp == NULL ||
       (pos >= 0 && ! io_file_seek(f, pos, true))) {
        return 0;
    }
    
    int limbs = v->limbs;
    
    ulong wr = fwrite(v->vals, 8, v->count * limbs, f);
    if(wr < v->count * limbs) {
        return 0;
    }
    
    ulong wrs = fwrite(v->sexp, 4, v->count, f);
    if(wrs < v->count) {
        return 0;
    }
    
    return 8 * wr + 4 * wrs;
}

bool mpv_write_csv(mpv v, char *fn, bool complex, int digits, ulong start, ulong count, bool append) {
    int cr = complex ? 2 : 1;
    if(fn == NULL || v == NULL || digits < 15 || count == 0 || start + count > v->count / cr) {
        return false;
    }
    
    uint maxLen = digits + 50;
    char *dx = malloc(sizeof(char) * maxLen);
    if(dx == NULL) {
        return false;
    }
    
    char *dy = NULL;
    if(complex) {
        dy = malloc(sizeof(char) * maxLen);
        if(dy == NULL) {
            free(dx);
            
            return false;
        }
    }
    
    FILE *f = fopen(fn, append ? "a" : "w");
    if(f == NULL) {
        free(dx);
        if(complex) {
            free(dy);
        }
        
        return false;
    }
    
    long end = start + count;
    long prec = v->prec;
    mpfr_t x, y;
    mpfr_init2(x, prec);
    mpfr_init2(y, prec);
    bool err = false;
    
    char rf[100];
    sprintf(rf, "%% .%dRg", digits);
    
    for (long i = start; i < end && ! err; i++) {
        mpv_get(x, v, i * cr);
        mpfr_snprintf(dx, maxLen, rf, x);
        
        if(complex) {
            mpv_get(y, v, i * cr + 1);
            mpfr_snprintf(dy, maxLen, rf, y);
            
            fprintf(f, "%s, %s\n", dx, dy);
        } else {
            fprintf(f, "%s\n", dx);
        }
    }
    
    fclose(f);
    
    free(dx);
    if(complex) {
        free(dy);
    }
    
    mpfr_clear(x);
    mpfr_clear(y);
    
    return ! err;
}

mpv mpv_read_csv(char *fn, int prec, bool complex, long max) {
    FILE *f = fopen(fn, "r");
    if(f == NULL) {
        return false;
    }
    
    mpfr_t x, y;
    mpfr_init2(x, prec);
    mpfr_init2(y, prec);
    
    int cr = complex ? 2 : 1;
    mpv v = mpv_new(prec, cr * max);
    if(v == NULL) {
        fclose(f);
        
        return NULL;
    }
    
    char sx[20051], sy[20051];
    long pos = 0;
    
    if(complex) {
        while(! feof(f) && pos < 2 * max) {
            if(fscanf(f, "%20050s%20050s", sx, sy) == 2) {
                long lx = strlen(sx);
                if(lx > 1 && sx[lx - 1] == ',') {
                    sx[lx - 1] = 0;
                }
                
                bool ok = 0 == mpfr_set_str(x, sx, 10, MPFR_RNDN);
                ok = ok && 0 == mpfr_set_str(y, sy, 10, MPFR_RNDN);
                
                if(ok) {
                    // set the values into the vector
                    mpv_set(v, pos++, x);
                    mpv_set(v, pos++, y);
                }
            }
        }
    } else {
        while(! feof(f) && pos < max) {
            if(fscanf(f, "%20050s", sx) == 1) {
                if(0 == mpfr_set_str(x, sx, 10, MPFR_RNDN)) {
                    // set the value into the vector
                    mpv_set(v, pos++, x);
                }
            }
        }
    }
    
    fclose(f);
    
    mpfr_clear(x);
    mpfr_clear(y);
    
    if(pos == 0 || (pos < cr * max && ! mpv_resize(v, pos))) {
        mpv_free(v);
        
        return NULL;
    }
    
    return v;
}

bool mpv_update_md5_partial(mpv vect, MD5_CTX *md5, ulong pos, ulong count) {
    if(vect == NULL || md5 == NULL || pos + count > vect->count) {
        return false;
    }
    
    if(count == 0) {
        return true;
    }
    
    uint limbs = vect->limbs;
    
    bool ok = true;
    ulong end = pos + count;
    for (ulong i = pos; ok && i < end; i++) {
        ok = ok && MD5_Update(md5, vect->vals + (i * limbs), limbs << 3);
    }
    
    for (ulong i = pos; ok && i < end; i++) {
        ok = ok && MD5_Update(md5, vect->sexp + i, 4);
    }
    
    return ok;
}

bool mpv_update_md5(mpv vect, MD5_CTX *md5) {
    return mpv_update_md5_partial(vect, md5, 0, vect->count);
}

ulong fwrite_mpfr(mpfr_t x, FILE *f) {
    if(f == NULL || x == NULL) {
        return 0;
    }
    
    ulong prec = x->_mpfr_prec;
    if(fwrite(&prec, 8, 1, f) != 1) {
        return 0;
    }
    
    int exp = mpfr_regular_p(x) ? (int) x->_mpfr_exp :
        mpfr_zero_p(x) ? MPV_ZERO_EXP :
        mpfr_nan_p(x) ? MPV_NAN_EXP : MPV_INF_EXP;
    exp = (exp << 1) | (x->_mpfr_sign < 0 ? 1 : 0);
    if(fwrite(&exp, 4, 1, f) != 1) {
        return 0;
    }
    
    ulong limbs = mpv_limbs(prec);
    if(fwrite(x->_mpfr_d, 8, limbs, f) != limbs) {
        return 0;
    }
    
    return 12 + 8 * limbs;
}

bool fread_mpfr(mpfr_t x, FILE *f) {
    if(f == NULL || x == NULL) {
        return false;
    }
    
    ulong prec;
    if(fread(&prec, 8, 1, f) != 1) {
        return false;
    }
    
    int exp;
    if(fread(&exp, 4, 1, f) != 1) {
        return false;
    }
    
    mpfr_init2(x, prec);
    
    x->_mpfr_sign = exp & 1 ? -1 : 1;
    exp >>= 1;
    
    x->_mpfr_exp = exp >= MPV_MIN_EXP && exp <= MPV_MAX_EXP ? exp :
        exp == MPV_ZERO_EXP ? __MPFR_EXP_ZERO :
        exp == MPV_INF_EXP ? __MPFR_EXP_INF : __MPFR_EXP_NAN;
    
    ulong limbs = mpv_limbs(prec);
    if(fread(x->_mpfr_d, 8, limbs, f) != limbs) {
        mpfr_clear(x);
        
        return false;
    }
    
    return true;
}

ulong mwrite_mpfr(mpfr_t x, mfile m) {
    if(m == NULL || x == NULL) {
        return 0;
    }
    
    ulong prec = x->_mpfr_prec;
    if(mwrite(&prec, 8, 1, m) != 1) {
        return 0;
    }
    
    int exp = mpfr_regular_p(x) ? (int) x->_mpfr_exp :
        mpfr_zero_p(x) ? MPV_ZERO_EXP :
        mpfr_nan_p(x) ? MPV_NAN_EXP : MPV_INF_EXP;
    exp = (exp << 1) | (x->_mpfr_sign < 0 ? 1 : 0);
    if(mwrite(&exp, 4, 1, m) != 1) {
        return 0;
    }
    
    ulong limbs = mpv_limbs(prec);
    if(mwrite(x->_mpfr_d, 8, limbs, m) != limbs) {
        return 0;
    }
    
    return 12 + 8 * limbs;
}

bool mread_mpfr(mpfr_t x, mfile m) {
    if(m == NULL || x == NULL) {
        return false;
    }
    
    ulong prec;
    if(mread(&prec, 8, 1, m) != 1) {
        return false;
    }
    
    int exp;
    if(mread(&exp, 4, 1, m) != 1) {
        return false;
    }
    
    mpfr_init2(x, prec);
    
    x->_mpfr_sign = exp & 1 ? -1 : 1;
    exp >>= 1;
    
    x->_mpfr_exp = exp >= MPV_MIN_EXP && exp <= MPV_MAX_EXP ? exp :
        exp == MPV_ZERO_EXP ? __MPFR_EXP_ZERO :
        exp == MPV_INF_EXP ? __MPFR_EXP_INF : __MPFR_EXP_NAN;
    
    ulong limbs = mpv_limbs(prec);
    if(mread(x->_mpfr_d, 8, limbs, m) != limbs) {
        mpfr_clear(x);
        
        return false;
    }
    
    return true;
}

bool dot_prod(mpfr_t p, mpv x, mpv y, uint prec) {
    if(p == NULL || x == NULL || y == NULL || x->count != y->count || prec < 32) {
        return false;
    }
    
    mpfr_t a, b, c, d;
    mpfr_inits2(prec, a, b, NULL);
    mpfr_set_zero(a, 1);
    
    ulong len = x->count;
    for (ulong i = 0; i < len; i++) {
        mpv_use(c, x, i);
        mpv_use(d, y, i);
        mpfr_mul(b, c, d, MPFR_RNDN);
        mpfr_add(a, a, b, MPFR_RNDN);
    }
    
    mpfr_set(p, a, MPFR_RNDN);
    
    mpfr_clears(a, b, NULL);
    
    return true;
}

bool dot_prodc(mpc p, mpv x, mpv y, uint prec) {
    if(p == NULL || x == NULL || y == NULL || x->count != y->count || x->count & 1 || prec < 32) {
        return false;
    }
    
    mpc a, b, c, d;
    mpc_inits(prec, a, b, NULL);
    mpc_set0(a);
    
    ulong len = x->count >> 1;
    for (ulong i = 0; i < len; i++) {
        mpv_usec(c, x, i);
        mpv_usec(d, y, i);
        mpc_mul(b, c, d);
        mpc_add(a, a, b);
    }
    
    mpc_set(p, a);
    
    mpc_clears(a, b, NULL);
    
    return true;
}

bool mpv_norm(mpfr_t x, mpv v, uint prec) {
    if(v == NULL || x == NULL || prec < 32) {
        return false;
    }
    
    mpfr_t a, s, b;
    mpfr_inits2(prec, a, s, NULL);
    mpfr_set_zero(s, 1);
    
    for (ulong i = 0; i < v->count; i++) {
        mpv_fuse(b, v, i);
        mpfr_sqr(a, b, MPFR_RNDU);
        mpfr_add(s, s, a, MPFR_RNDU);
    }
    
    mpfr_sqrt(x, s, MPFR_RNDU);
    mpfr_clears(a, s, NULL);
    
    return true;
}

bool mpv_norm2(mpfr_t x, mpv v, uint prec) {
    if(v == NULL || x == NULL || prec < 32) {
        return false;
    }
    
    mpfr_t a, s, b;
    mpfr_inits2(prec, a, s, NULL);
    mpfr_set_zero(s, 1);
    
    for (ulong i = 0; i < v->count; i++) {
        mpv_fuse(b, v, i);
        mpfr_sqr(a, b, MPFR_RNDU);
        mpfr_add(s, s, a, MPFR_RNDU);
    }
    
    mpfr_set(x, s, MPFR_RNDU);
    mpfr_clears(a, s, NULL);
    
    return true;
}

bool mpv_dist(mpfr_t x, mpv v, mpv u, uint prec) {
    if(v == NULL || u == NULL || v->count != u->count || x == NULL || prec < 32) {
        return false;
    }
    
    mpfr_t a, s, b, c;
    mpfr_inits2(prec, a, s, NULL);
    mpfr_set_zero(s, 1);
    
    for (ulong i = 0; i < v->count; i++) {
        mpv_fuse(b, v, i);
        mpv_fuse(c, u, i);
        mpfr_sub(a, b, c, MPFR_RNDA);
        mpfr_sqr(a, a, MPFR_RNDU);
        mpfr_add(s, s, a, MPFR_RNDU);
    }
    
    mpfr_sqrt(x, s, MPFR_RNDU);
    mpfr_clears(a, s, NULL);
    
    return true;
}

bool mpv_dist2(mpfr_t x, mpv v, mpv u, uint prec) {
    if(v == NULL || u == NULL || v->count != u->count || x == NULL || prec < 32) {
        return false;
    }
    
    mpfr_t a, s, b, c;
    mpfr_inits2(prec, a, s, NULL);
    mpfr_set_zero(s, 1);
    
    for (ulong i = 0; i < v->count; i++) {
        mpv_fuse(b, v, i);
        mpv_fuse(c, u, i);
        mpfr_sub(a, b, c, MPFR_RNDA);
        mpfr_sqr(a, a, MPFR_RNDU);
        mpfr_add(s, s, a, MPFR_RNDU);
    }
    
    mpfr_set(x, s, MPFR_RNDU);
    mpfr_clears(a, s, NULL);
    
    return true;
}

bool mpv_normalize(mpv v) {
    if(v == NULL || v->count == 0) {
        return false;
    }
    
    mpfr_t a;
    mpfr_init2(a, v->prec);
    
    if(! mpv_norm(a, v, v->prec) || mpfr_zero_p(a)) {
        return false;
    }
    
    mpfr_t x;
    for (ulong i = 0; i < v->count; i++) {
        mpv_fuse(x, v, i);
        mpfr_div(x, x, a, MPFR_RNDN);
        mpv_qsync(v, i, x);
    }
    
    mpfr_clear(a);
    
    return true;
}

bool mpv_print(mpv v, ulong pos, int digits) {
    if(v == NULL || digits < 2 || v->count <= pos) {
        return false;
    }
    
    char fmt[20];
    snprintf(fmt, 20, "%%.%dRg", digits);
    
    mpfr_t x;
    mpv_use(x, v, pos);
    
    mpfr_printf(fmt, x);
    
    return true;
}

bool mpv_printc(mpv v, ulong pos, int digits) {
    if(v == NULL || digits < 2 || v->count <= 2 * pos + 1) {
        return false;
    }
        
    char fmt[20];
    snprintf(fmt, 20, "(%%.%dRg, %%.%dRg)", digits, digits);
    
    mpc c;
    mpv_usec(c, v, pos);
    
    mpfr_printf(fmt, c->x, c->y);
    
    return true;
}

bool mpv_snprint(char *str, int len, mpv v, ulong pos, int digits) {
    if(v == NULL || digits < 2 || str == NULL || len < 2 * digits + 12 || v->count <= pos) {
        return false;
    }
    
    char fmt[20];
    snprintf(fmt, 20, "%%.%dRg", digits);
    
    mpfr_t x;
    mpv_use(x, v, pos);
    
    mpfr_snprintf(str, len, fmt, x);
    
    return true;
}

bool mpv_snprintc(char *str, int len, mpv v, ulong pos, int digits) {
    if(v == NULL || digits < 2 || v->count <= 2 * pos + 1 || len < 2 * digits + 12 || str == NULL) {
        return false;
    }
        
    char fmt[20];
    snprintf(fmt, 20, "(%%.%dRg, %%.%dRg)", digits, digits);
    
    mpc c;
    mpv_usec(c, v, pos);
    
    mpfr_snprintf(str, len, fmt, c->x, c->y);
    
    return true;
}

