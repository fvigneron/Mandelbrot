//
//  bitmap64.c
//  Mandel_v0.6
//
//  Created by MIHALACHE Nicolae on 1/23/21.
//  Copyright © 2021 MIHALACHE Nicolae. All rights reserved.
//

#include "bitmap.h"
#include "memFile.h"


bmap bmap_new(drect r) {
    if(r == NULL || r->w == 0 || r->h == 0) {
        return NULL;
    }
    
    ulong wh = r->w; // w, h are uint
    wh *= r->h;
    
    ulong size = (wh << 3) + sizeof(bitMap);
    bmap bm = malloc(size);
    
    bm->type = BMAP_TYPE_GENERIC;
    bm->subType = BMAP_SUB_TYPE_GENERIC;
    bm->r = *r;
    
    bm->pixelType = BMAP_PIXEL_TYPE_LINEAR;
    bm->sgn = 0;
    bm->typeBits = 0;
    bm->zeroTransp = 0;
    bm->colorMap = BMAP_COLOR_MAP_LINEAR;
    
    bm->a = 1;
    bm->b = 0;
    bm->power = 1;
    bm->mapA = 1;
    bm->mapB = 0;
    
    bm->colorLow = BMAP_COLOR_WHITE;
    bm->colorFirst = BMAP_COLOR_WHITE;
    bm->colorHigh = BMAP_COLOR_BLACK;
    bm->colorLast = BMAP_COLOR_BLACK;
    
    bm->useHD = false;
    
    ulong *d = bm->pix;
    for (ulong i = 0; i < wh; i ++) {
        d[i] = 0;
    }
    
    return bm;
}

bool bmap_fill(bmap m, drect r, ulong pix) {
    if(m == NULL || r == NULL || r->tpow != m->r.tpow || ! drect_contains_rect(&m->r, r)) {
        return false;
    }
    
    uint sy = (uint) (r->y - m->r.y), ey = sy + r->h;
    uint sx = (uint) (r->x - m->r.x), ex = sx + r->w;
    for (uint y = sy; y < ey; y++) {
        ulong *l = m->pix + (y * (ulong) m->r.w);
        
        for (uint x = sx; x < ex; x++) {
            l[x] = pix;
        }
    }
    
    return true;
}

bool bmap_set_pixel(bmap bm, ulong val, int x, int y) {
    if(bm == NULL || x < 0 || y < 0 || x >= bm->r.w || y >= bm->r.h) {
        return false;
    }
    
    long ind = y;
    ind *= bm->r.w;
    ind += x;
    
    bm->pix[ind] = val;
    
    return true;
}

inline bool bmap_increment_pixel(bmap bm, int x, int y) {
    if(bm == NULL || x < 0 || y < 0 || x >= bm->r.w || y >= bm->r.h) {
        return false;
    }
    
    long ind = y;
    ind *= bm->r.w;
    ind += x;
    
    bm->pix[ind] ++;
    
    return true;
}

inline bool bmap_increment_coords(bmap bm, fp80 p) {
    if(bm == NULL || p == NULL || ! drect_contains80(&bm->r, p)) {
        return false;
    }
    
    long x = drect_abs_to_rel_x80(&bm->r, p->x);
    long y = drect_abs_to_rel_y80(&bm->r, p->y);
    
    return bmap_increment_pixel(bm, (int) x, (int) y);
}

bool bmap_filter(bmap bm, ulong min, ulong max, ulong val) {
    if(bm == NULL || max < min) {
        return false;
    }
    
    ulong wh = bm->r.w; // w, h are uint
    wh *= bm->r.h;
    
    ulong *d = bm->pix, v;
    for (ulong i = 0; i < wh; i ++) {
        v = d[i];
        d[i] = v >= min && v <= max ? val : v;
    }
    
    return true;
}

bool bmap_norm_lin(bmap bm, ulong maxVal) {
    if(bm == NULL || maxVal == 0) {
        return false;
    }
    
    ulong wh = bm->r.w; // w, h are uint
    wh *= bm->r.h;
    
    ulong *d = bm->pix;
    ulong max = 0, v;
    for (ulong i = 0; i < wh; i ++) {
        v = d[i];
        max = v > max ? v : max;
    }
    
    if(max == 0) {
        return true;
    }
    
    ldbl t;
    for (ulong i = 0; i < wh; i ++) {
        t = d[i];
        t *= maxVal;
        t /= max;
        
        d[i] = t;
    }
    
    bm->a = 1;
    bm->a /= maxVal - 1;
    bm->b = 0;
    
    return true;
}

bool bmap_norm_pow(bmap bm, ulong maxVal, ldbl pow, ulong min, ulong floor) {
    if(bm == NULL || maxVal == 0 || pow < 0) {
        return false;
    }
    
    ulong wh = bm->r.w; // w, h are uint
    wh *= bm->r.h;
    
    ulong *d = bm->pix;
    ulong max = 0, v;
    for (ulong i = 0; i < wh; i ++) {
        v = d[i];
        max = v > max ? v : max;
    }
    
    if(max <= min) {
        return false;
    }
    
    ldbl minp = powl(min, pow);
    ldbl t, c = 1 / (powl(max, pow) - minp);
    for (ulong i = 0; i < wh; i ++) {        
        v = d[i];
        if(v <= min) {
            d[i] = v == 0 ? 0 : floor;
            
            continue;
        }
        
        t = (powl(v, pow) - minp) * c; // maps to [0, 1]
        t *= maxVal - floor;
        
        d[i] = t + floor;
    }
    
    bm->a = 1;
    bm->a /= maxVal - 1;
    bm->b = 0;
    
    return true;
}

bool bmap_norm_log(bmap bm, ulong maxVal) {
    if(bm == NULL || maxVal == 0 || pow < 0) {
        return false;
    }
    
    ulong wh = bm->r.w; // w, h are uint
    wh *= bm->r.h;
    
    ulong *d = bm->pix;
    ulong max = 0, v;
    for (ulong i = 0; i < wh; i ++) {
        v = d[i];
        max = v > max ? v : max;
    }
    
    if(max <= 1) {
        return true;
    }
    
    ldbl t, c = 1 / logl(max);
    for (ulong i = 0; i < wh; i ++) {
        t = logl(d[i]) * c;
        t *= maxVal - 1;
        
        d[i] = t + 1;
    }
    
    bm->a = 1;
    bm->a /= maxVal - 1;
    bm->b = 0;
    
    return true;
}

bool bmap_normalize_column(bmap bm, int x, ulong maxVal) {
    if(bm == NULL || maxVal == 0 || x < 0 || x >= bm->r.w) {
        return false;
    }
    
    int w = bm->r.w;
    ulong wh = w; // w, h are uint
    wh *= bm->r.h;
    
    ulong *d = bm->pix;
    ulong max = 0, v;
    for (ulong i = x; i < wh; i += w) {
        v = d[i];
        max = v > max ? v : max;
    }
    
    if(max == 0) {
        return true;
    }
    
    ldbl t;
    for (ulong i = x; i < wh; i += w) {
        t = d[i];
        t *= maxVal - 1;
        t /= max;
        
        d[i] = t;
    }
    
    return true;
}

bool bmap_normalize_line(bmap bm, int y, ulong maxVal) {
    if(bm == NULL || maxVal == 0 || y < 0 || y >= bm->r.h) {
        return false;
    }
    
    ulong w = bm->r.w;
    ulong st = w * y;
    
    ulong *d = bm->pix;
    ulong max = 0, v;
    for (ulong i = 0; i < w; i ++) {
        v = d[st + i];
        max = v > max ? v : max;
    }
    
    if(max == 0) {
        return true;
    }
    
    ldbl t;
    for (ulong i = 0; i < w; i ++) {
        t = d[st + i];
        t *= maxVal - 1;
        t /= max;
        
        d[st + i] = t;
    }
    
    return true;
}

ulong bmap_get_pixel(bmap bm, int x, int y) {
    if(bm == NULL || x < 0 || y < 0 || x >= bm->r.w || y >= bm->r.h) {
        return 0;
    }
    
    long ind = y;
    ind *= bm->r.w;
    ind += x;
    
    return bm->pix[ind];
}

bool bmap_draw(bmap dst, bmap src) {
    if(dst == NULL || src == NULL) {
        return false;
    }
    
    if(dst == src) {
        return true; // nothing to do
    }
    
    dyadic_rect r;
    if(! drect_intersection(&r, &src->r, &dst->r)) {
        return false;
    }
    
    ulong *srcl = src->pix + (r.x - src->r.x + src->r.w * (r.y - src->r.y));
    ulong *dstl = dst->pix + (r.x - dst->r.x + dst->r.w * (r.y - dst->r.y));
    ulong len = r.w << 3;
    for (int dy = 0; dy < r.h; dy++) {
        memcpy(dstl + dy * dst->r.w, srcl + dy * src->r.w, len);
    }
    
    return true;
}

ulong *bmap_histogram(bmap bm, int boxes) {
    if(bm == NULL || boxes < 2) {
        return NULL;
    }
    
    ulong *h = malloc(sizeof(ulong) * boxes);
    for (int i = 0; i < boxes; i++) {
        h[i] = 0;
    }
    
    ldbl a = bm->a * boxes, b = bm->b * boxes;
    ulong *p = bm->pix;
    ulong c = bm->r.w;
    c *= bm->r.h;
    
    int box, max = boxes - 1;
    ldbl v;
    for (ulong i = 0; i < c; i++) {
        v = a * p[i] + b;
        box = v <= 0 ? 0 : v >= max ? max : (int) v;
        
        h[box] ++;
    }
    
    return h;
}

ulong bmap_count(bmap bm, ulong v) {
    if(bm == NULL) {
        return 0;
    }
    
    ulong h = 0;
    ulong *p = bm->pix;
    ulong c = bm->r.w;
    c *= bm->r.h;    
    
    for (ulong i = 0; i < c; i++) {
        h += p[i] == v ? 1 : 0;
    }
    
    return h;
}

bmap bmap_read(char *fileName, bool ignoreMD5) {
    mfile h = mfile_load_unknown_header(fileName);
    if(h == NULL || mfile_header_id(h, BMAP_FILE_ID)) {
        return NULL;
    }
    
    h->pos = 16; // ignore id and total file length
    uint hlen = mfile_geti(h);
    
    byte type = mfile_getb(h);
    byte subType = mfile_getb(h);
    
    dyadic_rect r = {mfile_getl(h), mfile_getl(h), mfile_geti(h), mfile_geti(h), mfile_geti(h)};
    
    bmap m = bmap_new(&r);
    if(m == NULL) {
        return NULL;
    }
    
    m->type = type;
    m->subType = subType;
    
    m->pixelType = mfile_getb(h);
    m->sgn = mfile_getb(h);
    m->typeBits = mfile_getb(h);
    m->zeroTransp = mfile_getb(h);
    m->colorMap = mfile_getb(h);
    
    m->a = mfile_getd(h);
    m->b = mfile_getd(h);
    m->power = mfile_getd(h);
    m->mapA = mfile_getd(h);
    m->mapB = mfile_getd(h);
    
    m->colorLow = mfile_geti(h);
    m->colorFirst = mfile_geti(h);
    m->colorLast = mfile_geti(h);
    m->colorHigh = mfile_geti(h);
    
    m->useHD = mfile_getb(h);
    bool hasMD5 = mfile_getb(h);
    
    if(m->useHD) {
        char *dx = mfile_get_str(h);
        char *dy = mfile_get_str(h);
        
        if(dx == NULL || dy == NULL) {
            printf("Could not retrieve the delta from %s !\n", fileName);
            free(m);
            mfile_free(h);
            
            return NULL;
        }
        
        uint precx = strlen(dx) * 3.3;
        uint precy = strlen(dy) * 3.3;
        
        uint prec = precx > precy ? precx : precy;
        mpfr_inits2(prec, m->hdx, m->hdy, NULL);
        
        if(0 != mpfr_set_str(m->hdx, dx, 10, MPFR_RNDN) ||
           0 != mpfr_set_str(m->hdy, dy, 10, MPFR_RNDN)) {
            printf("Could not retrieve the delta from %s !\n", fileName);
            bmap_free(m);
            mfile_free(h);
            
            return NULL;
        }
    }
    
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        bmap_free(m);
        mfile_free(h);
        
        return NULL;
    }
    
    ulong wh = r.w; // w, h are uint
    wh *= r.h;
    
    if(fseek(f, hlen, SEEK_SET) != 0 || fread(m->pix, 8, wh, f) != wh) {
        bmap_free(m);
        mfile_free(h);
        fclose(f);
        
        return NULL;
    }
    
    fclose(f);
    
    if(hasMD5 && ! ignoreMD5) {
        if(! mfile_header_check_md5(h, m->pix, wh * 8)) {
            bmap_free(m);
            m = NULL;
        }
    }
    
    mfile_free(h);
    
    return m;
}

void bmap_free(bmap bitmap) {
    if(bitmap == NULL) {
        return;
    }
    
    if(bitmap->useHD) {
        mpfr_clears(bitmap->hdx, bitmap->hdy, NULL);
    }
    
    free(bitmap);
}

bool bmap_write(bmap bm, char *fileName, bool md5) {
    if(fileName == NULL || bm == NULL) {
        return false;
    }
    
    ulong wh = bm->r.w; // w, h are uint
    wh *= bm->r.h;
    
    uint hlen = BMAP_HEADER_MIN_LEN + (md5 ? 16 : 0);
    
    char *dx = NULL, *dy = NULL;
    uint plx = 0, ply = 0;
    if(bm->useHD) {
        uint lenx = (uint) (bm->hdx->_mpfr_prec / 3 + 25);
        dx = malloc(lenx);
        char fmt[20];
        snprintf(fmt, 20, "%%.%dRg", (uint) (bm->hdx->_mpfr_prec / 3));
        plx = mpfr_snprintf(dx, lenx, fmt, bm->hdx);
        
        uint leny = (uint) (bm->hdy->_mpfr_prec / 3 + 25);
        dy = malloc(leny);
        snprintf(fmt, 20, "%%.%dRg", (uint) (bm->hdy->_mpfr_prec / 3));
        ply = mpfr_snprintf(dy, leny, fmt, bm->hdy);
        
        plx ++;
        ply ++;
        
        if(plx == 0 || plx >= lenx || ply == 0 || ply >= leny) {
            printf("Could not write high precision delta to %s !!!\n", fileName);
            free(dx);
            free(dy);
            
            return false;
        }
        
        hlen += plx + ply;
    }
    
    mfile h = mfile_header(BMAP_FILE_ID, hlen + wh * 8, hlen);
    if(h == NULL) {
        return false;
    }
    
    mfile_putb(h, bm->type);
    mfile_putb(h, bm->subType);
    
    mfile_putl(h, bm->r.x);
    mfile_putl(h, bm->r.y);
    mfile_puti(h, bm->r.w);
    mfile_puti(h, bm->r.h);
    mfile_puti(h, bm->r.tpow);
    
    mfile_putb(h, bm->pixelType);
    mfile_putb(h, bm->sgn);
    mfile_putb(h, bm->typeBits);
    mfile_putb(h, bm->zeroTransp);
    mfile_putb(h, bm->colorMap);
    
    mfile_putd(h, bm->a);
    mfile_putd(h, bm->b);
    mfile_putd(h, bm->power);
    mfile_putd(h, bm->mapA);
    mfile_putd(h, bm->mapB);
    
    mfile_puti(h, bm->colorLow);
    mfile_puti(h, bm->colorFirst);
    mfile_puti(h, bm->colorLast);
    mfile_puti(h, bm->colorHigh);
    
    mfile_putb(h, bm->useHD);
    mfile_putb(h, md5);
    
    if(bm->useHD) {
        mfile_putbs(h, (byte *) dx, plx);
        mfile_putbs(h, (byte *) dy, ply);
        
        free(dx);
        free(dy);
    }
    
    bool ok = true;
    if(md5 && ! mfile_header_add_md5(h, bm->pix, wh * 8)) {
        printf("Could not compute the MD5 checksum for %s !!!\n", fileName);
        
        ok = false;
    }
    
    ok = ok && mfile_save(fileName, h, false);
    mfile_free(h);
    
    if(! ok) {
        return false;
    }
    
    FILE *f = fopen(fileName, "a");
    ok = f != NULL && fwrite(bm->pix, 8, wh, f) == wh;
    fclose(f);
    
    return ok;
}

static FILE* bmap_write_bmp_header(drect r, char *fileName) {
    if(fileName == NULL || r == NULL) {
        return NULL;
    }
    
    ulong line = ((((ulong) r->w) * 3 + 3) / 4) * 4;
    ulong fs = 54 + r->h * line;
    if(r->w == 0 || r->h == 0 || line >> 32 || fs >> 32) { // file too  large
        return NULL;
    }
    
    byte *bfs = (byte *) &fs;
    byte *bw = (byte *) &r->w;
    byte *bh = (byte *) &r->h;
    
    byte bmpFileHeader[14] = {'B','M', bfs[0], bfs[1], bfs[2], bfs[3], 0, 0, 0, 0, 54, 0, 0, 0};
    byte bmpInfoHeader[40] = {40, 0, 0, 0, bw[0], bw[1], bw[2], bw[3], bh[0], bh[1], bh[2], bh[3],
        1, 0, 24, 0};
    
    for (int i = 18; i < 40; i++) {
        bmpInfoHeader[i] = 0;
    }
    
    FILE *f = fopen(fileName, "w");
    bool ok = f != NULL && fwrite(bmpFileHeader, 1, 14, f) == 14;
    ok = ok && fwrite(bmpInfoHeader, 1, 40, f) == 40;
    
    if(ok) {
        return f;
    }
    
    fclose(f);
    
    return NULL;
}

bool bmp_write(drect r, char *fileName, rgb_pix colors, void *context) {
    if(fileName == NULL || r == NULL || colors == NULL) {
        return false;
    }
    
    FILE *f = bmap_write_bmp_header(r, fileName);
    if(f == NULL) {
        return false;
    }
    
    ulong line = ((((ulong) r->w) * 3 + 3) / 4) * 4;
    byte rgb[line];
    for (uint i = (uint) line - 3; i < line; i++) {
        rgb[i] = 0;
    }
    uint irgb;
    byte *brgb = (byte *) &irgb;
    
    bool ok = true;
    for (ulong y = 0; y < r->h && ok; y++) {
        for (uint p = 0, x = 0; x < r->w; x++) {
            irgb = colors(x, (int) y, context);
            
            rgb[p ++] = brgb[0];
            rgb[p ++] = brgb[1];
            rgb[p ++] = brgb[2];
        }
        
        ok = ok && fwrite(rgb, 1, line, f) == line;
    }
    
    fclose(f);
    
    return ok;
}

bool bmap_write_bmp(bmap bm, char *fileName, rgb_map colors) {
    if(fileName == NULL || bm == NULL || colors == NULL) {
        return false;
    }
    
    FILE *f = bmap_write_bmp_header(&bm->r, fileName);
    if(f == NULL) {
        return false;
    }
    
    drect r = &bm->r;
    ulong line = ((((ulong) r->w) * 3 + 3) / 4) * 4;
    byte rgb[line];
    for (uint i = (uint) line - 3; i < line; i++) {
        rgb[i] = 0;
    }
    uint irgb;
    byte *brgb = (byte *) &irgb;
    
    bool ok = true;
    for (ulong y = 0; y < r->h && ok; y++) {
        uint p = 0;
        ulong st = y * r->w;
        
        for (uint x = 0; x < r->w; x++) {
            irgb = colors(bm, bm->pix[st + x]);
            
            rgb[p ++] = brgb[0];
            rgb[p ++] = brgb[1];
            rgb[p ++] = brgb[2];
        }
        
        ok = ok && fwrite(rgb, 1, line, f) == line;
    }
    
    fclose(f);
    
    return ok;
}

bool bmap_write_bw_bmp(bmap bm, char *fileName, bool inv) {
    if(fileName == NULL || bm == NULL) {
        return false;
    }
    
    FILE *f = bmap_write_bmp_header(&bm->r, fileName);
    if(f == NULL) {
        return false;
    }
    
    drect r = &bm->r;
    ulong line = ((((ulong) r->w) * 3 + 3) / 4) * 4;
    byte rgb[line];
    for (uint i = (uint) line - 3; i < line; i++) {
        rgb[i] = 0;
    }
    
    bool ok = true;
    for (ulong y = 0; y < r->h && ok; y++) {
        uint p = 0;
        ulong st = y * r->w;
        
        for (uint x = 0; x < r->w; x++) {
            double c = bm->pix[st + x] * bm->a + bm->b;
            c = inv ? 1 - c : c;
            byte col = c <= 0 ? 255 : c >= 1 ? 0 : 255 * (1 - c);
            
            rgb[p ++] = col;
            rgb[p ++] = col;
            rgb[p ++] = col;
        }
        
        ok = ok && fwrite(rgb, 1, line, f) == line;
    }
    
    fclose(f);
    
    return ok;
}
