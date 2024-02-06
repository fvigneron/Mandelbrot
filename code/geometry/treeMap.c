//
//  treeMap.c
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

#include "treeMap.h"


// MARK: creation and deletion of nodes

tmap tmap_new(int x, int y, char level, char maxLevel, char levelStep, char minLevel) {
    if(levelStep < 1 || levelStep > TMAP_MAX_LEVEL_STEP || maxLevel > TMAP_MAX_LEVEL || level > maxLevel) {
        return NULL;
    }
    
    int chCount = 1 << (levelStep << 1);
    
    tmap tm = malloc(sizeof(treeMap_struct) + chCount * sizeof(mapCount_union));
    
    tm->count = 0;
    tm->level = level;
    tm->maxLevel = maxLevel;
    tm->levelStep = levelStep;
    tm->minLevel = minLevel;
    tm->children = 0;
    tm->x = x;
    tm->y = y;
    
    if(level < (int) maxLevel) {
        for (int i = 0; i < chCount; i++) {
            tm->chld[i].ptr = NULL;
        }
    } else {
        for (int i = 0; i < chCount; i++) {
            tm->chld[i].count = 0;
        }
    }
    
    return tm;
}

tmap tmap_new_mandel(char maxLevel, char levelStep) {
    int ml = TMAP_MANDEL_MIN_LEVEL;
    
    return tmap_new(TMAP_MANDEL_X, TMAP_MANDEL_Y, ml, maxLevel, levelStep, ml);
}

tmap tmap_map(nset ps, char maxLevel, char levelStep) {
    if(ps == NULL) {
        return NULL;
    }
    
    tmap t = tmap_new(0, 0, -2, maxLevel, levelStep, -2);
    if(t == NULL) {
        return NULL;
    }
    
    u128 u;
    for (long i = 0; i < ps->count; i++) {
        nset_point(u, ps, i);
        tmap_add(t, u);
    }
    
    return t;
}

void tmap_free(tmap tm) {
    if(tm == NULL) {
        return;
    }
    
    if(! tmap_isLeaf(tm)) {
        int chCount = tmap_maxChildren(tm);
        tmap ch;
        mapCount_union *chld = tm->chld;
        for (int i = 0; i < chCount; i++) {
            ch = chld[i].ptr;
            
            if(ch != NULL) {
                tmap_free(ch);
            }
        }
    }
    
    free(tm);
}

tmap tmap_clone(tmap tree) {
    if(tree == NULL) {
        return NULL;
    }
    
    int maxCh = tmap_maxChildren(tree);
    ulong size = sizeof(treeMap_struct) + maxCh * sizeof(mapCount_union);
    tmap nt = malloc(size);
    
    *nt = *tree;
    
    mapCount_union *nch = nt->chld;
    mapCount_union *tch = tree->chld;
    tmap child;
    
    if(! tmap_isLeaf(tree)) {
        for (int i = 0; i < maxCh; i++) {
            child = tch[i].ptr;
            
            nch[i].ptr = child != NULL ? tmap_clone(child) : NULL;
        }
    } else {
        for (int i = 0; i < maxCh; i ++) {
            nch[i].count = tch[i].count;
        }
    }
    
    return nt;
}

bool tmap_new_child(tmap tm, int i, int j) {
    if(tm == NULL || i < 0 || j < 0 || tmap_isLeaf(tm)) {
        return false;
    }
    
    int cols = tmap_columns(tm);
    int ind = i * cols + j;
    if(i >= cols || j >= cols || tm->chld[ind].ptr != NULL) {
        return false;
    }
    
    byte step = tm->levelStep;
    int left = (tm->x << step) + i;
    int bot = (tm->y << step) + j;
    
    tm->chld[ind].ptr = tmap_new(left, bot, tm->level + step, tm->maxLevel, step, tm->minLevel);
    if(tm->chld[ind].ptr != NULL) {
        tm->children ++;
        
        return true;
    }
    
    return false;
}

bool tmap_add(tmap tree, u128 p) {
    if(p == NULL) {
        return false;
    }
    
    return tmap_addl(tree, p->xh, p->yh);
}

bool tmap_addl(tmap tree, ulong x, ulong y) {
    if(tree == NULL || tree->level < -2) {
        return false;
    }
    
    int lev = tree->level;
    byte sh = 62 - lev;
    
    if(lev > -2 && (x >> sh != tree->x || y >> sh != tree->y)) {
        return false;
    }
    
    char ls = tree->levelStep;
    int colsMask = (1 << ls) - 1;
    sh -= ls;
    int cx = (x >> sh) & colsMask;
    int cy = (y >> sh) & colsMask;
    int ind = cx * (colsMask + 1) + cy;
    
    mapCount_union *pt = tree->chld + ind;
    if(lev + ls > tree->maxLevel) { // no other children in this leaf
        if(pt->count == 0) {
            tree->children ++;
        }
        
        pt->count ++;
        tree->count ++;
        
        return true;
    } else {
        if(pt->ptr == NULL && ! tmap_new_child(tree, cx, cy)) {
            return false;
        }
        
        if(tmap_addl(pt->ptr, x, y)) {
            tree->count ++;
            
            return true;
        }
        
        return false;
    }
}

ulong tmap_adds(tmap tree, nset ps) {
    if(tree == NULL || tree->level < -2 || ps == NULL) {
        return 0;
    }
    
    u128 p;
    ulong pa = 0;
    for (long i = 0; i < ps->count; i++) {
        nset_point(p, ps, i);
        
        pa += tmap_add(tree, p) != 0 ? 1 : 0;
    }
    
    return pa;
}

bool tmap_union(tmap tree, tmap op) {
    if(tree == NULL || op == NULL || tree->level != op->level || tree->x != op->x ||
       tree->y != op->y || tree->maxLevel != op->maxLevel || tree->levelStep != op->levelStep) {
        return false;
    }
    
    tree->count += op->count;
    
    mapCount_union *tch = tree->chld;
    mapCount_union *och = op->chld;
    tmap child, chldo;
    bool ok = true, maxCh = tmap_maxChildren(tree);
    if(tmap_isLeaf(tree)) {
        for (int i = 0; i < maxCh; i++) {
            if(tch[i].count == 0 && och[i].count > 0) {
                tree->children ++;
            }

            tch[i].count += och[i].count;
        }
    } else {
        for (int i = 0; ok && i < maxCh; i++) {
            chldo = och[i].ptr;
            
            if(chldo == NULL) {
                continue;
            }
            
            child = tch[i].ptr;
            if(child == NULL) {
                tch[i].ptr = tmap_clone(chldo);
                tree->children ++;
            } else {
                ok = ok && tmap_union(child, chldo);
            }
        }
    }
    
    return ok;
}

bool tmap_reduce(tmap tree, int maxLevel) {
    if(tree == NULL || maxLevel < tree->level || maxLevel > tree->maxLevel || tmap_isLeaf(tree)) {
        return false;
    }
    
    if(maxLevel == tree->maxLevel) {
        return true;
    }
    
    int lev = tree->level;
    int st = tree->levelStep;
    bool leaf = lev + st > maxLevel;
    
    mapCount_union *tch = tree->chld;
    tmap child;
    bool ok = true, maxCh = tmap_maxChildren(tree);
    for (int i = 0; ok && i < maxCh; i++) {
        child = tch[i].ptr;
        
        if(child == NULL) {
            if(leaf) { // redundant on most systems as NULL == 0, here only for safety
                tch[i].count = 0;
            }
            
            continue;
        }
        
        // if this node becomes a leaf
        if(leaf) {
            tch[i].count = child->count;
            tmap_free(child);
            
            continue;
        }
        
        // otherwise delegate to children
        ok = ok && tmap_reduce(child, maxLevel);
    }
    
    if(ok) {
        tree->maxLevel = maxLevel;
    }
    
    return ok;
}

// MARK: comparisons and node counts

bool tmap_eq(tmap t1, tmap t2) {
    if(t1 == NULL || t2 == NULL || t1->count != t2->count || t1->level != t2->level ||
       t1->children != t2->children || t1->x != t2->x || t1->y != t2->y) {
        return false;
    }

    int eq = 1;
    if(tmap_isLeaf(t1)) {
        for (int i = 0; eq && i < tmap_maxChildren(t1); i++) {
            eq = eq && t1->chld[i].count == t2->chld[i].count;
        }
        
        return eq;
    }
    
    for (int i = 0; eq && i < tmap_maxChildren(t1); i++) {
        int nn1 = t1->chld[i].ptr != NULL;
        int nn2 = t2->chld[i].ptr != NULL;
        
        eq = eq && nn1 == nn2;
        
        if(! nn1) {
            continue;
        }
        
        eq = eq && tmap_eq(t1->chld[i].ptr, t2->chld[i].ptr);
    }
    
    return eq;
}

bool tmap_leq(tmap t1, tmap t2) {
    if(t1 == NULL || t2 == NULL || t1->count > t2->count || t1->level != t2->level ||
       t1->children > t2->children || t1->x != t2->x || t1->y != t2->y) {
        return false;
    }

    int leq = 1;
    if(tmap_isLeaf(t1)) {
        for (int i = 0; leq && i < tmap_maxChildren(t1); i++) {
            leq = leq && t1->chld[i].count <= t2->chld[i].count;
        }
        
        return leq;
    }
    
    for (int i = 0; leq && i < tmap_maxChildren(t1); i++) {
        int nn1 = t1->chld[i].ptr != NULL;
        int nn2 = t2->chld[i].ptr != NULL;
        
        leq = leq && nn1 <= nn2;
        
        if(! leq || ! nn1) {
            continue;
        }
        
        leq = leq && tmap_leq(t1->chld[i].ptr, t2->chld[i].ptr);
    }
    
    return leq;
}

ulong tmap_nodes(tmap tm) {
    if(tm == NULL) {
        return 0;
    }
    
    if(tm->level + 2 * tm->levelStep > tm->maxLevel) { // pre-leaf
        return tm->children + 1;
    } else {
        ulong count = 1;

        int maxCh = tmap_maxChildren(tm);
        int ech = 0, ok = 1;
        tmap child;
        mapCount_union *chld = tm->chld;
        for (int i = 0; ok && i < maxCh; i++) {
            child = chld[i].ptr;
            if(child == NULL) {
                continue;
            }
            
            ech ++;
            count += tmap_nodes(child);
        }
        
        if(ech != tm->children) {
            printf("treeMap has bad structure !\n");
            
            ok = 0;
        }
        
        return count;
    }
}

ulong tmap_memory(tmap tm) {
    ulong count = tmap_nodes(tm);
    ulong node = sizeof(treeMap_struct) + tmap_maxChildren(tm) * sizeof(mapCount_union);
    
    return count * node;
}

ulong tmap_count_abs(tmap tm, drect r) {
    if(r == NULL) {
        return 0;
    }
    
    drect_translate(r, 2, 0, 0);
    ulong c = tmap_count(tm, r);
    drect_translate(r, -2, 0, 0);
    
    return c;
}

ulong tmap_count(tmap tm, drect r) {
    if(tm == NULL || r == NULL || r->w <= 0 || r->h <= 0 || r->tpow > tm->maxLevel ||
       r->x < 0 || r->y < 0) {
        return 0;
    }
    
    int lv = tm->level;
    
    long tx = tm->x;
    long ty = tm->y;
    int level = r->tpow;
    long x = r->x, y = r->y, w = r->w, h = r->h;
    if(level <= lv) { // the tree is a [sub]-pixel of the rectangle, quick final decision
        int dl = lv - level;
        
        tx >>= dl;
        ty >>= dl;
        
        int in = tx >= x && ty >= y && tx < x + w && ty < y + h;
        
        return in ? tm->count : 0;
    } else { // the rectangle has finer pixels than the tree, study the intersection
        int dl = level - lv;
        tx <<= dl;
        ty <<= dl;
        int size = 1 << dl; // the tree is the rectangle [square]: (tx, ty, size, size) x 2^{-level}
        
        int in = tx >= x && ty >= y && tx + size <= x + w && ty + size <= y + h;
        if(in) {
            return tm->count;
        }
        
        int out = tx >= x + w || ty >= y + h || x >= tx + size || y >= ty + size;
        if(out) {
            return 0;
        }
        
        // if partial overlap, then ask the children
        ulong count = 0;
        
        int maxCh = tmap_maxChildren(tm);
        int ech = 0, ok = 1;
        tmap child;
        mapCount_union *chld = tm->chld;
        for (int i = 0; ok && i < maxCh; i++) {
            child = chld[i].ptr;
            if(child == NULL) {
                continue;
            }
            
            ech ++;
            count += tmap_count(child, r);
        }
        
        if(ech != tm->children) {
            printf("treeMap has bad structure !\n");
        }
        
        return count;
    }
}

// MARK: bitmap operations

bool tmap_add_counts(bmap bm, tmap tm) {
    if(tm == NULL || bm == NULL || bm->r.w <= 0 || bm->r.h <= 0 || bm->r.tpow > tm->maxLevel) {
        return false;
    }
    
    int lv = tm->level;
    
    long tx = tm->x;
    long ty = tm->y;
    int level = bm->r.tpow;
    long x = bm->r.x, y = bm->r.y, w = bm->r.w, h = bm->r.h;
    if(level <= lv) { // the tree is a [sub]-pixel of the rectangle, quick final decision
        int dl = lv - level;
        
        tx >>= dl;
        ty >>= dl;
        
        int in = tx >= x && ty >= y && tx < x + w && ty < y + h;
        if(in) {
            ulong ind = (tx - x) + w * (ty - y);
            
            bm->pix[ind] += tm->count;
        }
        
        return true;
    } else { // the rectangle has finer pixels than the tree, study the intersection
        int dl = level - lv;
        tx <<= dl;
        ty <<= dl;
        int size = 1 << dl; // the tree is the rectangle [square]: (tx, ty, size, size) x 2^{-level}
        
        int out = tx >= x + w || ty >= y + h || x >= tx + size || y >= ty + size;
        if(out) { // nothing to do
            return true;
        }
        
        // if partial overlap, then ask the children
        int maxCh = tmap_maxChildren(tm);
        int ech = 0, ok = 1;
        tmap child;
        mapCount_union *chld = tm->chld;
        for (int i = 0; ok && i < maxCh; i++) {
            child = chld[i].ptr;
            if(child == NULL) {
                continue;
            }
            
            ech ++;
            ok = ok && tmap_add_counts(bm, child);
        }
        
        if(ech != tm->children) {
            printf("treeMap has bad structure !\n");
            
            ok = 0;
        }
        
        return  ok;
    }
}

bool tmap_bounds(drect r, tmap tm) {
    if(r == NULL || tm == NULL) {
        return false;
    }
    
    dyadic_rect u = {0}, c;
    
    int maxCh = tmap_maxChildren(tm);
    int ech = 0;
    mapCount_union *chld = tm->chld;
    
    if(tmap_isLeaf(tm)) {
        ulong c;
        int lv = tm->level;
        int st = tm->levelStep;
        int tx = tm->x << st;
        int ty = tm->y << st;
        int colMask = (1 << st) - 1;
        
        for (int i = 0; i < maxCh; i++) {
            c = chld[i].count;
            if(c == 0) {
                continue;
            }
            
            if(ech == 0) {
                u.tpow = lv + st;
                u.w = 1;
                u.h = 1;
                u.x = tx + (i & colMask);
                u.y = ty + (i >> st);
            } else {
                drect_add(&u, tx + (i & colMask), ty + (i >> st));
            }
            
            ech ++;
        }
    } else {
        tmap child;
        for (int i = 0; i < maxCh; i++) {
            child = chld[i].ptr;
            if(child == NULL) {
                continue;
            }
            
            if(ech == 0) {
                tmap_bounds(&u, child);
            } else {
                tmap_bounds(&c, child);
                
                drect_union(&u, &c);
            }
            
            ech ++;
        }
    }
    
    *r = u;
    
    if(ech != tm->children) {
        printf("treeMap has bad structure !\n");
        
        return false;
    }
    
    return true;
}

bmap tmap_bitmap(tmap tm, drect r) {
    bmap bm = bmap_new(r);
    if(bm == NULL) {
        return NULL;
    }
    
    if(tmap_add_counts(bm, tm)) {
        return bm;
    } else {
        free(bm);
        
        return NULL;
    }
}

// MARK: IO operations

/// @brief Reads the root node from the header of a tmap memory file @c h.
///
/// @param h the memory file containing the header of the tmap file
/// @param bounds pointer to store the bounds of the tree map from the file, @c NULL to ignore
///
/// @return the root node of the tree, @c NULL is some error occurred
static tmap tree_from_header(mfile h, drect bounds) {
    if(h == NULL) {
        return NULL;
    }
    
    h->pos = 0;
    ulong id = mfile_getl(h); // already checked
    ulong fileLen = mfile_getl(h);
    int headLen = mfile_geti(h);
    
    if(headLen < TMAP_HEADER_LEN || fileLen <= TMAP_HEADER_LEN || id == 0) {
        return NULL;
    }
    
    long count = mfile_getl(h);
    
    char lv[4];
    ulong rd = mread(lv, 1, 4, h);
    
    int xy[2];
    rd += mread(xy, 4, 2, h);
    
    // read the bounds, if needed
    if(bounds != NULL) {
        bounds->x = mfile_getl(h);
        bounds->y = mfile_getl(h);
        bounds->w = mfile_geti(h);
        bounds->h = mfile_geti(h);
        bounds->tpow = mfile_geti(h);
    }
    
    if(count == 0 || rd < 6 || lv[0] != lv[3]) { // not root
        return NULL;
    }
    
    tmap tm = tmap_new(xy[0], xy[1], lv[0], lv[1], lv[2], lv[3]);
    if(tm != NULL) {
        tm->count = count;
                
        if(! mseek(h, headLen)) {
            tmap_free(tm);
            
            tm = NULL;
        }
    }
    
    return tm;
}

/// @brief Reads a packed unsigned long from a memory file @c m.
///
/// Big endian encoding is used for packing reasons, slightly cumbersome.
/// Not a big performance penalty, as most packs are 1 or 2 bytes long.
///
/// @param m the memory file
static inline ulong read_packed_count(mfile m) {
    byte *mb = m->data + m->pos;
    byte t = mb[0];

    ulong c = 0;
    byte *buf = (byte *) &c;
    int len = 0;
    if((t & TMAP_MAX_CHAR) != 0) {
        c = t & (TMAP_MAX_CHAR - 1);
        
        len = 1;
    } else if((t & (TMAP_MAX_SHORT >> 8)) != 0) {
        buf[1] = t;
        buf[0] = mb[1];
        c &= (TMAP_MAX_SHORT - 1);
        
        len = 2;
    } else if((t & (TMAP_MAX_INT >> 24)) != 0) {
        buf[3] = t;
        buf[2] = mb[1];
        buf[1] = mb[2];
        buf[0] = mb[3];
        c &= (TMAP_MAX_INT - 1);
        
        len = 4;
    } else { // very, very rarely, if ever, we will be here
        buf[7] = t;     // looks ugly, but faster than a loop :)
        buf[6] = mb[1];
        buf[5] = mb[2];
        buf[4] = mb[3];
        buf[3] = mb[4];
        buf[2] = mb[5];
        buf[1] = mb[6];
        buf[0] = mb[7];
        
        len = 8;
    }
    
    if(m->pos + len <= m->len) {
        m->pos += len;
        
        return c;
    }
    
    return 0;
}

/// @brief Reads the children of the given @c parent from the memory file @c m.
///
/// @param parent the parent node
/// @param px the x position in the parent
/// @param py the y position in the parent
/// @param m the memory file
///
/// @return @ref true if successfull, @ref false otherwise
static bool read_tree(tmap parent, int px, int py, mfile m) {
    int n = tmap_bitDataLen(parent);
    if(mfile_space(m) < n) {
        return false;
    }
    
    byte *buf = m->data + m->pos;
    m->pos += n;
    
    int st = parent->levelStep;
    int ch = 1 << (st << 1);
    
    // create the new tree
    int cols = 1 << st;
    int ind = px * cols + py;
        
    tmap tree = malloc(sizeof(treeMap_struct) + ch * sizeof(mapCount_union));
    tree->level = parent->level + st;
    tree->maxLevel = parent->maxLevel;
    tree->levelStep = st;
    tree->minLevel = parent->minLevel;
    tree->children = 0;
    tree->x = (parent->x << st) + px;
    tree->y = (parent->y << st) + py;
    tree->count = 0;
    
    bool ok = true;
    mapCount_union *chld = tree->chld;
    if(tmap_isLeaf(tree)) {
        for (int i = 0; ok && i < ch; i++) {
            if(((buf[i >> 3] >> (i & 7)) & 1) != 0) {
                ulong c = read_packed_count(m);
                chld[i].count = c;
                tree->children ++;
                
                ok = c > 0;
                tree->count += c;
            } else {
                chld[i].count = 0;
            }
        }
    } else {
        int col = cols - 1;
        
        for (int i = 0; ok && i < ch; i++) {
            if(((buf[i >> 3] >> (i & 7)) & 1) != 0) {
                ok = ok && read_tree(tree, i >> st, i & col, m);
            } else {
                chld[i].ptr = NULL;
            }
        }
    }
    
    parent->chld[ind].ptr = tree;
    parent->children ++;
    parent->count += tree->count;
    
    return ok;
}

/// @brief Reads the tree with the given @c root from the memory file @c m.
///
/// @param root the root of the tree
/// @param m the memory file
///
/// @return @ref true if successfull, @ref false otherwise
static inline bool read_root(tmap root, mfile m) {
    if(root == NULL || m == NULL) {
        return false;
    }
    
    long count = root->count;
    root->count = 0;
  
    int n = tmap_bitDataLen(root);
    byte buf[n];
    
    bool ok = n == mread(buf, 1, n, m);
    
    int st = root->levelStep;
    int ch = 1 << (st << 1);
    int col = (1 << st) - 1;
  
    for (int i = 0; ok && i < ch; i++) {
        if(((buf[i / 8] >> (i & 7)) & 1) != 0) {
            ok = ok && read_tree(root, i >> st, i & col, m);
        } else {
            root->chld[i].ptr = NULL;
        }
    }
    
    if(! ok || count != root->count) {
        tmap_free(root);
        
        return false;
    }
    
    return ok;
}

tmap tmap_mread(mfile m, bool checkMD5, drect bounds) {
    tmap tree = tree_from_header(m, bounds);
    bool ok = read_root(tree, m); // reads the entire tree
    
    if(ok && checkMD5) {
        long flen = *(ulong *) (m->data + 8);
        uint hlen = *(uint *) (m->data + 16);
        
        if(hlen >= TMAP_HEADER_LEN + 16) {
            byte *fmd5sum = m->data + (hlen - 16); // MD5 checksum at the end of the header

            // compute MD5 checksum
            MD5_CTX md5;
            byte md5sum[16];
            
            MD5_Init(&md5);
            MD5_Update(&md5, m->data, hlen - 16); // the header except the MD5

            if(flen - hlen < 1L << 32) {
                MD5_Update(&md5, m->data + hlen, (uint) (flen - hlen)); // the rest of the file
            } else { // every 512 bytes
                long i = hlen;
                for (; i <= flen - 512; i += 512) {
                    MD5_Update(&md5, m->data + i, 512);
                }
                
                if(i < flen) {
                    MD5_Update(&md5, m->data + i, (uint) (flen - i));
                }
            }
            
            MD5_Final(md5sum, &md5);
            
            for (int i = 0; i < 16 && ok; i++) {
                ok = ok && fmd5sum[i] == md5sum[i];
            }
        } else {
            ok = false;
        }
    }
    
    if(! ok && tree != NULL) {
        tmap_free(tree);
    }
    
    return ok ? tree : NULL;
}

tmap tmap_read(FILE *f, long pos, bool checkMD5, drect bounds) {
    mfile m = mfile_read_id(f, pos, TMAP_FILE_ID);
    tmap t = tmap_mread(m, checkMD5, bounds);
    mfile_free(m);
    
    return t;
}

tmap tmap_load(char *fileName, ulong pos, bool checkMD5, drect bounds) {
    mfile m = mfile_load_id(fileName, pos, TMAP_FILE_ID);
    tmap t = tmap_mread(m, checkMD5, bounds);
    mfile_free(m);
    
    return t;
}

/// @brief Writes the structure of the node @c tree (that is not a leaf) to the memory file @c m.
///
/// @param m the memory file
/// @param tree the node
///
/// @return @ref true if successfull, @ref false otherwise
static inline bool bit_data(mfile m, tmap tree) {
    int n = tmap_bitDataLen(tree);
    if(! mfile_ensure_cap(m, n)) {
        return false;
    }
    
    byte *buf = m->data + m->pos;
    m->pos += n;
    m->len = m->pos > m->len ? m->pos : m->len;
    
    int mc = tmap_maxChildren(tree) - 1;
    int ch = mc >> 3; // that is, here mcm == 8 * ch
    
    byte b = 0;
    tmap *bufChld = (tmap *) tree->chld + (ch << 3);
    for (int j = mc & 7; j >= 0; j--) {
        b <<= 1;
        b |= bufChld[j] == NULL ? 0 : 1;
    }
    buf[ch] = b;
    
    for (ch --; ch >= 0; ch --) {
        bufChld = (tmap *) tree->chld + (ch << 3);
        
        b = 0;
        for (int j = mc & 7; j >= 0; j--) {
            b <<= 1;
            b |= bufChld[j] == NULL ? 0 : 1;
        }
        buf[ch] = b;
    }
    
    return true;
}

/// @brief Writes the content of the leaf @c tree to the memory file @c m.
///
/// @param m the memory file
/// @param tree the leaf
///
/// @return @ref true if successfull, @ref false otherwise
static inline bool leaf_data(mfile m, tmap tree) {
    int maxCh = tmap_maxChildren(tree);
    if(! mfile_ensure_cap(m, maxCh * 8)) {
        return false;
    }
    
    int bp = 0;
    
    byte *buf = m->data + m->pos;
    ulong c;
    byte *bc = (byte *) &c;
    int ech = 0;
    mapCount_union *chld = tree->chld;
    for (int i = 0; i < maxCh; i++) {
        c = chld[i].count;
        if(c == 0) {
            continue;
        }
        
        ech ++;
        
        if(c < TMAP_MAX_CHAR) {
            c |= TMAP_MAX_CHAR;
            
            buf[bp ++] = bc[0];
        } else if(c < TMAP_MAX_SHORT) {
            c |= TMAP_MAX_SHORT;
            
            buf[bp ++] = bc[1];
            buf[bp ++] = bc[0];
        } else if(c < TMAP_MAX_INT) {
            c |= TMAP_MAX_INT;
            
            buf[bp ++] = bc[3];
            buf[bp ++] = bc[2];
            buf[bp ++] = bc[1];
            buf[bp ++] = bc[0];
        } else {
            buf[bp ++] = bc[7];
            buf[bp ++] = bc[6];
            buf[bp ++] = bc[5];
            buf[bp ++] = bc[4];
            buf[bp ++] = bc[3];
            buf[bp ++] = bc[2];
            buf[bp ++] = bc[1];
            buf[bp ++] = bc[0];
        }
    }
    
    m->pos += bp;
    m->len = m->pos > m->len ? m->pos : m->len;
    
    if(ech != tree->children) {
        printf("treeMap has bad structure !\n");
        
        return false;
    }
    
    return true;
}

/// @brief Writes the content of the tree into the memory file @c m at the current position.
///
/// @param tree the tree
/// @param m the memory file
///
/// @return @ref true if successfull, @ref false otherwise
static bool write_tree(tmap tree, mfile m) {
    bool ok = bit_data(m, tree);
    
    if(tmap_isLeaf(tree)) {
        ok = ok && leaf_data(m, tree);
    } else {
        int ech = 0;
        int maxCh = tmap_maxChildren(tree);
        mapCount_union *chld = tree->chld;
        tmap child;
        for (int i = 0; ok && i < maxCh; i++) {
            child = chld[i].ptr;
            if(child == NULL) {
                continue;
            }
            
            ech ++;
            
            ok = ok && write_tree(child, m);
        }
        
        if(ech != tree->children) {
            printf("treeMap has bad structure !\n");
            
            return false;
        }
    }
    
    return ok;
}

/// @brief Writes the content of the tree into the memory file @c m at the current position.
///
/// @param tree the root of the tree
/// @param m the memory file
///
/// @return @ref true if successfull, @ref false otherwise
static inline bool write_root(tmap tree, mfile m) {
    if(tree == NULL || m == NULL) {
        return false;
    }
    
    bool ok = write_tree(tree, m);
    
    return ok;
}

static long fileStartPos;

/// @brief Writes the header of the @c tree into the memory file @c m.
///
/// @param tree the tree
/// @param m the memory file
/// @param bounds the precomputed bounds, @c NULL to recompute
///
/// @return @ref true if successfull, @ref false otherwise
static bool write_header(tmap tree, mfile m, drect bounds) {
    if(tree == NULL || m == NULL) {
        return false;
    }
    
    fileStartPos = m->pos;
    uint hlen = TMAP_HEADER_LEN + 16; // MD5 checksum
    mfile_write_header(TMAP_FILE_ID, hlen, hlen, m);
    
    mfile_putl(m, tree->count);
    
    char lv[4] = {tree->level, tree->maxLevel, tree->levelStep, tree->minLevel};
    mfile_putbs(m, (byte *) lv, 4);
    
    int xy[2] = {tree->x, tree->y};
    mfile_putis(m, (uint *) xy, 2);
        
    // bounding rect
    dyadic_rect r;
    if(bounds == NULL) {
        tmap_bounds(&r, tree);
        
        drect_translate(&r, -2, 0, 0);
    } else {
        r = *bounds;
    }
    
    mfile_putl(m, r.x);
    mfile_putl(m, r.y);
    mfile_puti(m, r.w);
    mfile_puti(m, r.h);
    mfile_puti(m, r.tpow);
    
    // space for the MD5 checksum
    mfile_putl(m, 0);
    mfile_putl(m, 0);
    
    return true;
}

/// @brief Computes the MD5 sum of the file and write it to the end of its header.
///
/// @param m the memory file of type tmap.
///
/// @return @ref true if successfull, @ref false otherwise
static bool complete_header(mfile m) {
    if(m == NULL) {
        return false;
    }
    
    ulong tmLen = m->pos - fileStartPos;
    
    bool ok = mseek(m, fileStartPos + 8);
    ok = ok && mfile_putl(m, tmLen);

    // compute MD5 checksum
    MD5_CTX md5;
    
    long flen = m->len;
    uint hlen = *(uint *) (m->data + 16);
    byte *md5sum = m->data + (hlen - 16);

    MD5_Init(&md5);
    MD5_Update(&md5, m->data, hlen - 16); // the header except the MD5

    if(flen - hlen < 1L << 32) {
        MD5_Update(&md5, m->data + hlen, (uint) (flen - hlen)); // the rest of the file
    } else { // every 512 bytes
        long i = hlen;
        for (; i <= flen - 512; i += 512) {
            MD5_Update(&md5, m->data + i, 512);
        }
        
        if(i < flen) {
            MD5_Update(&md5, m->data + i, (uint) (flen - i));
        }
    }
    
    MD5_Final(md5sum, &md5);
    
    return ok;
}

mfile tmap_write(tmap tree, drect bounds) {
    if(tree == NULL || tree->level != tree->minLevel) {
        return false;
    }
    
    mfile m = mfile_new(500000);
    
    bool ok = write_header(tree, m, bounds);
    ok = ok && write_root(tree, m);
    ok = ok && complete_header(m);
    
    if(ok) {
        return m;
    }
    
    mfile_free(m);
    
    return NULL;
}

bool tmap_save(tmap tree, char *fileName) {
    mfile m = tmap_write(tree, NULL);
    bool ok = mfile_save(fileName, m, false);
    mfile_free(m);
    
    return ok;
}
