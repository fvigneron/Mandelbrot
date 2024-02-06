//
//  memFile.c
//  Mandel_v0.6
//
//  Created by MIHALACHE Nicolae on 1/19/21.
//  Copyright © 2021 MIHALACHE Nicolae. All rights reserved.
//

#include "io.h"
#include "memFile.h"

mfile mfile_new(ulong cap) {
    mfile m = malloc(sizeof(memFile_struct));
    m->data = malloc(cap);
    
    m->pos = 0;
    m->len = 0;
    m->cap = cap;
    
    return m;
}

bool mfile_free(mfile m) {
    if(m == NULL || m->data == NULL) {
        return false;
    }
    
    free(m->data);
    m->data = NULL;
    m->cap = 0;
    
    free(m);
    
    return true;
}

mfile mfile_read(FILE *f) {
    if(f == NULL) {
        return NULL;
    }
    
    ulong pos = ftell(f);
    if(fseek(f, 0, SEEK_END) != 0) {
        return NULL;
    }
    
    ulong len = ftell(f) - pos;
    mfile m = mfile_new(len);
    
    if(! io_file_seek(f, pos, false) || fread(m->data, 1, len, f) != len) {
        mfile_free(m);
        
        return NULL;
    }
    
    m->len = len;
    
    return m;
}

mfile mfile_read_from(FILE *f, ulong pos, ulong len) {
    if(f == NULL) {
        return NULL;
    }
    
    if(fseek(f, 0, SEEK_END) != 0) {
        return NULL;
    }
    
    ulong end = ftell(f);
    if(end < pos + len) {
        return NULL;
    }
    
    mfile m = mfile_new(len);
    
    if(! io_file_seek(f, pos, false) || fread(m->data, 1, len, f) != len) {
        mfile_free(m);
        
        return NULL;
    }
    
    m->len = len;
    
    return m;
}

mfile mfile_load(char *fileName) {
    FILE *f = fopen(fileName, "r");
    mfile m = mfile_read(f);
    fclose(f);
    
    return m;
}

bool mfile_write(FILE *f, mfile m) {
    if(f == NULL || m == NULL) {
        return false;
    }
    
    return fwrite(m->data, 1, m->len, f) == m->len;
}

bool mfile_write_to(FILE *f, long fpos, mfile m, ulong mpos, ulong len) {
    if(f == NULL || m == NULL || mpos + len > m->len ||
       (fpos >= 0 && ! io_file_seek(f, fpos, true))) {
        return false;
    }
    
    return fwrite(m->data + mpos, 1, len, f) == len;
}

bool mfile_save(char *fileName, mfile m, bool append) {
    if(fileName == NULL || m == NULL) {
        return false;
    }
    
    FILE *f = fopen(fileName, append ? "a" : "w");
    if(f == NULL) {
        return false;
    }
    
    bool ok = fwrite(m->data, 1, m->len, f) == m->len;
    
    return fclose(f) == 0 && ok;
}

inline bool mseek(mfile m, ulong pos) {
    if(m == NULL || pos > m->len) {
        return false;
    }
    
    m->pos = pos;
    
    return true;
}

ulong mread(void *buf, ulong block, ulong count, mfile m) {
    ulong len = block * count;
    if(m == NULL || buf == NULL || m->pos + len > m->len) {
        return false;
    }
    
    memcpy(buf, m->data + m->pos, len);
    m->pos += len;
    
    return count;
}

char *mfile_get_str(mfile m) {
    if(m == NULL || m->pos >= m->len) {
        return NULL;
    }
    
    ulong p = m->pos, len = 0;
    while(p + len < m->len && m->data[p + len] != 0) {
        len ++;
    }
    
    if(p + len >= m->len) {
        return NULL;
    }
    
    m->pos += len + 1;
    
    return (char *) m->data + p;
}

inline bool mfile_ensure_cap(mfile m, ulong extra) {
    if(m == NULL) {
        return false;
    }
    
    if(m->pos + extra > m->cap) {
        ulong ns = m->cap * MFILE_FACTOR + MFILE_ADD;
        ulong nd = m->pos + extra;
        ns = ns > nd ? ns : nd;
        
        m->data = realloc(m->data, ns);
        if(m->data == NULL) {
            m->cap = 0;
            
            return false;
        }
        
        m->cap = ns;
    }
    
    return true;
}

ulong mwrite(void *buf, ulong block, ulong count, mfile m) {
    ulong len = block * count;
    if(m == NULL || buf == NULL) {
        return 0;
    }
    
    if(! mfile_ensure_cap(m, len)) {
        return 0;
    }
    
    memcpy(m->data + m->pos, buf, len);
    m->pos += len;
    m->len = m->pos > m->len ? m->pos : m->len;
    
    return count;
}

inline bool mfile_putb(mfile m, byte c) {
    if(m == NULL || m->pos > m->len || (m->pos + 1 > m->len && ! mfile_ensure_cap(m, 1))) {
        return false;
    }
    
    m->data[m->pos ++] = c;
    m->len = m->pos > m->len ? m->pos : m->len;
    
    return true;
}

inline bool mfile_putw(mfile m, word w) {
    if(m == NULL || m->pos > m->len || (m->pos + 2 > m->len && ! mfile_ensure_cap(m, 2))) {
        return false;
    }
    
    *(word *) (m->data + m->pos) = w;
    m->pos += 2;
    m->len = m->pos > m->len ? m->pos : m->len;
    
    return true;
}

inline bool mfile_puti(mfile m, uint i) {
    if(m == NULL || m->pos > m->len || (m->pos + 4 > m->len && ! mfile_ensure_cap(m, 4))) {
        return false;
    }
    
    *(uint *) (m->data + m->pos) = i;
    m->pos += 4;
    m->len = m->pos > m->len ? m->pos : m->len;
    
    return true;
}

inline bool mfile_putl(mfile m, ulong l) {
    if(m == NULL || m->pos > m->len || (m->pos + 8 > m->len && ! mfile_ensure_cap(m, 8))) {
        return false;
    }
    
    *(ulong *) (m->data + m->pos) = l;
    m->pos += 8;
    m->len = m->pos > m->len ? m->pos : m->len;
    
    return true;
}

inline bool mfile_putp(mfile m, void *p) {
    if(m == NULL || m->pos > m->len || (m->pos + 8 > m->len && ! mfile_ensure_cap(m, 8))) {
        return false;
    }
    
    *(void **) (m->data + m->pos) = p;
    m->pos += 8;
    m->len = m->pos > m->len ? m->pos : m->len;
    
    return true;
}

bool mfile_putd(mfile m, double d) {
    if(m == NULL || m->pos > m->len || (m->pos + 8 > m->len && ! mfile_ensure_cap(m, 8))) {
        return false;
    }
    
    *(double *) (m->data + m->pos) = d;
    m->pos += 8;
    m->len = m->pos > m->len ? m->pos : m->len;
    
    return true;
}

bool mfile_putld(mfile m, ldbl ld) {
    if(m == NULL || m->pos > m->len || (m->pos + 8 > m->len && ! mfile_ensure_cap(m, 8))) {
        return false;
    }
    
    memcpy(m->data + m->pos, &ld, 10);
    m->pos += 10;
    m->len = m->pos > m->len ? m->pos : m->len;
    
    return true;
}

inline bool mfile_putbs(mfile m, byte c[], ulong len) {
    return mwrite(c, 1, len, m) == len;
}

inline bool mfile_putws(mfile m, word w[], ulong len) {
    return mwrite(w, 2, len, m) == len;
}

inline bool mfile_putis(mfile m, uint i[], ulong len) {
    return mwrite(i, 4, len, m) == len;
}

inline bool mfile_putls(mfile m, ulong l[], ulong len) {
    return mwrite(l, 8, len, m) == len;
}

inline byte mfile_getb(mfile m) {
    if(m == NULL || m->pos >= m->len) {
        return 0;
    }
    
    return m->data[m->pos ++];
}

inline word mfile_getw(mfile m) {
    if(m == NULL || m->pos + 2 > m->len) {
        return 0;
    }
    
    short w = *(word *) (m->data + m->pos);
    m->pos += 2;
    
    return w;
}

inline uint mfile_geti(mfile m) {
    if(m == NULL || m->pos + 4 > m->len) {
        return 0;
    }
    
    uint w = *(uint *) (m->data + m->pos);
    m->pos += 4;
    
    return w;
}

inline ulong mfile_getl(mfile m) {
    if(m == NULL || m->pos + 8 > m->len) {
        return 0;
    }
    
    ulong w = *(ulong *) (m->data + m->pos);
    m->pos += 8;
    
    return w;
}

inline void *mfile_getp(mfile m) {
    if(m == NULL || m->pos + 8 > m->len) {
        return NULL;
    }
    
    void *p = *(void **) (m->data + m->pos);
    m->pos += 8;
    
    return p;
}

double mfile_getd(mfile m) {
    if(m == NULL || m->pos + 8 > m->len) {
        return NAN;
    }
    
    double w = *(double *) (m->data + m->pos);
    m->pos += 8;
    
    return w;
}

ldbl mfile_getld(mfile m) {
    if(m == NULL || m->pos + 10 > m->len) {
        return NAN;
    }
    
    ldbl w;
    memcpy(&w, m->data + m->pos, 10);
    m->pos += 10;
    
    return w;
}

inline bool mfile_getbs(mfile m, byte c[], ulong len) {
    return mread(c, 1, len, m) == len;
}

inline bool mfile_getws(mfile m, word w[], ulong len) {
    return mread(w, 2, len, m) == len;
}

inline bool mfile_getis(mfile m, uint i[], ulong len) {
    return mread(i, 4, len, m) == len;
}

inline bool mfile_getls(mfile m, ulong l[], ulong len) {
    return mread(l, 8, len, m) == len;
}

mfile mfile_header(char *textId, ulong fileLen, uint headerLen) {
    if(fileLen <= 20 || textId == NULL || headerLen < 20) {
        return NULL;
    }
    
    mfile m = mfile_new(headerLen);
    
    fileHeaderId id = {.fileTypeId = 0};
    snprintf(id.type, 9, "%s", textId);
    
    mfile_putl(m, id.fileTypeId);
    mfile_putl(m, fileLen);
    mfile_puti(m, headerLen);
    
    return m;
}

bool mfile_header_id(mfile h, char *textId) {
    if(h == NULL || h->len < 20 || textId == NULL) {
        return false;
    }
    
    fileHeaderId id = {.fileTypeId = 0};
    snprintf(id.type, 9, "%s", textId);
    
    ulong w = *(ulong *) (h->data);
    
    return w == id.fileTypeId;
}

bool mfile_header_add_md5(mfile h, void *data, ulong dataLen) {
    if(h == NULL || h->len < 20 || data == NULL || dataLen == 0) {
        return false;
    }
    
    MD5_CTX md5;
    byte md5sum[16];
    
    MD5_Init(&md5);
    MD5_Update(&md5, h->data, (uint) h->len);
    
    long i = 0;
    for (; i <= ((long) dataLen) - 512; i += 512) {
        MD5_Update(&md5, data + i, 512);
    }
    
    if(i < dataLen) {
        MD5_Update(&md5, data + i, (uint) (dataLen - i));
    }
    
    MD5_Final(md5sum, &md5);
    mfile_putbs(h, md5sum, 16);
    
    return true;
}

bool mfile_header_check_md5(mfile h, void *data, ulong dataLen) {
    if(h == NULL || h->len < 20 || data == NULL || dataLen == 0) {
        return false;
    }
    
    MD5_CTX md5;
    byte md5sum[16];
    
    MD5_Init(&md5);
    MD5_Update(&md5, h->data, (uint) h->len - 16);
    
    long i = 0;
    for (; i <= ((long) dataLen) - 512; i += 512) {
        MD5_Update(&md5, data + i, 512);
    }
    
    if(i < dataLen) {
        MD5_Update(&md5, data + i, (uint) (dataLen - i));
    }
    
    MD5_Final(md5sum, &md5);
    
    for (int j = 0; j < 16; j++) {
        if(md5sum[j] != h->data[h->len - 16 + j]) {
            return false;
        }
    }
    
    return true;
}

mfile mfile_from_header(char *textId, ulong fileLen, uint headerLen) {
    if(fileLen <= 20 || textId == NULL || headerLen < 20) {
        return NULL;
    }
    
    mfile m = mfile_new(fileLen);
    
    fileHeaderId id = {.fileTypeId = 0};
    snprintf(id.type, 9, "%s", textId);
    
    mfile_putl(m, id.fileTypeId);
    mfile_putl(m, fileLen);
    mfile_puti(m, headerLen);
    
    return m;
}

mfile mfile_read_header(FILE *f, long pos, char *textId, int hlen) {
    if(f == NULL || textId == NULL || (pos >= 0 && ! io_file_seek(f, pos, false))) {
        return NULL;
    }
    
    long fp = ftell(f);
    if(fp < 0) {
        return NULL;
    }
    
    fileHeaderId id = {.fileTypeId = 0};
    snprintf(id.type, 9, "%s", textId);
    
    ulong fid = 0, fileLen = 0;
    bool ok = fread(&fid, 8, 1, f) == 1 && fid == id.fileTypeId;
    ok = ok && fread(&fileLen, 8, 1, f) == 1;
    
    // check the header length
    uint headerLen = 0;
    ok = ok && fread(&headerLen, 4, 1, f) == 1;
    ok = ok && headerLen >= hlen;
    
    // rewind and read if all OK
    mfile m = ok ? mfile_read_from(f, fp, headerLen) : NULL;
        
    return m;
}

mfile mfile_read_unknown_header(FILE *f, long pos) {
    if(f == NULL || (pos >= 0 && ! io_file_seek(f, pos, false))) {
        return NULL;
    }
    
    long fp = ftell(f);
    if(fp < 0) {
        return NULL;
    }
        
    ulong fid = 0, fileLen = 0;
    bool ok = fread(&fid, 8, 1, f) == 1;
    ok = ok && fread(&fileLen, 8, 1, f) == 1;
    
    // check the header length
    uint headerLen = 0;
    ok = ok && fread(&headerLen, 4, 1, f) == 1;
    ok = ok && headerLen <= fileLen;
    
    // rewind and read if all OK
    ok = ok && io_file_seek(f, fp, false);
    mfile m = ok ? mfile_read_from(f, fp, headerLen) : NULL;
        
    return m;
}

mfile mfile_load_header(char *fileName, ulong pos, char *textId, int hlen) {
    if(fileName == NULL || textId == NULL) {
        return NULL;
    }
    
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return NULL;
    }
    
    mfile m = mfile_read_header(f, pos, textId, hlen);
    
    fclose(f);
    
    return m;
}

mfile mfile_load_unknown_header(char *fileName) {
    if(fileName == NULL) {
        return NULL;
    }
    
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return NULL;
    }
    
    mfile m = mfile_read_unknown_header(f, 0);
    
    fclose(f);
    
    return m;
}

bool mfile_write_header(char *textId, ulong fileLen, uint headerLen, mfile m) {
    if(fileLen <= 20 || textId == NULL || headerLen < 20 || m == NULL) {
        return false;
    }
        
    fileHeaderId id = {.fileTypeId = 0};
    snprintf(id.type, 9, "%s", textId);
    
    mfile_putl(m, id.fileTypeId);
    mfile_putl(m, fileLen);
    mfile_puti(m, headerLen);
    
    return true;
}

mfile mfile_load_id(char *fileName, ulong pos, char *textId) {
    if(fileName == NULL || textId == NULL) {
        return NULL;
    }
    
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return NULL;
    }
    
    if(pos != 0 && io_file_seek(f, pos, false)) {
        fclose(f);
        
        return NULL;
    }
        
    mfile m = mfile_read_id(f, pos, textId);
    
    fclose(f);
    
    return m;
}

mfile mfile_read_id(FILE *f, long pos, char *textId) {
    mfile m = mfile_read_header(f, pos, textId, 0);
    if(m == NULL) {
        return NULL;
    }
    
    m->pos = 8;
    ulong flen = mfile_getl(m);
    uint hlen = mfile_geti(m);
    ulong dlen = flen - hlen;
    
    m->pos = m->len; // mfile_ensureCap from current position, not on top of the length
    
    bool ok = hlen == m->len && mfile_ensure_cap(m, dlen);
    ok = ok && fread(m->data + hlen, 1, dlen, f) == dlen;
    
    if(ok) {
        m->len += dlen;
        
        return m;
    }
    
    mfile_free(m);
    
    return NULL;
}
