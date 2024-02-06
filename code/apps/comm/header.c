//
//  header.c
//  Mandel_v0.6
//
//  Created by MIHALACHE Nicolae on 2/20/21.
//  Copyright © 2021 MIHALACHE Nicolae. All rights reserved.
//

#include <stdio.h>
#include <sys/timeb.h>
#include "stopWatch.h"

#include "header.h"
#include "memFile.h"
#include "nSet.h"
#include "treeMap.h"
#include "levCurve.h"
#include "misCurve.h"
#include "bitmap.h"
#include "io.h"
#include "mpv.h"

static bool check_md5(mfile h, FILE *f) {
    h->pos = 8;
    long fsize = mfile_getl(h);
    int hsize = mfile_geti(h);
    long dsize = fsize - hsize;
    
    // compute MD5 checksum
    MD5_CTX md5;
    byte md5sum[16];
    
    MD5_Init(&md5);
    uint pos = (uint) h->len - 16;
    bool ok = MD5_Update(&md5, h->data, pos);
    
    byte *fmd5sum = h->data + pos;
    
    ok = ok && io_update_md5(&md5, f, dsize);
    ok = ok && MD5_Final(md5sum, &md5);
    
    for (int i = 0; i < 16 && ok; i++) {
        ok = ok && fmd5sum[i] == md5sum[i];
    }
    
    return ok;
}

static bool check_content(char *fn, bool minis, bool checkMD5) {
    printf(checkMD5 ? minis ? "Checking the MD5 sums of the mini-files.\n" :
           "Checking the MD5 sum of the file.\n" : "Checking mini-files.\n");
    if(checkMD5) {
        io_print_start();
    }
    
    char time[100];
    struct timeb ts;
    ftime(&ts);
    
    FILE *f = fopen(fn, "r");
    if(f == NULL) {
        return false;
    }
    
    long size = io_file_size(f);
    
    long pos = 0;
    int find = 0;
    mfile m = mfile_read_unknown_header(f, pos);
    
    if(m == NULL) {
        printf("Could not re-read the first header !\n");
        fclose(f);
        
        return false;
    }
    
    long fid = mfile_getl(m);
    bool uniType = true;
    
    long fmfsize = mfile_getl(m);
    bool uniSize = true;
    
    find ++;
    pos = fmfsize;
    
    bool md5ok = true;
    if(checkMD5) {
        md5ok = check_md5(m, f);
        printf(".");
        fflush(stdout);
    }
    
    mfile_free(m);
    
    int count = 1;
    
    while(pos < size) {
        m = mfile_read_unknown_header(f, pos);
        
        if(m == NULL) {
            printf("Could not read the header with index %d.\n", find + 1);
            fclose(f);
            
            return false;
        }
        
        long mfid = mfile_getl(m);
        uniType = uniType && (fid == mfid);
        
        long mfsize = mfile_getl(m);
        uniSize = uniSize && (fmfsize == mfsize);
        
        find ++;
        pos += mfsize;
        
        if(mfid == 0 || mfsize <= HEADER_MIN_LEN || mfsize >= (1L << 50)) {
            printf("The header with index %d is malformed: type = %016lX, size = %ld.\n",
                   find + 1, mfid, mfsize);
            fclose(f);

            return false;
        }
        
        if(checkMD5) {
            bool ok = check_md5(m, f);
            md5ok = md5ok && ok;
            count ++;
            
            printf(ok ? "." : "X");
            printf(count % 1000 == 0 ? "\n\n" : count % 100 == 0 ? "\n" : count % 50 == 0 ?
                   "  " : count % 10 == 0 ? " " : "");
            fflush(stdout);
        }
        
        mfile_free(m);
    }
    
    fclose(f);
    
    if(checkMD5) {
        if(md5ok) {
            printf("\nAll MD5 checksums are correct.\n\n");
        } else {
            printf("\nSome MD5 checksum have failed !!!\n\n");
        }
    }
    
    if(pos > size) {
        printf("At least one header of a mini-file is corrupted !\n");
        
        return false;
    }
    
    char *msgt = uniType ? "the same type" : "different types";
    char *msgs = uniType ? "the same size" : "different sizes";
    printf("The file contains %d mini-files of %s (and of %s).\n", find, msgt, msgs);
    
    lapse(&ts, time);
    printf("\nAll tasks completed in %s.\n\n", time);
    fflush(stdout);
    
    return true;
}

static bool header(char *fn, bool checkMD5) {
    printf("Summary info on the binary file %s:\n\n", fn);
    
    mfile h = mfile_load_unknown_header(fn);
    if(h == NULL) {
        printf("The file %s is not accessible.\n", fn);
        
        return false;
    }
    
    fileHeaderId idNset, idTmap, idLevc, idBmap, idLevm, idMpv;
    snprintf(idNset.type, 9, NSET_FILE_ID);
    snprintf(idTmap.type, 9, TMAP_FILE_ID);
    snprintf(idLevc.type, 9, LEVC_FILE_ID);
    snprintf(idBmap.type, 9, BMAP_FILE_ID);
    snprintf(idLevm.type, 9, LEVM_FILE_ID);
    snprintf(idMpv.type, 9, MPV_FILE_ID);
    fileHeaderId ids[] = {idNset, idTmap, idLevc, idBmap, idLevm, idMpv};
    
    ulong id = mfile_getl(h);
    ulong fs = mfile_getl(h);
    uint hl = mfile_geti(h);
    
    long fl = io_size_of(fn);
    
    char size[30];
    file_size(fs, size);
    
    if(fl == fs) {
        printf("The file %s is of size %s (%ld bytes) and has a header of length %u bytes.\n",
               fn, size, fs, hl);
    } else {
        char rsize[30];
        file_size(fl, rsize);
        
        printf("The file %s is of size %s (%ld bytes) but reported as %s (%ld bytes) in the header, which is of of length %u bytes.\n",
               fn, rsize, fl, size, fs, hl);
        
        if(fs > fl) {
            printf("The file is corrupted !\n\n");
        }
    }
    
    int type = -1;
    for (int i = 0; i < 6 && type < 0; i++) {
        if(id == ids[i].fileTypeId) {
            type = i;
        }
    }
    
    ulong count, rcount, epsl, sta, ena;
    int lev, maxLev, levStep, minLev, ix, iy;
    int pp, per, a2p, tp;
    mpfr_t r, eps;
    long x, y, dx, dy;
    uint prec, wd, ht;
    
    switch(type) {
        case 0: // nset
            count = mfile_getl(h);
            rcount = mfile_getl(h);
            epsl = mfile_getl(h);
            
            printf("\nThe file contains an nSet with %ld points (%ld of them real), separated by %.3Lg.\n",
                   count, rcount, ldexpl(epsl, -126));
            
            break;
            
        case 1: // tmap
            count = mfile_getl(h);
            lev = (char) mfile_getb(h);
            maxLev = (char) mfile_getb(h);
            levStep = (char) mfile_getb(h);
            minLev = (char) mfile_getb(h);
            
            ix = mfile_geti(h);
            iy = mfile_geti(h);
            
            printf("\nThe file contains a treeMap with %ld points, from level %d to %d, with step %d.\n",
                   count, minLev, maxLev, levStep);
            
            if(lev != minLev) {
                printf("The level %d of the root is not the minumum level %d, an error !", lev, minLev);
            }
            
            printf("The root dyadic square is of side length %Lg, the bottom left corner is (%Lg, %Lg) in abstract (positive) coordinates.\n",
                   ldexpl(1, -lev), ldexpl(ix, -lev), ldexpl(iy, -lev));
            
            
            if(h->len >= TMAP_HEADER_LEN + 16) { // has bounding rect
                dyadic_rect r;
                mread(&r, 28, 1, h);
                
                char str[80];
                drect_print(&r, str, 80);
                
                printf("The bounding rectangle in the plane is %s.\n", str);
            }
            break;
            
        case 2: // levc
            per = mfile_geti(h);
            a2p = mfile_geti(h);
            prec = mfile_geti(h);
            
            sta = mfile_getl(h);
            ena = mfile_getl(h);
            
            mpfr_init2(r, prec);
            mpfr_init2(eps, prec);
            
            mread_mpfr(r, h);
            mread_mpfr(eps, h);
            
            bool proven = mfile_getb(h);
            double guard = mfile_getd(h);
            
            ulong segs = ena - sta;
            
            printf("\nThe file contains a levCurve of period %d with %ld segments, from unit root %ld to root %ld of order 2^%d.\n",
                   per, segs, sta, ena, a2p);
            if(proven) {
                mpfr_printf("The points are stored with %d bits of precision, the radius is %Rg, the proven max error is %Rg.\n",
                        prec, r, eps);
            } else {
                mpfr_printf("The points are stored with %d bits of precision, the radius is %Rg, the max error is %Rg, the univalent guard is %g.\n",
                prec, r, eps, guard);
            }
            
            mpfr_clear(eps);
            mpfr_clear(r);
            
            break;
            
        case 3: // bm64
            x = mfile_getl(h);
            y = mfile_getl(h);
            wd = mfile_geti(h);
            ht = mfile_geti(h);
            tp = mfile_geti(h);
            dx = mfile_getl(h);
            dy = mfile_getl(h);
            
            int tb = mfile_getb(h);
            double a = mfile_getd(h);
            double b = mfile_getd(h);
            bool sgn = mfile_getb(h);
            
            printf("\nThe file contains a 64-bit bitmap of the rectangle (%ld, %ld, %u, %u) / 2^%d.\n",
                   x + dx, y + dy, wd, ht, tp);
            
            if(tb > 0) {
                printf("%d bits of each pixel are reserved for the type of point\n", tb);
            }
            
            if(sgn) {
                printf("The 64-bits integers of the pixels are signed.\n");
            }
            
            if((a != 1 || b != 0) && ! isnan(a) && ! isnan(b)) {
                printf("The stored value v of a pixel has to be multiplied by %g and translated by %g.\n",
                       a, b);
            }
            
            break;
            
        case 4: // levm
            pp = mfile_geti(h);
            per = mfile_geti(h);
            a2p = mfile_geti(h);
            prec = mfile_geti(h);
            
            ulong sta = mfile_getl(h);
            ulong ena = mfile_getl(h);
            
            mpfr_init2(r, prec);
            
            mread_mpfr(r, h);
            guard = mfile_getd(h);
            
            segs = ena - sta;
            
            printf("\nThe file contains a misCurve of type (%d, %d) with %ld segments, from unit root %ld to root %ld of order 2^%d.\n",
                   pp, per, segs, sta, ena, a2p);
            
                mpfr_printf("The points are stored with %d bits of precision, the radius is %Rg, the univalent guard is %g.\n",
                            prec, r, guard);
            
            mpfr_clear(r);
            
            break;
            
        case 5: // mpv
            prec = mfile_geti(h);
            ulong count = mfile_getl(h);
            
            printf("\nThe file contains a vector of length %ld of numbers, stored with precision %d.\n",
                   count, prec);
            
            break;
            
        default: // unknown
            printf("The file is damaged or of an unknown format.\n");
            
            return false;
    }
    
    // check if there is an MD5 and print it
    bool hasMD5 = h->len >= h->pos + 16;
    if(hasMD5) {
        h->pos = h->len - 16;
        
        printf("\nMD5 checksum:");
        for (int i = 0; i < 16; i++) {
            int b = mfile_getb(h);
            printf(" %02X", b);
        }
        printf("\n\n");
    } else {
        printf("\nThe file does not have a MD5 checksum.\n\n");
    }
    
    bool md5 = hasMD5 && checkMD5;
    if(fl > fs || md5) { // check for mini-files
        check_content(fn, fl > fs, md5);
    }
    
    return true;
}

// MARK: the help system and the main function

static const char* before = "This task explains the header of a binary file.\n\n";
static const char* after = "\nThe supported types are nset, tmap, levc, levm, bm64 and mpv (with mini-files).\n\n";

static const char *parameters[] = {
    "inputFileName",
    "check MD5"
};

static const char *types[] = {
    "required",
    "optional"
};

static const char *defaults[] = {
    "",
    "0"
};

static const char *descriptions[] = {
    "the name of the input file",
    "1 to load the content of the file and check the MD5 sum"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 2;
static const int columnWidths[] = {18, 18, 18};

/// Prints instructions for usage and some details about the command line arguments.
static void help(void) {
    printf("%s", before);
    char format[50];
    snprintf(format, 45, "    %%-%ds %%-%ds %%-%ds %%s\n", columnWidths[0], columnWidths[1], columnWidths[2]);
    
    printf(format, headers[0], headers[1], headers[2], headers[3]);
    printf("\n");
    for(int i = 0; i < paramCount; i++) {
        printf(format, parameters[i], types[i], defaults[i], descriptions[i]);
    }
    
    printf("%s", after);
}

int header_main(int argc, const char *argv[]) {
    int check = 0;
    
    if(argc < 1) {
        help();
        
        return 1;
    }
    
    if(argc > 1 && sscanf(argv[1], "%d", &check) < 1) {
        help();
        
        return 1;
    }
    
    return ! header((char *) argv[0], check != 0);
}
