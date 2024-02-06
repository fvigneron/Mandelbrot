//
//  hypCount.c
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
#include <sys/timeb.h>

#include "stopWatch.h"
#include "hypRaw.h"
#include "hypRawCount.h"
#include "hypCount.h"
#include "mandel.h"
#include "io.h"
#include "memFile.h"


u128_ptr count_load_and_check_nset(char *fileName, ulong *count, ulong *realCount, ulong *eps, bool verbose) {
    if(fileName == NULL || count == NULL || realCount == NULL || eps == NULL) {
        return NULL;
    }
    
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        return NULL;
    }
    
    char time[80];
    struct timeb lts, md5ts;
    ftime(&lts);
    
    // read an analyse the header of the binary file
    mfile h = mfile_read_header(f, 0, NSET_FILE_ID, NSET_HEADER_LEN);
    if(h == NULL) {
        return NULL;
    }
    
    h->pos = 8; // skip the file ID
    ulong flen = mfile_getl(h);    // the file size
    uint hlen = mfile_geti(h);     // the header size
    
    *count = mfile_getl(h);        // the count of points in the set
    *realCount = mfile_getl(h);    // the count of real points in the set
    *eps = mfile_getl(h);          // the max error on each coordinate
    
    ulong pslen = hlen + *count * sizeof(u128c_struct);
    u128_ptr list = malloc(*count * sizeof(u128c_struct));
    
    if(flen != pslen || fread_block(list, sizeof(u128c_struct), *count, FILE_BLOCK, f) != *count) {
        free(list);
        mfile_free(h);
        fclose(f);
        
        return NULL;
    }
    
    fclose(f);
    
    lapse(&lts, time);
    printf("File %s loaded in %s, ", fileName, time);
    fflush(stdout);
    
    ftime(&md5ts);
    
    // check md5
    if(hlen < NSET_HEADER_LEN) {
        mfile_free(h);
        free(list);
        
        return false;
    }
    
    byte fmd5sum[16];
    mfile_getbs(h, fmd5sum, 16);
    
    // compute MD5 checksum
    MD5_CTX md5;
    byte md5sum[16];
    
    MD5_Init(&md5);
    MD5_Update(&md5, h->data, (uint) h->len - 16); // except the MD5
    
    if(*count > 0) {
        long i = 0;
        for (; i <= ((long) *count) - 16; i += 16) {
            MD5_Update(&md5, list + i, 16 * sizeof(u128c_struct));
        }
        
        if(i < *count) {
            MD5_Update(&md5, list + i, (uint) (*count - i) * sizeof(u128c_struct));
        }
    }
    
    MD5_Final(md5sum, &md5);
    
    bool ok = true;
    for (int i = 0; i < 16 && ok; i++) {
        ok = ok && fmd5sum[i] == md5sum[i];
    }
    
    mfile_free(h);
    lapse(&lts, time);
    
    if(! ok) {
        printf("MD5 checksum is NOT correct (computed in %s).\n", time);
        fflush(stdout);
        
        free(list);
        
        return NULL;
    }

    printf("MD5 checksum is correct (computed in %s).\n", time);
    fflush(stdout);
    
    return list;
}

bool count_check_increasing(u128_ptr list, ulong count) {
    for (ulong i = 0; i < count - 1; i++) {
        if(! u128_sless(list + i, list + (i + 1))) {
            return false;
        }
    }
    
    return true;
}

bool count_check_min_dist(u128_ptr list, ulong count, ulong eps) {
    for (ulong i = 0; i < count - 1; i++) {
        bool done = false;
        uint128 dx = 0;
        
        for (ulong j = i + 1; j < count && ! done; j++) {
            dx = list[j].x - list[i].x;
            if(dx > eps) {
                done = true;
                
                continue;
            }

            if(u128_eq(list + i, list + j, eps)) {
                return false;
            }
        }        
    }
    
    return true;
}

bool count_check_boundary(u128_ptr left, ulong lc, u128_ptr right, ulong rc, ulong eps) {
    u128_ptr l = left + (lc - 1);
    bool done = false;
    uint128 dx = 0, lx = l->x;
    
    for (ulong i = 0; i < rc && ! done; i++) {
        dx = right[i].x - lx;
        if(dx > eps) {
            done = true;
            
            continue;
        }

        if(u128_eq(l, right + i, eps)) {
            return false;
        }
    }
    
    done = false;
    dx = 0;
    uint128 rx = right->x;
    
    for (long i = lc - 1; i >= 0 && ! done; i--) {
        dx = rx - left[i].x;
        if(dx > eps) {
            done = true;
            
            continue;
        }

        if(u128_eq(left + i, right, eps)) {
            return false;
        }
    }    
    
    return true;
}

ulong count_real_count(u128_ptr list, ulong count, ulong eps) {
    ulong r = 0;
    ulong eps2 = eps >> 1;
    
    for (ulong i = 0; i < count; i++) {
        if(list[i].y <= eps2) {
            r ++;
        }
    }
    
    return r;
}

// MARK: static functions

static bool count(int per) {
    printf("Counting hyperbolic centers of period %d.", per);
    
    char time[80];
    struct timeb ats, ts;
    ftime(&ats);
    
    time_stamp(time, 80, true, true);
    printf(" Started at %s\n\n", time);
    fflush(stdout);
    
    int res = hyp_resultsCount(per);
    u128_ptr prev = NULL, list = NULL;
    ulong prevc = 0, count = 0, real = 0, eps = 0, totc = 0, totr = 0, prevEps = 0;
    char fn[100];
    
    
    int jobs = 1 << (per - 27);
    int digits = 0;
    while(jobs > 0) {
        jobs /= 10;
        digits ++;
    }
    
    char pf[20];
    snprintf(pf, 20, " %%0%dd", digits);
    
    for (int i = 0; i < res; i++) {
        list = NULL;
        
        if(hyp_fileName(fn, 99, per, i) > 10) {
            list = count_load_and_check_nset(fn, &count, &real, &eps, true);
        }
        
        ftime(&ts);
        
        if(list == NULL) {
            printf("\nCould not load results file with index %d, will stop here.\n", i);
            fflush(stdout);
            
            if(prev != NULL) {
                free(prev);
            }
            
            return false;
        }
        
        // announce next files for external buffering
        printf("        prepare:");
        for (int j = i + 1; j < i + 21 && j < res; j++) {
            printf(pf, j);
        }
        printf("\n");
        
        if(! count_check_increasing(list, count)) {
            printf("The points in %s are not in strictly increasing order, count failed.\n", fn);
            fflush(stdout);
            
            free(list);
            if(prev != NULL) {
                free(prev);
            }
            
            return false;
        }
        
        if(count > 1 && list[count - 1].x - list[0].x <= eps) { // no overflow, as list has been checked to be increasing
            // needs at least two points to be a reason for failure
            printf("The points in %s are contained in very thin vertical strip, count failed.\n", fn);
            fflush(stdout);
            
            free(list);
            if(prev != NULL) {
                free(prev);
            }
            
            return false;
        }
        
        if(i > 0 && ! count_check_boundary(prev, prevc, list, count, eps)) {
            printf("Boundary intersection of files of %d and %d, count failed.\n",
                   i -1, i);
            fflush(stdout);
            
            free(list);
            if(prev != NULL) {
                free(prev);
            }
            
            return false;
        }
        
        if(! count_check_min_dist(list, count, eps)) {
            printf("The points in %s are not separated by %lg on at least one coordinate, count failed.\n",
                   fn, ldexp(eps, -126));
            fflush(stdout);
            
            free(list);
            if(prev != NULL) {
                free(prev);
            }
            
            return false;
        }
        
        if(real != count_real_count(list, count, eps)) {
            printf("The number of real parameters stored in %s is incorrect, count failed.\n", fn);
            fflush(stdout);
            
            free(list);
            if(prev != NULL) {
                free(prev);
            }
            
            return false;
        }
        
        if(i > 0) {
            free(prev);
        }
        
        totc += count;
        totr += real;
        
        prev = list;
        prevc = count;
        prevEps = eps;
        
        lapse(&ts, time);
        printf("  - %ld points of file %s checked in %s\n", count, fn, time);
        fflush(stdout);
    }
    
    free(list);
    ulong tr = 2 * totc - totr;
    ulong hc = mandel_hyp_count(per);
    
    if(hc != tr) {
        printf("\nThe files contain a list of unique points, but the total count is %ld instead of %ld\n",
               tr, hc);
        fflush(stdout);
        
        return false;
    }
    
    printf("\nThe files contain all %ld hyperbolic centers (%ld stored) of period %d,\nin strictly increasing order. Any two points are at distance at least %lg from each other.\n\n",
          hc, totc, per, ldexp(eps, -126));
    
    lapse(&ats, time);
    printf("Total loading an computing time %s\n", time);
    fflush(stdout);
    
    return true;
}

// MARK: the help system and the main function

static const char* before = "This task counts the hyperbolic centers of given period, checking for unicity.\n\n";
static const char* after = "\nIt also checks that the nset files are correct with an MD5 checksum.\n\n";

static const char *parameters[] = {
    "period"
};

static const char *types[] = {
    "required"
};

static const char *defaults[] = {
    ""
};

static const char *descriptions[] = {
    "the period, integer, at least 1, at most 41"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 1;
static const int columnWidths[] = {18, 18, 18};

/// Prints instructions for usage and some details about the command line arguments.
static void help(void) {
    printf("%s", before);
    char format[50];
    snprintf(format, 45, "    %%-%ds %%-%ds %%-%ds %%s\n",
             columnWidths[0], columnWidths[1], columnWidths[2]);
    
    printf(format, headers[0], headers[1], headers[2], headers[3]);
    printf("\n");
    for(int i = 0; i < paramCount; i++) {
        printf(format, parameters[i], types[i], defaults[i], descriptions[i]);
    }
    
    printf("%s", after);
    fflush(stdout);
}

int hyp_count_main(int argc, const char * argv[]) {
    int per;
    
    if(argc < 1 || sscanf(argv[0], "%d", &per) < 1 || per < 1 || per > HYP_RAW_MAX_PER) {
        help();
        
        return 1;
    }
    
    return ! count(per);
}
