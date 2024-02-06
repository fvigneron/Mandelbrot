//
//  render.c
//  Mandel
//
//  Created by MIHALACHE Nicolae on 10/19/22.
//  Copyright © 2022 MIHALACHE Nicolae. All rights reserved.
//

#include <sys/timeb.h>
#include "stopWatch.h"
#include "io.h"
#include "mthread.h"

#include "render.h"
#include "iterates.h"

static bool bitmaps(char *fn, int threads, int start, int end) {
    FILE *f = fopen(fn, "r");
    if(f == NULL) {
        printf("Could not read %s !\n", fn);
        
        return false;
    }
    
    if(start == 1 && end == INT_MAX) {
        printf("Rendering all bitmaps described in %s\n", fn);
    } else {
        printf("Rendering bitmaps %d to %d described in %s\n", start, end, fn);
    }
    printf(threads == 1 ? "Using %d thread.\n" : "Using %d threads.\n", threads);
    io_print_start();
    
    bool allOk = true;
    int lc = 1, imgc = 0;
    strings line = io_read_csv_line(f);
    while(line != NULL && line->count < 10) {
        line = io_read_csv_line(f);
    }
    while(line != NULL && lc < start) {
        line = io_read_csv_line(f);
        while(line != NULL && line->count < 10) {
            line = io_read_csv_line(f);
        }
        lc ++;
    }
    
    while(line != NULL && lc <= end) {
        iterBmap it;
        it.usesDelta = false;
        
        // CSV structure: type, subtype, precision, max_iter, log2_zoom, h-res, v-res,
        // x_cent, y_cent, fileName, Re(c) (or n), Im(c) (or k)
        int args[] = {10, 12, 11, 12};
        bool ok = sscanf(line->fields[0], "%d", &it.type) == 1 && it.type >= ITER_TYPE_M
                    && it.type <= ITERT_COUNT && line->count == args[it.type - 1];
        ok = ok && sscanf(line->fields[1], "%d", &it.subType) == 1 && it.subType >= 1;
        ok = ok && sscanf(line->fields[2], "%d", &it.prec) == 1 && it.prec >= 64;
        ok = ok && sscanf(line->fields[3], "%lu", &it.maxIter) == 1 && it.maxIter >= 64;
        ok = ok && sscanf(line->fields[4], "%d", &it.r.tpow) == 1
                    && it.r.tpow >= -10 && it.r.tpow <= 1000000;
        ok = ok && sscanf(line->fields[5], "%d", &it.r.w) == 1 && it.r.w >= 64;
        ok = ok && sscanf(line->fields[6], "%d", &it.r.h) == 1 && it.r.h >= 64;
        
        fp80 c80;
        mpc c;
        if(ok && it.prec <= 64) { // fp80
            ok = sscanf(line->fields[7], "%Lf", &c80->x) == 1;
            ok = ok && sscanf(line->fields[8], "%Lf", &c80->y) == 1;
            ok = ok && drect_center80((drect) &it.r, c80);
        } else if(ok) {
            mpc_init(c, it.prec);
            mpc_init((mpc_ptr) it.delta, it.prec);
            it.usesDelta = true;
            
            ok = mpfr_set_str(c->x, line->fields[7], 10, MPFR_RNDN) == 0;
            ok = ok && mpfr_set_str(c->y, line->fields[8], 10, MPFR_RNDN) == 0;
            ok = ok && drect_center((drect) &it.r, (mpc_ptr) it.delta, c);
        }
        
        switch(it.type) {
            case ITER_TYPE_J:
                if(it.prec == 64) { // fp80
                    ok = ok && sscanf(line->fields[10], "%Lf", &it.c80->x) == 1;
                    ok = ok && sscanf(line->fields[11], "%Lf", &it.c80->y) == 1;
                } else {  // MPFR
                    mpc_init((mpc_ptr) it.c, it.prec);
                    ok = ok && mpfr_set_str((mpfr_ptr) it.c->x, line->fields[7], 10, MPFR_RNDN) == 0;
                    ok = ok && mpfr_set_str((mpfr_ptr) it.c->y, line->fields[8], 10, MPFR_RNDN) == 0;
                }
                break;
                
            case ITER_TYPE_NM:
                ok = ok && sscanf(line->fields[11], "%d", &it.prePer) == 1;
            case ITER_TYPE_NH:
                ok = ok && sscanf(line->fields[10], "%d", &it.per) == 1;
                break;
        }
        
        if(! ok) {
            printf("Could not interpret line %d !\n", lc);
            allOk = false;
            
            strings_free(line);
            if(it.prec > 64 && it.type == 2) {
                mpc_clear((mpc_ptr) it.c);
            }
            
            continue;
        }
        struct timeb te;
        char time[80];
        
        char outFile[1000];
        printf("Rendering bitmap(s) specified at line %d ... ", lc);
        fflush(stdout);
        
        ftime(&te);
        bmap *bmps = iter_bitmaps(&it, threads);
        lapse(&te, time);
        printf("done in %s\n", time);
        fflush(stdout);
        
        if(bmps != NULL) {
            int st = it.subType, t = 1, pos = 0;
            while(t <= st) {
                if(st & t) {
                    if(st == t) {
                        snprintf(outFile, 1000, "%s.bmap", line->fields[9]);
                    } else {
                        snprintf(outFile, 1000, "%s_%d.bmap", line->fields[9], t);
                    }
                    
                    if(bmap_write(bmps[pos], outFile, false)) {
                        printf("    Written bitmap of type %d, subType %d, to %s\n", it.type, t, outFile);
                        imgc ++;
                    } else {
                        printf("    Could not write a bitmap to %s !!!\n", outFile);
                        allOk = false;
                    }
                }
                
                t <<= 1;
            }
        } else {
            printf("Could not render the bitmaps !\n");
            allOk = false;
        }
        
        strings_free(line);
        if(it.prec > 64 && it.type == 2) {
            mpc_clear((mpc_ptr) it.c);
        }
        
        line = io_read_csv_line(f);
        while(line != NULL && line->count < 10) {
            line = io_read_csv_line(f);
        }
        lc ++;
    }
    
    fclose(f);
    
    if(allOk) {
        printf("\nComputed and written %d bitmap(s) with no errors.\n\n", imgc);
    } else {
        if(imgc == 0) {
            printf("\nCould not compure or write any bitmap !!!\n\n");
        } else {
            printf("\nComputed and written %d bitmap(s) with SOME errors !\n\n", imgc);
        }
    }
    
    return allOk;
}

// MARK: the help system and the main function

static const char* before = "This task renders the bitmaps described in the input file.\n\n";
static const char* after = "\nWrites the result in the binary .bmap output file.\n\n";

static const char *parameters[] = {
    "descFile",
    "threads",
    "start",
    "end"
};

static const char *types[] = {
    "required",
    "optional",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "1",
    "1",
    "max"
};

static const char *descriptions[] = {
    "the name of the bitmap descriptions input file",
    "the number of threads to use",
    "the first (valid) line in the CSV file to render",
    "the last (valid) line in the CSV file to render"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 4;
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
}

int render_main(int argc, const char * argv[]) {
    if(argc < 1) {
        help();
        
        return 1;
    }
    
    char *fn = (char *) argv[0];
    int threads = 1, start = 1, end = INT_MAX;
    
    if(argc > 1) {
        if(sscanf(argv[1], "%d", &threads) < 1 || threads < 1 || threads > MTH_MAX_THREADS) {
            help();
            
            return 1;
        }
    }
    
    if(argc > 2) {
        if(sscanf(argv[2], "%d", &start) < 1 || start < 1) {
            help();
            
            return 1;
        }
    }
    
    if(argc > 3) {
        if(sscanf(argv[3], "%d", &end) < 1 || end < start) {
            help();
            
            return 1;
        }
    }
    
    return ! bitmaps(fn, threads, start, end);
}
