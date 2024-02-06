//
//  mpv2mpv.c
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2022.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2021.
//
//  Copyright 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the GNU Public Licence v3.0
//

#include <sys/timeb.h>
#include <string.h>
#include "stopWatch.h"
#include "io.h"

#include "ntypes.h"
#include "mpv2mpv.h"

// MARK: export functions

static FILE *fo = NULL;
static mpv d = NULL;
static ulong dpos = 0;

static bool m2m_transfer(mpv src, long st, long len) {
    long left = len;
    long pos = st;
    
    bool ok = true;
    while(ok && dpos + left >= d->count) {
        long batch = d->count - dpos;
        ok = mpv_copy_partial(d, dpos, src, pos, batch);
        
        ok = ok && mpv_write_to(d, fo, -1);
        
        pos += batch;
        left -= batch;
        dpos = 0;
    }
    
    if(ok && left > 0) {
        ok = mpv_copy_partial(d, dpos, src, pos, left);
        dpos += left;
    }
    
    return ok;
}

static bool m2m_flush(void) {
    bool ok = true;
    
    if(dpos > 0) {
        d->count = dpos;
        
        ok = mpv_write_to(d, fo, -1);
    }
    
    fclose(fo);
    mpv_free(d);
    
    return ok;
}

static bool mpv2mpv(char * src, char *dst, int dig, long start, long count, int tpm) {
    // welcome message
    int prec = 1 + (int) ceill((dig + 1) * LOG10 / LOG2);
    if(start == 0 && count == 0) {
        printf("Exporting %s to %s, with %d digits, or %d bits.\n", src, dst, dig, prec);
    } else if(start == 0) {
        printf("Exporting the first %ld points of %s to %s,\nwith %d digits, or %d bits.\n",
               count, src, dst, dig, prec);
    } else if(count == 0) {
        printf("Exporting the %s to %s,\nwith %d digits, or %d bits, starting at position %ld.\n",
               src, dst, dig, prec, start);
    } else { // here both start, count > 0
        printf("Exporting %ld points from %s to %s,\nwith %d digits, or %d bits, starting at position %ld.\n",
               count, src, dst, dig, prec, start);
    }
    if(tpm > 0) {
        printf("The destination file will contain mini-files with %ld numbers.\n", 1L << tpm);
    }
    
    io_print_start();
    fflush(stdout);
    
    char time[100];
    struct timeb ats, ts;
    ftime(&ats);
    
    ftime(&ts);
    
    long fpos = 0;
    long pdone = 0;
    FILE *f = fopen(src, "r");
    mfile header = f == NULL ? NULL : mfile_read_header(f, fpos, MPV_FILE_ID, MPV_HEADER_LEN);
    if(header == NULL) {
        printf("Could not read points from %s.\n", src);
        if(f != NULL) {
            fclose(f);
        }
        
        return false;
    }
    
    header->pos = 8; // position of mini-file size
    long mfsize = mfile_getl(header);
    
    header->pos = 24; // position of points count
    ulong pcount = mfile_getl(header);
    ulong dcount = tpm == 0 ? pcount : 1L << tpm;
    
    int skipped = 0;
    bool ok = true;
    while(ok && pdone + pcount <= start) { // skip the mini-files that are not needed
        fpos += mfsize;
        pdone += pcount;
        
        mfile_free(header);
        header = mfile_read_header(f, fpos, MPV_FILE_ID, MPV_HEADER_LEN);
        ok = header != NULL;
        
        if(ok) {
            header->pos = 8; // position of mini-file size
            mfsize = mfile_getl(header);
            
            header->pos = 24; // position of points count
            pcount = mfile_getl(header);
            
            skipped ++;
        }
    }
    
    if(! ok) {
        fclose(f);
        printf("Could only read %ld points from %s, not enough to export starting at position %ld.\n",
               pdone, src, start);
        
        return false;
    }
    
    if(skipped > 0) {
        lapse(&ts, time);
        printf("Skipped %d mini-files containing %ld points in %s.\n\n", skipped, pdone, time);
        fflush(stdout);
    }
    
    fo = fopen(dst, "w");
    if(fo == NULL) {
        printf("Could not open the file %s for writing.\n", dst);
        
        fclose(f);
        
        return false;
    }
    
    d = mpv_new(prec, dcount);
    if(d == NULL) {
        printf("Not enough available memory for two mini-files.\n");
        
        fclose(f);
        fclose(fo);
        
        return false;
    }
    
    long exp = 0;
    bool nfirst = false;
    while(ok && (pdone < start + count || count == 0)) { // read interesting points and export
        mpv s = mpv_read_from(f, nfirst ? -1 : fpos);
        ok = s != NULL;
        
        if(! ok) {
            if(count > 0) {
                printf("Could only read %ld points from %s, not enough to export %ld points starting at position %ld.\n",
                       pdone, src, count, start);
                fflush(stdout);
            }
            
            continue;
        }
        
        pcount = s->count;
        lapse(&ts, time);
        printf("Read %ld points in %s.\n", pcount, time);
        fflush(stdout);
        
        long st = nfirst ? 0 : start - pdone;
        long len = pdone + pcount < start + count || count == 0 ?
                    pcount - st :
                    start + count - pdone - st;
        ok = m2m_transfer(s, st, len);
        
        if(! ok) {
            printf("Could not write %ld points to %s.\n", len, dst);
            fflush(stdout);
            mpv_free(s);
            
            continue;
        }
        
        exp += len;
        
        lapse(&ts, time);
        printf("Written %ld points in %s.\n\n", len, time);
        fflush(stdout);
        
        nfirst = true;
        pdone += pcount;
        
        mpv_free(s);
    }
    
    m2m_flush();
    
    fclose(f);
    
    lapse(&ats, time);
    printf("Exported a total of %ld points to %s in %s.\n\n", exp, dst, time);
    
    return true;
}

// MARK: the help system and the main function

static const char* before = "This task [partially] converts .mpv files to .csv files.\n\n";
static const char* after = "\nThe \"end\" position (optional parameter) is not included in the list.\n\n";

static const char *parameters[] = {
    "inputFileName",
    "outputFileName",
    "digits",
    "start",
    "count",
    "mini-size"
};

static const char *types[] = {
    "required",
    "required",
    "required",
    "optional",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "",
    "",
    "0",
    "0",
    "0"
};

static const char *descriptions[] = {
    "the name of the input file",
    "the name of the output file",
    "total digits, at least 10",
    "first index in the input file, inclusive",
    "the number of points to export, 0 to export the entire file",
    "the power of two > 9 that gives the size of the mini-file, 0 to keep as the input"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 6;
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

int mpv2mpv_main(int argc, const char * argv[]) {
    int dig = 0, tpm = 0;
    long start = 0, count = 0;
    
    bool ok = argc >= 3;
    ok = ok && (sscanf(argv[2], "%d", &dig) == 1 &&
                             (dig >= 15 || dig == 0) && dig <= 10000);
    ok = ok && (argc < 5 || (sscanf(argv[3], "%ld", &start) == 1 &&
                sscanf(argv[4], "%ld", &count) == 1 && start >= 0 && count >= 0));
    ok = ok && (argc < 6 || (sscanf(argv[5], "%d", &tpm) == 1 && (tpm > 9 || tpm == 0)));
    
    if(! ok) {
        help();
        
        return 1;
    }
    
    return ! mpv2mpv((char *) argv[0], (char *) argv[1], dig, start, count, tpm);
}
