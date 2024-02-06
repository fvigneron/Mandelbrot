//
//  mpv2csv.c
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
#include "mpv2csv.h"

// MARK: export functions

static bool mpv2csv(char * src, char *dst, bool complex, int dig, long start, long count, bool append) {
    char *type = complex ? "complex" : "real";
    int cr = complex ? 2 : 1;
    
    // welcome message
    if(start == 0 && count == 0) {
        printf("Exporting %s to %s,\nas %s numbers.\n", src, dst, type);
    } else if(start == 0) {
        printf("Exporting the first %ld points of %s to %s,\nas %s numbers.\n", count, src, dst, type);
    } else if(count == 0) {
        printf("Exporting the %s to %s,\nas %s numbers, starting at position %ld.\n", src, dst, type, start);
    } else { // here both start, count > 0
        printf("Exporting %ld points from %s to %s,\nas %s numbers, starting at position %ld.\n", count, src, dst, type, start);
    }
    if(append && io_size_of(dst) > 0) {
        printf("The points will be appended to %s.\n", dst);
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
    ulong pcount = mfile_getl(header) / cr;
    
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
            pcount = mfile_getl(header) / cr;
            
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
        
        pcount = s->count / cr;
        lapse(&ts, time);
        printf("Read %ld points in %s.\n", pcount, time);
        fflush(stdout);
        
        long st = nfirst ? 0 : start - pdone;
        long len = pdone + pcount < start + count || count == 0 ?
                    pcount - st :
                    start + count - pdone - st;
        int digits = dig == 0 ? ceill(s->prec * LOG2 / LOG10) : dig;
        ok = mpv_write_csv(s, dst, complex, digits, st, len, nfirst ? true : append);
        
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
    "complex",
    "digits",
    "start",
    "count",
    "append"
};

static const char *types[] = {
    "required",
    "required",
    "required",
    "optional",
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
    "0",
    "1"
};

static const char *descriptions[] = {
    "the name of the input file",
    "the name of the output file",
    "0 for a vector of real numbers, 1 for a vector of complex numbers",
    "total digits; 0 to use to the precision in the .mpv file",
    "first index in the input file, inclusive",
    "the number of points to export, 0 to export the entire file",
    "1 to append points to the CSV file, 0 to erase existing points"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 7;
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

int mpv2csv_main(int argc, const char * argv[]) {
    int dig = 0, append = 1, complex = 0;
    long start = 0, count = 0;
    
    bool ok = argc >= 3 && sscanf(argv[2], "%d", &complex) == 1;
    ok = ok && (argc < 4 || (sscanf(argv[3], "%d", &dig) == 1 &&
                             (dig >= 15 || dig == 0) && dig <= 10000));
    ok = ok && (argc < 6 || (sscanf(argv[4], "%ld", &start) == 1 &&
                sscanf(argv[5], "%ld", &count) == 1 && start >= 0 && count >= 0));
    ok = ok && (argc < 7 || sscanf(argv[6], "%d", &append) == 1);
    
    if(! ok) {
        help();
        
        return 1;
    }
    
    return ! mpv2csv((char *) argv[0], (char *) argv[1], complex, dig, start, count, append);
}
