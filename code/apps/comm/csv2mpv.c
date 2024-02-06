//
//  csv2mpv.c
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2022.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2022.
//
//  Copyright 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the GNU Public Licence v3.0
//

#include <sys/timeb.h>
#include <string.h>
#include "stopWatch.h"

#include "ntypes.h"
#include "io.h"

#include "csv2mpv.h"

// MARK: import functions

bool csv_read_complex(mpfr_t x, mpfr_t y, FILE *f, int maxDigits) {
    int dig = maxDigits + 10;
    char sx[dig + 1];
    char sy[dig + 1];
    char format[50];
    
    sprintf(format, "%%%ds%%%ds", dig, dig);
    
    if(x == NULL || y == NULL || f == NULL || fscanf(f, format, sx, sy) < 2) {
        return false;
    }
    
    long lx = strlen(sx);
    if(lx > 1 && sx[lx - 1] == ',') {
        sx[lx - 1] = 0;
    }
    
    bool ok = 0 == mpfr_set_str(x, sx, 10, MPFR_RNDN);
    ok = ok && 0 == mpfr_set_str(y, sy, 10, MPFR_RNDN);
    
    return ok;
}

bool csv_read_real(mpfr_t x, FILE *f, int maxDigits) {
    int dig = maxDigits + 10;
    char sx[dig + 1];
    char format[50];
    
    sprintf(format, "%%%ds", dig);
    
    if(x == NULL || f == NULL || fscanf(f, format, sx) < 1) {
        return false;
    }
    
    return 0 == mpfr_set_str(x, sx, 10, MPFR_RNDN);
}

static bool csv2mpv(char * src, char *dst, bool complex, int dig, int pack, bool append) {
    int prec = ceill(dig * LOG10 / LOG2);
    long size = 1L << pack;
    
    char *type = complex ? "complex" : "real";
    int cr = complex ? 2 : 1;
    
    // welcome message
    if(pack == 0) {
        printf("Importing %s to %s, as %s numbers with %d digits.\n", src, dst, type, dig);
    } else {
        printf("Importing %s to %s, as %s numbers with %d digits, in mini-files with %ld points.\n",
               src, dst, type, dig, size);
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
    
    FILE *f = fopen(src, "r");
    if(f == NULL) {
        printf("Could not read %s.\n", src);
        fflush(stdout);
        
        return false;
    }
    
    FILE *g = fopen(dst, append ? "a" : "w");
    if(f == NULL) {
        printf("Could not open %s for writing.\n", dst);
        fflush(stdout);
        fclose(f);
        
        return false;
    }
    
    mpfr_t x, y;
    mpfr_init2(x, prec);
    mpfr_init2(y, prec);
    
    long read = 0, lread = 0;
    long written = 0;
    long vlen = 0;
    long cp = (1L << (pack == 0 ? 20 : pack)) * cr;
    mpv v = mpv_new(prec, cp);
    bool ok = v != NULL;
    if(! ok) {
        printf("Not enough memory for %ld points !\n", cp);
        fflush(stdout);
    }
    
    while(ok && ! feof(f)) {
        ok = complex ? csv_read_complex(x, y, f, 10000) : csv_read_real(x, f, 10000);
        read += ok ? 1 : 0;
        
        if(ok) {
            mpv_set(v, vlen++, x);
            if(complex) {
                mpv_set(v, vlen++, y);
            }
            
            if(vlen == v->count) {
                lapse(&ts, time);
                printf("Read %ld points in %s.\n", read - lread, time);
                fflush(stdout);
                lread = read;
                
                if(pack == 0) {
                    ok = mpv_resize(v, v->count << 1);
                    if(! ok) {
                        printf("Not enough memory, please retry using mini-files (pack > 0).\n");
                    }
                } else {
                    ok = mpv_write_to(v, g, -1);
                    if(ok) {
                        lapse(&ts, time);
                        written += v->count / cr;
                        printf("%ld points written in %s (total %ld).\n\n", v->count / cr, time, written);
                        vlen = 0;
                    } else {
                        printf("Could not write to %s.\n", dst);
                    }
                }
            }
            fflush(stdout);
        }
    }
    
    if(v != NULL && vlen > 0) {
        v->count = vlen;
        
        ok = mpv_write_to(v, g, -1);
        if(ok) {
            lapse(&ts, time);
            printf("%ld points written in %s.\n", v->count / cr, time);
            written += v->count / cr;
            vlen = 0;
        } else {
            printf("Could not write to %s.\n", dst);
        }
    }
    fflush(stdout);
    
    mpfr_clear(x);
    mpfr_clear(y);
    
    mpv_free(v);
    
    fclose(f);
    fclose(g);
    
    lapse(&ats, time);
    printf("%ld points read from %s, %ld points written to %s in %s.\n", read, src, written, dst, time);
    fflush(stdout);
    
    return true;
}

// MARK: the help system and the main function

static const char* before = "This task converts .csv files to .mpv files.\n\n";
static const char* after = "\nFor quick access, there will be several mini-files in the binary .mpv file.\nEach mini-file has it own MD5 check sum.\n\n";

static const char *parameters[] = {
    "inputFileName",
    "outputFileName",
    "complex",
    "digits",
    "pack",
    "append"
};

static const char *types[] = {
    "required",
    "required",
    "required",
    "required",
    "required",
    "optional"
};

static const char *defaults[] = {
    "",
    "",
    "",
    "",
    "",
    "0"
};

static const char *descriptions[] = {
    "the name of the input file",
    "the name of the output file",
    "0 for a vector of real numbers, 1 for a vector of complex numbers",
    "total digits",
    "mini-file of size 2^pack; at least 12 or 0 to create only one mini-file",
    "1 to append points to the CSV file, 0 to erase existing points"
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

int csv2mpv_main(int argc, const char * argv[]) {
    int dig = 26, append = 1, complex = 0, pack = 0;
    
    bool ok = argc >= 5 && sscanf(argv[2], "%d", &complex) == 1;
    ok = ok && sscanf(argv[3], "%d", &dig) == 1 &&
                             (dig >= 15 || dig == 0) && dig <= 10000;
    ok = ok && sscanf(argv[4], "%d", &pack) == 1 &&
        (pack >= 12 || pack == 0);
    ok = ok && (argc < 6 || sscanf(argv[5], "%d", &append) == 1);
    
    if(! ok) {
        help();
        
        return 1;
    }
    
    return ! csv2mpv((char *) argv[0], (char *) argv[1], complex != 0, dig, pack, append != 0);
}
