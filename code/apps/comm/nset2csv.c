//
//  nset2csv.c
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

#include <sys/timeb.h>
#include <string.h>

#include "stopWatch.h"
#include "ntypes.h"
#include "nset2csv.h"
#include "hypRaw.h"
#include "nSet.h"
#include "mandel.h"

// MARK: static functions 

/// @brief Erases the content and writes the points of the set @c ps to a CSV file with the given precision in @c digits.
///
/// The points are assumed to be hyperbolic centers of given @c period. If the number of digits is @c <=30, the
/// function fails. Otherwise, it increases the precision of the results as necessary by the Newton method (at most
/// 1500 digits).
///
/// @param ps the set of hyperbolic centers of period @c period
/// @param fileName the name of the CSV file to write to
/// @param digits the number of digits @c >30 and @c <=1500
/// @param period the period of the hyperbolic centers
/// @param append @ref true to append points to the CSV file, @ref false to erase existing points
///
/// @return @ref true if successfull, @ref false otherwise
static bool hyp_refine_and_write_csv(nset_t ps, char *fileName, int digits, int period, bool append) {
    if(fileName == NULL || ! nset_valid(ps) || digits <= 30 || digits > 1500 ||
       period < 3 || period > 41) {
        return false;
    }
    
    FILE *f = fopen(fileName, append ? "a" : "w");
    if(f == NULL) {
        return false;
    }
    
    int maxPrec = 64 * (1 + ((digits - 1) / 15)) - 8;
    mpc c, sp;
    mpc_init(c, maxPrec);
    mpc_init(sp, 128);
    
    u128 p;
    char rf[100], cf[100];
    sprintf(rf, "%%.%dRf, 0\n", digits);
    sprintf(cf, "%%.%dRf, %%.%dRf\n", digits, digits);

    bool err = false;
    for (long i = 0; i < ps->count && ! err; i ++) {
        nset_point(p, ps, i);
        u128_get(sp, p);
        
        if(hyp_refine(c, sp, period, digits)) {
            if(p->y <= ps->eps) {
                err = mpfr_fprintf(f, rf, c->x) < digits + 5;
            } else {
                err = mpfr_fprintf(f, cf, c->x, c->y) < 2 * digits + 6;
            }
        }
    }

    fclose(f);
    
    mpc_clear(c);
    
    return ! err;
}

// MARK: refinement functions

bool hyp_refine(mpc res, mpc sp, int per, int digits) {
    if(res == NULL || sp == NULL || digits <= 30 || digits > 1500 || per < 3 || per > HYP_RAW_MAX_PER) {
        return false;
    }
    
    int maxPrec = 64 * (1 + ((digits - 1) / 15)) - 8, prec = 170;
    mpc c, nt, bc;
    mpc_init(c, prec);
    mpc_init(nt, prec);
    mpc_init(bc, prec);
    
    int scale = -(7 * digits / 2);
    bool err = false;
    
    int ref = 0, texp;
    do { // Newton method with adaptive precision
        ref ++;
        mandel_nt(nt, c, per);
        
        mpc_sub(c, c, nt);
        texp = (int) mpc_2exp(nt);
        
        if(texp >= scale) {
            mpc_set_prec(bc, prec);
            mpc_set(bc, c);
            
            prec = -2 * texp - 70;
            prec = prec > maxPrec ? maxPrec : prec < 240 ? 240 : prec;
            
            mpc_set_prec(c, prec);
            mpc_set_prec(nt, prec);
            
            mpc_set(c, bc);
        }
    } while(texp >= scale && ref < 20);
    
    if(texp >= scale) {
        err = true;
    } else {
        mpc_set(res, c);
    }
    
    mpc_clear(c);
    mpc_clear(nt);
    mpc_clear(bc);
    
    return ! err;
}

bool hyp_refine_nset(nset_t ps, int per, ulong *maxErr) {
    if(! nset_valid(ps) || per < 1 || per > HYP_RAW_MAX_PER) {
        return false;
    }
        
    long me = 0;
    u128 p, np;
    mpc c, nt;
    mpc_init(c, 170);
    mpc_init(nt, 170);
    
    int scale = -128;
    bool err = false;
    for (long i = 0; i < ps->count; i ++) {
        nset_point(p, ps, i);

        u128_get(c, p);

        int ref = 0, conv = 0;
        do {
            ref ++;
            mandel_nt(nt, c, per);
            conv = mpc_2exp(nt) < scale;
            
            mpc_sub(c, c, nt);
        } while(! conv && ref < 3);
        
        if(! conv) {
            err = true;
            
            continue;
        }
        
        u128_set(np, c);
        err = ! nset_replace(ps, np);
        
        if(maxErr != NULL && (np->ul[0] != p->ul[0] || np->ul[3] != p->ul[3])) {
            long ex = np->ul[0] - p->ul[0];
            long ey = np->ul[3] - p->ul[3];
            
            ex = ex < 0 ? -ex : ex;
            ey = ey < 0 ? -ey : ey;
            
            me = ex > me ? ex : me;
            me = ey > me ? ey : me;
        }
    }
    
    if(maxErr != NULL) {
        *maxErr = me;
    }
    
    return ! err;
}

bool mis_refine(mpc res, mpc sp, int pp, int per, int digits) {
    if(res == NULL || sp == NULL || digits <= 30 || digits > 1500 ||
       pp < 2 || per < 1 || pp + per > HYP_RAW_MAX_PER) {
        return false;
    }
    
    int maxPrec = 64 * (1 + ((digits - 1) / 15)) - 8, prec = 180;
    mpc c, nt, bc;
    mpc_init(c, prec);
    mpc_init(nt, prec);
    mpc_init(bc, prec);
    
    mpc_set(c, sp);

    int scale = -(7 * digits / 2);
    bool err = false;
    
    int ref = 0, texp;
    do { // Newton method with adaptive precision
        ref ++;
        mandel_mis_nt(nt, c, pp, per);
        
        texp = (int) mpc_2exp(nt);
        
        if(texp >= scale) {
            mpc_set_prec(bc, prec);
            mpc_sub(bc, c, nt);

            prec = -2 * texp - 70;
            prec = prec > maxPrec ? maxPrec : (prec < 240 ? 240 : prec);

            mpc_set_prec(c, prec);
            mpc_set_prec(nt, prec);

            mpc_set(c, bc);
        } else {
            mpc_sub(c, c, nt);
        }
    } while(texp >= scale && ref < 20);
    
    if(texp >= scale) {
        err = true;
    } else {        
        mpc_set(res, c);
    }
    
    mpc_clear(c);
    mpc_clear(nt);
    mpc_clear(bc);
    
    return ! err;
}

bool mis_refine_nset(nset_t ps, int pp, int per, ulong *maxErr) {
    if(! nset_valid(ps) || pp < 2 || per < 1 || pp + per > HYP_RAW_MAX_PER) {
        return false;
    }
        
    long me = 0;
    u128 p, np;
    mpc c, nt;
    mpc_init(c, 170);
    mpc_init(nt, 170);
    
    int scale = -128;
    bool err = false;
    for (long i = 0; i < ps->count; i ++) {
        nset_point(p, ps, i);

        u128_get(c, p);

        int ref = 0, conv = 0;
        do {
            ref ++;
            mandel_mis_nt(nt, c, pp, per);
            conv = mpc_2exp(nt) < scale;
            
            mpc_sub(c, c, nt);
        } while(! conv && ref < 3);
        
        if(! conv) {
            err = true;
            
            continue;
        }
        
        u128_set(np, c);
        err = ! nset_replace(ps, np);
        
        if(maxErr != NULL && (np->ul[0] != p->ul[0] || np->ul[3] != p->ul[3])) {
            long ex = np->ul[0] - p->ul[0];
            long ey = np->ul[3] - p->ul[3];
            
            ex = ex < 0 ? -ex : ex;
            ey = ey < 0 ? -ey : ey;
            
            me = ex > me ? ex : me;
            me = ey > me ? ey : me;
        }
    }
    
    if(maxErr != NULL) {
        *maxErr = me;
    }
    
    return ! err;
}

// MARK: export and general text file functions

bool hyp_export2csv(char *in, char *out, int digits, long st, long count, int per, bool append) {
    struct timeb ts;
    ftime(&ts);

    nset_t ps;
    nset_init(ps, HYP_RAW_SET_EPS);
    
    if(count == 0) { // read the entire file
        if(! nset_read(ps, in, true)) {
            printf("Could not read %s\n", in);
            
            return 0;
        }
    } else if(! nset_read_partial(ps, in, 0, st, count)) {
        printf("Could not read %s\n", in);
        
        return 0;
    }
    
    char time[100];
    lapse(&ts, time);
    printf("Read %lu points from %s in %s\n", ps->count, in, time);
    
    ftime(&ts);
    
    long en = st + (count <= 0 || count > ps->count ? ps->count : count);
    
    if(digits == 0) {
        if(nset_write(ps, out)) {
            lapse(&ts, time);
            printf("%ld points written to %s int %s\n", ps->count, out, time);
        } else {
            printf("Could not write to %s\n", out);
        }
    } else if(digits > U128_DIGITS) {
        if(hyp_refine_and_write_csv(ps, out, digits, per, append)) {
            lapse(&ts, time);
            printf("%ld refined points written to %s in %s\n", en - st, out, time);
        } else {
            printf("Could not refine roots and write to %s\n", out);
        }
    } else if(nset_write_csv(ps, out, digits, append)) {
        lapse(&ts, time);
        printf("%ld points written to %s in %s\n", en - st, out, time);
    } else {
        printf("Could not write to %s\n", out);
    }
    
    nset_clear(ps);
    
    return 1;
}

bool txt_export_lines(char *input, char *output, long st, long len, int lineLen, bool append, bool verbose) {
    if(input == NULL || output == NULL || len < 1 || st < 0) {
        if(verbose) {
            printf("Invalid parameters.\n");
            
            return false;
        }
    }
    
    FILE *f1 = fopen(input, "r");
    if(f1 == NULL) {
        if(verbose) {
            printf("Cannot read file %s.\n", input);
            
            return false;
        }
    }
    
    FILE *f2 = fopen(output, append ? "a" : "w");
    if(f2 == NULL) {
        if(verbose) {
            printf("Cannot write to file %s.\n", output);
            
            return false;
        }
    }
    
    char time[80];
    struct timeb ts;
    if(verbose) {
        ftime(&ts);
    }
    
    char line[lineLen + 3];
    bool ok = true;
    
    long li = 0;
    while(ok && ! feof(f1) && li < st + len - 1) {
        ok = fgets(line, lineLen + 3, f1) == line;
        ok = ok && strlen(line) <= lineLen + 1;
        
        li ++;
        if(ok && li >= st && li < st + len) {
            ok = fputs(line, f2) != EOF;
        }
    }
    
    fclose(f1);
    fclose(f2);
    
    if(verbose) {
        if(ok) {
            lapse(&ts, time);
            printf("%ld lines written to %s in %s\n", len, output, time);
        } else {
            printf("An IO error occurred !\n");
        }
    }
    
    return ok;
}

ulong txt_line_count(char *input, int lineLen, bool verbose) {
    if(input == NULL ) {
        if(verbose) {
            printf("Invalid parameters.\n");
            
            return false;
        }
    }
    
    FILE *f1 = fopen(input, "r");
    if(f1 == NULL) {
        if(verbose) {
            printf("Cannot read file %s.\n", input);
            
            return false;
        }
    }
    
    char time[80];
    struct timeb ts;
    if(verbose) {
        ftime(&ts);
    }
    
    char line[lineLen + 3], *fg;
    bool ok = true;
    
    long li = 0;
    while(ok && ! feof(f1)) {
        fg = fgets(line, lineLen + 3, f1);
        
        if(fg == line) { // otherwise probably EOF
            ok = ok && strlen(line) <= lineLen + 1; // line too long, could not read
            
            li ++;
        }
    }
    
    fclose(f1);
    
    if(verbose) {
        if(ok) {
            lapse(&ts, time);
            printf("The text file %s contains %ld lines, counted in %s\n", input, li, time);
        } else {
            printf("An IO error occurred !\n");
        }
    }
    
    return li;
}


// MARK: the help system and the main function

static const char* before = "This task [partially] converts .nset files to .csv files.\nIf the requested precision is more than 36 digits, results are refined to comply.\nIn this case, the period should be provided as the last optional parameter.\n\n";
static const char* after = "\nThe \"end\" position (optional parameter) is not included in the list.\n\n";

static const char *parameters[] = {
    "inputFileName",
    "outputFileName",
    "digits",
    "start",
    "count",
    "period",
    "append"
};

static const char *types[] = {
    "required",
    "required",
    "optional",
    "optional",
    "optional",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "",
    "26",
    "0",
    "0",
    "0",
    "1"
};

static const char *descriptions[] = {
    "the name of the input file",
    "the name of the output file",
    "digits after the decimal point, from 15 to 36 (37 to 1500 to refine); 0 to output a .nset file",
    "first index in the input file, inclusive",
    "the number of points to export, 0 to export the entire file",
    "the period of the hyperbolic centers is necessary to refine the results",
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

int nset2csv_main(int argc, const char *argv[]) {
    int dig = 26, per = 0, append = 1;
    long start = 0, count = 0;
    
    if(argc < 2) {
        help();
        
        return 1;
    }
    
    if(argc >= 3 && (sscanf(argv[2], "%d", &dig) < 1 || (dig < 15 && dig != 0) || dig > 1500)) {
        help();
        
        return 1;
    }
    
    if(argc >= 5 && (sscanf(argv[3], "%ld", &start) < 1 || sscanf(argv[4], "%ld", &count) < 1 ||
                     start < 0 || count < 0)) {
        help();
        
        return 1;
    }
        
    if((argc < 6 && dig > U128_DIGITS) ||
       (argc >= 6 && (sscanf(argv[5], "%d", &per) < 1 || per < 3 || per > HYP_RAW_MAX_PER))) {
        help();
        
        return 1;
    }
    
    if(argc >= 7 && (sscanf(argv[6], "%d", &append) < 1)) {
        help();
        
        return 1;
    }
    
    return ! hyp_export2csv((char *) argv[0], (char *) argv[1], dig, start, count, per, append);
}
