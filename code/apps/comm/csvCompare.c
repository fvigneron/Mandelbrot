//
//  csvCompare.c
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
#include <string.h>
#include <sys/timeb.h>

#include "fp80.h"
#include "stopWatch.h"
#include "csvCompare.h"

static mpfr_t ldBuf, vl, vr, dif;
static bool ldBufIni = false, bufPrec = 0;

static mpfr_ptr getLdBuf(void) {
    if(! ldBufIni) {
        ldBufIni = true;
        
        mpfr_init2(ldBuf, 66);
    }
    
    return ldBuf;
}

static void clearLdBuf(void) {
    if(ldBufIni) {
        mpfr_clear(ldBuf);
        
        ldBufIni = false;
    }
}

static void initBuffers(int prec) {
    if(bufPrec != prec) {
        if(bufPrec == 0) {
            bufPrec = prec;
            
            mpfr_init2(vl, prec);
            mpfr_init2(vr, prec);
            mpfr_init2(dif, prec);
        } else {
            bufPrec = prec;
            
            mpfr_set_prec(vl, prec);
            mpfr_set_prec(vr, prec);
            mpfr_set_prec(dif, prec);
        }
    }
}

static void clearBuffers(void) {
    if(bufPrec > 0) {
        bufPrec = 0;
        
        mpfr_clear(vl);
        mpfr_clear(vr);
        mpfr_clear(dif);
    }
}

csvl csvl_read(FILE *f, uint prec) {
    if(f == NULL) {
        return NULL;
    }
    
    long fp = ftell(f);
    if(fp < 0) {
        return NULL;
    }
    
    uint len = 10000;
    char buf[len + 1];
    if(feof(f) || fgets(buf, len, f) != buf) {
        return NULL;
    }
    
    bool alloc = false;
    char *line = buf;
    if(strlen(line) >= len - 1) { // longer line
        bool all = false;
        
        do {
            len *= 2;
            line = malloc(len + 1);
            
            if(line == NULL) {
                return NULL;
            }
            
            alloc = true;
            
            if(fseek(f, fp, SEEK_SET) != 0 || fgets(line, len, f) != line) {
                free(line);
                
                return NULL;
            }
            
            all = strlen(line) < len - 1;
            if(! all) {
                free(line);
            }
        } while(! all && len < 1000000000);
        
        if(! all) {
            if(alloc) {
                free(line);
            }
            
            return NULL;
        }
    }
    
    // cut the line and count the numbers it may contain
    bool sep = true, done = false;
    int n = 0;
    for (int i = 0; i < len && ! done; i++) {
        char c = line[i];
        if(c == 0) {
            done = true;
            len = i;
            
            continue;
        }
        
        if(c == ' ' || c == ',' || c == '\n') {
            n += sep ? 0 : 1;
            sep = true;
            
            line[i] = 0;
        } else {
            sep = false;
        }
    }
    
    // read the numbers
    bool ld = prec <= LDBL_PREC;
    
    csvl cl = malloc(sizeof(csvLine) + n * sizeof(number_uni));
    cl->prec = prec;
    cl->count = n;
    
    if(! ld) {
        for (int i = 0; i < n; i++) {
            mpfr_init2(cl->vals[i].mp, prec);
        }
    }
    
    sep = true;
    bool ok = true;
    int i = 0, k = 0;
    for (; k < n && i < len && ok; k++) {
        while(line[i] == 0 && i < len) {
            i ++;
        }
        
        ok = ok && i < len;
        if(! ok) {
            continue;
        }
        
        // beginning of a new number
        if(ld) { // reading a long double with mpfr is quicker than scanf !
            mpfr_ptr ldb = getLdBuf();
            ok = ok && mpfr_set_str(ldb, line + i, 10, MPFR_RNDN) == 0;
            cl->vals[k].ld = mpfr_get_ld(ldb, MPFR_RNDN);
        } else {
            ok = ok && mpfr_set_str(cl->vals[k].mp, line + i, 10, MPFR_RNDN) == 0;
        }
        
        while(line[i] != 0 && i < len) {
            i ++;
        }
    }

    if(alloc) {
        free(line);
    }
    
    if(! ok || k < n) {
        csvl_free(cl);
        
        return NULL;
    }
    
    return cl;
}

bool csv_read_fp80(fp80 c, FILE *f) {
    if(c == NULL || f == NULL) {
        return false;
    }
    
    uint len = 10000;
    char buf[len + 1];
    if(feof(f) || fgets(buf, len, f) != buf || strlen(buf) >= len - 1) {
        return false;
    }
    
    // cut the line and count the numbers it may contain
    bool sep = true, done = false;
    int n = 0;
    for (int i = 0; i < len && ! done && n < 2; i++) {
        char c = buf[i];
        if(c == 0) {
            done = true;
            len = i;
            
            continue;
        }
        
        if(c == ' ' || c == ',' || c == '\n') {
            n += sep ? 0 : 1;
            sep = true;
            
            buf[i] = 0;
        } else {
            sep = false;
        }
    }
        
    ldbl xy[2];
    
    sep = true;
    bool ok = true;
    int i = 0, k = 0;
    mpfr_ptr ldb = getLdBuf();
    for (; k < n && i < len && ok; k++) {
        while(buf[i] == 0 && i < len) {
            i ++;
        }
        
        ok = ok && i < len;
        if(! ok) {
            continue;
        }
        
        // beginning of a new number, reading a long double with mpfr is quicker than scanf !
        ok = ok && mpfr_set_str(ldb, buf + i, 10, MPFR_RNDN) == 0;
        xy[k] = mpfr_get_ld(ldb, MPFR_RNDN);
        
        while(buf[i] != 0 && i < len) {
            i ++;
        }
    }
    
    if(! ok) {
        return false;
    }
    
    c->x = xy[0];
    c->y = xy[1];
    
    return true;
}

void csvl_free(csvl line) {
    if(line == NULL) {
        return;
    }
    
    if(line->prec > LDBL_PREC) {
        for (uint i = 0; i < line->count; i ++) {
            mpfr_clear(line->vals[i].mp);
        }
    }

    free(line);
}

bool csvl_eq(csvl l, csvl r, ldbl error, ldbl *maxError) {
    if(l == NULL || r == NULL || error < 0 || l->count != r->count) {
        return false;
    }
    
    uint count = l->count;
    bool eq = true;
    bool lld = l->prec <= LDBL_PREC;
    bool rld = r->prec <= LDBL_PREC;
    
    ldbl me = 0, df;
    
    if(lld || rld) {
        ldbl xl, xr;
        for (uint i = 0; i < count && eq; i++) {
            xl = lld ? l->vals[i].ld : mpfr_get_ld(l->vals[i].mp, MPFR_RNDN);
            xr = rld ? r->vals[i].ld : mpfr_get_ld(r->vals[i].mp, MPFR_RNDN);
            df = xl - xr;
            df = df < 0 ? -df : df;
            
            eq = eq && df <= error;
            me = df > me ? df : me;
        }
        
        if(eq && maxError != NULL) {
            *maxError = me;
        }
        
        return eq;
    }
    
    uint prec = l->prec > r->prec ? l->prec : r->prec;
    initBuffers(prec);
    
    for (uint i = 0; i < count && eq; i++) {
        if(lld) {
            mpfr_set_ld(vl, l->vals[i].ld, MPFR_RNDN);
        } else {
            mpfr_set(vl, l->vals[i].mp, MPFR_RNDN);
        }
        
        if(rld) {
            mpfr_set_ld(vr, r->vals[i].ld, MPFR_RNDN);
        } else {
            mpfr_set(vr, r->vals[i].mp, MPFR_RNDN);
        }
        
        mpfr_sub(dif, vl, vr, MPFR_RNDN);
        mpfr_abs(dif, dif, MPFR_RNDN);
                
        eq = eq && mpfr_cmp_ld(dif, error) <= 0;
        
        if(eq) {
            df = mpfr_get_ld(dif, MPFR_RNDN);
            me = df > me ? df : me;
        }
    }

    if(eq && maxError != NULL) {
        *maxError = me;
    }
    
    return eq;
}

static bool deleteAndAdd(csvl v[], int pos, int len, FILE *f, int prec, ulong *lineCount) {
    if(v[pos] == NULL || pos < 0 || pos >= len) {
        return false;
    }
    
    ulong line = v[pos]->line;
    csvl_free(v[pos]);
    
    for (int i = pos; i < len - 1; i++) {
        v[i] = v[i + 1];
    }
    
    v[len - 1] = csvl_read(f, prec);
    if(lineCount != NULL && v[len - 1] != NULL) {
        lineCount[0] ++;
    }
    
    if(v[len - 1] != NULL) {
        if(len == 1) {
            v[len - 1]->line = line + 1;
        } else if(v[len - 2] != NULL) {
            v[len - 1]->line = v[len - 2]->line + 1;
        } else {
            return false;
        }
    }
    
    return true;
}

static void printMissing(csvl l, char *format, char *fileName) {
    printf("%s L %ld: ", fileName, l->line + 1);
    
    bool ld = l->prec <= LDBL_PREC;
    for (uint i = 0; i < l->count; i++) {
        if(ld) {
            printf(format, l->vals[i].ld);
        } else {
            mpfr_printf(format, l->vals[i].mp);
        }
        
        if(i < l->count - 1) {
            printf(", ");
        }
    }
    
    printf("\n");
}

bool csv_compare(char *fn1, char *fn2, ldbl err, int prec, int lines, ldbl *maxError,
                 bool verbose, int diffCount) {
    if(fn1 == NULL || fn2 == NULL || err < 0 || prec < 0 || prec > 10000 || lines < 0) {
        if(verbose) {
            printf("Invalid parameters.\n");
            
            return false;
        }
    }
    
    FILE *f1 = fopen(fn1, "r");
    if(f1 == NULL) {
        if(verbose) {
            printf("Cannot read file %s.\n", fn1);
            
            return false;
        }
    }
    
    FILE *f2 = fopen(fn2, "r");
    if(f2 == NULL) {
        if(verbose) {
            printf("Cannot read file %s.\n", fn2);
            
            return false;
        }
    }
    
    char time[80];
    struct timeb ts;
    
    if(verbose) {
        ftime(&ts);
        
        printf("Comparing the following CSV files, up to error %Lg using precision %d bits\n", err, prec);
        printf("F1: %s\n", fn1);
        printf("F2: %s\n\n", fn2);
    }
    
    ldbl me = 0, le;
    ulong l1 = 0, l2 = 0;
    csvl ll[lines + 1], rl[lines + 1];
    bool ok = true;
    for (int i = 0; i <= lines && ok; i++) {
        ll[i] = csvl_read(f1, prec);
        if(ll[i] != NULL) {
            ll[i]->line = i;
            l1 ++;
            
            ok = ok && (i == 0 || ll[i - 1] != NULL);
        }
        
        rl[i] = csvl_read(f2, prec);
        if(rl[i] != NULL) {
            rl[i]->line = i;
            l2 ++;
            
            ok = ok && (i == 0 || rl[i - 1] != NULL);
        }
    }
       
    ok = ok && ll[0] != NULL && rl[0] != NULL;
    
    bool ld = prec <= LDBL_PREC;
    uint digits = ld ? 19 : prec / 3.2;
    char format[20];
    snprintf(format, 19, ld ? "%%.%dLf" : "%%.%dRf", digits);
    
    bool eq = true;
    int diffs = 0;
    
    while(ok && ll[0] != NULL && rl[0] != NULL && (verbose || eq)) {
        if(csvl_eq(ll[0], rl[0], err, &le)) {
            ok = ok && deleteAndAdd(ll, 0, lines + 1, f1, prec, &l1);
            ok = ok && deleteAndAdd(rl, 0, lines + 1, f2, prec, &l2);
            me = le > me ? le : me;
            
            continue;
        }
        
        bool foundLeft = false;
        for (int i = 1; i <= lines && ! foundLeft; i++) {
            if(csvl_eq(ll[0], rl[i], err, &le)) {
                ok = ok && deleteAndAdd(ll, 0, lines + 1, f1, prec, &l1);
                ok = ok && deleteAndAdd(rl, i, lines + 1, f2, prec, &l2);
                me = le > me ? le : me;
                
                foundLeft = true;
            }
        }
        
        bool foundRight = false;
        for (int i = foundLeft ? 0 : 1; i <= lines && ! foundRight; i++) {
            if(csvl_eq(ll[i], rl[0], err, &le)) {
                ok = ok && deleteAndAdd(ll, i, lines + 1, f1, prec, &l1);
                ok = ok && deleteAndAdd(rl, 0, lines + 1, f2, prec, &l2);
                me = le > me ? le : me;
                
                foundRight = true;
            }
        }

        if(foundLeft && foundRight) {
            continue;
        }
        
        // none of ll[0] or rl[0] has a pair
        if(verbose) {
            if(! foundLeft) {
                printMissing(ll[0], format, "F1");
            }
            
            if(! foundRight) {
                printMissing(rl[0], format, "F2");
            }
            
            printf("\n");
        }
        
        ok = ok && (foundLeft || deleteAndAdd(ll, 0, lines + 1, f1, prec, &l1));
        ok = ok && (foundRight || deleteAndAdd(rl, 0, lines + 1, f2, prec, &l2));
        
        eq = false;
        diffs ++;
        
        ok = ok && diffs < diffCount;
    }
    
    // if some lines are left without pair, the files are not equal
    eq = eq && ll[0] == NULL && rl[0] == NULL;
    
    bool print = false;
    while(ok && verbose && ll[0] != NULL) {
        printMissing(ll[0], format, "F1");
        ok = ok && deleteAndAdd(ll, 0, lines + 1, f1, prec, &l1);
        diffs ++;
        
        ok = ok && diffs < diffCount;
        print = true;
    }
    
    if(print) {
        printf("\n");
    }
    
    print = false;
    while(ok && verbose && rl[0] != NULL) {
        printMissing(rl[0], format, "F2");
        ok = ok && deleteAndAdd(rl, 0, lines + 1, f2, prec, &l2);
        diffs ++;
        
        ok = ok && diffs < diffCount;
        print = true;
    }
    
    if(print) {
        printf("\n");
    }
    
    fclose(f1);
    fclose(f2);
    
    // if ! ok, do some housekeeping
    for (int i = 0; i <= lines; i++) {
        if(ll[i] != NULL) {
            free(ll[i]);
        }
        
        if(rl[i] != NULL) {
            free(rl[i]);
        }
    }
    
    if(verbose) {
        if(ok) {
            if(eq) {
                printf("The files contain the same numbers, max error %Lg.\n", me);
            } else {
                if(diffs == 1) {
                    printf("There is one difference between the two files.\n");
                } else {
                    printf("There are %d differences between the two files.\n", diffs);
                }
                
            }
        } else {
            if(diffs == diffCount) {
                printf("There are at least %d differences between the two files !\n", diffs);
            } else {
                printf("\nIO errors have occurred, the files comparison failed !\n");
            }
        }
        
        if(l1 == l2) {
            printf("Both files contain %lu lines.\n", l1);
        } else {
            printf("The files contain %lu and respectively %lu lines.\n", l1, l2);
        }
        
        lapse(&ts, time);
        printf("The comparison took %s\n\n", time);
    }
    
    if(maxError != NULL) {
        *maxError = me;
    }
    
    clearLdBuf();
    clearBuffers();
    
    return ok && eq;
}

// MARK: the help system and the main function

static const char* before = "This app compares two CSV files containing numerical values.\n\n";
static const char* after = "\nValues that are closer than #error are considered equal.\n\n";

static const char *parameters[] = {
    "inputFileName",
    "outputFileName",
    "error",
    "precision",
    "lines",
    "max errors"
};

static const char *types[] = {
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
    "1E-17",
    "64",
    "10",
    "100"
};

static const char *descriptions[] = {
    "the name of the input file",
    "the name of the output file",
    "the max error to consider two numbers equal",
    "the precision in bits, 64 to use long double, more to use mpfr",
    "the max absolute distance of out of order lines to compare",
    "the max number of reported errors"
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

int csv_compare_main(int argc, const char *argv[]) {
    ldbl err = 1E-17;
    int prec = 64, lines = 10, errCount = 100;
    
    if(argc < 2) {
        help();
        
        return 1;
    }
    
    if(argc >= 3 && (sscanf(argv[2], "%Lf", &err) < 1 || err < 0 || err > 1)) {
        help();
        
        return 1;
    }
    
    if(argc >= 4 && (sscanf(argv[3], "%d", &prec) < 1 || prec > 10000 || prec < 0)) {
        help();
        
        return 1;
    }
        
    prec = prec < 64 ? 64 : prec;
    
    if(argc >= 5 && (sscanf(argv[4], "%d", &lines) < 1 || lines > 1000 || lines < 0)) {
        help();
        
        return 1;
    }
    
    if(argc >= 5 && (sscanf(argv[4], "%d", &lines) < 1 || lines > 1000 || lines < 0)) {
        help();
        
        return 1;
    }
    
    if(argc >= 6 && (sscanf(argv[5], "%d", &errCount) < 1 || errCount > 10000 || errCount < 1)) {
        help();
        
        return 1;
    }
    
    return ! csv_compare((char *) argv[0], (char *) argv[1], err, prec, lines, NULL, true, errCount);
}
