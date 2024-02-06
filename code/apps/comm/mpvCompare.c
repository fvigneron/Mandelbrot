//
//  mpvCompare.c
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
#include "stopWatch.h"

#include "io.h"
#include "mpvCompare.h"

// MARK: compare functions

static mpfr_t maxErr, nextErr, ulpMax, ulpNext, x, y, m, u;
static mpc a, b;
static bool bufInit = false;
static long posMax = -1, posNext = -1, ulpPos = -1, ulpNextPos = -1;

static void mpvc_init_buff(int prec) {
    if(bufInit) {
        if(prec > mpfr_get_prec(maxErr)) {
            mpfr_set_prec(maxErr, prec);
            mpfr_set_prec(nextErr, prec);
            mpfr_set_prec(ulpMax, prec);
            mpfr_set_prec(ulpNext, prec);
            mpfr_set_prec(x, prec);
            mpfr_set_prec(y, prec);
            mpfr_set_prec(m, prec);
            mpfr_set_prec(u, prec);
            mpc_set_prec(a, prec);
            mpc_set_prec(b, prec);
        }
        
        return;
    }
    
    mpfr_inits2(prec, maxErr, nextErr, ulpMax, ulpNext, x, y, m, u, NULL);
    mpfr_set_zero(maxErr, 1);
    mpfr_set_zero(nextErr, 1);
    mpfr_set_zero(ulpMax, 1);
    mpfr_set_zero(ulpNext, 1);
    
    mpc_inits(prec, a, b, NULL);
    
    bufInit = true;
}

static void mpvc_clear_buff(void) {
    mpfr_clears(maxErr, nextErr, ulpMax, ulpNext, x, y, m, u, NULL);
    mpc_clears(a, b, NULL);
    
    bufInit = false;
}

// called when the two numbers to compare are x and y
static void mpvc_comp(long i) {
    long prec = mpfr_get_prec(x);
    
    mpfr_sub(m, x, y, MPFR_RNDA);
    mpfr_abs(m, m, MPFR_RNDU);
    
    if(mpfr_zero_p(m)) {
        return;
    }
    
    if(mpfr_cmp(m, maxErr) > 0) {
        posNext = posMax;
        mpfr_set(nextErr, maxErr, MPFR_RNDU);
        
        posMax = i;
        mpfr_set(maxErr, m, MPFR_RNDU);
    } else if(mpfr_cmp(m, nextErr) > 0) {
        posNext = i;
        mpfr_set(nextErr, m, MPFR_RNDU);
    }
    long ex = mpfr_get_exp(x);
    long ey = mpfr_get_exp(y);
    long sc = ex < ey ? ex : ey;
    mpfr_div_2si(m, m, sc - prec, MPFR_RNDU);
    
    if(mpfr_cmp(m, ulpMax) > 0) {
        ulpNextPos = ulpPos;
        mpfr_set(ulpNext, ulpMax, MPFR_RNDU);
        
        ulpPos = i;
        mpfr_set(ulpMax, m, MPFR_RNDN);
    } else if(mpfr_cmp(m, ulpNext) > 0) {
        ulpNextPos = i;
        mpfr_set(ulpNext, m, MPFR_RNDN);
    }
}

// called when the two numbers to compare are a and b
static void mpvc_compc(long i) {
    long prec = mpfr_get_prec(a->x);
    
    mpfr_sub(m, a->x, b->x, MPFR_RNDA);
    mpfr_abs(m, m, MPFR_RNDU);
    
    mpfr_sub(u, a->y, b->y, MPFR_RNDA);
    mpfr_abs(u, u, MPFR_RNDU);
    
    mpfr_max(m, m, u, MPFR_RNDU);
    if(mpfr_zero_p(m)) {
        return;
    }
    
    mpfr_mul_2si(m, m, 1, MPFR_RNDU);
    
    long ex = mpfr_get_exp(a->x);
    long ey = mpfr_get_exp(a->y);
    long ea = ex < ey ? ey : ex;
    
    ex = mpfr_get_exp(b->x);
    ey = mpfr_get_exp(b->y);
    long eb = ex < ey ? ey : ex;
    
    long sc = ea < eb ? ea : eb;
    
    // dist(a, b) < m here
    if(mpfr_cmp(m, nextErr) > 0) {
        mpc_dist(m, a, b);
        
        if(mpfr_cmp(m, maxErr) > 0) {
            posNext = posMax;
            mpfr_set(nextErr, maxErr, MPFR_RNDN);
            
            posMax = i;
            mpfr_set(maxErr, m, MPFR_RNDN);
        } else if(mpfr_cmp(m, nextErr) > 0) {
            posNext = i;
            mpfr_set(nextErr, m, MPFR_RNDN);
        }
        
        mpfr_div_2si(m, m, sc - prec, MPFR_RNDU);
        
        if(mpfr_cmp(m, ulpMax) > 0) {
            ulpNextPos = ulpPos;
            mpfr_set(ulpNext, ulpMax, MPFR_RNDU);
            
            ulpPos = i;
            mpfr_set(ulpMax, m, MPFR_RNDN);
        } else if(mpfr_cmp(m, ulpNext) > 0) {
            ulpNextPos = i;
            mpfr_set(ulpNext, m, MPFR_RNDN);
        }
        
        return;
    }
    
    mpfr_div_2si(m, m, sc - prec, MPFR_RNDU);
    
    if(mpfr_cmp(m, nextErr) > 0) {
        mpc_dist(m, a, b);
        mpfr_div_2si(m, m, sc - prec, MPFR_RNDU);
        
        if(mpfr_cmp(m, ulpMax) > 0) {
            ulpNextPos = ulpPos;
            mpfr_set(ulpNext, ulpMax, MPFR_RNDU);
            
            ulpPos = i;
            mpfr_set(ulpMax, m, MPFR_RNDN);
        } else if(mpfr_cmp(m, ulpNext) > 0) {
            ulpNextPos = i;
            mpfr_set(ulpNext, m, MPFR_RNDN);
        }
    }
}

static bool mpv_comp_complex(mpv s, long sst, mpv d, long dst, long count, long done) {
    if(sst < 0 || dst < 0 || count <= 0 || done < 0) {
        return false;
    }
    
    int prec = s->prec > d->prec ? s->prec : d->prec;
    mpvc_init_buff(prec);
            
    long sc = s->count - (sst << 1);
    long dc = d->count - (dst << 1);
    
    long len = sc > dc ? dc : sc;
    
    if(len <= 0 || len > (count << 1)) {
        return false;
    }
    
    len /= 2;
    for (long i = 0; i < len; i++) {
        mpv_getc(a, s, i + sst);
        mpv_getc(b, d, i + dst);
        
        mpvc_compc(i + done);
    }
    
    return true;
}

static bool mpv_comp_real(mpv s, long sst, mpv d, long dst, long count, long done) {
    if(sst < 0 || dst < 0 || count <= 0 || done < 0) {
        return false;
    }
    
    int prec = s->prec > d->prec ? s->prec : d->prec;
    mpvc_init_buff(prec);
        
    long sc = s->count - sst;
    long dc = d->count - dst;
    
    long len = sc > dc ? dc : sc;
    
    if(len <= 0 || len > count) {
        return false;
    }
    
    for (long i = 0; i < len; i++) {
        mpv_get(x, s, i + sst);
        mpv_get(y, d, i + dst);
        
        mpvc_comp(i + done);
    }
        
    return true;
}

static bool mpv_compare(char * src, char *dst, bool complex) {
    char *type = complex ? "complex" : "real";
    int cr = complex ? 2 : 1;
    
    // welcome message
    printf("Comparing \n\t%s \nto \n\t%s, \nas vectors of %s numbers.\n\n", src, dst, type);
    fflush(stdout);
    
    char time[100];
    struct timeb ats, ts;
    ftime(&ats);
    
    ftime(&ts);
    
    FILE *f = fopen(src, "r");
    FILE *g = fopen(dst, "r");
    bool ok = f != NULL && g != NULL;
    if(! ok) {
        if(f == NULL) {
            printf("Could not open %s.\n", src);
        }
        if(g == NULL) {
            printf("Could not open %s.\n", dst);
        }
        
        return false;
    }
    
    long fs = io_file_size(f);
    long gs = io_file_size(g);
    ok = ok && (fs > HEADER_MIN_LEN) && (gs > HEADER_MIN_LEN);
    
    long compared = 0, srcLoaded = 0, dstLoaded = 0; // counts in real values
    
    mpv s = NULL, d = NULL;
    while(ok) {
        bool loads = srcLoaded <= compared;
        if(loads) {
            mpv_free(s);
            s = mpv_read_from(f, -1);
            lapse(&ts, time);
            if(s == NULL) {
                ok = false;
                
                if(compared == 0) {
                    printf("Could not read any points from %s.\n", src);
                } else {
                    long count = srcLoaded / cr;
                    char *msgt = count == 1 ? complex ? "point" : "number" :
                                              complex ? "points" : "numbers";
                    printf("A total of %ld %s have been read from %s.\n", count, msgt, src);
                    
                    long fp = ftell(f);
                    if(fp < fs) {
                        printf("The file %s is corrupted, there are %ld bytes that could not be read.\n",
                               src, fs - fp);
                    }
                }
            } else {
                long count = s->count / cr;
                char *msgt = count == 1 ? complex ? "point" : "number" :
                                          complex ? "points" : "numbers";
                printf("%ld %s loaded from %s in %s.\n", count, msgt, src, time);
                srcLoaded += s->count;
                
                if(s->count % cr != 0) {
                    printf("The file %s is corrupted, %ld real coordinates have been read.\n",
                           src, s->count);
                    
                    ok = false;
                }
            }
            fflush(stdout);
        }
        
        bool loadd = dstLoaded <= compared;
        if(loadd && srcLoaded > compared) {
            mpv_free(d);
            d = mpv_read_from(g, -1);
            lapse(&ts, time);
            if(d == NULL || d->count / cr == 0) {
                ok = false;
                
                if(compared == 0) {
                    printf("Could not read any points from %s.\n", src);
                } else {
                    long count = dstLoaded / cr;
                    char *msgt = count == 1 ? complex ? "point" : "number" :
                                              complex ? "points" : "numbers";
                    printf("A total of %ld %s have been read from %s.\n", count, msgt, dst);
                    
                    long gp = ftell(g);
                    if(gp < gs) {
                        printf("The file %s is probably corrupted, there are %ld bytes that could not be read.\n",
                               dst, gs - gp);
                    }
                }
            } else {
                long count = d->count / cr;
                char *msgt = count == 1 ? complex ? "point" : "number" :
                                          complex ? "points" : "numbers";
                printf("%ld %s loaded from %s in %s.\n", count, msgt, dst, time);
                dstLoaded += d->count;
                
                if(d->count % cr != 0) {
                    printf("The file %s is corrupted, %ld real coordinates have been read.\n",
                           dst, d->count);
                    
                    ok = false;
                }
            }
            fflush(stdout);
        }
        
        // all previously read points have been compared; if new points cannot be read, done !
        if(! ok) {
            continue;
        }
        
        long sts = loads ? 0 : s->count - (srcLoaded - compared);
        long std = loadd ? 0 : d->count - (dstLoaded - compared);
        long ls = s->count - sts;
        long ld = d->count - std;
        long count = ls <= ld ? ls : ld;
        count /= cr;
        
        if(complex) {
            ok = mpv_comp_complex(s, sts, d, std, count, compared);
        } else {
            ok = mpv_comp_real(s, sts, d, std, count, compared);
        }
        compared += ok ? count * cr : 0;
        
        lapse(&ts, time);
        if(ok) {
            char *msgt = count == 1 ? complex ? "point" : "number" :
                                      complex ? "points" : "numbers";
            printf("%ld %s compared in %s (total %ld).\n\n", count, msgt, time, compared / cr);
        } else {
            printf("Internal error, sorry ! sts = %ld, std = %ld, count = %ld\n", sts, std, count);
        }
        fflush(stdout);
    }
    
    if(mpfr_zero_p(maxErr)) {
        printf("\nThe two vectors are identical up to this position.\n");
    } else {
        char *msgt = complex ? "distance" : "difference";
        mpfr_printf("\nThe largest %s is %Rg, on position %ld.\n", msgt, maxErr, posMax);
        mpfr_printf("The second largest %s is %Rg, on position %ld.\n", msgt, nextErr, posNext);
        mpfr_printf("\nThe largest relative %s is %Rg ulp, on position %ld.\n",
                    msgt, ulpMax, ulpPos);
        mpfr_printf("The second largest relative %s is %Rg ulp, on position %ld.\n",
                    msgt, ulpNext, ulpNextPos);
    }
        
    mpv_free(s);
    mpv_free(d);
    mpvc_clear_buff();
    
    lapse(&ats, time);
    printf("\nAll operations completed in %s.\n\n", time);
    
    return true;
}

// MARK: the help system and the main function

static const char* before = "This task compares two .mpv files that may be assembled as several mini-files.\n\n";
static const char* after = "\nThe numbers may be seen as real or complex. The largest distance and its position will be reported.\n\n";

static const char *parameters[] = {
    "1stFileName",
    "2ndFileName",
    "complex"
};

static const char *types[] = {
    "required",
    "required",
    "optional"
};

static const char *defaults[] = {
    "",
    "",
    "0"
};

static const char *descriptions[] = {
    "the name of the input file",
    "the name of the output file",
    "0 for a vector of real numbers, 1 for a vector of complex numbers"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 3;
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

int mpv_compare_main(int argc, const char * argv[]) {
    int complex = 0;
    
    if(argc < 2 || (argc >= 3 && sscanf(argv[2], "%d", &complex) < 1)) {
        help();
        
        return 1;
    }
    
    return ! mpv_compare((char *) argv[0], (char *) argv[1], complex != 0);
}
