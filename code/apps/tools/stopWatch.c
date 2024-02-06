//
//  stopWatch.c
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

#include <stdio.h>
#include <time.h>
#include <sys/timeb.h>

#include "stopWatch.h"

void lapse(rtime *ts, char *str) {
    struct timeb te;
    ftime(&te);
    
    if(str != NULL) {
        long df = (te.time - ts->time) * 1000 + (te.millitm - ts->millitm); // time increment in milliseconds
        
        millis(df, str);
    }
    
    *ts = te;
}

long secs(void) {
    struct timeb te;
    ftime(&te);
    
    return te.time;
}

long msecs(void) {
    struct timeb te;
    ftime(&te);
    
    return te.time * 1000 + te.millitm;
}

void millis(long df, char *str) {
    if(df > 100) {    // pretty print durations that exceed 100 milliseconds
        double s = df - (df / MINUTE) * MINUTE;
        s = s / 1000;
        int ss = 2 * s + 1;
        ss = ss / 2;
        
        long dd = df / DAY;
        long hh = (df - dd * DAY) / HOUR;
        long r = df - dd * DAY - hh * HOUR;
        long mm = r / MINUTE;
        long mmm = r / MINUTE + (r % MINUTE > MINUTE / 2 ? 1 : 0);
        
        if(df < MINUTE) {
            snprintf(str, 30, "%.2lf s", s); // shorter than a minute : seconds, rounded up to 10 milliseconds
        } else if(df < HOUR) {
                snprintf(str, 30, "%ldm %.1lfs", mm, s); // minute long : rounded up to 0.1 seconds
        } else if(df < DAY) {
                snprintf(str, 30, "%ldh %ldm %ds", hh, mm, ss); // hour long : rounded up to 1 second
        } else {
            snprintf(str, 30, "%ldd %ldh %ldm", dd, hh, mmm); // day long : rounded up to 1 minute
        }
    } else {
        snprintf(str, 30, "%ld ms", df);  // raw if shorter than 100 milliseconds
    }
}

void time_stamp(char *str, int len, bool human, bool precise) {
    time_t now;
    time(&now);
    struct tm *local = localtime(&now);
    
    if(human) {
        if(precise) {
            ptime now;
            clock_gettime(CLOCK_MONOTONIC, &now);
            
            snprintf(str, 60, "%d-%d-%d %2d:%02d:%02d.%03lu", local->tm_year + 1900, local->tm_mon + 1,
                     local->tm_mday, local->tm_hour, local->tm_min, local->tm_sec, now.tv_nsec / 1000000);
        } else {
            snprintf(str, 60, "%d-%d-%d %2d:%02d:%02d", local->tm_year + 1900, local->tm_mon + 1,
                     local->tm_mday, local->tm_hour, local->tm_min, local->tm_sec);
        }
    } else {
        if(precise) {
            ptime now;
            clock_gettime(CLOCK_MONOTONIC, &now);
            
            snprintf(str, 60, "%02d%02d%02d%02d%02d%02d%03lu", local->tm_year % 100, local->tm_mon + 1,
                     local->tm_mday, local->tm_hour, local->tm_min, local->tm_sec, now.tv_nsec / 1000000);
        } else {
            snprintf(str, 60, "%02d%02d%02d%02d%02d%02d", local->tm_year % 100, local->tm_mon + 1,
                     local->tm_mday, local->tm_hour, local->tm_min, local->tm_sec);
        }
    }
}

void nanos(long ns, char *str) {
    if(ns > 100 * MILLION) {
        millis((ns + MILLION / 2) / MILLION, str);
        
        return;
    } else if(ns >= 10 * MILLION) {
        double dns = ns;
        dns /= MILLION;
        
        snprintf(str, 30, "%.1lf ms", dns);
    } else if(ns >= MILLION) {
        double dns = ns;
        dns /= MILLION;
        
        snprintf(str, 30, "%.2lf ms", dns);
    } else if(ns >= MILLION / 100) {
        double dns = ns;
        dns /= 1000;
        
        snprintf(str, 30, "%.1lf us", dns);
    } else if(ns > 1000) {
        double dns = ns;
        dns /= 1000;
        
        snprintf(str, 30, "%.2lf us", dns);
    } else if(ns == 1000) {
        snprintf(str, 30, "1 us");
    } else {
        snprintf(str, 30, "%ld ns", ns);
    }
}

long lap(clockid_t clock, ptime *ts, char *str) {
    ptime now;
    clock_gettime(clock, &now);
    
    long ns = (now.tv_sec - ts->tv_sec) * BILLION + now.tv_nsec - ts->tv_nsec;
    
    *ts = now;
    
    if(str != NULL) {
        nanos(ns, str);
    }
    
    return ns;
}

int date(char *now, int len) {
    time_t dt;
    time(&dt);
    
    int ch = snprintf(now, len, "%s", ctime(&dt));
    now[ch - 1] = 0;
    
    return ch - 1;
}
