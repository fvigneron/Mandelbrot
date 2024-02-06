//
//  arithmetic.c
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2022.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2022.
//
//  Copyright © 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the Lesser GNU Public Licence v3.0
//

#include "arithmetic.h"

// MARK: constants and buffers definitions

#define SMALL_PRIMES_LEN 100
static uint small_primes[SMALL_PRIMES_LEN];
static int primesInited = 0;

#define SMALL_MU_LEN 100
static uint small_mu[SMALL_MU_LEN];
static int muInited = 0;

// MARK: functions

/// Inits the buffer with the first @c SMALL_PRIMES_LEN prime numbers.
static void initPrimes(void) {
    int p = 2, i = 0;
    
    do {
        int pp = 1;
        for (int j = 2; pp && j * j <= p; j++) {
            pp &= p % j != 0;
        }
        
        if(pp) { // p is prime
            small_primes[i ++] = p;
        }
        
        p ++;
    } while(i < SMALL_PRIMES_LEN);
    
    primesInited = 1;
}

ulong factor(ulong n) {
    if(n < 4 || n == 5 || n == 7 || n == 11 || n == 13) {
        return n;
    }
    
    if(n % 2 == 0) {
        return 2;
    }
    
    if(n % 3 == 0) {
        return 3;
    }
    
    if(n % 5 == 0) {
        return 5;
    }
    
    // here n > 16
    ulong sqp = ceill(sqrtl(n));
    if(n % sqp == 0) {
        return factor(sqp);
    }
    
    sqp --;
    if(n % sqp == 0) { // just to protect from rounding quirks
        return factor(sqp);
    }
    
    if(! primesInited) {
        initPrimes(); // not optimal, but short and executed only once
    }
    
    ulong d = 7;
    for (int i = 3; i < SMALL_PRIMES_LEN && d < sqp; i ++) {
        d = small_primes[i];
        if(n % d == 0) {
            return d;
        }
    }
    
    if(d >= sqp) {
        return n;
    }
    
    d = 6 * (d / 6) + 1;
    for (; d < sqp; d += 6) {
        if(n % d == 0) {
            return d;
        }
        
        if(n % (d + 4) == 0) {
            return d + 4;
        }
    }
    
    return n;
}

bool is_prime(ulong p) {
    if(p < 2) {
        return false;
    }
        
    return factor(p) == p;
}

/// Inits the buffer with the first @c SMALL_MU_LEN values of the function @c mu()
static void initMu(void) {
    small_mu[0] = 0;
    small_mu[1] = 1;
    
    ulong d;
    for (int i = 2; i < SMALL_MU_LEN; i++) {
        d = factor(i);
        if(i == d) {
            small_mu[i] = -1; // prime
        } else {
            small_mu[i] = i % (d * d) == 0 ? 0 : -small_mu[i / d];
        }
    }
    
    muInited = 1;
}

int mu(ulong n) {
    if(n == 1) {
        return 1;
    }
    
    if(n % 4 == 0 || n % 9 == 0 || n % 25 == 0) {
        return 0;
    }
    
    if(n % 2 == 0) {
        return -mu(n / 2);
    }
    
    if(n % 3 == 0) {
        return -mu(n / 3);
    }
    
    if(n % 5 == 0) {
        return -mu(n / 5);
    }
    
    if(! muInited) {
        initMu();
    }
    
    if(n < SMALL_MU_LEN) {
        return small_mu[n];
    }
    
    ulong d = factor(n);
    if(d == n) {
        return -1;
    }
    
    if(n % (d * d) == 0) {
        return 0;
    }
    
    return -mu(n / d);
}

void print_mu(int count) {
    for (int i = 1; i < count; i++) {
        printf(i % 100 == 0 ? "%3d,\n\n" : i % 10 == 0 ? "%3d,\n" : "%3d,", mu(i));
    }
    
    printf("%3d\n", mu(count));
}

void print_phi(int count) {
    for (int i = 1; i < count; i++) {
        printf(i % 100 == 0 ? "%6lu,\n\n" : i % 10 == 0 ? "%6lu,\n" : "%6lu,", phi(i));
    }
    
    printf("%6lu\n", phi(count));
}

/// Formatted printf for the number @c p.
///
/// @param p the number
/// @param ln the index
/// @param last @c 1 if it is the last in the list, @c 0 otherwise
static void pp(ulong p, int ln, int last) {
    if(last) {
        printf("%6lu\n", p);
        
        return;
    }

    printf(ln % 100 == 0 ? "%6lu,\n\n" : ln % 10 == 0 ? "%6lu,\n" : "%6lu,", p);
}

void print_primes(int count) {
    ulong p = 2;
    for (int i = 0; i < count; i++) {
        while(! is_prime(p)) { // not optimal, bau printf dominates
            p ++;
        }
        
        pp(p ++, i + 1, i + 1 == count);
    }
}

ulong phi(ulong n) {
    if(n < 2) {
        return n;
    }
    
    ulong p = 1, d = 2;
    if((n & 1) == 0) { // 2 is a divisor
        n >>= 1;
        while((n & 1) == 0) {
            d <<= 1;
            n >>= 1;
        }
        
        p = d - (d >> 1);
    }
    
    if(! primesInited) {
        initPrimes();
    }
    
    d = 1;
    for (int i = 1; i < SMALL_PRIMES_LEN && n > d * d; i++) {
        d = small_primes[i];
        if(n % d != 0) {
            continue;
        }
        
        ulong dk = d;
        n /= d;
        
        while(n % d == 0) {
            dk *= d;
            n /= d;
        }
        
        p *= dk - dk / d;
    }
    
    if(n == 1) {
        return p;
    }
    
    if(n <= d * d) { // n prime
        return p * (n - 1);
    }
    
    // search for larger divisors (the smallest one will be prime)
    d = 6 * (d / 6) + 1;
    ulong dk;
    while(d * d < n) { // primes > 3 are only of forms 6k + 1 and 6k + 5
        if(n % d == 0) {
            dk = d;
            n /= d;
            
            while(n % d == 0) {
                dk *= d;
                n /= d;
            }
            
            p *= dk - dk / d;
        }
        
        d += 4;
        if(n % d == 0) {
            dk = d;
            n /= d;
            
            while(n % d == 0) {
                dk *= d;
                n /= d;
            }
            
            p *= dk - dk / d;
        }
        
        d += 2;
    }
    
    if(n == 1) {
        return p;
    }
    
    // n prime
    return p * (n - 1);
}

/// Assuming @c a>b>a, it is a quick version of @c gcd().
///
/// @param a a number
/// @param b another number, @c a>b>1
static ulong fgcd(ulong a, ulong b) {
    // here a > b > 1
    ulong x = b, t = 0;
    ulong r = a % b;
    while(r > 0) {
        t = r;
        r = x % r;
        x = t;
    }
    
    return x;
}

ulong gcd(ulong a, ulong b) {
    if(a * b == 0) {
        return a + b;
    }
    
    if(a == 1 || b == 1) {
        return 1;
    }
    
    if(a == b) {
        return a;
    }
    
    return fgcd(a < b ? b : a, a < b ? a : b);
}
