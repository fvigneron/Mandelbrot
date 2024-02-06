//
//  polynomial.c
//  Mandelbrot
//
//  Created by MIHALACHE Nicolae on 12/9/20.
//  Copyright © 2020 UPEC. All rights reserved.
//

#include <stdlib.h>
#include "polynomial.h"

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Constructors
// //////////////////////////////////////////////////////////////////////////////////////////

static int poly_init2(poly p, int deg, int prec) {
    if(deg < 0 || p == NULL) {
        return 0;
    }
    
    p->deg = deg;
    p->prec = prec;
    
    p->a = malloc((deg + 1) * sizeof(mpc_struct));
    p->x = malloc((deg) * sizeof(mpc_struct));
       
    mpc_init(p->a, prec);
    mpc_init(p->s, prec);
    mpc_init(p->t, prec);
    mpc_init(p->u, prec);
    mpc_init(p->n, prec);
    
    for (int i = 1; i <= deg; i++) {
        mpc_init(p->a + i, prec);
        mpc_init(p->x + (i - 1), prec);
    }
    
    return 1;
}

int poly_init(poly p, int deg, mpc_struct *a, int prec) {
    if(! poly_init2(p, deg, prec)) {
        return 0;
    }
    
    for (int i = 0; i <= deg; i++) {
        mpc_set(p->a + i, a + i);
    }
    
    return 1;
}

int poly_init_vect(poly p, int deg, mpv v, int prec) {
    if(v == NULL || v->count < deg * 2 + 2 || ! poly_init2(p, deg, prec)) {
        return 0;
    }
    
    for (int i = 0; i <= deg; i++) {
        mpv_getc(p->a + i, v, i);
    }
    
    return 1;
}

int poly_init_real(poly p, int deg, __mpfr_struct *a, int prec) {
   if(! poly_init2(p, deg, prec)) {
        return 0;
    }
        
    for (int i = 0; i <= deg; i++) {
        mpc_setr(p->a + i, a + i);
    }
    
    return 1;
}

int poly_init_root_powers(poly p, int deg, mpc_struct *pk, int prec) {
    // pk : sum of the powers of the roots, starting with power 1
    // p : polynomal that admits those roots
    if(! poly_init2(p, deg, prec)) {
        return 0;
    }
    
    mpc_setl(p->a + deg, 1, 0);
    
    for (int k = 1; k <= deg; k++) {
        mpc_setl(p->s, 0, 0);
        
        for (int i = 1; i <= k; i++) {
            mpc_mul(p->t, &pk[i - 1], &p->a[deg - k + i]);
            mpc_add(p->s, p->s, p->t);
        }
        
        mpc_divi(&p->a[deg - k], p->s, -k);
    }
    
    return 1;
}

bool poly_init_root_powers_real(poly p, int deg, __mpfr_struct *pk, int prec, ldbl eps) {
    if(! poly_init2(p, deg, prec)) {
        return false;
    }
    
    mpc_seti(p->a + deg, 1, 0);
    int firstNonZero = deg;
    
    for (int k = 1; k <= deg; k++) {
        mpfr_set_zero(p->s->x, 1); // s - sum
        
        for (int i = 1; i <= k; i++) {
            mpfr_mul(p->s->y, &pk[i - 1], p->a[deg - k + i].x, MPFR_RNDN);
            mpfr_add(p->s->x, p->s->x, p->s->y, MPFR_RNDN);
        }
        
        mpfr_div_si(p->a[deg - k].x, p->s->x, -k, MPFR_RNDN);
        mpfr_set_zero(p->a[deg - k].y, 1);
        
        if(mpfr_cmp_ld(p->a[deg - k].x, eps) > 0 || mpfr_cmp_ld(p->a[deg - k].x, -eps) < 0) {
            firstNonZero = deg - k;
        }
    }
    
    if(firstNonZero > 0) {
        p->deg -= firstNonZero;
        for (int i = 0; i <= p->deg; i++) {
            mpfr_set(p->a[i].x, p->a[i + firstNonZero].x, MPFR_RNDN);
        }
        
        for (int i = p->deg + 1; i <= deg; i++) {
            mpfr_set_zero(p->a[i].x, 1);
        }
    }
    
    return true;
}

int poly_clear(poly p) {
    if(p == NULL || p->deg == -1 || p->a == NULL || p->x == NULL) {
        return 0;
    }
    
    mpc_clear(p->a);
    mpc_clear(p->s);
    mpc_clear(p->t);
    mpc_clear(p->u);
    mpc_clear(p->n);
    
    for (int i = 1; i <= p->deg; i++) {
        mpc_clear(p->a + i);
        mpc_clear(p->x + (i - 1));
    }
    
    p->deg = -1;
    
    free(p->a);
    p->a = NULL;
    
    free(p->x);
    p->x = NULL;
    
    return 1;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Evaluations
// //////////////////////////////////////////////////////////////////////////////////////////

int poly_val(poly p, mpc v, mpc z) {
    if(p->deg < 0 || p == NULL) {
        return 0;
    }
    
    if(p->deg == 0) { // constant polynomial
        mpc_set(v, p->a);
        
        return 1;
    }
    
    mpc_set(p->x, z);
    mpc_set(p->s, p->a); // p->s == a[0]
    
    mpc_mul(p->t, p->x, p->a + 1);
    mpc_add(p->s, p->s, p->t); // p->s == a[0] + a[1]*z
    
    for (int i = 1; i < p->deg; i++) { // compute power i + 1
        if((i + 1) % 2 == 0) {
            mpc_sqr(p->x + i, p->x + (i / 2)); // because x[i] = z^{i + 1}
        } else {
            mpc_mul(p->x + i, z, p->x + (i - 1));
        }

        mpc_mul(p->t, p->x + i, p->a + (i + 1));
        mpc_add(p->s, p->s, p->t);
    }
    
    mpc_set(v, p->s);
    
    return 1;
}

int poly_der(poly p, mpc d, mpc z) {
    if(p->deg < 0 || p == NULL) {
        return 0;
    }
    
    if(p->deg == 0) { // constant polynomial
        mpc_set0(d);
        
        return 1;
    }
    
    mpc_set(p->x, z);
    mpc_set(p->u, p->a + 1); // p->u == a[1]
    
    for (int i = 1; i < p->deg; i++) { // compute power i + 1
        if((i + 1) % 2 == 0) {
            mpc_sqr(p->x + i, p->x + (i / 2)); // because x[i] = z^{i + 1}
        } else {
            mpc_mul(p->x + i, z, p->x + (i - 1));
        }
                
        mpc_mul(p->t, p->x + (i - 1), p->a + (i + 1));
        mpc_muli(p->t, p->t, i + 1);
        mpc_add(p->u, p->u, p->t);
    }
    
    mpc_set(d, p->u);
    
    return 1;
}

int poly_val_der(poly p, mpc v, mpc d, mpc z) {
    if(p->deg < 0 || p == NULL) {
        return 0;
    }
    
    if(p->deg == 0) {
        mpc_set(v, p->a);
        mpc_set0(d);
        
        return 1;
    }
    
    mpc_set(p->x, z);
    
    mpc_set(p->s, p->a); // p->s == a[0]

    mpc_mul(p->t, p->x, p->a + 1);
    mpc_add(p->s, p->s, p->t); // p->s == a[0] + a[1]*z
    
    mpc_set(p->u, p->a + 1); // p->u == a[1]
    
    for (int i = 1; i < p->deg; i++) { // compute power i + 1
        if((i + 1) % 2 == 0) {
            mpc_sqr(p->x + i, p->x + (i / 2)); // because x[i] = z^{i + 1}
        } else {
            mpc_mul(p->x + i, z, p->x + (i - 1));
        }
        
        mpc_mul(p->t, p->x + i, p->a + (i + 1));
        mpc_add(p->s, p->s, p->t);
        
        mpc_mul(p->t, p->x + (i - 1), p->a + (i + 1));
        mpc_muli(p->t, p->t, i + 1);
        mpc_add(p->u, p->u, p->t);
    }
    
    mpc_set(v, p->s);
    mpc_set(d, p->u);
    
    return 1;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Root finding and level sets using Newton's method
// //////////////////////////////////////////////////////////////////////////////////////////

int poly_newton_term(poly p, mpc nt, mpc z) {
    int ok = poly_val_der(p, p->s, p->u, z);
    
    return ok && mpc_div(nt, p->s, p->u);
}

int poly_newton_term_sol(poly p, mpc nt, mpc z, mpc t) {
    int ok = poly_val_der(p, p->s, p->u, z);
    mpc_sub(p->s, p->s, t);
    
    return ok && mpc_div(nt, p->s, p->u);
}

bool poly_root(poly p, mpc r, mpc z, int iter, long double err, int refine) {
    mpc iz;
    mpc_init(iz, p->prec);
    mpc_set(iz, z);
    
    int div = 0;
    int conv = 0;
    
    for (int i = 0; i < iter && ! div && ! conv; i++) {
        poly_newton_term(p, p->n, iz);
        mpc_sub(iz, iz, p->n);
        
        // divergence if z escapes too far
        ldbl mod = mpc_modl(iz);
        div = mod >= 1E500L; // FIXME: really ?
        
        // convergence if value is small enough
        mod = mpc_modl(p->n);
        conv = mod < err;
    }
    
    if(div || ! conv) {
        mpc_clear(iz);
        
        return false;
    }
    
    for (int i = 0; i < refine; i++) {
        poly_newton_term(p, p->n, iz);
        mpc_sub(iz, iz, p->n);
    }
    
    mpc_set(r, iz);
    mpc_clear(iz);
    
    return true;
}

mpv poly_roots(poly p, ldbl r, ldbl eps) {
    int n = p->deg;
    int prec = p->prec;
    
    bool real = true;
    for (int i = 0; i <= n && real; i++) {
        real = mpfr_zero_p(p->a[i].y) != 0;
    }
    ldbl tang = real ? PI : 2 * PI;
    int div = real ? n / 2 + 1 : n;
    
    mpv roots = mpv_new(prec, 2 * n);
    int found = 0, rootsCount = 0;
    
    mpc root, sp;
    mpc_inits(prec, root, sp, NULL);
    
    int iter = 2 * n + 15 + 2 * prec;
    for (int step = 0; step < 3 && rootsCount < n; step ++) {
        ldbl sta = step == 0 ? 0 : tang / (2 * div);
        
        int steps = real && step == 0 ? div + 1 : div;
        for (int i = 0; i < steps && rootsCount < n; i++) {
            ldbl a = i * tang / div + sta;
            
            fp80 sp80 = {cosl(a), sinl(a)};
            fp80_mull(sp80, sp80, 2.5 * r);
            
            mpc_set80(sp, sp80);
            if(poly_root(p, root, sp, iter, eps, 1)) {
                mpfr_abs(root->y, root->y, MPFR_RNDN);
                bool exists = false;
                
                for (int j = 0; j < found && ! exists; j++) {
                    mpv_getc(sp, roots, j);
                    ldbl dst = mpc_distl(root, sp);
                    
                    exists = dst <= eps;
                }
                
                if(! exists) {
                    mpv_setc(roots, found++, root);
                    
                    rootsCount += mpfr_cmp_ld(root->y, eps / 2) <= 0 ? 1 : 2;
                }
            }
        }
        
        if(step > 0) {
            div *= div;
        }
    }
    
    mpc_clears(root, sp, NULL);
    
    roots->count = 2 * found;
    
    return roots;
}

int poly_sol(poly p, mpc r, mpc z, mpc t, int iter, long double err, int refine) {
    // watch out for buffer override
    mpc_set(p->x, z);
    
    int div = 0;
    int conv = 0;
    
    for (int i = 0; i < iter && ! div && ! conv; i++) {
        
        poly_newton_term_sol(p, p->n, p->x, t);
        mpc_sub(p->x, p->x, p->n);
        
        // divergence if z escapes too far
        mpc_mod(p->t->x, p->x);
        div = mpfr_cmp_ld(p->t->x, 1E5) >= 0; // FIXME: really ?
        
        // convergence if value small enough
        mpc_mod(p->t->y, p->s);
        conv = mpfr_cmp_ld(p->t->y, err) < 0;
    }
    
    if(div || ! conv)
        return 0;
    
    for (int i = 0; i < refine; i++) {
        
        poly_newton_term_sol(p, p->n, p->x, t);
        mpc_sub(p->x, p->x, p->n);

    }
    
    mpc_set(r, p->x);
    
    return 1;
}

