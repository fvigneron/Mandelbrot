//
//  ntypes.c
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
#include "ntypes.h"

bool ntypes_check(void) {
    if(sizeof(ldbl) <= 8 || sizeof(byte) != 1 || sizeof(word) != 2 ||
       sizeof(uint) != 4 || sizeof(ulong) != 8) {
        return false;
    }
    
    return (1.0L + 1E-18L) != 1.0L;
}

// takes care of quircks of intel / xeon + gcc that produce -nan
inline bool is_number(ldbl x) {
    long *b = (long *) &x;
    
    return (b[1] & 0x7FFF) != 0x7FFF;
}
