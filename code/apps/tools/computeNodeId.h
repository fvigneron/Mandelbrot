//
//  computeNodeId.h
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

/**
 \dir tools
 \brief This folder contains tools related to parallelization, debuging and profiling of the apps.
 */

 /**
  \file computeNodeId.h
  \brief A mini-collection of functions for parallelization of super-computers, as
  [ROMEO] (https://romeo.univ-reims.fr/pages/aboutUs) of Univ. de Reims Champagne-Ardenne.
 */

#ifndef computeNodeId_h
#define computeNodeId_h

#include <stdio.h>

/// @brief Gets CPU identifier.
///
/// On a UNIX style system, querries @c /proc/self/stat and collects entry 39 (valid since Linux 2.2.8).
/// See @c PROC(5) manual page.
///
/// Returns a unique CPU identifier, or -1 if the querry failed.
///
/// @return cpu indentifier (@c int)
int get_cpu_id(void);

#endif /* computeNodeId_h */
