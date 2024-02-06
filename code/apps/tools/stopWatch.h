//
//  stopWatch.h
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
  \file stopWatch.h
  \brief A mini-collection of basic functions to time slow operations and provide human readable information to users.
 */

#ifndef STOPWATCH_H_
#define STOPWATCH_H_

#include <sys/timeb.h>
#include <time.h>
#include "ntypes.h"

/// One billion (nanoseconds in one second).
#define BILLION    1000000000L

/// One billion (miliseconds in one second).
#define MILLION       1000000L

/// Miliseconds in one minute.
#define MINUTE          60000L

/// Miliseconds in one hour.
#define HOUR              (60 * MINUTE)

/// Miliseconds in one day.
#define DAY               (24 * HOUR)

/// Short type name for real time.
typedef struct timeb rtime;

/// Short type name for precise time.
typedef struct timespec ptime;

/// @brief Formats the time duration in miliseconds @c ms into human readable string @c str.
///
/// @warning The string @c str should be at least @c 30 bytes long.
///
/// @param ms the time lapse in miliseconds
/// @param str the human readable print
void millis(long ms, char *str);

/// @brief Formats the time duration in nanoseconds @c ns into human readable string @c str.
///
/// @warning The string @c str should be at least @c 30 bytes long.
///
/// @param ns the time lapse in miliseconds
/// @param str the human readable print
void nanos(long ns, char *str);

/// @brief Computes how much time has passed since the time stamp @c ts.
///
/// Pretty prints the elapsed time into the string @c str and updates @c ts to the current time.
/// Does not fail if @c str==NULL, can also be used as a starter for the chronometer.
///
/// @warning The string @c str should be at least @c 30 bytes long.
///
/// @param ts old time stamp
/// @param str human readable string for the time increment
void lapse(rtime *ts, char *str);

/// @brief Pretty prints the current date & time into the string @c str.
///
/// The format for the time stamp is : "Year-Month-Day Hour:Minutes:Seconds".
///
/// The string should be at least @c 60 bytes long.
/// @param str human readable string that will be stamped
void time_stamp(char *str, int len, bool human, bool precise);

/// @brief Computes and prints the time lapse since @c ts.
///
/// Pretty prints the elapsed time into the string @c str and updates @c ts to the current time.
/// Does not fail if @c str==NULL, can also be used as a starter for the chronometer.
///
/// @warning The string @c str should be at least @c 30 bytes long.
///
/// @param clock the ID of the system clock to use
/// @param ts old time stamp
/// @param str human readable string for the time increment
///
/// @retrun the number of @c nano-seconds elapsed since the last call for @c ts
long lap(clockid_t clock, ptime *ts, char *str);

/// Return the number of seconds since [Epoch] (https://en.wikipedia.org/wiki/Epoch).
long secs(void);

/// Return the number of milli-seconds since [Epoch] (https://en.wikipedia.org/wiki/Epoch).
long msecs(void);

int date(char *now, int len);

#endif /* STOPWATCH_H_ */
