/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _NBODY_UTIL_H_
#define _NBODY_UTIL_H_

#include "nbody_boinc.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include "nbody_types.h"

/* Coordinate conversion */
void cartesianToLbr(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r);
void cartesianToLbr_rad(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r);
void lbrToCartesian(const NBodyCtx* ctx, vectorptr cart, const vectorptr lbr);
void lbrToCartesian_rad(const NBodyCtx* ctx, vectorptr cart, const vectorptr lbr);

/* xrandom: generate floating-point random number */
#define xrandom(st, xl, xh) ((real) (xl) + (real) ((xh) - (xl)) * dsfmt_genrand_open_open((st)))


void* callocSafe(size_t count, size_t size);
void* mallocSafe(size_t size);

#define warn(msg, ...) fprintf(stderr, msg, ##__VA_ARGS__)
#define fail(msg, ...) { fprintf(stderr, msg, ##__VA_ARGS__);  \
                         nbody_finish(EXIT_FAILURE); }


/*  ABS: returns the absolute value of its argument
 *  MAX: returns the argument with the highest value
 *  MIN: returns the argument with the lowest value
 */

#define   ABS(x)       (((x) < 0) ? -(x) : (x))
#define   MAX(x,y)     (((x) > (y)) ? (x) : (y))
#define   MIN(x,y)     (((x) < (y)) ? (x) : (y))


/* degrees to radians */
#define d2r(x) ((x) * M_PI / 180.0)

/* radians to degrees */
#define r2d(x) ((x) * 180.0 / M_PI)

/* simple math macros */
#define cube(x) ((x) * (x) * (x))
#define sqr(x)  ((x) * (x))
#define inv(x)  (1.0 / (x))

/* other useful nonstandard constants */

/* (4 * pi) / 3 Should be enough for 128-bit long double should we choose to use
 * that */
#define PI_4_3 (4.1887902047863909846168578443727)


__attribute__ ((visibility("default"))) double get_time();

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_UTIL_H_ */

