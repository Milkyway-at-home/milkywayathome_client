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

/* FIXME: This is just an arbitrary threshold I made up. What should it be? */
#define REQ(a, b) (rabs((a) - (b)) < 0.00001)

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

/* Taken from glibc */
#ifndef M_PI
# define M_E		2.7182818284590452354	/* e */
# define M_LOG2E	1.4426950408889634074	/* log_2 e */
# define M_LOG10E	0.43429448190325182765	/* log_10 e */
# define M_LN2		0.69314718055994530942	/* log_e 2 */
# define M_LN10		2.30258509299404568402	/* log_e 10 */
# define M_PI		3.14159265358979323846	/* pi */
# define M_PI_2		1.57079632679489661923	/* pi/2 */
# define M_PI_4		0.78539816339744830962	/* pi/4 */
# define M_1_PI		0.31830988618379067154	/* 1/pi */
# define M_2_PI		0.63661977236758134308	/* 2/pi */
# define M_2_SQRTPI	1.12837916709551257390	/* 2/sqrt(pi) */
# define M_SQRT2	1.41421356237309504880	/* sqrt(2) */
# define M_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */
#endif /* M_PI */

/* When we aren't using math.h (using fdlibm), the various constants
 * and NAN aren't defined, so just take this definition. */
#ifndef NAN
  #include "nan.h"
#endif /* NAN */

#ifndef isnan
  /* f != f only true for NaNs*/
  #define isnan(x) ( (x) != (x) )
#endif

__attribute__ ((visibility("default"))) double get_time();

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_UTIL_H_ */

