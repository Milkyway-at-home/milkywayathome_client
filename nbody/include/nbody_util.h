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

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

#include "stdinc.h"
#include "real.h"
#include "vectmath.h"
#include "nbody_types.h"

/* TODO: Naming, and sharing with separation */
#if BOINC_APPLICATION
  #include <boinc_api.h>
  #include <filesys.h>
  #define nbody_finish(x) boinc_finish(x)
  #define nbody_fopen(x,y) boinc_fopen((x),(y))
  #define nbody_remove(x) boinc_delete_file((x))
#else
  #define nbody_finish(x) exit(x)
  #define nbody_fopen(x,y) fopen((x),(y))
  #define nbody_remove(x) remove((x))
#endif /* BOINC_APPLICATION */

/* Coordinate conversion */
void cartesianToLbr(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r);
void cartesianToLbr_rad(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r);
void lbrToCartesian(const NBodyCtx* ctx, vectorptr cart, const vectorptr lbr);
void lbrToCartesian_rad(const NBodyCtx* ctx, vectorptr cart, const vectorptr lbr);

/* FIXME: New random source. drand48() doesn't work on Windows. */
/* xrandom: generate floating-point random number. */
#ifndef _WIN32
  #define xrandom(xl, xh) ((real) (xl) + (real) ((xh) - (xl)) * drand48())
  #define SET_SEED(x) (srand48((x)))
#else
  #define xrandom(xl, xh) ((real) (xl) + (real) ((xh) - (xl)) * ((double) rand() / (((double) RAND_MAX) + 1.0)))
  #define SET_SEED(x) ((void) 0)
#endif /* _WIN32 */

void* callocSafe(size_t count, size_t size);
void* mallocSafe(size_t size);

__attribute__ ((visibility("default"))) double get_time();

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_UTIL_H_ */

