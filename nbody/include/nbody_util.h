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

#include "stdinc.h"
#include "real.h"
#include "vectmath.h"
#include "nbody_types.h"

/* Coordinate conversion */
void cartesianToLbr(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r);
void cartesianToLbr_rad(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r);

/* xrandom: generate floating-point random number. */
#define xrandom(xl, xh) ((real) (xl) + (real) ((xh) - (xl)) * drand48())

void* callocSafe(size_t count, size_t size);
void* mallocSafe(size_t size);

__attribute__ ((visibility("default"))) double get_time();

#endif /* _NBODY_UTIL_H_ */

