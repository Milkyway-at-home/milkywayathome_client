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
#include "nbody_types.h"
#include "nbody_boinc.h"
#include "milkyway_util.h"


/* FIXME: This is just an arbitrary threshold I made up. What should it be? */
#define REQ(a, b) (rabs((a) - (b)) < 0.00001)

/* Coordinate conversion */
void cartesianToLbr(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r);
void cartesianToLbr_rad(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r);
void lbrToCartesian(const NBodyCtx* ctx, vectorptr cart, const vectorptr lbr);
void lbrToCartesian_rad(const NBodyCtx* ctx, vectorptr cart, const vectorptr lbr);

/* xrandom: generate floating-point random number */
#define xrandom(st, xl, xh) ((real) (xl) + (real) ((xh) - (xl)) * dsfmt_genrand_open_open((st)))


#ifdef __cplusplus
}
#endif

#endif /* _NBODY_UTIL_H_ */

