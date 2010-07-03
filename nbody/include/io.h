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

#ifndef _IO_H_
#define _IO_H_

#include "nbody_types.h"
#include "vectmath.h"

/* Basic IO and initialization */
void initOutput(NBodyCtx*);             /* open files for output */
void output(const NBodyCtx* ctx, const NBodyState* st);
void boincOutput(const NBodyCtx* ctx, const NBodyState* st); /* Output with xml tags */

/* Checkpointing */
void openCheckpoint(NBodyCtx* ctx);
void closeCheckpoint(NBodyCtx* ctx);
void nbodyCheckpoint(const NBodyCtx* ctx, const NBodyState* st);
int thawState(const NBodyCtx* ctx, NBodyState* st);

/* Close output and free things */
void nbodyCtxDestroy(NBodyCtx* ctx);
void nbodyStateDestroy(NBodyState* st);

/* Coordinate conversion */
void cartesianToLbr(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r);
void cartesianToLbr_rad(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r);

#endif /* _IO_H_ */

