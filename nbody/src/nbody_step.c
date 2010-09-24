/*
  Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
  Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
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

#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */

#include "nbody_step.h"
#include "grav.h"


static inline void bodyAdvancePosVel(bodyptr p, const vectorptr a, const real dt)
{
    vector dvel;
    vector dpos;

    MULVS(dvel, (vectorptr) a, 0.5 * dt);      /* get velocity increment */
    INCADDV(Vel(p), dvel);                     /* advance v by 1/2 step */
    MULVS(dpos, Vel(p), dt);                   /* get positon increment */
    INCADDV(Pos(p), dpos);                     /* advance r by 1 step */
}

static inline void bodyAdvanceVelocity(bodyptr p, const vectorptr a, const real dt)
{
    vector dvel;

    MULVS(dvel, (vectorptr) a, 0.5 * dt);   /* get velocity increment */
    INCADDV(Vel(p), dvel);                  /* advance v by 1/2 step */
}

#ifndef _OPENMP

static inline void advancePosVel(NBodyState* st, const unsigned int nbody, const real dt)
{
    bodyptr p;
    vector* a;
    const bodyptr endp = st->bodytab + nbody;

    for (p = st->bodytab, a = st->acctab; p < endp; ++p, ++a)    /* loop over all bodies */
        bodyAdvancePosVel(p, (vectorptr) a, dt);
}


static inline void advanceVelocities(NBodyState* st, const unsigned int nbody, const real dt)
{
    bodyptr p;
    vector* a;
    const bodyptr endp = st->bodytab + nbody;

    for (p = st->bodytab, a = st->acctab; p < endp; ++p, ++a)
        bodyAdvanceVelocity(p, (vectorptr) a, dt);
}

#else

static inline void advancePosVel(NBodyState* st, const unsigned int nbody, const real dt)
{
    unsigned int i;

  #pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < nbody; ++i)
        bodyAdvancePosVel(&st->bodytab[i], &st->acctab[i], dt);

}

static inline void advanceVelocities(NBodyState* st, const unsigned int nbody, const real dt)
{
    unsigned int i;

    #pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < nbody; ++i)      /* loop over all bodies */
        bodyAdvanceVelocity(&st->bodytab[i], (vectorptr) &st->acctab[i], dt);
}

#endif /* _OPENMP */

/* stepSystem: advance N-body system one time-step. */
void stepSystem(const NBodyCtx* ctx, NBodyState* st)
{
    const real dt = ctx->model.timestep;

    advancePosVel(st, ctx->model.nbody, dt);

  #if !NBODY_OPENCL
    gravMap(ctx, st);
  #else
    gravMapCL(ctx, st);
  #endif /* !NBODY_OPENCL */

    advanceVelocities(st, ctx->model.nbody, dt);

    st->tnow += dt;                           /* finally, advance time */
}

