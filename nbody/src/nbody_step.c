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

/* Advance velocity by half a timestep */
ALWAYS_INLINE
static inline void bodyAdvanceVel(bodyptr p, const mwvector a, const real dt)
{
    mwvector dv;

    dv = mw_mulvs(a, 0.5 * dt);   /* get velocity increment */
    mw_incaddv(Vel(p), dv);       /* advance v by 1/2 step */
}


/* Advance body position by 1 timestep */
static inline void bodyAdvancePos(bodyptr p, const real dt)
{
    mwvector dr;

    dr = mw_mulvs(Vel(p), dt);  /* get position increment */
    mw_incaddv(Pos(p), dr);     /* advance r by 1 step */
}

#ifndef _OPENMP

/* Advance positions of all bodies by 1 timestep, and the velocities half a timestep. */
ALWAYS_INLINE
static inline void advancePosVel(NBodyState* st, const unsigned int nbody, const real dt)
{
    bodyptr p;
    mwvector* a;
    const bodyptr endp = st->bodytab + nbody;

    for (p = st->bodytab, a = st->acctab; p < endp; ++p, ++a)    /* loop over all bodies */
    {
        bodyAdvanceVel(p, *a, dt);
        bodyAdvancePos(p, dt);
    }
}

/* Advance velocities of all bodies by 1 timestep. */
ALWAYS_INLINE
static inline void advanceVelocities(NBodyState* st, const unsigned int nbody, const real dt)
{
    bodyptr p;
    mwvector* a;
    const bodyptr endp = st->bodytab + nbody;

    for (p = st->bodytab, a = st->acctab; p < endp; ++p, ++a)
        bodyAdvanceVel(p, *a, dt);
}

#else

ALWAYS_INLINE
static inline void advancePosVel(NBodyState* st, const unsigned int nbody, const real dt)
{
    unsigned int i;

  #pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < nbody; ++i)
    {
        bodyAdvanceVel(&st->bodytab[i], st->acctab[i], dt);
        bodyAdvancePos(&st->bodytab[i], dt);
    }

}

ALWAYS_INLINE
static inline void advanceVelocities(NBodyState* st, const unsigned int nbody, const real dt)
{
    unsigned int i;

    #pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < nbody; ++i)      /* loop over all bodies */
        bodyAdvanceVel(&st->bodytab[i], st->acctab[i], dt);
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

