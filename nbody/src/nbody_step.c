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

#include "milkyway_util.h"
#include "nbody_step.h"
#include "nbody_grav.h"

/* Advance velocity by half a timestep */
static inline void bodyAdvanceVel(Body* p, const mwvector a, const real dtHalf)
{
    mwvector dv;

    dv = mw_mulvs(a, dtHalf);   /* get velocity increment */
    mw_incaddv(Vel(p), dv);     /* advance v by 1/2 step */
}

/* Advance body position by 1 timestep */
static inline void bodyAdvancePos(Body* p, const real dt)
{
    mwvector dr;

    dr = mw_mulvs(Vel(p), dt);  /* get position increment */
    mw_incaddv(Pos(p), dr);     /* advance r by 1 step */
}

static inline void advancePosVel(NBodyState* st, const unsigned int nbody, const real dt)
{
    unsigned int i;
    real dtHalf = 0.5 * dt;

  #ifdef _OPENMP
    #pragma omp parallel for private(i) schedule(static)
  #endif
    for (i = 0; i < nbody; ++i)
    {
        bodyAdvanceVel(&st->bodytab[i], st->acctab[i], dtHalf);
        bodyAdvancePos(&st->bodytab[i], dt);
    }
}

static inline void advanceVelocities(NBodyState* st, const unsigned int nbody, const real dt)
{
    unsigned int i;
    real dtHalf = 0.5 * dt;

  #ifdef _OPENMP
    #pragma omp parallel for private(i) schedule(static)
  #endif
    for (i = 0; i < nbody; ++i)      /* loop over all bodies */
        bodyAdvanceVel(&st->bodytab[i], st->acctab[i], dtHalf);
}

/* stepSystem: advance N-body system one time-step. */
NBodyStatus stepSystem(const NBodyCtx* ctx, NBodyState* st)
{
    NBodyStatus rc;
    const real dt = ctx->timestep;

    advancePosVel(st, st->nbody, dt);

    rc = gravMap(ctx, st);
    advanceVelocities(st, st->nbody, dt);

    st->tnow += dt;                           /* finally, advance time */

    return rc;
}

