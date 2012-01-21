/*
 * Copyright (c) 2010-2011 Matthew Arsenault
 * Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody_plain.h"
#include "nbody_shmem.h"
#include "nbody_curses.h"
#include "nbody_defaults.h"
#include "nbody_util.h"
#include "nbody_checkpoint.h"
#include "nbody_grav.h"

/* If enough time has passed, record the next center of mass position */
static void nbAddTracePoint(const NBodyCtx* ctx, NBodyState* st)
{
    int i = st->step * N_ORBIT_TRACE_POINTS / ctx->nStep;

    if (st->usesExact) /* FIXME?. We don't get the CM without the tree */
        return;

    if (i >= N_ORBIT_TRACE_POINTS) /* Just in case */
        return;

    if (X(st->orbitTrace[i]) < REAL_MAX)
        return;

    st->orbitTrace[i] = Pos(st->tree.root);
}

static void nbReportProgress(const NBodyCtx* ctx, NBodyState* st)
{
    double frac = (double) st->step / (double) ctx->nStep;

    mw_fraction_done(frac);

    if (st->reportProgress)
    {
        mw_mvprintw(0, 0,
                    "Running: %f / %f (%f%%)\n",
                    frac * ctx->timeEvolve,
                    ctx->timeEvolve,
                    100.0 * frac
            );

        mw_refresh();
    }
}

static inline NBodyStatus nbCheckpoint(const NBodyCtx* ctx, NBodyState* st)
{
    if (nbTimeToCheckpoint(ctx, st))
    {
        if (nbWriteCheckpoint(ctx, st))
        {
            return NBODY_CHECKPOINT_ERROR;
        }

        mw_checkpoint_completed();
    }

    return NBODY_SUCCESS;
}

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

static inline void advancePosVel(NBodyState* st, const int nbody, const real dt)
{
    int i;
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

static inline void advanceVelocities(NBodyState* st, const int nbody, const real dt)
{
    int i;
    real dtHalf = 0.5 * dt;

  #ifdef _OPENMP
    #pragma omp parallel for private(i) schedule(static)
  #endif
    for (i = 0; i < nbody; ++i)      /* loop over all bodies */
        bodyAdvanceVel(&st->bodytab[i], st->acctab[i], dtHalf);
}

/* stepSystem: advance N-body system one time-step. */
NBodyStatus nbStepSystemPlain(const NBodyCtx* ctx, NBodyState* st)
{
    NBodyStatus rc;
    const real dt = ctx->timestep;

    advancePosVel(st, st->nbody, dt);

    rc = nbGravMap(ctx, st);
    advanceVelocities(st, st->nbody, dt);

    st->step++;

    return rc;
}

NBodyStatus nbRunSystemPlain(const NBodyCtx* ctx, NBodyState* st)
{
    NBodyStatus rc = NBODY_SUCCESS;

    rc |= nbGravMap(ctx, st); /* Calculate accelerations for 1st step this episode */
    if (nbStatusIsFatal(rc))
        return rc;

    while (st->step < ctx->nStep)
    {
        nbAddTracePoint(ctx, st);
        nbUpdateDisplayedBodies(ctx, st);
        rc |= nbStepSystemPlain(ctx, st);
        if (nbStatusIsFatal(rc))   /* advance N-body system */
            return rc;

        rc |= nbCheckpoint(ctx, st);
        if (nbStatusIsFatal(rc))
            return rc;

        nbReportProgress(ctx, st);
    }

    if (BOINC_APPLICATION || ctx->checkpointT >= 0)
    {
        mw_report("Making final checkpoint\n");
        if (nbWriteCheckpoint(ctx, st))
        {
            mw_printf("Failed to write final checkpoint\n");
            return NBODY_CHECKPOINT_ERROR;
        }
    }

    return rc;
}


