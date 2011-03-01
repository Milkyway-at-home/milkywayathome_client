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

#include <stdlib.h>
#include "nbody.h"
#include "nbody_priv.h"
#include "milkyway_util.h"
#include "nbody_step.h"
#include "grav.h"
#include "orbitintegrator.h"
#include "show.h"
#include "nbody_lua.h"

#if BOINC_APPLICATION
static inline void nbodyCheckpoint(const NBodyCtx* ctx, const NBodyState* st)
{
    if (boinc_time_to_checkpoint())
    {
        if (writeCheckpoint(ctx, st))
            fail("Failed to write checkpoint\n");

        boinc_checkpoint_completed();
    }

    boinc_fraction_done(st->tnow / ctx->timeEvolve);
}
#endif /* BOINC_APPLICATION */

static void runSystem(const NBodyCtx* ctx, NBodyState* st)
{
    const real tstop = ctx->timeEvolve - ctx->timestep / 1024.0;

    while (st->tnow < tstop)
    {
        stepSystem(ctx, st);   /* advance N-body system */
      #if BOINC_APPLICATION
        nbodyCheckpoint(ctx, st);
      #endif

      #if PERIODIC_OUTPUT
        st->outputTime = (st->outputTime + 1) % ctx->freqOut;
        if (st->outputTime == 0)
            outputBodyPositionBin(ctx, st);
      #endif /* PERIODIC_OUTPUT */
    }

  #if BOINC_APPLICATION
    mw_report("Making final checkpoint\n");
    if (writeCheckpoint(ctx, st))
        fail("Failed to write final checkpoint\n");
  #endif

}

static void endRun(NBodyCtx* ctx, NBodyState* st, const real chisq)
{
    finalOutput(ctx, st, chisq);
    nbodyCtxDestroy(ctx);     /* finish up output */
    nbodyStateDestroy(st);
}

static int setupRun(NBodyCtx* ctx, NBodyState* st)
{
    /* If the checkpoint exists, try to use it */
    if (resolvedCheckpointExists())
    {
        mw_report("Checkpoint exists. Attempting to resume from it.\n");
        if (readCheckpoint(ctx, st))
        {
            mw_report("Failed to read checkpoint\n");
            nbodyStateDestroy(st);
            return 1;
        }
        else
        {
            mw_report("Successfully read checkpoint\n");
        }
    }

    /* Accelerations are scratch space */
    st->acctab = (mwvector*) mwMallocA(ctx->nbody * sizeof(mwvector));

    gravMap(ctx, st); /* Start 1st step */

    return 0;
}

/* Set context fields read from command line flags */
static inline void nbodySetCtxFromFlags(NBodyCtx* ctx, const NBodyFlags* nbf)
{
    ctx->outputCartesian = nbf->outputCartesian;
    ctx->outputBodies    = nbf->printBodies;
    ctx->outputHistogram = nbf->printHistogram;
}

static int verifyFile(const NBodyCtx* ctx, int rc)
{
    printNBodyCtx(ctx);
    printHistogramParams(&ctx->histogramParams);

    if (rc)
        warn("File failed\n");
    else
        warn("File is OK\n");

    return rc;
}

/* Takes parsed json and run the simulation, using outFileName for
 * output. */
int runNBodySimulation(const NBodyFlags* nbf)       /* Misc. parameters to control output */
{
    NBodyCtx ctx  = EMPTY_NBODYCTX;
    NBodyState st = EMPTY_STATE;

    real chisq;
    double ts = 0.0, te = 0.0;
    int rc = 0;

    rc |= setupNBody(nbf->inputFile, &ctx, &st);
    if (rc && !nbf->verifyOnly)   /* Fail right away, unless we are diagnosing file problems */
        return warn1("Failed to read input parameters file\n");

    nbodySetCtxFromFlags(&ctx, nbf); /* These will get nuked by the checkpoint */

    if (nbf->verifyOnly)
        return verifyFile(&ctx, rc);
    if (rc)
        return warn1("Failed to read input parameters file\n");

    if (resolveCheckpoint(nbf->checkpointFileName))
        return warn1("Failed to resolve checkpoint\n");

    if (setupRun(&ctx, &st))
        return warn1("Failed to setup run\n");

    if (initOutput(&ctx, nbf))
        return warn1("Failed to open output files\n");

    if (nbf->printTiming)     /* Time the body of the calculation */
        ts = mwGetTime();

    runSystem(&ctx, &st);
    mw_report("Simulation complete\n");

    if (nbf->printTiming)
    {
        te = mwGetTime();
        printf("<run_time> %g </run_time>\n", te - ts);
    }

    /* Get the likelihood */
    chisq = nbodyChisq(&ctx, &st, nbf, &ctx.histogramParams);
    if (isnan(chisq))
        warn("Failed to calculate chisq\n");

    endRun(&ctx, &st, chisq);

    return 0;
}

