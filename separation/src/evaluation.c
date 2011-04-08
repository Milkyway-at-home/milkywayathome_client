/*
Copyright 2008-2010 Travis Desell, Dave Przybylo, Nathan Cole, Matthew
Arsenault, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

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

#include "separation.h"
#include "separation_utils.h"

#if SEPARATION_OPENCL
  #include "run_cl.h"
  #include "setup_cl.h"
#elif SEPARATION_CAL
  #include "separation_cal_types.h"
  #include "separation_cal_setup.h"
  #include "separation_cal_run.h"
#endif /* SEPARATION_OPENCL */

#if SEPARATION_GRAPHICS
  #include "separation_graphics.h"
#endif /* SEPARATION_GRAPHICS */

#include <stdlib.h>
#include <stdio.h>


static void getFinalIntegrals(SeparationResults* results,
                              const EvaluationState* es,
                              const unsigned int number_streams,
                              const unsigned int number_integrals)
{
    unsigned int i, j;

    results->backgroundIntegral = es->cuts[0].bgIntegral;
    for (i = 0; i < number_streams; ++i)
        results->streamIntegrals[i] = es->cuts[0].streamIntegrals[i];

    for (i = 1; i < number_integrals; ++i)
    {
        results->backgroundIntegral -= es->cuts[i].bgIntegral;
        for (j = 0; j < number_streams; j++)
            results->streamIntegrals[j] -= es->cuts[i].streamIntegrals[j];
    }
}

#if 0
static void printStreamIntegrals(const FinalStreamIntegrals* fsi, const unsigned int number_streams)
{
    unsigned int i;
    fprintf(stderr, "<background_integral> %.20lf </background_integral>\n", fsi->background_integral);
    fprintf(stderr, "<stream_integrals>");
    for (i = 0; i < number_streams; i++)
        fprintf(stderr, " %.20lf", fsi->streamIntegrals[i]);
    fprintf(stderr, " </stream_integrals>\n");
}
#endif

/* Add up completed integrals for progress reporting */
static inline unsigned int completedIntegralProgress(const IntegralArea* ias, const EvaluationState* es)
{
    const IntegralArea* ia;
    unsigned int i, current_calc_probs = 0;

    for (i = 0; i < es->currentCut; ++i)
    {
        ia = &ias[i];
        current_calc_probs += ia->r_steps * ia->mu_steps * ia->nu_steps;
    }

    return current_calc_probs;
}

/* Zero insignificant streams */
static void cleanStreamIntegrals(real* stream_integrals,
                                 const StreamConstants* sc,
                                 const unsigned int number_streams)
{
    unsigned int i;

    for (i = 0; i < number_streams; ++i)
    {
        /* Rather than not adding up these streams, let them add and then
         * ignore them. They would have ended up being zero anyway */
        if (!sc[i].large_sigma)
            stream_integrals[i] = 0.0;
    }
}


static void calculateIntegrals(const AstronomyParameters* ap,
                               const IntegralArea* ias,
                               const StreamConstants* sc,
                               const StreamGauss sg,
                               EvaluationState* es,
                               const CLRequest* clr)
{
    const IntegralArea* ia;
    double t1, t2;
    int rc;

  #if SEPARATION_OPENCL
    DevInfo di;
    CLInfo ci = EMPTY_CL_INFO;
    cl_bool useImages = CL_TRUE;
  #elif SEPARATION_CAL
    MWCALInfo ci;
  #endif

  #if SEPARATION_OPENCL
    if (setupSeparationCL(&ci, &di, ap, clr, useImages) != CL_SUCCESS)
        fail("Failed to setup CL\n");

    useImages = useImages && di.imgSupport;

  #elif SEPARATION_CAL
    memset(&ci, 0, sizeof(MWCALInfo));
    if (separationCALInit(&ci, ap, sc, clr) != CAL_RESULT_OK)
        fail("Failed to setup CAL\n");

  #endif /* SEPARATION_OPENCL */

    for (; es->currentCut < es->numberCuts; es->currentCut++)
    {
        es->cut = &es->cuts[es->currentCut];
        ia = &ias[es->currentCut];
        es->current_calc_probs = completedIntegralProgress(ias, es);

        t1 = mwGetTime();
      #if SEPARATION_OPENCL
        rc = integrateCL(ap, ia, sc, sg, es, clr, &ci, &di, useImages);
      #elif SEPARATION_CAL
        rc = integrateCAL(ap, ia, sg, es, clr, &ci);
      #else
        rc = integrate(ap, ia, sc, sg, es);
      #endif /* SEPARATION_OPENCL */

        t2 = mwGetTime();
        warn("Integral %u time = %f s\n", es->currentCut, t2 - t1);

        if (rc || isnan(es->cut->bgIntegral))
            fail("Failed to calculate integral %u\n", es->currentCut);

        cleanStreamIntegrals(es->cut->streamIntegrals, sc, ap->number_streams);
        clearEvaluationStateTmpSums(es);
    }

  #if SEPARATION_OPENCL
    mwDestroyCLInfo(&ci);
  #elif SEPARATION_CAL
    mwCALShutdown(&ci);
  #endif
}

int evaluate(SeparationResults* results,
             const AstronomyParameters* ap,
             const IntegralArea* ias,
             const Streams* streams,
             const StreamConstants* sc,
             const char* star_points_file,
             const CLRequest* clr,
             const int do_separation,
             const char* separation_outfile)
{
    int rc = 0;
    EvaluationState* es;
    StreamGauss sg;
    StarPoints sp = EMPTY_STAR_POINTS;

    es = newEvaluationState(ap);
    sg = getStreamGauss(ap->convolve);

  #if SEPARATION_GRAPHICS
    if (separationInitSharedEvaluationState(es))
        warn("Failed to initialize shared evaluation state\n");
  #endif /* SEPARATION_GRAPHICS */

    if (resolveCheckpoint())
        fail("Failed to resolve checkpoint file '%s'\n", CHECKPOINT_FILE);

    if (maybeResume(es))
        fail("Failed to resume checkpoint\n");

    calculateIntegrals(ap, ias, sc, sg, es, clr);

  #if BOINC_APPLICATION && !SEPARATION_OPENCL && !SEPARATION_CAL
    /* Final checkpoint. */
    if (writeCheckpoint(es))
        fail("Failed to write final checkpoint\n");
  #endif

    getFinalIntegrals(results, es, ap->number_streams, ap->number_integrals);
    freeEvaluationState(es);

    if (readStarPoints(&sp, star_points_file))
    {
        rc = 1;
        warn("Failed to read star points file\n");
    }
    else
    {
        rc = likelihood(results, ap, &sp, sc, streams, sg, do_separation, separation_outfile);
        rc |= checkSeparationResults(results, ap->number_streams);
    }

    freeStarPoints(&sp);
    freeStreamGauss(sg);

    return rc;
}


