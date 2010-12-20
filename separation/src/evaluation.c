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

#if SEPARATION_OPENCL
  #include "run_cl.h"
  #include "setup_cl.h"
#endif /* SEPARATION_OPENCL */

#if SEPARATION_GRAPHICS
  #include "separation_graphics.h"
#endif /* SEPARATION_GRAPHICS */

#include <stdlib.h>
#include <stdio.h>


static void finalStreamIntegrals(FinalStreamIntegrals* fsi,
                                 const EvaluationState* es,
                                 const unsigned int number_streams,
                                 const unsigned int number_integrals)
{
    unsigned int i, j;

    fsi->stream_integrals = (real*) mwCalloc(number_streams, sizeof(real));

    fsi->background_integral = es->integrals[0].background_integral;
    for (i = 0; i < number_streams; ++i)
        fsi->stream_integrals[i] = es->integrals[0].stream_integrals[i];

    for (i = 1; i < number_integrals; ++i)
    {
        fsi->background_integral -= es->integrals[i].background_integral;
        for (j = 0; j < number_streams; j++)
            fsi->stream_integrals[j] -= es->integrals[i].stream_integrals[j];
    }

}

static void freeFinalStreamIntegrals(FinalStreamIntegrals* fsi)
{
    free(fsi->stream_integrals);
}

static void printStreamIntegrals(const FinalStreamIntegrals* fsi, const unsigned int number_streams)
{
    unsigned int i;
    fprintf(stderr, "<background_integral> %.20lf </background_integral>\n", fsi->background_integral);
    fprintf(stderr, "<stream_integrals>");
    for (i = 0; i < number_streams; i++)
        fprintf(stderr, " %.20lf", fsi->stream_integrals[i]);
    fprintf(stderr, " </stream_integrals>\n");
}

/* Add up completed integrals for progress reporting */
static inline unsigned int completedIntegralProgress(const IntegralArea* ias, const EvaluationState* es)
{
    const IntegralArea* ia;
    unsigned int i, current_calc_probs = 0;

    for (i = 0; i < es->current_integral; ++i)
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
    Integral* integral;
    const IntegralArea* ia;
    double t1, t2;

  #if SEPARATION_OPENCL
    DevInfo di;
    CLInfo ci = EMPTY_CL_INFO;
    cl_bool useImages = CL_TRUE;
  #endif

  #if SEPARATION_OPENCL
    if (setupSeparationCL(&ci, &di, ap, sc, sg, clr, useImages) != CL_SUCCESS)
        fail("Failed to setup up CL\n");

    useImages = useImages && di.imgSupport;
  #endif /* SEPARATION_OPENCL */

    for (; es->current_integral < ap->number_integrals; es->current_integral++)
    {
        integral = &es->integrals[es->current_integral];
        ia = &ias[es->current_integral];
        es->current_calc_probs = completedIntegralProgress(ias, es);

        t1 = mwGetTime();
      #if SEPARATION_OPENCL
        integral->background_integral = integrateCL(ap, ia, sc, sg,
                                                    integral->stream_integrals, es,
                                                    clr, &ci, &di, useImages);
      #elif SEPARATION_CAL
        integral->background_integral = integrateCAL(ap, ia, sc, sg,
                                                     integral->stream_integrals, es,
                                                     clr);
      #else
        integral->background_integral = integrate(ap, ia, sc, sg,
                                                  integral->stream_integrals, integral->probs, es);
      #endif /* SEPARATION_CL */

        t2 = mwGetTime();
        printf("Time = %.20g\n", t2 - t1);

        if (isnan(integral->background_integral))
            fail("Failed to calculate integral %u\n", es->current_integral);

        cleanStreamIntegrals(integral->stream_integrals, sc, ap->number_streams);

        CLEAR_KAHAN(es->sum);
    }


  #if SEPARATION_OPENCL
    mwDestroyCLInfo(&ci);
  #endif
}

real evaluate(const AstronomyParameters* ap,
              const IntegralArea* ias,
              const Streams* streams,
              const StreamConstants* sc,
              const char* star_points_file,
              const CLRequest* clr,
              const int do_separation,
              const char* separation_outfile)
{
    real likelihood_val = NAN;
    EvaluationState* es;
    StreamGauss sg;
    FinalStreamIntegrals fsi;
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

  #if BOINC_APPLICATION && !SEPARATION_OPENCL
    /* Final checkpoint. */
    if (writeCheckpoint(es))
        fail("Failed to write final checkpoint\n");
  #endif

    finalStreamIntegrals(&fsi, es, ap->number_streams, ap->number_integrals);
    printStreamIntegrals(&fsi, ap->number_streams);
    freeEvaluationState(es);

    if (readStarPoints(&sp, star_points_file))
        warn("Failed to read star points file\n");
    else
    {
        likelihood_val = likelihood(ap, &sp, sc, streams, &fsi, sg, do_separation, separation_outfile);
    }

    freeStarPoints(&sp);
    freeFinalStreamIntegrals(&fsi);
    freeStreamGauss(sg);

    return likelihood_val;
}


