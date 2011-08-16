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
#include "probabilities.h"
#include "mw_cpuid.h"

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

/* FIXME: deal with this correctly */
#if SEPARATION_CAL
typedef MWCALInfo GPUInfo;
#elif SEPARATION_OPENCL
typedef CLInfo GPUInfo;
#else
typedef int GPUInfo;
#endif


ProbabilityFunc probabilityFunc = NULL;

/* MSVC can't do weak imports */
#if !HAVE_SSE41
#define initProbabilities_SSE41 NULL
#endif

#if !HAVE_SSE3
#define initProbabilities_SSE3 NULL
#endif

#if !HAVE_SSE2
#define initProbabilities_SSE2 NULL
#endif

#if !HAVE_OTHER_MATH_X87
  #define initProbabilities NULL
#endif

/* Can't use the functions themselves if defined to NULL */
static ProbInitFunc initSSE41 = initProbabilities_SSE41;
static ProbInitFunc initSSE3 = initProbabilities_SSE3;
static ProbInitFunc initSSE2 = initProbabilities_SSE2;
static ProbInitFunc initOther = initProbabilities;

#if MW_IS_X86

/* Use one of the faster functions if available, or use something forced */
static int probabilityFunctionDispatch(const AstronomyParameters* ap, const CLRequest* clr)
{
    int hasSSE2, hasSSE3, hasSSE41;
    int forcingInstructions = clr->forceSSE41 || clr->forceSSE3 || clr->forceSSE2 || clr->forceX87;
    int abcd[4];

    mw_cpuid(abcd, 1, 0);

    hasSSE41 = mwHasSSE41(abcd);
    hasSSE3 = mwHasSSE3(abcd);
    hasSSE2 = mwHasSSE2(abcd);

    if (clr->verbose)
    {
        mw_printf("CPU features:        SSE2 = %d, SSE3 = %d, SSE4.1 = %d\n"
                  "Available functions: SSE2 = %d, SSE3 = %d, SSE4.1 = %d, x87 = %d\n"
                  "Forcing:             SSE2 = %d, SSE3 = %d, SSE4.1 = %d, x87 = %d\n",
                  hasSSE2, hasSSE3, hasSSE41,
                  initSSE2 != NULL, initSSE3 != NULL, initSSE41 != NULL, initProbabilities != NULL,
                  clr->forceSSE2, clr->forceSSE3, clr->forceSSE41, clr->forceX87);
    }

    /* If multiple instructions are forced, the highest will take precedence */
    if (forcingInstructions)
    {
        if (clr->forceSSE41 && hasSSE41 && initSSE41)
        {
            mw_printf("Using SSE4.1 path\n");
            probabilityFunc = initSSE41(ap, clr->forceNoIntrinsics);
        }
        else if (clr->forceSSE3 && hasSSE3 && initSSE3)
        {
            mw_printf("Using SSE3 path\n");
            probabilityFunc = initSSE3(ap, clr->forceNoIntrinsics);
        }
        else if (clr->forceSSE2 && hasSSE2 && initSSE2)
        {
            mw_printf("Using SSE2 path\n");
            probabilityFunc = initSSE2(ap, clr->forceNoIntrinsics);
        }
        else if (clr->forceX87 && initOther)
        {
            mw_printf("Using other path\n");
            probabilityFunc = initOther(ap, clr->forceNoIntrinsics);
        }
        else
        {
            mw_printf("Tried to force an unusable path\n");
            return 1;
        }
    }
    else
    {
        /* Choose the highest level with available function and instructions */
        if (hasSSE41 && initSSE41)
        {
            mw_printf("Using SSE4.1 path\n");
            probabilityFunc = initSSE41(ap, clr->forceNoIntrinsics);
        }
        else if (hasSSE3 && initSSE3)
        {
            mw_printf("Using SSE3 path\n");
            probabilityFunc = initSSE3(ap, clr->forceNoIntrinsics);
        }
        else if (hasSSE2 && initSSE2)
        {
            mw_printf("Using SSE2 path\n");
            probabilityFunc = initSSE2(ap, clr->forceNoIntrinsics);
        }
        else if (initOther)
        {
            mw_printf("Using other path\n");
            probabilityFunc = initOther(ap, clr->forceNoIntrinsics);
        }
        else
        {
            mw_printf("No paths usable\n");
            return 1;
        }
    }

    if (!probabilityFunc)
    {
        mw_panic("Probability function not set!:\n"
                 "  Has SSE4.1           = %d\n"
                 "  Has SSE3             = %d\n"
                 "  Has SSE2             = %d\n"
                 "  Forced SSE4.1        = %d\n"
                 "  Forced SSE3          = %d\n"
                 "  Forced SSE2          = %d\n"
                 "  Forced x87           = %d\n"
                 "  Forced no intrinsics = %d\n"
                 "  Arch                 = %s\n",
                 hasSSE41, hasSSE3, hasSSE2,
                 clr->forceSSE41, clr->forceSSE3, clr->forceSSE2,
                 clr->forceX87, clr->forceNoIntrinsics,
                 ARCH_STRING);
    }

    return 0;
}


#else
static int probabilityFunctionDispatch(const AstronomyParameters* ap, const CLRequest* clr)
{
    probabilityFunc = initProbabilities(ap, clr->forceNoIntrinsics);
    return 0;
}
#endif /* MW_IS_X86 */


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

static void finalCheckpoint(EvaluationState* es)
{
  #if BOINC_APPLICATION
    boinc_begin_critical_section();
  #endif

    if (writeCheckpoint(es))
        mw_fail("Failed to write final checkpoint\n");

  #if BOINC_APPLICATION
    boinc_end_critical_section();
  #endif
}

static void calculateIntegrals(const AstronomyParameters* ap,
                               const IntegralArea* ias,
                               const StreamConstants* sc,
                               const StreamGauss sg,
                               EvaluationState* es,
                               const CLRequest* clr,
                               GPUInfo* ci,
                               int useImages)
{
    const IntegralArea* ia;
    double t1, t2;
    int rc;

  #if SEPARATION_CAL
    if (separationLoadKernel(ci, ap, sc) != CAL_RESULT_OK)
        mw_fail("Failed to load integral kernel");
  #endif /* SEPARATION_OPENCL */

    for (; es->currentCut < es->numberCuts; es->currentCut++)
    {
        es->cut = &es->cuts[es->currentCut];
        ia = &ias[es->currentCut];
        es->current_calc_probs = completedIntegralProgress(ias, es);

        t1 = mwGetTime();
      #if SEPARATION_OPENCL
        rc = integrateCL(ap, ia, sc, sg, es, clr, ci, (cl_bool) useImages);
      #elif SEPARATION_CAL
        rc = integrateCAL(ap, ia, sg, es, clr, ci);
      #else
        rc = integrate(ap, ia, sc, sg, es, clr);
      #endif /* SEPARATION_OPENCL */

        t2 = mwGetTime();
        mw_printf("Integral %u time = %f s\n", es->currentCut, t2 - t1);

        if (rc || isnan(es->cut->bgIntegral))
            mw_fail("Failed to calculate integral %u\n", es->currentCut);

        cleanStreamIntegrals(es->cut->streamIntegrals, sc, ap->number_streams);
        clearEvaluationStateTmpSums(es);
    }

  #if SEPARATION_CAL
    mwUnloadKernel(ci);
  #endif /* SEPARATION_CAL */
}

int evaluate(SeparationResults* results,
             const AstronomyParameters* ap,
             const IntegralArea* ias,
             const Streams* streams,
             const StreamConstants* sc,
             const char* star_points_file,
             const CLRequest* clr,
             int do_separation,
             int ignoreCheckpoint,
             const char* separation_outfile)
{
    int rc = 0;
    EvaluationState* es;
    StreamGauss sg;
    GPUInfo ci;
    StarPoints sp = EMPTY_STAR_POINTS;
    int useImages = FALSE; /* Only applies to CL version */

    memset(&ci, 0, sizeof(ci));

    if (probabilityFunctionDispatch(ap, clr))
        return 1;

    es = newEvaluationState(ap);
    sg = getStreamGauss(ap->convolve);

  #if SEPARATION_GRAPHICS
    if (separationInitSharedEvaluationState(es))
        mw_printf("Failed to initialize shared evaluation state\n");
  #endif /* SEPARATION_GRAPHICS */

    if (!ignoreCheckpoint)
    {
        if (resolveCheckpoint())
            mw_fail("Failed to resolve checkpoint file '%s'\n", CHECKPOINT_FILE);

        if (maybeResume(es))
            mw_fail("Failed to resume checkpoint\n");
    }

  #if SEPARATION_OPENCL
    if (setupSeparationCL(&ci, ap, ias, clr, &useImages) != CL_SUCCESS)
        mw_fail("Failed to setup CL\n");
  #elif SEPARATION_CAL
    if (separationCALInit(&ci, clr) != CAL_RESULT_OK)
        mw_fail("Failed to setup CAL\n");
  #endif

    calculateIntegrals(ap, ias, sc, sg, es, clr, &ci, useImages);

    if (!ignoreCheckpoint)
    {
        finalCheckpoint(es);
    }

    getFinalIntegrals(results, es, ap->number_streams, ap->number_integrals);
    freeEvaluationState(es);

    if (readStarPoints(&sp, star_points_file))
    {
        rc = 1;
        mw_printf("Failed to read star points file\n");
    }
    else
    {
        /* TODO: likelihood on GPU with OpenCL. Make this less of a
         * mess. The different versions should appear to be the
         * same. */

      #if SEPARATION_CAL
        if (do_separation)
        {
            /* No separation on GPU */
            rc = likelihood(results, ap, &sp, sc, streams, sg, do_separation, separation_outfile);
        }
        else
        {
            //rc = likelihoodCAL(results, ap, &sp, sc, streams, sg, clr, &ci);
            rc = likelihood(results, ap, &sp, sc, streams, sg, do_separation, separation_outfile);
        }
      #else
        rc = likelihood(results, ap, &sp, sc, streams, sg, do_separation, separation_outfile);
      #endif /* SEPARATION_CAL */

        rc |= checkSeparationResults(results, ap->number_streams);
    }

    freeStarPoints(&sp);
    freeStreamGauss(sg);

  #if SEPARATION_OPENCL
    mwDestroyCLInfo(&ci);
  #elif SEPARATION_CAL
    mwCALShutdown(&ci);
  #endif

    return rc;
}

