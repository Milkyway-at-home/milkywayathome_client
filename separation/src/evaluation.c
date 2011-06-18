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


#if MW_IS_X86

#define bit_CMPXCHG8B (1 << 8)
#define bit_CMOV (1 << 15)
#define bit_MMX (1 << 23)
#define bit_SSE (1 << 25)
#define bit_SSE2 (1 << 26)
#define bit_SSE3 (1 << 0)
#define bit_CMPXCHG16B (1 << 13)
#define bit_3DNOW (1 << 31)
#define bit_3DNOWP (1 << 30)
#define bit_LM (1 << 29)


#ifndef _WIN32
  #if defined(__i386__) && defined(__PIC__)
    /* %ebx may be the PIC register.  */
    #define cpuid(level, a, b, c, d) __asm__ ("xchglt%%ebx, %1nt"                        \
                                              "cpuid"                                    \
                                              "xchglt%%ebx, %1nt"                        \
                                              : "=a" (a), "=r" (b), "=c" (c), "=d" (d)   \
                                              : "0" (level))
  #else
    #define cpuid(level, a, b, c, d) __asm__ ("cpuid"                                    \
                                              : "=a" (a), "=b" (b), "=c" (c), "=d" (d)   \
                                              : "0" (level))
  #endif /* defined(__i386__) && defined(__PIC__) */

static void getSSELevelSupport(int* hasSSE2, int* hasSSE3)
{
#ifndef __APPLE__
    /* http://peter.kuscsik.com/drupal/?q=node/10 */
    unsigned int eax, ebx, ecx, edx;
    cpuid(1, eax, ebx, ecx, edx);

    *hasSSE2 = !!(edx & bit_SSE2);
    *hasSSE3 = !!(ecx & bit_SSE3);
#else
    /* Something weird happens with this inline asm in Apple GCC that I
     * don't feel like figuring out now */
    *hasSSE2 = TRUE;
    *hasSSE3 = TRUE;
#endif /* __APPLE__ */
}


#else

static void getSSELevelSupport(int* hasSSE2_out, int* hasSSE3_out)
{
    int nIds = 0;
    int cpuInfo[4] = { 0, 0, 0, 0 };
    int hasSSE2 = FALSE;
    int hasSSE3 = FALSE;

    __cpuid(cpuInfo, 0);
    nIds = cpuInfo[0];

    if (nIds >= 1)
    {
        __cpuid(cpuInfo, 1);
        hasSSE2 = !!(cpuInfo[3] & bit_SSE2);
        hasSSE3 = !!(cpuInfo[2] & bit_SSE3);
    }

    *hasSSE2_out = hasSSE2;
    *hasSSE3_out = hasSSE3;
}
#endif /* _WIN32 */


/* Use one of the faster functions if available */
static void probabilityFunctionDispatch(const AstronomyParameters* ap, const CLRequest* clr)
{
    int hasSSE2 = FALSE, hasSSE3 = FALSE;
    int useSSE2 = FALSE, useSSE3 = FALSE;

    getSSELevelSupport(&hasSSE2, &hasSSE3);

    if (clr->verbose)
    {
        warn("CPU features: SSE2 = %d, SSE3 = %d\n"
             "Forcing: SSE2 = %d, SSE3 = %d, x87 = %d\n",
             hasSSE2, hasSSE3,
             clr->forceSSE2, clr->forceSSE3, clr->forceX87);
    }

    if (!clr->forceSSE2 && !clr->forceSSE3 && !clr->forceX87)
    {
        /* Default to using highest capability if not forcing anything */
        useSSE3 = hasSSE3;
        useSSE2 = hasSSE2;
    }
    else if (clr->forceSSE2)
    {
        useSSE2 = TRUE;
    }
    else if (clr->forceSSE3)
    {
        useSSE3 = TRUE;
    }
    else if (clr->forceX87)
    {
      #if MW_IS_X86_32
        useSSE2 = useSSE3 = FALSE;
      #elif MW_IS_X86_64
        useSSE2 = TRUE;  /* Ignore flag */
      #endif
    }

    if (useSSE3)  /* Precedence to higher level */
    {
        warn("Using SSE3 path\n");
        probabilityFunc = initProbabilities_SSE3(ap, clr->forceNoIntrinsics);
    }
    else if (useSSE2)
    {
        warn("Using SSE2 path\n");
        probabilityFunc = initProbabilities_SSE2(ap, clr->forceNoIntrinsics);
    }
    else
    {
      #if !MW_NO_X87_EVER
        warn("Using x87 path\n");
        probabilityFunc = initProbabilities(ap, clr->forceNoIntrinsics);
      #else
        mw_unreachable();
      #endif
    }

    if (!probabilityFunc)
    {
        mw_panic("Probability function not set!:\n"
                 "  Has SSE2             = %d\n"
                 "  Has SSE3             = %d\n"
                 "  Forced SSE3          = %d\n"
                 "  Forced SSE2          = %d\n"
                 "  Forced x87           = %d\n"
                 "  Forced no intrinsics = %d\n"
                 "  Arch                 = %s\n",
                 hasSSE2, hasSSE3,
                 clr->forceSSE3, clr->forceSSE2, clr->forceX87, clr->forceNoIntrinsics,
                 ARCH_STRING);
    }
}


#else
static void probabilityFunctionDispatch(const AstronomyParameters* ap, const CLRequest* clr)
{
    probabilityFunc = initProbabilities(ap, clr->forceNoIntrinsics);
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
        fail("Failed to write final checkpoint\n");

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
        fail("Failed to load integral kernel");
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
        warn("Integral %u time = %f s\n", es->currentCut, t2 - t1);

        if (rc || isnan(es->cut->bgIntegral))
            fail("Failed to calculate integral %u\n", es->currentCut);

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

    probabilityFunctionDispatch(ap, clr);

    es = newEvaluationState(ap);
    sg = getStreamGauss(ap->convolve);

  #if SEPARATION_GRAPHICS
    if (separationInitSharedEvaluationState(es))
        warn("Failed to initialize shared evaluation state\n");
  #endif /* SEPARATION_GRAPHICS */

    if (!ignoreCheckpoint)
    {
        if (resolveCheckpoint())
            fail("Failed to resolve checkpoint file '%s'\n", CHECKPOINT_FILE);

        if (maybeResume(es))
            fail("Failed to resume checkpoint\n");
    }

  #if SEPARATION_OPENCL
    if (setupSeparationCL(&ci, ap, ias, clr, &useImages) != CL_SUCCESS)
        fail("Failed to setup CL\n");
  #elif SEPARATION_CAL
    if (separationCALInit(&ci, clr) != CAL_RESULT_OK)
        fail("Failed to setup CAL\n");
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
        warn("Failed to read star points file\n");
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

