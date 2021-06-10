/*
 *  Copyright (c) 2008-2010 Travis Desell, Nathan Cole, Dave Przybylo
 *  Copyright (c) 2008-2010 Boleslaw Szymanski, Heidi Newberg
 *  Copyright (c) 2008-2010 Carlos Varela, Malik Magdon-Ismail
 *  Copyright (c) 2008-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "milkyway_cpuid.h"
#include "milkyway_util.h"
#include "probabilities_dispatch.h"
#include "probabilities.h"


ProbabilityFunc probabilityFunc = NULL;


/* MSVC can't do weak imports. Using dlsym()/GetProcAddress() etc. would be better */
#if !HAVE_AVX || !DOUBLEPREC || defined(MSVC32_AVX_WORKAROUND)
  #define initProbabilities_AVX NULL
#endif

#if !HAVE_SSE41 || !DOUBLEPREC
#define initProbabilities_SSE41 NULL
#endif

#if !HAVE_SSE3 || !DOUBLEPREC
#define initProbabilities_SSE3 NULL
#endif

#if !HAVE_SSE2 || !DOUBLEPREC
#define initProbabilities_SSE2 NULL
#endif

/* Can't use the functions themselves if defined to NULL */
static ProbInitFunc initAVX = initProbabilities_AVX;
static ProbInitFunc initSSE41 = initProbabilities_SSE41;
static ProbInitFunc initSSE3 = initProbabilities_SSE3;
static ProbInitFunc initSSE2 = initProbabilities_SSE2;


static int usingIntrinsicsIsAcceptable(const AstronomyParameters* ap, int forceNoIntrinsics)
{
    if (!DOUBLEPREC)
    {
        return FALSE;
    }

    if (forceNoIntrinsics)
    {
        mw_printf("Forced to not use intrinsics functions\n");
        return FALSE;
    }

    if (ap->background_profile == SLOW_HERNQUIST || ap->aux_bg_profile)
    {
        mw_printf("Intrinsics are not acceptable for workunit\n");
        return FALSE;
    }

    return TRUE;
}

static ProbabilityFunc selectStandardFunction(const AstronomyParameters* ap)
{
	switch(ap->background_profile)
	{
	case SLOW_HERNQUIST:
		return probabilities_slow_hprob;
		break;
	case FAST_HERNQUIST:
		return probabilities_fast_hprob;
		break;
	case BROKEN_POWER_LAW:
		return probabilities_broken_power_law;
		break;
	default:
		return probabilities_slow_hprob;
		break;
	}
	return probabilities_slow_hprob;
}

#if MW_IS_X86 || defined(__e2k__)

/* Use one of the faster functions if available, or use something forced */
int probabilityFunctionDispatch(const AstronomyParameters* ap, const CLRequest* clr)
{
    int hasSSE2, hasSSE3, hasSSE41, hasAVX;
    int forcingInstructions = clr->forceAVX || clr->forceSSE41 || clr->forceSSE3 || clr->forceSSE2 || clr->forceX87;
    int abcd[4];

    if (!usingIntrinsicsIsAcceptable(ap, clr->forceNoIntrinsics))
    {
        probabilityFunc = selectStandardFunction(ap);
        return 0;
    }

    mw_cpuid(abcd, 1, 0);

#if defined(__e2k__)
   hasAVX = hasSSE41 = hasSSE3 = hasSSE2 = 1;
#else
    hasAVX = mwHasAVX(abcd) && mwOSHasAVXSupport();
    hasSSE41 = mwHasSSE41(abcd);
    hasSSE3 = mwHasSSE3(abcd);
    hasSSE2 = mwHasSSE2(abcd);
#endif

    if (clr->verbose)
    {
        mw_printf("CPU features:        SSE2 = %d, SSE3 = %d, SSE4.1 = %d, AVX = %d\n"
                  "Available functions: SSE2 = %d, SSE3 = %d, SSE4.1 = %d, AVX = %d\n"
                  "Forcing:             SSE2 = %d, SSE3 = %d, SSE4.1 = %d, AVX = %d\n",
                  hasSSE2, hasSSE3, hasSSE41, hasAVX,
                  initSSE2 != NULL, initSSE3 != NULL, initSSE41 != NULL, initAVX != NULL,
                  clr->forceSSE2, clr->forceSSE3, clr->forceSSE41, clr->forceAVX);
    }

    /* If multiple instructions are forced, the highest will take precedence */
    if (forcingInstructions)
    {
        if (clr->forceAVX && hasAVX && initAVX)
        {
            mw_printf("Using AVX path\n");
            probabilityFunc = initAVX(ap);
        }
        else if (clr->forceSSE41 && hasSSE41 && initSSE41)
        {
            mw_printf("Using SSE4.1 path\n");
            probabilityFunc = initSSE41(ap);
        }
        else if (clr->forceSSE3 && hasSSE3 && initSSE3)
        {
            mw_printf("Using SSE3 path\n");
            probabilityFunc = initSSE3(ap);
        }
        else if (clr->forceSSE2 && hasSSE2 && initSSE2)
        {
            mw_printf("Using SSE2 path\n");
            probabilityFunc = initSSE2(ap);
        }
        else if (clr->forceX87)
        {
            mw_printf("Using other path\n");
            probabilityFunc = selectStandardFunction(ap);
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
        if (hasAVX && initAVX)
        {
            mw_printf("Using AVX path\n");
            probabilityFunc = initAVX(ap);
        }
        else if (hasSSE41 && initSSE41)
        {
            mw_printf("Using SSE4.1 path\n");
            probabilityFunc = initSSE41(ap);
        }
        else if (hasSSE3 && initSSE3)
        {
            mw_printf("Using SSE3 path\n");
            probabilityFunc = initSSE3(ap);
        }
        else if (hasSSE2 && initSSE2)
        {
            mw_printf("Using SSE2 path\n");
            probabilityFunc = initSSE2(ap);
        }
        else
        {
            mw_printf("Using other path\n");
            probabilityFunc = selectStandardFunction(ap);
        }
    }

    if (!probabilityFunc)
    {
        mw_panic("Probability function not set!:\n"
                 "  Has AVX              = %d\n"
                 "  Has SSE4.1           = %d\n"
                 "  Has SSE3             = %d\n"
                 "  Has SSE2             = %d\n"
                 "  Forced AVX           = %d\n"
                 "  Forced SSE4.1        = %d\n"
                 "  Forced SSE3          = %d\n"
                 "  Forced SSE2          = %d\n"
                 "  Forced x87           = %d\n"
                 "  Forced no intrinsics = %d\n"
                 "  Arch                 = %s\n",
                 hasSSE41, hasSSE3, hasSSE2, hasAVX,
                 clr->forceAVX, clr->forceSSE41, clr->forceSSE3, clr->forceSSE2,
                 clr->forceX87, clr->forceNoIntrinsics,
                 ARCH_STRING);
    }

    return 0;
}

#else

int probabilityFunctionDispatch(const AstronomyParameters* ap, const CLRequest* clr)
{
    probabilityFunc = selectStandardFunction(ap);
    return 0;
}

#endif /* MW_IS_X86 */


