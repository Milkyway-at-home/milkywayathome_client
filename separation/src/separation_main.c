/*
 *  Copyright (c) 2008-2010 Travis Desell, Nathan Cole, Dave Przybylo
 *  Copyright (c) 2008-2010 Boleslaw Szymanski, Heidi Newberg
 *  Copyright (c) 2008-2010 Carlos Varela, Malik Magdon-Ismail
 *  Copyright (c) 2008-2012 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2012 Matthew Arsenault
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

#ifdef _WIN32
  #if MW_IS_X86_32
    #pragma comment(linker, "/SUBSYSTEM:CONSOLE,5.0")
  #else
    #pragma comment(linker, "/SUBSYSTEM:CONSOLE,5.02")
  #endif
#endif

#include "separation.h"
#include "separation_lua.h"
#include "milkyway_util.h"
#include "milkyway_boinc_util.h"
#include "milkyway_git_version.h"
#include "io_util.h"
#include <popt.h>


#define DEFAULT_ASTRONOMY_PARAMETERS "astronomy_parameters.txt"
#define DEFAULT_STAR_POINTS "stars.txt"

#define SEED_ARGUMENT (1 << 1)

const char* separationCommitID = MILKYWAY_GIT_COMMIT_ID;
const char* separationCommitDescribe = MILKYWAY_GIT_DESCRIBE;

static void printCopyright()
{
    mw_printf(
        "Milkyway@Home Separation client %d.%d\n\n"
        "Copyright (c) 2008-2011 Travis Desell, Nathan Cole, Boleslaw Szymanski\n"
        "Copyright (c) 2008-2011 Heidi Newberg, Carlos Varela, Malik Magdon-Ismail\n"
        "Copyright (c) 2008-2012 Rensselaer Polytechnic Institute.\n"
        "Copyright (c) 2010-2012 Matthew Arsenault\n"
        "Copyright (c) 1991-2000 University of Groningen, The Netherlands.\n"
        "Copyright (c) 2001-2009 The GROMACS Development Team\n"
        "\n"
        "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        "This is free software: you are free to change and redistribute it.\n"
        "There is NO WARRANTY, to the extent permitted by law.\n"
        "\n"
        " Incorporates works covered by the following copyright and\n"
        " permission notices:\n"
        "\n"
        "Copyright (C) 2007, 2008 Mutsuo Saito, Makoto Matsumoto and Hiroshima University\n"
        "Copyright (c) 2010, Naoaki Okazaki\n"
        "\n"
        " Redistribution and use in source and binary forms, with or without\n"
        " modification, are permitted provided that the following conditions are met:\n"
        "     * Redistributions of source code must retain the above copyright\n"
        "       notice, this list of conditions and the following disclaimer.\n"
        "     * Redistributions in binary form must reproduce the above copyright\n"
        "       notice, this list of conditions and the following disclaimer in the\n"
        "       documentation and/or other materials provided with the distribution.\n"
        "     * Neither the names of the authors nor the names of its contributors\n"
        "       may be used to endorse or promote products derived from this\n"
        "       software without specific prior written permission.\n"
        "\n"
        " THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
        " \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
        " LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n"
        " A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER\n"
        " OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,\n"
        " EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,\n"
        " PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR\n"
        " PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF\n"
        " LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING\n"
        " NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\n"
        " SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n"
        "\n"
        "\n"
        "Copyright (C) 1994-2008 Lua.org, PUC-Rio.\n"
        "\n"
        " Permission is hereby granted, free of charge, to any person obtaining a copy\n"
        " of this software and associated documentation files (the \"Software\"), to deal\n"
        " in the Software without restriction, including without limitation the rights\n"
        " to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
        " copies of the Software, and to permit persons to whom the Software is\n"
        " furnished to do so, subject to the following conditions:\n"
        "\n"
        " The above copyright notice and this permission notice shall be included in\n"
        " all copies or substantial portions of the Software.\n"
        "\n"
        " THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
        " IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
        " FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE\n"
        " AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
        " LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
        " OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
        " THE SOFTWARE.\n"
        "\n",
        SEPARATION_VERSION_MAJOR,
        SEPARATION_VERSION_MINOR
        );
}


static void printVersion(int boincTag, int verbose)
{
    char versionStr[2048];

    snprintf(versionStr, sizeof(versionStr), "%s %u.%02u %s %s %s%s%s%s",
             SEPARATION_PROJECT_NAME,
             SEPARATION_VERSION_MAJOR, SEPARATION_VERSION_MINOR,
             MILKYWAY_SYSTEM_NAME,
             ARCH_STRING,
             PRECSTRING,
             DENORMAL_STRING,
             SEPARATION_SPECIAL_STR,
             SEPARATION_SPECIAL_LIBM_STR);

    if (boincTag)
    {
        mw_printf("<search_application> %s </search_application>\n", versionStr);
    }
    else
    {
        mw_printf("%s %s\n",
                  versionStr,
                  BOINC_APPLICATION ? "BOINC" : "");
    }

    if (verbose)
    {
        mw_printf("Commit %s\n", MILKYWAY_GIT_COMMIT_ID);
    }
}

static void freeSeparationFlags(SeparationFlags* sf)
{
    free(sf->star_points_file);
    free(sf->ap_file);
    free(sf->separation_outfile);
    free(sf->forwardedArgs);
    free(sf->numArgs);
    free(sf->preferredPlatformVendor);
}

/* Use hardcoded names if files not specified for compatability */
static void setDefaults(SeparationFlags* sf)
{
    mwStringDefault(sf->star_points_file, DEFAULT_STAR_POINTS);
    mwStringDefault(sf->ap_file, DEFAULT_ASTRONOMY_PARAMETERS);
}

static void setCLReqFlags(CLRequest* clr, const SeparationFlags* sf)
{
    clr->forceNoIntrinsics = sf->forceNoIntrinsics;
    clr->forceX87 = sf->forceX87;
    clr->forceSSE2 = sf->forceSSE2;
    clr->forceSSE3 = sf->forceSSE3;
    clr->forceSSE41 = sf->forceSSE41;
    clr->forceAVX = sf->forceAVX;
    clr->verbose = sf->verbose;
    clr->nonResponsive = sf->nonResponsive;
    clr->enableCheckpointing = !sf->disableGPUCheckpointing;

    clr->devNum = sf->useDevNumber;
    clr->platform = sf->usePlatform;
    clr->preferredPlatformVendor = sf->preferredPlatformVendor;

    clr->targetFrequency = (sf->targetFrequency <= 0.01) ? DEFAULT_TARGET_FREQUENCY : sf->targetFrequency;
    clr->gpuWaitFactor = (sf->waitFactor <= 0.01 || sf->waitFactor > 10.0) ? DEFAULT_WAIT_FACTOR : sf->waitFactor;
    clr->pollingMode = sf->pollingMode;

    clr->forceNoILKernel = sf->forceNoILKernel;
    clr->forceNoOpenCL = sf->forceNoOpenCL;

    clr->enableProfiling = FALSE;
}

typedef struct
{
    double gpuTargetFrequency;
    double gpuWaitFactor;
    int gpuNonResponsive;
    int gpuProcessPriority;
    int gpuDisableCheckpoint;
} SeparationPrefs;


/* If using BOINC try reading a few of the settings from the project
 * preferences. If command line arguments are used, those will
 * override the preferences. The command line arguments will also
 * still work without BOINC */
static void separationReadPreferences(SeparationPrefs* prefsOut)
{
    static SeparationPrefs prefs;
    static MWProjectPrefs sepPrefs[] =
        {
            { "gpu_target_frequency", MW_PREF_DOUBLE, FALSE, &prefs.gpuTargetFrequency   },
            { "gpu_wait_factor",      MW_PREF_DOUBLE, FALSE, &prefs.gpuWaitFactor        },
            { "gpu_non_responsive",   MW_PREF_BOOL,   FALSE, &prefs.gpuNonResponsive     },
            { "gpu_process_priority", MW_PREF_INT,    FALSE, &prefs.gpuProcessPriority   },
            { "no_gpu_checkpoint",    MW_PREF_BOOL,   FALSE, &prefs.gpuDisableCheckpoint },
            END_MW_PROJECT_PREFS
        };

    prefs.gpuTargetFrequency   = DEFAULT_TARGET_FREQUENCY;
    prefs.gpuWaitFactor        = DEFAULT_WAIT_FACTOR;
    prefs.gpuNonResponsive     = DEFAULT_NON_RESPONSIVE;
    prefs.gpuProcessPriority   = DEFAULT_GPU_PRIORITY;
    prefs.gpuDisableCheckpoint = DEFAULT_DISABLE_GPU_CHECKPOINTING;

    if (BOINC_APPLICATION)
    {
        if (mwGetAppInitData())
        {
            mw_printf("Error reading app init data. Project preferences will not be used\n");
        }
        else
        {
            mwReadProjectPrefs(sepPrefs, mwGetProjectPrefs());
        }
    }

    *prefsOut = prefs;
}

static void setInitialFlags(SeparationFlags* sf)
{
    memset(sf, 0, sizeof(*sf));

    /* Preferences or BOINC settings */
    sf->useDevNumber = -1;
    sf->targetFrequency = -1.0;
    sf->processPriority = MW_PRIORITY_INVALID;
    sf->nonResponsive = -1;

    sf->usePlatform = -1;

    sf->nonResponsive = DEFAULT_NON_RESPONSIVE;
    sf->pollingMode = DEFAULT_POLLING_MODE;
    sf->disableGPUCheckpointing = DEFAULT_DISABLE_GPU_CHECKPOINTING;
    sf->forceNoOpenCL = DEFAULT_DISABLE_OPENCL;
    sf->forceNoILKernel = DEFAULT_DISABLE_IL_KERNEL;
    sf->background = 0;
}

/* Set any flags based on project preferences that weren't specified
 * on the command line.
 *
 * This is a bit convoluted since we need to boinc_init before we read
 * the preferences, but that needs to be delayed until after argument
 * reading in case we want to disable output redirection, and then we
 * still want the command line to supersede the project prefs / the
 * device specified by app_init_data
 */
static void setFlagsFromPreferences(SeparationFlags* flags, const SeparationPrefs* prefs, const char* progName)
{
    if (flags->useDevNumber < 0)
    {
        /* Try to use BOINC's suggestion from app_init_data stuff;
           We might not get it so just use the first device. */
        flags->useDevNumber = mwGetBoincOpenCLDeviceIndex();
        if (flags->useDevNumber < 0)
        {
            flags->useDevNumber = 0;
        }
    }

    if (!flags->preferredPlatformVendor)
    {
        const char* vendor = mwGetBoincOpenCLPlatformVendor();
        if (vendor)
        {
            mw_printf("BOINC GPU type suggests using OpenCL vendor '%s'\n", vendor);
        }
        else
        {
            /* If BOINC doesn't tell us, guess based on the binary name */
            vendor = mwGuessPreferredPlatform(progName);
            if (vendor)
            {
                mw_printf("Guessing preferred OpenCL vendor '%s'\n", vendor);
            }
        }

        flags->preferredPlatformVendor = vendor ? strdup(vendor) : NULL;
    }

    if (flags->targetFrequency <= 0.0)
    {
        flags->targetFrequency = prefs->gpuTargetFrequency;
    }

    if (flags->nonResponsive < 0)
    {
        flags->nonResponsive = prefs->gpuNonResponsive;
    }

    if (flags->processPriority == MW_PRIORITY_INVALID)
    {
        /* For GPU versions, default to using a higher process priority if not set */
        if (SEPARATION_OPENCL && !flags->forceNoOpenCL)
        {
            flags->processPriority = prefs->gpuProcessPriority;
        }
    }
}

/* Read project preferences and command line arguments */
static int parseParameters(int argc, const char** argv, SeparationFlags* sfOut)
{
    poptContext context;
    int argRead;
    static int version = FALSE;
    static int copyright = FALSE;
    static unsigned int numParams = 0;
    static int serverParams = 0;
    static const char** rest = NULL;
    static SeparationFlags sf;

    static const struct poptOption options[] =
        {
            {
                "astronomy-parameter-file", 'a',
                POPT_ARG_STRING, &sf.ap_file,
                0, "Astronomy parameter file", NULL
            },

            {
                "star-points-file", 's',
                POPT_ARG_STRING, &sf.star_points_file,
                0, "Star points files", NULL
            },

            {
                "output", 'o',
                POPT_ARG_STRING, &sf.separation_outfile,
                0, "Output file for separation (enables separation)", NULL
            },

            {
                "seed", 'e',
                POPT_ARG_INT, &sf.separationSeed,
                SEED_ARGUMENT, "Seed for random number generator", NULL
            },

            {
                "modfit", 'f',
                POPT_ARG_NONE, &sf.modfit,
                0, "Modified fit from Newby 2011", NULL
            },

            {
				"newbg", 'y',
				POPT_ARG_NONE, &sf.background,
				0, "Uses broken power law as background fit", NULL
			},

            {
                "ignore-checkpoint", 'i',
                POPT_ARG_NONE, &sf.ignoreCheckpoint,
                0, "Ignore the checkpoint file", NULL
            },

            {
                "cleanup-checkpoint", 'c',
                POPT_ARG_NONE, &sf.cleanupCheckpoint,
                0, "Delete checkpoint on successful", NULL
            },

            {
                "print-likelihood-text", 't',
                POPT_ARG_NONE, &sf.LikelihoodToText,
                0, "Create text file with likelihood for use in local MLE", NULL
            },

            {
                "debug-boinc", 'g',
                POPT_ARG_NONE, &sf.debugBOINC,
                0, "Init BOINC with debugging. No effect if not built with BOINC_APPLICATION", NULL
            },

            {
                "process-priority", 'b',
                POPT_ARG_INT, &sf.processPriority,
                0, "Set process priority. Set priority 0 (lowest) to 4 (highest)", NULL
            },

            {
                "device", 'd',
                POPT_ARG_INT, &sf.useDevNumber,
                0, "Device number passed by BOINC to use", NULL
            },

            {
                "non-responsive", 'r',
                POPT_ARG_NONE, &sf.nonResponsive,
                0, "Do not care about display responsiveness (use with caution)", NULL
            },

            {
                "gpu-target-frequency", 'q',
                POPT_ARG_DOUBLE, &sf.targetFrequency,
                0, "Target frequency for GPU tasks" , NULL
            },

            {
                "gpu-wait-factor", 'w',
                POPT_ARG_DOUBLE, &sf.waitFactor,
                0, "Wait correction factor when using high CPU workarounds" , NULL
            },

            {
                "gpu-polling-mode", 'm',
                POPT_ARG_INT, &sf.pollingMode,
                0, "Interval for polling GPU: (-2 (default): Use mode -1 unless working around high CPU driver issue.  -1: use clWaitForEvents(). 0: Use clWaitForEvents() with initial wait, >= 1: sets manual interval polling in ms)" , NULL
            },

            {
                "gpu-disable-checkpointing", 'k',
                POPT_ARG_NONE, &sf.disableGPUCheckpointing,
                0, "Disable checkpointing with GPUs" , NULL
            },

            {
                "platform", 'l',
                POPT_ARG_INT, &sf.usePlatform,
                0, "CL platform index to use", NULL
            },

            {
                "platform-vendor", '\0',
                POPT_ARG_STRING, &sf.preferredPlatformVendor,
                0, "CL Platform vendor name to try to use", NULL
            },

            {
                "verbose", '\0',
                POPT_ARG_NONE, &sf.verbose,
                0, "Print some extra debugging information", NULL
            },

            {
                "force-no-opencl", '\0',
                POPT_ARG_NONE, &sf.forceNoOpenCL,
                0, "Use regular CPU path instead of OpenCL if available", NULL
            },

            {
                "force-no-il-kernel", '\0',
                POPT_ARG_NONE, &sf.forceNoILKernel,
                0, "Do not use AMD IL replacement kernels if available", NULL
            },

            {
                "force-no-intrinsics", '\0',
                POPT_ARG_NONE, &sf.forceNoIntrinsics,
                0, "Use old default path", NULL
            },

            {
                "force-x87", '\0',
                POPT_ARG_NONE, &sf.forceX87,
                0, "Force to use x87 path (ignored if x86_64)", NULL
            },

            {
                "force-sse2", '\0',
                POPT_ARG_NONE, &sf.forceSSE2,
                0, "Force to use SSE2 path", NULL
            },

            {
                "force-sse3", '\0',
                POPT_ARG_NONE, &sf.forceSSE3,
                0, "Force to use SSE3 path", NULL
            },

            {
                "force-sse4.1", '\0',
                POPT_ARG_NONE, &sf.forceSSE41,
                0, "Force to use SSE4.1 path", NULL
            },

            {
                "force-avx", '\0',
                POPT_ARG_NONE, &sf.forceAVX,
                0, "Force to use AVX path", NULL
            },

            {
                "p", 'p',
                POPT_ARG_NONE, &serverParams,
                0, "Unused dummy argument to satisfy primitive arguments the server sends", NULL
            },

            {
                "np", '\0',
                POPT_ARG_INT | POPT_ARGFLAG_ONEDASH, &numParams,
                0, "Unused dummy argument to satisfy primitive arguments the server sends", NULL
            },

            {
                "version", 'v',
                POPT_ARG_NONE, &version,
                0, "Print version information", NULL
            },

            {
                "copyright", '\0',
                POPT_ARG_NONE, &copyright,
                0, "Print copyright information and exit", NULL
            },

            POPT_AUTOHELP
            POPT_TABLEEND
        };

    setInitialFlags(&sf);

    context = poptGetContext(argv[0], argc, argv, options, POPT_CONTEXT_POSIXMEHARDER);
    if (!context)
    {
        mw_printf("Failed to get popt context\n");
        exit(EXIT_FAILURE);
    }

    if (argc < 2)
    {
        poptPrintUsage(context, stderr, 0);
        poptFreeContext(context);
        exit(EXIT_FAILURE);
    }

    argRead = mwReadArguments(context);
    if (argRead < 0)
    {
        poptFreeContext(context);
        freeSeparationFlags(&sf);
        exit(EXIT_FAILURE);
    }

    if (version)
    {
        printVersion(FALSE, sf.verbose);
    }

    if (copyright)
    {
        printCopyright();
    }

    if (version || copyright)
    {
        exit(EXIT_SUCCESS);
    }

    sf.setSeed = !!(argRead & SEED_ARGUMENT); /* Check if these flags were used */

    sf.do_separation = (sf.separation_outfile && strcmp(sf.separation_outfile, ""));
    if (sf.do_separation)
        prob_ok_init(sf.separationSeed, sf.setSeed);

    rest = poptGetArgs(context);
    sf.forwardedArgs = mwGetForwardedArguments(rest, &sf.nForwardedArgs);
    sf.numArgs = mwReadRestArgs(rest, sf.nForwardedArgs); /* Temporary */

    poptFreeContext(context);
    setDefaults(&sf);
    *sfOut = sf;

    return 0;
}

static IntegralArea* prepareParameters(const SeparationFlags* sf,
                                       AstronomyParameters* ap,
                                       BackgroundParameters* bgp,
                                       Streams* streams)
{
    IntegralArea* ias;

    ias = setupSeparation(ap, bgp, streams, sf);
    /* Try the new file first. If that doesn't work, try the old one. */
    if (!ias)
    {
        mw_printf("Switching to Parameter File '%s'\n", sf->ap_file);
        ias = readParameters(sf->ap_file, ap, bgp, streams);
    }

    if (!ias)
    {
        mw_printf("Failed to read parameters file\n");
        return NULL;
    }


    return ias;
}
//Needs to loop to account for number of WUs being crunched
static int worker(const SeparationFlags* sf)
{
    AstronomyParameters ap;
    BackgroundParameters bgp = EMPTY_BACKGROUND_PARAMETERS;
    Streams streams = EMPTY_STREAMS;
    IntegralArea* ias = NULL;
    StreamConstants* sc = NULL;
    SeparationResults* results = NULL;
    int rc;
    CLRequest clr;

    memset(&ap, 0, sizeof(ap));
    memset(&clr, 0, sizeof(clr));

    ap.modfit = sf->modfit;

    if(sf->background)
    {
    	ap.background_profile = BROKEN_POWER_LAW;
    }
    else
    {
    	ap.background_profile = FAST_HERNQUIST;
    }

    setCLReqFlags(&clr, sf);
    /*Assume we are crunching at least 1 work unit (These numbers will be properly set in prepareParameters when the parameter file is read)*/

    ias = prepareParameters(sf, &ap, &bgp, &streams);
    if (!ias)
        return 1;

    if(sf->nForwardedArgs)
    {
        ap.totalWUs = sf->nForwardedArgs/ap.params_per_workunit;
    }
    else
    {
        ap.totalWUs = 1;
    }
    mw_printf("<number_WUs> %d </number_WUs>\n", ap.totalWUs);
    mw_printf("<number_params_per_WU> %d </number_params_per_WU>\n", ap.params_per_workunit);
    int ignoreCheckpoint = sf->ignoreCheckpoint;
    for(ap.currentWU = 0; ap.currentWU < ap.totalWUs; ap.currentWU++)
    {

        if (sf->numArgs && setParameters(&ap, &bgp, &streams, &(sf->numArgs[ap.params_per_workunit * ap.currentWU]), ap.params_per_workunit))
        {
            mwFreeA(ias);
            freeStreams(&streams);
            return 1;
        }

        rc = setAstronomyParameters(&ap, &bgp);
        if (rc)
        {
            mwFreeA(ias);
            freeStreams(&streams);
            return 1;
        }

        setExpStreamWeights(&ap, &streams);
        sc = getStreamConstants(&ap, &streams);
        if (!sc)
        {
            mw_printf("Failed to get stream constants\n");
            mwFreeA(ias);
            freeStreams(&streams);
            return 1;
        }

        results = newSeparationResults(ap.number_streams);
        int currentWU = ap.currentWU;
        rc = evaluate(results, &ap, ias, &streams, sc, sf->LikelihoodToText, sf->star_points_file,
                  &clr, sf->do_separation, &ignoreCheckpoint, sf->separation_outfile);
        if (rc)
            mw_printf("Failed to calculate likelihood\n");
    }
    
    mwFreeA(ias);
    mwFreeA(sc);
    freeStreams(&streams);
    if(results) freeSeparationResults(results);

    return rc;
}

static int separationInit(int debugBOINC)
{
    int rc;
    MWInitType initType = MW_PLAIN;

  #if DISABLE_DENORMALS
    mwDisableDenormalsSSE();
  #endif

    mwFixFPUPrecision();

    if (debugBOINC)
        initType |= MW_DEBUG;

    if (SEPARATION_OPENCL)
        initType |= MW_OPENCL;

    rc = mwBoincInit(initType);
    if (rc)
        return rc;

    if (BOINC_APPLICATION && mwIsFirstRun())
    {
        /* Print the version, but only once for the workunit */
        printVersion(TRUE, FALSE);
    }

  #if (SEPARATION_OPENCL) && defined(_WIN32)
    /* We need to increase timer resolution to prevent big slowdown on windows when CPU is loaded. */
    mwSetTimerMinResolution();
  #endif /* defined(_WIN32) */

    return 0;
}


#ifdef MILKYWAY_IPHONE_APP
  #define main _iphone_main
#endif

int main(int argc, const char* argv[])
{
    int rc;
    SeparationFlags sf;
    SeparationPrefs preferences;
    const char** argvCopy = NULL;

  #ifdef NDEBUG
    mwDisableErrorBoxes();
  #endif /* NDEBUG */

    argvCopy = mwFixArgv(argc, argv);
    rc = parseParameters(argc, argvCopy ? argvCopy : argv, &sf);
    if (rc)
    {
        if (BOINC_APPLICATION)
        {
            freeSeparationFlags(&sf);
            mwBoincInit(MW_PLAIN);
            parseParameters(argc, argvCopy, &sf);
            printVersion(TRUE, FALSE);
        }

        mw_printf("Failed to parse parameters\n");
        free(argvCopy);
        mw_finish(EXIT_FAILURE);
    }


    rc = separationInit(sf.debugBOINC);
    free(argvCopy);
    if (rc)
        return rc;

    separationReadPreferences(&preferences);
    setFlagsFromPreferences(&sf, &preferences, argv[0]);

    if (sf.processPriority != MW_PRIORITY_INVALID)
    {
        mwSetProcessPriority(sf.processPriority);
    }

    rc = worker(&sf);

    freeSeparationFlags(&sf);

    if (!sf.ignoreCheckpoint && sf.cleanupCheckpoint && rc == 0)
    {
        mw_report("Removing checkpoint file '%s'\n", CHECKPOINT_FILE);
        mw_remove(CHECKPOINT_FILE);
    }

    mw_finish(rc);
    return rc;
}

