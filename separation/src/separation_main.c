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

#include "separation.h"
#include "separation_lua.h"
#include "milkyway_util.h"
#include "milkyway_cpp_util.h"
#include "io_util.h"
#include <popt.h>


#define DEFAULT_ASTRONOMY_PARAMETERS "astronomy_parameters.txt"
#define DEFAULT_STAR_POINTS "stars.txt"

#if NVIDIA_OPENCL
  #define DEFAULT_PREFERRED_PLATFORM_VENDOR "NVIDIA Corporation"
#elif AMD_OPENCL
  #define DEFAULT_PREFERRED_PLATFORM_VENDOR "Advanced Micro Devices, Inc."
#elif defined(__APPLE__)
  #define DEFAULT_PREFERRED_PLATFORM_VENDOR "Apple"
#else
  #define DEFAULT_PREFERRED_PLATFORM_VENDOR ""
#endif

#define SEED_ARGUMENT (1 << 1)
#define PRIORITY_ARGUMENT (1 << 2)

static void printCopyright()
{
    warn(
        "Milkyway@Home Separation client %d.%d\n\n"
        "Copyright (c) 2008-2011 Travis Desell, Nathan Cole, Boleslaw Szymanski\n"
        "Copyright (c) 2008-2011 Heidi Newberg, Carlos Varela, Malik Magdon-Ismail\n"
        "Copyright (c) 2008-2011 Rensselaer Polytechnic Institute.\n"
        "Copyright (c) 2010-2011 Matthew Arsenault\n"
        "Copyright (c) 1991-2000 University of Groningen, The Netherlands.\n"
        "Copyright (c) 2001-2009 The GROMACS Development Team\n"
        "\n"
        "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        "This is free software: you are free to change and redistribute it.\n"
        "There is NO WARRANTY, to the extent permitted by law.\n"
        "\n"
        " Incorporates works covered by the following copyright and\n"
        " permission notice:\n"
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
        "\n",
        SEPARATION_VERSION_MAJOR,
        SEPARATION_VERSION_MINOR
        );
}


static void printVersion(int boincTag)
{
    char versionStr[2048];

    snprintf(versionStr, sizeof(versionStr), "%s %u.%u %s %s %s%s%s%s",
             SEPARATION_APP_NAME,
             SEPARATION_VERSION_MAJOR, SEPARATION_VERSION_MINOR,
             SEPARATION_SYSTEM_NAME,
             ARCH_STRING,
             PRECSTRING,
             DENORMAL_STRING,
             SEPARATION_SPECIAL_STR,
             SEPARATION_SPECIAL_LIBM_STR);

    if (boincTag)
    {
        warn("<search_application> %s </search_application>\n", versionStr);
    }
    else
    {
        warn("%s %s\n",
             versionStr,
             BOINC_APPLICATION ? "BOINC" : "");
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

/* Use hardcoded names if files not specified */
static void setDefaultFiles(SeparationFlags* sf)
{
    stringDefault(sf->star_points_file, DEFAULT_STAR_POINTS);
    stringDefault(sf->ap_file, DEFAULT_ASTRONOMY_PARAMETERS);
    stringDefault(sf->preferredPlatformVendor, DEFAULT_PREFERRED_PLATFORM_VENDOR);
}

static void setCLReqFlags(CLRequest* clr, const SeparationFlags* sf)
{
    clr->forceNoIntrinsics = sf->forceNoIntrinsics;
    clr->forceX87 = sf->forceX87;
    clr->forceSSE2 = sf->forceSSE2;
    clr->forceSSE3 = sf->forceSSE3;
    clr->forceSSE41 = sf->forceSSE41;
    clr->verbose = sf->verbose;
    clr->nonResponsive = sf->nonResponsive;
    clr->enableCheckpointing = !sf->disableGPUCheckpointing;

    clr->devNum = sf->useDevNumber;
    clr->platform = sf->usePlatform;
    clr->preferredPlatformVendor = sf->preferredPlatformVendor;

    clr->targetFrequency = sf->targetFrequency <= 0.01 ? DEFAULT_TARGET_FREQUENCY : sf->targetFrequency;
    clr->magicFactor = sf->magicFactor;
}

#if BOINC_APPLICATION

/* If someone has mixed Windows and non-Windows they can set for each*/
#ifdef _WIN32
  #define GPU_PRIORITY_SETTING "gpu_process_priority"
#else
  #define GPU_PRIORITY_SETTING "gpu_process_nice"
#endif /* _WIN32 */

/* If using BOINC try reading a few of the settings from the project
 * preferences. If command line arguments are used, those will
 * override the preferences. The command line arguments will also
 * still work without BOINC */
static void separationReadPreferences(SeparationFlags* sf)
{
    MWAppInitData aid;

    static struct
    {
        double gpuTargetFrequency;
        int gpuNonResponsive;
        int gpuProcessPriority;
        int gpuDisableCheckpoint;
    } prefs;

    static MWProjectPrefs sepPrefs[] =
        {
            { "gpu_target_frequency", MW_PREF_DOUBLE, FALSE, &prefs.gpuTargetFrequency   },
            { "gpu_non_responsive",   MW_PREF_BOOL,   FALSE, &prefs.gpuNonResponsive     },
            { GPU_PRIORITY_SETTING,   MW_PREF_INT,    FALSE, &prefs.gpuProcessPriority   },
            { "no_gpu_checkpoint",    MW_PREF_BOOL,   FALSE, &prefs.gpuDisableCheckpoint },
            END_MW_PROJECT_PREFS
        };

    prefs.gpuTargetFrequency   = DEFAULT_TARGET_FREQUENCY;
    prefs.gpuNonResponsive     = DEFAULT_NON_RESPONSIVE;
    prefs.gpuProcessPriority   = DEFAULT_GPU_PRIORITY;
    prefs.gpuDisableCheckpoint = DEFAULT_DISABLE_GPU_CHECKPOINTING;

    if (mwGetMWAppInitData(&aid))
    {
        warn("Error reading app init data. Project preferences will not be used\n");
    }
    else
    {
        if (aid.projectPrefs)
        {
            mwReadProjectPrefs(sepPrefs, aid.projectPrefs);
        }
    }

    /* Any successfully found setting will be used; otherwise it will get the default */
    sf->targetFrequency = prefs.gpuTargetFrequency;
    sf->waitFactor = prefs.gpuWaitFactor;
    sf->nonResponsive = prefs.gpuNonResponsive;
    sf->processPriority = prefs.gpuProcessPriority;
    sf->pollingMode = prefs.gpuPollingMode;
    sf->disableGPUCheckpointing = prefs.gpuDisableCheckpoint;
}

#endif /* BOINC_APPLICATION */


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
    static SeparationFlags sf = EMPTY_SEPARATION_FLAGS;

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
                "debug-boinc", 'g',
                POPT_ARG_NONE, &sf.debugBOINC,
                0, "Init BOINC with debugging. No effect if not built with BOINC_APPLICATION", NULL
            },

            {
                "process-priority", 'b',
                POPT_ARG_INT, &sf.processPriority,
              #ifndef _WIN32
                PRIORITY_ARGUMENT, "Set process nice value (-20 to 20)", NULL
              #else
                PRIORITY_ARGUMENT, "Set process priority class. Set priority class 0 (lowest) to 4 (highest)", NULL
              #endif /* _WIN32 */
            },

            {
                "device", 'd',
                POPT_ARG_INT, &sf.useDevNumber,
                0, "Device number passed by BOINC to use", NULL
            },

            {
                "magic-factor", 'm',
                POPT_ARG_INT, &sf.magicFactor,
                0, "Number of blocks to run on GPU at once.", NULL
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

  #if BOINC_APPLICATION
    separationReadPreferences(&sf);
  #endif /* BOINC_APPLICATION */

    context = poptGetContext(argv[0], argc, argv, options, POPT_CONTEXT_POSIXMEHARDER);
    if (argc < 2)
    {
        poptPrintUsage(context, stderr, 0);
        poptFreeContext(context);
        exit(EXIT_FAILURE);
    }

    argRead = mwReadArguments(context);
    if (argRead < 0)
    {
        poptPrintHelp(context, stderr, 0);
        poptFreeContext(context);
        freeSeparationFlags(&sf);
        exit(EXIT_FAILURE);
    }

    if (version)
    {
        printVersion(FALSE);
    }

    if (copyright)
    {
        printCopyright();
    }

    if (version || copyright)
    {
        exit(EXIT_SUCCESS);
    }

    sf.setSeed = argRead & SEED_ARGUMENT; /* Check if these flags were used */
    sf.setPriority = argRead & PRIORITY_ARGUMENT;

    sf.do_separation = (sf.separation_outfile && strcmp(sf.separation_outfile, ""));
    if (sf.do_separation)
        prob_ok_init(sf.separationSeed, sf.setSeed);

    rest = poptGetArgs(context);
    sf.forwardedArgs = mwGetForwardedArguments(rest, &sf.nForwardedArgs);
    sf.numArgs = mwReadRestArgs(rest, sf.nForwardedArgs); /* Temporary */

    poptFreeContext(context);
    setDefaultFiles(&sf);
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
        warn("Error reading astronomy parameters from file '%s'\n"
             "  Trying old parameters file\n", sf->ap_file);
        ias = readParameters(sf->ap_file, ap, bgp, streams);
    }

    if (!ias)
    {
        warn("Failed to read parameters file\n");
        return NULL;
    }

    if (sf->numArgs && setParameters(ap, bgp, streams, sf->numArgs, sf->nForwardedArgs))
    {
        mwFreeA(ias);
        freeStreams(streams);
        return NULL;
    }

    return ias;
}

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

    setCLReqFlags(&clr, sf);
    ias = prepareParameters(sf, &ap, &bgp, &streams);
    if (!ias)
        return 1;

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
        warn("Failed to get stream constants\n");
        mwFreeA(ias);
        freeStreams(&streams);
        return 1;
    }

    results = newSeparationResults(ap.number_streams);

    rc = evaluate(results, &ap, ias, &streams, sc, sf->star_points_file,
                  &clr, sf->do_separation, sf->ignoreCheckpoint, sf->separation_outfile);
    if (rc)
        warn("Failed to calculate likelihood\n");

    printSeparationResults(results, ap.number_streams);

    mwFreeA(ias);
    mwFreeA(sc);
    freeStreams(&streams);
    freeSeparationResults(results);

    return rc;
}

static int separationInit(int debugBOINC, MWPriority priority, int setPriority)
{
    int rc;
    MWInitType initType = 0;

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

    /* For GPU versions, default to using a higher process priority if not set */

    /* If a  priority was specified, use that */
    if (setPriority)
    {
        mwSetProcessPriority(priority);
    }
    else if (SEPARATION_OPENCL)
    {
        mwSetProcessPriority(DEFAULT_GPU_PRIORITY);
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
    SeparationFlags sf = EMPTY_SEPARATION_FLAGS;
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
        }

        warn("Failed to parse parameters\n");
        free(argvCopy);
        mw_finish(EXIT_FAILURE);
    }

    rc = separationInit(sf.debugBOINC, sf.processPriority, sf.setPriority);
    free(argvCopy);
    if (rc)
        return rc;

    rc = worker(&sf);

    freeSeparationFlags(&sf);

  #if !SEPARATION_OPENCL
    if (!sf.ignoreCheckpoint && sf.cleanupCheckpoint && rc == 0)
    {
        mw_report("Removing checkpoint file '%s'\n", CHECKPOINT_FILE);
        mw_remove(CHECKPOINT_FILE);
    }
  #endif

  #if BOINC_APPLICATION
    printVersion(TRUE);
    mw_finish(rc);
  #endif

    return rc;
}

