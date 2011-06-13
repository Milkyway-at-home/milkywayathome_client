/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
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
#include <stdio.h>
#include <popt.h>

#include "milkyway_util.h"
#include "nbody.h"
#include "nbody_defaults.h"

#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */

#if NBODY_CRLIBM
  #include <crlibm.h>
#endif /* NBODY_CRLIBM */


#if !BOINC_APPLICATION
static void nbodyPrintVersion() { }

static int nbodyInit(const NBodyFlags* nbf) { return 0; }

#else

static void nbodyPrintVersion()
{
    warn("<search_application>" BOINC_NBODY_APP_VERSION "</search_application>\n");
}

static int nbodyInit(const NBodyFlags* nbf)
{
    MWInitType initType = 0;

  #ifdef _OPENMP
    initType |= MW_MULTITHREAD;
  #endif

    if (nbf->debugBOINC)
        initType |= MW_DEBUG;

    return mwBoincInit(initType);
}

#endif


/* Maybe set up some platform specific issues */
static void specialSetup()
{
    mwDisableErrorBoxes();

  #if ENABLE_CRLIBM
    /* Try to handle inconsistencies with x87. We shouldn't use
     * this. This helps, but there can still be some problems for some
     * values. Sticking with SSE2 is the way to go. */
    crlibm_init();
  #endif

  #ifdef _MSC_VER
    /* FIXME: Also for mingw, but seems to be missing */
    /* Make windows printing be more consistent. For some reason it
     * defaults to printing 3 digits in the exponent. There are still
     * issues where the rounding of the last digit by printf on
     * windows in a small number of cases. */
    _set_output_format(_TWO_DIGIT_EXPONENT);
  #endif /* _WIN32 */
}

/* For automated testing, pass extra arguments on to Lua script */
static void setForwardedArguments(NBodyFlags* nbf, const char** args)
{
    nbf->forwardedArgs = mwGetForwardedArguments(args, &nbf->numForwardedArgs);
}

/* Read the command line arguments, and do the inital parsing of the parameter file. */
static mwbool readParameters(const int argc, const char* argv[], NBodyFlags* nbf)
{
    int argRead;
    poptContext context;
    const char** rest = NULL;   /* Leftover arguments */
    mwbool failed = FALSE;

    unsigned int numParams = 0, params = 0;

    /* FIXME: There's a small leak of the inputFile from use of
       poptGetNextOpt(). Some mailing list post suggestst that this
       is some kind of semi-intended bug to work around something or other */
    const struct poptOption options[] =
    {
        {
            "input-file", 'f',
            POPT_ARG_STRING, &nbf->inputFile,
            0, "Input Lua file to read", NULL
        },

        {
            "histoout-file", 'z',
            POPT_ARG_STRING, &nbf->histoutFileName,
            0, "Output histogram file", NULL
        },

        {
            "histogram-file", 'h',
            POPT_ARG_STRING, &nbf->histogramFileName,
            0, "Histogram file", NULL
        },

        {
            "output-file", 'o',
            POPT_ARG_STRING, &nbf->outFileName,
            0, "Output file", NULL
        },

        {
            "output-cartesian", 'x',
            POPT_ARG_NONE, &nbf->outputCartesian,
            0, "Output Cartesian coordinates instead of lbR", NULL
        },

        {
            "timing", 't',
            POPT_ARG_NONE, &nbf->printTiming,
            0, "Print timing of actual run", NULL
        },

        {
            "verify-file", 'v',
            POPT_ARG_NONE, &nbf->verifyOnly,
            0, "Check that the input file is valid only; perform no calculation.", NULL
        },

        {
            "checkpoint", 'c',
            POPT_ARG_STRING, &nbf->checkpointFileName,
            0, "Checkpoint file to use", NULL
        },

        {
            "clean-checkpoint", 'k',
            POPT_ARG_NONE, &nbf->cleanCheckpoint,
            0, "Cleanup checkpoint after finishing run", NULL
        },

        {
            "checkpoint-interval", 'w',
            POPT_ARG_INT, &nbf->checkpointPeriod,
            0, "Period (in seconds) to checkpoint. -1 to disable", NULL
        },

        {
            "debug-boinc", 'g',
            POPT_ARG_NONE, &nbf->debugBOINC,
            0, "Init BOINC with debugging. No effect if not built with BOINC_APPLICATION", NULL
        },

        {
            "lua-debug-libraries", 'a',
            POPT_ARG_NONE, &nbf->debugLuaLibs,
            0, "Load extra Lua libraries not normally allowed (e.g. io) ", NULL
        },

        {
            "visualizer", 'u',
            POPT_ARG_NONE, &nbf->visualizer,
            0, "Try to run N-body visualization", NULL
        },

        {
            "visualizer-args", '\0',
            POPT_ARG_STRING, &nbf->visArgs,
            0, "Command line to pass on to visualizer", NULL
        },

        {
            "ignore-checkpoint", 'i',
            POPT_ARG_NONE, &nbf->ignoreCheckpoint,
            0, "Ignore the checkpoint file", NULL
        },

        {
            "print-bodies", 'b',
            POPT_ARG_NONE, &nbf->printBodies,
            0, "Print bodies", NULL
        },

        {
            "print-histogram", 'm',
            POPT_ARG_NONE, &nbf->printHistogram,
            0, "Print histogram", NULL
        },

      #ifdef _OPENMP
        {
            "nthreads", 'n',
            POPT_ARG_INT, &nbf->numThreads,
            0, "BOINC argument for number of threads", NULL
        },
      #endif /* _OPENMP */

        {
            "p", 'p',
            POPT_ARG_NONE, &params,
            0, "Unused dummy argument to satisfy primitive arguments the server sends", NULL
        },

        {
            "np", '\0',
            POPT_ARG_INT | POPT_ARGFLAG_ONEDASH, &numParams,
            0, "Unused dummy argument to satisfy primitive arguments the server sends", NULL
        },

        {
            "seed", 'e',
            POPT_ARG_INT, &nbf->setSeed,
            'e', "seed for PRNG", NULL
        },

        POPT_AUTOHELP
        POPT_TABLEEND
    };

    context = poptGetContext(argv[0], argc, argv, options, 0);

    if (argc < 2)
    {
        poptPrintUsage(context, stderr, 0);
        poptFreeContext(context);
        return TRUE;
    }

    /* Check for invalid options, and must have the input file or a
     * checkpoint to resume from */
    argRead = mwReadArguments(context);
    if (argRead < 0 || (!nbf->inputFile && !nbf->checkpointFileName))
    {
        poptPrintHelp(context, stderr, 0);
        failed = TRUE;
    }

    rest = poptGetArgs(context);
    if ((params || numParams) && !rest)
    {
        warn("Expected arguments to follow, got 0\n");
        failed = TRUE;
    }
    else
    {
        setForwardedArguments(nbf, rest);
    }

    poptFreeContext(context);

    return failed;
}

static void setDefaultFlags(NBodyFlags* nbf)
{
   /* Use default if checkpoint file not specified */
    stringDefault(nbf->checkpointFileName, DEFAULT_CHECKPOINT_FILE);
    stringDefault(nbf->histogramFileName,  DEFAULT_HISTOGRAM_FILE);

    /* Specifying output files implies using them */
    if (!nbf->printBodies)
        nbf->printBodies = (nbf->outFileName != NULL);
    if (!nbf->printHistogram)
        nbf->printHistogram = (nbf->histoutFileName != NULL);
    if (nbf->checkpointPeriod == 0)
        nbf->checkpointPeriod = NOBOINC_DEFAULT_CHECKPOINT_PERIOD;
}

static void freeNBodyFlags(NBodyFlags* nbf)
{
    free(nbf->inputFile);
    free(nbf->outFileName);
    free(nbf->checkpointFileName);
    free(nbf->histogramFileName);
    free(nbf->histoutFileName);
    free(nbf->forwardedArgs);
    free(nbf->visArgs);
}

#ifdef _OPENMP
static void setNumThreads(int numThreads)
{
    if (numThreads != 0)
    {
        omp_set_num_threads(numThreads);
        mw_report("Using OpenMP %d max threads on a system with %d processors\n",
                  omp_get_max_threads(),
                  omp_get_num_procs());
    }
}
#else

static void setNumThreads(int numThreads) { }

#endif /* _OPENMP */

int main(int argc, const char* argv[])
{
    NBodyFlags nbf = EMPTY_NBODY_FLAGS;
    int rc = 0;

    specialSetup();

    if (readParameters(argc, argv, &nbf))
        exit(EXIT_FAILURE);
    free(argvCopy);

    if (nbodyInit(&nbf))
    {
        exit(EXIT_FAILURE);
    }

    nbodyPrintVersion();
    setDefaultFlags(&nbf);
    setNumThreads(nbf.numThreads);

    if (nbf.verifyOnly)
    {
        rc = verifyFile(&nbf);
    }
    else
    {
        rc = runNBodySimulation(&nbf);
        if (nbf.cleanCheckpoint)
        {
            mw_report("Removing checkpoint file '%s'\n", nbf.checkpointFileName);
            mw_remove(nbf.checkpointFileName);
        }
    }

    freeNBodyFlags(&nbf);
    mw_finish(rc);

    return rc;
}

