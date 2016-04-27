/*
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
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

#ifdef _WIN32
  #if MW_IS_X86_32
    #pragma comment(linker, "/SUBSYSTEM:CONSOLE,5.0")
  #else
    #pragma comment(linker, "/SUBSYSTEM:CONSOLE,5.02")
  #endif
#endif

#include <popt.h>

#include "milkyway_util.h"
#include "nbody.h"
#include "nbody_likelihood.h"
#include "nbody_defaults.h"
#include "milkyway_git_version.h"

#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */

#if NBODY_CRLIBM
  #include <crlibm.h>
#endif /* NBODY_CRLIBM */

#define SEED_ARGUMENT (1 << 1)


const char* nbCommitID = MILKYWAY_GIT_COMMIT_ID;
const char* nbCommitDescribe = MILKYWAY_GIT_DESCRIBE;


static void nbPrintCopyright(void)
{
    mw_printf(
        "Milkyway@Home N-body client %d.%d\n\n"
        "Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.\n"
        "Copyright (c) 2010 Ben Willett\n"
        "Copyright (c) 2010-2012 Matthew Arsenault\n"
        "Copyright (c) 2010-2011 Rensselaer Polytechnic Institute.\n"
        "Copyright (c) 2010 The University of Texas at Austin\n"
        "Copyright (c) 2010 Dr. Martin Burtscher\n"
        "\n"
        "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        "This is free software: you are free to change and redistribute it.\n"
        "There is NO WARRANTY, to the extent permitted by law.\n"
        "\n"
        " Incorporates works covered by the following copyright and\n"
        " permission notices:\n"
        "\n"
        "Copyright (C) 2007, 2008 Mutsuo Saito, Makoto Matsumoto and Hiroshima University\n"
        "Copyright (C) 2000, Intel Corporation\n"
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
        NBODY_VERSION_MAJOR,
        NBODY_VERSION_MINOR
        );
}

static void nbPrintVersion(int boincTag, int verbose)
{
    char versionStr[2048];

    snprintf(versionStr, sizeof(versionStr),
             "%s %u.%02u %s %s %s %s %s, %s",
             NBODY_PROJECT_NAME,
             NBODY_VERSION_MAJOR, NBODY_VERSION_MINOR,
             MILKYWAY_SYSTEM_NAME,
             ARCH_STRING,
             PRECSTRING,
             DENORMAL_STRING,
             NBODY_EXTRAVER,
             NBODY_EXTRALIB);

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

static int nbInit(const NBodyFlags* nbf)
{
    MWInitType initType = MW_PLAIN;

  #ifdef _OPENMP
    initType |= MW_MULTITHREAD;
  #endif

    if (nbf->debugBOINC)
        initType |= MW_DEBUG;

    return mwBoincInit(initType);
}


/* Maybe set up some platform specific issues */
static void nbSpecialSetup()
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
static void nbSetForwardedArguments(NBodyFlags* nbf, const char** args)
{
    nbf->forwardedArgs = mwGetForwardedArguments(args, &nbf->numForwardedArgs);
}

/* Read the command line arguments, and do the inital parsing of the parameter file. */
static mwbool nbReadParameters(const int argc, const char* argv[], NBodyFlags* nbfOut)
{
    int argRead;
    poptContext context;
    const char** rest = NULL;   /* Leftover arguments */
    static int version = FALSE;
    static int copyright = FALSE;
    static NBodyFlags nbf = EMPTY_NBODY_FLAGS;
    static unsigned int numParams = 0, params = 0;

    /* FIXME: There's a small leak of the inputFile from use of
       poptGetNextOpt(). Some mailing list post suggestst that this
       is some kind of semi-intended bug to work around something or other */
    static const struct poptOption options[] =
    {
        {
            "input-file", 'f',
            POPT_ARG_STRING, &nbf.inputFile,
            0, "Input Lua file to read", NULL
        },

        {
            "histoout-file", 'z',
            POPT_ARG_STRING, &nbf.histoutFileName,
            0, "Output histogram file", NULL
        },

        {
            "histogram-file", 'h',
            POPT_ARG_STRING, &nbf.histogramFileName,
            0, "Histogram file", NULL
        },

        {
            "match-histogram", 's',
            POPT_ARG_STRING, &nbf.matchHistogram,
            0, "Only match this histogram against other histogram (requires histogram argument)", NULL
        },

        {
            "output-file", 'o',
            POPT_ARG_STRING, &nbf.outFileName,
            0, "Output file", NULL
        },

      #if 0
        {
            "binary-output", 'B',
            POPT_ARG_NONE, &nbf.outputBinary,
            0, "Write output dump as a binary", NULL
        },
      #endif

        {
            "output-cartesian", 'x',
            POPT_ARG_NONE, &nbf.outputCartesian,
            0, "Output Cartesian coordinates instead of lbR", NULL
        },
        
        {
            "output-lbrcartesian", 'b',
            POPT_ARG_NONE, &nbf.outputlbrCartesian,
            0, "Output both lbr and Cartesian coordinates", NULL
        },

        {
            "timing", 't',
            POPT_ARG_NONE, &nbf.printTiming,
            0, "Print timing of actual run", NULL
        },

        {
            "verify-file", 'v',
            POPT_ARG_NONE, &nbf.verifyOnly,
            0, "Check that the input file is valid only; perform no calculation.", NULL
        },

        {
            "checkpoint", 'c',
            POPT_ARG_STRING, &nbf.checkpointFileName,
            0, "Checkpoint file to use", NULL
        },

        {
            "checkpoint-interval", 'w',
            POPT_ARG_INT, &nbf.checkpointPeriod,
            0, "Period (in seconds) to checkpoint. -1 to disable", NULL
        },

        {
            "gpu-disable-checkpointing", 'k',
            POPT_ARG_NONE, &nbf.disableGPUCheckpointing,
            0, "Disable checkpointing with GPUs" , NULL
        },

        {
            "debug-boinc", 'g',
            POPT_ARG_NONE, &nbf.debugBOINC,
            0, "Init BOINC with debugging. No effect if not built with BOINC_APPLICATION", NULL
        },

        {
            "lua-debug-libraries", 'a',
            POPT_ARG_NONE, &nbf.debugLuaLibs,
            0, "Load extra Lua libraries not normally allowed (e.g. io) ", NULL
        },

        {
            "visualizer", 'u',
            POPT_ARG_NONE, &nbf.visualizer,
            0, "Try to run N-body visualization", NULL
        },

        {
            "visualizer-args", '\0',
            POPT_ARG_STRING, &nbf.visArgs,
            0, "Command line to pass on to visualizer", NULL
        },

        {
            "visualizer-bin", '\0',
            POPT_ARG_STRING, &nbf.graphicsBin,
            0, "Path to visualize", NULL
        },

        {
            "ignore-checkpoint", 'i',
            POPT_ARG_NONE, &nbf.ignoreCheckpoint,
            0, "Ignore the checkpoint file", NULL
        },

        {
            "print-histogram", 'm',
            POPT_ARG_NONE, &nbf.printHistogram,
            0, "Print generated histogram to stderr", NULL
        },

        {
            "nthreads", 'n',
            POPT_ARG_INT, &nbf.numThreads,
            0, "BOINC argument for number of threads. No effect if built without OpenMP", NULL
        },

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
            POPT_ARG_INT, &nbf.seed,
            SEED_ARGUMENT, "seed for PRNG", NULL
        },

        {
            "device", 'd',
            POPT_ARG_INT, &nbf.devNum,
            0, "OpenCL device number", NULL
        },

        {
            "platform", 'p',
            POPT_ARG_INT, &nbf.platform,
            0, "OpenCL platform", NULL
        },

        {
            "disable-opencl", '\0',
            POPT_ARG_NONE, &nbf.noCL,
            0, "Use normal CPU path instead of OpenCL. No effect if not built with OpenCL", NULL
        },

        {
            "non-responsive", 'r',
            POPT_ARG_NONE, &nbf.ignoreResponsive,
            0, "Do not care about display responsiveness (use with caution)", NULL
        },

        {
            "progress", 'P',
            POPT_ARG_NONE, &nbf.reportProgress,
            0, "Print verbose progress information, possibly with curses", NULL
        },

        {
            "no-clean-checkpoint", 'k',
            POPT_ARG_NONE, &nbf.noCleanCheckpoint,
            0, "Do not delete checkpoint on finish", NULL
        },

        {
            "verbose", '\0',
            POPT_ARG_NONE, &nbf.verbose,
            0, "Print some extra debugging information", NULL
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
        return TRUE;
    }

    /* Check for invalid options, and must have the input file or a
     * checkpoint to resume from */
    argRead = mwReadArguments(context);
    if (argRead < 0)
    {
        mw_printf("Failed to read arguments\n");
        poptFreeContext(context);
        return TRUE;
    }

    if (version)
    {
        nbPrintVersion(FALSE, nbf.verbose);
    }

    if (copyright)
    {
        nbPrintCopyright();
    }

    if (version || copyright)
    {
        poptFreeContext(context);
        exit(EXIT_SUCCESS);
    }

    if (!nbf.inputFile && !nbf.checkpointFileName && !nbf.matchHistogram)
    {
        mw_printf("An input file, checkpoint, or matching histogram argument is required\n");
        poptFreeContext(context);
        return TRUE;
    }

    if (nbf.matchHistogram && !nbf.histogramFileName)
    {
        mw_printf("--match-histogram argument requires --histogram-file\n");
        poptFreeContext(context);
        return TRUE;
    }

    nbf.setSeed = !!(argRead & SEED_ARGUMENT);

    rest = poptGetArgs(context);
    if ((params || numParams) && !rest)
    {
        mw_printf("Expected arguments to follow, got 0\n");
    }
    else
    {
        nbSetForwardedArguments(&nbf, rest);
    }

    poptFreeContext(context);

    *nbfOut = nbf;

    return FALSE;
}

static void nbSetDefaultFlags(NBodyFlags* nbf)
{
    /* Use default if checkpoint file not specified */
    mwStringDefault(nbf->checkpointFileName, DEFAULT_CHECKPOINT_FILE);

    mwStringDefault(nbf->graphicsBin, NBODY_GRAPHICS_NAME);

    /* Use a specified seed or time seeding */
    nbf->seed = nbf->setSeed ? nbf->seed : (uint32_t) time(NULL);

    if (nbf->checkpointPeriod == 0)
    {
        nbf->checkpointPeriod = NOBOINC_DEFAULT_CHECKPOINT_PERIOD;
    }

    if (BOINC_APPLICATION && nbf->debugLuaLibs)
    {
        mw_printf("Warning: disabling --lua-debug-libraries\n");
        nbf->debugLuaLibs = FALSE;
    }
}

static void freeNBodyFlags(NBodyFlags* nbf)
{
    free(nbf->inputFile);
    free(nbf->outFileName);
    free(nbf->checkpointFileName);
    free(nbf->histogramFileName);
    free(nbf->histoutFileName);
    free(nbf->matchHistogram);
    free(nbf->forwardedArgs);
    free(nbf->graphicsBin);
    free(nbf->visArgs);
}

static int nbSetNumThreads(int numThreads)
{
  #ifdef _OPENMP
    int nProc = omp_get_num_procs();
    int nBoinc = mwGetBoincNumCPU();

    if (nProc <= 0) /* It's happened before... */
    {
        mw_printf("Number of processors %d is crazy\n", nProc);
        return 1;
    }

    /* If command line argument not given, and BOINC gives us a value use that */
    if (numThreads <= 0 && nBoinc > 0)
    {
        numThreads = nBoinc;
    }

    if (numThreads != 0)
    {
        omp_set_num_threads(numThreads);
        mw_printf("Using OpenMP %d max threads on a system with %d processors\n",
                  omp_get_max_threads(),
                  nProc);
    }
  #endif

    return 0;
}

/* Maximum exit code is 255 which ruins everything even though we want
 * to have or'able errors. */
static int nbStatusToRC(NBodyStatus rc)
{
    unsigned int n = (unsigned int) rc;
    unsigned int shift = 0;

    if (rc == NBODY_SUCCESS || (nbStatusIsWarning(rc) && !nbStatusIsFatal(rc)))
    {
        return 0;
    }

    while (n >> shift)
    {
        ++shift;
    }

    return (int) shift - 1;
}

int main(int argc, const char* argv[])
{
    NBodyFlags nbf;
    int rc = 0;
    const char** argvCopy = mwFixArgv(argc, argv);
    nbSpecialSetup();

    if (nbReadParameters(argc, argvCopy ? argvCopy : argv, &nbf))
    {
        if (BOINC_APPLICATION)
        {
            mwBoincInit(MW_PLAIN);
            nbReadParameters(argc, argvCopy ? argvCopy : argv, &nbf);
            nbPrintVersion(TRUE, FALSE);
        }

        mw_finish(EXIT_FAILURE);
    }

    if (nbInit(&nbf))
    {
        exit(EXIT_FAILURE);
    }

    if (BOINC_APPLICATION && mwIsFirstRun())
    {
        nbPrintVersion(TRUE, FALSE);
    }

    nbSetDefaultFlags(&nbf);
    if (nbSetNumThreads(nbf.numThreads))
    {
        mw_finish(EXIT_FAILURE);
    }

    if (nbf.verifyOnly)
    {
        rc = nbVerifyFile(&nbf);
    }
    else if (nbf.matchHistogram)
    {
        real emd;

        emd = nbMatchHistogramFiles(nbf.histogramFileName, nbf.matchHistogram);
        mw_printf("<search_likelihood>%.15f</search_likelihood>\n", -emd);
        rc = isnan(emd);
    }
    else
    {
        rc = nbMain(&nbf);
        rc = nbStatusToRC(rc);

        if (!nbf.noCleanCheckpoint)
        {
            mw_report("Removing checkpoint file '%s'\n", nbf.checkpointFileName);
            mw_remove(nbf.checkpointFileName);
        }
    }

    fflush(stderr);
    fflush(stdout); /* Odd things happen with the OpenCL one where stdout starts disappearing */


    freeNBodyFlags(&nbf);

    if (BOINC_APPLICATION)
    {
        mw_finish(rc);
    }

    return rc;
}

