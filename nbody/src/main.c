/* ************************************************************************** */
/* CODE.C: hierarchical N-body code. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

/* TODO: wuh wuh windows */
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <popt.h>
#include <errno.h>

#include <json/json.h>

#include "nbody_config.h"

#if ENABLE_CRLIBM
  #include <crlibm.h>
#endif /* ENABLE_CRLIBM */

#include "nbody.h"

#ifdef _WIN32
  #define R_OK 4 /* FIXME: Windows */
#endif


#define DEFAULT_CHECKPOINT_FILE "nbody_checkpoint"
#define DEFAULT_HISTOGRAM_FILE  "histogram"
#define DEFAULT_HISTOUT_FILE    "histout"

/* If one of these options is null, use the default. */
#define stringDefault(s, d) ((s) = (s) ? (s) : strdup((d)))


#if !BOINC_APPLICATION
static void nbodyBoincInit() { }
#else

#if BOINC_DEBUG
/* Use BOINC, but prevent it from redirecting stderr to a file, which
 * is really annoying for debugging */
static void nbodyBoincInit()
{
    int rc =  boinc_init_diagnostics(  BOINC_DIAG_DUMPCALLSTACKENABLED
                                     | BOINC_DIAG_HEAPCHECKENABLED
                                     | BOINC_DIAG_MEMORYLEAKCHECKENABLED);
    if (rc)
    {
        warn("boinc_init failed: %d\n", rc);
        exit(EXIT_FAILURE);
    }
}

#else
/* For BOINC releases */
static void nbodyBoincInit()
{
    int rc = boinc_init();
    if (rc)
    {
        warn("boinc_init failed: %d\n", rc);
        exit(EXIT_FAILURE);
    }
}
#endif /* BOINC_DEBUG */
#endif /* !BOINC_APPLICATION */


/* Maybe set up some platform specific issues */
static void specialSetup()
{
  #if !defined(__SSE2__) && ENABLE_CRLIBM
    /* Try to handle inconsistencies with x87. We shouldn't use
     * this. This helps, but there can still be some problems for some
     * values. Sticking with SSE2 is the way to go. */
    crlibm_init();
  #endif

  #ifdef _WIN32
    /* Make windows printing be more consistent. For some reason it
     * defaults to printing 3 digits in the exponent. There are still
     * issues where the rounding of the last digit by printf on
     * windows in a small number of cases. */
    _set_output_format(_TWO_DIGIT_EXPONENT);
  #endif /* _WIN32 */
}

#if BOINC_APPLICATION

/* Read JSON from a file using BOINC file functions */
static json_object* nbodyJSONObjectFromFile(char* inputFile)
{
    char resolvedPath[1024];
    int ret;
    char* buf;
    json_object* obj;

    ret = boinc_resolve_filename(inputFile, resolvedPath, sizeof(resolvedPath));
    if (ret)
        fail("Error resolving file '%s': %d\n", inputFile, ret);

    buf = nbodyReadFile(resolvedPath);
    obj = json_tokener_parse(buf);
    free(buf);

    return obj;
}

#else

static json_object* nbodyJSONObjectFromFile(char* inputFile)
{
    return json_object_from_file(inputFile);
}

#endif /* BOINC_APPLICATION */

/* Read the command line arguments, and do the inital parsing of the parameter file. */
static json_object* readParameters(const int argc,
                                   const char** argv,
                                   FitParams* fitParams,
                                   char** outFileName,
                                   char** checkpointFileName,
                                   char** histogramFileName,
                                   char** histoutFileName,
                                   int* setSeed,
                                   int* ignoreCheckpoint,
                                   int* outputCartesian,
                                   int* printTiming,
                                   int* verifyOnly)
{
  #if !BOINC_APPLICATION
    #pragma unused(checkpointFileName)
    #pragma unused(ignoreCheckpoint)
  #endif

    poptContext context;
    int o;
    static char* inputFile   = NULL;   /* input JSON file */
    static char* inputStr    = NULL;   /* a string of JSON to use directly */
    json_object* obj         = NULL;
    static const char** rest = NULL;   /* Leftover arguments */

    unsigned int numParams = 0, params = 0, paramCount = 0;

    /* FIXME: There's a small leak of the inputFile from use of
       poptGetNextOpt(). Some mailing list post suggestst that this
       is some kind of semi-intended bug to work around something or
       other while maintaining ABI compatability */
    const struct poptOption options[] =
    {
        {
            "input-file", 'f',
            POPT_ARG_STRING, &inputFile,
            0, "Input file to read", NULL
        },

        {
            "histoout-file", 'z',
            POPT_ARG_STRING, histoutFileName,
            0, "Output histogram file", NULL
        },

        {
            "histogram-file", 'h',
            POPT_ARG_STRING, histogramFileName,
            0, "Histogram file", NULL
        },

        {
            "output-file", 'o',
            POPT_ARG_STRING, outFileName,
            0, "Output file", NULL
        },

        {
            "input-string", 's',
            POPT_ARG_STRING, &inputStr,
            0, "Input given as string", NULL
        },

        {
            "output-cartesian", 'x',
            POPT_ARG_NONE, outputCartesian,
            0, "Output Cartesian coordinates instead of lbR", NULL
        },

        {
            "timing", 't',
            POPT_ARG_NONE, printTiming,
            0, "Print timing of actual run", NULL
        },

        {
            "check-file", 'g',
            POPT_ARG_NONE, verifyOnly,
            0, "Check that the input file is valid only; perform no calculation.", NULL
        },

      #if BOINC_APPLICATION
        {
            "checkpoint", 'c',
            POPT_ARG_STRING, checkpointFileName,
            0, "Checkpoint file to use", NULL
        },

        {
            "ignore-checkpoint", 'i',
            POPT_ARG_NONE, ignoreCheckpoint,
            0, "Ignore the checkpoint file", NULL
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

        {  /* FIXME: Only used when using the server arguments. */
            "seed", 'e',
            POPT_ARG_INT, setSeed,
            0, "seed for PRNG", NULL
        },

      #endif /* BOINC_APPLICATION */

        POPT_AUTOHELP

        { NULL, 0, 0, NULL, 0, NULL, NULL }
    };

    context = poptGetContext(argv[0],
                             argc,
                             argv,
                             options,
                             POPT_CONTEXT_POSIXMEHARDER);

    if (argc < 2)
    {
        poptPrintUsage(context, stderr, 0);
        poptFreeContext(context);
        nbody_finish(EXIT_FAILURE);
    }

    while ( ( o = poptGetNextOpt(context)) >= 0 );

    /* Check for invalid options, and must have one of input file or input string */
    if (     o < -1
         ||  (inputFile && inputStr)
         || !(inputFile || inputStr))
    {
        poptPrintHelp(context, stderr, 0);
        free(inputFile);
        free(inputStr);
        poptFreeContext(context);
        nbody_finish(EXIT_FAILURE);
    }

    if (params) /* Using primitive server arguments */
    {
        unsigned int i;
        real* parameters = NULL;

        fitParams->useFitParams = TRUE;

        /* Read through all the server arguments, and make sure we can
         * read everything and have the right number before trying to
         * do anything with them */
        rest = poptGetArgs(context);
        if (!rest)
        {
            poptFreeContext(context);
            fail("Expected arguments to follow, got 0\n");
        }

        while (rest[++paramCount]);  /* Count number of parameters */

        if (numParams == 0)
        {
            poptFreeContext(context);
            fail("numParams = 0 makes no sense\n");
        }

        /* Make sure the number of extra parameters matches the number
         * we were told to expect. */
        if (numParams != paramCount)
        {
            poptFreeContext(context);
            fail("Parameter count mismatch: Expected %u, got %u\n", numParams, paramCount);
        }

        parameters = (real*) mallocSafe(sizeof(real) * numParams);

        errno = 0;
        for ( i = 0; i < numParams; ++i )
        {
            parameters[i] = (real) strtod(rest[i], NULL);
            if (errno)
            {
                perror("Error parsing command line fit parameters");
                poptPrintHelp(context, stderr, 0);
                free(parameters);
                poptFreeContext(context);
                nbody_finish(EXIT_FAILURE);
            }
        }

        fitParams->modelMass        = parameters[0];
        fitParams->modelRadius      = parameters[1];
        fitParams->reverseOrbitTime = parameters[2];
        fitParams->simulationTime   = parameters[3];

        free(parameters);
    }

    poptFreeContext(context);

    if (inputFile)
    {
        /* Check if we can read the file, so we can fail saying that
         * and not be left to guessing if it's that or a parse
         * error */
        if (access(inputFile, R_OK) < 0)
        {
            perror("Failed to read input file");
            free(inputFile);
            nbody_finish(EXIT_FAILURE);
        }

        /* The lack of parse errors from json-c is unfortunate.
           TODO: If we use the tokener directly, can get them.
         */

        obj = nbodyJSONObjectFromFile(inputFile);
        if (is_error(obj))
        {
            warn("Parse error in file '%s'\n", inputFile);
            free(inputFile);
            nbody_finish(EXIT_FAILURE);
        }
        free(inputFile);
    }
    else
    {
        obj = json_tokener_parse(inputStr);
        free(inputStr);
        if (is_error(obj))
            fail("Failed to parse given string\n");
    }

    return obj;
}

/* main: toplevel routine for hierarchical N-body code. */
int main(int argc, const char* argv[])
{
    char* outFile        = NULL;
    json_object* obj     = NULL;
    int outputCartesian  = FALSE;
    int ignoreCheckpoint = FALSE;
    int printTiming      = FALSE;
    int verifyOnly       = FALSE;
    char* checkpointFile = NULL;
    char* histogramFile  = NULL;
    char* histoutFile    = NULL;
    FitParams fitParams  = EMPTY_FIT_PARAMS;
    int setSeed;

    specialSetup();
    nbodyBoincInit();

    obj = readParameters(argc,
                         argv,
                         &fitParams,
                         &outFile,
                         &checkpointFile,
                         &histogramFile,
                         &histoutFile,
                         &setSeed,
                         &ignoreCheckpoint,
                         &outputCartesian,
                         &printTiming,
                         &verifyOnly);

    /* Use default if checkpoint file not specified */
    stringDefault(checkpointFile, DEFAULT_CHECKPOINT_FILE);
    stringDefault(histogramFile,  DEFAULT_HISTOGRAM_FILE);
    stringDefault(histoutFile,    DEFAULT_HISTOUT_FILE);

    /* the PRNG uses a long for a seed, but int is more portable. */
    runNBodySimulation(obj, &fitParams,
                       outFile, checkpointFile, histogramFile, histoutFile, (long) setSeed,
                       outputCartesian, printTiming, verifyOnly);

    free(outFile);
    free(checkpointFile);
    free(histogramFile);
    free(histoutFile);

    nbody_finish(EXIT_SUCCESS);
}

