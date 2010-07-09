/* ************************************************************************** */
/* CODE.C: hierarchical N-body code. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

/* TODO: wuh wuh windows */
#include <unistd.h>
#include <fcntl.h>

#include <stdlib.h>
#include <string.h>
#include <popt.h>
#include <errno.h>
#include "nbody.h"
#include "nbody_priv.h"

#ifdef _WIN32
  #define R_OK 2 /* FIXME: Windows */
#endif

#define DEFAULT_CHECKPOINT_FILE "nbody_checkpoint"

/* Read the command line arguments, and do the inital parsing of the parameter file */
static json_object* readParameters(const int argc,
                                   const char** argv,
                                   char** outFileName,
                                   int* useDouble,
                                   char** checkpointFileName,
                                   int* ignoreCheckpoint,
                                   int* outputCartesian,
                                   int* printTiming,
                                   int* verifyOnly)
{
  #if !defined(DYNAMIC_PRECISION)
    #pragma unused(useDouble)
  #endif

  #if !BOINC_APPLICATION
    #pragma unused(checkpointFileName)
    #pragma unused(ignoreCheckpoint)
  #endif

    poptContext context;
    int o;
    static char* inputFile = NULL;        /* input JSON file */
    static char* inputStr  = NULL;        /* a string of JSON to use directly */
    json_object* obj;
    static const char** rest;

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
      #endif /* BOINC_APPLICATION */


      #if DYNAMIC_PRECISION
        {
            "double", 'd',
            POPT_ARG_NONE, useDouble,
            0, "Use double precision", NULL
        },
      #endif

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

        rest = poptGetArgs(context);
        if (rest)
        {
            while (rest[++paramCount]);  /* Count number of parameters */

            assert(numParams == paramCount);

            parameters = (real*) mallocSafe(sizeof(real) * numParams);

            errno = 0;
            for ( i = 0; i < numParams; ++i )
            {
                parameters[i] = (real) strtod(rest[i], NULL);
                if (errno)
                {
                    perror("Error parsing command line parameters");
                    poptPrintHelp(context, stderr, 0);
                    free(parameters);
                    poptFreeContext(context);
                    nbody_finish(EXIT_FAILURE);
                }
            }
        }

        /* TODO: Stuff with parameters */
        free(parameters);
    }

    poptFreeContext(context);

    if (inputFile)
    {
        /* check if we can read the file, so we can fail saying that
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
        obj = json_object_from_file(inputFile);
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

/* Run with double, float, or whatever we have */
static void runSimulationWrapper(json_object* obj,
                                 const char* outFileName,
                                 const char* checkpointFileName,
                                 const int useDouble,
                                 const int outputCartesian,
                                 const int printTiming,
                                 const int verifyOnly)
{
  #if DYNAMIC_PRECISION
    if (useDouble)
    {
        printf("Using double precision\n");
        runNBodySimulation_double(obj, outFileName, checkpointFileName, outputCartesian, printTiming, verifyOnly);
        printf("Done with double\n");
    }
    else
    {
        printf("Using float precision\n");
        runNBodySimulation_float(obj, outFileName, checkpointFileName, outputCartesian, printTiming, verifyOnly);
        printf("Done with float\n");
    }
  #else
    #pragma unused(useDouble)

    runNBodySimulation(obj, outFileName, checkpointFileName, outputCartesian, printTiming, verifyOnly);
  #endif /* DYNAMIC_PRECISION */
}

/* FIXME: Clean up the separation between boinc and nonboinc. Right
 * now it's absolutely disgusting. */

/* main: toplevel routine for hierarchical N-body code. */
int main(int argc, const char* argv[])
{
    char* outFileName = NULL;
    json_object* obj = NULL;
    int useDouble = FALSE;
    int outputCartesian = FALSE;
    int ignoreCheckpoint = FALSE;
    int printTiming = FALSE;
    int verifyOnly = FALSE;
    char* checkpointFileName = NULL;

#if BOINC_APPLICATION
    int boincInitStatus = 0;
  #if !BOINC_DEBUG
    boincInitStatus = boinc_init();
  #else
    boincInitStatus = boinc_init_diagnostics(  BOINC_DIAG_DUMPCALLSTACKENABLED
                                             | BOINC_DIAG_HEAPCHECKENABLED
                                             | BOINC_DIAG_MEMORYLEAKCHECKENABLED);
  #endif /* !BOINC_DEBUG */

    if (boincInitStatus)
    {
        fprintf(stderr, "boinc_init failed: %d\n", boincInitStatus);
        exit(EXIT_FAILURE);
    }
#endif /* BOINC_APPLICATION */

    obj = readParameters(argc,
                         argv,
                         &outFileName,
                         &useDouble,
                         &checkpointFileName,
                         &ignoreCheckpoint,
                         &outputCartesian,
                         &printTiming,
                         &verifyOnly);

    /* Use default if checkpoint file not specified */
    checkpointFileName = checkpointFileName ? checkpointFileName : strdup(DEFAULT_CHECKPOINT_FILE);

    runSimulationWrapper(obj,
                         outFileName,
                         checkpointFileName,
                         useDouble,
                         outputCartesian,
                         printTiming,
                         verifyOnly);

    free(outFileName);
    free(checkpointFileName);

    nbody_finish(EXIT_SUCCESS);
}

