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

#include <popt.h>
#include <errno.h>
#include <stdio.h>

#include "separation.h"

static char* boinc_graphics = NULL;
static char* star_points_file = NULL;
static char* astronomy_parameter_file = NULL;
static char* output_file = NULL;


#ifdef _WIN32
void AppInvalidParameterHandler(const wchar_t* expression,
                                const wchar_t* function,
                                const wchar_t* file,
                                unsigned int line,
                                uintptr_t pReserved )
{
    fprintf(stderr, "Invalid parameter detected in function %s. File: %s Line: %d\n", function, file, line);
    fprintf(stderr, "Expression: %s\n", expression);
    // Cause a Debug Breakpoint.
    DebugBreak();
}
#endif /* _WIN32 */


/* Returns the newly allocated array of parameters */
static real* parse_parameters(int argc, const char** argv, int* paramnOut)
{
    poptContext context;
    int o;
    unsigned int i, paramn = 0;
    real* parameters = NULL;
    static unsigned int numParams;
    static int server_params = 0;
    static const char** rest;
    static const struct poptOption options[] =
    {
        {
            "boinc-init-graphics", 'b',
            POPT_ARG_STRING, &boinc_graphics,
            'b', "Argument to boinc_init_graphics", NULL
        },

        {
            "star-points-file", 's',
            POPT_ARG_STRING, &star_points_file,
            's', "Star points files", NULL
        },

        {
            "astronomy-parameter-file", 'a',
            POPT_ARG_STRING, &astronomy_parameter_file,
            'a', "Astronomy parameter file", NULL
        },

        {
            "output", 'o',
            POPT_ARG_STRING, &output_file,
            'o', "Output file", NULL
        },

        {
            "p", 'p',
            POPT_ARG_NONE, &server_params,
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

    while ( ( o = poptGetNextOpt(context)) >= 0 );

    if ( o < -1 )
    {
        poptPrintHelp(context, stderr, 0);
        mw_finish(EXIT_FAILURE);
    }

    MW_DEBUG("Got arguments: "
             "boinc_graphics = '%s' "
             "star_points_file = '%s' "
             "astronomy_parameter_file = '%s' "
             "output_file = '%s'\n",
             boinc_graphics,
             star_points_file,
             astronomy_parameter_file,
             output_file);

    rest = poptGetArgs(context);
    if (rest)
    {
        while (rest[++paramn]);  /* Count number of parameters */

        MW_DEBUG("%u arguments leftover\n", paramn);

        if (server_params && (paramn != numParams))
        {
            poptFreeContext(context);
            fprintf(stderr, "Parameter count mismatch: Expected %u, got %u\n", numParams, paramn);
            mw_finish(EXIT_FAILURE);
        }

        parameters = (real*) mallocSafe(sizeof(real) * paramn);
        errno = 0;
        for (i = 0; i < paramn; ++i)
        {
            parameters[i] = (real) strtod(rest[i], NULL);

            if (errno)
            {
                perror("error parsing command line parameters");
                poptPrintHelp(context, stderr, 0);
                free(parameters);
                poptFreeContext(context);
                mw_finish(EXIT_FAILURE);
            }
        }
    }

    poptFreeContext(context);

    *paramnOut = paramn;
    return parameters;
}

static void cleanup_worker()
{
	free(boinc_graphics);
	free(star_points_file);
	free(astronomy_parameter_file);
	free(output_file);
}

static void worker(int argc, const char** argv)
{
    real* parameters;
    int number_parameters, ap_number_parameters;
    ASTRONOMY_PARAMETERS ap = EMPTY_ASTRONOMY_PARAMETERS;
    BACKGROUND_PARAMETERS bgp = EMPTY_BACKGROUND_PARAMETERS;
    STREAMS streams = EMPTY_STREAMS;

    parameters = parse_parameters(argc, argv, &number_parameters);
    if (!parameters)
    {
        fprintf(stderr, "Could not parse parameters from the command line\n");
        mw_finish(EXIT_FAILURE);
    }

    if (read_parameters(astronomy_parameter_file, &ap, &bgp, &streams))
    {
        fprintf(stderr,
                "Error reading astronomy parameters from file '%s'\n",
                astronomy_parameter_file);
        free(parameters);
        cleanup_worker();
		mw_finish(EXIT_FAILURE);
    }

    ap_number_parameters = get_optimized_parameter_count(&ap, &bgp, &streams);

    if (number_parameters < 1 || number_parameters != ap_number_parameters)
    {
        fprintf(stderr,
                "Error reading parameters: number of parameters from the "
                "command line (%d) does not match the number of parameters "
                "to be optimized in %s (%d)\n",
                number_parameters,
                astronomy_parameter_file,
                ap_number_parameters);

        free(parameters);
        cleanup_worker();
        mw_finish(EXIT_FAILURE);
    }

    set_parameters(&ap, &bgp, &streams, parameters);
    free(parameters);

    real likelihood;
    STREAM_CONSTANTS* sc = init_constants(&ap, &bgp, &streams);
    free_background_parameters(&bgp);

    likelihood = evaluate(&ap, &streams, sc, star_points_file);

    fprintf(stderr, "<search_likelihood> %0.20f </search_likelihood>\n", likelihood);
    fprintf(stderr, "<search_application> %s </search_application>\n", BOINC_APP_VERSION);

    free(sc);
    free_astronomy_parameters(&ap);
    free_streams(&streams);

	cleanup_worker();

    mw_finish(EXIT_SUCCESS);
}

#if BOINC_APPLICATION

static int separation_init(int argc, char** argv)
{
    int rc;

  #if BOINC_DEBUG
    rc = boinc_init_diagnostics(  BOINC_DIAG_DUMPCALLSTACKENABLED
                                | BOINC_DIAG_HEAPCHECKENABLED
                                | BOINC_DIAG_MEMORYLEAKCHECKENABLED);
  #else
    rc = boinc_init();
  #endif /* BOINC_DEBUG */


  #if BOINC_APP_GRAPHICS
    #if defined(_WIN32) || defined(__APPLE__)
    rc = boinc_init_graphics(worker);
    #else
    rc = boinc_init_graphics_lib(worker, argv[0]);
    #endif /*  defined(_WIN32) || defined(__APPLE__) */
  #endif /* BOINC_APP_GRAPHICS */

  #if defined(_WIN32) && COMPUTE_ON_GPU
    //make the windows GPU app have a higher priority
    BOINC_OPTIONS options;
    boinc_options_defaults(options);
    options.normal_thread_priority = 1; // higher priority (normal instead of idle)
    rc = boinc_init_options(&options);
  #endif /* defined(_WIN32) && COMPUTE_ON_GPU */

    return rc;
}

#else

static int separation_init(int argc, char** argv)
{
  #pragma unused(argc)
  #pragma unused(argv)
    return 0;
}

#endif /* BOINC_APPLICATION */



int main(int argc, char** argv)
{
    int rc = separation_init(argc, argv);
    if (rc)
        exit(rc);

    worker(argc, (const char**) argv);

    return rc;
}

#if defined(_WIN32) && !defined(__MINGW32__)
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR Args, int WinMode)
{
    LPSTR command_line;
    char* argv[100];
    int argc;

    command_line = GetCommandLine();
    argc = parse_command_line( command_line, argv );
    return main(argc, argv);
}
#endif /* defined(_WIN32) && !defined(__MINGW__) */

const char BOINC_RCSID_33ac47a071[] = "$Id: boinc_astronomy.C,v 1.24 2010/05/04 04:21:24 deselt Exp $";

