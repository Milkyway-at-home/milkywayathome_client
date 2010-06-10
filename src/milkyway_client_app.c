/*
Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

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

#include <boinc_api.h>

#if BOINC_APP_GRAPHICS
	#include <graphics_api.h>
	#include <graphics_lib.h>
#endif

/* I'm not sure what the MSVC macro is.
 This only needs the Windows API for the stuff to deal with the
 truly awful Windows API / WinMain, which you only need to deal with
 for visual studio and should be avoided as much as possible. */
#if defined(_WIN32) && !defined(__MINGW32__)
    #include <str_util.h>
#endif


#ifdef __cplusplus  /* Workaround for Windows / boinc issue */
 extern "C" {
#endif

#include "milkyway.h"
#include "milkyway_priv.h"
#include "star_points.h"
#include "evaluation_state.h"

#if USE_CUDA
    #include "evaluation_gpu.h"
#elif USE_OCL
    #include "evaluation_ocl.h"
#endif


#ifdef __cplusplus
 }
#endif

static char* boinc_graphics = NULL;
static char* search_parameter_file = NULL;
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
static double* parse_parameters(int argc, const char** argv, int* paramnOut)
{
    poptContext context;
    int o;
    unsigned int i, paramn = 0;
    double* parameters = NULL;
    static const char** rest;

    static const struct poptOption options[] =
    {
        {
            "boinc-init-graphics", 'b',
            POPT_ARG_STRING, &boinc_graphics,
            'b', "Argument to boinc_init_graphics", NULL
        },

        {
            "search-parameter-file", 's',
            POPT_ARG_STRING, &search_parameter_file,
            's', "Search parameter file name", NULL
        },

        {
            "star-points-file", 'p',
            POPT_ARG_STRING, &star_points_file,
            'p', "Star points files", NULL
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
             "search_parameter_file = '%s' "
             "star_points_file = '%s' "
             "astronomy_parameter_file = '%s' "
             "output_file = '%s'\n",
             boinc_graphics,
             search_parameter_file,
             star_points_file,
             astronomy_parameter_file,
             output_file);

    rest = poptGetArgs(context);
    if (rest)
    {
        while (rest[++paramn]);  /* Count number of parameters */

        MW_DEBUG("%u arguments leftover\n", paramn);

        parameters = (double*) malloc(sizeof(double) * paramn);

        errno = 0;
        for ( i = 0; i < paramn; ++i )
        {
            parameters[i] = strtod(rest[i], NULL);

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
	free(search_parameter_file);
	free(star_points_file);
	free(astronomy_parameter_file);
	free(output_file);
}

static void worker(int argc, const char** argv)
{
    double* parameters;
    int ret1, ret2;
    int number_parameters, ap_number_parameters;
    ASTRONOMY_PARAMETERS ap = { 0 };
    STAR_POINTS sp = { 0 };
    EVALUATION_STATE es = { 0 };

    parameters = parse_parameters(argc, argv, &number_parameters);

    if (!parameters)
    {
        fprintf(stderr, "Could not parse parameters from the command line\n");
        mw_finish(EXIT_FAILURE);
    }

    ret1 = read_astronomy_parameters(astronomy_parameter_file, &ap);
    ret2 = read_star_points(star_points_file, &sp);

    MW_DEBUG("ap.number_stream_parameters = %d\n", ap.number_stream_parameters);

    if (ret1)
    {
        fprintf(stderr,
                "APP: error reading astronomy parameters from file %s: %d\n",
                astronomy_parameter_file,
                ret1);
    }

    if (ret2)
    {
        fprintf(stderr,
                "APP: error reading star points from file %s: %d\n",
                star_points_file,
                ret2);
    }

    if (ret1 | ret2)
    {
        free(parameters);
        cleanup_worker();
		mw_finish(EXIT_FAILURE);
    }

    initialize_state(&ap, &sp, &es);

    ap_number_parameters = get_optimized_parameter_count(&ap);

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

    set_astronomy_parameters(&ap, parameters);

#if COMPUTE_ON_CPU
    init_constants(&ap);
    init_simple_evaluator(cpu_evaluate);
#elif USE_CUDA
    init_constants(&ap);
    init_simple_evaluator(cuda_evaluate);
#elif USE_OCL
    init_constants(&ap);
    init_simple_evaluator(ocl_evaluate);
#else
    #error "Must choose CUDA, OpenCL or CPU"
#endif /* COMPUTE_ON_CPU */

    /* CHECKME: What is this magic 3.0, and why was it being
     * subtracted from CPU and CUDA result, but not OpenCL? */
    double likelihood = evaluate(parameters, &ap, &es, &sp) - 3.0;

    fprintf(stderr, "<search_likelihood> %0.20f </search_likelihood>\n", likelihood);
    fprintf(stderr, "<search_application> %s %s </search_application>\n", BOINC_APP_VERSION, PRECISION);

    free(parameters);
	cleanup_worker();

    mw_finish(EXIT_SUCCESS);

}

int main(int argc, char** argv)
{
    int retval = 0;

#if BOINC_APP_GRAPHICS
  #if defined(_WIN32) || defined(__APPLE__)
      retval = boinc_init_graphics(worker);
  #else
      retval = boinc_init_graphics_lib(worker, argv[0]);
  #endif /*  defined(_WIN32) || defined(__APPLE__) */

  if (retval)
      exit(retval);
#endif /* BOINC_APP_GRAPHICS */

#if defined(_WIN32) && COMPUTE_ON_GPU
    //make the windows GPU app have a higher priority
    BOINC_OPTIONS options;
    boinc_options_defaults(options);
    options.normal_thread_priority = 1; // higher priority (normal instead of idle)
    retval = boinc_init_options(&options);
#else
    retval = boinc_init();
#endif /* defined(_WIN32) && COMPUTE_ON_GPU */
    if (retval)
        exit(retval);


#if COMPUTE_ON_GPU
    //Choose the GPU to execute on, first look
    //at the command line argument for a
    //--device 0..n string, then enumerate all CUDA
    //devices on the system and choose the one
    //with double precision support and the most
    //GFLOPS
    //APP_INIT_DATA init_data;
    //boinc_get_init_data_p(&init_data);
    char* project_prefs = NULL;  //init_data.project_preferences;
    if (choose_gpu(argc, argv) == -1)
    {
        fprintf(stderr, "Unable to find a capable GPU\n");
        mw_finish(EXIT_FAILURE);
    }
    MW_DEBUGMSG("got here\n");
    parse_prefs(project_prefs);
#endif /* COMPUTE_ON_GPU */

    worker(argc, (const char**) argv);

    return retval;
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

