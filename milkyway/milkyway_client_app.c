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

#define BREAKPOINT() (__asm__ __volatile__ ("int  $03"))

#include "milkyway.h"
#include "parameters.h"
#include "star_points.h"
#include "evaluation_state.h"
//#include "../searches/search_parameters.h"
//#include "../searches/result.h"
//#include "../evaluation/simple_evaluator.h"
//#include "../evaluation/evaluator.h"

#if MILKYWAY_GPU
	#if COMPUTE_ON_GPU
        void init_constants(ASTRONOMY_PARAMETERS *ap);
        #include "evaluation_gpu.h"
    #endif
	#if COMPUTE_ON_CPU
		#include "evaluation_optimized.h"
	#endif

//#include "../searches/hessian.h"
//	#include "../searches/gradient.h"
//	#include "../searches/newton_method.h"
//	#include "../searches/line_search.h"
	#include "../util/matrix.h"
	#include "../util/io_util.h"

	#define hessian_checkpoint_file "hessian_checkpoint"
	#define gradient_checkpoint_file "gradient_checkpoint"
#else
	#include "evaluation_optimized.h"
#endif

#if USE_OCL
#include "evaluation_ocl.h"
#endif

static char* boinc_graphics = NULL;
static char* search_parameter_file = NULL;
static char* star_points_file = NULL;
static char* astronomy_parameter_file = NULL;
static char* output_file = NULL;


#ifdef _WIN32
	void AppInvalidParameterHandler(const wchar_t* expression, const wchar_t* function, const wchar_t* file, unsigned int line, uintptr_t pReserved ) {
		fprintf(stderr, "Invalid parameter detected in function %s. File: %s Line: %d\n", function, file, line);
		fprintf(stderr, "Expression: %s\n", expression);
		// Cause a Debug Breakpoint.
		DebugBreak();
	}
#endif


extern ASTRONOMY_PARAMETERS *ap;
extern STAR_POINTS *sp;
extern EVALUATION_STATE *es;

void init_simple_evaluator(double (*likelihood_function)(double*));
extern double (*evaluate)(double*);


#if COMPUTE_ON_CPU
double astronomy_evaluate(double *parameters) {
	int retval;
	/**
	 * 	Reset the evaluation state
	 */
	set_astronomy_parameters(ap, parameters);
	reset_evaluation_state(es);
	/********
		*	CALCULATE THE INTEGRALS
	 ********/
	retval = calculate_integrals(ap, es, sp);
	if (retval) {
		fprintf(stderr, "APP: error calculating integrals: %d\n", retval);
		boinc_finish(retval);
	}
//	printf("calculated integrals: %lf, %lf\n", es->background_integral, es->stream_integrals[0]);

	/********
		*	CALCULATE THE LIKELIHOOD
	 ********/
	retval = calculate_likelihood(ap, es, sp);
	if (retval) {
		fprintf(stderr, "APP: error calculating likelihood: %d\n", retval);
		boinc_finish(retval);
	}
	return es->prob_sum / (sp->number_stars - es->bad_jacobians);
}
#endif /* COMPUTE_ON_CPU */


/* Returns the newly allocated array of parameters */
static double* parse_parameters(int argc, const char** argv, int* paramnOut)
{
    poptContext context;
    int i, o, paramn;
    double* parameters = NULL;
    static const char** rest;

    static const struct poptOption options[] =
    {
        { "boinc-init-graphics", 'b',
          POPT_ARG_STRING, &boinc_graphics,
          'b', "Agrument to boinc_init_graphics", NULL
        },

        { "search-parameter-file", 's',
          POPT_ARG_STRING, &search_parameter_file,
          's', "Search parameter file name", NULL
        },

        { "star-points-file", 'p',
          POPT_ARG_STRING, &star_points_file,
          'p', "Star points files", NULL
        },

        { "astronomy-parameter-file", 'a',
          POPT_ARG_STRING, &astronomy_parameter_file,
          'a', "Astronomy parameter file", NULL
        },

        { "output", 'o',
          POPT_ARG_STRING, &output_file,
          'o', "Output file", NULL
        },

        POPT_AUTOHELP

        { NULL, 0, 0, NULL, 0, NULL, NULL }
    };

    context = poptGetContext(argv[0],
                             argc,
                             (const char**) argv,
                             options,
                             POPT_CONTEXT_POSIXMEHARDER);

    while ( ( o = poptGetNextOpt(context)) >= 0 );

    if ( o < -1 )
    {
        poptPrintHelp(context, stderr, 0);
        return NULL;
    }

    rest = poptGetArgs(context);
    if (rest)
    {
        paramn = 0;
        while (rest[++paramn]) ;

        parameters = (double*) malloc(sizeof(double) * paramn);

        errno = 0;
        for ( i = 0; i < paramn; ++i )
        {
            parameters[i] = strtod(rest[i], NULL);

            if (errno)
            {
                free(parameters);
                perror("error parsing command line parameters");
                return NULL;
            }
        }
    }

    poptFreeContext(context);

    *paramnOut = paramn;
    return parameters;
}

static void cleanup_worker()
{
    free_state(es);
    free(es);
    free_parameters(ap);
    free(ap);
    free_star_points(sp);
    free(sp);
}

static void fail_worker()
{
    cleanup_worker();
#ifdef _WIN32
    _set_printf_count_output( 1 );
    _get_printf_count_output();
#endif

}

static void worker(int argc, const char** argv)
{
    double* parameters;
    int ret1, ret2;
    int number_parameters, ap_number_parameters;

    parameters = parse_parameters(argc, argv, &number_parameters);

    if (!parameters)
    {
        fprintf(stderr, "Could not parse parameters from the command line\n");
        boinc_finish(0);
        return;
    }

    ap = (ASTRONOMY_PARAMETERS*) malloc(sizeof(ASTRONOMY_PARAMETERS));
    sp = (STAR_POINTS*) malloc(sizeof(STAR_POINTS));
    es = (EVALUATION_STATE*) malloc(sizeof(EVALUATION_STATE));

    ret1 = read_astronomy_parameters(astronomy_parameter_file, ap);
    ret2 = read_star_points(star_points_file, sp);

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
        fail_worker();
        boinc_finish(1);
    }

    initialize_state(ap, sp, es);

    ap_number_parameters = get_optimized_parameter_count(ap);

    if (number_parameters < 1 || number_parameters != ap_number_parameters)
    {
        fprintf(stderr,
                "Error reading parameters: number of parameters from the "
                 "command line (%d) does not match the number of parameters "
                 "to be optimized in astronomy_parameters.txt (%d)\n",
                number_parameters,
                ap_number_parameters);

        fail_worker();
        boinc_finish(0);
        return;
    }

    set_astronomy_parameters(ap, parameters);

    #if COMPUTE_ON_CPU
        init_simple_evaluator(astronomy_evaluate);
	#endif

	#if COMPUTE_ON_GPU
		int *r_steps = (int*)malloc(ap->number_integrals * sizeof(int));
 		int *mu_steps = (int*)malloc(ap->number_integrals * sizeof(int));
		int *nu_steps = (int*)malloc(ap->number_integrals * sizeof(int));
		double *r_min = (double*)malloc(ap->number_integrals * sizeof(double));
		double *mu_min = (double*)malloc(ap->number_integrals * sizeof(double));
		double *nu_min = (double*)malloc(ap->number_integrals * sizeof(double));
		double *r_step_size = (double*)malloc(ap->number_integrals * sizeof(double));
		double *mu_step_size = (double*)malloc(ap->number_integrals * sizeof(double));
		double *nu_step_size = (double*)malloc(ap->number_integrals * sizeof(double));
		for (int i = 0; i < ap->number_integrals; ++i)
        {
			r_steps[i] = ap->integral[i]->r_steps;
			mu_steps[i] = ap->integral[i]->mu_steps;
			nu_steps[i] = ap->integral[i]->nu_steps;
			r_min[i] = ap->integral[i]->r_min;
			mu_min[i] = ap->integral[i]->mu_min;
			nu_min[i] = ap->integral[i]->nu_min;
			r_step_size[i] = ap->integral[i]->r_step_size;
			mu_step_size[i] = ap->integral[i]->mu_step_size;
			nu_step_size[i] = ap->integral[i]->nu_step_size;
		}
		init_constants(ap);
		gpu__initialize();

		init_simple_evaluator(gpu__likelihood);

    #endif /* COMPUTE_ON_GPU */

#if USE_OCL
    printf("Using OpenCL CodePath\n");
    init_constants(ap);
    ocl_mem_t *ocl_mem = setup_ocl(ap, sp);
	double likelihood = ocl_likelihood(parameters,
					   ap,
					   sp,
					   ocl_mem);
	destruct_ocl(ocl_mem);
#else
	double likelihood = evaluate(parameters) - 3.0;
#endif /* USE_OCL */

	fprintf(stderr,"<search_likelihood> %0.20f </search_likelihood>\n", likelihood);
	fprintf(stderr,"<search_application> %s %s </search_application>\n", BOINC_APP_VERSION, PRECISION);

	free(parameters);
	free_state(es);
	free(es);
	free_parameters(ap);
	free(ap);
	free_star_points(sp);
	free(sp);

#if COMPUTE_ON_GPU
    gpu__free_constants();

    free(r_steps);
    free(mu_steps);
    free(nu_steps);
    free(r_min);
    free(mu_min);
    free(nu_min);
    free(r_step_size);
    free(mu_step_size);
    free(nu_step_size);
#endif /* COMPUTE_ON_GPU */

#ifdef _WIN32
	_set_printf_count_output( 1 );
	_get_printf_count_output();
#endif

	boinc_finish(0);

}

int main(int argc, char **argv)
{
    int retval = 0;

#if BOINC_APP_GRAPHICS
#if defined(_WIN32) || defined(__APPLE__)
    retval = boinc_init_graphics(worker);
#else
    retval = boinc_init_graphics_lib(worker, argv[0]);
#endif
    if (retval) exit(retval);
#endif /* BOINC_APP_GRAPHICS */

#if defined(_WIN32) && COMPUTE_ON_GPU
	//make the windows GPU app have a higher priority
	BOINC_OPTIONS options;
	boinc_options_defaults(options);
	options.normal_thread_priority = 1; // higher priority (normal instead of idle)
	retval = boinc_init_options(&options);
#else
	retval = boinc_init();
#endif
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
	char *project_prefs = NULL;  //init_data.project_preferences;
	if (choose_gpu(argc, argv) == -1)
	  {
	    fprintf(stderr, "Unable to find a capable GPU\n");
	    exit(EXIT_FAILURE);
	  }
	printf("got here\n");
	parse_prefs(project_prefs);
#endif /* COMPUTE_ON_GPU */

    worker(argc, (const char**) argv);

    return retval;
}

#ifdef _WIN32
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR Args, int WinMode){
        LPSTR command_line;
        char* argv[100];
	int argc;

	command_line = GetCommandLine();
	argc = parse_command_line( command_line, argv );
	return main(argc, argv);
}
#endif

const char *BOINC_RCSID_33ac47a071 = "$Id: boinc_astronomy.C,v 1.24 2010/05/04 04:21:24 deselt Exp $";

