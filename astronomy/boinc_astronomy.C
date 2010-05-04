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

/********
        *       Includes for BOINC
********/
#include "util.h"
#ifdef _WIN32
	#include "boinc_win.h"
	#include "str_util.h"
#else
	#include "config.h"
#endif

#ifndef _WIN32
	#include <cstdio>
	#include <cctype>
	#include <ctime>
	#include <cstring>
	#include <cstdlib>
	#include <csignal>
	#include <unistd.h>
#endif

/* #define BOINC_APP_GRAPHICS */

#ifdef BOINC_APP_GRAPHICS
	#include "graphics_api.h"
	#include "graphics_lib.h"
#endif

#include "diagnostics.h"
#include "filesys.h"
#include "boinc_api.h"
#include "mfile.h"
#include "stdio.h"

using std::string;

#define OUTPUT_FILENAME "out"
#define SEARCH_PARAMETER_FILENAME "search_parameters.txt"
#define ASTRONOMY_PARAMETER_FILENAME "astronomy_parameters.txt"
#define STAR_POINTS_FILENAME "stars.txt"

/********
	*	Includes for astronomy
 ********/
#include "parameters.h"
#include "star_points.h"
#include "evaluation_state.h"
#include "../searches/search_parameters.h"
#include "../searches/result.h"
#include "../evaluation/simple_evaluator.h"
#include "../evaluation/evaluator.h"

#ifdef MILKYWAY_GPU
	#ifdef COMPUTE_ON_GPU
void init_constants(ASTRONOMY_PARAMETERS *ap);
		#include "../astronomy_gpu/evaluation_gpu.h"
	#endif
	#ifdef COMPUTE_ON_CPU
		#include "evaluation_optimized.h"
	#endif

	#include "../searches/hessian.h"
	#include "../searches/gradient.h"
	#include "../searches/newton_method.h"
	#include "../searches/line_search.h"
	#include "../util/matrix.h"
	#include "../util/io_util.h"

	#define hessian_checkpoint_file "hessian_checkpoint"
	#define gradient_checkpoint_file "gradient_checkpoint"
#else
	#include "evaluation_optimized.h"
#endif

#ifdef USE_OCL
#include "../astronomy_ocl/evaluation_ocl.h"
#endif

#ifdef _WIN32
	void AppInvalidParameterHandler(const wchar_t* expression, const wchar_t* function, const wchar_t* file, unsigned int line, uintptr_t pReserved ) {
		fprintf(stderr, "Invalid parameter detected in function %s. File: %s Line: %d\n", function, file, line);
		fprintf(stderr, "Expression: %s\n", expression);
		// Cause a Debug Breakpoint.
		DebugBreak();
	}
#endif


ASTRONOMY_PARAMETERS *ap;
STAR_POINTS *sp;
EVALUATION_STATE *es;

#ifdef COMPUTE_ON_CPU
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
#endif

int get_number_parameters(int argc, char** argv) {
	int i;
	for (i = 0; i < argc; i++) {
		if ( !strcmp(argv[i], "-np") ) {
			return atoi(argv[++i]);
		}
	}
	return -1;
}

int get_parameters(int argc, char** argv, int number_parameters, double* parameters) {
	int i, j;
	for (i = 0; i < argc; i++) {
//		fprintf(stderr, "np: %d, parsing parameter: %s\n", number_parameters, argv[i]);

		if ( !strcmp(argv[i], "-p") ) {
			for (j = 0; j < number_parameters; j++) {
//				fprintf(stderr, "i: %d, j: %d, argc: %d, arg: %s\n", i, j, argc, argv[1 + i + j]);
				if (1 + i + j >= argc) return 1;

				parameters[j] = atof(argv[1 + i + j]);
			}
			return 0;
		}
	}

	return 1;
}


void worker(int argc, char** argv) {
	int number_parameters, ap_number_parameters, retval;
	double* parameters;

	/********
		*	READ THE ASTRONOMY PARAMETERS
	 ********/
	ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	retval = read_astronomy_parameters(ASTRONOMY_PARAMETER_FILENAME, ap);
	if (retval) {
		fprintf(stderr, "APP: error reading astronomy parameters: %d\n", retval);
		boinc_finish(1);
	}
//	printf("read astronomy parameters\n");

	/********
		*	READ THE STAR POINTS
	 ********/
	sp = (STAR_POINTS*)malloc(sizeof(STAR_POINTS));
	retval = read_star_points(STAR_POINTS_FILENAME, sp);
	if (retval) {
		fprintf(stderr, "APP: error reading star points: %d\n", retval);
		boinc_finish(1);
	}
//	printf("read star points\n");

	/********
		*	INITIALIZE THE EVALUATION STATE
	 ********/
	es = (EVALUATION_STATE*)malloc(sizeof(EVALUATION_STATE));
	initialize_state(ap, sp, es);
//	printf("read evaluation state\n");

	number_parameters = get_number_parameters(argc, argv);
	ap_number_parameters = get_optimized_parameter_count(ap);

	if (number_parameters < 1 || number_parameters != ap_number_parameters) {
		fprintf(stderr, "error reading parameters, number of parameters from the command line (%d) does not match the number of parameters to be optimized in astronomy_parameters.txt (%d)\n", number_parameters, ap_number_parameters);

		free_state(es);
		free(es);
		free_parameters(ap);
		free(ap);
		free_star_points(sp);
		free(sp);

#ifdef _WIN32
		_set_printf_count_output( 1 );
		_get_printf_count_output();
#endif

		boinc_finish(0);
		return;
	}

	parameters = (double*)malloc(sizeof(double) * number_parameters);
	retval = get_parameters(argc, argv, number_parameters, parameters);

	if (retval) {
		fprintf(stderr, "could not parse parameters from the command line, retval: %d\n", retval);

		free_state(es);
		free(es);
		free_parameters(ap);
		free(ap);
		free_star_points(sp);
		free(sp);

#ifdef _WIN32
		_set_printf_count_output( 1 );
		_get_printf_count_output();
#endif

		boinc_finish(0);
		return;
	}

	set_astronomy_parameters(ap, parameters);

	#ifdef COMPUTE_ON_CPU
		init_simple_evaluator(astronomy_evaluate);
		char *precision = "double";
	#endif
	#ifdef COMPUTE_ON_GPU
		int *r_steps = (int*)malloc(ap->number_integrals * sizeof(int));
 		int *mu_steps = (int*)malloc(ap->number_integrals * sizeof(int));
		int *nu_steps = (int*)malloc(ap->number_integrals * sizeof(int));
		double *r_min = (double*)malloc(ap->number_integrals * sizeof(double));
		double *mu_min = (double*)malloc(ap->number_integrals * sizeof(double));
		double *nu_min = (double*)malloc(ap->number_integrals * sizeof(double));
		double *r_step_size = (double*)malloc(ap->number_integrals * sizeof(double));
		double *mu_step_size = (double*)malloc(ap->number_integrals * sizeof(double));
		double *nu_step_size = (double*)malloc(ap->number_integrals * sizeof(double));
		for (int i = 0; i < ap->number_integrals; i++) {
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

		#ifdef DOUBLE_PRECISION
			char *precision = "double";
		#else
			char *precision = "single";
		#endif
	#endif
#ifdef USE_OCL
			//use OpenCL
			printf("Using OpenCL CodePath\n");
			char precision[] = "single";
			init_constants(ap);
			ocl_mem_t *ocl_mem = setup_ocl(ap, sp);
#endif
	char app_version[256];
	sprintf(app_version, "%s: %1.2lf", BOINC_APP_NAME, BOINC_APP_VERSION);
       	
#ifdef USE_OCL
	double likelihood = ocl_likelihood(parameters, 
					   ap,
					   sp,
					   ocl_mem);
	destruct_ocl(ocl_mem);
#else
	double likelihood = evaluate(parameters) - 3.0;
#endif
	fprintf(stderr,"<search_likelihood> %0.20f </search_likelihood>\n", likelihood);
	fprintf(stderr,"<search_application> %s %s </search_application>\n", app_version, precision);


	free(parameters);

	free_state(es);
	free(es);
	free_parameters(ap);
	free(ap);
	free_star_points(sp);
	free(sp);

	#ifdef COMPUTE_ON_GPU
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
	#endif

#ifdef _WIN32
	_set_printf_count_output( 1 );
	_get_printf_count_output();
#endif

	boinc_finish(0);
}

int main(int argc, char **argv){
        int retval = 0;

#ifdef BOINC_APP_GRAPHICS
#if defined(_WIN32) || defined(__APPLE__)
	        retval = boinc_init_graphics(worker);
#else
                retval = boinc_init_graphics_lib(worker, argv[0]);
#endif
                if (retval) exit(retval);
#endif

#if defined(_WIN32) && defined(COMPUTE_ON_GPU)
	//make the windows GPU app have a higher priority
	BOINC_OPTIONS options; 
	boinc_options_defaults(options); 
	options.normal_thread_priority = 1; // higher priority (normal instead of idle) 
	retval = boinc_init_options(&options); 
	if (retval) exit(retval);
#else
	retval = boinc_init();
        if (retval) exit(retval);
#endif

#ifdef COMPUTE_ON_GPU
	//Choose the GPU to execute on, first look
	//at the command line argument for a 
	//--device 0..n string, then enumerate all CUDA
	//devices on the system and choose the one
	//with double precision support and the most
	//GFLOPS
	APP_INIT_DATA init_data;
	boinc_get_init_data_p(&init_data);
	char *project_prefs = init_data.project_preferences;
	if (choose_gpu(argc, argv) == -1)
	  {
	    fprintf(stderr, "Unable to find a capable GPU\n");
	    exit(1);
	  }
	printf("got here\n");
	parse_prefs(project_prefs);
#endif
        worker(argc, argv);
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
