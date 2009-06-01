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
#include "util.h"
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

#ifdef MILKYWAY_GPU
	#ifdef COMPUTE_ON_GPU
		#include "../astronomy_gpu/evaluation_gpu.h"
	#endif
	#ifdef COMPUTE_ON_CPU
		#include "evaluation_optimized.h"
	#endif

	#include "../evaluation/simple_evaluator.h"
	#include "../evaluation/evaluator.h"
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
SEARCH_PARAMETERS *s;

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


void worker() {
	/********
		*	READ THE ASTRONOMY PARAMETERS
	 ********/
	ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	int retval = read_astronomy_parameters(ASTRONOMY_PARAMETER_FILENAME, ap);
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

	/********
		*	READ AND SET THE SEARCH PARAMETERS
	 ********/
	init_search_parameters(&s, get_optimized_parameter_count(ap));
	retval = read_search_parameters(SEARCH_PARAMETER_FILENAME, s);
//	fwrite_search_parameters(stdout, s);
	set_astronomy_parameters(ap, s->parameters);
//	printf("read search parameters\n");

#ifdef MILKYWAY_GPU
	double **hessian, *gradient, *direction, *step;
	double minimum_fitness, *minimum;
	double likelihood;
	int result, evaluations_done, i;

	#ifdef COMPUTE_ON_CPU
		init_simple_evaluator(astronomy_evaluate);
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
		for (i = 0; i < ap->number_integrals; i++) {
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

		gpu__initialize(ap->wedge, ap->convolve, ap->number_streams, ap->number_integrals,
				r_steps, r_min, r_step_size,
				mu_steps, mu_min, mu_step_size,
				nu_steps, nu_min, nu_step_size,
				sp->number_stars, sp->stars);

		init_simple_evaluator(gpu__likelihood);
	#endif
	get_step(ap, &step);

	likelihood = evaluate(s->parameters);
	fprintf(stdout, "initial likelihood: %.20lf\n\n", likelihood);

	fwrite_double_array(stdout, "point", s->number_parameters, s->parameters);
	fwrite_double_array(stdout, "step", s->number_parameters, step);

	fprintf(stdout, "\n");

	new_matrix(&hessian, s->number_parameters, s->number_parameters);
	get_hessian__checkpointed(s->number_parameters, s->parameters, step, hessian, hessian_checkpoint_file);

	gradient = (double*)malloc(sizeof(double) * s->number_parameters);
	get_gradient__checkpointed(s->number_parameters, s->parameters, step, gradient, gradient_checkpoint_file);

	newton_step__alloc(s->number_parameters, hessian, gradient, &direction);

	minimum = (double*)malloc(sizeof(double) * s->number_parameters);
	result = line_search(s->parameters, likelihood, direction, s->number_parameters, minimum, &minimum_fitness, &evaluations_done);

	//initial likelihood calculation
	evaluations_done++;
	//gradient calculation
	evaluations_done += 2 * s->number_parameters;
	//hessian calculation
	evaluations_done += (3 * s->number_parameters) + 2 * ((s->number_parameters * s->number_parameters) - s->number_parameters);

	write_gpu_result(OUTPUT_FILENAME, s->number_parameters, hessian, gradient, likelihood, s->parameters, minimum_fitness, minimum, evaluations_done, s->metadata);

	#ifdef COMPUTE_ON_GPU
		gpu__free_constants();
	#endif

	free_matrix(&hessian, s->number_parameters, s->number_parameters);
	free(gradient);
	free(direction);
	free(minimum);
#else
	double likelihood = astronomy_evaluate(s->parameters);
//	printf("calculated likelihood: %lf\n", likelihood);

	/********
		*	RESOLVE THE OUTPUT FILE & WRITE THE RESULT
	 ********/
	write_cpu_result(OUTPUT_FILENAME, s->number_parameters, s->parameters, likelihood, s->metadata);
#endif

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

	retval = boinc_init();
        if (retval) exit(retval);
        worker();
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

const char *BOINC_RCSID_33ac47a071 = "$Id: boinc_astronomy.C,v 1.15 2009/06/01 07:28:17 deselt Exp $";
