#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "asynchronous_search.h"
#include "asynchronous_newton_method.h"
#include "regression.h"
#include "newton_method.h"
#include "outliers.h"
#include "population.h"
#include "recombination.h"
#include "search_log.h"
#include "search_parameters.h"

#include "../evaluation/search_manager.h"
#include "../util/settings.h"
#include "../util/matrix.h"
#include "../util/io_util.h"

ASYNCHRONOUS_SEARCH* get_asynchronous_genetic_search() {
	ASYNCHRONOUS_SEARCH *as = (ASYNCHRONOUS_SEARCH*)malloc(sizeof(ASYNCHRONOUS_SEARCH));
	as->search_qualifier = (char*)malloc(sizeof(char) * SEARCH_QUALIFIER_SIZE);
	strcpy(as->search_qualifier, "gs");
	as->create_search = create_genetic_search;
	as->read_search = read_genetic_search;
	as->checkpoint_search = checkpoint_genetic_search;
	as->generate_parameters = genetic_search_generate_parameters;
	as->insert_parameters = genetic_search_insert_parameters;
	return as;
}

int create_genetic_search(char* search_name, int number_arguments, char** arguments, int number_parameters, double *point, double *range, double *min_bound, double *max_bound) {
	int type, remove_outliers, maximum_iteration, evaluations_per_iteration, ls_iterations;
	char search_directory[FILENAME_SIZE];
	NEWTON_METHOD_SEARCH *nms;
	int i, simplex_points, population_size;
	double simplex_min, simplex_max;

	remove_outliers = 0;
	type = GENETIC_AVERAGE;
	int max_evaluations = 10000;
	for (i = 0; i < number_arguments; i++) {
		if (!strcmp(arguments[i], "-gs_type")) {
			i++;
			if (!strcmp(arguments[i], "average")) {
				type = GENETIC_AVERAGE;
			} else if (!strcmp(arguments[i], "double_shot")) {
				type = GENETIC_DOUBLE_SHOT;
			} else if (!strcmp(arguments[i], "simplex")) {
				type = GENETIC_SIMPLEX;
				simplex_points = atoi[++i];
				simplex_min = atof[++i];
				simplex_max = atof[++i];
			} else {
				sprintf(AS_MSG, "unknown search type: %s", arguments[i]);
				return AS_CREATE_FAIL;
			}
		} else if (!strcmp(arguments[i], "-gs_remove_outliers")) {
			remove_outliers = 1;
		} else if (!strcmp(arguments[i], "-gs_population_size")) {
			population_size = atoi(arguments[++i]);
		} else if (!strcmp(arguments[i], "-gs_evaluations")) {
			max_evaluations = atoi(arguments[++i]);
		}
	}

	sprintf(search_directory, "%s/%s", get_working_directory(), search_name);
	printf("making directory: %s\n", search_directory);
	mkdir(search_directory, 0777);

	gs = (GENETIC_SEARCH*)malloc(sizeof(GENETIC_SEARCH));

	gs->type = type;
	gs->remove_outliers = remove_outliers;
	gs->current_evaluation = 0;
	gs->population_size = population_size;
	gs->maximum_iteration = max_evaluation;
	gs->number_parameters = number_parameters;

	ss->min_bound = (double*)malloc(sizeof(double) * nms->number_parameters);
	ss->max_bound = (double*)malloc(sizeof(double) * nms->number_parameters);

	memcpy(gs->min_bound, min_bound, sizeof(double) * nms->number_parameters);
	memcpy(gs->max_bound, max_bound, sizeof(double) * nms->number_parameters);

	new_population(gs->population_size, nms->number_parameters, &(nms->population));

	return checkpoint_genetic_search(search_name, gs);
}

int fwrite_genetic_search(FILE* file, GENETIC_SEARCH* gs) {
}

int fread_genetic_search(FILE* file, GENETIC_SEARCH* gs) {
}

int checkpoint_genetic_search(char* search_name, void* search_data) {
	NEWTON_METHOD_SEARCH *nms = (NEWTON_METHOD_SEARCH*)search_data;

	sprintf(AS_MSG, "evaluation: %d/%d, iteration: %d/%d", nms->current_evaluation, nms->evaluations_per_iteration, nms->current_iteration, nms->maximum_iteration);

	if (nms->current_iteration >= nms->maximum_iteration) return AS_CP_OVER;
	else return write_newton_method(search_name, search_data);
}

int bound_parameters(int number_parameters, double *parameters, double *min_bound, double *max_bound) {
	int j;
	for (j = 0; j < number_parameters; j++) {
		if (parameters[j] < min_bound[j]) parameters[j] = min_bound[j];
		if (parameters[j] > max_bound[j]) parameters[j] = max_bound[j];
	}
	return 0;
}

int newton_generate_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS* sp) {
	POPULATION *p;
	NEWTON_METHOD_SEARCH *nms = (NEWTON_METHOD_SEARCH*)search_data;

	nms = (NEWTON_METHOD_SEARCH*)(search_data);
	p = nms->population;

	if (nms->current_iteration >= nms->maximum_iteration) return AS_GEN_OVER;

	if (nms->line_search == NULL || nms->mode == NMS_HESSIAN) {
		range_recombination(p->number_parameters, nms->current_point, nms->parameter_range, sp->parameters);
		sprintf(sp->metadata, "it: %d, ev: %d", nms->current_iteration, nms->current_evaluation);
	} else {
		double point = random_linear_recombination(nms->number_parameters, nms->line_search->min_range, nms->line_search->max_range, nms->current_point, nms->line_search->direction, sp->parameters);
		sprintf(sp->metadata, "ls, point: %.20lf, it: %d, ev: %d", point, nms->current_iteration, nms->current_evaluation);
	}
	if (bound_parameters(sp->number_parameters, sp->parameters, nms->min_bound, nms->max_bound)) return AS_GEN_FAIL;

	return AS_GEN_SUCCESS;
}

int verify(NEWTON_METHOD_SEARCH *nms, SEARCH_PARAMETERS *sp) {
	int i;

	if (isnan(sp->fitness)) {
		return AS_INSERT_FITNESS_NAN;
	} else if (sp->fitness >= -2.8) {
		return AS_INSERT_FITNESS_INVALID;
	}

	for (i = 0; i < nms->number_parameters; i++) {
		if (isnan(sp->parameters[i])) return AS_INSERT_PARAMETERS_NAN;
		if (sp->parameters[i] < nms->min_bound[i] || sp->parameters[i] > nms->max_bound[i]) return AS_INSERT_OUT_OF_BOUNDS;
	}

	if (nms->line_search != NULL && nms->mode == NMS_LINE_SEARCH) {
		int current_iteration, ls_evaluation;
		double point;

		if (3 != sscanf(sp->metadata, "ls, point: %lf, it: %d, ev: %d\n", &point, &current_iteration, &ls_evaluation)) return AS_INSERT_BAD_METADATA;

		if (current_iteration != nms->current_iteration) return AS_INSERT_OUT_OF_ITERATION;
		sp->number_parameters = 1;
		free(sp->parameters);
		sp->parameters = (double*)malloc(sizeof(double));
		sp->parameters[0] = point;
	} else {
		int current_iteration, current_evaluation;

		if (2 != sscanf(sp->metadata, "it: %d, ev: %d", &current_iteration, &current_evaluation)) return AS_INSERT_BAD_METADATA;

		for (i = 0; i < nms->number_parameters; i++) {
			if (sp->parameters[i] < (nms->current_point[i] - nms->parameter_range[i]) || sp->parameters[i] > (nms->current_point[i] + nms->parameter_range[i])) return AS_INSERT_OUT_OF_RANGE;
		}
	}

	if (population_contains(nms->population, sp->fitness, sp->parameters)) return AS_INSERT_NOT_UNIQUE;
	return 0;
}


void get_newton_step(NEWTON_METHOD_SEARCH *nms, double *step, double *step_error) {
	double **hessian, **hessian_error, *gradient, *gradient_error, c, c_error;
	POPULATION *p = nms->population;
	int j, k;

	printf("modifying individuals\n");
	for (j = 0; j < p->size; j++) {
		for (k = 0; k < p->number_parameters; k++) p->individuals[j][k] -= nms->current_point[k];
	}
	printf("modified\n");

	init_matrix(&hessian, nms->number_parameters, nms->number_parameters);
	init_matrix(&hessian_error, nms->number_parameters, nms->number_parameters);
	gradient = (double*)malloc(sizeof(double) * nms->number_parameters);
	gradient_error = (double*)malloc(sizeof(double) * nms->number_parameters);

	printf("initialized gradient and hessian\n");

	parabolic_2d_regression_error(p->size, p->number_parameters, p->individuals, p->fitness, hessian, hessian_error, gradient, gradient_error, &c, &c_error);

	printf("did regression\n");

	for (j = 0; j < p->size; j++) {
		for (k = 0; k < p->number_parameters; k++) p->individuals[j][k] += nms->current_point[k];
	}

	printf("fixed individuals\n");

	newton_step(nms->number_parameters, hessian, gradient, step);
	newton_step(nms->number_parameters, hessian_error, gradient_error, step_error);

	printf("did newton steps\n");
	
	free_matrix(&hessian, nms->number_parameters, nms->number_parameters);
	free_matrix(&hessian_error, nms->number_parameters, nms->number_parameters);
	free(gradient);
	free(gradient_error);
}

void get_line_search_step(char* search_name, NEWTON_METHOD_SEARCH *nms, double *new_center, double *new_range) {
        double a, b, c, a_error, b_error, c_error, error;
        double *individuals;
        int j;
        LINE_SEARCH *ls = nms->line_search;
        POPULATION *p = nms->population;
        
        individuals = (double*)malloc(sizeof(double) * p->size);
        for (j = 0; j < p->size; j++) {
		individuals[j] = p->individuals[j][0];
		printf("\tindividuals[%d]: %.20lf\n", j, individuals[j]);
        }

        parabolic_regression(p->number_parameters, individuals, p->fitness, &a, &b, &c, &a_error, &b_error, &c_error);
	log_printf(search_name, "   y = %.20lfx^2 + %.20lfx + %.20lf\n", a, b, c);
	log_printf(search_name, "maxy = %.20lfx^2 + %.20lfx + %.20lf\n", a+a_error, b+b_error, c+c_error);
	log_printf(search_name, "miny = %.20lfx^2 + %.20lfx + %.20lf\n", a-a_error, b-b_error, c-c_error);

        ls->center = parabolic_center(a, b, c);
	if (ls->center > ls->max_range) ls->center = ls->max_range;
	if (ls->center < ls->min_range) ls->center = ls->min_range;

	error = parabolic_center(a_error, b_error, c_error);

        ls->min_range = ls->center - error;
	ls->max_range = ls->center + error;
         
        free(individuals);
}

int newton_insert_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS* sp) {
	NEWTON_METHOD_SEARCH *nms = (NEWTON_METHOD_SEARCH*)search_data;
	POPULATION *p = nms->population;
	int result, i;

	/********
		*	Insert parameters into population.  If cutoff reached, calculate hessian
		*	and generate new population.
	 ********/

	if (nms->line_search == NULL) {
		sprintf(AS_MSG, "ev: %*d/%*d, it: %*d/%*d, f: %.15lf", 3, nms->current_evaluation, 3, nms->evaluations_per_iteration, 3, nms->current_iteration, 3, nms->maximum_iteration, sp->fitness);
	} else {
		if (nms->mode == NMS_HESSIAN) {
		sprintf(AS_MSG, "ev: %*d/%*d, it: %*d/%*d, f: %.15lf, hessian", 3, nms->current_evaluation, 3, nms->evaluations_per_iteration, 3, nms->current_iteration, 3, nms->maximum_iteration, sp->fitness);
		} else {
			sprintf(AS_MSG, "ev: %*d/%*d, it: %*d/%*d, f: %.15lf, ls", 3, nms->current_evaluation, 3, nms->evaluations_per_iteration, 3, nms->current_iteration, 3, nms->maximum_iteration, sp->fitness);
		}
	}

	if (nms->current_iteration < nms->maximum_iteration) {
		result = verify(nms, sp);
		if (result != 0) return result;

		insert_incremental_info(p, sp->parameters, sp->fitness, sp->host_os, sp->app_version);
		nms->current_evaluation++;

		if (nms->current_evaluation >= nms->evaluations_per_iteration) {
			nms->current_evaluation = 0;

			if (nms->current_iteration < nms->maximum_iteration) {
				double *best_point, best_fitness, average_fitness, worst_fitness, standard_deviation;
				best_point = (double*)malloc(sizeof(double) * nms->number_parameters);

				write_newton_method(search_name, nms);

				if (nms->line_search == NULL || nms->mode == NMS_HESSIAN) {
					double *step, *step_error;
					double scaling_factor = 0.4;

					log_printf(search_name, "iteration: %d/%d\n", nms->current_iteration, nms->maximum_iteration);
					get_population_statistics(p, best_point, &best_fitness, &average_fitness, &worst_fitness, &standard_deviation);
					log_printf(search_name, "[before outliers] best_fitness: %.20lf, average_fitness: %.20lf, worst_fitness: %.20lf, st_dev: %.20lf\n", best_fitness, average_fitness, worst_fitness, standard_deviation);

					remove_outliers(p, 25.0);

					get_population_statistics(p, best_point, &best_fitness, &average_fitness, &worst_fitness, &standard_deviation);
					log_printf(search_name, "[after outliers ] best_fitness: %.20lf, average_fitness: %.20lf, worst_fitness: %.20lf, st_dev: %.20lf\n", best_fitness, average_fitness, worst_fitness, standard_deviation);

					step = (double*)malloc(sizeof(double) * nms->number_parameters);
					step_error = (double*)malloc(sizeof(double) * nms->number_parameters);

					/********
						*	Remove outliers here
					 ********/

					get_newton_step(nms, step, step_error);
					for (i = 0; i < nms->number_parameters; i++) {
						if (nms->line_search == NULL) {
							nms->current_point[i] -= step[i] * scaling_factor;
							if (nms->type == NEWTON_ERROR_RANGE) {
								nms->parameter_range[i] = fabs(step_error[i]);
							} else {
								nms->parameter_range[i] = fabs(step[i] * scaling_factor);
							}
						} else {
							nms->line_search->direction[i] = step[i];
						}
					}

					if (nms->line_search != NULL) {
						nms->mode = NMS_LINE_SEARCH;
						nms->line_search->center = 0;
						nms->line_search->min_range = -2;
						nms->line_search->max_range = 1;
					}

					if (bound_newton_method(nms)) {
						nms->current_iteration = nms->maximum_iteration;
						write_newton_method(search_name, nms);
						sprintf(AS_MSG, "Search completed because of NaN");
						log_printf(search_name, "Search completed because of NaN");
						return AS_INSERT_OVER;
					}

					log_print_double_array(search_name, "current_point", nms->number_parameters, nms->current_point);

					free_population(nms->population);
					free(nms->population);
					if (nms->line_search == NULL) {
						log_print_double_array(search_name, "parameter_range", nms->number_parameters, nms->parameter_range);
						nms->current_iteration++;
						new_population(nms->evaluations_per_iteration, nms->number_parameters, &(nms->population));
					} else {
						log_print_double_array(search_name, "direction", nms->number_parameters, nms->line_search->direction);
						log_printf(search_name, "min_range: %.20lf, center: %.20lf, max_range: %.20lf\n", nms->line_search->min_range, nms->line_search->center, nms->line_search->max_range);
						new_population(nms->evaluations_per_iteration, 1, &(nms->population));
					}

					free(step);
					free(step_error);
				} else {
					log_printf(search_name, "iteration: %d/%d, line_search\n", nms->current_iteration, nms->maximum_iteration);
					get_population_statistics(p, best_point, &best_fitness, &average_fitness, &worst_fitness, &standard_deviation);
					log_printf(search_name, "[before outliers] best_fitness: %.20lf, average_fitness: %.20lf, worst_fitness: %.20lf, st_dev: %.20lf\n", best_fitness, average_fitness, worst_fitness, standard_deviation);

					remove_outliers(p, 25.0);

					get_population_statistics(p, best_point, &best_fitness, &average_fitness, &worst_fitness, &standard_deviation);
					log_printf(search_name, "[after outliers ] best_fitness: %.20lf, average_fitness: %.20lf, worst_fitness: %.20lf, st_dev: %.20lf\n", best_fitness, average_fitness, worst_fitness, standard_deviation);

					for (i = 0; i < nms->number_parameters; i++) {
						nms->current_point[i] = nms->current_point[i] + (nms->line_search->direction[i] * best_point[0]);
						nms->parameter_range[i] = fmin(fabs(nms->line_search->direction[i] * best_point[0]), nms->initial_range[i]);
					}

					log_print_double_array(search_name, "current_point", nms->number_parameters, nms->current_point);
					log_print_double_array(search_name, "parameter_range", nms->number_parameters, nms->parameter_range);

					free_population(nms->population);
					free(nms->population);
					new_population(nms->evaluations_per_iteration, nms->number_parameters, &(nms->population));

					nms->mode = NMS_HESSIAN;
					nms->current_iteration++;

					if (bound_newton_method(nms)) {
						nms->current_iteration = nms->maximum_iteration;
						write_newton_method(search_name, nms);
						sprintf(AS_MSG, "Search completed because of NaN");
						log_printf(search_name, "Search completed because of NaN");
						return AS_INSERT_OVER;
					}
				}
				free(best_point);
				write_newton_method(search_name, nms);
			}
		}
	} else {
		return AS_INSERT_OVER;
	}
	return AS_INSERT_SUCCESS;
}
