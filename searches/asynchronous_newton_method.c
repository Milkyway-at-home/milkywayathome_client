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

ASYNCHRONOUS_SEARCH* get_asynchronous_newton_method() {
	ASYNCHRONOUS_SEARCH *as = (ASYNCHRONOUS_SEARCH*)malloc(sizeof(ASYNCHRONOUS_SEARCH));
	as->search_qualifier = (char*)malloc(sizeof(char) * SEARCH_QUALIFIER_SIZE);
	strcpy(as->search_qualifier, "nm");
	as->create_search = create_newton_method;
	as->read_search = read_newton_method;
	as->checkpoint_search = checkpoint_newton_method;
	as->generate_parameters = newton_generate_parameters;
	as->insert_parameters = newton_insert_parameters;
	return as;
}

int create_newton_method(char* search_name, int number_arguments, char** arguments, int number_parameters, double *point, double *range, double *min_bound, double *max_bound) {
	int type, remove_outliers, maximum_iteration, evaluations_per_iteration, ls_iterations;
	char search_directory[FILENAME_SIZE];
	NEWTON_METHOD_SEARCH *nms;
	int i;

	remove_outliers = 0;
	type = NEWTON_UPDATE_RANGE;
	evaluations_per_iteration = 400;
	maximum_iteration = 10;
	ls_iterations = 2;
	for (i = 0; i < number_arguments; i++) {
		if (!strcmp(arguments[i], "-nm_type")) {
			i++;
			if (!strcmp(arguments[i], "line_search")) {
				type = NEWTON_LINE_SEARCH;
			} else if (!strcmp(arguments[i], "error_range")) {
				type = NEWTON_ERROR_RANGE;
			} else if (!strcmp(arguments[i], "update_range")) {
				type = NEWTON_UPDATE_RANGE;
			} else {
				sprintf(AS_MSG, "unknown search type: %s", arguments[i]);
				return AS_CREATE_FAIL;
			}
		} else if (!strcmp(arguments[i], "-nm_remove_outliers")) {
			remove_outliers = 1;
		} else if (!strcmp(arguments[i], "-nm_evaluations")) {
			evaluations_per_iteration = atoi(arguments[++i]);
		} else if (!strcmp(arguments[i], "-nm_iterations")) {
			maximum_iteration = atoi(arguments[++i]);
		}
	}

	sprintf(search_directory, "%s/%s", get_working_directory(), search_name);
	printf("making directory: %s\n", search_directory);
	mkdir(search_directory, 0777);

	nms = (NEWTON_METHOD_SEARCH*)malloc(sizeof(NEWTON_METHOD_SEARCH));

	nms->type = type;
	nms->remove_outliers = remove_outliers;
	nms->current_iteration = 0;
	nms->maximum_iteration = maximum_iteration;
	nms->current_evaluation = 0;
	nms->evaluations_per_iteration = evaluations_per_iteration;
	nms->number_parameters = number_parameters;
	nms->mode = NMS_HESSIAN;

	if (type == NEWTON_LINE_SEARCH) {
		nms->line_search = (LINE_SEARCH*)malloc(sizeof(LINE_SEARCH));
		nms->line_search->center = 0;
		nms->line_search->max_range = 0;
		nms->line_search->min_range = 0;
		nms->line_search->direction = (double*)malloc(sizeof(double) * nms->number_parameters);
	} else {
		nms->line_search = NULL;
	}

	printf("nms->type: %d\nnms->maximum_iteration: %d\nnms->evaluations_per_iteration: %d\nnms->number_parameters: %d\n", nms->type, nms->maximum_iteration, nms->evaluations_per_iteration, nms->number_parameters);

	nms->current_point = (double*)malloc(sizeof(double) * nms->number_parameters);
	nms->initial_range = (double*)malloc(sizeof(double) * nms->number_parameters);
	nms->parameter_range = (double*)malloc(sizeof(double) * nms->number_parameters);
	nms->min_bound = (double*)malloc(sizeof(double) * nms->number_parameters);
	nms->max_bound = (double*)malloc(sizeof(double) * nms->number_parameters);

	memcpy(nms->current_point, point, sizeof(double) * nms->number_parameters);
	memcpy(nms->parameter_range, range, sizeof(double) * nms->number_parameters);
	memcpy(nms->initial_range, range, sizeof(double) * nms->number_parameters);
	memcpy(nms->min_bound, min_bound, sizeof(double) * nms->number_parameters);
	memcpy(nms->max_bound, max_bound, sizeof(double) * nms->number_parameters);

	for (i = 0; i < nms->number_parameters; i++) {
		if (type == NEWTON_LINE_SEARCH) nms->line_search->direction[i] = 0;
		if (nms->parameter_range[i] < 0) nms->parameter_range[i] *= -1.0;
	}

	new_population(nms->evaluations_per_iteration, nms->number_parameters, &(nms->population));

	return checkpoint_newton_method(search_name, nms);	
}

int write_newton_method(char* search_name, void* search_data) {
	char search_filename[FILENAME_SIZE], population_filename[FILENAME_SIZE];
	FILE *search_file;
	int result;
	NEWTON_METHOD_SEARCH *nms = (NEWTON_METHOD_SEARCH*)search_data;

	sprintf(search_filename, "%s/%s/search", get_working_directory(), search_name);
	search_file = fopen(search_filename, "w+");
	if (search_file == NULL) return -1;

	fprintf(search_file, "type: %d\n", nms->type);
	fprintf(search_file, "mode: %d\n", nms->mode);

	print_double_array(search_file, "current_point", nms->number_parameters, nms->current_point);
	print_double_array(search_file, "initial_range", nms->number_parameters, nms->parameter_range);
	print_double_array(search_file, "parameter_range", nms->number_parameters, nms->parameter_range);
	print_double_array(search_file, "min_bound", nms->number_parameters, nms->min_bound);
	print_double_array(search_file, "max_bound", nms->number_parameters, nms->max_bound);

	fprintf(search_file, "current_iteration: %d, maximum_iteration: %d\n", nms->current_iteration, nms->maximum_iteration);
	fprintf(search_file, "current_evaluation: %d, evaluations_per_iteration: %d\n", nms->current_evaluation, nms->evaluations_per_iteration);

	fprintf(search_file, "remove_outliers: %d\n", nms->remove_outliers);
	fprintf(search_file, "line_search: %d\n", (nms->line_search != NULL));
	if (nms->line_search != NULL) {
		fprintf(search_file, "center: %lf, min_range: %lf, max_range: %lf\n", nms->line_search->center, nms->line_search->min_range, nms->line_search->max_range);
		print_double_array(search_file, "direction", nms->number_parameters, nms->line_search->direction);
	}
	fclose(search_file);

	if (nms->line_search != NULL && nms->mode == NMS_LINE_SEARCH) {
		sprintf(population_filename, "%s/%s/population_%d_ls", get_working_directory(), search_name, nms->current_iteration);
		result = write_population(population_filename, nms->population);
		if (result < 0) return AS_CP_ERROR;
	} else {
		sprintf(population_filename, "%s/%s/population_%d", get_working_directory(), search_name, nms->current_iteration);
		result = write_population(population_filename, nms->population);
		if (result < 0) return AS_CP_ERROR;
	}

	return AS_CP_SUCCESS;
}

int bound_newton_method(NEWTON_METHOD_SEARCH *nms) {
	int i;
	for (i = 0; i < nms->number_parameters; i++) {
		if (isnan(nms->current_point[i]) || isnan(nms->parameter_range[i])) return 1;

		if (nms->current_point[i] > nms->max_bound[i]) nms->current_point[i] = nms->max_bound[i];
		if (nms->current_point[i] < nms->min_bound[i]) nms->current_point[i] = nms->min_bound[i];

		if (nms->parameter_range[i] < 0) nms->parameter_range[i] *= -1.0;
	}

	if (nms->line_search != NULL && nms->mode == NMS_LINE_SEARCH) {
		double min_bound, max_bound;
		if (nms->line_search->min_range < nms->line_search->max_range) {
			min_bound = nms->line_search->min_range;
			max_bound = nms->line_search->max_range;
		} else {
			min_bound = nms->line_search->max_range;
			max_bound = nms->line_search->min_range;
		}

		if (isnan(min_bound) || isnan(max_bound)) return 1;
		for (i = 0; i < nms->number_parameters; i++) {
			double temp_min = min_bound;
			double temp_max = max_bound;

			if (isnan(nms->line_search->direction[i])) return 1;

			if (nms->line_search->direction[i] > 0) {
				temp_max = (nms->max_bound[i] - nms->current_point[i]) / nms->line_search->direction[i];
				temp_min = (nms->min_bound[i] - nms->current_point[i]) / nms->line_search->direction[i];
			} else if (nms->line_search->direction[i] < 0) {
				temp_min = (nms->max_bound[i] - nms->current_point[i]) / nms->line_search->direction[i];
				temp_max = (nms->min_bound[i] - nms->current_point[i]) / nms->line_search->direction[i];
			}
			if (temp_min > min_bound) min_bound = temp_min;
			if (temp_max < max_bound) max_bound = temp_max;
			printf("min_bound: %.20lf, temp_min: %.20lf, max_bound: %.20lf, temp_max: %.20lf\n", min_bound, temp_min, max_bound, temp_max);
		}
		nms->line_search->min_range = min_bound * .98;
		nms->line_search->max_range = max_bound * .98;
		printf("new min bound: %.20lf, max bound: %.20lf\n", min_bound, max_bound);
	}
	return 0;
}


int read_newton_method(char* search_name, void** search_data) {
	char search_filename[FILENAME_SIZE], population_filename[FILENAME_SIZE];
	FILE *search_file;
	int result, line_search;
	NEWTON_METHOD_SEARCH **nms;
	nms = (NEWTON_METHOD_SEARCH**)search_data;
	(*nms) = (NEWTON_METHOD_SEARCH*)malloc(sizeof(NEWTON_METHOD_SEARCH));

	sprintf(search_filename, "%s/%s/search", get_working_directory(), search_name);
	search_file = fopen(search_filename, "r");
	if (search_file == NULL) return -1;

	fscanf(search_file, "type: %d\n", &((*nms)->type));
	fscanf(search_file, "mode: %d\n", &((*nms)->mode));

	(*nms)->number_parameters = read_double_array(search_file, "current_point", &((*nms)->current_point));
	read_double_array(search_file, "initial_range", &((*nms)->initial_range));
	read_double_array(search_file, "parameter_range", &((*nms)->parameter_range));
	read_double_array(search_file, "min_bound", &((*nms)->min_bound));
	read_double_array(search_file, "max_bound", &((*nms)->max_bound));

	fscanf(search_file, "current_iteration: %d, maximum_iteration: %d\n", &((*nms)->current_iteration), &((*nms)->maximum_iteration));
	fscanf(search_file, "current_evaluation: %d, evaluations_per_iteration: %d\n", &((*nms)->current_evaluation), &((*nms)->evaluations_per_iteration));

	fscanf(search_file, "remove_outliers: %d\n", &((*nms)->remove_outliers));
	fscanf(search_file, "line_search: %d\n", &line_search);
	if (line_search) {
		LINE_SEARCH *ls;

		printf("initializing line search\n");

		(*nms)->line_search = (LINE_SEARCH*)malloc(sizeof(LINE_SEARCH));
		ls = (*nms)->line_search;
		fscanf(search_file, "center: %lf, min_range: %lf, max_range: %lf\n", &(ls->center), &(ls->min_range), &(ls->max_range));
		read_double_array(search_file, "direction", &((*nms)->line_search->direction));
	} else {
		(*nms)->line_search = NULL;
	}
	fclose(search_file);

	printf("line search: %d, mode: %d\n", line_search, (*nms)->mode);

	if (line_search && (*nms)->mode == NMS_LINE_SEARCH) {
		printf("reading line search population\n");
		sprintf(population_filename, "%s/%s/population_%d_ls", get_working_directory(), search_name, (*nms)->current_iteration);
		result = read_population(population_filename, &((*nms)->population));
		if (result < 0) return result;
	} else {
		printf("reading population\n");
		sprintf(population_filename, "%s/%s/population_%d", get_working_directory(), search_name, (*nms)->current_iteration);
		result = read_population(population_filename, &((*nms)->population));
		if (result < 0) return result;
	}
	(*nms)->current_evaluation = (*nms)->population->size;

	if (bound_newton_method((*nms))) {
		(*nms)->current_iteration = (*nms)->maximum_iteration;
		write_newton_method(search_name, (*nms));
		sprintf(AS_MSG, "Search completed because of NaN");
		log_printf(search_name, "Search completed because of NaN");
		return AS_INSERT_OVER;
	}

	return 1;
}


int checkpoint_newton_method(char* search_name, void* search_data) {
	NEWTON_METHOD_SEARCH *nms = (NEWTON_METHOD_SEARCH*)search_data;

	sprintf(AS_MSG, "evaluation: %d/%d, iteration: %d/%d", nms->current_evaluation, nms->evaluations_per_iteration, nms->current_iteration, nms->maximum_iteration);
	if (nms->current_iteration >= nms->maximum_iteration) return AS_CP_OVER;
	else return write_newton_method(search_name, search_data);
}

int newton_bound_parameters(int number_parameters, double *parameters, double *min_bound, double *max_bound) {
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
	if (newton_bound_parameters(sp->number_parameters, sp->parameters, nms->min_bound, nms->max_bound)) return AS_GEN_FAIL;

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

double fmin(double f1, double f2) {
	if (f1 < f2) return f1;
	else return f2;
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

					printf("removing outliers for: %s\n", search_name);
					remove_outliers_incremental(p, 25.0);
					printf("removed outliers.\n");

					get_population_statistics(p, best_point, &best_fitness, &average_fitness, &worst_fitness, &standard_deviation);
					log_printf(search_name, "[after outliers ] best_fitness: %.20lf, average_fitness: %.20lf, worst_fitness: %.20lf, st_dev: %.20lf\n", best_fitness, average_fitness, worst_fitness, standard_deviation);

					step = (double*)malloc(sizeof(double) * nms->number_parameters);
					step_error = (double*)malloc(sizeof(double) * nms->number_parameters);

					printf("getting newton step\n");
					get_newton_step(nms, step, step_error);
					printf("got newton step\n");
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

					remove_outliers_incremental(p, 25.0);

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
