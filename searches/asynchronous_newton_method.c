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
	as->read_search = read_newton_method;
	as->checkpoint_search = checkpoint_newton_method;
	as->generate_parameters = newton_generate_parameters;
	as->insert_parameters = newton_insert_parameters;
	return as;
}

int create_newton_method(char* search_name, int type, int line_search, int remove_outliers, int maximum_iteration, int evaluations_per_iteration, int number_parameters, double *point, double *range, double *min_bound, double *max_bound) {
	char search_directory[FILENAME_SIZE];
	NEWTON_METHOD_SEARCH *nms;
	int i = 0;

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

	if (line_search) {
		nms->line_search = (LINE_SEARCH*)malloc(sizeof(LINE_SEARCH));
		nms->line_search->iteration = -1;
		nms->line_search->max_iteration = 1;
		nms->line_search->evaluation = 0;
		nms->line_search->center = 0;
		nms->line_search->max_range = 0;
		nms->line_search->min_range = 0;
	} else {
		nms->line_search = NULL;
	}

	printf("nms->type: %d\nnms->maximum_iteration: %d\nnms->evaluations_per_iteration: %d\nnms->number_parameters: %d\n", nms->type, nms->maximum_iteration, nms->evaluations_per_iteration, nms->number_parameters);

	nms->current_point = (double*)malloc(sizeof(double) * nms->number_parameters);
	nms->previous_point = (double*)malloc(sizeof(double) * nms->number_parameters);
	nms->direction = (double*)malloc(sizeof(double) * nms->number_parameters);
	nms->parameter_range = (double*)malloc(sizeof(double) * nms->number_parameters);
	nms->min_bound = (double*)malloc(sizeof(double) * nms->number_parameters);
	nms->max_bound = (double*)malloc(sizeof(double) * nms->number_parameters);

	memcpy(nms->current_point, point, sizeof(double) * nms->number_parameters);
	memcpy(nms->parameter_range, range, sizeof(double) * nms->number_parameters);
	memcpy(nms->min_bound, min_bound, sizeof(double) * nms->number_parameters);
	memcpy(nms->max_bound, max_bound, sizeof(double) * nms->number_parameters);

	for (i = 0; i < nms->number_parameters; i++) {
		nms->previous_point[i] = nms->current_point[i];
		nms->direction[i] = 0;
		if (nms->parameter_range[i] < 0) nms->parameter_range[i] *= -1.0;
	}

	new_population(nms->evaluations_per_iteration, nms->number_parameters, &(nms->population));

	return checkpoint_newton_method(search_name, nms);	
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
	read_double_array(search_file, "previous_point", &((*nms)->previous_point));
	read_double_array(search_file, "direction", &((*nms)->direction));
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
		fscanf(search_file, "iteration: %d, max_iteration: %d\n", &(ls->iteration), &(ls->max_iteration));
		fscanf(search_file, "evaluation: %d\n", &(ls->evaluation));
		fscanf(search_file, "center: %lf, min_range: %lf, max_range: %lf\n", &(ls->center), &(ls->min_range), &(ls->max_range));
	} else {
		(*nms)->line_search = NULL;
	}
	fclose(search_file);

	printf("line search: %d, mode: %d\n", line_search, (*nms)->mode);

	if (line_search && (*nms)->mode == NMS_LINE_SEARCH) {
		printf("reading line search population\n");
		sprintf(population_filename, "%s/%s/population_%d_ls_%d", get_working_directory(), search_name, (*nms)->current_iteration, (*nms)->line_search->iteration);
		result = read_population(population_filename, &((*nms)->population));
		if (result < 0) return result;
	} else {
		printf("reading population\n");
		sprintf(population_filename, "%s/%s/population_%d", get_working_directory(), search_name, (*nms)->current_iteration);
		result = read_population(population_filename, &((*nms)->population));
		if (result < 0) return result;
	}
	(*nms)->current_evaluation = (*nms)->population->size;

	return 1;
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
	print_double_array(search_file, "previous_point", nms->number_parameters, nms->previous_point);
	print_double_array(search_file, "direction", nms->number_parameters, nms->direction);
	print_double_array(search_file, "parameter_range", nms->number_parameters, nms->parameter_range);
	print_double_array(search_file, "min_bound", nms->number_parameters, nms->min_bound);
	print_double_array(search_file, "max_bound", nms->number_parameters, nms->max_bound);

	fprintf(search_file, "current_iteration: %d, maximum_iteration: %d\n", nms->current_iteration, nms->maximum_iteration);
	fprintf(search_file, "current_evaluation: %d, evaluations_per_iteration: %d\n", nms->current_evaluation, nms->evaluations_per_iteration);

	fprintf(search_file, "remove_outliers: %d\n", nms->remove_outliers);
	fprintf(search_file, "line_search: %d\n", (nms->line_search != NULL));
	if (nms->line_search != NULL) {
		fprintf(search_file, "iteration: %d, max_iteration: %d\n", nms->line_search->iteration, nms->line_search->max_iteration);
		fprintf(search_file, "evaluation: %d\n", nms->line_search->evaluation);
		fprintf(search_file, "center: %lf, min_range: %lf, max_range: %lf\n", nms->line_search->center, nms->line_search->min_range, nms->line_search->max_range);
	}
	fclose(search_file);

	if (nms->line_search != NULL && nms->mode == NMS_LINE_SEARCH) {
		sprintf(population_filename, "%s/%s/population_%d_ls_%d", get_working_directory(), search_name, nms->current_iteration, nms->line_search->iteration);
		result = write_population(population_filename, nms->population);
		if (result < 0) return AS_CP_ERROR;
	} else {
		sprintf(population_filename, "%s/%s/population_%d", get_working_directory(), search_name, nms->current_iteration);
		result = write_population(population_filename, nms->population);
		if (result < 0) return AS_CP_ERROR;
	}

	return AS_CP_SUCCESS;
}


int checkpoint_newton_method(char* search_name, void* search_data) {
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

	if (nms->line_search != NULL && nms->mode == NMS_LINE_SEARCH) {
		double point; 
		point = random_linear_recombination(p->number_parameters, nms->line_search->min_range, nms->line_search->max_range, nms->current_point, nms->direction, sp->parameters);
		sprintf(sp->metadata, "ls, point: %.20lf, it: %d, ls_it: %d, ev: %d\n", point, nms->current_iteration, nms->line_search->iteration, nms->line_search->evaluation);
	} else {
		sprintf(sp->metadata, "it: %d, ev: %d", nms->current_iteration, nms->current_evaluation);
		range_recombination(p->number_parameters, nms->current_point, nms->parameter_range, sp->parameters);
	}
	if (bound_parameters(sp->number_parameters, sp->parameters, nms->min_bound, nms->max_bound)) return AS_GEN_FAIL;

	return AS_GEN_SUCCESS;
}

int bound_newton_method(NEWTON_METHOD_SEARCH *nms) {
	int j;
	for (j = 0; j < nms->number_parameters; j++) {
		if (isnan(nms->current_point[j]) || isnan(nms->parameter_range[j])) return 1;

		if (nms->current_point[j] > nms->max_bound[j]) nms->current_point[j] = nms->max_bound[j];
		if (nms->current_point[j] < nms->min_bound[j]) nms->current_point[j] = nms->min_bound[j];

		if (nms->parameter_range[j] < 0) nms->parameter_range[j] *= -1.0;
	}

	if (nms->mode == NMS_LINE_SEARCH) {
		if (nms->line_search->min_range < nms->line_search->max_range) {
			double temp = nms->line_search->min_range;
			nms->line_search->min_range = nms->line_search->max_range;
			nms->line_search->max_range = temp;
		}
	}

	return 0;
}

int verify(NEWTON_METHOD_SEARCH* nms, SEARCH_PARAMETERS* sp) {
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

	if (nms->mode == NMS_LINE_SEARCH) {
		int current_iteration, ls_iteration, ls_evaluation;
		double point;
		LINE_SEARCH* ls = nms->line_search;

		if (4 != sscanf(sp->metadata, "ls, point: %lf, it: %d, ls_it: %d, ev: %d\n", &point, &current_iteration, &ls_iteration, &ls_evaluation)) return AS_INSERT_BAD_METADATA;

		if (current_iteration != nms->current_iteration || ls_iteration != ls->iteration) return AS_INSERT_OUT_OF_ITERATION;
		sp->number_parameters = 1;
		free(sp->parameters);
		sp->parameters = (double*)malloc(sizeof(double));
		sp->parameters[0] = point;
	} else {
		for (i = 0; i < nms->number_parameters; i++) {
			if (sp->parameters[i] < (nms->current_point[i] - nms->parameter_range[i]) || sp->parameters[i] > (nms->current_point[i] + nms->parameter_range[i])) return AS_INSERT_OUT_OF_RANGE;
		}
	}

	if (population_contains(nms->population, sp->fitness, sp->parameters)) return AS_INSERT_NOT_UNIQUE;

	return 0;
}

void iterate_newton_method(NEWTON_METHOD_SEARCH *nms) {
	double **hessian, **hessian_error;
	double *gradient, *gradient_error;
	double c, c_error;
	POPULATION *p = nms->population; 
	int j, k;

	for (j = 0; j < p->size; j++) {
		for (k = 0; k < p->number_parameters; k++) p->individuals[j][k] -= nms->current_point[k];
	}

	printf("updated individuals\n");

	init_matrix(&hessian, nms->number_parameters, nms->number_parameters);
	init_matrix(&hessian_error, nms->number_parameters, nms->number_parameters);
	gradient = (double*)malloc(sizeof(double) * nms->number_parameters);
	gradient_error = (double*)malloc(sizeof(double) * nms->number_parameters);

	printf("doing regression\n");

	parabolic_2d_regression_error(p->size, p->number_parameters, p->individuals, p->fitness, hessian, hessian_error, gradient, gradient_error, &c, &c_error);

	printf("did regression\n");

	for (j = 0; j < p->size; j++) {
		for (k = 0; k < p->number_parameters; k++) p->individuals[j][k] += nms->current_point[k];
	}

	printf("doing newton step\n");

	newton_step(nms->number_parameters, hessian, gradient, nms->direction);

	printf("doing update/error\n");
	if (nms->type == NEWTON_ERROR_RANGE) {
		newton_step(nms->number_parameters, hessian_error, gradient_error, nms->parameter_range);
		for (j = 0; j < nms->number_parameters; j++) {
			nms->parameter_range[j] = fabs(nms->parameter_range[j]);
		}
	} else if (nms->type == NEWTON_UPDATE_RANGE) {
		for (j = 0; j < nms->number_parameters; j++) {
			nms->parameter_range[j] = fabs(nms->direction[j]);
		}
	}

	printf("did newton step\n");

	for (j = 0; j < nms->number_parameters; j++) {
		nms->previous_point[j] = nms->current_point[j];
		nms->current_point[j] -= nms->direction[j];
	}

	printf("freeing\n");

	free_matrix(&hessian, nms->number_parameters, nms->number_parameters);
	free_matrix(&hessian_error, nms->number_parameters, nms->number_parameters);
	free(gradient);
	free(gradient_error);

	printf("finished\n");
}

void iterate_line_search(NEWTON_METHOD_SEARCH *nms) {
	double a, b, c, a_error, b_error, c_error;
	double *individuals;
	int j;
	LINE_SEARCH *ls = nms->line_search;
	POPULATION *p = nms->population;

	individuals = (double*)malloc(sizeof(double) * p->size);
	for (j = 0; j < p->size; j++) individuals[j] = p->individuals[j][0];

	parabolic_regression(p->number_parameters, individuals, p->fitness, &a, &b, &c, &a_error, &b_error, &c_error);
	ls->center = parabolic_center(a, b, c);
	ls->min_range = parabolic_center(a-a_error, b-b_error, c-c_error);
	ls->max_range = parabolic_center(a+a_error, b+b_error, c+c_error);

	free(individuals);
}

int newton_insert_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS* sp) {
	NEWTON_METHOD_SEARCH *nms = (NEWTON_METHOD_SEARCH*)search_data;
	POPULATION *p = nms->population;
	int result, j, k;

	/********
		*	Insert parameters into population.  If cutoff reached, calculate hessian
		*	and generate new population.
	 ********/
	
	if (nms->line_search == NULL) {
		sprintf(AS_MSG, "ev: %*d/%*d, it: %*d/%*d, f: %.15lf, m: %s", 3, nms->current_evaluation, 3, nms->evaluations_per_iteration, 3, nms->current_iteration, 3, nms->maximum_iteration, sp->fitness, sp->metadata);
	} else {
		sprintf(AS_MSG, "ev: %*d/%*d, it: %*d/%*d, f: %.15lf, ls: %d/%d, m: %s", 3, nms->current_evaluation, 3, nms->evaluations_per_iteration, 3, nms->current_iteration, 3, nms->maximum_iteration, sp->fitness, nms->line_search->iteration, nms->line_search->max_iteration, sp->metadata);
	}

	if (nms->current_iteration < nms->maximum_iteration) {
		result = verify(nms, sp);
		if (result != 0) return result;

		insert_incremental_info(p, sp->parameters, sp->fitness, sp->host_os, sp->app_version);
//		insert_incremental(p, sp->parameters, sp->fitness);
		nms->current_evaluation++;

		if (nms->current_evaluation >= nms->evaluations_per_iteration) {
			nms->current_evaluation = 0;

			if (nms->current_iteration < nms->maximum_iteration) {
				/********
					*	Remove any outliers.
				 ********/
				double *best_point, best_fitness, average_fitness, worst_fitness, standard_deviation;
				best_point = (double*)malloc(sizeof(double) * nms->number_parameters);

				write_newton_method(search_name, nms);

				get_population_statistics(p, best_point, &best_fitness, &average_fitness, &worst_fitness, &standard_deviation);
				log_printf(search_name, "best_fitness: %.20lf, average_fitness: %.20lf, worst_fitness: %.20lf, st_dev: %.20lf\n", best_fitness, average_fitness, worst_fitness, standard_deviation);

				if (nms->remove_outliers > 0) {
					for (j = 0; j < p->size; j++) {
						for (k = 0; k < p->number_parameters; k++) p->individuals[j][k] -= nms->current_point[k];
					}
					remove_outliers(p, 2.0);
					for (j = 0; j < p->size; j++) {
						for (k = 0; k < p->number_parameters; k++) p->individuals[j][k] += nms->current_point[k];
					}
					get_population_statistics(p, best_point, &best_fitness, &average_fitness, &worst_fitness, &standard_deviation);
					log_printf(search_name, "OUTLIERS REMOVED: best_fitness: %.20lf, average_fitness: %.20lf, worst_fitness: %.20lf, st_dev: %.20lf\n", best_fitness, average_fitness, worst_fitness, standard_deviation);
				}
				free(best_point);

				if (nms->line_search != NULL) {
					if (nms->mode == NMS_LINE_SEARCH) {
						iterate_line_search(nms);
						nms->line_search->iteration++;
						if (nms->line_search->iteration >= nms->line_search->max_iteration) {
							for (j = 0; j < nms->number_parameters; j++) {
								nms->previous_point[j] = nms->current_point[j];
								nms->current_point[j] -= nms->line_search->center * nms->direction[j];
							}
							nms->line_search->iteration = -1;
							nms->mode = NMS_HESSIAN;
							nms->current_iteration++;

							free_population(nms->population);
							free(nms->population);
							new_population(nms->evaluations_per_iteration, nms->number_parameters, &(nms->population));
						} else {
							free_population(nms->population);
							free(nms->population);
							new_population(nms->evaluations_per_iteration, 1, &(nms->population));
						}
					} else {
						iterate_newton_method(nms);

						nms->line_search->iteration = 0;
						nms->line_search->center = 0;
						nms->line_search->min_range = -1;
						nms->line_search->max_range = 2;
						nms->mode = NMS_LINE_SEARCH;

						free_population(nms->population);
						free(nms->population);
						new_population(nms->evaluations_per_iteration, 1, &(nms->population));
					}
				} else {
					iterate_newton_method(nms);
					nms->current_iteration++;

					free_population(nms->population);
					free(nms->population);
					new_population(nms->evaluations_per_iteration, nms->number_parameters, &(nms->population));
				}
				write_newton_method(search_name, nms);

				if (bound_newton_method(nms)) {
					nms->current_iteration = nms->maximum_iteration;
					write_newton_method(search_name, nms);
					sprintf(AS_MSG, "Search completed because of NaN");
					return AS_INSERT_OVER;
				}

				log_print_double_array(search_name, "current_point", nms->number_parameters, nms->current_point);
				log_print_double_array(search_name, "direction", nms->number_parameters, nms->direction);
				log_printf(search_name, "\n");
			}
		}
	} else {
		return AS_INSERT_OVER;
	}
	return AS_INSERT_SUCCESS;
}
