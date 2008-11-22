#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "asynchronous_search.h"
#include "asynchronous_newton_method.h"
#include "gradient.h"
#include "hessian.h"
#include "population.h"
#include "recombination.h"
#include "search_parameters.h"

#include "../evaluation/search_manager.h"
#include "../util/settings.h"
#include "../util/matrix.h"
#include "../util/io_util.h"

const ASYNCHRONOUS_SEARCH asynchronous_newton_method = { "nm", create_newton_method, read_newton_method, checkpoint_newton_method, newton_generate_parameters, newton_insert_parameters };

/********
	*	Parameters: int maximum_iteration, int evaluations_per_iteration, int number_parameters, double* point, double* range
 ********/
int create_newton_method(char* search_name, ...) {
	char search_directory[FILENAME_SIZE];
	NEWTON_METHOD_SEARCH *nms;
	va_list vl;
	double *point, *range;

	sprintf(search_directory, "%s/%s", get_working_directory(), search_name);
	printf("making directory: %s\n", search_directory);
	mkdir(search_directory, 0777);

	nms = (NEWTON_METHOD_SEARCH*)malloc(sizeof(NEWTON_METHOD_SEARCH));
	va_start(vl, search_name);
	nms->current_iteration = 0;
	nms->maximum_iteration = va_arg(vl, int);
	nms->current_evaluation = 0;
	nms->evaluations_per_iteration = va_arg(vl, int);
	nms->number_parameters = va_arg(vl, int);
	point = va_arg(vl, double*);
	range = va_arg(vl, double*);

	nms->parameters = (double*)malloc(sizeof(double) * nms->number_parameters);
	nms->parameter_range = (double*)malloc(sizeof(double) * nms->number_parameters);
	memcpy(nms->parameters, point, sizeof(double) * nms->number_parameters);
	memcpy(nms->parameter_range, range, sizeof(double) * nms->number_parameters);

	new_population(nms->evaluations_per_iteration, nms->number_parameters, &(nms->population));

	return checkpoint_newton_method(search_name, nms);	
}

int read_newton_method(char* search_name, void** search_data) {
	char search_filename[FILENAME_SIZE], population_filename[FILENAME_SIZE];
	FILE *search_file;
	int result, np, i;
	NEWTON_METHOD_SEARCH **nms;
	nms = (NEWTON_METHOD_SEARCH**)search_data;
	(*nms) = (NEWTON_METHOD_SEARCH*)malloc(sizeof(NEWTON_METHOD_SEARCH));

	sprintf(search_filename, "%s/%s/search", get_working_directory(), search_name);
	search_file = fopen(search_filename, "r");
	if (search_file == NULL) return -1;
	(*nms)->number_parameters = read_double_array(search_file, "parameters", &((*nms)->parameters));
	read_double_array(search_file, "parameter_range", &((*nms)->parameter_range));
	fscanf(search_file, "current_iteration: %d, maximum_iteration: %d\n", &((*nms)->current_iteration), &((*nms)->maximum_iteration));
	fscanf(search_file, "current_evaluation: %d, evaluations_per_iteration: %d\n", &((*nms)->current_evaluation), &((*nms)->evaluations_per_iteration));

	np = (*nms)->number_parameters;
	(*nms)->min_parameters = (double*)malloc(sizeof(double) * np);
	(*nms)->max_parameters = (double*)malloc(sizeof(double) * np);
	for (i = 0; i < np; i++) {
		(*nms)->min_parameters[i] = (*nms)->parameters[i] - (*nms)->parameter_range[i];
		(*nms)->max_parameters[i] = (*nms)->parameters[i] + (*nms)->parameter_range[i];
	}
	fclose(search_file);

	sprintf(population_filename, "%s/%s/population_%d", get_working_directory(), search_name, (*nms)->current_iteration);
	result = read_population(population_filename, &((*nms)->population));
	if (result < 0) return result;
	(*nms)->current_evaluation = (*nms)->population->size;

	return 1;
}

int checkpoint_newton_method(char* search_name, void* search_data) {
	char search_filename[FILENAME_SIZE], population_filename[FILENAME_SIZE];
	FILE *search_file;
	int result;
	NEWTON_METHOD_SEARCH *nms = (NEWTON_METHOD_SEARCH*)search_data;

	sprintf(search_filename, "%s/%s/search", get_working_directory(), search_name);
	search_file = fopen(search_filename, "w+");
	if (search_file == NULL) return -1;
	print_double_array(search_file, "parameters", nms->number_parameters, nms->parameters);
	print_double_array(search_file, "parameter_range", nms->number_parameters, nms->parameter_range);
	fprintf(search_file, "current_iteration: %d, maximum_iteration: %d\n", nms->current_iteration, nms->maximum_iteration);
	fprintf(search_file, "current_evaluation: %d, evaluations_per_iteration: %d\n", nms->current_evaluation, nms->evaluations_per_iteration);
	fclose(search_file);

	sprintf(population_filename, "%s/%s/population_%d", get_working_directory(), search_name, nms->current_iteration);
	result = write_population(population_filename, nms->population);
	if (result < 0) return result;

	return 1;
}


int newton_generate_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS** sp) {
	POPULATION *p;
	NEWTON_METHOD_SEARCH *nms = (NEWTON_METHOD_SEARCH*)search_data;

	nms = (NEWTON_METHOD_SEARCH*)(search_data);
	p = nms->population;

	if (nms->current_iteration < nms->maximum_iteration) {
		char metadata[METADATA_SIZE];
		sprintf(metadata, "iteration: %d, evaluation: %d", nms->current_iteration, nms->current_evaluation);
		new_search_parameters(sp, search_name, p->number_parameters, random_recombination(nms->min_parameters, nms->max_parameters, p->number_parameters), metadata);
	}
	return 1;
}

int newton_insert_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS* sp) {
	NEWTON_METHOD_SEARCH *nms = (NEWTON_METHOD_SEARCH*)search_data;
	POPULATION *p = nms->population;

	/********
		*	Insert parameters into population.  If cutoff reached, calculate hessian
		*	and generate new population.
	 ********/
	
	if (nms->current_iteration < nms->maximum_iteration) {
		replace(p, nms->current_evaluation, sp->parameters, sp->fitness);
		nms->current_evaluation++;
		if (nms->current_evaluation >= nms->evaluations_per_iteration) {
			nms->current_evaluation = 0;
			nms->current_iteration++;
			if (nms->current_iteration < nms->maximum_iteration) {
				double **hessian, **inverse_hessian;
				double *gradient;
				int j, k;
				char filename[FILENAME_SIZE];

				randomized_hessian(p->individuals, nms->parameters, p->fitness, p->size, p->number_parameters, &hessian, &gradient);

				printf("\n");
				matrix_print(stdout, "hessian", hessian, p->number_parameters, p->number_parameters);
				printf("\n");
				print_double_array(stdout, "gradient: ", p->number_parameters, gradient);
				printf("\n");

				/********
					*	Take the newton step:  y = -hessian^-1 * gradient
				 ********/
				matrix_invert(hessian, p->number_parameters, p->number_parameters, &inverse_hessian);
				for (j = 0; j < p->number_parameters; j++) {
					nms->parameters[j] = 0;
					for (k = 0; k < p->number_parameters; k++) nms->parameters[j] -= inverse_hessian[j][k] * gradient[j];
					nms->min_parameters[j] = nms->parameters[j] - nms->parameter_range[j];
					nms->max_parameters[j] = nms->parameters[j] + nms->parameter_range[j];
				}
				print_double_array(stdout, "current point", nms->number_parameters, nms->parameters);
				print_double_array(stdout, "min range", nms->number_parameters, nms->min_parameters);
				print_double_array(stdout, "max range", nms->number_parameters, nms->max_parameters);

				for (j = 0; j < p->number_parameters; j++) {
					free(inverse_hessian[j]);
					free(hessian[j]);
				}
				free(hessian);
				free(gradient);

				/********
					*	TODO
					*	Optional: use previous parameters with new point as a line search direction
				 ********/
				sprintf(filename, "%s/%s/population_%d", get_working_directory(), search_name, nms->current_iteration-1);
				write_population(filename, p);
				free_population(p);
				free(p);

				new_population(nms->evaluations_per_iteration, nms->number_parameters, &(nms->population));
				sprintf(filename, "%s/%s/population_%d", get_working_directory(), search_name, nms->current_iteration);
				write_population(filename, p);
			}
		}
	}
	return 1;
}
