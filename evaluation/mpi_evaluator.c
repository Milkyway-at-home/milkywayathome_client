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

#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

#include "mpi_evaluator.h"
#include "evaluator.h"

#include "../util/io_util.h"

char    	hostname[MPI_MAX_PROCESSOR_NAME];
int     	hostname_length;
int     	rank, number_processes;

int		integral_defined = 0;

void		(*__integral_function)(double*, double**) = NULL;
void		(*__integral_combinator)(double*, int, double**) = NULL;
int		integral_parameter_length, integral_results_length;

void		(*__likelihood_function)(double*, double**) = NULL;
double		(*__likelihood_combinator)(double*, int) = NULL;
int		likelihood_parameter_length, likelihood_results_length;

double (*evaluate)(double*);

double* 	likelihood_parameters;
double* 	likelihood_results_send;
double* 	likelihood_results_recv;
double* 	integral_parameters;
double* 	integral_results_send;
double*		integral_results_recv;

int get_mpi_rank() {
	return rank;
}

double mpi_evaluate(double* likelihood_parameters) {
	double result = -1.0;

	if (rank == 0) {
		MPI_Bcast(likelihood_parameters, likelihood_parameter_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		__likelihood_function(likelihood_parameters, &likelihood_results_send);
		MPI_Gather(likelihood_results_send, likelihood_results_length, MPI_DOUBLE, likelihood_results_recv, likelihood_results_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		result = __likelihood_combinator(likelihood_results_recv, number_processes);
	} else {
		fprintf(stderr, "ERROR: calling evaluate from non-master process.");
	}
	return result;
}

double mpi_integral_evaluate(double* integral_parameters) {
	double result;

	if (rank == 0) {
		MPI_Bcast(integral_parameters, integral_parameter_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		__integral_function(integral_parameters, &integral_results_send);

		MPI_Gather(integral_results_send, integral_results_length, MPI_DOUBLE, integral_results_recv, integral_results_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		__integral_combinator(integral_results_recv, number_processes, &likelihood_parameters);

		MPI_Bcast(likelihood_parameters, likelihood_parameter_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		__likelihood_function(likelihood_parameters, &likelihood_results_send);

		MPI_Gather(likelihood_results_send, likelihood_results_length, MPI_DOUBLE, likelihood_results_recv, likelihood_results_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		result = __likelihood_combinator(likelihood_results_recv, number_processes);

		free(integral_results_send);
		free(likelihood_results_send);

		return result;
	} else {
		fprintf(stderr, "ERROR: calling evaluate from non-master process.");
		return -1.0;
	}
}

void evaluator__init_likelihood(void (*l_f)(double*, double**), int l_p_l, double (*l_c)(double*, int), int l_r_l) {
	__likelihood_function = (*l_f);
	__likelihood_combinator = (*l_c);
	likelihood_parameter_length = l_p_l;
	likelihood_results_length = l_r_l;

	if (evaluate == NULL) evaluate = mpi_evaluate;
	likelihood_results_recv = (double*)malloc(sizeof(double) * likelihood_results_length * number_processes);
}

void evaluator__init_integral(void (*i_f)(double*, double**), int i_p_l, void (*i_c)(double*, int, double**), int i_r_l) {
	__integral_function = (*i_f);
	__integral_combinator = (*i_c);
	integral_parameter_length = i_p_l;
	integral_results_length = i_r_l;

	integral_defined = 1;
	evaluate = mpi_integral_evaluate;
	integral_results_recv = (double*)malloc(sizeof(double) * integral_results_length * number_processes);
}

void mpi_evaluator__init(int *number_arguments, char*** arguments) {
	MPI_Init(number_arguments, arguments);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
	MPI_Get_processor_name(hostname, &hostname_length);
}

void mpi_evaluator__read_data(void (*r_d)(int, int)) {
	r_d(rank, number_processes);
}

void mpi_evaluator__start() {
	int completed, evaluations;

	evaluations = 0;
	completed = 0;
	if (rank != 0) {
		//malloc likelihood and integral arrays
		likelihood_parameters = (double*)malloc(sizeof(double) * likelihood_parameter_length);
		integral_parameters = (double*)malloc(sizeof(double) * integral_parameter_length);

		while (!completed) {
			if (integral_defined > 0) {
				MPI_Bcast(integral_parameters, integral_parameter_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				__integral_function(integral_parameters, &integral_results_send);
				MPI_Gather(integral_results_send, integral_results_length, MPI_DOUBLE, integral_results_recv, integral_results_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
			MPI_Bcast(likelihood_parameters, likelihood_parameter_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			__likelihood_function(likelihood_parameters, &likelihood_results_send);
			MPI_Gather(likelihood_results_send, likelihood_results_length, MPI_DOUBLE, likelihood_results_recv, likelihood_results_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			evaluations++;

			free(integral_results_send);
			free(likelihood_results_send);
		}
	}
}
