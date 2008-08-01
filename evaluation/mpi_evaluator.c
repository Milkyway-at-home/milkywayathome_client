#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

#include "mpi_evaluator.h"
#include "evaluator.h"

char    	hostname[MPI_MAX_PROCESSOR_NAME];
int     	hostname_length;
int     	rank, number_processes;

int		integral_defined = 0;
double*		(*integral_function)(double*) = NULL;
double*		(*integral_combinator)(double*, int) = NULL;
int		integral_parameter_length;
int		integral_results_length;

double*		(*likelihood_function)(double*) = NULL;
double		(*likelihood_combinator)(double*, int) = NULL;
int		likelihood_parameter_length;
int		likelihood_results_length;

double mpi_evaluate(double* likelihood_parameters) {
	double* likelihood_results_send;
	double* likelihood_results_recv;
	double result;

	likelihood_results_recv = (double*)malloc(sizeof(double) * likelihood_results_length);
	if (rank == 0) {

		MPI_Bcast(likelihood_parameters, likelihood_parameter_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		likelihood_results_send = (*likelihood_function)(likelihood_parameters);
		MPI_Gather(likelihood_results_send, likelihood_results_length, MPI_DOUBLE, likelihood_results_recv, likelihood_results_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		result = (*likelihood_combinator)(likelihood_results_recv, number_processes);

		free(likelihood_results_send);
		free(likelihood_results_recv);
		return result;
	} else {
		fprintf(stderr, "ERROR: calling evaluate from non-master process.");
		return -1.0;
	}
}

double mpi_integral_evaluate(double* integral_parameters) {
	double* integral_results_send;
	double* integral_results_recv;
	double* likelihood_parameters;
	double* likelihood_results_send;
	double* likelihood_results_recv;
	double result;
	int i;

	likelihood_parameters = (double*)malloc(sizeof(double) * likelihood_parameter_length);
	likelihood_results_recv = (double*)malloc(sizeof(double) * likelihood_results_length);
	integral_results_recv = (double*)malloc(sizeof(double) * integral_results_length);
	if (rank == 0) {
		MPI_Bcast(integral_parameters, integral_parameter_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		integral_results_send = (*integral_function)(integral_parameters);
		MPI_Gather(integral_results_send, integral_results_length, MPI_DOUBLE, integral_results_recv, integral_results_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		likelihood_parameters = (*integral_combinator)(integral_results_recv, number_processes);

		MPI_Bcast(likelihood_parameters, likelihood_parameter_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		likelihood_results_send = (*likelihood_function)(likelihood_parameters);
		MPI_Gather(likelihood_results_send, likelihood_results_length, MPI_DOUBLE, likelihood_results_recv, likelihood_results_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		result = (*likelihood_combinator)(likelihood_results_recv, number_processes);

		free(integral_results_send);
		free(integral_results_recv);
		free(likelihood_parameters);
		free(likelihood_results_send);
		free(likelihood_results_recv);
		return result;
	} else {
		fprintf(stderr, "ERROR: calling evaluate from non-master process.");
		return -1.0;
	}
}

void evaluator__init_likelihood(double* (*l_f)(double*), int l_p_l, double (*l_c)(double*, int), int l_r_l) {
	likelihood_function = (*l_f);
	likelihood_combinator = (*l_c);
	likelihood_parameter_length = l_p_l;
	likelihood_results_length = l_r_l;

	if (integral_defined == 0) evaluate = mpi_evaluate;
}

void evaluator__init_integral(double* (*i_f)(double*), int i_p_l, double* (*i_c)(double*, int), int i_r_l) {
	integral_function = (*i_f);
	integral_combinator = (*i_c);
	integral_parameter_length = i_p_l;
	integral_results_length = i_r_l;

	integral_defined = 1;
	evaluate = mpi_integral_evaluate;
}

void evaluator__init_data(void (*read_data)(int, int)) {
        read_data(rank, number_processes);
}

void mpi_evaluator__start(int number_arguments, char** arguments) {
	int completed, evaluations;
	double* likelihood_parameters;
	double* likelihood_results_send;
	double* likelihood_results_recv;
	double* integral_parameters;
	double* integral_results_send;
	double* integral_results_recv;

	MPI_Init(&number_arguments, &arguments);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
	MPI_Get_processor_name(hostname, &hostname_length);

	evaluations = 0;
	completed = 0;
	if (rank != 0) {
		//malloc likelihood and integral arrays
		likelihood_parameters = (double*)malloc(sizeof(double) * likelihood_parameter_length);
		likelihood_results_send = (double*)malloc(sizeof(double) * likelihood_results_length);
		likelihood_results_recv = (double*)malloc(sizeof(double) * likelihood_results_length);
		integral_parameters = (double*)malloc(sizeof(double) * integral_parameter_length);
		integral_results_send = (double*)malloc(sizeof(double) * integral_results_length);
		integral_results_recv = (double*)malloc(sizeof(double) * integral_results_length);

		while (!completed) {
			if (integral_defined > 0) {
				MPI_Bcast(integral_parameters, integral_parameter_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				integral_results_send = (*integral_function)(integral_parameters);
				MPI_Gather(integral_results_send, integral_results_length, MPI_DOUBLE, integral_results_recv, integral_results_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
			MPI_Bcast(likelihood_parameters, likelihood_parameter_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			likelihood_results_send = (*likelihood_function)(likelihood_parameters);
			MPI_Gather(likelihood_results_send, likelihood_results_length, MPI_DOUBLE, likelihood_results_recv, likelihood_results_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			evaluations++;
		}
	}
}
