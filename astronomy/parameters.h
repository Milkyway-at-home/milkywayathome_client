#ifndef ASTRONOMY_PARAMETERS_H
#define ASTRONOMY_PARAMETERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct astronomy_parameters {
	int number_parameters;

	int number_background_parameters;
	double background_weight;
	double *background_parameters;
	double *background_step;
	double *background_min;
	double *background_max;
	int *background_optimize;

	int number_streams;
	int number_stream_parameters;

	double *stream_weights;
	double *stream_weight_step;
	double *stream_weight_min;
	double *stream_weight_max;
	int *stream_weight_optimize;

	double **stream_parameters;
	double **stream_step;
	double **stream_min;
	double **stream_max;
	int **stream_optimize;

	int convolve;
	int sgr_coordinates;
	int wedge;
	int r_steps;
	int mu_steps;
	int nu_steps;
	double r_min, r_max, r_step_size;
	double mu_min, mu_max, mu_step_size;
	double nu_min, nu_max, nu_step_size;
	
	int r_cut_steps;
	int mu_cut_steps;
	int nu_cut_steps;
	double r_cut_min, r_cut_max, r_cut_step_size;
	double mu_cut_min, mu_cut_max, mu_cut_step_size;
	double nu_cut_min, nu_cut_max, nu_cut_step_size;

} ASTRONOMY_PARAMETERS;

int	read_astronomy_parameters(const char* file, ASTRONOMY_PARAMETERS *ap);
void	fread_astronomy_parameters(FILE* file, ASTRONOMY_PARAMETERS *ap);
int	write_astronomy_parameters(const char* file, ASTRONOMY_PARAMETERS *ap);
void	fwrite_astronomy_parameters(FILE* file, ASTRONOMY_PARAMETERS *ap);
void	free_parameters(ASTRONOMY_PARAMETERS* ap);

void	split_astronomy_parameters(ASTRONOMY_PARAMETERS* ap, int rank, int max_rank);

void	set_astronomy_parameters(ASTRONOMY_PARAMETERS* ap, double* parameters);
void	get_search_parameters(ASTRONOMY_PARAMETERS* ap, double** parameters);
void	get_min_parameters(ASTRONOMY_PARAMETERS* ap, double** parameters);
void	get_max_parameters(ASTRONOMY_PARAMETERS* ap, double** parameters);
void	get_step(ASTRONOMY_PARAMETERS* ap, double** step);

#ifdef GMLE_BOINC
	int	boinc_read_astronomy_parameters(const char* file, ASTRONOMY_PARAMETERS *ap);
	int	boinc_write_astronomy_parameters(const char* file, ASTRONOMY_PARAMETERS *ap);
#endif
#endif
