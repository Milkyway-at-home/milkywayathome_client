#ifndef ASTRONOMY_EVALUATION_STATE_H
#define ASTRONOMY_EVALUATION_STATE_H

#include "parameters.h"
#include "star_points.h"

typedef struct integral_area {
	int mu_steps, nu_steps, r_steps;
	double mu_min, nu_min, r_min;
	double mu_max, nu_max, r_max;
	double mu_step_size, nu_step_size, r_step_size;
	int mu_step_current, nu_step_current, r_step_current;	

	int number_streams;
	double background_integral, *stream_integrals;
} INTEGRAL_AREA;

typedef struct evaluation_state {
	/********
		*	State for integral calculation.
	 ********/
	INTEGRAL_AREA *main_integral, **cuts;
	int current_cut, number_streams, number_cuts;

	double background_integral;
	double* stream_integrals;

	/********
		*	State for likelihood calculation.
	 ********/
	int current_star_point, total_stars;
	int num_zero;
	int bad_jacobians;
	double prob_sum;
} EVALUATION_STATE;

void	initialize_state(ASTRONOMY_PARAMETERS *ap, STAR_POINTS *sp, EVALUATION_STATE* es);
void	free_state(EVALUATION_STATE* es);
void	reset_evaluation_state(EVALUATION_STATE *es);

int	write_checkpoint(EVALUATION_STATE* es);
int	read_checkpoint(EVALUATION_STATE* es);

#endif
