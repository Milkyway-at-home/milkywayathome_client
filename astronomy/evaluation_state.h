#ifndef ASTRONOMY_EVALUATION_STATE_H
#define ASTRONOMY_EVALUATION_STATE_H

#include "parameters.h"

typedef struct integral_area {
	int mu_steps, nu_steps, r_steps;
	double mu_min, nu_min, r_min;
	double mu_step_size, nu_step_size, r_step_size;
	int mu_step_current, nu_step_current, r_step_current;	

	double background_integral, *stream_integrals;
} INTEGRAL_AREA;

typedef struct evaluation_state {
	/********
		*	State for integral calculation.
	 ********/
	INTEGRAL_AREA *main_integral, **cuts;
	int current_cut;

	double background_integral;
	double* stream_integrals;

	/********
		*	State for likelihood calculation.
	 ********/
	int current_star_point;
	int num_zero;
	int bad_jacobians;
	double prob_sum;
} EVALUATION_STATE;

void	initialize_state(ASTRONOMY_PARAMETERS *ap, EVALUATION_STATE* es);
void	free_state(ASTRONOMY_PARAMETERS *ap, EVALUATION_STATE* es);
void	reset_evaluation_state(ASTRONOMY_PARAMETERS *ap, EVALUATION_STATE *es);

int	write_checkpoint(EVALUATION_STATE* es);
int	read_checkpoint(EVALUATION_STATE* es);

#endif
