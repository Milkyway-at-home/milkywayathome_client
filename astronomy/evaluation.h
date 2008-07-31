#ifndef ASTRONOMY_EVALUATION_H
#define ASTRONOMY_EVALUATION_H

#include "parameters.h"
#include "star_points.h"

typedef struct evaluation_state {
	/********
		*	State for integral calculation.
	 ********/
	int r_step_current;
	int mu_step_current;
	int nu_step_current;
	int number_streams;
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

void	initialize_state(EVALUATION_STATE* es, int number_streams);
void	free_state(EVALUATION_STATE* es);

int	write_checkpoint(EVALUATION_STATE* es);
int	read_checkpoint(EVALUATION_STATE* es);

int	calculate_integrals(ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, STAR_POINTS* sp);
int	calculate_likelihood(ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, STAR_POINTS* sp);

#endif


