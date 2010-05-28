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
	int r_cut_step_current;
	int mu_cut_step_current;
	int nu_cut_step_current;
	int number_streams;
	int current_cut;
	int main_integral_calculated;
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
void	reset_evaluation_state(EVALUATION_STATE *es);

int	write_checkpoint(EVALUATION_STATE* es);
int	read_checkpoint(EVALUATION_STATE* es);

int	calculate_integrals(ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, STAR_POINTS* sp);
int	calculate_likelihood(ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, STAR_POINTS* sp);

#endif


