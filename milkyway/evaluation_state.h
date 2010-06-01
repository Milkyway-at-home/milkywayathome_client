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


#ifndef ASTRONOMY_EVALUATION_STATE_H
#define ASTRONOMY_EVALUATION_STATE_H

#include "parameters.h"
#include "star_points.h"

typedef struct integral_area
{
    int mu_steps, nu_steps, r_steps;
    double mu_min, nu_min, r_min;
    double mu_max, nu_max, r_max;
    double mu_step_size, nu_step_size, r_step_size;

    long min_calculation, max_calculation, current_calculation;

    int number_streams;
    double background_integral, *stream_integrals;
} INTEGRAL_AREA;

typedef struct evaluation_state
{
    /********
        *   State for integral calculation.
     ********/
    INTEGRAL_AREA** integral;
    int current_integral, number_streams, number_integrals;

    double background_integral;
    double* stream_integrals;

    /********
        *   State for likelihood calculation.
     ********/
    int current_star_point, total_stars;
    int num_zero;
    int bad_jacobians;
    double prob_sum;
} EVALUATION_STATE;

void    get_steps(INTEGRAL_AREA* ia, int* mu_step_current, int* nu_step_current, int* r_step_current);

void    initialize_state(ASTRONOMY_PARAMETERS* ap, STAR_POINTS* sp, EVALUATION_STATE* es);
void    free_state(EVALUATION_STATE* es);
void    reset_evaluation_state(EVALUATION_STATE* es);

void    fwrite_integral_area(FILE* file, INTEGRAL_AREA* ia);

int write_checkpoint(EVALUATION_STATE* es);
int read_checkpoint(EVALUATION_STATE* es);

#endif
