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

#include "evaluation_state.h"
#include "parameters.h"
#include "star_points.h"

void init_constants(ASTRONOMY_PARAMETERS* ap);
void free_constants(ASTRONOMY_PARAMETERS* ap);
void set_probability_constants(int n_convolve,
                               double coords,
                               double* r_point,
                               double* r_in_mag,
                               double* r_in_mag2,
                               double* qw_r3_N,
                               double* reff_xr_rp3);

void calculate_probabilities(double* r_point,
                             double* r_in_mag,
                             double* r_in_mag2,
                             double* qw_r3_N,
                             double reff_xr_rp3,
                             double* integral_point,
                             const ASTRONOMY_PARAMETERS* ap,
                             double* bg_prob,
                             double* st_prob);

int calculate_integrals(const ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, const STAR_POINTS* sp);
int calculate_likelihood(const ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, const STAR_POINTS* sp);

void cpu__r_constants(  int n_convolve,
                        int r_steps, double r_min, double r_step_size,
                        int mu_steps, double mu_min, double mu_step_size,
                        int nu_steps, double nu_min, double nu_step_size,
                        double* irv, double** r_point, double** r_in_mag, double** r_in_mag2, double** qw_r3_N,
                        double* reff_xr_rp3, double* nus, double* ids);

double cpu_evaluate(double* parameters,
                    ASTRONOMY_PARAMETERS* ap,
                    EVALUATION_STATE* es,
                    STAR_POINTS* sp);

#endif

