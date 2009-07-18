/*
 * Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
 * Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
 * and Rensselaer Polytechnic Institute.
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 * */

#ifndef EVALUATION_GPU_H
#define EVALUATION_GPU_H

#include <cuda.h>
#include <cuda_runtime.h>

#include "../astronomy/star_points.h"
#include "../astronomy/parameters.h"

void choose_gpu();
void gpu__initialize(	int ap_sgr_coordinates, int ap_wedge, int ap_convolve, int ap_number_streams, int ap_number_integrals, 
			int *in__r_steps, double *r_min, double *r_step_size,
			int *in__mu_steps, double *mu_min, double *mu_step_size,
			int *in__nu_steps, double *nu_min, double *nu_step_size,
			int in__number_stars, double **stars);

double gpu__likelihood(double *parameters);
void gpu__free_constants();

#endif
