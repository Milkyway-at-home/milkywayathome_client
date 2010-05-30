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

#ifndef FGDO_MPI_WORKER_ASTRONOMY_H
#define FGDO_MPI_WORKER_ASTRONOMY_H

/********
	*	Includes for astronomy
 ********/
#include "parameters.h"
#include "star_points.h"
#include "evaluation_optimized.h"

extern char astronomy_parameters_file[1024];
extern char star_points_file[1024];

extern ASTRONOMY_PARAMETERS *ap;
extern STAR_POINTS *sp;
extern EVALUATION_STATE *es;

void read_data(int rank, int max_rank);

void integral_f(double* parameters, double* results);
void integral_compose(double* integral_results, int num_results, double* results);

void likelihood_f(double* integrals, double* results);
double likelihood_compose(double* results, int num_results);

#endif

