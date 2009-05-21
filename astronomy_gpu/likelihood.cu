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

extern "C" {
#include "../astronomy/parameters.h"
#include "../astronomy/star_points.h"

#include "cpu_integrals.h"
#include "cpu_coords.h"
#include "cpu_r_constants.h"
}

#include <cuda.h>
#include <evaluation_gpu.cu>

void compare_arrays(const char *array_name, int length, float *gpu, double *cpu) {
	int i;
	double min_difference, max_difference, sum_difference, avg_difference, diff;
	double min_percentage, max_percentage, sum_percentage, avg_percentage, perc;

	diff = (double)gpu[0] - cpu[0];
	perc = fabs(diff)/fabs(cpu[0]);
	min_difference = diff;
	min_percentage = perc;
	max_difference = diff;
	max_percentage = perc;
	sum_difference = fabs(diff);
	sum_percentage = perc;

	for (i = 1; i < length; i++) {
		diff = (double)gpu[i] - cpu[i];
		perc = fabs(diff)/fabs(cpu[i]);

		if (diff < min_difference) min_difference = diff;
		if (perc < min_percentage) min_percentage = perc;
		if (diff > max_difference) max_difference = diff;
		if (perc > max_percentage) max_percentage = perc;
		sum_difference += fabs(diff);
		sum_percentage += perc;
	}
	avg_difference = sum_difference / length;
	avg_percentage = sum_percentage / length;

	printf("[%s] min_difference: %.20lf, max_difference: %.20lf, avg_difference: %.20lf\n", array_name, min_difference, max_difference, avg_difference);
	printf("[%s] min_percentage: %.20lf, max_percentage: %.20lf, avg_percentage: %.20lf\n", array_name, min_percentage, max_percentage, avg_percentage);
}


int main (int argc, char **argv) {
	ASTRONOMY_PARAMETERS *ap;
	STAR_POINTS *sp;
	double *parameters;

	ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	read_astronomy_parameters(argv[1], ap);

	sp = (STAR_POINTS*)malloc(sizeof(STAR_POINTS));
	read_star_points(argv[2], sp);

	cudaSetDevice( 0 );

	unsigned int free_memory;
	unsigned int total_memory;
	cuMemGetInfo(&free_memory, &total_memory);
	printf("free memory: %d, total memory: %d (before initialize)\n", free_memory, total_memory);

	gpu__initialize(ap, sp);
	get_search_parameters(ap, &parameters);

	cuMemGetInfo(&free_memory, &total_memory);
	printf("free memory: %d, total memory: %d (after initialize)\n", free_memory, total_memory);

	double likelihood = gpu__likelihood(parameters);

	cuMemGetInfo(&free_memory, &total_memory);
	printf("free memory: %d, total memory: %d (after integrals)\n", free_memory, total_memory);

	gpu__free_constants();
	cuMemGetInfo(&free_memory, &total_memory);
	printf("free memory: %d, total memory: %d (after free)\n", free_memory, total_memory);

	printf("gpu completed\n");
}

