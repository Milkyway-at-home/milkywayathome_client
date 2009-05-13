extern "C" {
#include "../astronomy/parameters.h"

#include "cpu_integrals.h"
#include "cpu_coords.h"
#include "cpu_r_constants.h"
}

#include <gpu_integrals.cu>
#include <gpu_coords.cu>
#include <gpu_r_constants.cu>

#include <cuda.h>

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
	int i, j;
	ASTRONOMY_PARAMETERS *ap;
	INTEGRAL *integral;

	ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	read_astronomy_parameters(argv[1], ap);
	integral = ap->integral[0];

//	ap->convolve = 120;
//	integral->r_steps = 350;
//	integral->nu_steps = 160;
//	integral->mu_steps = 200;

	cudaSetDevice( 0 );

	double *cpu__sinb, *cpu__sinl, *cpu__cosb, *cpu__cosl;
	cpu__gc_to_gal(ap->wedge, integral, &cpu__sinb, &cpu__sinl, &cpu__cosb, &cpu__cosl);

	double *cpu__V, *cpu__r_point, *cpu__qw_r3_N;
	cpu__r_constants(ap->convolve, integral, &cpu__V, &cpu__r_point, &cpu__qw_r3_N);

	printf("calculating integrals\n");

	unsigned int free_memory;
	unsigned int total_memory;
	cuMemGetInfo(&free_memory, &total_memory);
	printf("free memory: %d, total memory: %d (before initialize)\n", free_memory, total_memory);

	gpu__initialize_constants(ap, integral, cpu__V, cpu__r_point, cpu__qw_r3_N, cpu__sinb, cpu__sinl, cpu__cosb, cpu__cosl);

	cuMemGetInfo(&free_memory, &total_memory);
	printf("free memory: %d, total memory: %d (after initialize)\n", free_memory, total_memory);

	double gpu_bg_int, *gpu_st_int;
	gpu__integrals(ap, integral, &gpu_bg_int, &gpu_st_int);

	cuMemGetInfo(&free_memory, &total_memory);
	printf("free memory: %d, total memory: %d (after integrals)\n", free_memory, total_memory);

	gpu__free_constants();
	cuMemGetInfo(&free_memory, &total_memory);
	printf("free memory: %d, total memory: %d (after free)\n", free_memory, total_memory);

	printf("gpu completed\n");

	double *background_integrals, *stream_integrals;

	double bg_int = 0.0;
	double *st_int = (double*)malloc(ap->number_streams * sizeof(double));

	for (i = 0; i < ap->number_streams; i++) {
		st_int[i] = 0.0;
	}

	for (i = 0; i < integral->r_steps * integral->mu_steps * integral->nu_steps; i++) {
		bg_int += background_integrals[i];
		for (j = 0; j < ap->number_streams; j++) {
			st_int[j] += stream_integrals[(i * ap->number_streams) + j];
//			printf("    st_int[%d]: %.15lf\n", (i * ap->number_streams) + j, stream_integrals[(i * ap->number_streams) + j]);
		}
	}
	printf("     background_integral_sum: %.15lf\n", bg_int);
	for (i = 0; i < ap->number_streams; i++) {
		printf("     stream_integral_sum[%d]: %.15lf\n", i, st_int[i]);
	}
}


