
#include <r_constants__kernel.cu>
#include <coords__kernel.cu>

__device__ float **device__reff_xr_rp3;
__device__ float **device__r_point;
__device__ float **device__qw_r3_N;

__device__ float **device__glong;
__device__ float **device__glat;


double **gold__reff_xr_rp3;
double **gold__r_point;
double **gold__qw_r3_N;

double **gold__glong;
double **gold__glat;

void gold__initialize_constants(ASTRONOMY_PARAMETERS *ap) {
	int i;
	INTEGRAL *integral;

	gold__reff_xr_rp3 = (double**)malloc(ap->number_integrals * sizeof(double*));
	gold__r_point = (double**)malloc(ap->number_integrals * sizeof(double*));
	gold__qw_r3_N = (double**)malloc(ap->number_integrals * sizeof(double*));

	gold__glong = (double**)malloc(ap->number_integrals * sizeof(double*));
	gold__glat = (double**)malloc(ap->number_integrals * sizeof(double*));

	for (i = 0; i < ap->number_integrals; i++) {
		integral = ap->integrals[i];

		gold__r_constants(ap->convolve, integral->r_steps, integral->r_min, integral->r_step_size, &(gold__reff_xr_rp3[i]), &(gold__r_point[i]), &(gold__qw_r3_N[i]));
		gold__gc_to_gal(ap->wedge, integral->mu_steps, integral->mu_min, integral->mu_step_size, integral->nu_steps, integral->nu_min, integral->nu_step_size, &(gold__glong[i]), &(gold__glat[i]));
	}
}

double gold__milkyway_evaluate(double *parameters) {
}

void gpu__initialize_constants(ASTRONOMY_PARAMETERS *ap) {
	int i;
	INTEGRAL *integral;

	cudaMalloc((void**) &(device__reff_xr_rp3), ap->number_integrals * sizeof(float*));
	cudaMalloc((void**) &(device__r_point), ap->number_integrals * sizeof(float*));
	cudaMalloc((void**) &(device__qw_r3_N), ap->number_integrals * sizeof(float*));

	cudaMalloc((void**) &(device__glong), ap->number_integrals * sizeof(float*));
	cudaMalloc((void**) &(device__glat), ap->number_integrals * sizeof(float*));

	for (i = 0; i < ap->number_integrals; i++) {
		integral = ap->integrals[i];

		gpu__r_constants(ap->convolve, integral->r_steps, integral->r_min, integral->r_step_size, &(device__reff_xr_rp3[i]), &(device__r_point[i]), &(device__qw_r3_N[i]));
		gpu__gc_to_gal(ap->wedge, integral->mu_steps, integral->mu_min, integral->mu_step_size, integral->nu_steps, integral->nu_min, integral->nu_step_size, &(device__glong[i]), &(device__glat[i]));
	}
}

double gpu__milkyway_evaluate(double *parameters) {
}

