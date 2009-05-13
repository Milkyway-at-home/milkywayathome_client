extern "C" {
#include "../astronomy/parameters.h"
#include "coords.h"
#include "cpu_coords.h"
#include "r_constants.h"
#include "pi_constants.h"
#include "time.h"
}

#include <cutil_inline.h>

#include "/Developer/CUDA/projects/reduction/reduction_kernel_sm10.cu"

float *device__V;
float *device__r_point, *device__qw_r3_N;
float *device__sinb, *device__sinl, *device__cosb, *device__cosl;

float *device__background_integrals, *device__stream_integrals;
float *host__background_integrals;
float *host__stream_integrals;

float *device__reduce;
float *host__reduce;


void gpu__initialize_constants(ASTRONOMY_PARAMETERS *ap, INTEGRAL *integral, double *cpu__V, double *cpu__r_point, double *cpu__qw_r3_N, double *cpu__sinb, double *cpu__sinl, double *cpu__cosb, double *cpu__cosl) {
	int i; 
	int sizeof_V = integral->nu_steps * integral->r_steps;
	int sizeof_r_constants = integral->r_steps * ap->convolve;
	int sizeof_g_constants = integral->mu_steps * integral->nu_steps;
	float *host__V			= (float*)malloc(sizeof_V * sizeof(float));
	float *host__r_point		= (float*)malloc(sizeof_r_constants * sizeof(float));
	float *host__qw_r3_N		= (float*)malloc(sizeof_r_constants * sizeof(float));
	float *host__sinb		= (float*)malloc(sizeof_g_constants * sizeof(float));
	float *host__sinl		= (float*)malloc(sizeof_g_constants * sizeof(float));
	float *host__cosb		= (float*)malloc(sizeof_g_constants * sizeof(float));
	float *host__cosl		= (float*)malloc(sizeof_g_constants * sizeof(float));

	long constants_size = 0;
	constants_size += sizeof_V * sizeof(float);
	constants_size += sizeof_r_constants * sizeof(float) * 2; 
	constants_size += sizeof_g_constants * sizeof(float) * 4;

	printf("Allocating %ld bytes for constants on GPU.\n", constants_size);

	for (i = 0; i < sizeof_V; i++) {
		host__V[i] = (float)cpu__V[i];
	}
	for (i = 0; i < sizeof_r_constants; i++) {
		host__r_point[i] = (float)cpu__r_point[i];
		host__qw_r3_N[i] = (float)cpu__qw_r3_N[i];
	}
	for (i = 0; i < sizeof_g_constants; i++) {
		host__sinb[i] = (float)cpu__sinb[i];
		host__sinl[i] = (float)cpu__sinl[i];
		host__cosb[i] = (float)cpu__sinb[i];
		host__cosl[i] = (float)cpu__sinl[i];
	}

	cutilSafeCall( cudaMalloc((void**) &device__V, sizeof_V * sizeof(float)) );
	cutilSafeCall( cudaMalloc((void**) &device__r_point, sizeof_r_constants * sizeof(float)) );
	cutilSafeCall( cudaMalloc((void**) &device__qw_r3_N, sizeof_r_constants * sizeof(float)) );

	cutilSafeCall( cudaMalloc((void**) &device__sinb, sizeof_g_constants * sizeof(float)) );
	cutilSafeCall( cudaMalloc((void**) &device__sinl, sizeof_g_constants * sizeof(float)) );
	cutilSafeCall( cudaMalloc((void**) &device__cosb, sizeof_g_constants * sizeof(float)) );
	cutilSafeCall( cudaMalloc((void**) &device__cosl, sizeof_g_constants * sizeof(float)) );

	cutilSafeCall( cudaMemcpy(device__V, host__V, sizeof_V * sizeof(float), cudaMemcpyHostToDevice) );
	cutilSafeCall( cudaMemcpy(device__r_point, host__r_point, sizeof_r_constants * sizeof(float), cudaMemcpyHostToDevice) );
	cutilSafeCall( cudaMemcpy(device__qw_r3_N, host__qw_r3_N, sizeof_r_constants * sizeof(float), cudaMemcpyHostToDevice) );

	cutilSafeCall( cudaMemcpy(device__sinb, host__sinb, sizeof_g_constants * sizeof(float), cudaMemcpyHostToDevice) );
	cutilSafeCall( cudaMemcpy(device__sinl, host__sinl, sizeof_g_constants * sizeof(float), cudaMemcpyHostToDevice) );
	cutilSafeCall( cudaMemcpy(device__cosb, host__cosb, sizeof_g_constants * sizeof(float), cudaMemcpyHostToDevice) );
	cutilSafeCall( cudaMemcpy(device__cosl, host__cosl, sizeof_g_constants * sizeof(float), cudaMemcpyHostToDevice) );

	int integral_size = ap->convolve * integral->r_steps * integral->mu_steps; 
	printf("Allocating %d bytes for integral data on GPU\n", integral_size * sizeof(float));

	cutilSafeCall( cudaMalloc((void**) &device__background_integrals, integral_size * sizeof(float)) );
	cutilSafeCall( cudaMalloc((void**) &device__stream_integrals, ap->number_streams * integral_size * sizeof(float)) );
	host__background_integrals = (float*)malloc(integral_size * sizeof(float));
	host__stream_integrals = (float*)malloc(ap->number_streams * integral_size * sizeof(float));

	cutilSafeCall( cudaMalloc((void**) &device__reduce, 64 * sizeof(float)) );
	host__reduce = (float*)malloc(64 * sizeof(float));
}

void gpu__free_constants() {
	cutilSafeCall( cudaFree(device__V) );
	cutilSafeCall( cudaFree(device__r_point) );
	cutilSafeCall( cudaFree(device__qw_r3_N) );
	cutilSafeCall( cudaFree(device__sinb) );
	cutilSafeCall( cudaFree(device__sinl) );
	cutilSafeCall( cudaFree(device__cosb) );
	cutilSafeCall( cudaFree(device__cosl) );
	cutilSafeCall( cudaFree(device__background_integrals) );
	cutilSafeCall( cudaFree(device__stream_integrals) );
	cutilSafeCall( cudaFree(device__reduce) );
	free(host__background_integrals);
	free(host__stream_integrals);
	free(host__reduce);
}

template <unsigned int convolve>
__global__ void gpu__integrals_kernel	(
					float alpha, float q, float r0, float delta,										//background parameters
					float stream_sigma_sq2, float stream_a[3], float stream_c[3],		//stream parameters
					float *device__sinb, float *device__sinl, float *device__cosb, float *device__cosl,					//mu and nu constants
					float *device__V, float *device__r_point, float *device__qw_r3_N,							//r constants
					float *background_integrals, float *stream_integrals									//output
					) {

	__shared__ float bg_integrals[convolve];
	__shared__ float st_integrals[convolve];

//	int rc_position = (r_step * n_convolve) + convolve;
	int rc_position = (blockIdx.x * blockDim.x) + threadIdx.x;

	float xyz2 = device__r_point[rc_position] * device__sinb[blockIdx.y];
	float zp = device__r_point[rc_position] * device__cosb[blockIdx.y];
	float xyz0 = zp * device__cosl[blockIdx.y] - f_lbr_r;
	float xyz1 = zp * device__sinl[blockIdx.y];

	float rg = sqrt( (xyz0 * xyz0) + (xyz1 * xyz1) + (xyz2 * xyz2)/(q * q) );
	float rs = rg + r0;

	bg_integrals[threadIdx.x] = device__qw_r3_N[rc_position] / (rg * rs * rs * rs);

	float sxyz0 = xyz0 - stream_c[0];
	float sxyz1 = xyz1 - stream_c[1];
	float sxyz2 = xyz2 - stream_c[2];

	float dotted = stream_a[0] * sxyz0 + stream_a[1] * sxyz1 + stream_a[2] * sxyz2;

	sxyz0 -= dotted * stream_a[0];
	sxyz1 -= dotted * stream_a[1];
	sxyz2 -= dotted * stream_a[2];

	float xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);

	st_integrals[threadIdx.x] = device__qw_r3_N[rc_position] * exp(-xyz_norm / stream_sigma_sq2);

	//Perform a reduce on the convolve data
	__syncthreads();
	if (threadIdx.x < 64 && (threadIdx.x + 64 < blockDim.x)) {
		bg_integrals[threadIdx.x] += bg_integrals[threadIdx.x + 64];
		st_integrals[threadIdx.x] += st_integrals[threadIdx.x + 64];
	}
	if (threadIdx.x < 32) {
		bg_integrals[threadIdx.x] += bg_integrals[threadIdx.x + 32];
		st_integrals[threadIdx.x] += st_integrals[threadIdx.x + 32];
	}
	if (threadIdx.x < 16) {
		bg_integrals[threadIdx.x] += bg_integrals[threadIdx.x + 16];
		st_integrals[threadIdx.x] += st_integrals[threadIdx.x + 16];
	}
	if (threadIdx.x < 8) {
		bg_integrals[threadIdx.x] += bg_integrals[threadIdx.x + 8];
		st_integrals[threadIdx.x] += st_integrals[threadIdx.x + 8];
	}
	if (threadIdx.x < 4) {
		bg_integrals[threadIdx.x] += bg_integrals[threadIdx.x + 4];
		st_integrals[threadIdx.x] += st_integrals[threadIdx.x + 4];
	}
	if (threadIdx.x < 2) {
		bg_integrals[threadIdx.x] += bg_integrals[threadIdx.x + 2];
		st_integrals[threadIdx.x] += st_integrals[threadIdx.x + 2];
	}
	if (threadIdx.x < 1) {
		bg_integrals[threadIdx.x] += bg_integrals[threadIdx.x + 1];
		st_integrals[threadIdx.x] += st_integrals[threadIdx.x + 1];
	}
	//multiple the result by V
	if (threadIdx.x == 0) {
		background_integrals[blockIdx.x + (blockIdx.y * blockDim.x)] = bg_integrals[0] * device__V[blockIdx.x];
		stream_integrals[blockIdx.x + (blockIdx.y * blockDim.x)] = st_integrals[0] * device__V[blockIdx.x];
	}
}

void gpu__integrals(ASTRONOMY_PARAMETERS *ap, INTEGRAL *integral, double *background_integral, double **stream_integrals) {
	int i;
	double alpha = ap->background_parameters[0];
	double q = ap->background_parameters[1];
	double r0 = ap->background_parameters[2];
	double delta = ap->background_parameters[3];

	float fstream_sigma_sq2[ap->number_streams], fstream_a[ap->number_streams][3], fstream_c[ap->number_streams][3];
	double stream_c[3], lbr[3];
	for (i = 0; i < ap->number_streams; i++) {
		fstream_sigma_sq2[i] = (float)(2.0 * ap->stream_parameters[i][4] * ap->stream_parameters[i][4]);

		fstream_a[i][0] = (float)( sin(ap->stream_parameters[i][2]) * cos(ap->stream_parameters[i][3]) );
		fstream_a[i][1] = (float)( sin(ap->stream_parameters[i][2]) * sin(ap->stream_parameters[i][3]) );
		fstream_a[i][2] = (float)( cos(ap->stream_parameters[i][2]) );

		gc_to_gal(ap->wedge, ap->stream_parameters[i][0] * D_DEG2RAD, 0 * D_DEG2RAD, &(lbr[0]), &(lbr[1]));
		lbr[2] = ap->stream_parameters[i][1];
		d_lbr2xyz(lbr, stream_c);
	
		fstream_c[i][0] = (float)stream_c[0];
		fstream_c[i][1] = (float)stream_c[1];
		fstream_c[i][2] = (float)stream_c[2];
	}

	dim3 dimGrid(integral->r_steps, integral->mu_steps);

	switch (ap->convolve) {
		case 120:	gpu__integrals_kernel<120><<<dimGrid, ap->convolve>>>(	(float)alpha, (float)q, (float)r0, (float)delta,
											fstream_sigma_sq2[0], fstream_a[0], fstream_c[0],
											device__sinb, device__sinl, device__cosb, device__cosl,
											device__V, device__r_point, device__qw_r3_N,
											device__background_integrals, device__stream_integrals);
				break;
	}
	*stream_integrals = (double*)malloc(ap->number_streams * sizeof(double*));
}

/*
__global__ void gpu__background_integrals(int convolve, int r_step, int r_steps, float alpha, float q, float r0, float delta, 
		float *device__glong, float *device__glat, float *device__V, float *device__r_point, float *device__qw_r3_N,
		float *background_integral) {
//	int r_step, r_steps;
	int nu_step, nu_steps, mu_step;
	float V, *qw_r3_N, *r_point, glong, glat;

	mu_step = blockIdx.x;
//	r_step = blockIdx.x;
//	r_steps = gridDim.x;
	nu_steps = blockDim.x;
	nu_step = threadIdx.x;

	glong = device__glong[(mu_step * nu_steps) + nu_step];
	glat = device__glat[(mu_step * nu_steps) + nu_step];

	float sinb = sin(glat);
	float sinl = sin(glong);
	float cosb = cos(glat);
	float cosl = cos(glong);

	V = device__V[(nu_step * r_steps) + r_step];
	qw_r3_N = &( device__qw_r3_N[r_step * convolve] );
	r_point = &( device__r_point[r_step * convolve] );

	float zp, rg, rs;
	float xyz0, xyz1, xyz2;

	int bg_pos = r_step + (nu_step * r_steps) + (mu_step * nu_steps * r_steps);

	float bg_int = 0.0;
//	float pow_rg, pow_rg_r0;
	for (int i = 0; i < convolve; i++) {
		xyz2 = r_point[i] * sinb;
		zp = r_point[i] * cosb;
		xyz0 = zp * cosl - f_lbr_r;
		xyz1 = zp * sinl;

		rg = sqrt(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2)/(q*q));
		rs = rg + r0;
		bg_int += qw_r3_N[i] / (rg * rs * rs * rs);

//		pow_rg = powf(rg, alpha);
//		pow_rg_r0 = powf(rg + r0, 3.0f - alpha + delta);
//		bg_int += pow_rg * pow_rg_r0;

//		bg_int += qw_r3_N[i] / (powf(rg, alpha) * powf(rg + r0, 3.0f - alpha + delta));
	}
	background_integral[bg_pos] = bg_int * V;
}

__global__ void gpu__stream_integrals(int convolve, int r_step, int r_steps, float stream_sigma_sq2, float stream_a0, float stream_a1, float stream_a2, float stream_c0, float stream_c1, float stream_c2, 
		float *device__glong, float *device__glat, float *device__V, float *device__r_point, float *device__qw_r3_N,
		float *stream_integrals) {
	int i;
//	int r_step, r_steps;
	int nu_step, nu_steps, mu_step;
	float V, *qw_r3_N, *r_point, glong, glat;

	mu_step = blockIdx.x;
//	r_step = blockIdx.x;
//	r_steps = gridDim.x;
	nu_steps = blockDim.x;
	nu_step = threadIdx.x;

	glong = device__glong[(mu_step * nu_steps) + nu_step];
	glat = device__glat[(mu_step * nu_steps) + nu_step];

	V = device__V[(nu_step * r_steps) + r_step];
	qw_r3_N = &( device__qw_r3_N[r_step * convolve] );
	r_point = &( device__r_point[r_step * convolve] );

	float sinb = sin(glat);
	float sinl = sin(glong);
	float cosb = cos(glat);
	float cosl = cos(glong);

	float xyz0, xyz1, xyz2;
	float sxyz0, sxyz1, sxyz2;
	float zp, dotted, xyz_norm;

	int st_pos = r_step + (nu_step * r_steps) + (mu_step * nu_steps * r_steps);

	float st_int = 0.0;
	for (i = 0; i < convolve; i++) {
		xyz2 = r_point[i] * sinb;
		zp = r_point[i] * cosb;
		xyz0 = zp * cosl - f_lbr_r;
		xyz1 = zp * sinl;

		sxyz0 = xyz0 - stream_c0;
		sxyz1 = xyz1 - stream_c1;
		sxyz2 = xyz2 - stream_c2;

		dotted = stream_a0 * sxyz0 + stream_a1 * sxyz1 + stream_a2 * sxyz2;

		sxyz0 -= dotted * stream_a0;
		sxyz1 -= dotted * stream_a1;
		sxyz2 -= dotted * stream_a2;

		xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);

		st_int += qw_r3_N[i] * exp(-xyz_norm / stream_sigma_sq2);
	}
	stream_integrals[st_pos] = st_int * V;
}
*/

/*
double sum_integrals(int integral_size) {
	int i;
	time_t start_time, finish_time;

	printf("starting reduce\n");

	time(&start_time);
	reduce_sm10<float>(integral_size, 128, 64, 6, device__integrals, device__reduce);

	printf("did reduce\n");

	cudaMemcpy(host__reduce, device__reduce, 64 * sizeof(float), cudaMemcpyDeviceToHost);
	time(&finish_time);
	printf("reduce took %ld seconds\n", (long)(finish_time-start_time));
	double sum = 0.0;
	for (i = 0; i < 64; i++) sum += host__reduce[i];
	printf("reduce sum: %lf\n", sum);
	return sum;
//	time(&start_time);
//	cudaMemcpy(host__integrals, device__integrals, integral_size * sizeof(float), cudaMemcpyDeviceToHost);
//	time(&finish_time);
//	printf("transfer took %ld seconds\n", (long)(finish_time-start_time));

//	double integral = 0.0;
//	for (i = 0; i < integral_size; i++) {
//		integral += (double)host__integrals[i];
//	}
//	return integral;
}
*/

/*
void gpu__integrals(ASTRONOMY_PARAMETERS *ap, INTEGRAL *integral, double *background_integral, double **stream_integrals) {
	int i, j;
	double alpha = ap->background_parameters[0];
	double q = ap->background_parameters[1];
	double r0 = ap->background_parameters[2];
	double delta = ap->background_parameters[3];
	time_t start_time, finish_time;

	int integral_size = integral->r_steps * integral->mu_steps * integral->nu_steps;

	dim3 dimGrid(integral->r_steps, integral->mu_steps);
	time(&start_time);
	for (i = 0; i < integral->r_steps; i++) {
		gpu__background_integrals<<<integral->mu_steps, integral->nu_steps>>>(ap->convolve, i, integral->r_steps, alpha, q, r0, delta, device__glong, device__glat, device__V, device__r_point, device__qw_r3_N, device__integrals);
	}
	time(&finish_time);
	printf("calculate took %ld seconds\n", (long)(finish_time - start_time));

	*background_integral = sum_integrals(integral_size);

	printf("gpu__background_integral_sum: %.15lf\n", *background_integral);

	double stream_sigma_sq2, stream_a[3], stream_c[3], lbr[3];
	*stream_integrals = (double*)malloc(ap->number_streams * sizeof(double*));
	for (i = 0; i < ap->number_streams; i++) {
		stream_sigma_sq2 = 2.0 * ap->stream_parameters[i][4] * ap->stream_parameters[i][4];

		stream_a[0] = sin(ap->stream_parameters[i][2]) * cos(ap->stream_parameters[i][3]);
		stream_a[1] = sin(ap->stream_parameters[i][2]) * sin(ap->stream_parameters[i][3]);
		stream_a[2] = cos(ap->stream_parameters[i][2]);

		gc_to_gal(ap->wedge, ap->stream_parameters[i][0] * D_DEG2RAD, 0 * D_DEG2RAD, &(lbr[0]), &(lbr[1]));
		lbr[2] = ap->stream_parameters[i][1];
		d_lbr2xyz(lbr, stream_c);

		time(&start_time);
		for (j = 0; j < integral->r_steps; j++) {
			gpu__stream_integrals<<<integral->mu_steps, integral->nu_steps>>>(ap->convolve, j, integral->r_steps, (float)stream_sigma_sq2, (float)stream_a[0], (float)stream_a[1], (float)stream_a[2], (float)stream_c[0], (float)stream_c[1], (float)stream_c[2], device__glong, device__glat, device__V, device__r_point, device__qw_r3_N, device__integrals);
		}
		time(&finish_time);
		printf("calculate took %ld seconds\n", (long)(finish_time - start_time));

		(*stream_integrals)[i] = sum_integrals(integral_size);
		printf("gpu__stream_integrals[%d]: %.15lf\n", i, (*stream_integrals)[i]);
	}
}
*/
