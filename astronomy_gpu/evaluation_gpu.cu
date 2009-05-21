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
#include "coords.h"
#include "cpu_coords.h"
#include "cpu_r_constants.h"
#include "r_constants.h"
#include "pi_constants.h"
#include "gauss_legendre.h"
}

#include <cuda.h>
#include <cuda_runtime.h>
#include <cutil_inline.h>

#include "/Developer/CUDA/projects/reduction/reduction_kernel_sm10.cu"

int	number_threads = 256;

int	mu_increment = 1;
int	wedge;
int	convolve;
int	number_streams;
int	number_integrals;

int	*r_steps;
//double	*r_min;
//double	*r_step_size;

int	*mu_steps;
int	*nu_steps;

int	*sizeof_V;
float	**device__V;		//V				-- float[nu][r]

int	*sizeof_r_constants;
float	**device__r_constants;	//r_point, qw_r3_N		-- float[r][convolve][2]

int	*sizeof_lb;
float	**device__lb;		//sinb, sinl, cosb, cosl	-- float[nu][mu][4]

int	*integral_size;
float	**host__background_integrals;
float	**host__stream_integrals;
float	**device__background_integrals;
float	**device__stream_integrals;

__device__ __constant__ float device__fstream_sigma_sq2[8];
__device__ __constant__ float device__fstream_a[12];
__device__ __constant__ float device__fstream_c[12];

__device__ __constant__ float device__dx[120];
__device__ __constant__ float device__qgaus_W[120];

__device__ __constant__ float device__background_weight[1];
__device__ __constant__ float device__stream_weight[4];

int	number_stars;
float	*device__stars;

int	probability_size;
float	*device__probability;
float	*host__probability;

float	*device__reduce;
float	*host__reduce;

//extern "C" void gpu__initialize(ASTRONOMY_PARAMETERS *ap, STAR_POINTS *sp);
//extern "C" double gpu__likelihood(double *parameters);
//extern "C" void gpu__free_constants();


void gpu__initialize(ASTRONOMY_PARAMETERS *ap, STAR_POINTS *sp) { 
	int i, j, pos;

	printf("astronomy parameters?\n");
	fwrite_astronomy_parameters(stdout, ap);

	printf("wedge: %d, convolve: %d, number_streams: %d, number_integrals: %d\n", wedge, convolve, number_streams, number_integrals);


	wedge = ap->wedge;
	convolve = ap->convolve;
	number_streams = ap->number_streams;
	number_integrals = ap->number_integrals;

	sizeof_V = (int*)malloc(number_integrals * sizeof(int));
	sizeof_r_constants = (int*)malloc(number_integrals * sizeof(int));
	sizeof_lb = (int*)malloc(number_integrals * sizeof(int));

	integral_size = (int*)malloc(number_integrals * sizeof(int));

	device__background_integrals = (float**)malloc(number_integrals * sizeof(float*));
	device__stream_integrals = (float**)malloc(number_integrals * sizeof(float*));
	host__background_integrals = (float**)malloc(number_integrals * sizeof(float*));
	host__stream_integrals = (float**)malloc(number_integrals * sizeof(float*));

	device__V = (float**)malloc(number_integrals * sizeof(float*));
	device__lb = (float**)malloc(number_integrals * sizeof(float*));
	device__r_constants = (float**)malloc(number_integrals * sizeof(float*));

	printf("calculating integral constants\n");
	r_steps = (int*)malloc(number_integrals * sizeof(int));
//	r_min = (double*)malloc(number_integrals * sizeof(double));
//	r_step_size = (double*)malloc(number_integrals * sizeof(double));

	mu_steps = (int*)malloc(number_integrals * sizeof(int));
	nu_steps = (int*)malloc(number_integrals * sizeof(int));
	for (i = 0; i < number_integrals; i++) {
		printf("first integral loop\n");
		if (ap->integral == NULL) printf("INTEGRAL == NULL!\n");
		printf("number integrals: %d\n", ap->number_integrals);
		INTEGRAL *integral = ap->integral[i];
		printf("got integral\n");


		r_steps[i] = integral->r_steps;
//		r_min[i] = integral->r_min;
//		r_step_size[i] = integral->r_step_size;
		mu_steps[i] = integral->mu_steps;
		nu_steps[i] = integral->nu_steps;

		sizeof_V[i] = integral->nu_steps * integral->r_steps;
		sizeof_r_constants[i] = integral->r_steps * ap->convolve * 2;
		sizeof_lb[i] = integral->mu_steps * integral->nu_steps * 4;

		double *cpu__V, *cpu__r_const, *cpu__lb;
		cpu__gc_to_lb(ap->wedge, integral, &cpu__lb);
		cpu__r_constants(ap->convolve, integral, &cpu__V, &cpu__r_const);

		float *host__V			= (float*)malloc(sizeof_V[i] * sizeof(float));
		float *host__r_constants	= (float*)malloc(sizeof_r_constants[i] * sizeof(float));
		float *host__lb			= (float*)malloc(sizeof_lb[i] * sizeof(float));

		long constants_size = 0;
		constants_size += sizeof_V[i] * sizeof(float);
		constants_size += sizeof_r_constants[i] * sizeof(float); 
		constants_size += sizeof_lb[i] * sizeof(float);

		printf("Allocating %ld bytes for constants on GPU.\n", constants_size);

		for (j = 0; j < sizeof_V[i]; j++) {
			host__V[j] = (float)cpu__V[j];
		}
		for (j = 0; j < sizeof_r_constants[i]; j++) {
			host__r_constants[j] = (float)cpu__r_const[j];
		}
		for (j = 0; j < sizeof_lb[i]; j++) {
			host__lb[j] = (float)cpu__lb[j];
		}

		printf("freeing cpu constants\n");
		free(cpu__V);
		free(cpu__r_const);
		free(cpu__lb);

		printf("device malloc\n");

		cutilSafeCall( cudaMalloc((void**) &(device__V[i]), sizeof_V[i] * sizeof(float)) );
		cutilSafeCall( cudaMalloc((void**) &(device__r_constants[i]), sizeof_r_constants[i] * sizeof(float)) );
		cutilSafeCall( cudaMalloc((void**) &(device__lb[i]), sizeof_lb[i] * sizeof(float)) );

		printf("device memcpy\n");

		cutilSafeCall( cudaMemcpy(device__V[i], host__V, sizeof_V[i] * sizeof(float), cudaMemcpyHostToDevice) );
		cutilSafeCall( cudaMemcpy(device__r_constants[i], host__r_constants, sizeof_r_constants[i] * sizeof(float), cudaMemcpyHostToDevice) );
		cutilSafeCall( cudaMemcpy(device__lb[i], host__lb, sizeof_lb[i] * sizeof(float), cudaMemcpyHostToDevice) );

		free(host__V);
		free(host__r_constants);
		free(host__lb);

		integral_size[i] = integral->r_steps * integral->nu_steps * mu_increment;
		printf("Allocating %d bytes for integral data on GPU\n", (number_streams + 1) * integral_size[i] * sizeof(float));

		cutilSafeCall( cudaMalloc((void**) &device__background_integrals[i], integral_size[i] * sizeof(float)) );
		cutilSafeCall( cudaMalloc((void**) &device__stream_integrals[i], number_streams * integral_size[i] * sizeof(float)) );
		host__background_integrals[i] = (float*)malloc(integral_size[i] * sizeof(float));
		host__stream_integrals[i] = (float*)malloc(number_streams * integral_size[i] * sizeof(float));
	}

	cutilSafeCall( cudaMalloc((void**) &device__reduce, 64 * sizeof(float)) );
	host__reduce = (float*)malloc(64 * sizeof(float));

	printf("initializing constants for %d stars\n", sp->number_stars);

	number_stars = sp->number_stars;
	float *host__stars = (float*)malloc(number_stars * 5 * sizeof(float));
	for (i = 0; i < number_stars; i++) {
		pos = i * 5;
		host__stars[pos] = (float)sin(sp->stars[i][1] * D_DEG2RAD);
		host__stars[pos + 1] = (float)sin(sp->stars[i][0] * D_DEG2RAD);
		host__stars[pos + 2] = (float)cos(sp->stars[i][1] * D_DEG2RAD);
		host__stars[pos + 3] = (float)cos(sp->stars[i][0] * D_DEG2RAD);
		host__stars[pos + 4] = (float)sp->stars[i][2];
	}
	printf("allocating %d bytes for device__stars\n", number_stars * 5 * sizeof(float));
	cutilSafeCall( cudaMalloc((void**) &device__stars, number_stars * 5 * sizeof(float)) );
	cutilSafeCall( cudaMemcpy(device__stars, host__stars, number_stars * 5 * sizeof(float), cudaMemcpyHostToDevice) );

	free(host__stars);

	printf("initializing qgaus constants\n");

	double *d_qgaus_W = (double*)malloc(convolve * sizeof(double));
	double *d_qgaus_X = (double*)malloc(convolve * sizeof(double));

	d_gauss_legendre(-1.0, 1.0, d_qgaus_X, d_qgaus_W, convolve);
	float *host__dx = (float*)malloc(convolve * sizeof(float));
	float *host__qgaus_W = (float*)malloc(convolve * sizeof(float));
	for (i = 0; i < convolve; i++) {
		host__dx[i] = (float)(3.0 * d_stdev * d_qgaus_X[i]);
		host__qgaus_W[i] = (float)d_qgaus_W[i];
	}
	free(d_qgaus_W);
	free(d_qgaus_X);

	cutilSafeCall( cudaMemcpyToSymbol(device__dx, host__dx, convolve * sizeof(float), 0, cudaMemcpyHostToDevice) );
	cutilSafeCall( cudaMemcpyToSymbol(device__qgaus_W, host__qgaus_W, convolve * sizeof(float), 0, cudaMemcpyHostToDevice) );

	//	free(host__dx);
	//	free(host__qgaus_W);

	printf("mallocing GPU bg and stream probability: %d bytes\n", number_threads * sizeof(float));

	probability_size = number_threads;
	cutilSafeCall( cudaMalloc((void**) &device__probability, probability_size * sizeof(float)) );

	printf("mallocing host bg and stream probability\n");

	host__probability = (float*)malloc(probability_size * sizeof(float));
}

void gpu__free_constants() {
	int i;
	for (i = 0; i < number_integrals; i++) {
		cutilSafeCall( cudaFree(device__V[i]) );
		cutilSafeCall( cudaFree(device__r_constants[i]) );
		cutilSafeCall( cudaFree(device__lb[i]) );
		cutilSafeCall( cudaFree(device__background_integrals[i]) );
		cutilSafeCall( cudaFree(device__stream_integrals[i]) );
		free(host__background_integrals[i]);
		free(host__stream_integrals[i]);
	}

	cutilSafeCall( cudaFree(device__stars) );

	free(host__background_integrals);
	free(host__stream_integrals);
	free(device__V);
	free(device__r_constants);
	free(device__lb);
	free(device__background_integrals);
	free(device__stream_integrals);

	free(r_steps);
	free(mu_steps);
	free(nu_steps);

	free(sizeof_V);
	free(sizeof_lb);
	free(sizeof_r_constants);
	free(integral_size);

	cutilSafeCall( cudaFree(device__probability) );
	free(host__probability);

	cutilSafeCall( cudaFree(device__reduce) );
	free(host__reduce);
}

template <unsigned int number_streams>
__global__ void gpu__zero_integrals(float *background_integrals, float *stream_integrals) {
	int pos = threadIdx.x + (blockIdx.x * blockDim.x) + (blockIdx.y * blockDim.x * gridDim.x);

	background_integrals[pos] = 0;
	for (int i = 0; i < number_streams; i++) stream_integrals[(i * gridDim.y * gridDim.x * blockDim.x) + pos] = 0;
}

template <unsigned int number_streams> 
__global__ void gpu__integral_kernel(	int convolve, int mu_step, int mu_steps,
		float q, float r0,
		float *device__lb, float *device__r_constants, float *device__V,
		float *background_integrals, float *stream_integrals) {
	int i, j;
	int pos = ((threadIdx.x * mu_steps) + mu_step + blockIdx.y) * 4;

	float sinb = device__lb[pos];
	float sinl = device__lb[pos + 1];
	float cosb = device__lb[pos + 2];
	float cosl = device__lb[pos + 3];

	float rg, xyz0, xyz1, xyz2;
	float dotted, sxyz0, sxyz1, sxyz2;

	float bg_int = 0.0;
	float st_int[number_streams];
	for (i = 0; i < number_streams; i++) st_int[i] = 0;
	float r_point, qw_r3_N;
	float zp, rs;

	for (i = 0; i < convolve; i++) {
		pos = 2 * (i + (blockIdx.x * convolve));
		r_point = device__r_constants[pos];
		qw_r3_N = device__r_constants[pos + 1];

		xyz2 = r_point * sinb;
		zp = r_point * cosb;
		xyz0 = zp * cosl - f_lbr_r;
		xyz1 = zp * sinl;

		rg = sqrt(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2)/(q*q));
		rs = rg + r0;
		bg_int += qw_r3_N / (rg * rs * rs * rs);

		/********
			*	Not sure why this doesn't work.
		 ********/
		//		pow_rg = powf(rg, alpha);
		//		pow_rg_r0 = powf(rg + r0, 3.0f - alpha + delta);
		//		bg_int += pow_rg * pow_rg_r0;

		//		bg_int += qw_r3_N[i] / (powf(rg, alpha) * powf(rg + r0, 3.0f - alpha + delta));

		for (j = 0; j < number_streams; j++) {
			pos = (j * 3);
			sxyz0 = xyz0 - device__fstream_c[pos];
			sxyz1 = xyz1 - device__fstream_c[pos + 1];
			sxyz2 = xyz2 - device__fstream_c[pos + 2];

			dotted = device__fstream_a[pos] * sxyz0 + device__fstream_a[pos + 1] * sxyz1 + device__fstream_a[pos + 2] * sxyz2;

			sxyz0 -= dotted * device__fstream_a[pos];
			sxyz1 -= dotted * device__fstream_a[pos + 1];
			sxyz2 -= dotted * device__fstream_a[pos + 2];

			st_int[j] += qw_r3_N * exp(-((sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2)) / device__fstream_sigma_sq2[j]);
		}
	}
	float V = device__V[blockIdx.x + (gridDim.x * threadIdx.x)];

	pos = threadIdx.x + (blockIdx.x * blockDim.x) + (blockIdx.y * blockDim.x * gridDim.x);

	background_integrals[pos] += bg_int * V;
	for (i = 0; i < number_streams; i++) stream_integrals[(i * gridDim.y * gridDim.x * blockDim.x) + pos] += st_int[i] * V;
}

/*
template <unsigned int number_streams>
__global__ void gpu__integral_kernel2(	int convolve, int mu_step, int mu_steps,
					float r_min, float r_step_size,
					float q, float r0,
					float coeff, 
					float *device__lb, float *device__V,
					float *background_integrals, float *stream_integrals) {
	int i;
	int pos = ((threadIdx.x * mu_steps) + mu_step + blockIdx.y) * 4;
	float sinb = device__lb[pos];
	float sinl = device__lb[pos + 1];
	float cosb = device__lb[pos + 2];
	float cosl = device__lb[pos + 3];

	float coords = r_min + (blockIdx.x * r_step_size);

	float rg, xyz0, xyz1, xyz2;
	float dotted, sxyz0, sxyz1, sxyz2;

	float bg_int = 0.0;
	float st_int[number_streams];
	for (i = 0; i < number_streams; i++) st_int[i] = 0.0;

	float gPrime = 5.0f * (log10(coords * 1000.0f) - 1.0f) + f_absm;
	float exponent = exp(sigmoid_curve_1 * (gPrime - sigmoid_curve_2));
	float reff_value = sigmoid_curve_0 / (exponent + 1);
	float rPrime3 = coords * coords * coords;

	float reff_xr_rp3 = reff_value * f_xr / rPrime3;

	float r_point, qw_r3_N;
	float zp, rs, g;

	for (i = 0; i < convolve; i++) {
		g = gPrime + device__dx[i];

		r_point = pow(10.0f, (g - f_absm)/5.0f + 1.0f) / 1000.0f;
		rPrime3 = r_point * r_point * r_point;

		qw_r3_N = device__qgaus_W[i] * rPrime3 * coeff * exp( -((g - gPrime) * (g - gPrime) / (2 * f_stdev * f_stdev)) );

		xyz2 = r_point * sinb;
		zp = r_point * cosb;
		xyz0 = zp * cosl - f_lbr_r;
		xyz1 = zp * sinl;

		rg = sqrt(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2)/(q*q));
		rs = rg + r0;

		bg_int += qw_r3_N / (rg * rs * rs * rs);

		for (int j = 0; j < number_streams; j++) {
			pos = (j * 3);
			sxyz0 = xyz0 - device__fstream_c[pos];
			sxyz1 = xyz1 - device__fstream_c[pos + 1];
			sxyz2 = xyz2 - device__fstream_c[pos + 2];

			dotted = device__fstream_a[pos] * sxyz0 + device__fstream_a[pos + 1] * sxyz1 + device__fstream_a[pos + 2] * sxyz2;

			sxyz0 -= dotted * device__fstream_a[pos];
			sxyz1 -= dotted * device__fstream_a[pos + 1];
			sxyz2 -= dotted * device__fstream_a[pos + 2];

			st_int[j] += qw_r3_N * exp(-((sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2)) / device__fstream_sigma_sq2[j]);
		}
	}
	float V = device__V[blockIdx.x + (gridDim.x * threadIdx.x)];

	pos = threadIdx.x + (blockIdx.x * blockDim.x) + (blockIdx.y * blockDim.x * gridDim.x);

	background_integrals[pos] += bg_int * V * reff_xr_rp3;
	for (i = 0; i < number_streams; i++) stream_integrals[(i * gridDim.y * gridDim.x * blockDim.x) + pos] += st_int[i] * V * reff_xr_rp3;
}
*/

void cpu__sum_integrals(int iteration, double *background_integral, double *stream_integrals) {
	int i, j;

	cutilSafeCall( cudaMemcpy(host__background_integrals[iteration], device__background_integrals[iteration], integral_size[iteration] * sizeof(float), cudaMemcpyDeviceToHost) );

	double sum = 0.0;
	for (i = 0; i < integral_size[iteration]; i++) {
		sum += (double)(host__background_integrals[iteration][i]);
	}
	if (iteration == 0) *background_integral = sum;
	else *background_integral -= sum;

	cutilSafeCall( cudaMemcpy(host__stream_integrals[iteration], device__stream_integrals[iteration], number_streams * integral_size[iteration] * sizeof(float), cudaMemcpyDeviceToHost) );
	for (i = 0; i < number_streams; i++) {
		sum = 0.0;
		for (j = 0; j < integral_size[iteration]; j++) {
			sum += (double)(host__stream_integrals[iteration][j + (i * integral_size[iteration])]);
		}
		if (iteration == 0) stream_integrals[i] = sum;
		else stream_integrals[i] -= sum;
	}
}

void gpu__sum_integrals(int iteration, double *background_integral, double *stream_integrals) {
	int i, j;
	reduce_sm10<float>(integral_size[iteration], 128, 64, 6, device__background_integrals[iteration], device__reduce);
	cutilSafeCall( cudaMemcpy(host__reduce, device__reduce, 64 * sizeof(float), cudaMemcpyDeviceToHost) );

	double reduce_sum = 0.0;
	for (i = 0; i < 64; i++) {
		reduce_sum += (double)(host__reduce[i]);
	}
	if (iteration == 0) *background_integral = reduce_sum;
	else *background_integral -= reduce_sum;

	for (i = 0; i < number_streams; i++) {
		reduce_sm10<float>(integral_size[iteration], 128, 64, 6, &(device__stream_integrals[iteration][i * integral_size[iteration]]), device__reduce);
		cutilSafeCall( cudaMemcpy(host__reduce, device__reduce, 64 * sizeof(float), cudaMemcpyDeviceToHost) );

		double reduce_sum = 0.0;
		for (j = 0; j < 64; j++) {
			reduce_sum += (double)(host__reduce[j]);
		}
		if (iteration == 0) stream_integrals[i] = reduce_sum;
		else stream_integrals[i] -= reduce_sum;
	}

}

/********
 *	Likelihood calculation
 ********/

template <unsigned int number_streams>
__global__ void gpu__zero_likelihood(int block_size, float *device__probability) {
	device__probability[threadIdx.x] = 0;
}

void cpu__sum_likelihood(int block_size, double *probability) {
	int i;

	cutilSafeCall( cudaMemcpy(host__probability, device__probability, probability_size * sizeof(float), cudaMemcpyDeviceToHost) );

//	*probability = 0.0;

	for (i = 0; i < block_size; i++) {
		*probability += host__probability[i];
//		printf("bg_prob: %.15f\n", host__probability[i]);
	}
}

template <unsigned int number_streams>
__global__ void gpu__likelihood_kernel(	int offset, int convolve,
		float q, float r0,
		float coeff, 
		float *device__stars,
		float *probability) {
	int i;
	int pos = (offset + threadIdx.x) * 5;
	float sinb = device__stars[pos];
	float sinl = device__stars[pos + 1];
	float cosb = device__stars[pos + 2];
	float cosl = device__stars[pos + 3];
	float coords = device__stars[pos + 4];

	float rg, xyz0, xyz1, xyz2;
	float dotted, sxyz0, sxyz1, sxyz2;

	float bg_int = 0.0;
	float st_int[number_streams];
	for (i = 0; i < number_streams; i++) st_int[i] = 0.0;

	float gPrime = 5.0f * (log10(coords * 1000.0f) - 1.0f) + f_absm;
	float exponent = exp(sigmoid_curve_1 * (gPrime - sigmoid_curve_2));
	float reff_value = sigmoid_curve_0 / (exponent + 1);
	float rPrime3 = coords * coords * coords;

	float reff_xr_rp3 = reff_value * f_xr / rPrime3;

	float r_point, qw_r3_N;
	float zp, rs, g;

	for (i = 0; i < convolve; i++) {
		g = gPrime + device__dx[i];

		r_point = pow(10.0f, (g - f_absm)/5.0f + 1.0f) / 1000.0f;
		rPrime3 = r_point * r_point * r_point;

		qw_r3_N = device__qgaus_W[i] * rPrime3 * coeff * exp( -((g - gPrime) * (g - gPrime) / (2 * f_stdev * f_stdev)) );

		xyz2 = r_point * sinb;
		zp = r_point * cosb;
		xyz0 = zp * cosl - f_lbr_r;
		xyz1 = zp * sinl;

		rg = sqrt(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2)/(q*q));
		rs = rg + r0;

		bg_int += qw_r3_N / (rg * rs * rs * rs);

		for (int j = 0; j < number_streams; j++) {
			pos = (j * 3);
			sxyz0 = xyz0 - device__fstream_c[pos];
			sxyz1 = xyz1 - device__fstream_c[pos + 1];
			sxyz2 = xyz2 - device__fstream_c[pos + 2];

			dotted = device__fstream_a[pos] * sxyz0 + device__fstream_a[pos + 1] * sxyz1 + device__fstream_a[pos + 2] * sxyz2;

			sxyz0 -= dotted * device__fstream_a[pos];
			sxyz1 -= dotted * device__fstream_a[pos + 1];
			sxyz2 -= dotted * device__fstream_a[pos + 2];

			st_int[j] += qw_r3_N * exp(-((sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2)) / device__fstream_sigma_sq2[j]);
		}
	}
	float probability_sum = 0.0;
	probability_sum += bg_int * reff_xr_rp3 * device__background_weight[0];
	for (i = 0; i < number_streams; i++) {
		probability_sum += st_int[i] * reff_xr_rp3 * device__stream_weight[i];
	}
//	printf("bg_prob %.15f st_prob[0]: %.15f st_prob[1]: %.15f, prob_sum: %.15f\n", (bg_int * reff_xr_rp3), (st_int[0] * reff_xr_rp3), (st_int[1] * reff_xr_rp3), probability_sum);

	if (probability_sum == 0.0) probability[threadIdx.x] += -238.0;
	else probability[threadIdx.x] += log(probability_sum)/log(10.0);
}

/********
 *	Run the GPU kernels and get the probability
 ********/

#define stream_parameters(x, y) parameters[(x * 6) + y + 3]
#define stream_weights(x) parameters[(x * 6) + 2]
//#define background_weight parameters[0]
#define background_weight 0.0
//#define alpha parameters[1]
#define q parameters[0]
#define r0 parameters[1]
//#define delta parameters[4]

double gpu__likelihood(double *parameters) {
	int i, j;

	double stream_c[3], lbr[3];
	float fstream_a[number_streams * 3], fstream_c[number_streams * 3], fstream_sigma_sq2[number_streams];

	for (i = 0; i < number_streams; i++) {
		fstream_sigma_sq2[i] = (float)(2.0 * stream_parameters(i,4) * stream_parameters(i,4));

		fstream_a[(i * 3)] = (float)( sin(stream_parameters(i,2)) * cos(stream_parameters(i,3)) );
		fstream_a[(i * 3) + 1] = (float)( sin(stream_parameters(i,2)) * sin(stream_parameters(i,3)) );
		fstream_a[(i * 3) + 2] = (float)( cos(stream_parameters(i,2)) );

		gc_to_gal(wedge, stream_parameters(i,0) * D_DEG2RAD, 0 * D_DEG2RAD, &(lbr[0]), &(lbr[1]));
		lbr[2] = stream_parameters(i,1);
		d_lbr2xyz(lbr, stream_c);

		fstream_c[(i * 3)] = (float)stream_c[0]; 
		fstream_c[(i * 3) + 1] = (float)stream_c[1];
		fstream_c[(i * 3) + 2] = (float)stream_c[2];
	}

	cutilSafeCall( cudaMemcpyToSymbol(device__fstream_sigma_sq2, fstream_sigma_sq2, number_streams * sizeof(float), 0, cudaMemcpyHostToDevice) ); 
	cutilSafeCall( cudaMemcpyToSymbol(device__fstream_a, fstream_a, number_streams * 3 * sizeof(float), 0, cudaMemcpyHostToDevice) );
	cutilSafeCall( cudaMemcpyToSymbol(device__fstream_c, fstream_c, number_streams * 3 * sizeof(float), 0, cudaMemcpyHostToDevice) );

	double background_integral = 0.0;
	double *stream_integrals = (double*)malloc(number_streams * sizeof(double));
	for (i = 0; i < number_streams; i++) stream_integrals[i] = 0.0;

	double coeff = 1.0 / (d_stdev * sqrt(2.0 * D_PI));

	for (i = 0; i < number_integrals; i++) {
		dim3 dimGrid(r_steps[i], mu_increment);

		gpu__zero_integrals<2><<<dimGrid, nu_steps[i]>>>(device__background_integrals[i], device__stream_integrals[i]);
		printf("zeroed integrals\n");
		for (j = 0; j < mu_steps[i]; j += mu_increment) {
/*			gpu__integral_kernel2<2><<<dimGrid, nu_steps[i]>>>(	convolve, j, mu_steps[i],
										r_min[i], r_step_size[i],
										q, r0,
										(float)coeff,
										device__lb[i], device__V[i],
										device__background_integrals[i], device__stream_integrals[i]);
*/			switch(number_streams) {
				case 1:	gpu__integral_kernel<1><<<dimGrid, nu_steps[i]>>>(	convolve, j, mu_steps[i], 
							q, r0,
							device__lb[i], device__r_constants[i], device__V[i], 
							device__background_integrals[i],
							device__stream_integrals[i]);
					break;
				case 2:	gpu__integral_kernel<2><<<dimGrid, nu_steps[i]>>>(	convolve, j, mu_steps[i], 
							q, r0,
							device__lb[i], device__r_constants[i], device__V[i], 
							device__background_integrals[i],
							device__stream_integrals[i]);
					break;
				case 3:	gpu__integral_kernel<3><<<dimGrid, nu_steps[i]>>>(	convolve, j, mu_steps[i], 
							q, r0,
							device__lb[i], device__r_constants[i], device__V[i], 
							device__background_integrals[i],
							device__stream_integrals[i]);
					break;
				case 4:	gpu__integral_kernel<4><<<dimGrid, nu_steps[i]>>>(	convolve, j, mu_steps[i], 
							q, r0,
							device__lb[i], device__r_constants[i], device__V[i], 
							device__background_integrals[i],
							device__stream_integrals[i]);
					break;
			}
		}
		printf("completed integral kernal evaluations\n");
		cpu__sum_integrals(i, &background_integral, stream_integrals);
		printf("bg_integral[%d]   : %.15lf\n", i, background_integral);
		printf("st_integral[%d][0]: %.15lf\n", i, stream_integrals[0]);
		printf("st_integral[%d][1]: %.15lf\n", i, stream_integrals[1]);
		printf("\n");
	}

	int block_size;

	double *stream_weight = (double*)malloc(number_streams * sizeof(double));
	double exp_weight = exp(background_weight);
	double sum_exp_weights = exp_weight; 
	double bg_weight = exp_weight/background_integral;
	for (i = 0; i < number_streams; i++) {
		exp_weight = exp(stream_weights(i));
		sum_exp_weights += exp_weight;
		stream_weight[i] = exp_weight/stream_integrals[i];
	}

	float f_background_weight[1];
	float *f_stream_weight = (float*)malloc(number_streams * sizeof(float));
	f_background_weight[0] = (float)( bg_weight / sum_exp_weights );
	for (i = 0; i < number_streams; i++) {
		f_stream_weight[i] = (float)( stream_weight[i] / sum_exp_weights );
	}

	cutilSafeCall( cudaMemcpyToSymbol(device__background_weight, f_background_weight, 1 * sizeof(float), 0, cudaMemcpyHostToDevice) );
	cutilSafeCall( cudaMemcpyToSymbol(device__stream_weight, f_stream_weight, number_streams * sizeof(float), 0, cudaMemcpyHostToDevice) );

	double likelihood = 0.0;
	gpu__zero_likelihood<2><<<1, number_threads>>>(number_threads, device__probability);
	for (i = 0; i < number_stars; i += number_threads) {
		block_size = min(number_threads, number_stars - i);
		switch (number_streams) {
			case 1:	gpu__likelihood_kernel<1><<<1, block_size>>>(	i, convolve,
										q, r0,
										(float)coeff,
										device__stars,
										device__probability);
			break;
			case 2:	gpu__likelihood_kernel<2><<<1, block_size>>>(	i, convolve,
										q, r0,
										(float)coeff,
										device__stars,
										device__probability);
			break;
			case 3:	gpu__likelihood_kernel<3><<<1, block_size>>>(	i, convolve,
										q, r0,
										(float)coeff,
										device__stars,
										device__probability);
			break;
			case 4:	gpu__likelihood_kernel<4><<<1, block_size>>>(	i, convolve,
										q, r0,
										(float)coeff,
										device__stars,
										device__probability);
			break;
		}
	}
	cpu__sum_likelihood(number_threads, &likelihood);
	likelihood /= number_stars;
	printf("likelihood: %.15lf\n", likelihood);
	return likelihood;
}
