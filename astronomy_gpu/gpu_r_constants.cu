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

#ifndef GPU__R_CONSTANTS
#define GPU__R_CONSTANTS

extern "C" {
#include "pi_constants.h"
#include "r_constants.h"
#include "gauss_legendre.h"
#include <stdio.h>
}


__global__ void gpu__reff_gPrime(float r_min, float r_step_size, float mu_step_size, float *gPrime, float *reff_xr_rp3_irv) {
        int r_step;
        float log_r, r, next_r, coords, gp, exp_result, reff_value, rPrime3, irv;

	r_step = threadIdx.x;

	log_r = r_min + (r_step * r_step_size);
	r = pow(10.0, (log_r - 14.2)/5.0);
	next_r = pow(10.0, (log_r - 14.2 + r_step_size)/5.0);
	coords = (next_r + r)/2.0;
	gp = 5.0 * (log10(coords * 1000) - 1.0) + f_absm;

	irv = (((next_r * next_r * next_r) - (r * r * r))/3.0) * mu_step_size * D_DEG2RAD;

	gPrime[r_step] = gp;

	exp_result = exp(sigmoid_curve_1 * (gp - sigmoid_curve_2));
	reff_value = sigmoid_curve_0 / (exp_result + 1);
	rPrime3 = coords * coords * coords;
	reff_xr_rp3_irv[r_step] = irv * reff_value * f_xr / rPrime3;
}

__global__ void gpu__r_qw(float coeff, float *dx, float *qgaus_W, float *gPrime, float *r_point, float *qw_r3_N) {
	int r_step, convolve_step, n_convolve, position;
	float g, rp, r3, exponent, N;

	r_step = threadIdx.x;
	convolve_step = blockIdx.x;
	n_convolve  = gridDim.x;

	position = (r_step * n_convolve) + convolve_step;

	g = gPrime[r_step] + dx[convolve_step];
	rp = pow(10.0, (g - f_absm)/5.0 + 1.0) / 1000.0;

	r_point[position] = rp;

	r3 = rp * rp * rp;
	exponent = (dx[convolve_step] * dx[convolve_step]) / (2 * f_stdev * f_stdev);
	N = coeff * exp(-exponent);

	qw_r3_N[position] = qgaus_W[convolve_step] * r3 * N;
}


__global__ void gpu__reff_V(float nu_min, float nu_step_size, float *reff_xr_rp3_irv, float *V) {
	int nu_step, r_step, r_steps;
	float nu, ids;

	r_step = threadIdx.x;
	r_steps = blockDim.x;
	nu_step = threadIdx.y;

	nu = nu_min + (nu_step * nu_step_size);
	ids = cos((90 - nu - nu_step_size) * D_DEG2RAD) - cos((90 - nu) * D_DEG2RAD);

	V[(nu_step * r_steps) + r_step] = reff_xr_rp3_irv[r_step] * ids;
}

void gpu__r_constants(int n_convolve, INTEGRAL *integral, float **gpu__V, float **gpu__r_point, float **gpu__qw_r3_N) {
	int i;
	float *device__gPrime, *device__reff_xr_rp3_irv;

	cudaMalloc((void**) &device__gPrime, integral->r_steps * sizeof(float));
	cudaMalloc((void**) &device__reff_xr_rp3_irv, integral->r_steps * sizeof(float));

	gpu__reff_gPrime<<<1,integral->r_steps>>>(integral->r_min, integral->r_step_size, integral->mu_step_size, device__gPrime, device__reff_xr_rp3_irv);
	printf("reff_gPrime\n");

	cudaMalloc((void**) gpu__V, integral->nu_steps * integral->r_steps * sizeof(float));

	dim3 dimBlock(integral->r_steps, integral->nu_steps);
	gpu__reff_V<<<1, dimBlock>>>(integral->nu_min, integral->nu_step_size, device__reff_xr_rp3_irv, *gpu__V);
	printf("reff_V\n");

	float *host__qgaus_W = (float*)malloc(n_convolve * sizeof(float));
	float *host__qgaus_X = (float*)malloc(n_convolve * sizeof(float));
	float *host__dx = (float*)malloc(n_convolve * sizeof(float));
	printf("did qgaus mallocs\n");

	f_gauss_legendre(-1.0, 1.0, host__qgaus_X, host__qgaus_W, n_convolve);
	printf("did f_gauss_legendre\n");

	for (i = 0; i < n_convolve; i++) {
		host__dx[i] = 3 * d_stdev * host__qgaus_X[i];
	}

	printf("calculated inputs\n");

	float *device__qgaus_W, *device__dx;
	cudaMalloc((void**) &device__qgaus_W, n_convolve * sizeof(float));
	cudaMalloc((void**) &device__dx, n_convolve * sizeof(float));
	printf("did malloc\n");
	cudaMemcpy(device__qgaus_W, host__qgaus_W, n_convolve * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(device__dx, host__dx, n_convolve * sizeof(float), cudaMemcpyHostToDevice);
	printf("did memcpy\n");

	double coeff = 1.0 / (d_stdev * sqrt(2.0 * F_PI));

	printf("did coeff\n");
	cudaMalloc((void**) gpu__r_point, integral->r_steps * n_convolve * sizeof(float));
	cudaMalloc((void**) gpu__qw_r3_N, integral->r_steps * n_convolve * sizeof(float));
	printf("did last mallocs\n");

	gpu__r_qw<<<n_convolve, integral->r_steps>>>((float)coeff, device__dx, device__qgaus_W, device__gPrime, *gpu__r_point, *gpu__qw_r3_N);
	printf("r_qw\n");

	cudaFree(device__reff_xr_rp3_irv);
	cudaFree(device__gPrime);
	cudaFree(device__qgaus_W);
	cudaFree(device__dx);
	free(host__qgaus_X);
	free(host__qgaus_W);
	free(host__dx);
}

void gpu__copy_r_constants(int n_convolve, INTEGRAL *integral, float *device__V, float *device__r_point, float *device__qw_r3_N, float **host__V, float **host__r_point, float **host__qw_r3_N) {
	*host__V = (float*)malloc(integral->nu_steps * integral->r_steps * sizeof(float));
	*host__r_point = (float*)malloc(integral->r_steps * n_convolve * sizeof(float));
	*host__qw_r3_N = (float*)malloc(integral->r_steps * n_convolve * sizeof(float));

	cudaMemcpy(*host__V, device__V, integral->nu_steps * integral->r_steps * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(*host__r_point, device__r_point, integral->r_steps * n_convolve * sizeof(float), cudaMemcpyDeviceToHost); 
	cudaMemcpy(*host__qw_r3_N, device__qw_r3_N, integral->r_steps * n_convolve * sizeof(float), cudaMemcpyDeviceToHost); 
}

#endif
