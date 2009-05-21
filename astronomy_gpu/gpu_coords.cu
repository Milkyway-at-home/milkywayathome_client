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

#ifndef GPU__COORDS
#define GPU__COORDS

extern "C" {
#include <stdio.h>
#include <math.h>
#include "pi_constants.h"
#include "coords.h"

#include "../astronomy/parameters.h"
}


__global__ void calculate_gc_to_gal(float a_inc_rad, float mu_min_rad, float mu_step_rad, float nu_min_rad, float nu_step_rad, float *glong_out, float *glat_out) {
	int nu = threadIdx.x;
	int nu_steps = blockDim.x;
	int mu = blockIdx.x;
	int position = mu * nu_steps + nu;

	float x1, y1, z1, x2, y2, z2;

	float amu = (mu_min_rad + ((mu + 0.5) * mu_step_rad));
	float anu = (nu_min_rad + ((nu + 0.5) * nu_step_rad));
	float cos_ainc = cos(a_inc_rad);
	float sin_ainc = sin(a_inc_rad);
	x2 = cos(amu - F_A_NODE_RAD) * cos(anu);
	y2 = sin(amu - F_A_NODE_RAD) * cos(anu);
	z2 = sin(anu);
	x1 = x2;
	y1 = (y2 * cos_ainc) - (z2 * sin_ainc);
	z1 = (y2 * sin_ainc) + (z2 * cos_ainc);

	float ra = atan2(y1, x1) + F_A_NODE_RAD;	//a_node = 95.0
	float dec = asin(z1);

//	printf("amu: %.15f, anu: %.15f, anode: %.15f, ainc: %.15f, x2: %.15f, y2: %.15f, z2: %.15f, y1: %.15f, z1: %.15f\n", amu, anu, F_A_NODE_RAD, ainc, x2, y2, z2, y1, z1);

//`	printf("ra[%d][%d]: %.15f, dec[%d][%d]: %.15f\n", mu, nu, ra, mu, nu, dec);
	
//	ra /= F_DEG2RAD;
//	dec /= F_DEG2RAD;

	//slaDcs2c(ra, dec, v1)
	float cosb = cos(dec);
	float v0 = cos(ra) * cosb;
	float v1 = sin(ra) * cosb;
	float v2 = sin(dec);

	//slaDmxv(rmat, v1, v2)
	float v0_temp = (rmat00 * v0) + (rmat01 * v1) + (rmat02 * v2);
	float v1_temp = (rmat10 * v0) + (rmat11 * v1) + (rmat12 * v2);
	float v2_temp = (rmat20 * v0) + (rmat21 * v1) + (rmat22 * v2);

	v0 = v0_temp;
	v1 = v1_temp;
	v2 = v2_temp;

	//slaDcc2s(v2, dl, db)
	float r = sqrt( (v0 * v0) + (v1 * v1) );
	/*
	 *	Might need to check if r is 0, then glong and glat are both 0
	**/
	float glong = atan2(v1, v0); 
	float glat = atan2(v2, r);

	//glong = slaDranrm( glong )
	//	Might need to check if glong within 0 .. 2PI
	glong = fmod(glong, F_2PI);
	//glat = slaDrange( glat)
	//	Might need to check if glat within 0 .. 2PI
	glat = fmod(glat, F_2PI);

	glong_out[position] = glong;
	glat_out[position] = glat;
}

void gpu__gc_to_gal(int wedge, INTEGRAL *integral, float **device__glong, float **device__glat) {
	double ainc_rad = d_get_incl_rad(wedge);
	double mu_min_rad = integral->mu_min * D_DEG2RAD;
	double mu_step_rad = integral->mu_step_size * D_DEG2RAD;
	double nu_min_rad = integral->nu_min * D_DEG2RAD;
	double nu_step_rad = integral->nu_step_size * D_DEG2RAD;

        cudaMalloc((void**) device__glong, integral->mu_steps * integral->nu_steps * sizeof(float));
        cudaMalloc((void**) device__glat, integral->mu_steps * integral->nu_steps * sizeof(float));

	calculate_gc_to_gal<<<integral->mu_steps, integral->nu_steps>>>((float)ainc_rad, (float)mu_min_rad, (float)mu_step_rad, (float)nu_min_rad, (float)nu_step_rad, *device__glong, *device__glat);
}

void gpu__copy_gal_to_host(INTEGRAL *integral, float *device__glong, float *device__glat, float **host__glong, float **host__glat) {
	int size = integral->mu_steps * integral->nu_steps;
	*host__glong = (float*)malloc(size * sizeof(float));
	*host__glat = (float*)malloc(size * sizeof(float));

	cudaMemcpy(*host__glong, device__glong, size * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(*host__glat, device__glat, size * sizeof(float), cudaMemcpyDeviceToHost);
}

#endif
