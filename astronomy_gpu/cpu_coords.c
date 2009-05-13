#include <stdlib.h>
#include <math.h>

#include "pi_constants.h"
#include "coords.h"

#include "../astronomy/parameters.h"


void gc_to_gal(int wedge, double amu_rad, double anu_rad, double *cpu__sinb, double *cpu__sinl, double *cpu__cosb, double *cpu__cosl) {
        double a_inc_rad = d_get_incl_rad(wedge);
	double cos_ainc = cos(a_inc_rad);
	double sin_ainc = sin(a_inc_rad);
	double x2 = cos(amu_rad - D_A_NODE_RAD) * cos(anu_rad);
	double y2 = sin(amu_rad - D_A_NODE_RAD) * cos(anu_rad);
	double z2 = sin(anu_rad);
	double x1 = x2;
	double y1 = (y2 * cos_ainc) - (z2 * sin_ainc);
	double z1 = (y2 * sin_ainc) + (z2 * cos_ainc);

	double ra = atan2(y1, x1) + D_A_NODE_RAD;
	double dec = asin(z1);

	double cosb = cos(dec);
	double v0 = cos(ra) * cosb;
	double v1 = sin(ra) * cosb;
	double v2 = sin(dec);

	double v0_temp = (rmat00 * v0) + (rmat01 * v1) + (rmat02 * v2);
	double v1_temp = (rmat10 * v0) + (rmat11 * v1) + (rmat12 * v2);
	double v2_temp = (rmat20 * v0) + (rmat21 * v1) + (rmat22 * v2);

	v0 = v0_temp;
	v1 = v1_temp;
	v2 = v2_temp;

	//slaDcc2s(v2, dl, db)
	double r = sqrt( (v0 * v0) + (v1 * v1) );
	/*
	 *      Might need to check if r is 0, then glong and glat are both 0
	**/
	double glong = fmod(atan2(v1, v0), D_2PI);
	double glat = fmod(atan2(v2, r), D_2PI);

	(*cpu__sinb) = sin(glat);
	(*cpu__sinl) = sin(glong);
	(*cpu__cosb) = cos(glat);
	(*cpu__cosl) = cos(glong);
}

void cpu__gc_to_gal(int wedge, INTEGRAL *integral, double **cpu__sinb, double **cpu__sinl, double **cpu__cosb, double **cpu__cosl) {
	int i, j;
        double mu_min_rad = integral->mu_min * D_DEG2RAD;
        double mu_step_rad = integral->mu_step_size * D_DEG2RAD;
        double nu_min_rad = integral->nu_min * D_DEG2RAD;
        double nu_step_rad = integral->nu_step_size * D_DEG2RAD;

	*cpu__sinb = (double*)malloc(integral->mu_steps * integral->nu_steps * sizeof(double));
	*cpu__sinl = (double*)malloc(integral->mu_steps * integral->nu_steps * sizeof(double));
	*cpu__cosb = (double*)malloc(integral->mu_steps * integral->nu_steps * sizeof(double));
	*cpu__cosl = (double*)malloc(integral->mu_steps * integral->nu_steps * sizeof(double));

	double anu, amu;
	int pos;
	for (i = 0; i < integral->mu_steps; i++) {
		amu = mu_min_rad + ((i + 0.5) * mu_step_rad);
		for (j = 0; j < integral->nu_steps; j++) {
			anu = nu_min_rad + ((j + 0.5) * nu_step_rad);
       			pos = (i * integral->nu_steps) + j;
			gc_to_gal(wedge, amu, anu, &( (*cpu__sinb)[pos] ), &( (*cpu__sinl)[pos] ), &( (*cpu__cosb)[pos] ), &( (*cpu__cosl)[pos] ));
		}
	}
}
