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

#include <stdlib.h>
#include <math.h>

#include "pi_constants.h"
#include "coords.h"

#include "../astronomy/parameters.h"
#include "../astronomy/atSurveyGeometry.h"


void gc_eq_gal(int wedge, double amu_rad, double anu_rad, double *cpu__glong, double *cpu__glat) {
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
	*cpu__glong = fmod(atan2(v1, v0), D_2PI);
	*cpu__glat = fmod(atan2(v2, r), D_2PI);
}


void gc_eq_gal_lb(int wedge, double amu_rad, double anu_rad, double *cpu__lb) {
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
	double glong = dmod(atan2(v1, v0), D_2PI);
	double glat = dmod(atan2(v2, r), D_2PI);

	cpu__lb[0] = sin(glat);
	cpu__lb[1] = sin(glong);
	cpu__lb[2] = cos(glat);
	cpu__lb[3] = cos(glong);
}

void cpu__gc_eq_gal_lb(	int wedge,
			int mu_steps, double mu_min, double mu_step_size,
			int nu_steps, double nu_min, double nu_step_size,
			double **cpu__lb) {
	int i, j;
        double mu_min_rad = mu_min * D_DEG2RAD;
        double mu_step_rad = mu_step_size * D_DEG2RAD;
        double nu_min_rad = nu_min * D_DEG2RAD;
        double nu_step_rad = nu_step_size * D_DEG2RAD;

	*cpu__lb = (double*)malloc(4 * mu_steps * nu_steps * sizeof(double));

	double anu, amu;
	int pos;
	for (i = 0; i < mu_steps; i++) {
		amu = mu_min_rad + ((i + 0.5) * mu_step_rad);
 		for (j = 0; j < nu_steps; j++) {
			anu = nu_min_rad + ((j + 0.5) * nu_step_rad);
       			pos = ((i * nu_steps) + j) * 4;
			gc_eq_gal_lb(wedge, amu, anu, &( (*cpu__lb)[pos] ));
		}
	}
}


void gc_sgr_gal(int wedge, double amu_rad, double anu_rad, double *cpu__glong, double *cpu__glat) {
//	gcToSgr(ap->stream_parameters[i][0], 0, wedge, &lamda, &beta);
	double lamda, beta;

	double x = cos(amu_rad)*cos(anu_rad);
	double y = -sin(anu_rad);
	double z = sin(amu_rad)*cos(anu_rad);

	lamda = atan2(y, x);
	lamda = lamda * D_RAD2DEG; 
	lamda = lamda + 2.5 * wedge;
	if (lamda < 0) lamda = lamda + 360;

	beta = asin(z);
	beta *= D_RAD2DEG;

	sgrToGal(lamda, beta, cpu__glong, cpu__glat); 
	(*cpu__glong) *= D_DEG2RAD;
	(*cpu__glat) *= D_DEG2RAD;
}

void gc_sgr_gal_lb(int wedge, double amu_rad, double anu_rad, double *cpu__lb) {
	double glong, glat;
	gc_sgr_gal(wedge, amu_rad, anu_rad, &glong, &glat); 

	cpu__lb[0] = sin(glat);
	cpu__lb[1] = sin(glong);
	cpu__lb[2] = cos(glat);
	cpu__lb[3] = cos(glong);
}

void cpu__gc_sgr_gal_lb(	int wedge,
				int mu_steps, double mu_min, double mu_step_size,
				int nu_steps, double nu_min, double nu_step_size,
				double **cpu__lb) {
	int i, j;
        double mu_min_rad = mu_min * D_DEG2RAD;
        double mu_step_rad = mu_step_size * D_DEG2RAD;
        double nu_min_rad = nu_min * D_DEG2RAD;
        double nu_step_rad = nu_step_size * D_DEG2RAD;

	*cpu__lb = (double*)malloc(4 * mu_steps * nu_steps * sizeof(double));

	double anu, amu;
	int pos;
	for (i = 0; i < mu_steps; i++) {
		amu = mu_min_rad + ((i + 0.5) * mu_step_rad);
		for (j = 0; j < nu_steps; j++) {
			anu = nu_min_rad + ((j + 0.5) * nu_step_rad);
       			pos = ((i * nu_steps) + j) * 4;
			gc_sgr_gal_lb(wedge, amu, anu, &( (*cpu__lb)[pos] ));
		}
	}
}

void cpu__gc_to_lb(int wedge, INTEGRAL *integral, double **cpu__lb) {
	cpu__gc_eq_gal_lb(wedge, integral->mu_steps, integral->mu_min, integral->mu_step_size, integral->nu_steps, integral->nu_min, integral->nu_step_size, cpu__lb);
}



void populate_cpu__lb( int sgr_coordinates, int wedge,
		       int mu_steps, double mu_min, double mu_step_size,
		       int nu_steps, double nu_min, double nu_step_size,
		       double **cpu__lb)
{
  	int i, j;
	*cpu__lb = (double*)malloc(4 * mu_steps * nu_steps * sizeof(double));

	int pos;
	for (i = 0; i < mu_steps; i++) {
	  double amu_deg = mu_min + ((i + 0.5) * mu_step_size);
		for (j = 0; j < nu_steps; j++) {
		  double anu_deg = (nu_min + (j + 0.5) * nu_step_size);
			pos = ((i * nu_steps) + j) * 4;
			//gc_eq_gal_lb(wedge, amu, anu, &( (*cpu__lb)[pos] ));
			//use the same code as the CPU app to determine cpu__lb
			double glong, glat;
			if (sgr_coordinates == 0) {
			  double ra, dec;
			  atGCToEq(amu_deg, anu_deg, &ra, &dec, 95.0, d_get_incl(wedge));
			  atEqToGal(ra, dec, &glong, &glat);
			} else {
			  double lamda, beta;
			  gcToSgr(amu_deg, anu_deg, wedge, &lamda, &beta);
			  sgrToGal(lamda, beta, &glong, &glat);
			}
			(&(*cpu__lb)[pos])[0] = sin(glat * D_DEG2RAD);
			(&(*cpu__lb)[pos])[1] = sin(glong * D_DEG2RAD);
			(&(*cpu__lb)[pos])[2] = cos(glat * D_DEG2RAD);
			(&(*cpu__lb)[pos])[3] = cos(glong * D_DEG2RAD);
		}
	}
}
