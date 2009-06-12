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

#include "../astronomy/parameters.h"

#include "r_constants.h"
#include "pi_constants.h"
#include "coords.h"
#include "cpu_integrals.h"
#include "cpu_coords.h"

void cpu__integral_point(int convolve, double alpha, double q, double r0, double delta, int number_streams, double *stream_a, double *stream_c, double *stream_sigma_sq2, double V, double *r_point, double *qw_r3_N, double glong, double glat, double *background_integral, double *stream_integrals) {
	int i, j, pos;
	double sinb = sin(glat);
	double sinl = sin(glong);
	double cosb = cos(glat);
	double cosl = cos(glong);

	for (i = 0; i < number_streams; i++) stream_integrals[i] = 0.0;

	double xyz0, xyz1, xyz2;
	double sxyz0, sxyz1, sxyz2;
	double zp, rg, rs, dotted, xyz_norm;

	*background_integral = 0.0;
	for (i = 0; i < convolve; i++) {
		xyz2 = r_point[i] * sinb;
		zp = r_point[i] * cosb;
		xyz0 = zp * cosl - d_lbr_r;
		xyz1 = zp * sinl;

		rg = sqrt(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2)/(q*q));

		if (alpha == delta == 1) {
			rs = rg + r0;
			*background_integral += qw_r3_N[i] / (rg * rs * rs * rs);
		} else {
			*background_integral += qw_r3_N[i] / (pow(rg, alpha) * pow(rg + r0, 3 - alpha + delta));
		}

		for (j = 0; j < number_streams; j++) {
			pos = (j * 3);

			sxyz0 = xyz0 - stream_c[pos];
			sxyz1 = xyz1 - stream_c[pos + 1];
			sxyz2 = xyz2 - stream_c[pos + 2];

			dotted = stream_a[pos] * sxyz0 + stream_a[pos + 1] * sxyz1 + stream_a[pos + 2] * sxyz2;

			sxyz0 -= dotted * stream_a[pos];
			sxyz1 -= dotted * stream_a[pos + 1];
			sxyz2 -= dotted * stream_a[pos + 2];

			xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);

			stream_integrals[j] += qw_r3_N[i] * exp(-xyz_norm / stream_sigma_sq2[j]);
		}
	}

	(*background_integral) *= V;
	for (i = 0; i < number_streams; i++) stream_integrals[i] *= V;
}


void cpu__integrals(ASTRONOMY_PARAMETERS *ap, INTEGRAL *integral, double *V, double *r_point, double *qw_r3_N, double *glong, double *glat, double **background_integrals, double **stream_integrals) {
	int i, j, pos;
	int k;

	double alpha = ap->background_parameters[0];
	double q = ap->background_parameters[1];
	double r0 = ap->background_parameters[2];
	double delta = ap->background_parameters[3];

	double *stream_a = (double*)malloc(3 * ap->number_streams * sizeof(double));
	double *stream_c = (double*)malloc(3 * ap->number_streams * sizeof(double));
	double *stream_sigma_sq2 = (double*)malloc(ap->number_streams * sizeof(double));

	double lbr[3];
	for (i = 0; i < ap->number_streams; i++) {
		stream_sigma_sq2[i] = 2.0 * ap->stream_parameters[i][4] * ap->stream_parameters[i][4];

		stream_a[(i * 3)]	= sin(ap->stream_parameters[i][2]) * cos(ap->stream_parameters[i][3]);
		stream_a[(i * 3) + 1]	= sin(ap->stream_parameters[i][2]) * sin(ap->stream_parameters[i][3]);
		stream_a[(i * 3) + 2]	= cos(ap->stream_parameters[i][2]);

		if (ap->sgr_coordinates == 0) {
			gc_eq_gal(ap->wedge, ap->stream_parameters[i][0] * D_DEG2RAD, 0 * D_DEG2RAD, &(lbr[0]), &(lbr[1]));
		} else {
			gc_sgr_gal(ap->wedge, ap->stream_parameters[i][0] * D_DEG2RAD, 0 * D_DEG2RAD, &(lbr[0]), &(lbr[1]));
		}
		lbr[2] = ap->stream_parameters[i][1];
		d_lbr2xyz(lbr, &(stream_c[(i * 3)]));
	}

	//check if q = 0

	(*background_integrals) = (double*)malloc(integral->r_steps * integral->mu_steps * integral->nu_steps * sizeof(double));
	(*stream_integrals) = (double*)malloc(integral->r_steps * integral->mu_steps * integral->nu_steps * ap->number_streams * sizeof(double));

	int nu_step = 0;
	double *r_point__in, *qw_r3_N__in, *bg_integral__in, *st_integral__in;
	for (j = 0; j < integral->mu_steps * integral->nu_steps; j++) {
		nu_step = j % integral->nu_steps;
		for (i = 0; i < integral->r_steps; i++) {
			r_point__in = &(r_point[i * ap->convolve]);
			qw_r3_N__in = &(qw_r3_N[i * ap->convolve]);

			pos = i + (j * integral->r_steps);
			bg_integral__in = &((*background_integrals)[pos]);
			st_integral__in = &((*stream_integrals)[pos * ap->number_streams]);

			cpu__integral_point(ap->convolve, alpha, q, r0, delta, ap->number_streams, stream_a, stream_c, stream_sigma_sq2, V[(nu_step * integral->r_steps) + i], r_point__in, qw_r3_N__in, glong[j], glat[j], bg_integral__in, st_integral__in);
		}
	}
}


