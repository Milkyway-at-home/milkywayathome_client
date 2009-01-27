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

/****
        *       BOINC includes
*****/

#ifdef GMLE_BOINC
	#ifdef _WIN32
		#include "boinc_win.h"
	#else
		#include "config.h"
	#endif

	#ifndef _WIN32
		#include <cstdio>
		#include <cctype>
		#include <cstring>
		#include <cstdlib>
		#include <csignal>
	#endif

	#ifdef BOINC_APP_GRAPHICS
		#include "graphics_api.h"
		#include "graphics_lib.h"
	#endif

	#include "diagnostics.h"
	#include "util.h"
	#include "filesys.h"
	#include "boinc_api.h"
	#include "mfile.h"
#endif

#define CHECKPOINT_FILE "astronomy_checkpoint"

/****
	*	Astronomy includes
*****/
#include <math.h>
#include <time.h>

#include "evaluation_optimized.h"
#include "parameters.h"
#include "probability.h"
#include "stCoords.h"
#include "atSurveyGeometry.h"
#include "stVector.h"
#include "star_points.h"
#include "numericalIntegration.h"

#ifndef _WIN32
	#define pi M_PI
#else
	#define pi 3.1415926535897932384626433832795028841971693993751
#endif

#define NEW_FORMULA 0
#define WEDGE_ALLOW_ZERO 0

#define deg (180.0/pi)
#define stdev 0.6
#define xr 3 * stdev
#define lbr_r 8.5
#define absm 4.2

double sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };

double alpha, q, r0, delta, coeff, alpha_delta3;
double *qgaus_X, *qgaus_W, **xyz, *dx;
double **stream_a, **stream_c, *stream_sigma;

void init_constants(ASTRONOMY_PARAMETERS *ap) {
	if (ap->convolve > 0) {
		int i;

		qgaus_X		= (double*)malloc(sizeof(double) * ap->convolve);
		qgaus_W		= (double*)malloc(sizeof(double) * ap->convolve);
		dx		= (double*)malloc(sizeof(double) * ap->convolve);
		xyz		= (double**)malloc(sizeof(double*) * ap->convolve);
		stream_sigma	= (double*)malloc(sizeof(double) * ap->number_streams);
		stream_a	= (double**)malloc(sizeof(double*) * ap->number_streams);	
		stream_c	= (double**)malloc(sizeof(double*) * ap->number_streams);	

		alpha	= ap->background_parameters[0];
		q	= ap->background_parameters[1];
		r0	= ap->background_parameters[2];
		delta	= ap->background_parameters[3];
		alpha_delta3 = 3 - alpha + delta;
		coeff	= 1 / (stdev * sqrt(2*pi));

		gaussLegendre(-1.0, 1.0, qgaus_X, qgaus_W, ap->convolve);

		for (i = 0; i < ap->convolve; i++) {
			xyz[i] = (double*)malloc(sizeof(double) * 3);
			dx[i] = 3 * stdev * qgaus_X[i];
		}

		for (i = 0; i < ap->number_streams; i++) {
			double ra, dec, lamda, beta, l, b, lbr[3];

			stream_a[i] = (double*)malloc(sizeof(double) * 3);
			stream_c[i] = (double*)malloc(sizeof(double) * 3);
			stream_sigma[i] = ap->stream_parameters[i][4];
   
			if (ap->sgr_coordinates == 0) {
				atGCToEq(ap->stream_parameters[i][0], 0, &ra, &dec, get_node(), wedge_incl(ap->wedge));
				atEqToGal(ra, dec, &l, &b);
			} else if (ap->sgr_coordinates == 1) {
				gcToSgr(ap->stream_parameters[i][0], 0, ap->wedge, &lamda, &beta); //vickej2
				sgrToGal(lamda, beta, &l, &b); //vickej2
			} else {
				printf("Error: sgr_coordinates not valid");
			}

			lbr[0] = l;
			lbr[1] = b;
			lbr[2] = ap->stream_parameters[i][1];
			lbr2xyz(lbr, stream_c[i]);
    
			stream_a[i][0] = sin(ap->stream_parameters[i][2]) * cos(ap->stream_parameters[i][3]);
			stream_a[i][1] = sin(ap->stream_parameters[i][2]) * sin(ap->stream_parameters[i][3]);
			stream_a[i][2] = cos(ap->stream_parameters[i][2]);
		}
	}
}

void free_constants(ASTRONOMY_PARAMETERS *ap) {
	if (ap->convolve > 0) {
		int i;

		free(qgaus_X);
		free(qgaus_W);
		free(dx);
		free(stream_sigma);
		for (i = 0; i < ap->number_streams; i++) {
			free(stream_a[i]);
			free(stream_c[i]);
			free(xyz[i]);
		}
		free(stream_a);
		free(stream_c);
		free(xyz);
	}
}

void set_probability_constants(ASTRONOMY_PARAMETERS *ap, double coords, double *r_point, double *r3, double *N, double *rPrime3, double *reff_value) {
	double gPrime, exp_result, g, exponent;
	int i;

	//R2MAG
	gPrime = 5.0 * (log10(coords * 1000) - 1.0) + absm;

	//REFF
	exp_result = exp(sigmoid_curve_params[1] * (gPrime - sigmoid_curve_params[2]));
	(*reff_value) = sigmoid_curve_params[0] / (exp_result + 1);

	(*rPrime3) = coords * coords * coords;

	for (i = 0; i < ap->convolve; i++) {
		g = gPrime + dx[i];

		//MAG2R
		r_point[i] = pow(10.0, (g - absm)/5.0 + 1.0) / 1000.0;

		r3[i] = r_point[i] * r_point[i] * r_point[i];
		exponent = (g-gPrime) * (g-gPrime) / (2 * stdev * stdev);
		N[i] = coeff * exp(-exponent);
	}
}

void calculate_probabilities(double *r_point, double *r3, double *N, double reff_value, double rPrime3, double *integral_point, ASTRONOMY_PARAMETERS *ap, double *bg_prob, double *st_prob) {
	double sinb, sinl, cosb, cosl, zp;
	double psg, xyzs[3], dotted, xyz_norm;
	double rg, pbx;
	int i, j;

        sinb = sin(integral_point[1] / deg);
        sinl = sin(integral_point[0] / deg);
        cosb = cos(integral_point[1] / deg);
        cosl = cos(integral_point[0] / deg);

	(*bg_prob) = 0;
	/* if q is 0, there is no probability */
	if (q == 0) {
		(*bg_prob) -= 1;
	} else {
		for (i = 0; i < ap->convolve; i++) {
			xyz[i][2] = r_point[i] * sinb;
			zp = r_point[i] * cosb;
			xyz[i][0] = zp * cosl - lbr_r;
			xyz[i][1] = zp * sinl;

			/* background probability */
			rg = sqrt(xyz[i][0]*xyz[i][0] + xyz[i][1]*xyz[i][1] + (xyz[i][2]/q)*(xyz[i][2]/q));
			pbx = 1 / (pow(rg, alpha) * pow(rg + r0, 3 - alpha + delta));

			(*bg_prob) += qgaus_W[i] * pbx * r3[i] * N[i];
		}
	}
	(*bg_prob) *= reff_value * xr / rPrime3;

	for (i = 0; i < ap->number_streams; i++) {
		st_prob[i] = 0;
		if (stream_sigma[i] > -0.0001 && stream_sigma[i] < 0.0001) continue;
		for (j = 0; j < ap->convolve; j++) {
			xyzs[0] = xyz[j][0] - stream_c[i][0];
			xyzs[1] = xyz[j][1] - stream_c[i][1];
			xyzs[2] = xyz[j][2] - stream_c[i][2];

			dotted = stream_a[i][0] * xyzs[0] + stream_a[i][1] * xyzs[1] + stream_a[i][2] * xyzs[2];

			xyzs[0] = xyzs[0] - dotted * stream_a[i][0];
			xyzs[1] = xyzs[1] - dotted * stream_a[i][1];
			xyzs[2] = xyzs[2] - dotted * stream_a[i][2];

			xyz_norm = xyzs[0] * xyzs[0] + xyzs[1] * xyzs[1] + xyzs[2] * xyzs[2];

			psg = exp( -xyz_norm / 2 / (stream_sigma[i] * stream_sigma[i]) );

			st_prob[i] += qgaus_W[j] * r3[j] * N[j] * psg;
		}
		st_prob[i] *= reff_value * xr / rPrime3;
	}
}

double calculate_progress(EVALUATION_STATE *s) {
	long total_calc_probs, current_calc_probs, current_probs, i;
	int mu_step_current, nu_step_current, r_step_current;
	INTEGRAL_AREA *ia;

	total_calc_probs = 0;
	current_calc_probs = 0;

 	for (i = 0; i < s->number_integrals; i++) {
		ia = s->integral[i];

		get_steps(ia, &mu_step_current, &nu_step_current, &r_step_current);
		current_probs = ia->r_steps * ia->mu_steps * ia->nu_steps;
		total_calc_probs += current_probs;
		if (i < s->current_integral) {
			current_calc_probs += current_probs;
		} else if (i == s->current_integral) {
			current_calc_probs += r_step_current + (nu_step_current * ia->r_steps) + (mu_step_current * ia->nu_steps * ia->r_steps);
		}
	}

	total_calc_probs += s->total_stars;
	current_calc_probs += s->current_star_point;

	return (double)current_calc_probs / (double)total_calc_probs;
}

void calculate_integral_unconvolved(ASTRONOMY_PARAMETERS *ap, INTEGRAL_AREA *ia, EVALUATION_STATE *es) {
	int i;
	int mu_step_current, nu_step_current, r_step_current;
	double bg_prob, st_prob, V;
	double *irv, *r_unconvolved;
	double integral_point[3];

	irv = (double*)malloc(sizeof(double) * ia->r_steps);
	r_unconvolved	= (double*)malloc(sizeof(double) * ia->r_steps);
	
	double rPrime, log_r, r, next_r;
	for (i = 0; i < ia->r_steps; i++) {
		log_r	=	ia->r_min + (i * ia->r_step_size);
		r	=	pow(10.0, (log_r-14.2)/5.0);
		next_r	=	pow(10.0, (log_r + ia->r_step_size - 14.2)/5.0);

		irv[i]	= (((next_r * next_r * next_r) - (r * r * r))/3.0) * ia->mu_step_size / deg;
		rPrime	= (next_r+r)/2.0;

		r_unconvolved[i] = rPrime;
	}

	get_steps(ia, &mu_step_current, &nu_step_current, &r_step_current);
	for (; mu_step_current < ia->mu_steps; mu_step_current++) {
		double mu = ia->mu_min + (mu_step_current * ia->mu_step_size);

		for (; nu_step_current < ia->nu_steps; nu_step_current++) {
			double id = 0;
			double nu = ia->nu_min + (nu_step_current * ia->nu_step_size);

			#ifdef GMLE_BOINC
				if (boinc_time_to_checkpoint()) {
					int retval = write_checkpoint(es);
					if (retval) {
						fprintf(stderr,"APP: astronomy checkpoint failed %d\n",retval);
						return;
					}
					boinc_checkpoint_completed();
				}
				boinc_fraction_done(calculate_progress(es));
			#endif

			#if WEDGE_ALLOW_ZERO == 1
				if (ap->wedge > 0) {
					id = cos((90 - nu - ia->nu_step_size)/deg) - cos((90 - nu)/deg);
					if (ap->sgr_coordinates == 0) {
						double ra, dec;
						atGCToEq(mu + 0.5 * ia->mu_step_size, nu + 0.5 * ia->nu_step_size, &ra, &dec, get_node(), wedge_incl(ap->wedge));
						atEqToGal(ra, dec, &integral_point[0], &integral_point[1]);
					} else if (ap->sgr_coordinates == 1) {					
						double lamda, beta;
						gcToSgr(mu + 0.5 * ia->mu_step_size, nu + 0.5 * ia->nu_step_size, ap->wedge, &lamda, &beta);
						sgrToGal(lamda, beta, &integral_point[0], &integral_point[1]);
					} else { 
						printf("Error: ap->sgr_coordinates not valid");
					}
				}
			#else
				id = cos((90 - nu - ia->nu_step_size)/deg) - cos((90 - nu)/deg);
				if (ap->sgr_coordinates == 0) {
					double ra, dec;
					atGCToEq(mu + 0.5 * ia->mu_step_size, nu + 0.5 * ia->nu_step_size, &ra, &dec, get_node(), wedge_incl(ap->wedge));
					atEqToGal(ra, dec, &integral_point[0], &integral_point[1]);
				} else if (ap->sgr_coordinates == 1) {					
					double lamda, beta;
					gcToSgr(mu + 0.5 * ia->mu_step_size, nu + 0.5 * ia->nu_step_size, ap->wedge, &lamda, &beta);
					sgrToGal(lamda, beta, &integral_point[0], &integral_point[1]);
				} else { 
					printf("Error: ap->sgr_coordinates not valid");
				}
			#endif

			for (; r_step_current < ia->r_steps; r_step_current++) {
				#if WEDGE_ALLOW_ZERO == 1
					if (ap->wedge == 0) {
						double xyz[3];
						double log_r = ia->r_min + (r_step_current * ia->r_step_size);
						double r = pow(10.0, (log_r-14.2)/5.0);

						V = ia->mu_step_size * ia->nu_step_size * ia->r_step_size;
						xyz[0] = mu + (0.5 * ia->mu_step_size);
						xyz[1] = nu + (0.5 * ia->nu_step_size);
						xyz[2] = r + (0.5 * ia->r_step_size);
						xyz2lbr(xyz, integral_point);
					} else {
						V = irv[r_step_current] * id;
					}
				#else
					V = irv[r_step_current] * id;
				#endif

				integral_point[2] = r_unconvolved[r_step_current];
				bg_prob = stPbx(integral_point, ap->background_parameters);
				ia->background_integral += bg_prob * V;
				for (i = 0; i < ap->number_streams; i++) {
					st_prob = stPsg(integral_point, ap->stream_parameters[i], ap->wedge, ap->sgr_coordinates);
					ia->stream_integrals[i] += st_prob * V;
				}
			}
			r_step_current = 0;
		}
		nu_step_current = 0;
	}
	mu_step_current = 0;

	free(irv);
	free(r_unconvolved);
}

void calculate_integral_convolved(ASTRONOMY_PARAMETERS *ap, INTEGRAL_AREA *ia, EVALUATION_STATE *es) {
	int i;
	int mu_step_current, nu_step_current, r_step_current;
	double bg_prob, *st_probs, V;
	double *irv;
	double *rPrime3, *reff_value, **N, **r3, **r_point;
	double integral_point[3];

	irv		= (double*)malloc(sizeof(double) * ia->r_steps);
	st_probs	= (double*)malloc(sizeof(double) * ap->number_streams); 
	rPrime3		= (double*)malloc(sizeof(double) * ia->r_steps);
	reff_value	= (double*)malloc(sizeof(double) * ia->r_steps);
	N		= (double**)malloc(sizeof(double*) * ia->r_steps);
	r3		= (double**)malloc(sizeof(double*) * ia->r_steps);
	r_point		= (double**)malloc(sizeof(double*) * ia->r_steps);
	
	double rPrime, log_r, r, next_r;

	for (i = 0; i < ia->r_steps; i++) {
		log_r	=	ia->r_min + (i * ia->r_step_size);
		r	=	pow(10.0, (log_r-14.2)/5.0);
		next_r	=	pow(10.0, (log_r + ia->r_step_size - 14.2)/5.0);

		irv[i]	= (((next_r * next_r * next_r) - (r * r * r))/3.0) * ia->mu_step_size / deg;
		rPrime	= (next_r+r)/2.0;

		r_point[i] = (double*)malloc(sizeof(double) * ap->convolve);
		r3[i] = (double*)malloc(sizeof(double) * ap->convolve);
		N[i] = (double*)malloc(sizeof(double) * ap->convolve);

		set_probability_constants(ap, rPrime, r_point[i], r3[i], N[i], &(rPrime3[i]), &(reff_value[i]));
	}

	get_steps(ia, &mu_step_current, &nu_step_current, &r_step_current);
	for (; mu_step_current < ia->mu_steps; mu_step_current++) {
		double mu = ia->mu_min + (mu_step_current * ia->mu_step_size);

		for (; nu_step_current < ia->nu_steps; nu_step_current++) {
			double id = 0;
			double nu = ia->nu_min + (nu_step_current * ia->nu_step_size);

			#ifdef GMLE_BOINC
				if (boinc_time_to_checkpoint()) {
					int retval = write_checkpoint(es);
					if (retval) {
						fprintf(stderr,"APP: astronomy checkpoint failed %d\n",retval);
						return;
					}
					boinc_checkpoint_completed();
				}
				boinc_fraction_done(calculate_progress(es));
			#endif

			#if WEDGE_ALLOW_ZERO == 1
				if (ap->wedge > 0) {
					id = cos((90 - nu - ia->nu_step_size)/deg) - cos((90 - nu)/deg);
					if (ap->sgr_coordinates == 0) {
						double ra, dec;
						atGCToEq(mu + 0.5 * ia->mu_step_size, nu + 0.5 * ia->nu_step_size, &ra, &dec, get_node(), wedge_incl(ap->wedge));
						atEqToGal(ra, dec, &integral_point[0], &integral_point[1]);
					} else if (ap->sgr_coordinates == 1) {					
						double lamda, beta;
						gcToSgr(mu + 0.5 * ia->mu_step_size, nu + 0.5 * ia->nu_step_size, ap->wedge, &lamda, &beta);
						sgrToGal(lamda, beta, &integral_point[0], &integral_point[1]);
					} else { 
						printf("Error: ap->sgr_coordinates not valid");
					}
				}
			#else
				id = cos((90 - nu - ia->nu_step_size)/deg) - cos((90 - nu)/deg);
				if (ap->sgr_coordinates == 0) {
					double ra, dec;
					atGCToEq(mu + 0.5 * ia->mu_step_size, nu + 0.5 * ia->nu_step_size, &ra, &dec, get_node(), wedge_incl(ap->wedge));
					atEqToGal(ra, dec, &integral_point[0], &integral_point[1]);
				} else if (ap->sgr_coordinates == 1) {					
					double lamda, beta;
					gcToSgr(mu + 0.5 * ia->mu_step_size, nu + 0.5 * ia->nu_step_size, ap->wedge, &lamda, &beta);
					sgrToGal(lamda, beta, &integral_point[0], &integral_point[1]);
				} else { 
					printf("Error: ap->sgr_coordinates not valid");
				}
			#endif

			for (; r_step_current < ia->r_steps; r_step_current++) {
				#if WEDGE_ALLOW_ZERO == 1
					if (ap->wedge == 0) {
						double xyz[3];
						double log_r = ia->r_min + (r_step_current * ia->r_step_size);
						double r = pow(10.0, (log_r-14.2)/5.0);

						V = ia->mu_step_size * ia->nu_step_size * ia->r_step_size;
						xyz[0] = mu + (0.5 * ia->mu_step_size);
						xyz[1] = nu + (0.5 * ia->nu_step_size);
						xyz[2] = r + (0.5 * ia->r_step_size);
						xyz2lbr(xyz, integral_point);
					} else {
						V = irv[r_step_current] * id;
					}
				#else
					V = irv[r_step_current] * id;
				#endif

				calculate_probabilities(r_point[r_step_current], r3[r_step_current], N[r_step_current], reff_value[r_step_current], rPrime3[r_step_current], integral_point, ap, &bg_prob, st_probs);
				ia->background_integral += bg_prob * V;
				for (i = 0; i < ap->number_streams; i++) {
					ia->stream_integrals[i] += st_probs[i] * V;
				}
				ia->current_calculation++;
				if (ia->current_calculation >= ia->max_calculation) break;
			}
			if (ia->current_calculation >= ia->max_calculation) break;
			r_step_current = 0;
		}
		if (ia->current_calculation >= ia->max_calculation) break;
		nu_step_current = 0;
	}
	mu_step_current = 0;

	free(irv);
	free(st_probs);
	free(rPrime3);
	free(reff_value);
	for (i = 0; i < ia->r_steps; i++) {
		free(N[i]);
		free(r_point[i]);
		free(r3[i]);
	}
	free(N);
	free(r_point);
	free(r3);
}

int calculate_integrals(ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, STAR_POINTS* sp) {
	int i, j;
//	time_t start_time, finish_time;
//	time(&start_time);

	#ifdef GMLE_BOINC
		read_checkpoint(es);
	#endif

	init_constants(ap);
	for (; es->current_integral < ap->number_integrals; es->current_integral++) {
		if (ap->convolve > 0) {
			calculate_integral_convolved(ap, es->integral[es->current_integral], es);
		} else {
			calculate_integral_unconvolved(ap, es->integral[es->current_integral], es);
		}
//		printf("bg_int: %lf\n", es->integral[es->current_integral]->background_integral);
//		for (i = 0; i < ap->number_streams; i++) printf("st_int[%d]: %lf\n", i, es->integral[es->current_integral]->stream_integrals[i]);
	}

	es->background_integral = es->integral[0]->background_integral;
	for (i = 0; i < ap->number_streams; i++) {
		es->stream_integrals[i] = es->integral[0]->stream_integrals[i];
	}

	for (i = 1; i < ap->number_integrals; i++) {
		es->background_integral -= es->integral[i]->background_integral;
		for (j = 0; j < ap->number_streams; j++) es->stream_integrals[j] -= es->integral[i]->stream_integrals[j];
	}

//	time(&finish_time);
//	printf("background integral: %.20lf\n", es->background_integral);
//	for (i = 0; i < ap->number_streams; i++) printf("stream integral[%d]: %.20lf\n", i, es->stream_integrals[i]);
//	printf("integrals calculated in: %lf\n", (double)finish_time - (double)start_time);
	return 0;
}

int calculate_likelihood(ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, STAR_POINTS* sp) {
	int i, current_stream;
	double bg_prob, *st_prob;
	double background_weight, sum_exp_weights, *exp_stream_weights;
	double *r_point, *r3, *N, rPrime3, reff_value;
//	time_t start_time, finish_time;
//	time (&start_time);

	background_weight = exp(ap->background_weight);
	sum_exp_weights = 0.0;
	exp_stream_weights = (double*)malloc(sizeof(double) * ap->number_streams);
	for (i = 0; i < ap->number_streams; i++) {
		exp_stream_weights[i] = exp(ap->stream_weights[i]);
		sum_exp_weights += exp(ap->stream_weights[i]);
	}
	sum_exp_weights += background_weight;

	if (ap->convolve > 0) {
		st_prob = (double*)malloc(sizeof(double) * ap->number_streams);
		r_point = (double*)malloc(sizeof(double) * ap->convolve);
		r3 = (double*)malloc(sizeof(double) * ap->convolve);
		N = (double*)malloc(sizeof(double) * ap->convolve);

		for (; es->current_star_point < sp->number_stars; es->current_star_point++) {
			double star_prob = 0.0;
			double star_coords[3];
			#if NEW_FORUMLA == 1
				double sum_integrals = 0.0;
			#endif

			#ifdef GMLE_BOINC
				if (boinc_time_to_checkpoint()) {
					int retval = write_checkpoint(es);
					if (retval) {
						fprintf(stderr,"APP: astronomy checkpoint failed %d\n",retval);
						return retval;
					}
					boinc_checkpoint_completed();
				}
				boinc_fraction_done(calculate_progress(es));
			#endif

			star_coords[0] = sp->stars[es->current_star_point][0];
			star_coords[1] = sp->stars[es->current_star_point][1];
			star_coords[2] = sp->stars[es->current_star_point][2];

			set_probability_constants(ap, star_coords[2], r_point, r3, N, &rPrime3, &reff_value);
			calculate_probabilities(r_point, r3, N, reff_value, rPrime3, star_coords, ap, &bg_prob, st_prob);
			#if NEW_FORMULA == 1
				for (current_stream = 0; current_stream < ap->number_streams; current_stream++) {
					star_prob += st_prob[current_stream] * exp_stream_weights[current_stream];
					sum_integrals += es->stream_integrals[current_stream] * exp_stream_weights[current_stream];
				}
				star_prob += bg_prob * background_weight;
			#else
				for (current_stream = 0; current_stream < ap->number_streams; current_stream++) {
					star_prob += (st_prob[current_stream]/es->stream_integrals[current_stream]) * exp_stream_weights[current_stream];
				}
				star_prob += (bg_prob/es->background_integral) * background_weight;
			#endif

			//Calculate the probability for this star.
			star_prob /= sum_exp_weights;
			#if NEW_FORMULA == 1
				star_prob /= sum_integrals;		
			#endif		

			//update: check star_prob==0, not prob_sum
			if (star_prob != 0.0) {
				es->prob_sum += log(star_prob)/log(10.0);
			} else {
				es->num_zero++;
				es->prob_sum -= 238.0;
			}
		}

		free(st_prob);
		free(r_point);
		free(r3);
		free(N);
	} else {
		for (; es->current_star_point < sp->number_stars; es->current_star_point++) {
			double star_prob = 0.0;
			double star_coords[3];
			#if NEW_FORMULA == 1
				double sum_integrals = 0.0;
			#endif

			#ifdef GMLE_BOINC
				if (boinc_time_to_checkpoint()) {
					int retval = write_checkpoint(es);
					if (retval) {
						fprintf(stderr,"APP: astronomy checkpoint failed %d\n",retval);
						return retval;
					}
					boinc_checkpoint_completed();
				}
				boinc_fraction_done(calculate_progress(es));
			#endif

			star_coords[0] = sp->stars[es->current_star_point][0];
			star_coords[1] = sp->stars[es->current_star_point][1];
			star_coords[2] = sp->stars[es->current_star_point][2];

			#if NEW_FORMULA == 1
				for (current_stream = 0; current_stream < ap->number_streams; current_stream++) {
					star_prob += stPsg(star_coords, ap->stream_parameters[current_stream], ap->wedge, ap->sgr_coordinates) * exp_stream_weights[current_stream];
					sum_integrals += es->stream_integrals[current_stream] * exp_stream_weights[current_stream];
				}
				star_prob += stPbx(star_coords, ap->background_parameters) * background_weight;
			#else
				for (current_stream = 0; current_stream < ap->number_streams; current_stream++) {
					star_prob += (stPsg(star_coords, ap->stream_parameters[current_stream], ap->wedge, ap->sgr_coordinates)/es->stream_integrals[current_stream]) * exp_stream_weights[current_stream];
				}
				star_prob += (stPbx(star_coords, ap->background_parameters)/es->background_integral) * background_weight;
			#endif

			//Calculate the probability for this star.
			star_prob /= sum_exp_weights;
			#if NEW_FORMULA == 1
				star_prob /= sum_integrals;		
			#endif

			//update: check star_prob==0, not prob_sum
			if (star_prob != 0.0) {
				es->prob_sum += log(star_prob)/log(10.0);
			} else {
				es->num_zero++;
				es->prob_sum -= 238.0;
			}
		}
	}

	free_constants(ap);
	free(exp_stream_weights);

//	time(&finish_time);
//	printf("likelihood calculated in: %lf\n", (double)finish_time - (double)start_time);

	return 0;
}
