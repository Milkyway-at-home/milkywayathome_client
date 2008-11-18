/*
 *  integral_function.c
 *  Astronomy
 *
 *  Created by Travis Desell on 2/21/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
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

#include "evaluation.h"
#include "parameters.h"
#include "probability.h"
#include "stCoords.h"
#include "atSurveyGeometry.h"
#include "star_points.h"
#include "numericalIntegration.h"
#include "../evaluation/mpi_evaluator.h"

#ifndef _WIN32
	#define pi M_PI
#else
	#define pi 3.14159265358979323846
#endif

#define deg (180.0/pi)

void initialize_state(EVALUATION_STATE* es, int number_streams) {
	int i;
	es->r_step_current = 0;
	es->mu_step_current = 0;
	es->nu_step_current = 0;
	
	es->r_cut_step_current = 0;
	es->mu_cut_step_current = 0;
	es->nu_cut_step_current = 0;

	es->number_streams = number_streams;

	es->main_integral_calculated = 0;
	es->current_cut = 0;

	es->background_integral = 0.0;
	es->stream_integrals = (double*)malloc(sizeof(double) * number_streams);
	for (i = 0; i < number_streams; i++) {
		es->stream_integrals[i] = 0.0;
	}

	es->current_star_point = 0;
	es->num_zero = 0;
	es->bad_jacobians = 0;
	es->prob_sum = 0.0;
}

void reset_evaluation_state(EVALUATION_STATE *es) {
	int i;
	es->r_step_current = 0;
	es->mu_step_current = 0;
	es->nu_step_current = 0;
	
	es->r_cut_step_current = 0;
	es->mu_cut_step_current = 0;
	es->nu_cut_step_current = 0;

	es->main_integral_calculated = 0;
	es->current_cut = 0;

	es->background_integral = 0.0;
	for (i = 0; i < es->number_streams; i++) {
		es->stream_integrals[i] = 0.0;
	}

	es->current_star_point = 0;
	es->num_zero = 0;
	es->bad_jacobians = 0;
	es->prob_sum = 0.0;
}

void free_state(EVALUATION_STATE* es) {
	free( es->stream_integrals );
}

#ifdef GMLE_BOINC
	int write_checkpoint(EVALUATION_STATE* es) {
		int i;
		char output_path[512];
		boinc_resolve_filename(CHECKPOINT_FILE, output_path, sizeof(output_path));

		MFILE cp;
		int retval = cp.open(output_path, "w");
		if (retval) {
	                fprintf(stderr, "APP: error writing checkpoint (opening checkpoint file) %d\n", retval);
	                return retval;
		}

		cp.printf("r_step_current: %d\n", es->r_step_current);
		cp.printf("mu_step_current: %d\n", es->mu_step_current);
		cp.printf("nu_step_current: %d\n", es->nu_step_current);

		cp.printf("main_integral_calculated: %d\n", es->main_integral_calculated);
		cp.printf("current_cut: %d\n", es->current_cut);
		cp.printf("background_integral: %lf\n", es->background_integral);
		cp.printf("stream_integrals [%d]:", es->number_streams);
		for (i = 0; i < es->number_streams; i++) cp.printf(" %lf", es->stream_integrals[i]);
		cp.printf("\ncurrent_star_point: %d\n", es->current_star_point);
		cp.printf("num_zero: %d\n", es->num_zero);
		cp.printf("bad_jacobians: %d\n", es->bad_jacobians);
		cp.printf("prob_sum: %lf\n", es->prob_sum);

		retval = cp.flush();
		if (retval) {
	                fprintf(stderr, "APP: error writing checkpoint (flushing checkpoint file) %d\n", retval);
	                return retval;
		}
		retval = cp.close();
		if (retval) {
	                fprintf(stderr, "APP: error writing checkpoint (closing checkpoint file) %d\n", retval);
	                return retval;
		}
		return 0;
	}

	int read_checkpoint(EVALUATION_STATE* es) {
		int i;
		char input_path[512];
		int retval = boinc_resolve_filename(CHECKPOINT_FILE, input_path, sizeof(input_path));
		if (retval) {
			return 0;
		}

	        FILE* data_file = boinc_fopen(input_path, "r");
	        if (!data_file) {
	                fprintf(stderr, "APP: error reading checkpoint (opening file)\n");
	                return 1;
	        }
		fscanf(data_file, "r_step_current: %d\n", &es->r_step_current);
		fscanf(data_file, "mu_step_current: %d\n", &es->mu_step_current);
		fscanf(data_file, "nu_step_current: %d\n", &es->nu_step_current);
	
		fscanf(data_file, "main_integral_calculated: %d\n", &es->main_integral_calculated);
		fscanf(data_file, "current_cut: %d\n", &es->current_cut);
		fscanf(data_file, "background_integral: %lf\n", &es->background_integral);
		fscanf(data_file, "stream_integrals [%d]:", &es->number_streams);
		for (i = 0; i < es->number_streams; i++) fscanf(data_file, " %lf", &es->stream_integrals[i]);

		fscanf(data_file, "\ncurrent_star_point: %d\n", &es->current_star_point);
		fscanf(data_file, "num_zero: %d\n", &es->num_zero);
		fscanf(data_file, "bad_jacobians: %d\n", &es->bad_jacobians);
		fscanf(data_file, "prob_sum: %lf\n", &es->prob_sum);

	        fclose(data_file);
	        return 0;
	}
#endif

int calculate_integrals(ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, STAR_POINTS* sp) {
	int s;
	int first_run = 1;
	double V = 0;
	time_t start_time, finish_time;

	time(&start_time);
	
	setWeights(ap->convolve);

	es->background_integral = 0;
	for (s = 0; s < ap->number_streams; s++) {
		es->stream_integrals[s] = 0;
	}

	#ifdef GMLE_BOINC
		int retval = read_checkpoint(es);
//		if (retval) {
//			fprintf(stderr,"APP: failed reading checkpoint %d\n", retval);
//			return retval;
//		}
	#endif

	if (!es->main_integral_calculated) {
		for (; es->mu_step_current < ap->mu_steps; es->mu_step_current++) {
			double mu = ap->mu_min + (es->mu_step_current * ap->mu_step_size);

			for (; es->nu_step_current < ap->nu_steps; es->nu_step_current++) {
				double nu = ap->nu_min + (es->nu_step_current * ap->nu_step_size);

				for (; es->r_step_current < ap->r_steps; es->r_step_current++) {
					double integral_point[3], xyz[3];

					double log_r = ap->r_min + (es->r_step_current * ap->r_step_size);
					double r = pow(10.0, (log_r-14.2)/5.0);
					double next_r = pow(10.0, (log_r+ap->r_step_size-14.2)/5.0);

					if (ap->wedge > 0) {
						double ir = (pow(next_r,3.0) - pow(r, 3.0))/3.0;
						double id = cos((90 - nu - ap->nu_step_size)/deg) - cos((90 - nu)/deg);

						V = ir * id * ap->mu_step_size / deg;

						double ra = 0.0;
						double dec = 0.0;
						double point0 = 0.0;
						double point1 = 0.0;
						double lamda = 0.0;
						double beta = 0.0;
				
						if (ap->sgr_coordinates == 0) {
							atGCToEq(mu + 0.5 * ap->mu_step_size, nu + 0.5 * ap->nu_step_size, &ra, &dec, get_node(), wedge_incl(ap->wedge));
							atEqToGal(ra, dec, &point0, &point1);
						} else if (ap->sgr_coordinates == 1) {					
							gcToSgr(mu + 0.5 * ap->mu_step_size, nu + 0.5 * ap->nu_step_size, ap->wedge, &lamda, &beta); //vickej2
							sgrToGal(lamda, beta, &point0, &point1); //vickej2
						} else { 
							printf("Error: ap->sgr_coordinates not valid");
						}
	
						integral_point[0] = point0;
						integral_point[1] = point1;
						integral_point[2] = (next_r+r)/2.0;
					
						//vickej2 <<<testing if r stays within its bounds>>>
                		   		//printf("<<<%f>>>\n", r); //vickej2
						//vickej2 <<<end>>>
					} else {
						V = ap->mu_step_size * ap->nu_step_size * ap->r_step_size;
						xyz[0] = mu + (0.5 * ap->mu_step_size);
						xyz[1] = nu + (0.5 * ap->nu_step_size);
						xyz[2] = r + (0.5 * ap->r_step_size);
						xyz2lbr(xyz, integral_point);
					}

					double bg_prob = 0.0;
					double st_prob = 0.0;
							   
					if (ap->convolve > 0) {
						bg_prob = stPbxConvolved(integral_point, ap->background_parameters, ap->wedge, ap->convolve);
					} else {
						bg_prob = stPbx(integral_point, ap->background_parameters);
					}
					es->background_integral += bg_prob * V;
					for (s = 0; s < ap->number_streams; s++) {
						if (ap->convolve > 0) {
							st_prob = stPsgConvolved(integral_point, ap->stream_parameters[s], ap->wedge, ap->convolve, ap->sgr_coordinates);
						} else {
							st_prob = stPsg(integral_point, ap->stream_parameters[s], ap->wedge, ap->sgr_coordinates);
						}
						es->stream_integrals[s] += st_prob * V;
	 				}
					first_run = 0;

					#ifdef GMLE_BOINC
						if (boinc_time_to_checkpoint()) {
							retval = write_checkpoint(es);
							if (retval) {
								fprintf(stderr,"APP: astronomy checkpoint failed %d\n",retval);
								return retval;
							}
							boinc_checkpoint_completed();
						}
						double f = es->r_step_current + (es->nu_step_current * ap->r_steps) + (es->mu_step_current * ap->nu_steps * ap->r_steps);
						f /= (ap->mu_steps * ap->nu_steps * ap->r_steps);
						f *= 0.5;
						boinc_fraction_done(f);
					#endif
				}
				es->r_step_current = 0;
			}
			es->nu_step_current = 0;
		}
		es->mu_step_current = 0;
	}
	es->main_integral_calculated = 1;
	printf("background_integral: %.10lf\n", es->background_integral);
	for (s = 0; s < ap->number_streams; s++) {
		printf("stream_integral[s]: %.10lf\n", es->stream_integrals[s]);
	}



	/*** Begin volume removal ***/
	for (; es->current_cut < ap->number_cuts; es->current_cut++) {
        	for (; es->mu_cut_step_current < ap->mu_cut[es->current_cut][2]; es->mu_cut_step_current++) {
        	        double mu = ap->mu_cut[es->current_cut][0] + (es->mu_cut_step_current * ap->mu_cut_step_size[es->current_cut]);
	
	                for (; es->nu_cut_step_current < ap->nu_cut[es->current_cut][2]; es->nu_cut_step_current++) {
	                        double nu = ap->nu_cut[es->current_cut][0] + (es->nu_cut_step_current * ap->nu_cut_step_size[es->current_cut]);
	
	                        for (; es->r_cut_step_current < ap->r_cut[es->current_cut][2]; es->r_cut_step_current++) {
	                                double integral_point[3], xyz[3];

        	                        double log_r = ap->r_cut[es->current_cut][0] + (es->r_cut_step_current * ap->r_cut_step_size[es->current_cut]);
        	                        double r = pow(10.0, (log_r-14.2)/5.0);
        	                        double next_r = pow(10.0, (log_r+ap->r_cut_step_size[es->current_cut]-14.2)/5.0);
	
	                                if (ap->wedge > 0) {
	                                        double ir = (pow(next_r,3.0) - pow(r, 3.0))/3.0;
	                                        double id = cos((90 - nu - ap->nu_cut_step_size[es->current_cut])/deg) - cos((90 - nu)/deg);
	
	                                        V = ir * id * ap->mu_cut_step_size[es->current_cut] / deg;
	
	                                        double ra = 0.0;
	                                        double dec = 0.0;
	                                        double point0 = 0.0;
	                                        double point1 = 0.0;
	                                        double lamda = 0.0;
	                                        double beta = 0.0;
	                                        //vickej2 <<<make sure all my parameters are being taken correctly>>>
	                                        //printf("rmax=%f",ap->r_max); //vickej2
	                                        //printf("wedge=%i, r_steps=%i, mu_steps=%i, nu_steps=%i, nu_min=%f, nu_max=%f, r_min=%f, r_max=%f, mu_min=%f, mu_max=%f, nu_step_size=%f, r_step_size=%f, mu_step_size=%f", ap->wedge, ap->r_steps, ap->mu_steps, ap->nu_steps, ap->nu_min, ap->nu_max, ap->r_min, ap->r_max, ap->mu_min, ap->mu_max, ap->nu_step_size, ap->r_step_size, ap->mu_step_size);  //vickej2
	                                        //vickej2 <<<end>>>
	
	                                        if (ap->sgr_coordinates == 0) {
	                                                atGCToEq(mu + 0.5 * ap->mu_cut_step_size[es->current_cut], nu + 0.5 * ap->nu_cut_step_size[es->current_cut], &ra, &dec, get_node(), wedge_incl(ap->wedge));
	                                                atEqToGal(ra, dec, &point0, &point1);
       		                               	} else if (ap->sgr_coordinates == 1) {                                 
                                                	gcToSgr(mu + 0.5 * ap->mu_cut_step_size[es->current_cut], nu + 0.5 * ap->nu_cut_step_size[es->current_cut], ap->wedge, &lamda, &beta); //vickej2
                                                	sgrToGal(lamda, beta, &point0, &point1); //vickej2
	
	                                                //vickej2 <<<make sure the conversion is correct (check with conversiontester.vb)>>>
	                                                //printf(" mui=%f, nui=%f, lamda=%f, beta=%f, l=%f, b=%f", mu + 0.5 * ap->mu_step_size, nu + 0.5 * ap->nu_step_size, lamda, beta, point0, point1);  //vickej2
	                                                //vickej2 <<<end>>>
	                                        } else {
	                                                printf("Error: ap->sgr_coordinates not valid");
	                                        }
	
						integral_point[0] = point0;
						integral_point[1] = point1;
						integral_point[2] = (next_r+r)/2.0;
	
	                                        //vickej2 <<<testing if r stays within its bounds>>>
	                                        //printf("<<<%f>>>\n", r); //vickej2
	                                        //vickej2 <<<end>>>
	                                } else {
	                                        V = ap->mu_cut_step_size[es->current_cut] * ap->nu_cut_step_size[es->current_cut] * ap->r_cut_step_size[es->current_cut];
	                                        xyz[0] = mu + (0.5 * ap->mu_cut_step_size[es->current_cut]);
	                                        xyz[1] = nu + (0.5 * ap->nu_cut_step_size[es->current_cut]);
	                                        xyz[2] = r + (0.5 * ap->r_cut_step_size[es->current_cut]);
	                                        xyz2lbr(xyz, integral_point);
	                                }
	                                double bg_prob = 0.0;
	                                double st_prob = 0.0;
	
					if (ap->convolve > 0) {
						bg_prob = stPbxConvolved(integral_point, ap->background_parameters, ap->wedge, ap->convolve);
					} else {
						bg_prob = stPbx(integral_point, ap->background_parameters);
					}
					es->background_integral -= bg_prob * V;

					for (s = 0; s < ap->number_streams; s++) {
						if (ap->convolve > 0) {
							st_prob = stPsgConvolved(integral_point, ap->stream_parameters[s], ap->wedge, ap->convolve, ap->sgr_coordinates);
						} else {
							st_prob = stPsg(integral_point, ap->stream_parameters[s], ap->wedge, ap->sgr_coordinates);
						}
						es->stream_integrals[s] -= st_prob * V;
					}
	                                first_run = 0;
	
	                                #ifdef GMLE_BOINC
	                                        if (boinc_time_to_checkpoint()) {
	                                                retval = write_checkpoint(es);
	                                                if (retval) {
	                                                        fprintf(stderr,"APP: astronomy checkpoint failed %d\n",retval);
	                                                        return retval;
	                                                }
	                                                boinc_checkpoint_completed();
	                                        }
	                                        double f = es->r_cut_step_current + (es->nu_cut_step_current * ap->r_cut[es->current_cut][2]) + (es->mu_cut_step_current * ap->nu_cut[es->current_cut][2] * ap->r_cut[es->current_cut][2]);
	                                        f /= (ap->mu_cut[es->current_cut][2] * ap->nu_cut[es->current_cut][2] * ap->r_cut[es->current_cut][2]);
	                                        f *= 0.5;
	                                        boinc_fraction_done(f);
	                                #endif
	                        }
				es->r_cut_step_current = 0;
	                }
			es->nu_cut_step_current = 0;
	        }
		es->mu_cut_step_current = 0;
	}
	/*** End volume removal ***/

	#ifdef GMLE_BOINC
		retval = write_checkpoint(es);
		if (retval) {
			fprintf(stderr,"APP: astronomy checkpoint failed %d\n",retval);
			return retval;
		}
	#endif

	time(&finish_time);
	printf("integrals calculated in: %lf\n", (double)finish_time - (double)start_time);

	return 0;
}

int calculate_likelihood(ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, STAR_POINTS* sp) {
	int i, current_stream;
	int new_formula = 0;
	time_t start_time, finish_time;

	time(&start_time);

	double background_weight = exp(ap->background_weight);
	double sum_exp_weights = 0.0;
	double* exp_stream_weights = (double*)malloc(sizeof(double) * ap->number_streams);
	for (i = 0; i < ap->number_streams; i++) {
		exp_stream_weights[i] = exp(ap->stream_weights[i]);
		sum_exp_weights += exp(ap->stream_weights[i]);
	}
	sum_exp_weights += background_weight;

	int first_run = 1;
	for (; es->current_star_point < sp->number_stars; es->current_star_point++) {
		double star_prob = 0.0;
		double sum_integrals = 0.0;
		double* star_coords = sp->stars[es->current_star_point];

		for (current_stream = 0; current_stream < ap->number_streams; current_stream++) {
			if (new_formula) {
				if (ap->convolve > 0) {
					star_prob += stPsgConvolved(star_coords, ap->stream_parameters[current_stream], ap->wedge, ap->convolve, ap->sgr_coordinates) * exp_stream_weights[current_stream];
				} else {
					star_prob += stPsg(star_coords, ap->stream_parameters[current_stream], ap->wedge, ap->sgr_coordinates) * exp_stream_weights[current_stream];
				}
				sum_integrals += exp_stream_weights[current_stream] * es->stream_integrals[current_stream];
			} else {
				if (ap->convolve > 0) {
					star_prob += (stPsgConvolved(star_coords, ap->stream_parameters[current_stream], ap->wedge, ap->convolve, ap->sgr_coordinates)/es->stream_integrals[current_stream]) * exp_stream_weights[current_stream];
				} else {
					star_prob += (stPsg(star_coords, ap->stream_parameters[current_stream], ap->wedge, ap->sgr_coordinates)/es->stream_integrals[current_stream]) * exp_stream_weights[current_stream];
				}
			}
			first_run = 0;
		}

		//Find the background probability
		if (new_formula) {
			if (ap->convolve > 0) {
				star_prob += stPbxConvolved(star_coords, ap->background_parameters, ap->wedge, ap->convolve) * background_weight;
			} else {
				star_prob += stPbx(star_coords, ap->background_parameters) * background_weight;
			}
			sum_integrals += background_weight * es->background_integral;
		} else {
			if (ap->convolve > 0) {
				star_prob += (stPbxConvolved(star_coords, ap->background_parameters, ap->wedge, ap->convolve)/es->background_integral) * background_weight;
			} else {
				star_prob += (stPbx(star_coords, ap->background_parameters)/es->background_integral) * background_weight;
			}
		}

		//Calculate the probability for this star.
		star_prob /= sum_exp_weights;
		if (new_formula) star_prob /= sum_integrals;		
		
		//update: check star_prob==0, not prob_sum
		if (star_prob != 0.0) {
			es->prob_sum += log(star_prob)/log(10.0);
		} else {
			es->num_zero++;
			es->prob_sum -= 238.0;
		}

		#ifdef GMLE_BOINC
			if (boinc_time_to_checkpoint()) {
				int retval = write_checkpoint(es);
				if (retval) {
					fprintf(stderr,"APP: astronomy checkpoint failed %d\n",retval);
					return retval;
				}
				boinc_checkpoint_completed();
			}

			double f = (double)es->current_star_point/(double)sp->number_stars;
			f *= 0.5;
			f += 0.5;
			boinc_fraction_done(f);
		#endif
	}


	#ifdef GMLE_BOINC
		int retval = write_checkpoint(es);
		if (retval) {
			fprintf(stderr,"APP: astronomy checkpoint failed %d\n",retval);
			return retval;
		}
	#endif

	free(exp_stream_weights);
	freeWeights();

	time(&finish_time);
	printf("likelihood calculated in: %lf\n", (double)finish_time - (double)start_time);

	return 0;
}
