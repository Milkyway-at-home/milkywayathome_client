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

#ifndef ASTRONOMY_EVALUATION_H
#define ASTRONOMY_EVALUATION_H

#include "evaluation_state.h"
#include "parameters.h"
#include "star_points.h"

typedef struct
{
    vector stream_a;
    vector stream_c;
    double stream_sigma;
    double stream_sigma_sq2;
} STREAM_CONSTANTS;

/* FIXME: Better name */
typedef struct
{
    double alpha, q, r0, delta, coeff, alpha_delta3;
    double bg_a, bg_b, bg_c;
} STREAM_NUMS;

typedef struct
{
    double st_prob;
    double st_prob_int;    /* for kahan summation */
    double st_prob_int_c;
} ST_PROBS;

/* Replacement for mess of double** used before */
typedef struct
{
    double r_point;
    double r_in_mag;
    double r_in_mag2;
    double qw_r3_N;
} R_STEP_STATE;

/* Scratch space used by each integral */
typedef struct
{
    double* irv;
    double* reff_xr_rp3;
    double* ids;
    double* nus;
    ST_PROBS* probs;
    R_STEP_STATE* rss;
} INTEGRAL_STATE;


STREAM_CONSTANTS* init_constants(ASTRONOMY_PARAMETERS* ap, STREAM_NUMS* sn);
void free_constants(ASTRONOMY_PARAMETERS* ap, STREAM_CONSTANTS* sc);

void set_probability_constants(const STREAM_NUMS* sn,
                               unsigned int n_convolve,
                               double coords,
                               R_STEP_STATE* rss,
                               unsigned int r_step,
                               double* reff_xr_rp3);

void calculate_probabilities(const ASTRONOMY_PARAMETERS* ap,
                             const STREAM_CONSTANTS* sc,
                             const STREAM_NUMS* sn,
                             R_STEP_STATE* rss,
                             unsigned int r_step_current,
                             unsigned int r_steps,
                             double reff_xr_rp3,
                             double* integral_point,
                             double* bg_prob,
                             ST_PROBS* probs);

int calculate_integrals(const ASTRONOMY_PARAMETERS* ap,
                        const STREAM_CONSTANTS* sc,
                        const STREAM_NUMS* sn,
                        EVALUATION_STATE* es);

int calculate_likelihood(const ASTRONOMY_PARAMETERS* ap,
                         const STREAM_CONSTANTS* sc,
                         const STREAM_NUMS* sn,
                         EVALUATION_STATE* es,
                         const STAR_POINTS* sp);

void cpu__r_constants(const STREAM_NUMS* sn,
                      unsigned int n_convolve,
                      unsigned int r_steps, double r_min, double r_step_size,
                      double mu_step_size,
                      unsigned int nu_steps, double nu_min, double nu_step_size,
                      double* irv,
                      R_STEP_STATE* rss,
                      double* reff_xr_rp3, double* nus, double* ids);

double cpu_evaluate(const ASTRONOMY_PARAMETERS* ap,
                    const STAR_POINTS* sp,
                    const STREAM_CONSTANTS* sc,
                    const STREAM_NUMS* sn);

#endif

