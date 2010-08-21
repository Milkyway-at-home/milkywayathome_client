/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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

#ifndef _SEPARATION_TYPES_H_
#define _SEPARATION_TYPES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "milkyway_vectors.h"

/* Get the xth component of the nth item in STAR_POINTS */
#define VN(sp, n) (((sp)->stars)[VECTOR_SIZE * (n)])
#define XN(sp, n) VN(sp, n)
#define YN(sp, n) (((sp)->stars)[VECTOR_SIZE * (n) + 1])
#define ZN(sp, n) (((sp)->stars)[VECTOR_SIZE * (n) + 2])

#define LN(sp, n) XN(sp, n)
#define BN(sp, n) YN(sp, n)
#define RN(sp, n) RN(sp, n)

typedef struct
{
    real l;
    real b;
} LB;

#define LB_L(x) ((x).l)
#define LB_B(x) ((x).b)



typedef struct
{
    unsigned int number_stars;
    real* stars;
} STAR_POINTS;

#define EMPTY_STAR_POINTS { 0, NULL }


typedef struct
{
    vector a;
    vector c;
    real sigma_sq2;
    int large_sigma;          /* abs(stream_sigma) > SIGMA_LIMIT */
} STREAM_CONSTANTS;

typedef struct
{
    real r_point;
    real r_in_mag;
    real r_in_mag2;
    real qw_r3_N;
} R_POINTS;

typedef struct
{
    real nu;
    real id;
} NU_CONSTANTS;

typedef struct
{
    real dx;
    real qgaus_W;
} STREAM_GAUSS;


/* Parameter related types */

typedef struct
{
    real r_min, r_max, r_step_size;
    real nu_min, nu_max, nu_step_size;
    real mu_min, mu_max, mu_step_size;
    unsigned int r_steps, nu_steps, mu_steps;
} INTEGRAL_AREA;

typedef struct
{
    real weight;
    real step;
    real min;
    real max;
    int optimize;
} STREAM_WEIGHT;

#define EMPTY_STREAM_WEIGHT { 0.0, 0.0, 0.0, 0.0, 0 }

typedef struct
{
    real* stream_parameters;
    real* stream_step;
    real* stream_min;
    real* stream_max;
    int* stream_optimize;
} STREAM_PARAMETERS;

#define EMPTY_STREAM_PARAMETERS { NULL, NULL, NULL, NULL, NULL }

typedef struct
{
    STREAM_WEIGHT* stream_weight;
    STREAM_PARAMETERS* parameters;

    unsigned int number_streams;
    unsigned int number_stream_parameters;
} STREAMS;

#define EMPTY_STREAMS { NULL, NULL, 0, 0 }

typedef struct
{
    real* parameters;
    real* step;
    real* min;
    real* max;
    int* optimize;
} BACKGROUND_PARAMETERS;

#define EMPTY_BACKGROUND_PARAMETERS { NULL, NULL, NULL, NULL, NULL }


typedef struct
{
    real background_integral;
    real* stream_integrals;
} FINAL_STREAM_INTEGRALS;

#define EMPTY_FINAL_STREAM_INTEGRALS { 0.0, NULL }


typedef struct
{
    real parameters_version;
    real total_calc_probs;  /* sum of (r_steps * mu_steps * nu_steps) for all integrals */

    unsigned int number_background_parameters;
    real background_weight;

    unsigned int number_streams;
    unsigned int convolve;

    int sgr_coordinates;

    int aux_bg_profile;
    int wedge;

    unsigned int number_integrals;
    INTEGRAL_AREA* integral;

    /* Constants determined by other parameters */
    real alpha, q, sn, r0, delta, coeff, alpha_delta3;
    real bg_a, bg_b, bg_c;
} ASTRONOMY_PARAMETERS;


#define EMPTY_ASTRONOMY_PARAMETERS { 0.0, 0.0, \
                                     0, 0.0,   \
                                     0, 0,     \
                                     0, 0, 0, 0, NULL, \
                                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
                                     0.0, 0.0, 0.0 }

typedef struct
{
    real st_prob_int;    /* for Kahan summation */
    real st_prob_int_c;
} ST_PROBS;

#define ZERO_ST_PROBS = { 0.0, 0.0 }

typedef struct
{
    real bg_int;
    real correction;   /* Correction for Kahan summation */
} BG_PROB;

/* TODO: All these tuples of reals really serve the same
 * purpose. Fix having all of them. */
typedef struct
{
    real sum;
    real correction;
} PROB_SUM;

#define ZERO_PROB_SUM { 0.0, 0.0 }

#define CLEAR_BG_PROB(bgp) { (bgp).bg_int = 0.0; (bgp).correction = 0.0; }

/* Add b to a */
#define INCADD_BG_PROB(a, b) { (a).bg_int += (b).bg_int; (a).correction += (b).correction; }

#define ZERO_BG_PROB { 0.0, 0.0 }

typedef struct
{
    real st_only_sum;
    real st_only_sum_c;
} ST_SUM;


#ifdef __cplusplus
}
#endif

#endif /* _SEPARATION_TYPES_H_ */

