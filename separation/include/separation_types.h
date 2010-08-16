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

typedef struct
{
    unsigned int number_stars;
    double* stars;
} STAR_POINTS;

#define EMPTY_STAR_POINTS { 0, NULL }


typedef struct
{
    vector a;
    vector c;
    double sigma_sq2;
    int large_sigma;          /* abs(stream_sigma) > SIGMA_LIMIT */
} STREAM_CONSTANTS;

typedef struct
{
    double irv;
    double reff_xr_rp3;
} R_CONSTANTS;

typedef struct
{
    double r_point;
    double r_in_mag;
    double r_in_mag2;
    double qw_r3_N;
} R_POINTS;

typedef struct
{
    double nu;
    double id;
} NU_CONSTANTS;

typedef struct
{
    R_CONSTANTS* r_step_consts;
    R_POINTS* r_pts;
    NU_CONSTANTS* nu_consts;
} INTEGRAL_CONSTANTS;


typedef struct
{
    double* dx;
    double* qgaus_W;
} STREAM_GAUSS;


/* Parameter related types */

typedef struct
{
    double r_min, r_max, r_step_size;
    double nu_min, nu_max, nu_step_size;
    double mu_min, mu_max, mu_step_size;
    unsigned int r_steps, nu_steps, mu_steps;
} INTEGRAL_AREA;

typedef struct
{
    double weight;
    double step;
    double min;
    double max;
    int optimize;
} STREAM_WEIGHT;

#define EMPTY_STREAM_WEIGHT { 0.0, 0.0, 0.0, 0.0, 0 }

typedef struct
{
    double* stream_parameters;
    double* stream_step;
    double* stream_min;
    double* stream_max;
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
    double* parameters;
    double* step;
    double* min;
    double* max;
    int* optimize;
} BACKGROUND_PARAMETERS;

#define EMPTY_BACKGROUND_PARAMETERS { NULL, NULL, NULL, NULL, NULL }

#ifndef __OPENCL_VERSION__
/* No function pointers allowed in kernels, but we don't need it. */
/* wedge, mu, nu, l, b. gc2lb or gc2sgr  */
typedef void (*SGRConversion)(int, double, double, double*, double*);

#else

typedef void* SGRConversion;

#endif


typedef struct
{
    double parameters_version;
    double total_calc_probs;  /* sum of (r_steps * mu_steps * nu_steps) for all integrals */

    unsigned int number_background_parameters;
    double background_weight;

    unsigned int number_streams;
    unsigned int convolve;

    unsigned int sgr_coordinates;
    SGRConversion sgr_conversion;

    int aux_bg_profile;
    int wedge;

    unsigned int number_integrals;
    INTEGRAL_AREA* integral;

    /* Constants determined by other parameters */
    double alpha, q, sn, r0, delta, coeff, alpha_delta3;
    double bg_a, bg_b, bg_c;
} ASTRONOMY_PARAMETERS;


#define EMPTY_ASTRONOMY_PARAMETERS { 0.0, 0.0, \
                                     0, 0.0,   \
                                     0, 0,     \
                                     0, NULL, 0, 0, 0, NULL, \
                                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
                                     0.0, 0.0, 0.0 }



typedef struct
{
    double st_prob_int;    /* for Kahan summation */
    double st_prob_int_c;
} ST_PROBS;

typedef struct
{
    double bg_int;
    double correction;   /* Correction for Kahan summation */
} BG_PROB;

/* TODO: All these tuples of doubles really serve the same
 * purpose. Fix having all of them. */
typedef struct
{
    double sum;
    double correction;
} PROB_SUM;

#define ZERO_PROB_SUM { 0.0, 0.0 }

#define CLEAR_BG_PROB(bgp) { (bgp).bg_int = 0.0; (bgp).correction = 0.0; }

/* Add b to a */
#define INCADD_BG_PROB(a, b) { (a).bg_int += (b).bg_int; (a).correction += (b).correction; }

#define ZERO_BG_PROB { 0.0, 0.0 }

typedef struct
{
    double st_only_sum;
    double st_only_sum_c;
} ST_SUM;


#ifdef __cplusplus
}
#endif

#endif /* _SEPARATION_TYPES_H_ */

