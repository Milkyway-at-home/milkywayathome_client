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

#include "separation_config.h"
#include "milkyway_math.h"
#include "milkyway_types.h"

#ifndef _MSC_VER
  #define SEPARATION_ALIGN(x) __attribute__ ((packed, aligned(x)))
#else
  #define SEPARATION_ALIGN(x)
#endif /* _MSC_VER */

typedef struct SEPARATION_ALIGN(2 * sizeof(real))
{
    real l;
    real b;
} LB;

typedef struct SEPARATION_ALIGN(4 * sizeof(real))
{
    real lsin, lcos;
    real bsin, bcos;
} LB_TRIG;

#define LB_L(x) ((x).l)
#define LB_B(x) ((x).b)

typedef struct SEPARATION_ALIGN(4 * sizeof(real))
{
    real reff_xr_rp3;
    real irv;
    real gPrime;
} R_CONSTS;

typedef struct SEPARATION_ALIGN(2 * sizeof(real))
{
    real nu;
    real id;
} NU_ID;


typedef struct
{
    unsigned int number_stars;
    vector* stars;
} STAR_POINTS;

#define EMPTY_STAR_POINTS { 0, NULL }

typedef struct SEPARATION_ALIGN(64)
{
    vector a;
    vector c;
    real sigma_sq2_inv;
    mw_int large_sigma;          /* abs(stream_sigma) > SIGMA_LIMIT */
} STREAM_CONSTANTS;

typedef struct SEPARATION_ALIGN(4 * sizeof(real))
{
    real r_point;
    real r_in_mag;
    real r_in_mag2;
    real qw_r3_N;
} R_POINTS;

typedef struct SEPARATION_ALIGN(2 * sizeof(real))
{
    real nu;
    real id;
} NU_CONSTANTS;

typedef struct SEPARATION_ALIGN(2 * sizeof(real))
{
    real dx;
    real qgaus_W;
} STREAM_GAUSS;


/* Parameter related types */

typedef struct SEPARATION_ALIGN(128)
{
    real r_min, r_max, r_step_size;
    real nu_min, nu_max, nu_step_size;
    real mu_min, mu_max, mu_step_size;
    mw_uint r_steps, nu_steps, mu_steps;
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

typedef struct SEPARATION_ALIGN(2 * sizeof(real))
{
    real irv;
    real rPrime;
} R_PRIME;

typedef struct SEPARATION_ALIGN(128)
{
    real parameters_version;
    real total_calc_probs;  /* sum of (r_steps * mu_steps * nu_steps) for all integrals */
    mw_int sgr_coordinates;
    mw_uint number_integrals;
    mw_uint number_background_parameters;
    mw_uint number_streams;
    mw_uint convolve;

    mw_int aux_bg_profile;
    mw_int fast_h_prob;
    mw_int zero_q;
    mw_int wedge;
    real background_weight;

    /* Constants determined by other parameters */
    real alpha, sn, r0, delta, coeff, alpha_delta3;
    real q;
    real q_inv_sqr;  /* 1 / q^2 */
    real bg_a, bg_b, bg_c;
} ASTRONOMY_PARAMETERS;


#define EMPTY_ASTRONOMY_PARAMETERS { 0.0, 0.0, \
                                     0, 0,   \
                                     0, 0, 0, 0, \
                                     0, 0, 0, 0.0, 0.0,   \
                                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
                                     0.0, 0.0, 0.0 }

//    typedef struct __attribute__((packed, aligned(4 * sizeof(real)), vec_type_hint(float2)))  //SEPARATION_ALIGN(2 * sizeof(real))

typedef struct SEPARATION_ALIGN(2 * sizeof(real))
{
    real sum;
    real correction;
} KAHAN;

#define ZERO_KAHAN { 0.0, 0.0 }

#define CLEAR_KAHAN(k) { (k).sum = 0.0; (k).correction = 0.0; }

#ifdef __cplusplus
}
#endif

#endif /* _SEPARATION_TYPES_H_ */

