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

#ifndef _SEPARATION_KERNEL_TYPES_H_
#define _SEPARATION_KERNEL_TYPES_H_

/* Should only be directly included by CL kernel */

#ifndef _MSC_VER
  #define SEPARATION_ALIGN(x) __attribute__ ((aligned(x)))
#else
  #define SEPARATION_ALIGN(x) __declspec(align(x))
#endif /* _MSC_VER */


#ifdef __OPENCL_VERSION__

#ifdef cl_amd_fp64
  #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
  #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif /* cl_amd_fp64 */


#if DOUBLEPREC
typedef double real;
typedef double2 real2;
typedef double4 real4;
#else
typedef float real;
typedef float2 real2;
typedef float4 real4;
#endif /* DOUBLEPREC */

typedef struct
{
    real x, y, z, w;
} mwvector;

#define X(v) ((v).x)
#define Y(v) ((v).y)
#define Z(v) ((v).z)
#define W(v) ((v).w)
#endif /* __OPENCL_VERSION__ */



#ifndef __OPENCL_VERSION__
typedef struct
{
    real irv_reff_xr_rp3;
    real gPrime;
} RConsts;

  #define IRV_REFF_XR_RP3(r) ((r).irv_reff_xr_rp3)
  #define GPRIME(r) ((r).gPrime)
#else

typedef real2 RConsts;

  #define IRV_REFF_XR_RP3(r) ((r).x)
  #define GPRIME(r) ((r).y)
#endif



typedef struct SEPARATION_ALIGN(128)
{
    mwvector a;
    mwvector c;
    real sigma_sq2_inv;
    int large_sigma;          /* abs(stream_sigma) > SIGMA_LIMIT */
} StreamConstants;

#ifndef __OPENCL_VERSION__
typedef struct
{
    real r_point;
    real qw_r3_N;
} RPoints;

#define R_POINT(r) ((r).r_point)
#define QW_R3_N(r) ((r).qw_r3_N)

typedef struct
{
    real lCosBCos;
    real lSinBCos;
    real bSin;
    real _pad;
} LBTrig;

#define LCOS_BCOS(x) ((x).lCosBCos)
#define LSIN_BCOS(x) ((x).lSinBCos)
#define BSIN(x) ((x).bSin)

#else

/* x = r_point; y = qw_r3_N */
typedef real2 RPoints;

#define R_POINT(r) ((r).x)
#define QW_R3_N(r) ((r).y)


typedef real4 LBTrig;

#define LCOS_BCOS(l) ((l).x)
#define LSIN_BCOS(l) ((l).y)
#define BSIN(l) ((l).z)

#endif /* __OPENCL_VERSION__ */


/* Parameter related types */

typedef struct SEPARATION_ALIGN(128)
{
    real r_min, r_max, r_step_size;
    real nu_min, nu_max, nu_step_size;
    real mu_min, mu_max, mu_step_size;
    unsigned int r_steps, nu_steps, mu_steps;
} IntegralArea;


/* Kitchen sink of constants, etc. */
typedef struct SEPARATION_ALIGN(128)
{
    /* Constants determined by other parameters */
    real m_sun_r0;
    real q_inv;
    real q_inv_sqr;  /* 1 / q^2 */
    real r0;

    unsigned int convolve;
    unsigned int number_streams;
    int fast_h_prob;
    int aux_bg_profile;

    real alpha;
    real delta;
    real alpha_delta3;
    real bg_a, bg_b, bg_c;

    int wedge;
    real sun_r0;
    real q;
    real coeff;

    real total_calc_probs;  /* sum of (r_steps * mu_steps * nu_steps) for all integrals */
    unsigned int number_integrals;

    real exp_background_weight;
} AstronomyParameters;

typedef struct SEPARATION_ALIGN(16)
{
    real sum;
    real correction;
} Kahan;

#define ZERO_KAHAN { 0.0, 0.0 }

#define CLEAR_KAHAN(k) { (k).sum = 0.0; (k).correction = 0.0; }


#endif /* _SEPARATION_KERNEL_TYPES_H_ */

