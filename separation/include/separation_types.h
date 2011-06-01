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

#include "separation_config.h"
#include "milkyway_math.h"
#include "separation_kernel_types.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    real l;
    real b;
} LB;

#define LB_L(x) ((x).l)
#define LB_B(x) ((x).b)



typedef struct
{
    real nu;
    real id;
} NuId;


typedef struct
{
    unsigned int number_stars;
    mwvector* stars;
} StarPoints;

#define EMPTY_STAR_POINTS { 0, NULL }


/* Convenience structure for passing mess of LBTrig to CAL kernel in 2 parts */
typedef struct
{
    real lCosBCos, lSinBCos;
} LTrigPair;

typedef struct
{
    real nu;
    real id;
} NuConstants;

typedef struct
{
    real* dx;
    real* qgaus_W;
} StreamGauss;


/* Parameter related types */


typedef struct
{
    real epsilon;    /* Stream weight */
    real mu;         /* mu, angular position */
    real r;          /* r, radial distance from the sun in kpc */
    real theta;      /* theta, zenith angle from Z-axis */
    real phi;        /* phi, azimuthal angle from X-axis */
    real sigma;      /* sigma, standard deviation of the Gaussian defining the stream falloff in kpc */

    real epsilonExp; /* exp(epsilon) */
} StreamParameters;

#define EMPTY_STREAM_PARAMETERS { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }

typedef struct
{
    StreamParameters* parameters;
    real sumExpWeights;
    unsigned int number_streams;
} Streams;

#define EMPTY_STREAMS { NULL, 0.0, 0 }

typedef struct
{
    real alpha;
    real r0;
    real q;
    real delta;

    real epsilon;

    real a, b, c;
} BackgroundParameters;

#define EMPTY_BACKGROUND_PARAMETERS { 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 }

typedef struct
{
    real irv;
    real rPrime;
} RPrime;




typedef struct
{
    real likelihood;

    real backgroundIntegral;
    real backgroundLikelihood;

    real* streamIntegrals;

    real* streamLikelihoods;
} SeparationResults;


#ifdef __cplusplus
}
#endif

#endif /* _SEPARATION_TYPES_H_ */

