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
    real weight;
    real step;
    real min;
    real max;
    int optimize;
} StreamWeight;

#define EMPTY_STREAM_WEIGHT { 0.0, 0.0, 0.0, 0.0, 0 }

typedef struct
{
    real* stream_parameters;
    real* stream_step;
    real* stream_min;
    real* stream_max;
    int* stream_optimize;
} StreamParameters;

#define EMPTY_STREAM_PARAMETERS { NULL, NULL, NULL, NULL, NULL }

typedef struct
{
    StreamWeight* stream_weight;
    StreamParameters* parameters;

    real* expStreamWeights;
    real sumExpWeights;
    unsigned int number_streams;
    unsigned int number_stream_parameters;
} Streams;

#define EMPTY_STREAMS { NULL, NULL, NULL, 0.0, 0, 0 }

typedef struct
{
    real* parameters;
    real* step;
    real* min;
    real* max;
    int* optimize;
} BackgroundParameters;

#define EMPTY_BACKGROUND_PARAMETERS { NULL, NULL, NULL, NULL, NULL }

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

