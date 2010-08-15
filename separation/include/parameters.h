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

#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <stdio.h>
#include "coordinates.h"

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

#define EMPTY_STREAM_PARAMETERS { NULL, NULL, NULL, NULL, NULL }

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

void free_astronomy_parameters(ASTRONOMY_PARAMETERS* ap);
void free_background_parameters(BACKGROUND_PARAMETERS* bgp);
void free_streams(STREAMS* streams);

unsigned int get_optimized_parameter_count(ASTRONOMY_PARAMETERS* ap,
                                           BACKGROUND_PARAMETERS* bgp,
                                           STREAMS* streams);

int read_parameters(const char* file,
                    ASTRONOMY_PARAMETERS* ap,
                    BACKGROUND_PARAMETERS* bgp,
                    STREAMS* streams);

void fread_parameters(FILE* file,
                      ASTRONOMY_PARAMETERS* ap,
                      BACKGROUND_PARAMETERS* bgp,
                      STREAMS* streams);

int write_parameters(const char* file,
                     ASTRONOMY_PARAMETERS* ap,
                     BACKGROUND_PARAMETERS* bgp,
                     STREAMS* streams);

void fwrite_parameters(FILE* file,
                       ASTRONOMY_PARAMETERS* ap,
                       BACKGROUND_PARAMETERS* bgp,
                       STREAMS* streams);

void set_parameters(ASTRONOMY_PARAMETERS* ap,
                    BACKGROUND_PARAMETERS* bgp,
                    STREAMS* streams,
                    double* parameters);


#endif /* _PARAMETERS_H_ */

