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

#ifndef ASTRONOMY_PARAMETERS_H
#define ASTRONOMY_PARAMETERS_H

#include <stdio.h>

typedef struct
{
    double r_min, r_max, r_step_size;
    double nu_min, nu_max, nu_step_size;
    double mu_min, mu_max, mu_step_size;
    int r_steps, nu_steps, mu_steps;

    long min_calculation, max_calculation;
} INTEGRAL;

typedef struct
{
    double weights;
    double weight_step;
    double weight_min;
    double weight_max;
    int weight_optimize;
} STREAM;

#define EMPTY_STREAM { 0.0, 0.0, 0.0, 0.0, 0 }

typedef struct
{
    double* stream_parameters;
    double* stream_step;
    double* stream_min;
    double* stream_max;
    int* stream_optimize;
} STREAM_PARAMETERS;

#define EMPTY_STREAM_PARAMETERS { NULL, NULL, NULL, NULL, NULL }

#define STREAM_N(ap, n) ((ap)->stream[(n)])
#define STREAM_PARAM_N(ap, n) ((ap)->parameters[(n)])

typedef struct
{
    double parameters_version;

    int number_background_parameters;
    double background_weight;
    double* background_parameters;
    double* background_step;
    double* background_min;
    double* background_max;
    int* background_optimize;

    unsigned int number_streams;
    unsigned int number_stream_parameters;

    STREAM* stream;
    STREAM_PARAMETERS* parameters;

    int convolve;
    int sgr_coordinates;
    int aux_bg_profile;
    int wedge;

    int number_integrals;
    INTEGRAL** integral;
} ASTRONOMY_PARAMETERS;

#define EMPTY_ASTRONOMY_PARAMETERS { 0.0, 0, 0.0, NULL, NULL, NULL, NULL, NULL, \
                                     0, 0, NULL, NULL, \
                                     0, 0, 0, 0, 0, NULL }


int get_optimized_parameter_count(ASTRONOMY_PARAMETERS* ap);

int read_astronomy_parameters(const char* file, ASTRONOMY_PARAMETERS* ap);
void fread_astronomy_parameters(FILE* file, ASTRONOMY_PARAMETERS* ap);
int write_astronomy_parameters(const char* file, ASTRONOMY_PARAMETERS* ap);
void fwrite_astronomy_parameters(FILE* file, ASTRONOMY_PARAMETERS* ap);
void free_parameters(ASTRONOMY_PARAMETERS* ap);

void set_astronomy_parameters(ASTRONOMY_PARAMETERS* ap, double* parameters);
void get_step(ASTRONOMY_PARAMETERS* ap, double** step);

#endif

