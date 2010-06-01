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
#include <stdlib.h>

typedef struct integral
{
    double r_min, r_max, r_step_size;
    double nu_min, nu_max, nu_step_size;
    double mu_min, mu_max, mu_step_size;
    int r_steps, nu_steps, mu_steps;

    long min_calculation, max_calculation;
} INTEGRAL;

typedef struct astronomy_parameters
{
    double parameters_version;

    int number_background_parameters;
    double background_weight;
    double* background_parameters;
    double* background_step;
    double* background_min;
    double* background_max;
    int* background_optimize;

    int number_streams;
    int number_stream_parameters;

    double* stream_weights;
    double* stream_weight_step;
    double* stream_weight_min;
    double* stream_weight_max;
    int* stream_weight_optimize;

    double** stream_parameters;
    double** stream_step;
    double** stream_min;
    double** stream_max;
    int** stream_optimize;

    int convolve;
    int sgr_coordinates;
    int aux_bg_profile; //vickej2_bg
    int wedge;

    int number_integrals;
    INTEGRAL** integral;
} ASTRONOMY_PARAMETERS;

int get_optimized_parameter_count(ASTRONOMY_PARAMETERS* ap);

int read_astronomy_parameters(const char* file, ASTRONOMY_PARAMETERS* ap);
void    fread_astronomy_parameters(FILE* file, ASTRONOMY_PARAMETERS* ap);
int write_astronomy_parameters(const char* file, ASTRONOMY_PARAMETERS* ap);
void    fwrite_astronomy_parameters(FILE* file, ASTRONOMY_PARAMETERS* ap);
void    free_parameters(ASTRONOMY_PARAMETERS* ap);

void    split_astronomy_parameters(ASTRONOMY_PARAMETERS* ap, int rank, int max_rank);

void    set_astronomy_parameters(ASTRONOMY_PARAMETERS* ap, double* parameters);
void    get_search_parameters(ASTRONOMY_PARAMETERS* ap, double** parameters);
void    get_min_parameters(ASTRONOMY_PARAMETERS* ap, double** parameters);
void    get_max_parameters(ASTRONOMY_PARAMETERS* ap, double** parameters);
void    get_step(ASTRONOMY_PARAMETERS* ap, double** step);

#endif

