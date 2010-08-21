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

#include <stdio.h>
#include <stdlib.h>

#include "io_util.h"
#include "milkyway_util.h"
#include "separation_types.h"

/****
	*	Functions for printing parameters to files/output
*****/
void fwrite_double_array(FILE *file, const char *array_name, real* array_t, size_t size)
{
	size_t i;

	fprintf(file, "%s[%u]: ", array_name, size);
	for (i = 0; i < size; i++)
    {
		fprintf(file, "%.20lf", array_t[i]);
		if (i < size - 1)
            fprintf(file, ", ");
	}
	fprintf(file, "\n");
}

void fwrite_int_array(FILE* file, const char* array_name, int* array_t, size_t size)
{
	size_t i;
	fprintf(file, "%s[%u]: ", array_name, size);
	for (i = 0; i < size; i++)
    {
		if (i == 0)
            fprintf(file, " %d", array_t[i]);
		else
            fprintf(file, ", %d", array_t[i]);
	}
	fprintf(file, "\n");
}

#if DOUBLEPREC
  #define READ_DOUBLE_ARRAY_READ_STR "%lf"
#else
  #define READ_DOUBLE_ARRAY_READ_STR "%f"
#endif /* DOUBLEPREC */

/*Functions for reading parameters from files */
real* fread_double_array(FILE* file, const char* array_name, unsigned int* sizeOut)
{
    unsigned int i, size;
    int rc;
    real* arr;

	fscanf(file, array_name);
	fscanf(file, "[%u]: ", &size);

    arr = mallocSafe(sizeof(real) * size);

	for (i = 0; i < size; i++)
    {
        rc = fscanf(file, READ_DOUBLE_ARRAY_READ_STR, &arr[i]);
		if (rc != 1)
			fail("Error reading into %s\n", array_name);
		if (i < size - 1)
            fscanf(file, ", ");
	}
	fscanf(file, "\n");

    if (sizeOut)
        *sizeOut = size;

	return arr;
}

int* fread_int_array(FILE *file, const char *array_name, unsigned int* sizeOut)
{
	unsigned int i, size;
	fscanf(file, array_name);
	fscanf(file, "[%u]: ", &size);
    int* arr;

	arr = mallocSafe(sizeof(int) * size);

	for (i = 0; i < size; i++)
    {
		if (fscanf(file, "%d", &arr[i]) != 1)
        {
			fprintf(stderr,"Error reading into %s\n", array_name);
			exit(-1);
		}

		if (i < size - 1)
            fscanf(file, ", ");
	}
	fscanf(file, "\n");

    if (sizeOut)
        *sizeOut = size;

	return arr;
}

void printIntegralArea(INTEGRAL_AREA* ia)
{
    printf("integral-area {\n"
           "  r_min        = %g\n"
           "  r_max        = %g\n"
           "  r_step_size  = %g\n"
           "  nu_min       = %g\n"
           "  nu_max       = %g\n"
           "  nu_step_size = %g\n"
           "  mu_min       = %g\n"
           "  mu_max       = %g\n"
           "  mu_step_size = %g\n"
           "  r_steps      = %u\n"
           "  nu_steps     = %u\n"
           "  mu_steps     = %u\n",
           ia->r_min, ia->r_max, ia->r_step_size,
           ia->nu_min, ia->nu_max, ia->nu_step_size,
           ia->mu_min, ia->mu_max, ia->mu_step_size,
           ia->r_steps, ia->nu_steps, ia->mu_steps);
}

void printAstronomyParameters(ASTRONOMY_PARAMETERS* ap)
{
    unsigned int i;
    printf("astronomy-parameters {\n"
           "  parameters_version           = %g\n"
           "  total_calc_probs             = %g\n"
           "  number_background_parameters = %u\n"
           "  background_weight            = %g\n"
           "  number_streams               = %u\n"
           "  convolve                     = %u\n"
           "  sgr_coordinates              = %u\n"
           "  aux_bg_profile               = %d\n"
           "  wedge                        = %d\n"
           "  number_integrals             = %u\n"
           "  integral                     = { integrals } %p\n"
           "  alpha                        = %g\n"
           "  q                            = %g\n"
           "  sn                           = %g\n"
           "  r0                           = %g\n"
           "  delta                        = %g\n"
           "  coeff                        = %g\n"
           "  alpha_delta                  = %g\n"
           "  bg_a                         = %g\n"
           "  bg_b                         = %g\n"
           "  bg_c                         = %g\n",
           ap->parameters_version,
           ap->total_calc_probs,
           ap->number_background_parameters,
           ap->background_weight,
           ap->number_streams,
           ap->convolve,
           ap->sgr_coordinates,
           ap->aux_bg_profile,
           ap->wedge,
           ap->number_integrals,
           ap->integral,
           ap->alpha, ap->q, ap->sn, ap->r0, ap->delta, ap->coeff, ap->alpha_delta3,
           ap->bg_a, ap->bg_b, ap->bg_c);

    printf("ap integral: \n");
    for (i = 0; i < ap->number_integrals; ++i)
    {
        printf("integral area[%u]]\n", i);
        printIntegralArea(&ap->integral[i]);
    }
}

