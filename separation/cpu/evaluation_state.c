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

#define CHECKPOINT_FILE "astronomy_checkpoint"

#include "milkyway.h"
#include "milkyway_priv.h"


void initialize_integral(INTEGRAL* integral, unsigned int number_streams)
{
    integral->mu_step = 0;
    integral->nu_step = 0;
    integral->r_step = 0;

    integral->number_streams = number_streams;
    integral->background_integral = 0;
    integral->stream_integrals = calloc(number_streams, sizeof(double));
}

void initialize_state(const ASTRONOMY_PARAMETERS* ap, const STAR_POINTS* sp, EVALUATION_STATE* es)
{
    unsigned int i;

    es->current_integral = 0;
    es->background_integral = 0;
    es->stream_integrals = calloc(ap->number_streams, sizeof(double));

    es->number_streams = ap->number_streams;
    es->num_zero = 0;
    es->bad_jacobians = 0;

    es->number_integrals = ap->number_integrals;
    es->integrals = malloc(sizeof(INTEGRAL) * ap->number_integrals);

    for (i = 0; i < ap->number_integrals; i++)
        initialize_integral(&es->integrals[i], ap->number_streams);
}

void reset_evaluation_state(EVALUATION_STATE* es)
{
    unsigned int i, j;

    es->current_integral = 0;
    es->background_integral = 0;
    memset(es->stream_integrals, 0, sizeof(double) * es->number_streams);
    es->num_zero = 0;
    es->bad_jacobians = 0;

    for (i = 0; i < es->number_integrals; i++)
    {
        es->integrals[i].background_integral = 0;
        for (j = 0; j < es->integrals[i].number_streams; j++)
            es->integrals[i].stream_integrals[j] = 0;

        es->integrals[i].mu_step = 0;
        es->integrals[i].nu_step = 0;
        es->integrals[i].r_step = 0;
    }
}

void free_integral(INTEGRAL* i)
{
    free(i->stream_integrals);
}

void free_evaluation_state(EVALUATION_STATE* es)
{
    unsigned int i;

    free(es->stream_integrals);

    for (i = 0; i < es->number_integrals; ++i)
        free_integral(&es->integrals[i]);
    free(es->integrals);
}

#if BOINC_APPLICATION
static void fwrite_integral_area(FILE* file, INTEGRAL* ia)
{
    fprintf(file,
            "mu[min,max,steps]: %.3lf, %.3lf, %d\n",
            ia->mu_min,
            ia->mu_max,
            ia->mu_steps);

    fprintf(file,
            "nu[min,max,steps]: %.3lf, %.3lf, %d\n",
            ia->nu_min,
            ia->nu_max,
            ia->nu_steps);

    fprintf(file,
            " r[min,max,steps]: %.3lf, %.3lf, %d\n",
            ia->r_min,
            ia->r_max,
            ia->r_steps);

    fprintf(file,
            "mu_step: %d, nu_step: %d, r_step: %d\n",
	    ia->mu_step,
	    ia->nu_step,
	    ia->r_step);

    fprintf(file,
            "background_integral: %.20lf\n",
            ia->background_integral);

    fwrite_double_array(file,
                        "stream_integrals",
                        ia->stream_integrals,
                        ia->number_streams);
}

static void fread_integral_area(FILE* file, INTEGRAL* ia)
{
    unsigned int i;

    fscanf(file, "mu[min,max,steps]: %lf, %lf, %d\n", &ia->mu_min, &ia->mu_max, &ia->mu_steps);
    fscanf(file, "nu[min,max,steps]: %lf, %lf, %d\n", &ia->nu_min, &ia->nu_max, &ia->nu_steps);
    fscanf(file, " r[min,max,steps]: %lf, %lf, %d\n", &ia->r_min, &ia->r_max, &ia->r_steps);

    ia->mu_step_size = (ia->mu_max - ia->mu_min) / ia->mu_steps;
    ia->nu_step_size = (ia->nu_max - ia->nu_min) / ia->nu_steps;
    ia->r_step_size = (ia->r_max - ia->r_min) / ia->r_steps;


    fscanf(file, "mu_step: %d, nu_step: %d, r_step: %d\n", &(ia->mu_step), &(ia->nu_step), &(ia->r_step));

    fscanf(file, "background_integral: %lf\n", &ia->background_integral);
    fscanf(file, "stream_integrals[%d]: ", &ia->number_streams);

    for (i = 0; i < ia->number_streams; i++)
    {
        fscanf(file, "%lf", &ia->stream_integrals[i]);
        if (i != ia->number_streams - 1)
            fscanf(file, ", ");
    }
}

int read_checkpoint(EVALUATION_STATE* es)
{
    unsigned int i;
    char input_path[512];
    int retval = boinc_resolve_filename(CHECKPOINT_FILE, input_path, sizeof(input_path));
    if (retval)
        return 0;

    FILE* file = boinc_fopen(input_path, "r");
    if (file == NULL)
        return 0;

    if (1 > fscanf(file, "background_integral: %lf\n", &es->background_integral))
        return 1;

    es->stream_integrals = fread_double_array(file, "stream_integrals", &es->number_streams);

    fscanf(file, "prob_sum: %lf, num_zero: %d, bad_jacobians: %d\n",
           &es->prob_sum, &es->num_zero, &es->bad_jacobians);
    fscanf(file, "current_star_point: %d\n", &es->current_star_point);
    fscanf(file, "current_integral: %d\n", &es->current_integral);
    fscanf(file, "number_integrals: %d\n", &es->number_integrals);

    for (i = 0; i < es->number_integrals; i++)
        fread_integral_area(file, &es->integrals[i]);

    fclose(file);
    return 0;
}

int write_checkpoint(EVALUATION_STATE* es)
{
    int retval;
    unsigned int i;
    char output_path[512];
    FILE* file;

    boinc_resolve_filename(CHECKPOINT_FILE, output_path, sizeof(output_path));

    file = boinc_fopen(output_path, "w+");
    if (!file)
    {
        fprintf(stderr, "APP: error writing checkpoint (opening checkpoint file)\n");
        return 1;
    }

    fprintf(file, "background_integral: %.20lf\n", es->background_integral);
    fprintf(file, "stream_integrals[%d]: ", es->number_streams);

    for (i = 0; i < es->number_streams; i++)
    {
        fprintf(file, "%.20lf", es->stream_integrals[i]);
        if (i != (es->number_streams - 1))
            fprintf(file, ", ");
    }
    fprintf(file, "\n");

    fprintf(file, "prob_sum: %.20lf, num_zero: %d, bad_jacobians: %d\n", es->prob_sum, es->num_zero, es->bad_jacobians);
    fprintf(file, "current_star_point: %d\n", es->current_star_point);
    fprintf(file, "current_integral: %d\n", es->current_integral);
    fprintf(file, "number_integrals: %d\n", es->number_integrals);

    for (i = 0; i < es->number_integrals; i++)
    {
        fwrite_integral_area(file, &es->integrals[i]);
    }

    if ((retval = fclose(file)))
    {
        fprintf(stderr, "APP: error writing checkpoint (closing checkpoint file) %d\n", retval);
        return retval;
    }

    return 0;
}

#endif /* BOINC_APPLICATION */


