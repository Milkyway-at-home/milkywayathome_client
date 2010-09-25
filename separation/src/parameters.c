/*
Copyright 2008-2010 Travis Desell, Dave Przybylo, Nathan Cole,
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

#include <stdlib.h>
#include <stdio.h>

#include "separation_types.h"
#include "parameters.h"
#include "io_util.h"
#include "milkyway_util.h"

void free_background_parameters(BACKGROUND_PARAMETERS* bgp)
{
    free(bgp->parameters);
    free(bgp->step);
    free(bgp->min);
    free(bgp->max);
    free(bgp->optimize);
}

void free_stream_parameters(STREAM_PARAMETERS* p)
{
    free(p->stream_parameters);
    free(p->stream_step);
    free(p->stream_min);
    free(p->stream_max);
    free(p->stream_optimize);
}

void free_streams(STREAMS* streams)
{
    unsigned int i;

    for (i = 0; i < streams->number_streams; ++i)
        free_stream_parameters(&streams->parameters[i]);
    free(streams->parameters);
    free(streams->stream_weight);
}

static void calc_integral_step_sizes(INTEGRAL_AREA* i)
{
    i->r_step_size = (i->r_max - i->r_min) / (real)i->r_steps;
    i->mu_step_size = (i->mu_max - i->mu_min) / (real)i->mu_steps;
    i->nu_step_size = (i->nu_max - i->nu_min) / (real)i->nu_steps;
}

static INTEGRAL_AREA* fread_parameters(FILE* file,
                                       ASTRONOMY_PARAMETERS* ap,
                                       BACKGROUND_PARAMETERS* bgp,
                                       STREAMS* streams)
{
    unsigned int i, temp;
    int retval;
    double tmp1, tmp2;
    INTEGRAL_AREA* integral;
    const INTEGRAL_AREA* ia;
    unsigned int total_calc_probs;

    retval = fscanf(file, "parameters_version: %lf\n", &tmp1);
    ap->parameters_version = (real) tmp1;
    if (retval < 1)
    {
        ap->parameters_version = 0.01;
//      fprintf(stderr, "Error reading astronomy parameters file. Parameters version not specified\n");
//      return;
    }

    fscanf(file, "number_parameters: %u\n", &ap->number_background_parameters);
    fscanf(file, "background_weight: %lf\n", &tmp1);
    ap->background_weight = (real) tmp1;

    bgp->parameters = fread_double_array(file, "background_parameters", NULL);
    bgp->step       = fread_double_array(file, "background_step", NULL);
    bgp->min        = fread_double_array(file, "background_min", NULL);
    bgp->max        = fread_double_array(file, "background_max", NULL);
    bgp->optimize   = fread_int_array(file, "optimize_parameter", NULL);

    fscanf(file, "number_streams: %u, %u\n", &streams->number_streams, &streams->number_stream_parameters);

    ap->number_streams = streams->number_streams;

    streams->stream_weight = (STREAM_WEIGHT*) mallocSafe(sizeof(STREAM_WEIGHT) * streams->number_streams);
    streams->parameters = (STREAM_PARAMETERS*) mallocSafe(sizeof(STREAM_PARAMETERS) * streams->number_streams);

    for (i = 0; i < streams->number_streams; ++i)
    {
        fscanf(file, "stream_weight: %lf\n", &tmp1);
        streams->stream_weight[i].weight = (real) tmp1;
        fscanf(file, "stream_weight_step: %lf\n", &tmp1);
        streams->stream_weight[i].step = (real) tmp1;
        fscanf(file, "stream_weight_min: %lf\n", &tmp1);
        streams->stream_weight[i].min = (real) tmp1;
        fscanf(file, "stream_weight_max: %lf\n", &tmp1);
        streams->stream_weight[i].max = (real) tmp1;
        fscanf(file, "optimize_weight: %d\n", &streams->stream_weight[i].optimize);

        streams->parameters[i].stream_parameters = fread_double_array(file, "stream_parameters", NULL);
        streams->parameters[i].stream_step       = fread_double_array(file, "stream_step", NULL);
        streams->parameters[i].stream_min        = fread_double_array(file, "stream_min", NULL);
        streams->parameters[i].stream_max        = fread_double_array(file, "stream_max", NULL);
        streams->parameters[i].stream_optimize   = fread_int_array(file, "optimize_parameter", NULL);
    }

    if (fscanf(file, "convolve: %u\n", &ap->convolve) < 1)
        warn("Error reading convolve\n");

    fscanf(file, "sgr_coordinates: %d\n", &ap->sgr_coordinates);
    if (ap->parameters_version > 0.01)
    {
        if (fscanf(file, "aux_bg_profile: %d\n", &ap->aux_bg_profile) < 1)
            warn("Error reading aux_bg_profile\n");
    }

    if (fscanf(file, "wedge: %d\n", &ap->wedge) < 1)
        warn("Error reading wedge\n");

    //integral = (INTEGRAL_AREA*) mallocSafe(sizeof(INTEGRAL_AREA));
    integral = (INTEGRAL_AREA*) mwMallocAligned(sizeof(INTEGRAL_AREA), sizeof(INTEGRAL_AREA));

    fscanf(file,
           "r[min,max,steps]: %lf, %lf, %u\n",
           &tmp1, &tmp2, &integral[0].r_steps);

    integral[0].r_min = (real) tmp1;
    integral[0].r_max = (real) tmp2;

    fscanf(file,
           "mu[min,max,steps]: %lf, %lf, %u\n",
           &tmp1, &tmp2, &integral[0].mu_steps);

    integral[0].mu_min = (real) tmp1;
    integral[0].mu_max = (real) tmp2;


    fscanf(file,
           "nu[min,max,steps]: %lf, %lf, %u\n",
           &tmp1, &tmp2, &integral[0].nu_steps);

    integral[0].nu_min = (real) tmp1;
    integral[0].nu_max = (real) tmp2;

    calc_integral_step_sizes(&integral[0]);

    fscanf(file, "number_cuts: %u\n", &ap->number_integrals);
    ap->number_integrals++;
    if (ap->number_integrals > 1)
    {
        /* FIXME?: Alignment in case where copy happens */
        integral = (INTEGRAL_AREA*) reallocSafe(integral, sizeof(INTEGRAL_AREA) * ap->number_integrals);
        for (i = 1; i < ap->number_integrals; i++)
        {
            fscanf(file,
                   "r_cut[min,max,steps][%u]: %lf, %lf, %u\n",
                   &temp, &tmp1, &tmp2, &integral[i].r_steps);

            integral[i].r_min = (real) tmp1;
            integral[i].r_max = (real) tmp2;

            fscanf(file,
                   "mu_cut[min,max,steps][%u]: %lf, %lf, %u\n",
                   &temp, &tmp1, &tmp2, &integral[i].mu_steps);

            integral[i].mu_min = (real) tmp1;
            integral[i].mu_max = (real) tmp2;


            fscanf(file,
                   "nu_cut[min,max,steps][%u]: %lf, %lf, %u\n",
                   &temp, &tmp1, &tmp2, &integral[i].nu_steps);

            integral[i].nu_min = (real) tmp1;
            integral[i].nu_max = (real) tmp2;


            calc_integral_step_sizes(&integral[i]);
        }
    }

    /* Calculate total probability calculations for checkpointing */

    total_calc_probs = 0;
    for (i = 0; i < ap->number_integrals; ++i)
    {
        ia = &integral[i];
        total_calc_probs += ia->mu_steps * ia->nu_steps * ia->r_steps;
    }

    ap->total_calc_probs = (real) total_calc_probs;
    return integral;
}

INTEGRAL_AREA* read_parameters(const char* filename,
                               ASTRONOMY_PARAMETERS* ap,
                               BACKGROUND_PARAMETERS* bgp,
                               STREAMS* streams)
{

    FILE* f;
    INTEGRAL_AREA* integral;

    f = mwOpenResolved(filename, "r");
    if (!f)
    {
        perror("Opening astronomy parameters file");
        return NULL;
    }

    integral = fread_parameters(f, ap, bgp, streams);
    if (!integral)
        warn("Error reading parameters file\n");

    fclose(f);
    return integral;
}

static void fwrite_parameters(FILE* file,
                              ASTRONOMY_PARAMETERS* ap,
                              INTEGRAL_AREA* integral,
                              BACKGROUND_PARAMETERS* bgp,
                              STREAMS* streams)
{
    unsigned int i;

    fprintf(file, "parameters_version: %lf\n", ap->parameters_version);

    fprintf(file, "number_parameters: %u\n", ap->number_background_parameters);
    fprintf(file, "background_weight: %lf\n", ap->background_weight);
    fwrite_double_array(file, "background_parameters", bgp->parameters, ap->number_background_parameters);
    fwrite_double_array(file, "background_step", bgp->step, ap->number_background_parameters);
    fwrite_double_array(file, "background_min", bgp->min, ap->number_background_parameters);
    fwrite_double_array(file, "background_max", bgp->max, ap->number_background_parameters);
    fwrite_int_array(file, "optimize_parameter", bgp->optimize, ap->number_background_parameters);

    fprintf(file, "number_streams: %u, %u\n", streams->number_streams, streams->number_stream_parameters);
    for (i = 0; i < streams->number_streams; i++)
    {
        fprintf(file, "stream_weight: %lf\n", streams->stream_weight[i].weight);
        fprintf(file, "stream_weight_step: %lf\n", streams->stream_weight[i].step);
        fprintf(file, "stream_weight_min: %lf\n", streams->stream_weight[i].min);
        fprintf(file, "stream_weight_max: %lf\n", streams->stream_weight[i].max);
        fprintf(file, "optimize_weight: %d\n", streams->stream_weight[i].optimize);

        fwrite_double_array(file, "stream_parameters",
                            streams->parameters[i].stream_parameters, streams->number_stream_parameters);
        fwrite_double_array(file, "stream_step",
                            streams->parameters[i].stream_step, streams->number_stream_parameters);
        fwrite_double_array(file, "stream_min",
                            streams->parameters[i].stream_min, streams->number_stream_parameters);
        fwrite_double_array(file, "stream_max",
                            streams->parameters[i].stream_max, streams->number_stream_parameters);
        fwrite_int_array(file, "optimize_parameter",
                         streams->parameters[i].stream_optimize, streams->number_stream_parameters);
    }

    fprintf(file, "convolve: %d\n", ap->convolve);
    fprintf(file, "sgr_coordinates: %d\n", ap->sgr_coordinates);
    fprintf(file, "aux_bg_profile: %d\n", ap->aux_bg_profile); //vickej2_bg
    fprintf(file, "wedge: %d\n", ap->wedge);

    fprintf(file,
            "r[min,max,steps]: %lf, %lf, %u\n",
            integral[0].r_min,
            integral[0].r_max,
            integral[0].r_steps);

    fprintf(file,
            "mu[min,max,steps]: %lf, %lf, %u\n",
            integral[0].mu_min,
            integral[0].mu_max,
            integral[0].mu_steps);

    fprintf(file,
            "nu[min,max,steps]: %lf, %lf, %u\n",
            integral[0].nu_min,
            integral[0].nu_max,
            integral[0].nu_steps);

    fprintf(file, "number_cuts: %u\n", ap->number_integrals - 1);

    for (i = 1; i < ap->number_integrals; i++)
    {
        fprintf(file,
                "r_cut[min,max,steps][3]: %lf, %lf, %u\n",
                integral[i].r_min,
                integral[i].r_max,
                integral[i].r_steps);

        fprintf(file,
                "mu_cut[min,max,steps][3]: %lf, %lf, %u\n",
                integral[i].mu_min,
                integral[i].mu_max,
                integral[i].mu_steps);

        fprintf(file,
                "nu_cut[min,max,steps][3]: %lf, %lf, %u\n",
                integral[i].nu_min,
                integral[i].nu_max,
                integral[i].nu_steps);
    }
}

int write_parameters(const char* filename,
                     ASTRONOMY_PARAMETERS* ap,
                     INTEGRAL_AREA* ias,
                     BACKGROUND_PARAMETERS* bgp,
                     STREAMS* streams)
{
    FILE* f;

    f = mw_fopen(filename, "w");
    if (!f)
    {
        perror("Opening parameter file for writing");
        warn("Couldn't write output file '%s' to write astronomy parameters.\n", filename);
        return 1;
    }

    fwrite_parameters(f, ap, ias, bgp, streams);
    fclose(f);

    return 0;
}

unsigned int get_optimized_parameter_count(ASTRONOMY_PARAMETERS* ap,
                                           BACKGROUND_PARAMETERS* bgp,
                                           STREAMS* streams)
{
    unsigned int i, j, count = 0;

    for (i = 0; i < ap->number_background_parameters; i++)
    {
        if (bgp->optimize[i])
            ++count;
    }

    for (i = 0; i < streams->number_streams; i++)
    {
        if (streams->stream_weight[i].optimize)
            ++count;

        for (j = 0; j < streams->number_stream_parameters; j++)
        {
            if (streams->parameters[i].stream_optimize[j])
                ++count;
        }
    }

    return count;
}

void set_parameters(ASTRONOMY_PARAMETERS* ap,
                    BACKGROUND_PARAMETERS* bgp,
                    STREAMS* streams,
                    const real* parameters)
{
    unsigned int i, j;
    unsigned int current = 0;

    for (i = 0; i < ap->number_background_parameters; i++)
    {
        if (bgp->optimize[i])
        {
            bgp->parameters[i] = parameters[current];
            current++;
        }
    }

    for (i = 0; i < ap->number_streams; i++)
    {
        if (streams->stream_weight[i].optimize)
        {
            streams->stream_weight[i].weight = parameters[current];
            ++current;
        }

        for (j = 0; j < streams->number_stream_parameters; j++)
        {
            if (streams->parameters[i].stream_optimize[j])
            {
                streams->parameters[i].stream_parameters[j] = parameters[current];
                ++current;
            }
        }
    }
}

