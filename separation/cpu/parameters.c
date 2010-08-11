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

#include "milkyway.h"
#include "milkyway_priv.h"
#include "parameters.h"

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

void free_astronomy_parameters(ASTRONOMY_PARAMETERS* ap)
{
    free(ap->integral);
}

int read_parameters(const char* filename,
                    ASTRONOMY_PARAMETERS* ap,
                    BACKGROUND_PARAMETERS* bgp,
                    STREAMS* streams)
{
#if BOINC_APPLICATION
    char input_path[512];
    int retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));

    if (retval)
    {
        fprintf(stderr, "APP: error resolving parameters file [%s], %d\n", filename, retval);
        return retval;
    }

    FILE* data_file = boinc_fopen(input_path, "r");
#else
    FILE* data_file = fopen(filename, "r");
#endif
    if (!data_file)
    {
        fprintf(stderr, "Couldn't find input file [%s] to read astronomy parameters.\n", filename);
        return 1;
    }

    fread_parameters(data_file, ap, bgp, streams);
    if (ap->parameters_version < 0)
    {
        fprintf(stderr, "Input file [%s] did not specify parameter file version.\n", filename);
        return 1;
    }
    fclose(data_file);
    return 0;
}

int write_parameters(const char* filename,
                     ASTRONOMY_PARAMETERS* ap,
                     BACKGROUND_PARAMETERS* bgp,
                     STREAMS* streams)
{
    FILE* data_file = fopen(filename, "w");
    if (!data_file)
    {
        fprintf(stderr, "Couldn't find output file [%s] to write astronomy parameters.\n", filename);
        return 1;
    }

    fwrite_parameters(data_file, ap, bgp, streams);
    fclose(data_file);
    return 0;
}

static void calc_integral_step_sizes(INTEGRAL* i)
{
    i->r_step_size = (i->r_max - i->r_min) / (double)i->r_steps;
    i->mu_step_size = (i->mu_max - i->mu_min) / (double)i->mu_steps;
    i->nu_step_size = (i->nu_max - i->nu_min) / (double)i->nu_steps;
    i->min_calculation = 0;
    i->max_calculation = i->r_steps * i->mu_steps * i->nu_steps;
}

void fread_parameters(FILE* file,
                      ASTRONOMY_PARAMETERS* ap,
                      BACKGROUND_PARAMETERS* bgp,
                      STREAMS* streams)
{
    unsigned int i, temp;
    int retval;

    retval = fscanf(file, "parameters_version: %lf\n", &ap->parameters_version);
    if (retval < 1)
    {
        ap->parameters_version = 0.01;
//      fprintf(stderr, "Error reading astronomy parameters file. Parameters version not specified\n");
//      return;
    }

    fscanf(file, "number_parameters: %u\n", &ap->number_background_parameters);
    fscanf(file, "background_weight: %lf\n", &ap->background_weight);

    bgp->parameters = fread_double_array(file, "background_parameters", NULL);
    bgp->step       = fread_double_array(file, "background_step", NULL);
    bgp->min        = fread_double_array(file, "background_min", NULL);
    bgp->max        = fread_double_array(file, "background_max", NULL);
    bgp->optimize   = fread_int_array(file, "optimize_parameter", NULL);

    fscanf(file, "number_streams: %u, %u\n", &streams->number_streams, &streams->number_stream_parameters);

    ap->number_streams = streams->number_streams;

    streams->stream_weight = (STREAM_WEIGHT*) malloc(sizeof(STREAM_WEIGHT) * streams->number_streams);
    streams->parameters = (STREAM_PARAMETERS*) malloc(sizeof(STREAM_PARAMETERS) * streams->number_streams);

    for (i = 0; i < streams->number_streams; ++i)
    {
        fscanf(file, "stream_weight: %lf\n", &streams->stream_weight[i].weight);
        fscanf(file, "stream_weight_step: %lf\n", &streams->stream_weight[i].step);
        fscanf(file, "stream_weight_min: %lf\n", &streams->stream_weight[i].min);
        fscanf(file, "stream_weight_max: %lf\n", &streams->stream_weight[i].max);
        fscanf(file, "optimize_weight: %d\n", &streams->stream_weight[i].optimize);

        streams->parameters[i].stream_parameters = fread_double_array(file, "stream_parameters", NULL);
        streams->parameters[i].stream_step       = fread_double_array(file, "stream_step", NULL);
        streams->parameters[i].stream_min        = fread_double_array(file, "stream_min", NULL);
        streams->parameters[i].stream_max        = fread_double_array(file, "stream_max", NULL);
        streams->parameters[i].stream_optimize   = fread_int_array(file, "optimize_parameter", NULL);
    }

    fscanf(file, "convolve: %u\n", &ap->convolve);
    fscanf(file, "sgr_coordinates: %u\n", &ap->sgr_coordinates);
    if (ap->parameters_version > 0.01)
    {
        fscanf(file, "aux_bg_profile: %u\n", &ap->aux_bg_profile);
    }
    fscanf(file, "wedge: %u\n", &ap->wedge);

    ap->integral = (INTEGRAL*) malloc(sizeof(INTEGRAL));

    fscanf(file,
           "r[min,max,steps]: %lf, %lf, %u\n",
           &ap->integral[0].r_min,
           &ap->integral[0].r_max,
           &ap->integral[0].r_steps);

    fscanf(file,
           "mu[min,max,steps]: %lf, %lf, %u\n",
           &ap->integral[0].mu_min,
           &ap->integral[0].mu_max,
           &ap->integral[0].mu_steps);

    fscanf(file,
           "nu[min,max,steps]: %lf, %lf, %u\n",
           &ap->integral[0].nu_min,
           &ap->integral[0].nu_max,
           &ap->integral[0].nu_steps);

    calc_integral_step_sizes(&ap->integral[0]);

    fscanf(file, "number_cuts: %u\n", &ap->number_integrals);
    ap->number_integrals++;
    if (ap->number_integrals > 1)
    {
        ap->integral = (INTEGRAL*) realloc(ap->integral, sizeof(INTEGRAL) * ap->number_integrals);
        if (!ap->integral)
        {
            fprintf(stderr, "realloc failed\n");
            mw_finish(EXIT_FAILURE);
        }

        for (i = 1; i < ap->number_integrals; i++)
        {
            fscanf(file,
                   "r_cut[min,max,steps][%u]: %lf, %lf, %u\n",
                   &temp,
                   &ap->integral[i].r_min,
                   &ap->integral[i].r_max,
                   &ap->integral[i].r_steps);

            fscanf(file,
                   "mu_cut[min,max,steps][%u]: %lf, %lf, %u\n",
                   &temp,
                   &ap->integral[i].mu_min,
                   &ap->integral[i].mu_max,
                   &ap->integral[i].mu_steps);

            fscanf(file,
                   "nu_cut[min,max,steps][%u]: %lf, %lf, %u\n",
                   &temp,
                   &ap->integral[i].nu_min,
                   &ap->integral[i].nu_max,
                   &ap->integral[i].nu_steps);

            calc_integral_step_sizes(&ap->integral[i]);
        }
    }
}

void fwrite_parameters(FILE* file,
                       ASTRONOMY_PARAMETERS* ap,
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
            ap->integral[0].r_min,
            ap->integral[0].r_max,
            ap->integral[0].r_steps);

    fprintf(file,
            "mu[min,max,steps]: %lf, %lf, %u\n",
            ap->integral[0].mu_min,
            ap->integral[0].mu_max,
            ap->integral[0].mu_steps);

    fprintf(file,
            "nu[min,max,steps]: %lf, %lf, %u\n",
            ap->integral[0].nu_min,
            ap->integral[0].nu_max,
            ap->integral[0].nu_steps);

    fprintf(file, "number_cuts: %u\n", ap->number_integrals - 1);

    for (i = 1; i < ap->number_integrals; i++)
    {
        fprintf(file,
                "r_cut[min,max,steps][3]: %lf, %lf, %u\n",
                ap->integral[i].r_min,
                ap->integral[i].r_max,
                ap->integral[i].r_steps);

        fprintf(file,
                "mu_cut[min,max,steps][3]: %lf, %lf, %u\n",
                ap->integral[i].mu_min,
                ap->integral[i].mu_max,
                ap->integral[i].mu_steps);

        fprintf(file,
                "nu_cut[min,max,steps][3]: %lf, %lf, %u\n",
                ap->integral[i].nu_min,
                ap->integral[i].nu_max,
                ap->integral[i].nu_steps);
    }
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
                    double* parameters)
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

