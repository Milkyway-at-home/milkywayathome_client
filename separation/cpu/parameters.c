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
#include "../util/io_util.h"

void free_parameters(ASTRONOMY_PARAMETERS* ap)
{
    unsigned int i;
    free(ap->background_parameters);
    free(ap->background_step);
    free(ap->background_min);
    free(ap->background_max);
    free(ap->background_optimize);

    free(ap->stream);
    free(ap->parameters);

    for (i = 0; i < ap->number_integrals; i++)
        free(ap->integral[i]);
    free(ap->integral);
}

int read_astronomy_parameters(const char* filename, ASTRONOMY_PARAMETERS* ap)
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

    fread_astronomy_parameters(data_file, ap);
    if (ap->parameters_version < 0)
    {
        fprintf(stderr, "Input file [%s] did not specify parameter file version.\n", filename);
        return 1;
    }
    fclose(data_file);
    return 0;
}

int write_astronomy_parameters(const char* filename, ASTRONOMY_PARAMETERS* ap)
{
    FILE* data_file = fopen(filename, "w");
    if (!data_file)
    {
        fprintf(stderr, "Couldn't find output file [%s] to write astronomy parameters.\n", filename);
        return 1;
    }

    fwrite_astronomy_parameters(data_file, ap);
    fclose(data_file);
    return 0;
}

void fread_astronomy_parameters(FILE* file, ASTRONOMY_PARAMETERS* ap)
{
    unsigned int i;
    int retval;

    retval = fscanf(file, "parameters_version: %lf\n", &ap->parameters_version);
    if (retval < 1)
    {
        ap->parameters_version = 0.01;
//      fprintf(stderr, "Error reading astronomy parameters file. Parameters version not specified\n");
//      return;
    }

    fscanf(file, "number_parameters: %d\n", &ap->number_background_parameters);
    fscanf(file, "background_weight: %lf\n", &ap->background_weight);

    ap->background_parameters = fread_double_array(file, "background_parameters");
    ap->background_step       = fread_double_array(file, "background_step");
    ap->background_min        = fread_double_array(file, "background_min");
    ap->background_max        = fread_double_array(file, "background_max");
    ap->background_optimize   = fread_int_array(file, "optimize_parameter");

    fscanf(file, "number_streams: %d, %d\n", &ap->number_streams, &ap->number_stream_parameters);

    ap->stream = (STREAM*) malloc(sizeof(STREAM) * ap->number_streams);
    ap->parameters = (STREAM_PARAMETERS*) malloc(sizeof(STREAM_PARAMETERS) * ap->number_streams);

    for (i = 0; i < ap->number_streams; ++i)
    {
        fscanf(file, "stream_weight: %lf\n", &STREAM_N(ap,i).weights);
        fscanf(file, "stream_weight_step: %lf\n", &STREAM_N(ap, i).weight_step);
        fscanf(file, "stream_weight_min: %lf\n", &STREAM_N(ap, i).weight_min);
        fscanf(file, "stream_weight_max: %lf\n", &STREAM_N(ap, i).weight_max);
        fscanf(file, "optimize_weight: %d\n", &STREAM_N(ap, i).weight_optimize);

        STREAM_PARAM_N(ap, i).stream_parameters = fread_double_array(file, "stream_parameters");
        STREAM_PARAM_N(ap, i).stream_step       = fread_double_array(file, "stream_step");
        STREAM_PARAM_N(ap, i).stream_min        = fread_double_array(file, "stream_min");
        STREAM_PARAM_N(ap, i).stream_max        = fread_double_array(file, "stream_max");
        STREAM_PARAM_N(ap, i).stream_optimize   = fread_int_array(file, "optimize_parameter");
    }

    fscanf(file, "convolve: %d\n", &ap->convolve);
    fscanf(file, "sgr_coordinates: %d\n", &ap->sgr_coordinates);
    if (ap->parameters_version > 0.01)
    {
        fscanf(file, "aux_bg_profile: %d\n", &ap->aux_bg_profile);
    }
    fscanf(file, "wedge: %d\n", &ap->wedge);

    ap->integral = (INTEGRAL**)malloc(sizeof(INTEGRAL*));
    ap->integral[0] = (INTEGRAL*)malloc(sizeof(INTEGRAL));

    fscanf(file,
           "r[min,max,steps]: %lf, %lf, %d\n",
           &ap->integral[0]->r_min,
           &ap->integral[0]->r_max,
           &ap->integral[0]->r_steps);

    fscanf(file,
           "mu[min,max,steps]: %lf, %lf, %d\n",
           &ap->integral[0]->mu_min,
           &ap->integral[0]->mu_max,
           &ap->integral[0]->mu_steps);

    fscanf(file,
           "nu[min,max,steps]: %lf, %lf, %d\n",
           &ap->integral[0]->nu_min,
           &ap->integral[0]->nu_max,
           &ap->integral[0]->nu_steps);

    ap->integral[0]->r_step_size = (ap->integral[0]->r_max - ap->integral[0]->r_min) / (double)ap->integral[0]->r_steps;
    ap->integral[0]->mu_step_size = (ap->integral[0]->mu_max - ap->integral[0]->mu_min) / (double)ap->integral[0]->mu_steps;
    ap->integral[0]->nu_step_size = (ap->integral[0]->nu_max - ap->integral[0]->nu_min) / (double)ap->integral[0]->nu_steps;
    ap->integral[0]->min_calculation = 0;
    ap->integral[0]->max_calculation = ap->integral[0]->r_steps * ap->integral[0]->mu_steps * ap->integral[0]->nu_steps;

    fscanf(file, "number_cuts: %d\n", &ap->number_integrals);
    ap->number_integrals++;
    if (ap->number_integrals > 1)
    {
        ap->integral = (INTEGRAL**)realloc(ap->integral, sizeof(INTEGRAL*) * ap->number_integrals);
        for (i = 1; i < ap->number_integrals; i++)
        {
            int temp;
            ap->integral[i] = (INTEGRAL*)malloc(sizeof(INTEGRAL));
            fscanf(file,
                   "r_cut[min,max,steps][%d]: %lf, %lf, %d\n",
                   &temp,
                   &ap->integral[i]->r_min,
                   &ap->integral[i]->r_max,
                   &ap->integral[i]->r_steps);

            fscanf(file,
                   "mu_cut[min,max,steps][%d]: %lf, %lf, %d\n",
                   &temp,
                   &ap->integral[i]->mu_min,
                   &ap->integral[i]->mu_max,
                   &ap->integral[i]->mu_steps);

            fscanf(file,
                   "nu_cut[min,max,steps][%d]: %lf, %lf, %d\n",
                   &temp,
                   &ap->integral[i]->nu_min,
                   &ap->integral[i]->nu_max,
                   &ap->integral[i]->nu_steps);

            ap->integral[i]->r_step_size = (ap->integral[i]->r_max - ap->integral[i]->r_min) / (double)ap->integral[i]->r_steps;
            ap->integral[i]->mu_step_size = (ap->integral[i]->mu_max - ap->integral[i]->mu_min) / (double)ap->integral[i]->mu_steps;
            ap->integral[i]->nu_step_size = (ap->integral[i]->nu_max - ap->integral[i]->nu_min) / (double)ap->integral[i]->nu_steps;
            ap->integral[i]->min_calculation = 0;
            ap->integral[i]->max_calculation = ap->integral[i]->r_steps * ap->integral[i]->mu_steps * ap->integral[i]->nu_steps;
        }
    }
}

void fwrite_astronomy_parameters(FILE* file, ASTRONOMY_PARAMETERS* ap)
{
    unsigned int i;

    fprintf(file, "parameters_version: %lf\n", ap->parameters_version);

    fprintf(file, "number_parameters: %d\n", ap->number_background_parameters);
    fprintf(file, "background_weight: %lf\n", ap->background_weight);
    fwrite_double_array(file, "background_parameters", ap->background_parameters, ap->number_background_parameters);
    fwrite_double_array(file, "background_step", ap->background_step, ap->number_background_parameters);
    fwrite_double_array(file, "background_min", ap->background_min, ap->number_background_parameters);
    fwrite_double_array(file, "background_max", ap->background_max, ap->number_background_parameters);
    fwrite_int_array(file, "optimize_parameter", ap->background_optimize, ap->number_background_parameters);

    fprintf(file, "number_streams: %d, %d\n", ap->number_streams, ap->number_stream_parameters);
    for (i = 0; i < ap->number_streams; i++)
    {
        fprintf(file, "stream_weight: %lf\n", STREAM_N(ap, i).weights);
        fprintf(file, "stream_weight_step: %lf\n", STREAM_N(ap, i).weight_step);
        fprintf(file, "stream_weight_min: %lf\n", STREAM_N(ap, i).weight_min);
        fprintf(file, "stream_weight_max: %lf\n", STREAM_N(ap, i).weight_max);
        fprintf(file, "optimize_weight: %d\n", STREAM_N(ap, i).weight_optimize);

        fwrite_double_array(file, "stream_parameters",
                            STREAM_PARAM_N(ap, i).stream_parameters, ap->number_stream_parameters);
        fwrite_double_array(file, "stream_step",
                            STREAM_PARAM_N(ap, i).stream_step, ap->number_stream_parameters);
        fwrite_double_array(file, "stream_min",
                            STREAM_PARAM_N(ap, i).stream_min, ap->number_stream_parameters);
        fwrite_double_array(file, "stream_max",
                            STREAM_PARAM_N(ap, i).stream_max, ap->number_stream_parameters);
        fwrite_int_array(file, "optimize_parameter",
                         STREAM_PARAM_N(ap, i).stream_optimize, ap->number_stream_parameters);
    }

    fprintf(file, "convolve: %d\n", ap->convolve);
    fprintf(file, "sgr_coordinates: %d\n", ap->sgr_coordinates);
    fprintf(file, "aux_bg_profile: %d\n", ap->aux_bg_profile); //vickej2_bg
    fprintf(file, "wedge: %d\n", ap->wedge);

    fprintf(file,
            "r[min,max,steps]: %lf, %lf, %d\n",
            ap->integral[0]->r_min,
            ap->integral[0]->r_max,
            ap->integral[0]->r_steps);

    fprintf(file,
            "mu[min,max,steps]: %lf, %lf, %d\n",
            ap->integral[0]->mu_min,
            ap->integral[0]->mu_max,
            ap->integral[0]->mu_steps);

    fprintf(file,
            "nu[min,max,steps]: %lf, %lf, %d\n",
            ap->integral[0]->nu_min,
            ap->integral[0]->nu_max,
            ap->integral[0]->nu_steps);

    fprintf(file, "number_cuts: %d\n", ap->number_integrals - 1);

    for (i = 1; i < ap->number_integrals; i++)
    {
        fprintf(file,
                "r_cut[min,max,steps][3]: %lf, %lf, %d\n",
                ap->integral[i]->r_min,
                ap->integral[i]->r_max,
                ap->integral[i]->r_steps);

        fprintf(file,
                "mu_cut[min,max,steps][3]: %lf, %lf, %d\n",
                ap->integral[i]->mu_min,
                ap->integral[i]->mu_max,
                ap->integral[i]->mu_steps);

        fprintf(file,
                "nu_cut[min,max,steps][3]: %lf, %lf, %d\n",
                ap->integral[i]->nu_min,
                ap->integral[i]->nu_max,
                ap->integral[i]->nu_steps);
    }
}

unsigned int get_optimized_parameter_count(ASTRONOMY_PARAMETERS* ap)
{
    unsigned int i, j, count = 0;

    for (i = 0; i < ap->number_background_parameters; i++)
    {
        if (ap->background_optimize[i])
            ++count;
    }

    for (i = 0; i < ap->number_streams; i++)
    {
        if (STREAM_N(ap, i).weight_optimize)
            ++count;

        for (j = 0; j < ap->number_stream_parameters; j++)
        {
            if (STREAM_PARAM_N(ap, i).stream_optimize[j])
                ++count;
        }
    }

    return count;
}

void set_astronomy_parameters(ASTRONOMY_PARAMETERS* ap, double* parameters)
{
    unsigned int i, j;
    unsigned int current = 0;

    for (i = 0; i < ap->number_background_parameters; i++)
    {
        if (ap->background_optimize[i])
        {
            ap->background_parameters[i] = parameters[current];
            current++;
        }
    }

    for (i = 0; i < ap->number_streams; i++)
    {
        if (STREAM_N(ap, i).weight_optimize)
        {
            STREAM_N(ap, i).weights = parameters[current];
            ++current;
        }

        for (j = 0; j < ap->number_stream_parameters; j++)
        {
            if (STREAM_PARAM_N(ap, i).stream_optimize[j])
            {
                STREAM_PARAM_N(ap, i).stream_parameters[j] = parameters[current];
                ++current;
            }
        }
    }
}

