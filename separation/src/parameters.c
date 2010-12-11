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

#include "separation_types.h"
#include "parameters.h"
#include "io_util.h"
#include "milkyway_util.h"

void freeBackgroundParameters(BackgroundParameters* bgp)
{
    free(bgp->parameters);
    free(bgp->step);
    free(bgp->min);
    free(bgp->max);
    free(bgp->optimize);
}

void freeStreamParameters(StreamParameters* p)
{
    free(p->stream_parameters);
    free(p->stream_step);
    free(p->stream_min);
    free(p->stream_max);
    free(p->stream_optimize);
}

void freeStreams(Streams* streams)
{
    unsigned int i;

    for (i = 0; i < streams->number_streams; ++i)
        freeStreamParameters(&streams->parameters[i]);
    free(streams->parameters);
    free(streams->stream_weight);
}

static void calcIntegralStepSizes(IntegralArea* i)
{
    i->r_step_size = (i->r_max - i->r_min) / (real)i->r_steps;
    i->mu_step_size = (i->mu_max - i->mu_min) / (real)i->mu_steps;
    i->nu_step_size = (i->nu_max - i->nu_min) / (real)i->nu_steps;
}

static IntegralArea* freadParameters(FILE* file,
                                     AstronomyParameters* ap,
                                     BackgroundParameters* bgp,
                                     Streams* streams)
{
    unsigned int i, temp;
    double tmp1, tmp2;
    IntegralArea integralTmp;
    IntegralArea* integrals;
    const IntegralArea* ia;
    unsigned int total_calc_probs;
    unsigned int integralNumTmp;

    ap->parameters_version = (fscanf(file, "parameters_version: %lf\n", &tmp1) < 1) ? 0.01 : (real) tmp1;

    if (fscanf(file, "number_parameters: %u\n", &ap->number_background_parameters) < 1)
        warn("Error reading number_parameters\n");

    if (fscanf(file, "background_weight: %lf\n", &tmp1) < 1)
        warn("Error reading background_weight\n");

    ap->background_weight = (real) tmp1;

    bgp->parameters = fread_double_array(file, "background_parameters", NULL);
    bgp->step       = fread_double_array(file, "background_step", NULL);
    bgp->min        = fread_double_array(file, "background_min", NULL);
    bgp->max        = fread_double_array(file, "background_max", NULL);
    bgp->optimize   = fread_int_array(file, "optimize_parameter", NULL);

    if (fscanf(file, "number_streams: %u, %u\n", &streams->number_streams, &streams->number_stream_parameters) < 2)
        warn("Error reading number_streams\n");

    ap->number_streams = streams->number_streams;

    streams->stream_weight = (StreamWeight*) mwMalloc(sizeof(StreamWeight) * streams->number_streams);
    streams->parameters = (StreamParameters*) mwMalloc(sizeof(StreamParameters) * streams->number_streams);

    for (i = 0; i < streams->number_streams; ++i)
    {
        if (fscanf(file, "stream_weight: %lf\n", &tmp1) < 1)
            warn("Error reading stream_weight for stream %u\n", i);
        streams->stream_weight[i].weight = (real) tmp1;

        if (fscanf(file, "stream_weight_step: %lf\n", &tmp1) < 1)
            warn("Error reading stream_weight_step for stream %u\n", i);
        streams->stream_weight[i].step = (real) tmp1;

        if (fscanf(file, "stream_weight_min: %lf\n", &tmp1) < 1)
            warn("Error reading stream_weight_min for stream %u\n", i);
        streams->stream_weight[i].min = (real) tmp1;

        if (fscanf(file, "stream_weight_max: %lf\n", &tmp1) < 1)
            warn("Error reading stream_weight_max for stream %u\n", i);
        streams->stream_weight[i].max = (real) tmp1;

        if (fscanf(file, "optimize_weight: %d\n", &streams->stream_weight[i].optimize) < 1)
            warn("Error reading optimize_weight for stream %u\n", i);

        streams->parameters[i].stream_parameters = fread_double_array(file, "stream_parameters", NULL);
        streams->parameters[i].stream_step       = fread_double_array(file, "stream_step", NULL);
        streams->parameters[i].stream_min        = fread_double_array(file, "stream_min", NULL);
        streams->parameters[i].stream_max        = fread_double_array(file, "stream_max", NULL);
        streams->parameters[i].stream_optimize   = fread_int_array(file, "optimize_parameter", NULL);
    }

    if (fscanf(file, "convolve: %u\n", &ap->convolve) < 1)
        warn("Error reading convolve\n");

    if (fscanf(file, "sgr_coordinates: %d\n", &ap->sgr_coordinates) < 1)
        warn("Error reading sgr_coordinates\n");
    if (ap->parameters_version > 0.01)
    {
        if (fscanf(file, "aux_bg_profile: %d\n", &ap->aux_bg_profile) < 1)
            warn("Error reading aux_bg_profile\n");
    }

    if (fscanf(file, "wedge: %d\n", &ap->wedge) < 1)
        warn("Error reading wedge\n");

    if (fscanf(file,
               "r[min,max,steps]: %lf, %lf, %u\n",
               &tmp1, &tmp2, &integralTmp.r_steps) < 3)
        warn("Error reading r\n");

    integralTmp.r_min = (real) tmp1;
    integralTmp.r_max = (real) tmp2;

    if (fscanf(file,
               "mu[min,max,steps]: %lf, %lf, %u\n",
               &tmp1, &tmp2, &integralTmp.mu_steps) < 3)
        warn("Error reading mu\n");

    integralTmp.mu_min = (real) tmp1;
    integralTmp.mu_max = (real) tmp2;

    if (fscanf(file,
               "nu[min,max,steps]: %lf, %lf, %u\n",
               &tmp1, &tmp2, &integralTmp.nu_steps) < 3)
        warn("Error reading nu\n");

    integralTmp.nu_min = (real) tmp1;
    integralTmp.nu_max = (real) tmp2;

    calcIntegralStepSizes(&integralTmp);

    ap->number_integrals = 1;
    if (fscanf(file, "number_cuts: %u\n", &integralNumTmp) < 1)
        warn("Error reading number_cuts\n");
    ap->number_integrals += integralNumTmp;

    integrals = (IntegralArea*) mwMallocA(ap->number_integrals * sizeof(IntegralArea));

    integrals[0] = integralTmp;

    if (ap->number_integrals > 1)
    {
        for (i = 1; i < ap->number_integrals; i++)
        {
            if (fscanf(file,
                       "r_cut[min,max,steps][%u]: %lf, %lf, %u\n",
                       &temp, &tmp1, &tmp2, &integrals[i].r_steps) < 3)
            {
                warn("Error reading r for integral %u\n", i);
            }

            integrals[i].r_min = (real) tmp1;
            integrals[i].r_max = (real) tmp2;

            if (fscanf(file,
                       "mu_cut[min,max,steps][%u]: %lf, %lf, %u\n",
                       &temp, &tmp1, &tmp2, &integrals[i].mu_steps) < 3)
            {
                warn("Error reading mu for integral %u\n", i);
            }


            integrals[i].mu_min = (real) tmp1;
            integrals[i].mu_max = (real) tmp2;

            if (fscanf(file,
                       "nu_cut[min,max,steps][%u]: %lf, %lf, %u\n",
                       &temp, &tmp1, &tmp2, &integrals[i].nu_steps) < 3)
            {
                warn("Error reading nu for integral %u\n", i);
            }

            integrals[i].nu_min = (real) tmp1;
            integrals[i].nu_max = (real) tmp2;

            calcIntegralStepSizes(&integrals[i]);
        }
    }

    /* Calculate total probability calculations for checkpointing */
    total_calc_probs = 0;
    for (i = 0; i < ap->number_integrals; ++i)
    {
        ia = &integrals[i];
        total_calc_probs += ia->mu_steps * ia->nu_steps * ia->r_steps;
    }

    ap->total_calc_probs = (real) total_calc_probs;
    return integrals;
}

IntegralArea* readParameters(const char* filename,
                             AstronomyParameters* ap,
                             BackgroundParameters* bgp,
                             Streams* streams)
{
    FILE* f;
    IntegralArea* integral;

    f = mwOpenResolved(filename, "r");
    if (!f)
    {
        perror("Opening astronomy parameters file");
        return NULL;
    }

    integral = freadParameters(f, ap, bgp, streams);
    if (!integral)
        warn("Error reading parameters file\n");

    fclose(f);
    return integral;
}

static void fwriteParameters(FILE* file,
                             AstronomyParameters* ap,
                             IntegralArea* integral,
                             BackgroundParameters* bgp,
                             Streams* streams)
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
    fprintf(file, "aux_bg_profile: %d\n", ap->aux_bg_profile);
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

int writeParameters(const char* filename,
                    AstronomyParameters* ap,
                    IntegralArea* ias,
                    BackgroundParameters* bgp,
                    Streams* streams)
{
    FILE* f;

    f = mw_fopen(filename, "w");
    if (!f)
    {
        perror("Opening parameter file for writing");
        warn("Couldn't write output file '%s' to write astronomy parameters.\n", filename);
        return 1;
    }

    fwriteParameters(f, ap, ias, bgp, streams);
    fclose(f);

    return 0;
}

unsigned int getOptimizedParameterCount(AstronomyParameters* ap,
                                        BackgroundParameters* bgp,
                                        Streams* streams)
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

void setParameters(AstronomyParameters* ap,
                   BackgroundParameters* bgp,
                   Streams* streams,
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

