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

void freeStreams(Streams* streams)
{
    free(streams->parameters);
}

void calcIntegralStepSizes(IntegralArea* i)
{
    i->r_step_size = (i->r_max - i->r_min) / (real) i->r_steps;
    i->mu_step_size = (i->mu_max - i->mu_min) / (real) i->mu_steps;
    i->nu_step_size = (i->nu_max - i->nu_min) / (real) i->nu_steps;
}

static int checkIntegralAreasOK(const IntegralArea* ias, unsigned int n)
{
    unsigned int i;
    const IntegralArea* ia;

    for (i = 0; i < n; ++i)
    {
        ia = &ias[i];

        if (ia->nu_steps == 0 || ia->r_steps == 0 || ia->mu_steps == 0)
        {
            warn("Integral area dimensions cannot be 0\n");
            return 1;
        }

        if (!mwEven(ia->nu_steps) || !mwEven(ia->r_steps) || !mwEven(ia->mu_steps))
        {
            warn("Integral area dimensions must be even: cut %u: "
                 "{ nu_steps = %u, mu_steps = %u, r_steps = %u }\n",
                 i, ia->nu_steps, ia->mu_steps, ia->r_steps);
            return 1;
        }
    }

    return 0;
}

static IntegralArea* freadParameters(FILE* file,
                                     AstronomyParameters* ap,
                                     BackgroundParameters* bgp,
                                     Streams* streams)
{
    unsigned int i, temp;
    double tmp1, tmp2;
    double parametersVersion;
    IntegralArea integralTmp;
    IntegralArea* integrals;
    const IntegralArea* ia;
    unsigned int total_calc_probs;
    unsigned int integralNumTmp;
    int iTmp;
    int sgr_coordinates = 0;
    real* tmpArr = NULL;

    parametersVersion = (fscanf(file, "parameters_version: %lf\n", &tmp1) < 1) ? 0.01 : (real) tmp1;

    if (fscanf(file, "number_parameters: %u\n", &temp) < 1)
        fail("Error reading number_parameters\n");

    if (fscanf(file, "background_weight: %lf\n", &tmp1) < 1)
        fail("Error reading background_weight\n");

    bgp->epsilon = (real) tmp1;

    tmpArr = fread_double_array(file, "background_parameters", NULL);
    if (!tmpArr)
        return NULL;

    bgp->alpha = tmpArr[0];
    bgp->q     = tmpArr[1];
    bgp->r0    = tmpArr[2];
    bgp->delta = tmpArr[3];
    free(tmpArr);

    free(fread_double_array(file, "background_step", NULL));
    free(fread_double_array(file, "background_min", NULL));
    free(fread_double_array(file, "background_max", NULL));
    free(fread_int_array(file, "optimize_parameter", NULL));

    if (fscanf(file, "number_streams: %u, %u\n", &streams->number_streams, &temp) < 2)
        fail("Error reading number_streams\n");

    ap->number_streams = streams->number_streams;

    streams->parameters = (StreamParameters*) mwCalloc(streams->number_streams, sizeof(StreamParameters));

    for (i = 0; i < streams->number_streams; ++i)
    {
        if (fscanf(file, "stream_weight: %lf\n", &tmp1) < 1)
            fail("Error reading stream_weight for stream %u\n", i);
        streams->parameters[i].epsilon = (real) tmp1;

        if (fscanf(file, "stream_weight_step: %lf\n", &tmp1) < 1)
            fail("Error reading stream_weight_step for stream %u\n", i);

        if (fscanf(file, "stream_weight_min: %lf\n", &tmp1) < 1)
            fail("Error reading stream_weight_min for stream %u\n", i);

        if (fscanf(file, "stream_weight_max: %lf\n", &tmp1) < 1)
            fail("Error reading stream_weight_max for stream %u\n", i);

        if (fscanf(file, "optimize_weight: %d\n", &iTmp) < 1)
            fail("Error reading optimize_weight for stream %u\n", i);

        tmpArr = fread_double_array(file, "stream_parameters", NULL);
        if (!tmpArr)
            return NULL;
        streams->parameters[i].mu    = tmpArr[0];
        streams->parameters[i].r     = tmpArr[1];
        streams->parameters[i].theta = tmpArr[2];
        streams->parameters[i].phi   = tmpArr[3];
        streams->parameters[i].sigma = tmpArr[4];
        free(tmpArr);

        free(fread_double_array(file, "stream_step", NULL));
        free(fread_double_array(file, "stream_min", NULL));
        free(fread_double_array(file, "stream_max", NULL));
        free(fread_int_array(file, "optimize_parameter", NULL));
    }

    if (fscanf(file, "convolve: %u\n", &ap->convolve) < 1)
        fail("Error reading convolve\n");

    if (fscanf(file, "sgr_coordinates: %d\n", &sgr_coordinates) < 1)
        fail("Error reading sgr_coordinates\n");

    if (sgr_coordinates)
    {
        fail("sgr_coordinates unimplemented\n");
    }

    if (parametersVersion > 0.01)
    {
        if (fscanf(file, "aux_bg_profile: %d\n", &ap->aux_bg_profile) < 1)
            fail("Error reading aux_bg_profile\n");
    }

    if (fscanf(file, "wedge: %d\n", &ap->wedge) < 1)
        fail("Error reading wedge\n");

    if (fscanf(file,
               "r[min,max,steps]: %lf, %lf, %u\n",
               &tmp1, &tmp2, &integralTmp.r_steps) < 3)
        fail("Error reading r\n");

    integralTmp.r_min = (real) tmp1;
    integralTmp.r_max = (real) tmp2;

    if (fscanf(file,
               "mu[min,max,steps]: %lf, %lf, %u\n",
               &tmp1, &tmp2, &integralTmp.mu_steps) < 3)
        fail("Error reading mu\n");

    integralTmp.mu_min = (real) tmp1;
    integralTmp.mu_max = (real) tmp2;

    if (fscanf(file,
               "nu[min,max,steps]: %lf, %lf, %u\n",
               &tmp1, &tmp2, &integralTmp.nu_steps) < 3)
        fail("Error reading nu\n");

    integralTmp.nu_min = (real) tmp1;
    integralTmp.nu_max = (real) tmp2;

    calcIntegralStepSizes(&integralTmp);

    ap->number_integrals = 1;
    if (fscanf(file, "number_cuts: %u\n", &integralNumTmp) < 1)
        fail("Error reading number_cuts\n");
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
                fail("Error reading r for integral %u\n", i);
            }

            integrals[i].r_min = (real) tmp1;
            integrals[i].r_max = (real) tmp2;

            if (fscanf(file,
                       "mu_cut[min,max,steps][%u]: %lf, %lf, %u\n",
                       &temp, &tmp1, &tmp2, &integrals[i].mu_steps) < 3)
            {
                fail("Error reading mu for integral %u\n", i);
            }


            integrals[i].mu_min = (real) tmp1;
            integrals[i].mu_max = (real) tmp2;

            if (fscanf(file,
                       "nu_cut[min,max,steps][%u]: %lf, %lf, %u\n",
                       &temp, &tmp1, &tmp2, &integrals[i].nu_steps) < 3)
            {
                fail("Error reading nu for integral %u\n", i);
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

    if (checkIntegralAreasOK(integrals, ap->number_integrals))
    {
        free(integrals);
        return NULL;
    }

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


int setParameters(AstronomyParameters* ap,
                  BackgroundParameters* bgp,
                  Streams* streams,
                  const real* parameters,
                  unsigned int numberParameters)
{
    unsigned int i, idx;
    const unsigned int nBGParams = 2;
    const unsigned int nStreamParams = 6;
    unsigned int nStream = (numberParameters - nBGParams) / nStreamParams;

    if (nStream != ap->number_streams)
    {
        warn("Number of streams does not match\n");
        return 1;
    }

    if (!mwDivisible(numberParameters - nBGParams, nStreamParams))
    {
        warn("Number of parameters doesn't make sense\n");
        return 1;
    }

    bgp->q = parameters[0];
    bgp->r0 = parameters[1];

    for (i = 0; i < nStream; ++i)
    {
        idx = nBGParams + i * nStreamParams;

        streams->parameters[i].epsilon = parameters[idx + 0];
        streams->parameters[i].mu      = parameters[idx + 1];
        streams->parameters[i].r       = parameters[idx + 2];
        streams->parameters[i].theta   = parameters[idx + 3];
        streams->parameters[i].phi     = parameters[idx + 4];
        streams->parameters[i].sigma   = parameters[idx + 5];
    }

    return 0;
}

