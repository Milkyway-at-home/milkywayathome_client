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

#include "io_util.h"
#include "milkyway_util.h"
#include "separation_types.h"
#include "separation_utils.h"

/****
    *   Functions for printing parameters to files/output
*****/
void fwrite_double_array(FILE *file, const char *array_name, real* array_t, size_t size)
{
    size_t i;

    /* Should use %zu for size_t but Windows is broken */
    fprintf(file, "%s[%lu]: ", array_name, (unsigned long) size);
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
    fprintf(file, "%s[%lu]: ", array_name, (unsigned long) size);
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

    arr = (real*) mwMalloc(sizeof(real) * size);

    for (i = 0; i < size; i++)
    {
        rc = fscanf(file, READ_DOUBLE_ARRAY_READ_STR, &arr[i]);
        if (rc != 1)
        {
            warn("Error reading into %s\n", array_name);
            free(arr);
            return NULL;
        }

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
    int* arr;

    fscanf(file, array_name);
    fscanf(file, "[%u]: ", &size);

    arr = (int*) mwMalloc(sizeof(int) * size);

    for (i = 0; i < size; i++)
    {
        if (fscanf(file, "%d", &arr[i]) != 1)
        {
            warn("Error reading into %s\n", array_name);
            free(arr);
            return NULL;
        }

        if (i < size - 1)
            fscanf(file, ", ");
    }
    fscanf(file, "\n");

    if (sizeOut)
        *sizeOut = size;

    return arr;
}

void printIntegralArea(const IntegralArea* ia)
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
           "  mu_steps     = %u\n"
           "}\n",
           ia->r_min, ia->r_max, ia->r_step_size,
           ia->nu_min, ia->nu_max, ia->nu_step_size,
           ia->mu_min, ia->mu_max, ia->mu_step_size,
           ia->r_steps, ia->nu_steps, ia->mu_steps);
}

void printNuConstants(const NuConstants* c, unsigned int n)
{
    unsigned int i;
    printf("Nu constants:\n");
    for (i = 0; i < n; ++i)
        printf("[%u] { nu = %g, id = %g } \n", i, c[i].nu, c[i].id);
}

void printStreamGauss(const StreamGauss* c, unsigned int n)
{
    unsigned int i;
    printf("Stream gauss:\n");
    for (i = 0; i < n; ++i)
        printf("[%u] { dx = %g, qgaus_W = %g } \n", i, c->dx[i], c->qgaus_W[i]);
}

void printStreamConstants(const StreamConstants* c, unsigned int n)
{
    unsigned int i;
    printf("Stream constants:\n");
    for (i = 0; i < n; ++i)
    {
        printf("[%u] { a = { %g, %g, %g, %g }, c = { %g, %g, %g, %g }, sigma_sq2_inv = %g, large_sigma = %d } \n",
               i, X(c[i].a), Y(c[i].a), Z(c[i].a), W(c[i].a),
               X(c[i].c), Y(c[i].c), Z(c[i].c), W(c[i].c),
               c[i].sigma_sq2_inv, c[i].large_sigma);
    }
}

void printAstronomyParameters(const AstronomyParameters* ap)
{
    printf("astronomy-parameters {\n"
           "  m_sun_r0              = %f\n"
           "  q_inv                 = %f\n"
           "  q_inv_sqr             = %f\n"
           "  r0                    = %f\n"
           "  convolve              = %u\n"
           "  number_streams        = %u\n"
           "  fast_h_prob           = %d\n"
           "  aux_bg_profile        = %d\n"
           "  alpha                 = %f\n"
           "  delta                 = %f\n"
           "  alpha_delta           = %f\n"
           "  bg_a                  = %f\n"
           "  bg_b                  = %f\n"
           "  bg_c                  = %f\n"
           "  wedge                 = %d\n"
           "  sun_r0                = %f\n"
           "  q                     = %f\n"
           "  coeff                 = %f\n"
           "  total_calc_probs      = %f\n"
           "  number_integrals      = %u\n"
           "  exp_background_weight = %f\n",
           ap->m_sun_r0,
           ap->q_inv,
           ap->q_inv_sqr,
           ap->r0,
           ap->convolve,
           ap->number_streams,
           ap->fast_h_prob,
           ap->aux_bg_profile,
           ap->alpha,
           ap->delta,
           ap->alpha_delta3,
           ap->bg_a, ap->bg_b, ap->bg_c,
           ap->wedge,
           ap->sun_r0,
           ap->q,
           ap->coeff,
           ap->total_calc_probs,
           ap->number_integrals,
           ap->exp_background_weight);
}

void printSeparationResults(const SeparationResults* results, unsigned int numberStreams)
{
    unsigned int i;

    mw_begin_critical_section();

    /* Print integrals */
    warn("<background_integral> %.15f </background_integral>\n", results->backgroundIntegral);
    warn("<stream_integral> ");
    for (i = 0; i < numberStreams; ++i)
        warn(" %.15f ", results->streamIntegrals[i]);
    warn("</stream_integral>\n");

    /* Print individual likelihoods */
    warn("<background_likelihood> %.15f </background_likelihood>\n", results->backgroundLikelihood);
    warn("<stream_only_likelihood> ");
    for (i = 0; i < numberStreams; ++i)
        warn(" %.15f ", results->streamLikelihoods[i]);
    warn("</stream_only_likelihood>\n");

    /* Print overall likelihood */
    warn("<search_likelihood> %.15f </search_likelihood>\n", results->likelihood);

    mw_end_critical_section();
}

/* FIXME: Kill this with fire when we switch to JSON everything for separation */
static SeparationResults* freadReferenceResults(FILE* f, unsigned int nStream)
{
    SeparationResults* r;
    char buf[256] = "";
    unsigned int i;

    r = newSeparationResults(nStream);

    /* Can put comments at start of this junk */
    while (buf[0] == '#' || buf[0] == ' ' || !strncmp(buf, "", sizeof(buf)))
        fgets(buf, sizeof(buf), f);

    if (sscanf(buf, "likelihood: %lf\n", &r->likelihood) != 1)
        goto fail_read;

    if (fscanf(f, "background_integral: %lf\n", &r->backgroundIntegral) != 1)
        goto fail_read;

    if (fscanf(f, "background_likelihood: %lf\n", &r->backgroundLikelihood) != 1)
        goto fail_read;

    for (i = 0; i < nStream; ++i)
    {
        if (fscanf(f, "stream_integral: %lf\n", &r->streamIntegrals[i]) != 1)
            goto fail_read;

        if (fscanf(f, "stream_likelihood: %lf\n", &r->streamLikelihoods[i]) != 1)
            goto fail_read;
    }

    return r;

fail_read:
    freeSeparationResults(r);
    return NULL;
}

SeparationResults* readReferenceResults(const char* refFile, unsigned int nStream)
{
    FILE* f;
    SeparationResults* refResults;

    f = fopen(refFile, "r");
    if (!f)
    {
        perror("Opening reference results");
        return NULL;
    }

    refResults = freadReferenceResults(f, nStream);

    fclose(f);

    return refResults;
}


