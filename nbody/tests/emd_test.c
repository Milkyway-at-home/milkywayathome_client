/*
 * Copyright (c) 2011 Matthew Arsenault
 * Copyright (c) 2011 Rensselaer Polytechnic Institute
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "milkyway_util.h"
#include "nbody_emd.h"
#include "dSFMT.h"
#include <time.h>

static dsfmt_t _prng;

#define ZERO_THRESHOLD 1.0e-6

/* Function which assigns a sample distribution to arr1 and arr2 to
 * match of size n. Returns expected EMD for the distribution */
typedef float (*EMDTestDistribFunc)(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2, unsigned int n);

static void zeroWeights(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2, unsigned int n)
{
    unsigned int i;

    for (i = 0; i < n; ++i)
    {
        arr1[i].weight = arr2[i].weight = 0.0f;
    }
}

/* Returns expected value of dist */
static float oppositeSides(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2, unsigned int n)
{
    arr1[0].weight = 1.0f;
    arr2[n - 1].weight = 1.0f;

    return fabsf(arr1[0].pos - arr2[n-1].pos);
}

static float allInDifferentBins(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2, unsigned int n)
{
    unsigned int i, j;

    i = (unsigned int) mwXrandom(&_prng, 0, n);
    j = (unsigned int) mwXrandom(&_prng, 0, n);
    assert(i < n && j < n);

    arr1[i].weight = 1.0f;
    arr2[j].weight = 1.0f;

    return fabsf(arr1[i].pos - arr2[j].pos);
}

static float allInSameBin(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2, unsigned int n)
{
    unsigned int i;

    /* Pick one bin at random */
    i = (unsigned int) mwXrandom(&_prng, 0, n);
    assert(i < n);

    arr1[i].weight = arr2[i].weight = 1.0f;

    return 0.0f;
}

static float randomSelf(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2, unsigned int n)
{
    unsigned int i;
    float total = 0.0f;

    for (i = 0; i < n; ++i)
    {
        arr1[i].weight = (float) dsfmt_genrand_open_open(&_prng);
        total += arr1[i].weight;
    }

    /* Normalize it and copy */
    for (i = 0; i < n; ++i)
    {
        arr1[i].weight /= total;
        arr2[i].weight = arr1[i].weight;
    }

    return 0.0f;
}

static int testDistributionEMD(const char* distName, EMDTestDistribFunc distribf, unsigned int n)
{
    unsigned int i;
    WeightPos* arr1;
    WeightPos* arr2;
    float expected;
    float actual;
    int differs;

    arr1 = mwCalloc(n, sizeof(WeightPos));
    arr2 = mwCalloc(n, sizeof(WeightPos));

    for (i = 0; i < n; ++i)
    {
        arr2[i].pos = arr1[i].pos = (float) i;
    }

    expected = distribf(arr1, arr2, n);
    actual = emdCalc((const float*) arr1, (const float*) arr2, n, n, NULL);

    free(arr1);
    free(arr2);

    differs = (fabsf(expected - actual) >= ZERO_THRESHOLD);

    if (differs)
    {
        mw_printf("EMD different for test distribution '%s' with %u bins:\n"
                  "  Expected %f, Actual %f, |Diff| = %f\n",
                  distName,
                  n,
                  expected,
                  actual,
                  fabsf(actual - expected)
            );
    }
    else
    {
        printf("EMD test '%s'[%u] = %f, %f\n", distName, n, expected, actual);
    }

    return differs;
}

int main(int argc, const char* argv[])
{
    int fails = 0;

    dsfmt_init_gen_rand(&_prng, (uint32_t) time(NULL));


    fails += testDistributionEMD("randomSelf", randomSelf, 1);
    fails += testDistributionEMD("randomSelf", randomSelf, 7);
    fails += testDistributionEMD("randomSelf", randomSelf, 20);
    fails += testDistributionEMD("randomSelf", randomSelf, 33);
    fails += testDistributionEMD("randomSelf", randomSelf, 34);
    fails += testDistributionEMD("randomSelf", randomSelf, 50);

    fails += testDistributionEMD("allInSameBin", allInSameBin, 1);
    fails += testDistributionEMD("allInSameBin", allInSameBin, 7);
    fails += testDistributionEMD("allInSameBin", allInSameBin, 20);
    fails += testDistributionEMD("allInSameBin", allInSameBin, 33);
    fails += testDistributionEMD("allInSameBin", allInSameBin, 34);
    fails += testDistributionEMD("allInSameBin", allInSameBin, 50);

    fails += testDistributionEMD("oppositeSides", oppositeSides, 1);
    fails += testDistributionEMD("oppositeSides", oppositeSides, 7);
    fails += testDistributionEMD("oppositeSides", oppositeSides, 20);
    fails += testDistributionEMD("oppositeSides", oppositeSides, 33);
    fails += testDistributionEMD("oppositeSides", oppositeSides, 34);
    fails += testDistributionEMD("oppositeSides", oppositeSides, 50);

    fails += testDistributionEMD("allInDifferentBins", allInDifferentBins, 1);
    fails += testDistributionEMD("allInDifferentBins", allInDifferentBins, 7);
    fails += testDistributionEMD("allInDifferentBins", allInDifferentBins, 20);
    fails += testDistributionEMD("allInDifferentBins", allInDifferentBins, 33);
    fails += testDistributionEMD("allInDifferentBins", allInDifferentBins, 34);
    fails += testDistributionEMD("allInDifferentBins", allInDifferentBins, 50);

    if (fails != 0)
    {
        mw_printf("%d EMD test distributions failed\n", fails);
    }

    return fails;
}

