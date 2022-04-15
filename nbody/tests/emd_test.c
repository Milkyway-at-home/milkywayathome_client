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

#define ZERO_THRESHOLD 1.0e-4

/* Function which assigns a sample distribution to arr1 and arr2 to
 * match of size n. Returns expected EMD for the distribution */
typedef float (*EMDTestDistribFunc)(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2, unsigned int n);

/* Distance metric used by EMD code, currently uses L2 standard metric */
static inline float distMetric(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2, unsigned int i, unsigned int j)
{
    /* Calculate Euclidean Norm */
    float lambda = arr1[i].lambda - arr2[j].lambda;
    float beta = arr1[i].beta - arr2[j].beta;

    return sqrtf((lambda * lambda) + (beta * beta));
}

/* Compare two floats using ZERO_THRESHOLD */            
static inline int floatsDiffer(float value1, float value2)
{
    return (fabsf(value1 - value2) >= ZERO_THRESHOLD);
}

static void randomDist(WeightPos* RESTRICT arr, unsigned int n)
{
    unsigned int i;
    double total = 0.0f;
    double rand;

    /* Generate random histogram */
    for (i = 0; i < n; ++i)
    {
        rand = dsfmt_genrand_open_open(&_prng);
        arr[i].weight = (float) rand;
        total += rand;
    }

    /* Normalize it */
    for (i = 0; i < n; ++i) 
    {
        arr[i].weight = arr[i].weight / (float) total;
    }

    /* Check what the total of all the bins is, should be 1 */
    total = 0.0f;
    for (i = 0; i < n; ++i) 
    {
        total += arr[i].weight;
    }

    if(floatsDiffer(1.0f, total))
    {
        mw_printf("WARNING: Sum of all bins != 1,  %f\n", total);
    }
}

/* Two different histograms with opposite corners set to 1  */
static float oppositeSides(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2, unsigned int n)
{
    arr1[0].weight = 1.0f;
    arr2[n - 1].weight = 1.0f;

    /* Distance between opposite corners */
    return distMetric(arr1, arr2, 0, n-1);
}

// Two different histograms with a random bin in each set to 1 */
static float allInDifferentBins(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2, unsigned int n)
{
    unsigned int i, j;

    i = (unsigned int) mwXrandom(&_prng, 0, n);
    j = (unsigned int) mwXrandom(&_prng, 0, n);
    assert(i < n && j < n);

    arr1[i].weight = 1.0f;
    arr2[j].weight = 1.0f;

    /* Distance between the randomly selected bins */
    return distMetric(arr1, arr2, i, j);
}

/* Two identical histograms with a single bin set to 1 */
static float allInSameBin(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2, unsigned int n)
{
    unsigned int i;

    /* Pick one bin at random */
    i = (unsigned int) mwXrandom(&_prng, 0, n);
    assert(i < n);

    arr1[i].weight = arr2[i].weight = 1.0f;

    return 0.0f;
}

/* Two identical histograms with random weights*/
static float randomSelf(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2, unsigned int n)
{
    unsigned int i;

    /* Random normalized distribution */
    randomDist(arr1, n);

    /* Copy it to arr2 */
    for (i = 0; i < n; ++i)
    {
        arr2[i].weight = arr1[i].weight;
    }

    return 0.0f;
}


/* 
 * Write the position variables (lambda, beta) into the histogram
 * Increasing integers in each direction 
 *                    
 * e.g. (lambda,beta)
 * 
 * 0,0  1,0  2,0 ...
 * 0,1  1,1  2,1 ...
 * 0,2  1,2  2,2 ...
 * ...  ...  ... ... 
*/
static void generatePositions(WeightPos* RESTRICT arr1, WeightPos* RESTRICT arr2,
                              unsigned int dim1, unsigned int dim2)
{
    unsigned int i;
    unsigned int j;
    unsigned int k;

    for (i = 0; i < dim1; ++i)
    {
        for (j = 0; j < dim2; j++)
        {
            k = i * dim2 + j;
            arr1[k].lambda = (float) i; 
            arr2[k].lambda = (float) i;
            arr1[k].beta = (float) j;
            arr2[k].beta = (float) j;
        }
    }
}

/* Create two random histograms, compute the EMD value twice and compare the values */
/* Used to verify that MAX_ITERATIONS is large enough for large histograms */
static int testConsistentEMD(unsigned int dim1, unsigned int dim2) 
{
    unsigned int n = dim1 * dim2;
    WeightPos* arr1;
    WeightPos* arr2;
    float result1;
    float result2;
    int differs;

    arr1 = mwCalloc(n, sizeof(WeightPos));
    arr2 = mwCalloc(n, sizeof(WeightPos));

    generatePositions(arr1, arr2, dim1, dim2);

    randomDist(arr1, n);
    randomDist(arr2, n);
    
    result1 = emdCalc((const float*) arr1, (const float*) arr2, n, n, NULL);
    result2 = emdCalc((const float*) arr1, (const float*) arr2, n, n, NULL);

    differs = floatsDiffer(result1, result2);

    if(differs) 
    {
        mw_printf("ERROR: EMD returned inconsistent results with %u x %u bins:\n"
                  "  Result1 %f, Result2 %f, |Diff| = %f\n",
                  dim1, dim2,
                  result1, result2, fabsf(result1 - result2)
            );
    }
    else
    {
        mw_printf("EMD test [%u,%u] %-20s = %f, %f\n",
                  dim1, dim2, "consistent", result1, result2
            );
    }

    return differs;
}   
 
/* Test expected values for basic distributions */   
static int testDistributionEMD(const char* distName, EMDTestDistribFunc distribf,
                               unsigned int dim1, unsigned int dim2)
{
    unsigned int n = dim1 * dim2;
    WeightPos* arr1;
    WeightPos* arr2;
    float expected;
    float actual;
    int differs;

    arr1 = mwCalloc(n, sizeof(WeightPos));
    arr2 = mwCalloc(n, sizeof(WeightPos));

    generatePositions(arr1, arr2, dim1, dim2);

    expected = distribf(arr1, arr2, n);
    actual = emdCalc((const float*) arr1, (const float*) arr2, n, n, NULL);

    free(arr1);
    free(arr2);

    differs = floatsDiffer(expected, actual);

    if (differs)
    {
        mw_printf("ERROR: EMD different for test distribution '%s' with %u x %u bins:\n"
                  "  Expected %f, Actual %f, |Diff| = %f\n",
                  distName,
                  dim1, dim2,
                  expected,
                  actual,
                  fabsf(actual - expected)
            );
    }
    else
    { 
        mw_printf("EMD test [%u,%u] %-20s = %f, %f\n", 
                  dim1, dim2, distName, expected, actual);
    }

    return differs;
}

int runTestsEMD(unsigned int dim1, unsigned int dim2)
{
    int fails = 0;

    /* TODO randomSelf fails sometimes with a value slightly above zero */
    fails += testDistributionEMD("randomSelf", randomSelf, dim1, dim2);
    fails += testDistributionEMD("allInSameBin", allInSameBin, dim1, dim2);    
    fails += testDistributionEMD("oppositeSides", oppositeSides, dim1, dim2);
    fails += testDistributionEMD("allInDifferentBins", allInDifferentBins, dim1, dim2);

    fails += testConsistentEMD(dim1, dim2);

    return fails;
}

int main(int argc, const char* argv[])
{
    int fails = 0;

    dsfmt_init_gen_rand(&_prng, (uint32_t) time(NULL));

    fails += runTestsEMD(1, 1);
    fails += runTestsEMD(1, 7);
    fails += runTestsEMD(7, 1);
    fails += runTestsEMD(7, 7);

    fails += runTestsEMD(11, 11);
    fails += runTestsEMD(11, 34);
    fails += runTestsEMD(34, 11);

    if (fails != 0)
    {
        mw_printf("%d EMD test distributions failed\n", fails);
    }

    return fails;
}

