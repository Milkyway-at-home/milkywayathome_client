/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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

/* Quick test for comparing basic math library function output on
 * different platforms */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>
#include <sys/resource.h>

//#include <fdlibm.h>
#include <math.h>
//#include <crlibm.h>

/* Windows doesn't have drand48, so do stupid things */
#define RANDOM_DOUBLE (rand() / RAND_MAX)


double* fillRandoms(unsigned int n)
{
    unsigned int i;
    double* arr = malloc(sizeof(double) * n);

    for (i = 0; i < n; ++i)
        arr[i] = RANDOM_DOUBLE;
    /* arr[i] = drand48(); */

    return arr;
}

double get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}


#define MAKE_TEST_fdlibm(funcname, mfunc, strname)              \
    void (funcname)(const unsigned int n)                       \
    {                                                           \
        double* randoms;                                        \
        double* results;                                        \
        double ts, te;                                          \
        unsigned int i;                                         \
        randoms = fillRandoms(n);                               \
        results = malloc(sizeof(double) * n);                   \
        ts = get_time();                                        \
        for (i = 0; i < n; ++i)                                 \
            results[i] = (mfunc)(randoms[i]);                   \
        te = get_time();                                        \
        printf("Time for fdlibm "strname" = %g\n", te - ts);    \
        free(randoms);                                          \
        free(results);                                          \
    }


#if CRLIBM
#define MAKE_TEST_crlibm(funcname, mfunc, strname)              \
    void (funcname)(const unsigned int n)                       \
    {                                                           \
        double* randoms;                                        \
        double* results;                                        \
        double ts, te;                                          \
        unsigned int i;                                         \
        randoms = fillRandoms(n);                               \
        results = malloc(sizeof(double) * n);                   \
        ts = get_time();                                        \
        for (i = 0; i < n; ++i)                                 \
            results[i] = (mfunc)(randoms[i]);                   \
        te = get_time();                                        \
        printf("Time for crlibm "strname" = %g\n", te - ts);    \
        free(randoms);                                          \
        free(results);                                          \
    }
#else
  #define MAKE_TEST_crlibm(funcname, mfunc, strname)
#endif


MAKE_TEST_fdlibm(test_fdlibm_cbrt, cbrt, "cbrt")

#if CRLIBM
void test_crlibm_pow32(const unsigned int n)
{
    double* randoms;
    double* results;
    double ts, te;
    unsigned int i;
    randoms = fillRandoms(n);
    results = malloc(sizeof(double) * n);
    ts = get_time();
    for (i = 0; i < n; ++i)
        results[i] = (pow_rn)(randoms[i], 1.5);
    te = get_time();
    printf("Time for crlibm pow(x, 3/2) = %g\n", te - ts);
    free(randoms);
    free(results);
}

#else
  #define test_crlibm_pow32(n)
#endif


void test_fdlibm_pow32(const unsigned int n)
{
    double* randoms;
    double* results;
    double ts, te;
    unsigned int i;
    randoms = fillRandoms(n);
    results = malloc(sizeof(double) * n);
    ts = get_time();
    for (i = 0; i < n; ++i)
        results[i] = pow(randoms[i], 1.5);
    te = get_time();
    printf("Time for fdlibm pow(x, 3/2) = %g\n", te - ts);
    free(randoms);
    free(results);
}

#define cube(x) ((x) * (x) * (x))

void test_sqrt_cube(const unsigned int n)
{
    double* randoms;
    double* results;
    double ts, te;
    unsigned int i;
    randoms = fillRandoms(n);
    results = malloc(sizeof(double) * n);
    ts = get_time();
    for (i = 0; i < n; ++i)
    {
        results[i] = sqrt(cube(randoms[i]));
    }
    te = get_time();
    printf("Time for sqrt(x^3) = %g\n", te - ts);
    free(randoms);
    free(results);
}

#if CRLIBM
MAKE_TEST_crlibm(test_crlibm_log1p, log1p_rn, "log1p")
MAKE_TEST_crlibm(test_crlibm_exp, exp_rn, "exp")
MAKE_TEST_crlibm(test_crlibm_expm1, expm1_rn, "expm1")
MAKE_TEST_crlibm(test_crlibm_sin, sin_rn, "sin")
MAKE_TEST_crlibm(test_crlibm_cos, cos_rn, "cos")
#endif


MAKE_TEST_fdlibm(test_fdlibm_log1p, log1p,    "log1p")
MAKE_TEST_fdlibm(test_fdlibm_exp, exp,    "exp")
MAKE_TEST_fdlibm(test_fdlibm_expm1, expm1,    "expm1")
MAKE_TEST_fdlibm(test_fdlibm_sin, sin,    "sin")
MAKE_TEST_fdlibm(test_fdlibm_cos, cos,    "cos")






int main(int argc, char** argv)
{
    unsigned int n = 10000000;
    long seed = 0;

    if (argc >= 2)
        n = strtod(argv[1], NULL);

    /* srand48(seed); */
    srand(seed);


    //test_crlibm_exp(n);
    test_fdlibm_exp(n);
    printf("\n");

    //test_crlibm_log1p(n);
    test_fdlibm_log1p(n);
    printf("\n");

    //test_crlibm_sin(n);
    test_fdlibm_sin(n);
    printf("\n");

    //test_crlibm_cos(n);
    test_fdlibm_cos(n);
    printf("\n");

    //test_crlibm_expm1(n);
    test_fdlibm_expm1(n);
    printf("\n");

    //test_crlibm_pow32(n);
    test_fdlibm_pow32(n);
    test_sqrt_cube(n);
    printf("\n");

    //runTest_crlibm_sqrt(n);
    //runTest_fdlibm_sqrt(n);

    return 0;
}

