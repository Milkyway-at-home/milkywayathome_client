/* Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

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

#include "nbody_tests.h"
#include "nbody.h"

/* These are supposed be the max enum value. These will change as new
 * potentials are added */
#define RANDOM_SPHERICAL { RANDOM_INT(0, SphericalPotential), RANDOM_REAL, RANDOM_REAL }
#define RANDOM_DISK { RANDOM_INT(0, ExponentialDisk), RANDOM_REAL, RANDOM_REAL, RANDOM_REAL }
#define RANDOM_HALO { RANDOM_INT(0, TriaxialHalo), RANDOM_REAL, RANDOM_REAL, RANDOM_REAL }


/* Test the falloff property of the spherical potential */
int spherical_falloff(Spherical* s)
{
    vector a0, a1, scaleda0, rv1;

    const real r0 = s->scale;

    /* Take a random position */
    real r;
    vector rv = RANDOM_VECTOR;
    ABSV(r, rv);

    /* Find the acceleration there */
    sphericalAccel(a0, s, rv);

    /* If you scale the radius by an arbitrary factor k */
    real k = 10 * RANDOM_REAL;
    MULVS(rv1, rv, k);

    /* the new acceleration should be

             ( r^(1/2) + r0 )^2
          ----------------------- * (original acceleration)
            (k * r^(1/2) + r0)^2
    */

    /* Find the acceleration at the scaled vector */
    sphericalAccel(a1, &s, rv1);

    /* Calculate the factor */
    real numer = sqr( sqrt(r) + r0 );
    real denom = sqr( k * sqrt(r) + r0 );
    real factor = numer / denom;

    /* Calculate the expected potential */
    MULVS(scaleda0, a0, factor);

    if (!VECAPPROXEQ(scaleda0, a1))
    {
        char* scaledStr = showVector(scaleda0);
        char* a1Str     = showVector(a1);
        char* a0Str     = showVector(a0);
        char* rvStr     = showVector(rv);

        vector diff;
        real diffMag;
        SUBV(diff, scaleda0, a1);
        ABSV(diffMag, diff);

        char* diffStr = showVector(diff);

        fprintf(stderr,
                "Spherical potential falloff: Result differs significantly from expected value:\n"
                "\trv        = %s\n"
                "\tr         = %g\n"
                "\ta0        = %s\n"
                "\tr0        = %g\n"
                "\tk         = %g\n"
                "\tfactor    = %g/%g = %g\n"
                "\tExpected  = %s\n"
                "\tGot       = %s\n"
                "\tdiff      = %s\n"
                "\t|diff|    = %g\n"
                "\tTolerance = %g\n" ,
                rvStr, r, a0Str, r0, k, numer, denom, factor, scaledStr, a1Str, diffStr, diffMag, LIMIT);

        free(scaledStr);
        free(a1Str);
        free(a0Str);
        free(rvStr);
        free(diffStr);

        return 1;
    }

    return 0;
}


/* The acceleration should scale linearly with the mass */
int spherical_mass_scaling(Spherical* s)
{
    vector a0, a1, a0Scaled;
    int failed = 0;
    const real r0 = s->scale;
    const real m  = s->mass;

    /* Take a random position */
    real r;
    vector rv = RANDOM_VECTOR;
    ABSV(r, rv);

    /* Find the acceleration there */
    sphericalAccel(a0, s, rv);


    /* Take an arbitrary factor, and scale the mass of the bulge */
    real k = RANDOM_REAL;
    s->mass *= k;

    /* Find the new acceleration */
    sphericalAccel(a1, s, rv);


    /* The new acceleration should be the original scaled by k */
    MULVS(a0Scaled, a0, k);

    if (!VECAPPROXEQ(a0Scaled, a1))
    {
        char* scaledStr = showVector(a0Scaled);
        char* a1Str     = showVector(a1);
        char* rvStr     = showVector(rv);

        vector diff;
        real diffMag;
        SUBV(diff, a0Scaled, a1);
        ABSV(diffMag, diff);

        char* diffStr = showVector(diff);

        fprintf(stderr,
                "Spherical mass scaling: Result differs significantly from expected value:\n"
                "\trv        = %s\n"
                "\tr         = %g\n"
                "\tr0        = %g\n"
                "\tk         = %g\n"
                "\tExpected  = %s\n"
                "\tGot       = %s\n"
                "\tdiff      = %s\n"
                "\t|diff|    = %g\n"
                "\tTolerance = %g\n" ,
                rvStr, r, r0, k, scaledStr, a1Str, diffStr, diffMag, LIMIT);

        free(scaledStr);
        free(a1Str);
        free(rvStr);
        free(diffStr);

        failed = 1;
    }

    /* Restore the mass for future tests on this sample to use
     * unmodified */
    s->mass = m;

    return failed;
}


int spherical_r0_scaling(Spherical* s)
{
    vector a0, a1, a0Scaled;
    int failed = 0;
    const real r0 = s->scale;

    /* Take a random position */
    real r;
    vector rv = RANDOM_VECTOR;
    ABSV(r, rv);

    /* Find the acceleration there */
    sphericalAccel(a0, s, rv);

    /* Take an arbitrary factor, and scale the mass of the bulge */
    real k = RANDOM_REAL;
    s->scale *= k;

    /* Find the new acceleration */
    sphericalAccel(a1, s, rv);

    /*
                    (r^(1/2) + r0)^2
           a1 = ----------------------- * a0
                  (r^(1/2) + k * r0)^2
     */

    /* Calculate the factor */
    real numer = sqr( sqrt(r) + r0 );
    real denom = sqr( sqrt(r) + k * r0 );
    real factor = numer / denom;

    /* Calculate the expected potential */
    MULVS(a0Scaled, a0, factor);

    if (!VECAPPROXEQ(a0Scaled, a1))
    {
        char* scaledStr = showVector(a0Scaled);
        char* a1Str     = showVector(a1);
        char* rvStr     = showVector(rv);

        vector diff;
        real diffMag;
        SUBV(diff, a0Scaled, a1);
        ABSV(diffMag, diff);

        char* diffStr = showVector(diff);

        fprintf(stderr,
                "Spherical scale radius scaling: Result differs significantly from expected value:\n"
                "\trv        = %s\n"
                "\tr         = %g\n"
                "\tr0        = %g\n"
                "\tk         = %g\n"
                "\tfactor    = %g/%g = %g\n"
                "\tExpected  = %s\n"
                "\tGot       = %s\n"
                "\tdiff      = %s\n"
                "\t|diff|    = %g\n"
                "\tTolerance = %g\n" ,
                rvStr, r, r0, k, numer, denom, factor, scaledStr, a1Str, diffStr, diffMag, LIMIT);

        free(scaledStr);
        free(a1Str);
        free(rvStr);
        free(diffStr);

        failed = 1;
    }

    /* Restore the scale radius for future tests on this sample to use
     * unmodified */
    s->scale = r0;

    return failed;
}


/* One argument as number of tests to run */
int main(int argc, char** argv)
{
    size_t i;
    unsigned int fails = 0;
    unsigned int numTests;
    unsigned int totalTests;

    if (argc != 2)
        numTests = 1000000;
    else
        numTests = strtol(argv[1], NULL, 10);

    for (i = 0; i < numTests; ++i)
    {
        Spherical s = RANDOM_SPHERICAL;
        fails += spherical_falloff(&s);
        fails += spherical_mass_scaling(&s);
        fails += spherical_r0_scaling(&s);
    }

    totalTests = 3 * numTests;

    if (fails)
    {
        fprintf(stderr,
                "Acceleration tests: %u out of %u tests failed (%g%%)\n",
                fails,
                totalTests,
                100 * (double) fails / (double) totalTests);
    }

    return fails;
}


