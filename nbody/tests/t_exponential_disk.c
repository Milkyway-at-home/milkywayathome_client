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
#include "nbody_priv.h"

#define RANDOM_EXPONENTIAL_DISK { ExponentialDisk, RANDOM_REAL, RANDOM_REAL, RANDOM_REAL }


int exponential_disk_falloff(Disk* d)
{
    vector a0, a1, scaleda0, rv1;

    /* Take a random position */
    real r;
    vector rv = RANDOM_VECTOR;
    ABSV(r, rv);

    /* Find the acceleration there */
    exponentialDiskAccel(a0, d, rv);

    /* If you scale the radius by an arbitrary factor k */
    real k = RANDOM_REAL;
    MULVS(rv1, rv, k);

    exponentialDiskAccel(a1, d, rv1);

    /*
           let s = exp(- k * r / b) * (b + k * r) / b
               u = exp(-r / b) * (b + r) / b

                 s - 1
      a1 = ------------------ * a0
             k^2 * (u - 1)

     */


    /* Calculate the expected potential */
    const real m = d->mass;
    const real b = d->scale_length;

    const real s = rexp(-k * r / b) * (b + k * r) / b;
    const real u = rexp(-r / b) * (b + r) / b;

    const real numer  = s - 1;
    const real denom  = sqr(k) * (u - 1);
    const real factor = numer / denom;


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
                "Exponential disk potential falloff: Result differs significantly from expected value:\n"
                "\trv        = %s\n"
                "\tr         = %g\n"
                "\ta0        = %s\n"
                "\tmass      = %g\n"
                "\tb         = %g\n"
                "\tk         = %g\n"
                "\tfactor  = %g/%g = %g\n"
                "\ts         = %g\n"
                "\tu         = %g\n"
                "\tExpected  = %s\n"
                "\tGot       = %s\n"
                "\tdiff      = %s\n"
                "\t|diff|    = %g\n"
                "\tTolerance = %g\n" ,
                rvStr, r, a0Str, m, b, k, numer, denom, factor, s, u, scaledStr, a1Str, diffStr, diffMag, LIMIT);

        free(scaledStr);
        free(a1Str);
        free(a0Str);
        free(rvStr);
        free(diffStr);

        return 1;
    }

    return 0;
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
        Disk d = RANDOM_EXPONENTIAL_DISK;
        fails += exponential_disk_falloff(&d);
    }

    totalTests = numTests;

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


