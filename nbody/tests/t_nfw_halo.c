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
#include "nbody_priv.h"

#define RANDOM_NFW_HALO { NFWHalo, RANDOM_REAL, RANDOM_REAL, NAN, NAN, NAN, NAN }

int nfw_halo_falloff(Halo* h)
{
    vector a0, a1, scaleda0, rv1;

    /* Take a random position */
    real r;
    vector rv = RANDOM_VECTOR;
    ABSV(r, rv);

    /* Find the acceleration there */
    nfwHaloAccel(a0, h, rv);

    /* If you scale the radius by an arbitrary factor k */
    real k = RANDOM_REAL;
    MULVS(rv1, rv, k);

    nfwHaloAccel(a1, h, rv1);

    /*

               (a + r) * (kr - (a + kr) * log(1 + kr/a)
      a1 = ----------------------------------------------- * a0
               k^2 (a +kr) (r - (a+r) log( (a+r) / a))

     */


    /* Calculate the expected potential */
    const real v = h->vhalo;
    const real a = h->scale_length;

    const real ar  = a + r;
    const real akr = a + k * r;

    const real arst1 = ar / a;
    const real arst2 = k * r / a;

    const real numer = ar * ( k * r - akr * log1p(arst2));
    const real denom = sqr(k) * akr * ( r - ar * log(arst1));

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
                "NFW halo potential falloff: Result differs significantly from expected value:\n"
                "\trv        = %s\n"
                "\tr         = %g\n"
                "\ta0        = %s\n"
                "\tscale a   = %g\n"
                "\tvhalo     = %g\n"
                "\tk         = %g\n"
                "\tfactor  = %g/%g = %g\n"
                "\tExpected  = %s\n"
                "\tGot       = %s\n"
                "\tdiff      = %s\n"
                "\t|diff|    = %g\n"
                "\tTolerance = %g\n" ,
                rvStr, r, a0Str, a, v, k, numer, denom, factor, scaledStr, a1Str, diffStr, diffMag, LIMIT);

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
        Halo h = RANDOM_NFW_HALO;
        fails += nfw_halo_falloff(&h);
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


