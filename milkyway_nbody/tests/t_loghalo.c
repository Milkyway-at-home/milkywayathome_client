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

#define RANDOM_LOG_HALO { LogarithmicHalo, RANDOM_REAL, RANDOM_REAL, RANDOM_REAL, NAN, NAN, NAN }


/* Test the falloff property of the log halo potential */
int log_halo_falloff(Halo* h)
{
    vector a0, a1, scaleda0, rv1;

    const real v = h->vhalo;
    const real d = h->scale_length;
    const real q = h->flattenZ;

    /* Take a random position */
    real r;
    vector rv = RANDOM_VECTOR;
    ABSV(r, rv);

    /* Find the acceleration there */
    logHaloAccel(a0, h, rv);

    /* If you scale the radius by an arbitrary factor k */
    real k = RANDOM_REAL;
    MULVS(rv1, rv, k);

    logHaloAccel(a1, h, rv1);

    /* Each component is scaled differently */

    /*

               k (d^2 + x^2 +y^2 + (z/q)^2)
      a_1x = -------------------------------- * a_0x
             d^2 + k^2 (x^2 + y^2 + (z/q)^2)

             - same for y component

                   k(q^2 (d^2 + x^2 +y^2) + z^2)
      a_1z = ----------------------------------------- * a_0z
               d^2q^2 + k^2 (q^2 (x^2 + y^2) + z^2)

     */


    /* Calculate the expected potential */

    real zqsr = sqr(rv[2])/sqr(q);

    real xynum = k * ( sqr(d) + sqr(rv[0]) + sqr(rv[1]) + zqsr);
    real xyden = sqr(d) + sqr(k) * (sqr(rv[0]) + sqr(rv[1]) + zqsr);
    real xyfactor = xynum / xyden;

    real znum = k * (sqr(q) * (sqr(d) + sqr(rv[0]) + sqr(rv[1])) + sqr(rv[2]));
    real zden = sqr(d * q) + sqr(k) * ( sqr(q) * (sqr(rv[0]) + sqr(rv[1])) + sqr(rv[2]));
    real zfactor = znum / zden;


    scaleda0[0] = xyfactor * a0[0];
    scaleda0[1] = xyfactor * a0[1];
    scaleda0[2] = zfactor  * a0[2];

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
                "Logarithmic halo potential falloff: Result differs significantly from expected value:\n"
                "\trv        = %s\n"
                "\tr         = %g\n"
                "\ta0        = %s\n"
                "\tq         = %g\n"
                "\td         = %g\n"
                "\tvhalo     = %g\n"
                "\tk         = %g\n"
                "\txyfactor  = %g/%g = %g\n"
                "\tzfactor   = %g/%g = %g\n"
                "\tExpected  = %s\n"
                "\tGot       = %s\n"
                "\tdiff      = %s\n"
                "\t|diff|    = %g\n"
                "\tTolerance = %g\n" ,
                rvStr, r, a0Str, q, d, v, k, xynum, xyden, xyfactor, znum, zden, zfactor, scaledStr, a1Str, diffStr, diffMag, LIMIT);

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
        Halo h = RANDOM_LOG_HALO;
        fails += log_halo_falloff(&h);
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


