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

#define RANDOM_MN_DISK { MiyamotoNagaiDisk, RANDOM_REAL, RANDOM_REAL, RANDOM_REAL }


int mn_disk_falloff(Disk* d)
{
    vector a0, a1, scaleda0, rv1;

    const real m = d->mass;
    const real a = d->scale_length;
    const real b = d->scale_height;

    /* Take a random position */
    real r;
    vector rv = RANDOM_VECTOR;
    ABSV(r, rv);

    /* Find the acceleration there */
    miyamotoNagaiDiskAccel(a0, d, rv);

    /* If you scale the radius by an arbitrary factor k */
    real k = RANDOM_REAL;
    MULVS(rv1, rv, k);

    miyamotoNagaiDiskAccel(a1, d, rv1);

    /* Each component is scaled differently */

    /*
      let h = a + sqrt(b^2 + z^2)
          i = a + sqrt(b^2 + k^2 * z^2)
          w = x^2 + y^2 + h^2

                        k * w^(3/2)
      a_1x = ----------------------------------- * a_0x
               (k^2 * (x^2 + y^2) + i^2)^(3/2)

               same for y


                            k * sqrt(b^2 + z^2) * i * w^(3/2)
     a_1z = ------------------------------------------------------------------- * a_0z
               sqrt(b^2 + k^2 * z^2) * h * (k^2 * (x^2 + y^2) + i^2)^(3/2)

     */


    /* Calculate the expected potential */

    real h = a + sqrt( sqr(b) + sqr(rv[2]) );
    real i = a + sqrt( sqr(b) + sqr(k) * sqr(rv[2]) );
    real w = sqr(rv[0]) + sqr(rv[1]) + sqr(h);

    real xynum = k * pow(sqr(rv[0]) + sqr(rv[1]) + sqr(h), 1.5);
    real xyden = pow( sqr(k) * (sqr(rv[0]) + sqr(rv[1])) + sqr(i), 1.5);
    real xyfactor = xynum / xyden;

    real znum = k * sqrt( sqr(b) + sqr(rv[2]) ) * i * pow(w, 1.5);
    real zden = sqrt(sqr(b) + sqr(k) * sqr(rv[2])) * h * pow( sqr(k) * (sqr(rv[0]) + sqr(rv[1])) + sqr(i), 1.5);

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
                "Miyamoto-Nagai disk potential falloff: Result differs significantly from expected value:\n"
                "\trv        = %s\n"
                "\tr         = %g\n"
                "\ta0        = %s\n"
                "\tmass      = %g\n"
                "\ta         = %g\n"
                "\tb         = %g\n"
                "\tk         = %g\n"
                "\txyfactor  = %g/%g = %g\n"
                "\tzfactor   = %g/%g = %g\n"
                "\tExpected  = %s\n"
                "\tGot       = %s\n"
                "\tdiff      = %s\n"
                "\t|diff|    = %g\n"
                "\tTolerance = %g\n" ,
                rvStr, r, a0Str, m, a, b, k, xynum, xyden, xyfactor, znum, zden, zfactor, scaledStr, a1Str, diffStr, diffMag, LIMIT);

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
        Disk d = RANDOM_MN_DISK;
        fails += mn_disk_falloff(&d);
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


