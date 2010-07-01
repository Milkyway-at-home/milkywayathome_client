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

#define RANDOM_TRIAXIAL_HALO { TriaxialHalo, RANDOM_REAL, RANDOM_REAL, RANDOM_REAL, RANDOM_REAL, RANDOM_REAL, RANDOM_REAL, RANDOM_REAL, RANDOM_REAL, RANDOM_REAL }


/* Test the falloff property of the log halo potential */
static int triaxial_halo_falloff(Halo* h)
{
    vector a0, a1, scaleda0, rv1;

    /* Take a random position */
    real r;
    vector rv = RANDOM_VECTOR;
    ABSV(r, rv);

    /* Find the acceleration there */
    triaxialHaloAccel(a0, h, rv);

    /* If you scale the radius by an arbitrary factor k */
    real k = RANDOM_REAL;
    MULVS(rv1, rv, k);

    triaxialHaloAccel(a1, h, rv1);

    /*
            let l = (c1 * x^2) + (c3 * x * y) + (c2 * y^2)

                k * ( qz^2 * (rhalo^2 + l) + z^2 )
        a1 = ------------------------------------------ * a0
               qz^2 (rhalo^2 + k^2 * l) + (k^2 * z^2)

     */


    /* Calculate the expected potential */

    const real phi   = h->triaxAngle;
    const real vhalo = h->vhalo;
    const real rhalo = h->scale_length;
    const real qx    = h->flattenX;
    const real qy    = h->flattenY;
    const real qz    = h->flattenZ;

    const real cp  = cos(phi);
    const real cps = sqr(cp);
    const real sp  = sin(phi);
    const real sps = sqr(sp);

    const real qxs = sqr(qx);
    const real qys = sqr(qy);
    const real qzs = sqr(qz);

    const real c1 = (cps / qxs) + (sps / qys);
    const real c2 = (cps / qys) + (sps / qxs);

    const real c3 = 2 * sin(phi) * cos(phi) * (1/qxs - 1/qys);


    const real l = (c1 * sqr(rv[0])) + (c3 * rv[0] * rv[1]) + (c2 * sqr(rv[1]));

    const real numer  = k * ( qzs * ( sqr(rhalo) + l ) + sqr(rv[2]) );
    const real denom  = qzs * (sqr(rhalo) + sqr(k) * l) + (sqr(k) * sqr(rv[2]));
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
                "Triaxial halo potential falloff: Result differs significantly from expected value:\n"
                "\trv        = %s\n"
                "\tr         = %g\n"
                "\ta0        = %s\n"
                "\tqx        = %g\n"
                "\tqy        = %g\n"
                "\tqz        = %g\n"
                "\trhalo     = %g\n"
                "\tvhalo     = %g\n"
                "\tc1        = %g\n"
                "\tc2        = %g\n"
                "\tc3        = %g\n"
                "\tk         = %g\n"
                "\tfactor    = %g/%g = %g\n"
                "\tExpected  = %s\n"
                "\tGot       = %s\n"
                "\tdiff      = %s\n"
                "\t|diff|    = %g\n"
                "\tTolerance = %g\n" ,
                rvStr, r, a0Str, qx, qy, qz, rhalo, vhalo, c1, c2, c3, k, numer, denom, factor, scaledStr, a1Str, diffStr, diffMag, LIMIT);

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
        Halo h = RANDOM_TRIAXIAL_HALO;
        fails += triaxial_halo_falloff(&h);
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


