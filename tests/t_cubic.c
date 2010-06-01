/*
Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo, Nathan
Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
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

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "stMath.h"

#define LIMIT 0.00001  /* arbitrary */
#define PROP_TEST_NUM 100000
#define NUM_TESTS 5

#define APPROXEQ(x, y) (cabs(x - y) < LIMIT)

#define FSTRT 3
#define SNDRT 4
#define THDRT 5

#define RANDOM_DOUBLE (((double) rand()) / (((double) (RAND_MAX)) + (double) 1))
#define RANDOM_COMPLEX (RANDOM_DOUBLE + RANDOM_DOUBLE * I)


static const double complex tests[NUM_TESTS][7] =
{
    { 0.0 + 0.0 * I,  /* a */
      0.0 + 0.0 * I,  /* b */
      0.0 + 0.0 * I,  /* c */
      0.0 + 0.0 * I,  /* expected root 1 */
      0.0 + 0.0 * I,  /* expected root 2 */
      0.0 + 0.0 * I,  /* expected root 3 */
      0.0 + 0.0 * I   /* expected return value */
    },

    { 1.0 + 0.0 * I,
      0.0 + 0.0 * I,
      0.0 + 0.0 * I,
      0.0 + 0.0 * I,
      0.0 + 0.0 * I,
     -1.0 + 0.0 * I,
      0.0 + 0.0 * I
    },

    { 1.0 + 0.0 * I,
      1.0 + 0.0 * I,
      1.0 + 0.0 * I,
      0.0 + 1.0 * I,
      0.0 - 1.0 * I,
     -1.0 + 0.0 * I,
      0.0 + 0.0 * I
    },

    { 1.0 + 0.0 * I,
      2.0 + 0.0 * I,
      3.0 + 0.0 * I,
      -1.27568 + 0.0 * I,
      0.137841 + 1.52731 * I,
      0.137841 - 1.52731 * I,
      0.0 + 0.0 * I
    },


    { 3.0 + 0.0 * I,
      2.0 + 0.0 * I,
      1.0 + 0.0 * I,
      -2.32472 + 0.0 * I,
      -0.337641 - 0.56228 * I,
      -0.337641 + 0.56228 * I,
      0.0 + 0.0 * I
    }
};

/* for each root r, r^3 + a*r^2 + b*r + c == 0 */
unsigned int prop_tests()
{
    unsigned int i, j, failCount = 0;
    double complex roots[3];
    double complex a, b, c, x, x2, eval;
    int ret;

    for ( i = 0; i < PROP_TEST_NUM; ++i )
    {
        a = RANDOM_COMPLEX;
        b = RANDOM_COMPLEX;
        c = RANDOM_COMPLEX;
        ret = CnumCubic(a, b, c, roots, 0);

        for ( j = 0; j < 3; ++j )
        {
            x = roots[0];
            x2 = x * x;

            eval = (x * x2) + (a * x2) + (b * x) + c;

            if (!APPROXEQ(eval, COMPLEXZERO))
            {
                fprintf(stderr,
                        "CnumCubic failed property test:\n"
                        "  For a = (%g + %gI), b =  (%g + %gI), c = (%g + %gI)\n"
                        "  Got invalid root (%g + %gI), eval = %g\n",
                        creal(a), cimag(a),
                        creal(b), cimag(b),
                        creal(c), cimag(c),
                        creal(x), cimag(x),
                        eval);
                fprintf(stderr, "LOLOLO: %d\n", ret);
                ++failCount;
            }
        }

    }

    return failCount;
}

unsigned int run_tests()
{
    unsigned int i, j, failCount = 0;
    double complex roots[3];

    for ( i = 0; i < NUM_TESTS; ++i )
    {
        CnumCubic(tests[i][0], tests[i][1], tests[i][2], roots, 0);

        for ( j = 0; j < 3; ++j ) /* Look for each found root in expecteds */
        {

            if (    !APPROXEQ(roots[j], tests[i][FSTRT])
                 && !APPROXEQ(roots[j], tests[i][SNDRT])
                 && !APPROXEQ(roots[j], tests[i][THDRT])
               )
            {
                fprintf(stderr,
                        "CnumCubic failed test %d. Unexpected root (%g + %gI) found\n",
                        i,
                        creal(roots[j]),
                        cimag(roots[j]));
                ++failCount;
            }

        }
    }

    return failCount;

}

int main()
{
    unsigned int failCount = 0, propFails = 0, totalFails;

    failCount = run_tests();
    if ( failCount > 0 )
        fprintf(stderr, "t_cubic: %u out of %u tests failed.\n", failCount, NUM_TESTS);

    propFails = prop_tests();
    if ( propFails > 0 )
        fprintf(stderr, "t_cubic: %u out of %u property tests failed.\n", propFails, PROP_TEST_NUM);

    totalFails = failCount + propFails;

    if (totalFails == 0)
        printf("t_cubic: All tests passed.\n");
    else
        fprintf(stderr, "t_cubic: %u tests failed.\n", totalFails);

    return totalFails;
}

