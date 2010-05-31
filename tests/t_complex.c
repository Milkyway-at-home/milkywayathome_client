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
#include "milkyway_complex.h"

#define LIMIT 0.00001  /* arbitrary */
#define PROP_TEST_NUM 100000
#define NUM_TESTS 20

#define APPROXEQ(x, y) (cabs(x - y) < LIMIT)
                                     /* input , expected result */
static const double complex tests[NUM_TESTS][2] = { {0.0 + 0.0 * I, 0.0 + 0.0 * I},
                                                    {1.0 + 0.0 * I, 1.0 + 0.0 * I},
                                                    {0.0 + 1.0 * I, 0.866025 + 0.5 * I},
                                                    {-1.0 + 0.0 * I, 0.5 + 0.866025 * I},
                                                    {-1.0 - 1.0 * I, 0.793701 - 0.793701 * I},
                                                    {0.0 - 1.0 * I, 0.866025 - 0.5 * I},
                                                    {1 + 1 * I, 1.08422 + 0.290515 * I},
                                                    {3.0 + 4.0 * I, 1.62894 + 0.520175 * I},
                                                    {5.0 + 5.0 * I, 1.85398 + 0.496773 * I},
                                                    {7.0 + 2.0 * I, 1.92978 + 0.179534 * I},
                                                    {0.3 + 1.1 * I, 0.947472 + 0.440102 * I},
                                                    {1.1 + 2.2 * I, 1.25899 + 0.486938 * I},
                                                    {3.0 + 9.0 * I, 1.93609 + 0.856138 * I},
                                                    {0.452502 + 0.0410748 * I, 0.768428 + 0.0231942 * I},
                                                    {0.0523431 + 0.844328 * I, 0.828637 + 0.455877 * I},
                                                    {0.177803 + 0.228583 * I, 0.631415 + 0.197561 * I},
                                                    {0.0344606 + 0.345007 * I, 0.619723 + 0.330883 * I},
                                                    {0.660476 + 0.713974 * I, 0.953622 + 0.268824 * I},
                                                    {0.283584 + 0.114257 * I, 0.66819 + 0.0857739 * I},
                                                    {0.488534 + 0.576144 * I, 0.872918 + 0.259698 * I}
                                                  };


#define RANDOM_DOUBLE (((double) rand()) / (((double) (RAND_MAX)) + (double) 1))


/* CHECKME: This might not be true in general or something or other,
   something about -1 and ( (-1)^3 ) ^ (1/3) != -1, but this is the other way.
 */
/*run tests using the property that (x^1/3)^3 == x */
unsigned int prop_eq(unsigned int num)
{
    double complex cr, ret;
    unsigned int i;
    unsigned int fails = 0;

    for ( i = 0; i < num; ++i )
    {
        cr = RANDOM_DOUBLE + RANDOM_DOUBLE * I;
        ret = ccbrt(cr);
        ret = ret * ret * ret;

        if (cabs(ret - cr) > LIMIT)
            ++fails;

    }

    return fails;
}

unsigned int run_tests()
{
    double complex result;
    unsigned int i, failCount = 0;

    for ( i = 0; i < NUM_TESTS; ++i )
    {
        result = ccbrt(tests[i][0]);
        if (!APPROXEQ(result, tests[i][1]))
        {
            fprintf(stderr,
                    "ccbrt failed test %d: Expected (%g + %gI) got (%g + %gI). |difference| = %g\n",
                    i,
                    creal(tests[i][1]),
                    cimag(tests[i][1]),
                    creal(result),
                    cimag(result),
                    cabs(result - tests[i][1]));

            ++failCount;
        }
    }

    return failCount;

}

int main()
{
    unsigned int totalFails;
    unsigned int propFails = 0;
    unsigned int failCount = 0;

    failCount = run_tests();

    if ( failCount > 0 )
        fprintf(stderr, "t_complex: %u out of %u tests failed.\n", failCount, NUM_TESTS);

    propFails = prop_eq(PROP_TEST_NUM);

    if ( propFails > 0 )
        fprintf(stderr, "t_complex: %u out of %u property tests failed.\n", propFails, PROP_TEST_NUM);

    totalFails = failCount + propFails;

    if (totalFails == 0)
        fprintf(stderr, "t_complex: All tests passed.\n");
    else
        fprintf(stderr, "t_complex: %u tests failed.\n", totalFails);

    return totalFails;
}

