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
      0.0 + 0.0 * I,
      0.0 + 0.0 * I
    },

    { 1.0 + 0.0 * I,
      1.0 + 0.0 * I,
      1.0 + 0.0 * I,
      0.0 + 0.0 * I,
      -0.5 - 0.866025 * I,
      -0.5 + 0.866025 * I,
      0.0 + 0.0 * I
    },

    { 1.0 + 0.0 * I,
      2.0 + 0.0 * I,
      3.0 + 0.0 * I,
      0.0 + 0.0 * I,
      -1.0 - 1.41421 * I,
      -1.0 + 1.41421 * I,
      0.0 + 0.0 * I
    },


    { 3.0 + 0.0 * I,
      2.0 + 0.0 * I,
      1.0 + 0.0 * I,
      0.0 + 0.0 * I,
      -0.333333 - 0.471405 * I,
      -0.333333 + 0.471405 * I,
      0.0 + 0.0 * I
    }
};


#define RANDOM_DOUBLE (((double) rand()) / (((double) (RAND_MAX)) + (double) 1))



unsigned int run_tests()
{
    int ret;
    unsigned int i, failCount = 0;
    double complex roots[3];

    unsigned int j;

    for ( i = 0; i < NUM_TESTS; ++i )
    {
        ret = CnumCubic(tests[i][0], tests[i][1], tests[i][2], roots, 0);

        printf("test %u: ret = %d  ", i, ret);
        for ( j = 0; j < 3; ++j )
        {
            printf("root[%u] = %g + %gI ", j, creal(roots[j]), cimag(roots[j]));
        }
        printf("\n");
    }

    return failCount;

}

int main()
{
    unsigned int failCount = 0;

    failCount = run_tests();

    if ( failCount > 0 )
        fprintf(stderr, "t_cubic: %u out of %u tests failed.\n", failCount, NUM_TESTS);


    return failCount;
}

