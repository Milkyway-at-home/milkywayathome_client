/*
 * Copyright (c) 2020 Rensselaer Polytechnic Institute
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "milkyway_util.h"
#include "nbody_virial.h"
#include <stdio.h> 
#include <stdlib.h>
#include <time.h>

int main()
{
    srand(time(0));

    int test_fails=0;

    const int nTests = 100;
    const real ZERO_THRESHOLD = 1.0e-5;

    real a_b, a_d, M_b, M_d, flip, lhs, rhs, comp;
    real min_a = 0.1;
    real max_a = 5.0;
    real min_M = 1.0;
    real max_M = 20.0;

    for (int j=0; j<nTests; j++) {
        a_b = ((real)rand()/(real)RAND_MAX)*(max_a - min_a) + min_a;
        a_d = ((real)rand()/(real)RAND_MAX)*(max_a - min_a) + min_a;
        M_b = ((real)rand()/(real)RAND_MAX)*(max_M - min_M) + min_M;
        M_d = ((real)rand()/(real)RAND_MAX)*(max_M - min_M) + min_M;

        flip = ((real)rand()/(real)RAND_MAX);

        if (flip < 0.5) {
            M_b = 0.0;
            lhs = 16.0*a_d/3/3.1415926535897932;
        }
        else {
            M_d = 0.0;
            lhs = 16.0*a_b/3/3.1415926535897932;
        }

        rhs = nbCalculateVirial(a_b, a_d, M_b, M_d);
        comp = mw_abs_0(lhs-rhs);

        if (comp > ZERO_THRESHOLD) {
            mw_printf("LHS = %.15f\n", lhs);
            mw_printf("RHS = %.15f\n", rhs);
            test_fails += 1;
            break;
        }
    }

    return test_fails;
}
