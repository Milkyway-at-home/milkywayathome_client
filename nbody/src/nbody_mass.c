/*
 * Copyright (c) 2012 Rensselaer Polytechnic Institute
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

#include "nbody_mass.h"

#include "milkyway_math.h"
#include "nbody_types.h"


/*In order to decrease the size of the numbers
 * computed all these functions are
 * calculated in log space*/

static real factorial(int n)
{
    if (n <= 1)
    {
        return 0;
    }

    return mw_log(n) + factorial(n - 1);
}

static real choose(int n, int c)
{
    unsigned int i;
    real result = 0;

    /* This for loop calulates log(n!/(n-c)!) */
    for (i = n - c + 1; i <= (unsigned int) n; ++i)
    {
        result += mw_log(i);
    }

    result -= factorial(c);
    return result;
}

real probability_match(int n, int k, real pobs)
{
    real result;
    result = choose(n, k);
    result += mw_log(pobs) * (real) k;
    result += mw_log(1.0 - pobs) * (real)(n - k);
    return result;
}
