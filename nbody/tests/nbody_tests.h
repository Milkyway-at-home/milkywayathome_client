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

#ifndef _NBODY_TESTS_H_
#define _NBODY_TESTS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define LIMIT 0.00001         /* arbitrary */

#define test_xrandom(xl, xh) ((real) (xl) + (real) ((xh) - (xl)) * drand48())

#define APPROXEQ(x, y) (abs(x - y) < LIMIT)

#define RANDOM_REAL (test_xrandom(0.0, 1.0))
#define RANDOM_VECTOR { RANDOM_REAL, RANDOM_REAL, RANDOM_REAL }

#define RANDOM_INT(low, high) (rand() % (high - low + 1) + low)

#define REL_ERROR(expected, answer) (100 * ((answer) - (expected)) / (expected))

#define VECEQ(a,b) ((a)[0] == (b)[0] && (a)[1] == (b)[1] && (a)[2] == (b)[2])

#define VECAPPROXEQ(a,b) (APPROXEQ((a)[0], (b)[0]) && APPROXEQ((a)[1], (b)[1]) && APPROXEQ((a)[2], (b)[2]))

/* copy b to a */
#define COPYVECTOR(a,b) { (a)[0] = (b)[0]; (a)[1] = (b)[1]; (a)[2] = (b)[2]; }


#endif /* _NBODY_TESTS_H_ */

