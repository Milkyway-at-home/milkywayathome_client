/*
Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

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

#include "milkyway.h"
#include "milkyway_complex.h"

#define CUBE(x) x * x * x

double complex ccbrt(double complex z)
{
    unsigned int idx;
    double r, err, errp, theta;
    double complex roots[3], root;

    if (z == COMPLEXZERO)
        return COMPLEXZERO;

    r = cabs(z);
    r = pow(r, 1.0 / 3.0);
    theta = carg(z) / 3.0;

    roots[0] =  r * (cos(theta) + sin(theta) * I);
    root = CUBE(roots[0]);

    err = cabs(root - z);
    idx = 0;

    roots[1] = r * ((cos(theta + PI_2_3)) + (sin(theta + PI_2_3)) * I);
    root = CUBE(roots[1]);

    errp = cabs(root - z);
    if (err > errp)
    {
        err = errp;
        idx = 1;
    }

    roots[2] = r * ((cos(theta + PI_4_3)) + (sin(theta + PI_4_3)) * I);
    root = CUBE(roots[2]);

    errp = cabs(root - z);
    if (err > errp)
    {
        err = errp;
        idx = 2;
    }

    r = -r;

    roots[0] = r * (cos(theta) + sin(theta) * I);
    root = CUBE(roots[0]);

    errp = cabs(root - z);
    if (err > errp)
    {
        err = errp;
        idx = 3;
    }

    roots[1] = r * (cos(theta + PI_2_3) + sin(theta + PI_2_3) * I);
    root = CUBE(roots[1]);

    errp = cabs(root - z);
    if (err > errp)
    {
        err = errp;
        idx = 4;
    }

    roots[2] = r * (cos(theta + PI_4_3) + sin(theta + PI_4_3) * I);
    root = CUBE(roots[2]);

    errp = cabs(root - z);
    if (err > errp)
    {
        err = errp;
        idx = 5;
    }

    if (idx > 2)
        return roots[idx-3];

    return -roots[idx];
}

