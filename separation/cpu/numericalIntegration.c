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
#include "milkyway_priv.h"

/* Gauss-Legendre quadrature taken from Numerical Recipes in C */
void gaussLegendre(double x1, double x2, double x[], double w[], int n)
{
    int m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;

    m = (n + 1) / 2;
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);

    for (i = 1; i <= m; i++)
    {
        z = cos(PI * (i - 0.25) / (n + 0.5));
        do
        {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 1; j <= n; j++)
            {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        }
        while (fabs(z - z1) > EPS);

        x[i-1] = xm - xl * z;
        x[n-i] = xm + xl * z;
        w[i-1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w[n-i] = w[i-1];
    }
}

