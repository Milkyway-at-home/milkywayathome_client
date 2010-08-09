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

#ifdef _WIN32
#include <float.h>
#endif

#include "atSurveyGeometry.h"
#include "stCoords.h"

/* Convert sun-centered lbr into galactic xyz coordinates. */
void lbr2xyz(const double* lbr, vector xyz)
{
    double zp, d;

    const double bsin = sin(B(lbr) / deg);
    const double lsin = sin(L(lbr) / deg);
    const double bcos = cos(B(lbr) / deg);
    const double lcos = cos(L(lbr) / deg);

    Z(xyz) = R(lbr) * bsin;
    zp = R(lbr) * bcos;
    d = sqrt( sun_r0 * sun_r0 + zp * zp - 2.0 * sun_r0 * zp * lcos);
    X(xyz) = (zp * zp - sun_r0 * sun_r0 - d * d) / (2.0 * sun_r0);
    Y(xyz) = zp * lsin;
}

/* Get eta for the given wedge. */
double wedge_eta(int wedge)
{
    return wedge * at_stripeSeparation - 57.5 - (wedge > 46 ? 180.0 : 0.0);
}

/* Get inclination for the given wedge. */
double wedge_incl(int wedge)
{
    return wedge_eta(wedge) + at_surveyCenterDec;
}

/* Convert GC coordinates (mu, nu) into l and b for the given wedge. */
void gc2lb( int wedge, double mu, double nu, double* l, double* b )
{
    RA_DEC radec = atGCToEq( mu, nu, wedge_incl( wedge ) );
    atEqToGal( radec.ra, radec.dec, l, b );
}

void gc2sgr( int wedge, double mu, double nu, double* l, double* b )
{
    double lamda, beta;
    gcToSgr(wedge, mu, nu, &lamda, &beta);
    sgrToGal(lamda, beta, l, b);
}

