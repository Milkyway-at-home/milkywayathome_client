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

static const double r0 = 8.5;


/* Convert sun-centered lbr into galactic xyz coordinates. */
void lbr2xyz(const double* lbr, double* xyz)
{
    double bsin, lsin, bcos, lcos, zp, d;

    bsin = sin(lbr[1] / deg);
    lsin = sin(lbr[0] / deg);
    bcos = cos(lbr[1] / deg);
    lcos = cos(lbr[0] / deg);

    xyz[2] = lbr[2] * bsin;
    zp = lbr[2] * bcos;
    d = sqrt( r0 * r0 + zp * zp - 2.0 * r0 * zp * lcos);
    xyz[0] = (zp * zp - r0 * r0 - d * d) / (2 * r0);
    xyz[1] = zp * lsin;
}


/* Convert galactic xyz into sun-centered lbr coordinates. */
void xyz2lbr(const double* xyz, double* lbr)
{
    double temp, xsun;

    xsun = xyz[0] + r0;
    temp = xsun * xsun + xyz[1] * xyz[1];

    lbr[0] = atan2( xyz[1], xsun ) * deg;
    lbr[1] = atan2( xyz[2], sqrt( temp ) ) * deg;
    lbr[2] = sqrt( temp + xyz[2] * xyz[2] );

    if ( lbr[0] < 0 )
        lbr[0] += 360.0;
}

/* Get eta for the given wedge. */
double wedge_eta(int wedge)
{
    double d;
    d = at_stripeSeparation;
    return wedge * d - 57.5 - (wedge > 46 ? 180.0 : 0.0);
}

/* Get inclination for the given wedge. */
double wedge_incl(int wedge)
{
    double d;
    d = at_surveyCenterDec;
    return wedge_eta(wedge) + d;
}

/* Get the node of the GC coordinates used in the survey. */
double get_node()
{
    double d;
    d = at_surveyCenterRa;
    return d - 90.0;
}

/* Convert GC coordinates (mu, nu) into l and b for the given wedge. */
void gc2lb( int wedge, double mu, double nu, double* l, double* b )
{
    double ra, dec;

    atGCToEq( mu, nu, &ra, &dec, get_node(), wedge_incl( wedge ) );
    atEqToGal( ra, dec, l, b );
}

/* wrapper that converts a point into magnitude-space pseudo-xyz */
void xyz_mag(double* point, double offset, double* logPoint)
{
    double lbg[3];
    xyz2lbg(point, offset, lbg);
    lbr2xyz(lbg, logPoint);
}


void xyz2lbg(double* point, double offset, double* lbg)
{
    xyz2lbr(point, lbg);
    double g = 5.0 * (log(100.0 * lbg[2]) / log(10.0) ) + 4.2 - offset;

    lbg[2] = g;
}


