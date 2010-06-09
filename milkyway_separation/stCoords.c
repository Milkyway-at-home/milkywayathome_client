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
#include "stVector.h"

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
    d = sqrt( r0 * r0 + zp * zp - 2 * r0 * zp * lcos);
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
        lbr[0] += 360;
}

/* Convert stream coordinates to lbr given the specified stream parameters. */
void stream2lbr(const double* stream, const double* spars, double* lbr)
{
    double xyz[3];
    stream2xyz( stream, spars, xyz );
    xyz2lbr( xyz, lbr );
}


/* Same as stream2lbr, but with output in galactic-centered xyz. */
void stream2xyz(const double* stream, const double* spars, double* xyz)
{
    int i;
    const double* a;
    const double* b;
    const double* c;
    double anorm, bnorm, cost, sint;
    double ex[3], ey[3];

    c = &spars[0];
    a = &spars[3];
    b = &spars[6];

    anorm = norm(a);
    bnorm = norm(b);
    cost = cos(stream[0]);
    sint = sin(stream[0]);

    /* Get x-axis and y-axis normal vectors. */
    for (i = 0; i < 3; ++i)
        ex[i] = bnorm / anorm * a[i] * cost + anorm / bnorm * b[i] * sint;

    crossp(a, b, ey);

    normalize(ex);
    normalize(ey);

    for (i = 0; i < 3; ++i)
        xyz[i] = c[i] + a[i] * cost + b[i] * sint + ex[i] * stream[1] + ey[i] * stream[2];
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

/* Get normal vector of data slice from stripe number */
void stripe_normal( int wedge, double* xyz )
{
    double eta, ra, dec, l, b;

    eta = atEtaFromStripeNumber(wedge);
    atSurveyToEq(0, 90.0 + eta, &ra, &dec);
    atEqToGal(ra, dec, &l, &b);
    lbToXyz(l, b, xyz);
}

//vickej2 change made to account for sgr_Stripes, calculates normal vector by crossmultiplication
void sgr_stripe_normal(int wedge, double* xyz)
{
    double lamda1, beta1, lamda2, beta2, l1, b1, l2, b2, xyz1[3], xyz2[3];
    lamda1 = wedge * 2.5;
    lamda2 = wedge * 2.5;

    beta1 = 0;
    beta2 = 90;

    sgrToGal(lamda1, beta1, &l1, &b1);
    sgrToGal(lamda2, beta2, &l2, &b2);

    lbToXyz(l1, b1, xyz1);
    lbToXyz(l2, b2, xyz2);

    //printf("lamda=%f, beta=%f, l=%f, b=%f, x=%f, y=%f, z=%f", lamda1, beta1, l1, b1, xyz1[0], xyz1[1], xyz1[2]);

//crossmultiplication of the 2 vectors
    xyz[0] = xyz1[1] * xyz2[2] - xyz1[2] * xyz2[1];
    xyz[1] = xyz1[2] * xyz2[0] - xyz1[0] * xyz2[2];
    xyz[2] = xyz1[0] * xyz2[1] - xyz1[1] * xyz2[0];
}

/*convert galactic coordinates l,b into cartesian x,y,z*/
void lbToXyz(double l, double b, double* xyz)
{
    l = l / deg;
    b = b / deg;

    xyz[0] = cos(l) * cos(b);
    xyz[1] = sin(l) * cos(b);
    xyz[2] = sin(b);
}

/*wrapper that converts a point into magnitude-space pseudo-xyz*/
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

