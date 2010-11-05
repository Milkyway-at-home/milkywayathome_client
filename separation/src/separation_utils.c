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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "separation.h"
#include "integrals_common.h"
#include "milkyway_math.h"


#define at_rad2Deg (180.0 / M_PI)
#define at_deg2Rad (M_PI / 180.0)
#define deg (180.0 / M_PI)
#define r0 8.5

/* Cross product; stores result in prod */
static void crossp( const real* a, const real* b, real* prod )
{
    prod[0] = a[1] * b[2] - a[2] * b[1];
    prod[1] = a[2] * b[0] - a[0] * b[2];
    prod[2] = a[0] * b[1] - a[1] * b[0];
}

/* Dot product */
real dotp( const real* a, const real* b )
{
    return a[0] * b[0] +
           a[1] * b[1] +
           a[2] * b[2];
}

/* Get norm of input vector */
static real norm( const real* vec )
{
    return sqrt( vec[0] * vec[0] +
                 vec[1] * vec[1] +
                 vec[2] * vec[2] );
}

/* Normalize input vector */
static void normalize( real* vec )
{
    real vnorm = norm( vec );

    vec[0] /= vnorm;
    vec[1] /= vnorm;
    vec[2] /= vnorm;
}

/* Convert sun-centered lbr into galactic xyz coordinates. */
void lbr2xyz_old(const real* lbr, real* xyz)
{
    real bsin, lsin, bcos, lcos, zp, d;

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
static void xyz2lbr(const real* xyz, real* lbr)
{
    real temp, xsun;

    xsun = xyz[0] + 8.5;
    temp = xsun * xsun + xyz[1] * xyz[1];

    lbr[0] = atan2( xyz[1], xsun ) * deg;
    lbr[1] = atan2( xyz[2], sqrt( temp ) ) * deg;
    lbr[2] = sqrt( temp + xyz[2] * xyz[2] );

    if ( lbr[0] < 0 ) lbr[0] += 360;
}

static void xyz2lbg(real* point, real offset, real* lbg)
{
    xyz2lbr(point, lbg);
    real g = 5.0 * (log(100.0 * lbg[2]) / log(10.0) ) + 4.2 - offset;

    lbg[2] = g;
}

void marshal_mwvector_to_array(real* arr, mwvector v);
mwvector marshal_array_to_mwvector(real* arr);


/* wrapper that converts a point into magnitude-space pseudo-xyz */
static mwvector xyz_mag(const ASTRONOMY_PARAMETERS* ap, real* point, real offset)
{
    mwvector logPoint;
    real lbg[3];
    xyz2lbg(point, offset, lbg);

    mwvector arst = marshal_array_to_mwvector(lbg);

    logPoint = lbr2xyz(ap, arst);

    return logPoint;
}




static void slaDmxv(real dm[3][3], real va[3], real vb[3])
{
    int i, j;
    real w, vw[3];

    /* Matrix dm * vector va -> vector vw */
    for ( j = 0; j < 3; j++ )
    {
        w = 0.0;
        for ( i = 0; i < 3; i++ )
        {
            w += dm[j][i] * va[i];
        }
        vw[j] = w;
    }

    /* Vector vw -> vector vb */
    for ( j = 0; j < 3; j++ )
    {
        vb[j] = vw[j];
    }
}

static void atBound (
    real* angle,    /* MODIFIED -- the angle to bound in degrees*/
    real min,   /* IN -- inclusive minimum value */
    real max    /* IN -- exclusive maximum value */
)
{
    while (*angle < min)
    {
        *angle += 360.0;
    }
    while (*angle >= max)
    {
        *angle -= 360.0;
    }
    return;
}

static real slaDranrm ( real angle )
{
    real w;

    w = dmod ( angle, M_2PI );
    return ( w >= 0.0 ) ? w : w + M_2PI;
}

static real slaDrange ( real angle )
{
    real w;

    w = dmod ( angle, M_2PI );
    return ( fabs ( w ) < M_PI ) ? w : w - dsign ( M_2PI, angle );
}

static void atBound2(
    real* theta,    /* MODIFIED -- the -90 to 90 angle */
    real* phi   /* MODIFIED -- the 0 to 360 angle */
)
{
    atBound(theta, -180.0, 180.0);
    if (fabs(*theta) > 90.0)
    {
        *theta = 180.0 - *theta;
        *phi += 180;
    }
    atBound(theta, -180.0, 180.0);
    atBound(phi, 0.0, 360.0);
    if (fabs(*theta) == 90.0) *phi = 0.;
    return;
}

/* Return ra & dec from survey longitude and latitude */
static void atSurveyToEq (real slong, real slat, real* ra, real* dec)
{
    real anode, etaPole;
    real x1, y1, z1;

    /* Convert to radians */
    slong = slong * at_deg2Rad;
    slat = slat * at_deg2Rad;
    anode = surveyCenterRa - 90.0;
    anode = anode * at_deg2Rad;
    etaPole = surveyCenterDec * at_deg2Rad;

    /* Rotation */
    x1 = -sin(slong);
    y1 = cos(slat + etaPole) * cos(slong);
    z1 = sin(slat + etaPole) * cos(slong);
    *ra = atan2(y1, x1) + anode;
    *dec = asin(z1);
    *ra = *ra * at_rad2Deg;
    *dec = *dec * at_rad2Deg;
    atBound2(dec, ra);

    return;
}

static void slaDcc2s( real v[3], real* a, real* b )
{
    real x, y, z, r;

    x = v[0];
    y = v[1];
    z = v[2];
    r = sqrt ( x * x + y * y );

    *a = ( r != 0.0 ) ? atan2 ( y, x ) : 0.0;
    *b = ( z != 0.0 ) ? atan2 ( z, r ) : 0.0;
}


static void slaDcs2c( real a, real b, real v[3] )
{
    real cosb;

    cosb = cos ( b );
    v[0] = cos ( a ) * cosb;
    v[1] = sin ( a ) * cosb;
    v[2] = sin ( b );
}


static void slaEqgal( real dr, real dd, real* dl, real* db )
{
    real v1[3], v2[3];

    static real rmat[3][3];

    rmat[0][0] = -0.054875539726;
    rmat[0][1] = -0.873437108010;
    rmat[0][2] = -0.483834985808;
    rmat[1][0] =  0.494109453312;
    rmat[1][1] = -0.444829589425;
    rmat[1][2] =  0.746982251810;
    rmat[2][0] = -0.867666135858;
    rmat[2][1] = -0.198076386122;
    rmat[2][2] =  0.455983795705;

    /* Spherical to Cartesian */
    slaDcs2c ( dr, dd, v1 );

    /* Equatorial to Galactic */
    slaDmxv ( rmat, v1, v2 );

    /* Cartesian to spherical */
    slaDcc2s ( v2, dl, db );

    /* Express in conventional ranges */
    *dl = slaDranrm ( *dl );
    *db = slaDrange ( *db );
}


/* determines if star with prob p should be separrated into stream */
int prob_ok(int n, real* p)
{
    int ok;
    real r;
    real step1, step2, step3;

    r = drand48();

    switch (n)
    {
    case 1:
        if ( r > p[0] )
        {
            ok = 0;
        }
        else
        {
            ok = 1;
        }
        break;
    case 2:
        step1 = p[0] + p[1];
        if ( r > step1 )
        {
            ok = 0;
        }
        else if ( (r < p[0]) )
        {
            ok = 1;
        }
        else if ( (r > p[0]) && (r <= step1) )
        {
            ok = 2;
        }
        break;
    case 3:
        step1 = p[0] + p[1];
        step2 = p[0] + p[1] + p[2];
        if ( r > step2 )
        {
            ok = 0;
        }
        else if ( (r < p[0]) )
        {
            ok = 1;
        }
        else if ( (r > p[0]) && (r <= step1) )
        {
            ok = 2;
        }
        else if ( (r > step1) && (r <= step2) )
        {
            ok = 3;
        }
        break;
    case 4:
        step1 = p[0] + p[1];
        step2 = p[0] + p[1] + p[2];
        step3 = p[0] + p[1] + p[2] + p[3];
        if ( r > step3 )
        {
            ok = 0;
        }
        else if ( (r <= p[0]) )
        {
            ok = 1;
        }
        else if ( (r > p[0]) && (r <= step1) )
        {
            ok = 2;
        }
        else if ( (r > step1) && (r <= step2) )
        {
            ok = 3;
        }
        else if ( (r > step2) && (r <= step3) )
        {
            ok = 4;
        }
        break;
    default:
        fail("ERROR:  Too many streams to separate using current code; please update the switch statement in probability.c->prob_ok to handle %d streams", n);
    }
    return ok;
}

/* convert galactic coordinates l,b into cartesian x,y,z */
static void lbToXyz(real l, real b, real* xyz)
{
    l = l / deg;
    b = b / deg;

    xyz[0] = cos(l) * cos(b);
    xyz[1] = sin(l) * cos(b);
    xyz[2] = sin(b);
}

/* Return eta from stripe number */
static real atEtaFromStripeNumber(int wedge)
{
    if (wedge > 46)
        return wedge * stripeSeparation - 57.5 - 180.0;
    else
        return wedge * stripeSeparation - 57.5;
}

static void atEqToGal (
    real ra,  /* IN -- ra in degrees */
    real dec, /* IN -- dec in degrees */
    real* glong,  /* OUT -- Galactic longitude in degrees */
    real* glat    /* OUT -- Galactic latitude in degrees */
)
{
    /* Convert to radians */
    ra = ra * at_deg2Rad;
    dec = dec * at_deg2Rad;
    /* Use SLALIB to do the actual conversion */
    slaEqgal(ra, dec, glong, glat);
    /* Convert back to degrees */
    *glong = *glong * at_rad2Deg;
    *glat = *glat * at_rad2Deg;
    atBound2(glat, glong);
    return;
}


/* Determine the rotation matrix needed to transform f into t.  The result is an
   array of 9 elements representing a (flattened) 3x3 matrix.
   Adapted from information at http://www.flipcode.com/documents/matrfaq.html */
void get_transform(mwmatrix mat, const mwvector f, const mwvector t)
{
    real angle, sin_a;
    real x, y, z, w;
    mwvector axis;

    axis = mw_crossv(f, t);
    mw_normalize(axis);

    angle = mw_vecangle(f, t);
    sin_a = mw_sin(0.5 * angle);

    x = X(axis) * sin_a;
    y = Y(axis) * sin_a;
    z = Z(axis) * sin_a;
    w = mw_cos(0.5 * angle);

    X(mat[0]) = 1.0 - 2.0 * (y * y + z * z);
    Y(mat[0]) =       2.0 * (x * y - z * w);
    Z(mat[0]) =       2.0 * (x * z + y * w);

    X(mat[1]) =       2.0 * (x * y + z * w);
    Y(mat[1]) = 1.0 - 2.0 * (x * x + z * z);
    Z(mat[1]) =       2.0 * (y * z - x * w);

    X(mat[2]) =       2.0 * (x * z - y * w);
    Y(mat[2]) =       2.0 * (y * z + x * w);
    Z(mat[2]) = 1.0 - 2.0 * (x * x + y * y);
}


/* Transform v by applying the rotation matrix mat */
/* apply coordinate transformations to the given point */
mwvector transform_point(const ASTRONOMY_PARAMETERS* ap, real* point, const mwmatrix cmat, mwvector xsun)
{
    real mcutoff = 11.0;
    mwvector logPoint = xyz_mag(ap, point, mcutoff);

    mw_incsubv(logPoint, xsun);

    return mw_mulmv(cmat, logPoint); /* do transform */
}

/* Initialize seed for prob_ok */
void prob_ok_init()
{
    srand48(time(NULL));
}

/* Get normal vector of data slice from stripe number */
void stripe_normal(int wedge, real* xyz)
{
    real eta, ra, dec, l, b;

    eta = atEtaFromStripeNumber(wedge);
    atSurveyToEq(0, 90.0 + eta, &ra, &dec);
    atEqToGal(ra, dec, &l, &b);
    lbToXyz(l, b, xyz);
}

