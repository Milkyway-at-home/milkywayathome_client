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
static void crossp( const double* a, const double* b, double* prod )
{
    prod[0] = a[1] * b[2] - a[2] * b[1];
    prod[1] = a[2] * b[0] - a[0] * b[2];
    prod[2] = a[0] * b[1] - a[1] * b[0];
}

/* Dot product */
double dotp( const double* a, const double* b )
{
    return a[0] * b[0] +
           a[1] * b[1] +
           a[2] * b[2];
}

/* Get norm of input vector */
static double norm( const double* vec )
{
    return sqrt( vec[0] * vec[0] +
                 vec[1] * vec[1] +
                 vec[2] * vec[2] );
}

/* Normalize input vector */
static void normalize( double* vec )
{
    double vnorm = norm( vec );

    vec[0] /= vnorm;
    vec[1] /= vnorm;
    vec[2] /= vnorm;
}

/* Angle between two vectors, in the range [0,pi] */
static double vecangle( const double* a, const double* b )
{
    double anorm, bnorm, dprod;

    anorm = norm( a );
    bnorm = norm( b );
    dprod = dotp( a, b );

    return acos( dprod / (anorm * bnorm) );
}


/* Convert sun-centered lbr into galactic xyz coordinates. */
void lbr2xyz_old(const double* lbr, double* xyz)
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
static void xyz2lbr(const double* xyz, double* lbr)
{
    double temp, xsun;

    xsun = xyz[0] + 8.5;
    temp = xsun * xsun + xyz[1] * xyz[1];

    lbr[0] = atan2( xyz[1], xsun ) * deg;
    lbr[1] = atan2( xyz[2], sqrt( temp ) ) * deg;
    lbr[2] = sqrt( temp + xyz[2] * xyz[2] );

    if ( lbr[0] < 0 ) lbr[0] += 360;
}

static void xyz2lbg(double* point, double offset, double* lbg)
{
    xyz2lbr(point, lbg);
    double g = 5.0 * (log(100.0 * lbg[2]) / log(10.0) ) + 4.2 - offset;

    lbg[2] = g;
}


/* wrapper that converts a point into magnitude-space pseudo-xyz */
static void xyz_mag(double* point, double offset, double* logPoint)
{
    double lbg[3];
    xyz2lbg(point, offset, lbg);

    lbr2xyz_old(lbg, logPoint);
}




static void slaDmxv ( double dm[3][3], double va[3], double vb[3] )
{
    int i, j;
    double w, vw[3];

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
    double* angle,    /* MODIFIED -- the angle to bound in degrees*/
    double min,   /* IN -- inclusive minimum value */
    double max    /* IN -- exclusive maximum value */
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

static double slaDranrm ( double angle )
{
    double w;

    w = dmod ( angle, M_2PI );
    return ( w >= 0.0 ) ? w : w + M_2PI;
}

static double slaDrange ( double angle )
{
    double w;

    w = dmod ( angle, M_2PI );
    return ( fabs ( w ) < M_PI ) ? w : w - dsign ( M_2PI, angle );
}

static void atBound2(
    double* theta,    /* MODIFIED -- the -90 to 90 angle */
    double* phi   /* MODIFIED -- the 0 to 360 angle */
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
static void atSurveyToEq (double slong, double slat, double* ra, double* dec)
{
    double anode, etaPole;
    double x1, y1, z1;

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

static void slaDcc2s( double v[3], double* a, double* b )
{
    double x, y, z, r;

    x = v[0];
    y = v[1];
    z = v[2];
    r = sqrt ( x * x + y * y );

    *a = ( r != 0.0 ) ? atan2 ( y, x ) : 0.0;
    *b = ( z != 0.0 ) ? atan2 ( z, r ) : 0.0;
}


static void slaDcs2c( double a, double b, double v[3] )
{
    double cosb;

    cosb = cos ( b );
    v[0] = cos ( a ) * cosb;
    v[1] = sin ( a ) * cosb;
    v[2] = sin ( b );
}


static void slaEqgal( double dr, double dd, double* dl, double* db )
{
    double v1[3], v2[3];

    static double rmat[3][3];

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
int prob_ok(int n, double* p)
{
    int ok;
    double r;
    double step1, step2, step3;

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
static void lbToXyz(double l, double b, double* xyz)
{
    l = l / deg;
    b = b / deg;

    xyz[0] = cos(l) * cos(b);
    xyz[1] = sin(l) * cos(b);
    xyz[2] = sin(b);
}

/* Return eta from stripe number */
static double atEtaFromStripeNumber(int wedge)
{
    double eta;

    if (wedge <= 46)
    {
        eta = wedge * stripeSeparation - 57.5;
    }
    else
    {
        eta = wedge * stripeSeparation - 57.5 - 180.0;
    }

    return eta;
}

static void atEqToGal (
    double ra,  /* IN -- ra in degrees */
    double dec, /* IN -- dec in degrees */
    double* glong,  /* OUT -- Galactic longitude in degrees */
    double* glat    /* OUT -- Galactic latitude in degrees */
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
void get_transform( const double* f, const double* t, double** mat )
{
    double angle, sin_a;
    double x, y, z, w;
    double axis[3];

    crossp( f, t, axis );
    normalize( axis );

    angle = vecangle( f, t );
    sin_a = sin( angle / 2 );

    x = axis[0] * sin_a;
    y = axis[1] * sin_a;
    z = axis[2] * sin_a;
    w = cos( angle / 2 );

    mat[0][0] = 1 - 2 * (y * y + z * z);
    mat[0][1] =     2 * (x * y - z * w);
    mat[0][2] =     2 * (x * z + y * w);

    mat[1][0] =     2 * (x * y + z * w);
    mat[1][1] = 1 - 2 * (x * x + z * z);
    mat[1][2] =     2 * (y * z - x * w);

    mat[2][0] =     2 * (x * z - y * w);
    mat[2][1] =     2 * (y * z + x * w);
    mat[2][2] = 1 - 2 * (x * x + y * y);
}


/* Transform v by applying the rotation matrix mat */
static void do_transform(double* v, const double** mat)
{
    double newv[3];

    newv[0] = dotp( mat[0], v );
    newv[1] = dotp( mat[1], v );
    newv[2] = dotp( mat[2], v );

    v[0] = newv[0];
    v[1] = newv[1];
    v[2] = newv[2];
}

/* apply coordinate transformations to the given point */
void transform_point(double* point, const double** cmat, double* xsun, double* logPoint)
{
    double mcutoff = 11.0;

    xyz_mag(point, mcutoff, logPoint);

    double newx = logPoint[0] - xsun[0];
    double newy = logPoint[1] - xsun[1];
    double newz = logPoint[2] - xsun[2];
    logPoint[0] = newx;
    logPoint[1] = newy;
    logPoint[2] = newz;

    do_transform(logPoint, cmat);
}

/* Initialize seed for prob_ok */
void prob_ok_init()
{
    srand48(time(NULL));
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

