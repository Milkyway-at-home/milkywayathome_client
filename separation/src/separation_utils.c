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
#include "separation_utils.h"

/* Convert galactic xyz into sun-centered lbr coordinates. */
static mwvector xyz2lbr(const ASTRONOMY_PARAMETERS* ap, const mwvector xyz)
{
    mwvector lbr;
    real tmp, xsun;

    xsun = X(xyz) + ap->sun_r0;
    tmp = sqr(xsun) + sqr(Y(xyz));

    L(lbr) = mw_atan2(Y(xyz), xsun);
    B(lbr) = mw_atan2(Z(xyz), mw_sqrt(tmp));
    R(lbr) = mw_sqrt(tmp + sqr(Z(xyz)));

    L(lbr) = r2d(L(lbr));
    B(lbr) = r2d(B(lbr));

    if (L(lbr) < 0)
        L(lbr) += 360.0;

    return lbr;
}

static inline real calc_g(mwvector lbg)
{
    return 5.0 * (mw_log(100.0 * R(lbg)) / mw_log(10.0) ) + absm;
}

/* CHECKME: Isn't this supposed to be log10? */
static mwvector xyz2lbg(const ASTRONOMY_PARAMETERS* ap, mwvector point, real offset)
{
    mwvector lbg = xyz2lbr(ap, point);
    real g = calc_g(lbg) - offset;

    R(lbg) = g;

    return lbg;
}


/* wrapper that converts a point into magnitude-space pseudo-xyz */
static mwvector xyz_mag(const ASTRONOMY_PARAMETERS* ap, mwvector point, real offset)
{
    mwvector logPoint, lbg;

    lbg = xyz2lbg(ap, point, offset);
    logPoint = lbr2xyz(ap, lbg);

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
        *angle += 360.0;
    while (*angle >= max)
        *angle -= 360.0;
    return;
}

static real slaDranrm ( real angle )
{
    real w;

    w = dmod(angle, M_2PI);
    return (w >= 0.0) ? w : w + M_2PI;
}

static real slaDrange ( real angle )
{
    real w;

    w = dmod(angle, M_2PI);
    return (fabs(w) < M_PI) ? w : w - dsign(M_2PI, angle);
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
    if (fabs(*theta) == 90.0)
        *phi = 0.0;
    return;
}

/* Return ra & dec from survey longitude and latitude */
static void atSurveyToEq(real slong, real slat, real* ra, real* dec)
{
    real anode, etaPole;
    real x1, y1, z1;

    /* Convert to radians */
    slong = d2r(slong);
    slat = d2r(slat);
    anode = surveyCenterRa - 90.0;
    anode = d2r(anode);
    etaPole = d2r(surveyCenterDec);

    /* Rotation */
    x1 = -mw_sin(slong);
    y1 = mw_cos(slat + etaPole) * mw_cos(slong);
    z1 = mw_sin(slat + etaPole) * mw_cos(slong);
    *ra = mw_atan2(y1, x1) + anode;
    *dec = mw_asin(z1);
    *ra = r2d(*ra);
    *dec = r2d(*dec);
    atBound2(dec, ra);

    return;
}

static void slaDcc2s(real v[3], real* a, real* b)
{
    real x, y, z, r;

    x = v[0];
    y = v[1];
    z = v[2];
    r = mw_hypot(x, y);

    *a = (r != 0.0) ? mw_atan2(y, x) : 0.0;
    *b = (z != 0.0) ? mw_atan2(z, r) : 0.0;
}

static void slaDcs2c(real a, real b, real v[3])
{
    real cosb;

    cosb = mw_cos(b);
    v[0] = mw_cos(a) * cosb;
    v[1] = mw_sin(a) * cosb;
    v[2] = mw_sin(b);
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
    slaDcs2c( dr, dd, v1 );

    /* Equatorial to Galactic */
    slaDmxv( rmat, v1, v2 );

    /* Cartesian to spherical */
    slaDcc2s( v2, dl, db );

    /* Express in conventional ranges */
    *dl = slaDranrm(*dl);
    *db = slaDrange(*db);
}


/* convert galactic coordinates l,b into cartesian x,y,z */
static mwvector lbToXyz(real l, real b)
{
    mwvector xyz;
    l = d2r(l);
    b = d2r(b);

    X(xyz) = cos(l) * cos(b);
    Y(xyz) = sin(l) * cos(b);
    Z(xyz) = sin(b);

    return xyz;
}


static void atEqToGal(
    real ra,  /* IN -- ra in degrees */
    real dec, /* IN -- dec in degrees */
    real* glong,  /* OUT -- Galactic longitude in degrees */
    real* glat    /* OUT -- Galactic latitude in degrees */
)
{
    /* Convert to radians */
    ra = d2r(ra);
    dec = d2r(dec);

    /* Use SLALIB to do the actual conversion */
    slaEqgal(ra, dec, glong, glat);
    /* Convert back to degrees */
    *glong = r2d(*glong);
    *glat = r2d(*glat);
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
mwvector transform_point(const ASTRONOMY_PARAMETERS* ap,
                         mwvector point,
                         const mwmatrix cmat,
                         mwvector xsun)
{
    real mcutoff = 11.0;
    mwvector logPoint = xyz_mag(ap, point, mcutoff);

    mw_incsubv(logPoint, xsun);

    return mw_mulmv(cmat, logPoint); /* do transform */
}

/* Get normal vector of data slice from stripe number */
mwvector stripe_normal(int wedge)
{
    real eta, ra, dec, l, b;

    eta = atEtaFromStripeNumber_deg(wedge);
    atSurveyToEq(0, 90.0 + eta, &ra, &dec);
    atEqToGal(ra, dec, &l, &b);
    return lbToXyz(l, b);
}

/* Initialize seed for prob_ok; time based seed if 0 */
void prob_ok_init(long seed)
{
    if (seed)
        srand48(seed);
    else
        srand48(time(NULL));
}

/* FIXME: WTF? */
/* determines if star with prob p should be separrated into stream */
int prob_ok(StreamStats* ss, int n)
{
    int ok;
    real r;
    real step1, step2, step3;

    r = drand48();

    switch (n)
    {
        case 1:
            if (r > ss[0].sprob)
                ok = 0;
            else
                ok = 1;
            break;
        case 2:
            step1 = ss[0].sprob + ss[1].sprob;
            if (r > step1)
                ok = 0;
            else if (r < ss[0].sprob)
                ok = 1;
            else if (r > ss[0].sprob && r <= step1)
                ok = 2;
            break;
        case 3:
            step1 = ss[0].sprob + ss[1].sprob;
            step2 = ss[0].sprob + ss[1].sprob + ss[2].sprob;
            if (r > step2)
                ok = 0;
            else if (r < ss[0].sprob)
                ok = 1;
            else if (r > ss[0].sprob && r <= step1)
                ok = 2;
            else if (r > step1 && r <= step2)
                ok = 3;
            /* CHECME: else? */
            break;
        case 4:
            step1 = ss[0].sprob + ss[1].sprob;
            step2 = ss[0].sprob + ss[1].sprob + ss[2].sprob;
            step3 = ss[0].sprob + ss[1].sprob + ss[2].sprob + ss[3].sprob;
            if (r > step3)
                ok = 0;
            else if (r <= ss[0].sprob)
                ok = 1;
            else if (r > ss[0].sprob && r <= step1)
                ok = 2;
            else if (r > step1 && r <= step2)
                ok = 3;
            else if (r > step2 && r <= step3)
                ok = 4;
            break;
        default:
            fail("ERROR:  Too many streams to separate using current code; "
                 "please update the switch statement in prob_ok to handle %d streams", n);
    }
    return ok;
}

