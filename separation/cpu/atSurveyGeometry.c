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
#include "atSurveyGeometry.h"

void atGCToEq (
    double amu,  /* IN -- mu in degrees */
    double anu,  /* IN -- nu in degrees */
    double* ra,  /* OUT -- ra in degrees */
    double* dec, /* OUT -- dec in degrees */
    double anode,    /* IN -- node in degrees */
    double ainc  /* IN -- inclination in degrees */
)
{
    double x1, y1, z1, x2, y2, z2;
    /* Convert to radians */

    amu = amu * at_deg2Rad;
    anu = anu * at_deg2Rad;
    anode = anode * at_deg2Rad;
    ainc = ainc * at_deg2Rad;
    /* Rotation */
    x2 = cos(amu - anode) * cos(anu);
    y2 = sin(amu - anode) * cos(anu);
    z2 = sin(anu);
    x1 = x2;
    y1 = y2 * cos(ainc) - z2 * sin(ainc);
    z1 = y2 * sin(ainc) + z2 * cos(ainc);


    *ra = atan2 (y1, x1) + anode;
    *dec = asin(z1);
    /* Convert back to degrees */
    *ra = *ra * at_rad2Deg;
    *dec = *dec * at_rad2Deg;

//  printf("ra: %g, dec: %g\n", *ra, *dec);

    atBound2(dec, ra);
    return;
}

void atEqToGal (
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



void atBound (
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

void atBound2(
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

void slaDcc2s ( double v[3], double* a, double* b )
{
    double x, y, z, r;

    x = v[0];
    y = v[1];
    z = v[2];
    r = sqrt ( x * x + y * y );

    *a = ( r != 0.0 ) ? atan2 ( y, x ) : 0.0;
    *b = ( z != 0.0 ) ? atan2 ( z, r ) : 0.0;
}


void slaDcs2c ( double a, double b, vector v )
{
    double cosb;

    cosb = cos ( b );
    v[0] = cos ( a ) * cosb;
    v[1] = sin ( a ) * cosb;
    v[2] = sin ( b );
}


void slaDmxv ( double dm[3][3], vector va, vector vb )
{
    int i, j;
    double w;
    vector vw;

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


double slaDrange ( double angle )
{
    double w;

    w = dmod ( angle, D2PI );
    return ( fabs ( w ) < DPI ) ? w : w - dsign ( D2PI, angle );
}


double slaDranrm ( double angle )
{
    double w;

    w = dmod ( angle, D2PI );
    return ( w >= 0.0 ) ? w : w + D2PI;
}


void slaEqgal ( double dr, double dd, double* dl, double* db )
{
    vector v1;
    vector v2;

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

//vickej2 for sgr stripes, the great circles are defined thus:
//sgr stripes run parallel to sgr longitude lines, centered on lamda=2.5*wedge number
//with mu=0 at the sgr equator and increasing in the +z direction (increasing from the equator with beta)
//and nu=0 at the center and increasing in the -y direction (inversely to lamda)
//in this manner an equatorial stripe of standard coordinate conventions is created.
void gcToSgr ( double mu, double nu, int wedge, double* lamda, double* beta )
{
    double radpdeg = 3.141592653589793 / 180;
    mu = mu * radpdeg;
    nu = nu * radpdeg;

    double x = cos(mu) * cos(nu);
    double y = -sin(nu);
    double z = sin(mu) * cos(nu);

    *lamda = atan2(y, x);
    *lamda = *lamda / radpdeg;
    *lamda = *lamda + 2.5 * wedge;
    if (*lamda < 0)
    {
        *lamda = *lamda + 360;
    }

    *beta = asin(z);
    *beta = *beta / radpdeg;

    return;
}

#define sqr(x) ((x) * (x))

typedef struct
{
    double rot11, rot12, rot13;
    double rot21, rot22, rot23;
    double rot31, rot32, rot33;
} SGR_TO_GAL_CONSTANTS;

/* CHECKME: This mess needs testing, but I don't think it's actually used. */
void init_sgr_to_gal_constants(SGR_TO_GAL_CONSTANTS* sgc)
{
    const double radpdeg = M_PI / 180.0;
    const double phi = (180.0 + 3.75) * radpdeg;
    const double theta = (90.0 - 13.46) * radpdeg;
    const double psi = (180.0 + 14.111534) * radpdeg;

    const double sintsq = sqr(sin(theta));  /* sin^2(theta), cos^2(theta) */
    const double costsq = sqr(cos(theta));

    const double cosphisq = sqr(cos(phi));  /* sin(phi), cos(phi) */
    const double sinphisq = sqr(sin(phi));

    const double cospsisq = sqr(cos(psi));  /* sin^2(psi), cos^2(psi) */
    const double sinpsisq = sqr(sin(psi));

    const double sint = sin(theta);  /* sin(theta), cos(theta) */
    const double cost = cos(theta);

    const double sinphi = sin(phi);  /* sin(phi), cos(phi) */
    const double cosphi = cos(phi);

    const double sinpsi = sin(psi);  /* sin(psi), cos(psi) */
    const double cospsi = cos(psi);

    sgc->rot11 = -(  cost * sinphi * sinpsi
                   - costsq * cosphi * cospsi
                   - cospsi * sintsq * cosphi)
                            /
                  (  cospsisq * cosphisq * costsq
                   + cospsisq * cosphisq * sintsq
                   + costsq * sinpsisq * sinphisq
                   + sinpsisq * cosphisq * costsq
                   + sinpsisq * cosphisq * sintsq
                   + costsq * cospsisq * sinphisq
                   + sintsq * sinphisq * cospsisq
                   + sintsq * sinphisq * sinpsisq);

    sgc->rot12 = -(  cost * sinphi * cospsi
                   + costsq * cosphi * sinpsi
                   + sinpsi * sintsq * cosphi)
                           /
                 (  cospsisq * cosphisq * costsq
                   + cospsisq * cosphisq * sintsq
                   + costsq * sinpsisq * sinphisq
                   + sinpsisq * cosphisq * costsq
                   + sinpsisq * cosphisq * sintsq
                   + costsq * cospsisq * sinphisq
                   + sintsq * sinphisq * cospsisq
                   + sintsq * sinphisq * sinpsisq);

    sgc->rot13 = (sint * sinphi)
                      /
        (costsq * sinphisq + cosphisq * costsq + cosphisq * sintsq + sintsq * sinphisq);

    sgc->rot21 = (  cost * cosphi * sinpsi
                  + costsq * cospsi * sinphi
                  + cospsi * sintsq * sinphi)
                         /
                (  cospsisq * cosphisq * costsq
                  + cospsisq * cosphisq * sintsq
                  + costsq * sinpsisq * sinphisq
                  + sinpsisq * cosphisq * costsq
                  + sinpsisq * cosphisq * sintsq
                  + costsq * cospsisq * sinphisq
                  + sintsq * sinphisq * cospsisq
                  + sintsq * sinphisq * sinpsisq);

    sgc->rot22 = -(  -cost * cosphi * cospsi
                   + costsq * sinpsi * sinphi
                   + sinpsi * sintsq * sinphi)
                             /
                  (  cospsisq * cosphisq * costsq
                   + cospsisq * cosphisq * sintsq
                   + costsq * sinpsisq * sinphisq
                   + sinpsisq * cosphisq * costsq
                   + sinpsisq * cosphisq * sintsq
                   + costsq * cospsisq * sinphisq
                   + sintsq * sinphisq * cospsisq
                   + sintsq * sinphisq * sinpsisq);

    sgc->rot23 = -(sint * cosphi)
                        /
                  (  costsq * sinphisq
                   + cosphisq * costsq
                   + cosphisq * sintsq
                   + sintsq * sinphisq);

    sgc->rot31 = (sinpsi * sint)
                        /
                 (  cospsisq * costsq
                  + sinpsisq * sintsq
                  + cospsisq * sintsq
                  + sinpsisq * costsq);

    sgc->rot32 = (cospsi * sint)
                        /
                 (  cospsisq * costsq
                  + sinpsisq * sintsq
                  + cospsisq * sintsq
                  + sinpsisq * costsq);

    sgc->rot33 = cost / (costsq + sintsq);
}

//mathematic reversal of majewski's defined rotations for lbr->sgr conversion
void sgrToGal(double lamda, double beta, double* l, double* b)
{
    static const double radpdeg = M_PI / 180.0;
    double x2 = 0.0, y2 = 0.0;

    SGR_TO_GAL_CONSTANTS _sgc;  /* FIXME: Shouldn't be done each call */
    init_sgr_to_gal_constants(&_sgc);
    SGR_TO_GAL_CONSTANTS* sgc = &_sgc;

    if (beta > 90)
    {
        beta = 90 - (beta - 90);
        lamda = lamda + 180;
        if (lamda > 360)
        {
            lamda = lamda - 360;
        }
    }
    if (beta < -90)
    {
        beta = -90 - (beta + 90);
        lamda = lamda + 180;
        if (lamda > 360)
        {
            lamda = lamda - 360;
        }
    }
    if (lamda < 0)
    {
        lamda = lamda + 360;
    }

    beta = beta + 90;

    beta = beta * radpdeg;
    double z2 = cos(beta);

    if (lamda == 0)
    {
        lamda = lamda * radpdeg;
        x2 = sin(beta);
        y2 = 0;

    }
    else if (lamda < 90)
    {
        lamda = lamda * radpdeg;
        x2 = sqrt((1 - cos(beta) * cos(beta)) / (1 + tan(lamda) * tan(lamda)));
        y2 = x2 * tan(lamda);

    }
    else if (lamda == 90)
    {
        lamda = lamda * radpdeg;
        x2 = 0;
        y2 = sin(beta);

    }
    else  if (lamda < 180)
    {
        lamda = lamda * radpdeg;
        y2 = sqrt((1 - cos(beta) * cos(beta)) / (1 / (tan(lamda) * tan(lamda)) + 1));
        x2 = y2 / tan(lamda);

    }
    else if (lamda == 180)
    {
        lamda = lamda * radpdeg;
        x2 = -sin(beta);
        y2 = 0;

    }
    else if (lamda < 270)
    {
        lamda = lamda * radpdeg;
        x2 = sqrt((1 - cos(beta) * cos(beta)) / (1 + tan(lamda) * tan(lamda)));
        y2 = x2 * tan(lamda);
        x2 = -x2;
        y2 = -y2;

    }
    else if (lamda == 270)
    {
        lamda = lamda * radpdeg;
        x2 = 0;
        y2 = -sin(beta);

    }
    else if (lamda < 360)
    {
        lamda = lamda * radpdeg;
        x2 = sqrt((1.0 - cos(beta) * cos(beta)) / (1 + tan(lamda) * tan(lamda)));
        y2 = x2 * tan(lamda);

    }
    else if (lamda == 360)
    {
        lamda = lamda * radpdeg;
        x2 = sin(beta);
        y2 = 0;
    }

    double x1 = sgc->rot11 * x2 + sgc->rot12 * y2 + sgc->rot13 * z2;
    double y1 = sgc->rot21 * x2 + sgc->rot22 * y2 + sgc->rot23 * z2;
    double z1 = sgc->rot31 * x2 + sgc->rot32 * y2 + sgc->rot33 * z2;

    if (z1 > 1)
    {
        *l = 0;
        *b = 90;
    }
    else
    {
        *b = asin(z1);
        *b = *b / radpdeg;
        *l = atan2(y1, x1);
        *l = *l / radpdeg;
        if (*l < 0)
        {
            *l = *l + 360.0;
        }
    }

    return;
}

/* Return ra & dec from survey longitude and latitude */
void atSurveyToEq (double slong, double slat, double* ra, double* dec)
{
    double anode, etaPole;
    double x1, y1, z1;

    double surveyCenterRa = at_surveyCenterRa;
    double surveyCenterDec = at_surveyCenterDec;

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

/* Return eta from stripe number */
double atEtaFromStripeNumber(int wedge)
{
    double eta;
    double stripeSeparation = at_stripeSeparation;

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

