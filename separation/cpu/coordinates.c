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

#include "separation.h"
#include "coordinates.h"

/* Convert sun-centered lbr (degrees) into galactic xyz coordinates. */
void lbr2xyz(const double* lbr, vector xyz)
{
    double zp, d;
/*
    TODO: Use radians to begin with
    const double bsin = sin(B(lbr));
    const double lsin = sin(L(lbr));
    const double bcos = cos(B(lbr));
    const double lcos = cos(L(lbr));
*/

    double lsin, lcos;
    double bsin, bcos;

    sincos(d2r(B(lbr)), &bsin, &bcos);
    sincos(d2r(L(lbr)), &lsin, &lcos);

    Z(xyz) = R(lbr) * bsin;
    zp = R(lbr) * bcos;
    d = sqrt( sqr(sun_r0) + sqr(zp) - 2.0 * sun_r0 * zp * lcos);
    X(xyz) = (sqr(zp) - sqr(sun_r0) - sqr(d)) / (2.0 * sun_r0);
    Y(xyz) = zp * lsin;
}


//vickej2 for sgr stripes, the great circles are defined thus:
//sgr stripes run parallel to sgr longitude lines, centered on lamda=2.5*wedge number
//with mu=0 at the sgr equator and increasing in the +z direction (increasing from the equator with beta)
//and nu=0 at the center and increasing in the -y direction (inversely to lamda)
//in this manner an equatorial stripe of standard coordinate conventions is created.




/* (ra, dec) in degrees */
/* Get eta for the given wedge. */
inline static double wedge_eta(int wedge)
{
    return wedge * d2r(stripeSeparation) - d2r(57.5) - (wedge > 46 ? M_PI : 0.0);
}

/* Get inclination for the given wedge. */
inline static double wedge_incl(int wedge)
{
    return wedge_eta(wedge) + d2r(surveyCenterDec);
}

#define anode d2r(NODE_GC_COORDS)

/* Convert GC coordinates (mu, nu) into l and b for the given wedge. */
void gc2lb(int wedge, double mu, double nu, double* restrict l, double* restrict b)
{

    mu = d2r(mu);
    nu = d2r(nu);

    /* Rotation */
    double sinnu, cosnu;
    sincos(nu, &sinnu, &cosnu);

    double sinmunode, cosmunode;
    double munode = mu - anode;
    sincos(munode, &sinmunode, &cosmunode);

    const double x12 = cosmunode * cosnu;  /* x1 = x2 */
    const double y2 = sinmunode * cosnu;
    /* z2 = sin(nu) */

    double sininc, cosinc;
    sincos(wedge_incl(wedge), &sininc, &cosinc);

    const double y1 = y2 * cosinc - sinnu * sininc;
    const double z1 = y2 * sininc + sinnu * cosinc;

    const double ra = atan2(y1, x12) + anode;
    const double dec = asin(z1);

    /* Use SLALIB to do the actual conversion */
    vector v2;

    {
        unsigned int i, j;

        static const double rmat[3][3] =
            {
                { -0.054875539726, -0.873437108010, -0.483834985808 },
                {  0.494109453312, -0.444829589425,  0.746982251810 },
                { -0.867666135858, -0.198076386122,  0.455983795705 }
            };

        /* Spherical to Cartesian */
        double sinra, cosra;
        sincos(ra, &sinra, &cosra);

        const double cosdec = cos(dec);

        const vector v1 = { cosra * cosdec,
                            sinra * cosdec,
                            z1               /* sin(asin(z1)) = z1 */
                          };

        /* Equatorial to Galactic */

        /* Matrix rmat * vector v1 -> vector vb */
        for ( i = 0; i < 3; ++i )
        {
            v2[i] = 0.0;
            for ( j = 0; j < 3; ++j )
            {
                v2[i] += rmat[i][j] * v1[j];
            }
        }
    }

    /* Cartesian to spherical */
    {
        double r = sqrt(sqr(X(v2)) + sqr(Y(v2)));

        *l = ( r != 0.0 ) ? atan2( Y(v2), X(v2) ) : 0.0;
        *b = ( Z(v2) != 0.0 ) ? atan2( Z(v2), r ) : 0.0;
    }


    *l = r2d(*l);
    *b = r2d(*b);
}

