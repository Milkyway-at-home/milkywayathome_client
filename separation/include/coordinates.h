/*
Copyright 2008-2010 Travis Desell, Matthew Arsenault, Dave Przybylo,
Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

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

#ifndef _COORDINATES_H_
#define _COORDINATES_H_

#include "separation_constants.h"
#include "milkyway_vectors.h"
#include "milkyway_util.h"


/* Convert GC coordinates (mu, nu) into l and b for the given wedge. */
__attribute__ ((always_inline))
inline void gc2lb(const int wedge, double mu, double nu, double* restrict l, double* restrict b)
{
    mu = d2r(mu);
    nu = d2r(nu);

    /* Rotation */
    double sinnu, cosnu;
    sincos(nu, &sinnu, &cosnu);

    double sinmunode, cosmunode;
    double munode = mu - NODE_GC_COORDS_RAD;
    sincos(munode, &sinmunode, &cosmunode);

    const double x12 = cosmunode * cosnu;  /* x1 = x2 */
    const double y2 = sinmunode * cosnu;
    /* z2 = sin(nu) */

    const double wedge_eta = wedge * d2r(stripeSeparation) - d2r(57.5) - (wedge > 46 ? M_PI : 0.0);

    /* Get inclination for the given wedge. */
    const double wedge_incl = wedge_eta + d2r(surveyCenterDec);

    double sininc, cosinc;
    sincos(wedge_incl, &sininc, &cosinc);

    const double y1 = y2 * cosinc - sinnu * sininc;
    const double z1 = y2 * sininc + sinnu * cosinc;

    const double ra = atan2(y1, x12) + NODE_GC_COORDS_RAD;
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

        const vector v1 = VECTOR( cosra * cosdec,
                                  sinra * cosdec,
                                  z1         /* sin(asin(z1)) == z1 */
                                );

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


void lbr2xyz( const double* lbr, double* xyz );


#endif /* _COORDINATES_H_ */

