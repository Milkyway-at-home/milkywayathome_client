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
#include "separation_types.h"
#include "milkyway_cl.h"
#include "milkyway_extra.h"

//vickej2 for sgr stripes, the great circles are defined thus:
//sgr stripes run parallel to sgr longitude lines, centered on lamda=2.5*wedge number
//with mu=0 at the sgr equator and increasing in the +z direction (increasing from the equator with beta)
//and nu=0 at the center and increasing in the -y direction (inversely to lamda)
//in this manner an equatorial stripe of standard coordinate conventions is created.

/* Convert GC coordinates (mu, nu) into l and b for the given wedge. */
__attribute__ ((always_inline, hot, const))
inline LB gc2lb(const int wedge, const real mu, const real nu)
{
    LB lb;

    /* Rotation */
    real sinnu, cosnu;
    mw_sincos(d2r(nu), &sinnu, &cosnu);

    real sinmunode, cosmunode;
    real munode = mu - NODE_GC_COORDS;
    mw_sincos(d2r(munode), &sinmunode, &cosmunode);

    const real x12 = cosmunode * cosnu;  /* x1 = x2 */
    const real y2 = sinmunode * cosnu;
    /* z2 = sin(nu) */

    const real wedge_eta = wedge * d2r(stripeSeparation) - d2r((real) 57.5) - (wedge > 46 ? M_PI : 0.0);

    /* Get inclination for the given wedge. */
    const real wedge_incl = wedge_eta + d2r(surveyCenterDec);

    real sininc, cosinc;
    mw_sincos(wedge_incl, &sininc, &cosinc);

    const real y1 = y2 * cosinc - sinnu * sininc;
    const real z1 = y2 * sininc + sinnu * cosinc;

    const real ra = mw_atan2(y1, x12) + NODE_GC_COORDS_RAD;
    const real dec = mw_asin(z1);

    /* Use SLALIB to do the actual conversion */
    _MW_STATIC const matrix rmat =
        {
            VECTOR( -0.054875539726, -0.873437108010, -0.483834985808 ),
            VECTOR(  0.494109453312, -0.444829589425,  0.746982251810 ),
            VECTOR( -0.867666135858, -0.198076386122,  0.455983795705 )
        };

    /* Spherical to Cartesian */
    real sinra, cosra;
    mw_sincos(ra, &sinra, &cosra);

    const real cosdec = mw_cos(dec);
    vector v1 = VECTOR( cosra * cosdec,
                        sinra * cosdec,
                        z1         /* mw_sin(asin(z1)) == z1 */
                      );

    /* Equatorial to Galactic */

    /* Matrix rmat * vector v1 -> vector vb */
    vector v2;
    MULMV(v2, rmat, v1);

    /* Cartesian to spherical */
    const real r = mw_hypot(X(v2), Y(v2));

    LB_L(lb) = ( r != 0.0 ) ? mw_atan2( Y(v2), X(v2) ) : 0.0;
    LB_B(lb) = ( Z(v2) != 0.0 ) ? mw_atan2( Z(v2), r ) : 0.0;

    LB_L(lb) = r2d(LB_L(lb));
    LB_B(lb) = r2d(LB_B(lb));

    return lb;
}

#endif /* _COORDINATES_H_ */

