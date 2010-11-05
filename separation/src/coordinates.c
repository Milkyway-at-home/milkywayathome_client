/*
Copyright 2008-2010 Travis Desell, Dave Przybylo, Nathan Cole, Matthew
Arsenault, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
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

#include "coordinates.h"

/* Convert GC coordinates (mu, nu) into l and b for the given wedge. */
HOT CONST_F
LB gc2lb(const int wedge, const real mu, const real nu)
{
    LB lb;
    real sinmunode, cosmunode;
    real sinnu, cosnu;
    real munode;
    real sininc, cosinc;
    real sinra, cosra;

    real x12, y2, y1, z1;
    real ra, dec;
    real wedge_eta, wedge_incl;
    real cosdec;
    mwvector v1, v2;
    real r;

    /* Use SLALIB to do the actual conversion */
    static const mwmatrix rmat =
        {
            mw_vec( -0.054875539726, -0.873437108010, -0.483834985808 ),
            mw_vec(  0.494109453312, -0.444829589425,  0.746982251810 ),
            mw_vec( -0.867666135858, -0.198076386122,  0.455983795705 )
        };

    /* Rotation */
    mw_sincos(d2r(nu), &sinnu, &cosnu);

    munode = mu - NODE_GC_COORDS;
    mw_sincos(d2r(munode), &sinmunode, &cosmunode);

    x12 = cosmunode * cosnu;  /* x1 = x2 */
    y2 = sinmunode * cosnu;
    /* z2 = sin(nu) */

    wedge_eta = atEtaFromStripeNumber_rad(wedge);

    /* Get inclination for the given wedge. */
    wedge_incl = wedge_eta + d2r(surveyCenterDec);

    mw_sincos(wedge_incl, &sininc, &cosinc);

    y1 = y2 * cosinc - sinnu * sininc;
    z1 = y2 * sininc + sinnu * cosinc;

    ra = mw_atan2(y1, x12) + NODE_GC_COORDS_RAD;
    dec = mw_asin(z1);


    /* Spherical to Cartesian */
    mw_sincos(ra, &sinra, &cosra);

    cosdec = mw_cos(dec);
    SET_VECTOR(v1,
               cosra * cosdec,
               sinra * cosdec,
               z1         /* mw_sin(asin(z1)) == z1 */
              );

    /* Equatorial to Galactic */

    /* Matrix rmat * vector v1 -> vector vb */
    v2 = mw_mulmv(rmat, v1);

    /* Cartesian to spherical */
    r = mw_hypot(X(v2), Y(v2));

    LB_L(lb) = ( r != 0.0 ) ? mw_atan2( Y(v2), X(v2) ) : 0.0;
    LB_B(lb) = ( Z(v2) != 0.0 ) ? mw_atan2( Z(v2), r ) : 0.0;

    LB_L(lb) = r2d(LB_L(lb));
    LB_B(lb) = r2d(LB_B(lb));

    return lb;
}

