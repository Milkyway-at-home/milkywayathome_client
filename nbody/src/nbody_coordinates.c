/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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

#include "nbody_config.h"
#include "nbody_coordinates.h"

mwvector cartesianToLbr_rad(mwvector r, real sunGCDist)
{
    mwvector lbR;

    const real xp = X(r) + sunGCDist;

    L(lbR) = mw_atan2(Y(r), xp);
    B(lbR) = mw_atan2( Z(r), mw_sqrt( sqr(xp) + sqr(Y(r)) ) );
    R(lbR) = mw_sqrt(sqr(xp) + sqr(Y(r)) + sqr(Z(r)));
    W(lbR) = 0.0;

    if (L(lbR) < 0.0)
        L(lbR) += M_2PI;

    return lbR;
}

mwvector cartesianToLbr(mwvector r, real sunGCDist)
{
    mwvector lbR;
    lbR = cartesianToLbr_rad(r, sunGCDist);
    L(lbR) = r2d(L(lbR));
    B(lbR) = r2d(B(lbR));

    return lbR;
}

static inline mwvector _lbrToCartesian(const real l, const real b, const real r, const real sun)
{
    mwvector cart;

    X(cart) = r * mw_cos(l) * mw_cos(b) - sun;
    Y(cart) = r * mw_sin(l) * mw_cos(b);
    Z(cart) = r * mw_sin(b);
    W(cart) = 0.0;

    return cart;
}

mwvector lbrToCartesian_rad(mwvector lbr, real sunGCDist)
{
    return _lbrToCartesian(L(lbr), B(lbr), R(lbr), sunGCDist);
}

mwvector lbrToCartesian(mwvector lbr, real sunGCDist)
{
    return _lbrToCartesian(d2r(L(lbr)), d2r(B(lbr)), R(lbr), sunGCDist);
}

void nbGetHistTrig(NBHistTrig* ht, const HistogramParams* hp)
{
    real rphi = d2r(hp->phi);
    real rpsi = d2r(hp->psi);
    real rth  = d2r(hp->theta);

    ht->cosphi = mw_cos(rphi);
    ht->sinphi = mw_sin(rphi);
    ht->sinpsi = mw_sin(rpsi);
    ht->cospsi = mw_cos(rpsi);
    ht->costh  = mw_cos(rth);
    ht->sinth  = mw_sin(rth);
}

real nbXYZToLambda(const NBHistTrig* ht, mwvector xyz, real sunGCDist)
{
    real bcos, bsin, lsin, lcos;
    real lambda;
    mwvector lbr;

    real cosphi = ht->cosphi;
    real sinphi = ht->sinphi;
    real sinpsi = ht->sinpsi;
    real cospsi = ht->cospsi;
    real costh = ht->costh;
    real sinth = ht->sinth;

    /* Convert to (l,b) (involves convert x to Sun-centered)
       Leave in radians to make rotation easier */
    lbr = cartesianToLbr_rad(xyz, sunGCDist);

    /* Convert to (lambda, beta) (involves a rotation using the
       Newberg et al (2009) rotation matrices) which gets it from http://www.astro.virginia.edu/~srm4n/Sgr/SgrCoord.h**/
    bcos = mw_cos(B(lbr));
    bsin = mw_sin(B(lbr));
    lsin = mw_sin(L(lbr));
    lcos = mw_cos(L(lbr));

    lambda = r2d(mw_atan2(
                      -(sinpsi * cosphi + costh * sinphi * cospsi) * bcos * lcos
                     + (-sinpsi * sinphi + costh * cosphi * cospsi) * bcos * lsin
                     + cospsi * sinth * bsin,

                     (cospsi * cosphi - costh * sinphi * sinpsi) * bcos * lcos
                     + (cospsi * sinphi + costh * cosphi * sinpsi) * bcos * lsin
                     + sinpsi * sinth * bsin ));

    return lambda;
}

/* Transform positions from standard left handed Galactocentric XYZ to
/ the heliocentric Sgr system (lambda=0 at Sgr)
/ Input must be in kpc of the form X Y Z
/ Output is in kpc and degrees, of the form X_Sgr Y_Sgr Z_Sgr r lambda beta 
Adapted from http://www.astro.virginia.edu/~srm4n/Sgr/SgrCoord.h*/
/* Still Needs to be tested!!! */
mwvector nbXYZToLambdaBeta(const NBHistTrig* ht, mwvector xyz, real sunGCDist)  
{
    mwvector lambdabetar;
    real tempX, tempY, tempZ;
    real cosphi = ht->cosphi;
    real sinphi = ht->sinphi;
    real sinpsi = ht->sinpsi;
    real cospsi = ht->cospsi;
    real costh = ht->costh;
    real sinth = ht->sinth;

    /* Define the rotation matrix from the Euler angles */
    real rot11 = cospsi * cosphi - costh * sinphi * sinpsi;
    real rot12 = cospsi * sinphi + costh * cosphi * sinpsi;
    real rot13 = sinpsi * sinth;
    real rot21 = -sinpsi * cosphi - costh * sinphi * cospsi;
    real rot22 = -sinpsi * sinphi + costh * cosphi * cospsi;
    real rot23 = cospsi * sinth;
    real rot31 = sinth * sinphi;
    real rot32 = -sinth * cosphi;
    real rot33 = costh;

    X(xyz) = X(xyz) + sunGCDist;

    /* Calculate X,Y,Z,distance in the Sgr system */
    tempX = rot11 * X(xyz) + rot12 * Y(xyz) + rot13 * Z(xyz);
    tempY = rot21 * X(xyz) + rot22 * Y(xyz) + rot23 * Z(xyz);
    tempZ = rot31 * X(xyz) + rot32 * Y(xyz) + rot33 * Z(xyz);
    R(lambdabetar) = mw_sqrt(tempX * tempX + tempY * tempY + tempZ * tempZ);

    tempZ=-tempZ;
    /* Calculate the angular coordinates lambda,beta */
    L(lambdabetar) = r2d(mw_atan2(tempY,tempX));
    B(lambdabetar) = r2d(mw_asin(tempZ / R(lambdabetar)));
    return lambdabetar;
}
