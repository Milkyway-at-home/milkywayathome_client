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

mwvector cartesianToLbr_rad(mwvector r, real_0 sunGCDist)
{
    mwvector lbR;

    const real xp = mw_add(X(r), mw_real_const(sunGCDist));

    L(lbR) = mw_atan2(Y(r), xp);
    B(lbR) = mw_atan2( Z(r), mw_hypot(xp, Y(r)) );
    R(lbR) = mw_hypot(mw_hypot(xp, Y(r)), Z(r));
    W(lbR) = ZERO_REAL;

    if (showRealValue(L(lbR)) < 0.0)
        L(lbR) = mw_add(L(lbR), mw_real_const(M_2PI));

    return lbR;
}

mwvector cartesianToLbr(mwvector r, real_0 sunGCDist)
{
    mwvector lbR;
    lbR = cartesianToLbr_rad(r, sunGCDist);
    L(lbR) = r2d(L(lbR));
    B(lbR) = r2d(B(lbR));

    return lbR;
}

static inline mwvector _lbrToCartesian(const real l, const real b, const real r, const real_0 sun)
{
    mwvector cart;

    X(cart) = mw_sub(mw_mul(r, mw_mul(mw_cos(l), mw_cos(b))), mw_real_const(sun));
    Y(cart) = mw_mul(r, mw_mul(mw_sin(l), mw_cos(b)));
    Z(cart) = mw_mul(r, mw_sin(b));
    W(cart) = ZERO_REAL;

    return cart;
}

mwvector lbrToCartesian_rad(mwvector lbr, real_0 sunGCDist)
{
    return _lbrToCartesian(L(lbr), B(lbr), R(lbr), sunGCDist);
}

mwvector lbrToCartesian(mwvector lbr, real_0 sunGCDist)
{
    return _lbrToCartesian(d2r(L(lbr)), d2r(B(lbr)), R(lbr), sunGCDist);
}

void nbGetHistTrig(NBHistTrig* ht, const HistogramParams* hp)
{
    real_0 rphi = d2r_0(hp->phi);
    real_0 rpsi = d2r_0(hp->psi);
    real_0 rth  = d2r_0(hp->theta);

    ht->cosphi = mw_cos_0(rphi);
    ht->sinphi = mw_sin_0(rphi);
    ht->sinpsi = mw_sin_0(rpsi);
    ht->cospsi = mw_cos_0(rpsi);
    ht->costh  = mw_cos_0(rth);
    ht->sinth  = mw_sin_0(rth);
}

real nbXYZToLambda(const NBHistTrig* ht, mwvector xyz, real_0 sunGCDist)
{
    real bcos, bsin, lsin, lcos;
    real lambda;
    mwvector lbr;

    real_0 cosphi = ht->cosphi;
    real_0 sinphi = ht->sinphi;
    real_0 sinpsi = ht->sinpsi;
    real_0 cospsi = ht->cospsi;
    real_0 costh = ht->costh;
    real_0 sinth = ht->sinth;

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
         mw_add(mw_add(mw_mul_s(mw_mul(bcos, lcos), - (sinpsi * cosphi + costh * sinphi * cospsi)),
                       mw_mul_s(mw_mul(bcos, lsin), (-sinpsi * sinphi + costh * cosphi * cospsi))),
                       mw_mul_s(bsin, cospsi * sinth)),

         mw_add(mw_add(mw_mul_s(mw_mul(bcos, lcos), cospsi * cosphi - costh * sinphi * sinpsi),
                       mw_mul_s(mw_mul(bcos, lsin), cospsi * sinphi + costh * cosphi * sinpsi)),
                       mw_mul_s(bsin, sinpsi * sinth)) ));

    return lambda;
}

/* Transform positions from standard left handed Galactocentric XYZ to
/ the heliocentric Sgr system (lambda=0 at Sgr)
/ Input must be in kpc of the form X Y Z
/ Output is in kpc and degrees, of the form X_Sgr Y_Sgr Z_Sgr r lambda beta 
Adapted from http://www.astro.virginia.edu/~srm4n/Sgr/SgrCoord.h*/
/* Still Needs to be tested!!! */
mwvector nbXYZToLambdaBeta(const NBHistTrig* ht, mwvector xyz, real_0 sunGCDist)  
{
    mwvector lambdabetar;
    real tempX, tempY, tempZ;
    real_0 cosphi = ht->cosphi;
    real_0 sinphi = ht->sinphi;
    real_0 sinpsi = ht->sinpsi;
    real_0 cospsi = ht->cospsi;
    real_0 costh = ht->costh;
    real_0 sinth = ht->sinth;

    /* Define the rotation matrix from the Euler angles */
    real_0 rot11 = cospsi * cosphi - costh * sinphi * sinpsi;
    real_0 rot12 = cospsi * sinphi + costh * cosphi * sinpsi;
    real_0 rot13 = sinpsi * sinth;
    real_0 rot21 = -sinpsi * cosphi - costh * sinphi * cospsi;
    real_0 rot22 = -sinpsi * sinphi + costh * cosphi * cospsi;
    real_0 rot23 = cospsi * sinth;
    real_0 rot31 = sinth * sinphi;
    real_0 rot32 = -sinth * cosphi;
    real_0 rot33 = costh;

    X(xyz) = mw_add(X(xyz), sunGCDist);

    /* Calculate X,Y,Z,distance in the Sgr system */
    tempX = mw_add(mw_add(mw_mul_s(X(xyz), rot11), mw_mul_s(Y(xyz), rot12)), mw_mul_s(Z(xyz), rot13));
    tempY = mw_add(mw_add(mw_mul_s(X(xyz), rot21), mw_mul_s(Y(xyz), rot22)), mw_mul_s(Z(xyz), rot23));
    tempZ = mw_add(mw_add(mw_mul_s(X(xyz), rot31), mw_mul_s(Y(xyz), rot32)), mw_mul_s(Z(xyz), rot33));
    R(lambdabetar) = mw_hypot(mw_hypot(tempX, tempY), tempZ);

    tempZ = mw_neg(tempZ);
    /* Calculate the angular coordinates lambda,beta */
    L(lambdabetar) = r2d(mw_atan2(tempY,tempX));
    B(lambdabetar) = r2d(mw_asin(mw_div(tempZ, R(lambdabetar))));
    return lambdabetar;
}
