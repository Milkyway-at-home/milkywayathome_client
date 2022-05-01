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

mwvector cartesianToLbr_rad(mwvector* r, real_0 sunGCDist)
{
    mwvector lbR;

    const real xp = mw_add_s(&X(r), sunGCDist);
    const real rho = mw_hypot(&xp, &Y(r));

    lbR.x = mw_atan2(&Y(r), &xp);
    lbR.y = mw_atan2(&Z(r), &rho );
    lbR.z = mw_hypot(&rho, &Z(r));
    lbR.w = ZERO_REAL;

    if (showRealValue(&lbR.x) < 0.0)
        lbR.x = mw_add_s(&(lbR.x), M_2PI);

    return lbR;
}

mwvector cartesianToLbr(mwvector* r, real_0 sunGCDist)
{
    mwvector lbR;
    lbR = cartesianToLbr_rad(r, sunGCDist);
    lbR.x = r2d(&(lbR.x));
    lbR.y = r2d(&(lbR.y));

    return lbR;
}

static inline mwvector _lbrToCartesian(const real* l, const real* b, const real* r, const real_0 sun)
{
    mwvector cart;
    real tmp;
    real cosl = mw_cos(l);
    real cosb = mw_cos(b);
    real sinl = mw_sin(l);
    real sinb = mw_sin(b);

    tmp = mw_mul(&cosl, &cosb);
    tmp = mw_mul(r, &tmp);
    cart.x = mw_add_s(&tmp, -sun);

    tmp = mw_mul(&sinl, &cosb);
    cart.y = mw_mul(r, &tmp);
    cart.z = mw_mul(r, &sinb);
    cart.w = ZERO_REAL;

    return cart;
}

mwvector lbrToCartesian_rad(mwvector* lbr, real_0 sunGCDist)
{
    return _lbrToCartesian(&L(lbr), &B(lbr), &R(lbr), sunGCDist);
}

mwvector lbrToCartesian(mwvector* lbr, real_0 sunGCDist)
{
    real Lrad = d2r(&L(lbr));
    real Brad = d2r(&B(lbr));
    return _lbrToCartesian(&Lrad, &Brad, &R(lbr), sunGCDist);
}

void nbGetHistTrig(NBHistTrig* ht, const HistogramParams* hp, mwbool leftHanded)
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

    ht->leftHanded = leftHanded;
}

real nbXYZToLambda(const NBHistTrig* ht, mwvector* xyz, real_0 sunGCDist)
{
    real bcos, bsin, lsin, lcos, tmp;
    real lambda;
    mwvector lbr;
    real part1, part2, part3;

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
    bcos = mw_cos(&B(&lbr));
    bsin = mw_sin(&B(&lbr));
    lsin = mw_sin(&L(&lbr));
    lcos = mw_cos(&L(&lbr));

    part1 = mw_mul(&bcos, &lcos);
    part1 = mw_mul_s(&part1, - (sinpsi * cosphi + costh * sinphi * cospsi));
    part2 = mw_mul(&bcos, &lsin);
    part2 = mw_mul_s(&part2, (-sinpsi * sinphi + costh * cosphi * cospsi));
    part3 = mw_mul_s(&bsin, cospsi * sinth);
    real y_like = mw_add(&part1, &part2);
    y_like = mw_add(&y_like, &part3);

    part1 = mw_mul(&bcos, &lcos);
    part1 = mw_mul_s(&part1, cospsi * cosphi - costh * sinphi * sinpsi);
    part2 = mw_mul(&bcos, &lsin);
    part2 = mw_mul_s(&part2, cospsi * sinphi + costh * cosphi * sinpsi);
    part3 = mw_mul_s(&bsin, sinpsi * sinth);
    real x_like = mw_add(&part1, &part2);
    x_like = mw_add(&y_like, &part3);

    tmp = mw_atan2(&y_like, &x_like);
    lambda = r2d(&tmp);

    return lambda;
}

/* Transform positions from standard left handed Galactocentric XYZ to
/ the heliocentric Sgr system (lambda=0 at Sgr)
/ Input must be in kpc of the form X Y Z
/ Output is in kpc and degrees, of the form X_Sgr Y_Sgr Z_Sgr r lambda beta 
Adapted from http://www.astro.virginia.edu/~srm4n/Sgr/SgrCoord.h*/
/* Still Needs to be tested!!! */
mwvector nbXYZToLambdaBeta(const NBHistTrig* ht, mwvector* xyz, real_0 sunGCDist)  
{
    mwvector lambdabetar;
    real tempX, tempY, tempZ, tmp;
    real tempX1, tempX2, tempX3, tempY1, tempY2, tempY3, tempZ1,  tempZ2,  tempZ3; 
    real_0 cosphi = ht->cosphi;
    real_0 sinphi = ht->sinphi;
    real_0 sinpsi = ht->sinpsi;
    real_0 cospsi = ht->cospsi;
    real_0 costh = ht->costh;
    real_0 sinth = ht->sinth;

    mwbool LH = ht->leftHanded;

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

    real shiftX = mw_add_s(&X(xyz), sunGCDist);

    /* Calculate X,Y,Z,distance in the Sgr system */
    tempX1 = mw_mul_s(&shiftX, rot11);
    tempX2 = mw_mul_s(&Y(xyz), rot12);
    tempX3 = mw_mul_s(&Z(xyz), rot13);
    tempX  = mw_add(&tempX1, &tempX2);
    tempX  = mw_add(&tempX, &tempX3);

    tempY1 = mw_mul_s(&shiftX, rot21);
    tempY2 = mw_mul_s(&Y(xyz), rot22);
    tempY3 = mw_mul_s(&Z(xyz), rot23);
    tempY  = mw_add(&tempY1, &tempY2);
    tempY  = mw_add(&tempY, &tempY3);

    tempZ1 = mw_mul_s(&shiftX, rot31);
    tempZ2 = mw_mul_s(&Y(xyz), rot32);
    tempZ3 = mw_mul_s(&Z(xyz), rot33);
    tempZ = mw_add(&tempZ1, &tempZ2);
    tempZ = mw_add(&tempZ, &tempZ3);

    R(&lambdabetar) = mw_hypot(&tempX, &tempY);
    R(&lambdabetar) = mw_hypot(&R(&lambdabetar), &tempZ);

    if(LH) tempZ = mw_neg(&tempZ); //Flips Beta value if in left-handed system
    /* Calculate the angular coordinates lambda,beta */
    tmp = mw_atan2(&tempY, &tempX);
    L(&lambdabetar) = r2d(&tmp);
    tmp = mw_div(&tempZ, &R(&lambdabetar));
    tmp = mw_asin(&tmp);
    B(&lambdabetar) = r2d(&tmp);
    return lambdabetar;
}
