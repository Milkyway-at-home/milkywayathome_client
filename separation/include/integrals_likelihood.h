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

#ifndef _INTEGRALS_LIKELIHOOD_H_
#define _INTEGRALS_LIKELIHOOD_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "separation_constants.h"
#include "separation_types.h"
#include "milkyway_math.h"


ALWAYS_INLINE HOT
inline mwvector lbr2xyz_2(const real r_point, const LB_TRIG lbt)
{
    mwvector xyz;
    real zp = r_point * lbt.bcos;

    xyz.x = zp * lbt.lcos - sun_r0;
    xyz.y = zp * lbt.lsin;
    xyz.z = r_point * lbt.bsin;
    return xyz;
}

ALWAYS_INLINE HOT
inline real calc_st_prob_inc(__MW_CONSTANT STREAM_CONSTANTS* sc, const mwvector xyz, const real qw_r3_N)
{
    mwvector xyzs, tmp;
    real xyz_norm, dotted;

    xyzs = mw_subv(xyz, sc->c);
    dotted = mw_dotv(sc->a, xyzs);
    tmp = mw_mulvs(dotted, sc->a);
    mw_incsubv(xyzs, tmp);

    xyz_norm = mw_sqrv(xyzs);

    return qw_r3_N * mw_exp(-xyz_norm * sc->sigma_sq2_inv);
}

ALWAYS_INLINE HOT
inline real aux_prob(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                     const real qw_r3_N,
                     const real r_in_mag,
                     const real r_in_mag2)
{
    return qw_r3_N * (ap->bg_a * r_in_mag2 + ap->bg_b * r_in_mag + ap->bg_c);
}

ALWAYS_INLINE HOT CONST_F
inline real rg_calc(const mwvector xyz, const real q_inv_sqr)
{
    return mw_sqrt(sqr(X(xyz)) + sqr(Y(xyz)) + sqr(Z(xyz)) * q_inv_sqr);
}

ALWAYS_INLINE HOT CONST_F
inline real h_prob_fast(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap, const real qw_r3_N, const real rg)
{
    const real rs = rg + ap->r0;
    return qw_r3_N / (rg * cube(rs));
}

ALWAYS_INLINE HOT CONST_F
inline real h_prob_slow(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap, const real qw_r3_N, const real rg)
{
    const real rs = rg + ap->r0;
    return qw_r3_N / (mw_powr(rg, ap->alpha) * mw_powr(rs, ap->alpha_delta3));
}

ALWAYS_INLINE HOT CONST_F
inline LB_TRIG lb_trig(LB lb)
{
    LB_TRIG lbt;
    mw_sincos(d2r(LB_L(lb)), &lbt.lsin, &lbt.lcos);
    mw_sincos(d2r(LB_B(lb)), &lbt.bsin, &lbt.bcos);
    return lbt;
}


#ifdef __cplusplus
}
#endif

#endif /* _INTEGRALS_LIKELIHOOD_H_ */

