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

#define lbr2xyz_2(xyz, r_point, lbt)            \
    {                                           \
        real zp = r_point * lbt.bcos;           \
        X(xyz) = zp * lbt.lcos - sun_r0;        \
        Y(xyz) = zp * lbt.lsin;                 \
        Z(xyz) = r_point * lbt.bsin;            \
    }

ALWAYS_INLINE HOT
inline real calc_st_prob_inc(__MW_CONSTANT STREAM_CONSTANTS* sc, const vector xyz, const real qw_r3_N)
{
    vector xyzs;
    real xyz_norm, dotted;

    SUBV(xyzs, xyz, sc->c);
    DOTVP(dotted, sc->a, xyzs);
    INCSUBVMS(xyzs, dotted, sc->a);
    SQRV(xyz_norm, xyzs);

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
inline real rg_calc(const vector xyz, const real q_inv_sqr)
{
    return mw_sqrt(sqr(X(xyz)) + sqr(Y(xyz)) + sqr(Z(xyz)) * q_inv_sqr);
}

ALWAYS_INLINE HOT CONST_F
inline real h_prob_fast(const real qw_r3_N, const real rg, const real rs)
{
    return qw_r3_N / (rg * cube(rs));
}

ALWAYS_INLINE HOT CONST_F
inline real h_prob_slow(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap, const real qw_r3_N, const real rg)
{
    return qw_r3_N / (mw_powr(rg, ap->alpha) * mw_powr(rg + ap->r0, ap->alpha_delta3));
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

