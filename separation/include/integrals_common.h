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

#ifndef _INTEGRALS_COMMON_H_
#define _INTEGRALS_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "milkyway_cl.h"
#include "milkyway_math.h"
#include "separation_types.h"
#include "separation_constants.h"
#include "coordinates.h"
#include "r_points.h"



ALWAYS_INLINE HOT OLD_GCC_EXTERNINLINE
inline mwvector lbr2xyz_2(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                          const real r_point,
                          const LB_TRIG lbt)
{
    mwvector xyz;
    real zp = r_point * lbt.bcos;

    // This mad for some reason increases GPR usage by 1 pushing into next level of unhappy
    xyz.x = mw_mad(zp, lbt.lcos, ap->m_sun_r0);
    //xyz.x = zp * lbt.lcos - ap->sun_r0;
    xyz.y = zp * lbt.lsin;
    xyz.z = r_point * lbt.bsin;
    return xyz;
}

ALWAYS_INLINE HOT OLD_GCC_EXTERNINLINE
inline real calc_st_prob_inc(__MW_CONSTANT STREAM_CONSTANTS* sc, mwvector xyz, const real qw_r3_N)
{
    real xyz_norm, dotted;

    mw_incsubv(xyz, sc->c);
    dotted = mw_dotv(sc->a, xyz);
    mw_incsubv_s(xyz, sc->a, dotted);

    xyz_norm = mw_sqrv(xyz);

    return qw_r3_N * mw_exp(-xyz_norm * sc->sigma_sq2_inv);
}

ALWAYS_INLINE HOT OLD_GCC_EXTERNINLINE
inline real aux_prob(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                     const real qw_r3_N,
                     const real r_in_mag)
{
    //return qw_r3_N * (ap->bg_a * sqr(r_in_mag) + ap->bg_b * r_in_mag + ap->bg_c);
    real tmp;

    tmp = mw_mad(ap->bg_b, r_in_mag, ap->bg_c); /* bg_b * r_in_mag + bg_c */
    tmp = mw_mad(ap->bg_a, sqr(r_in_mag), tmp); /* bg_a * r_in_mag2 + (bg_b * r_in_mag + bg_c)*/

    return qw_r3_N * tmp;
}

ALWAYS_INLINE HOT CONST_F OLD_GCC_EXTERNINLINE
inline real rg_calc(const mwvector xyz, const real q_inv_sqr)
{
    /* sqrt(x^2 + y^2 + q_inv_sqr * z^2) */

    real tmp;

    tmp = sqr(X(xyz));
    tmp = mw_mad(Y(xyz), Y(xyz), tmp);           /* x^2 + y^2 */
    tmp = mw_mad(q_inv_sqr, sqr(Z(xyz)), tmp);   /* (q_invsqr * z^2) + (x^2 + y^2) */

    return mw_sqrt(tmp);
}

ALWAYS_INLINE HOT CONST_F OLD_GCC_EXTERNINLINE
inline real h_prob_fast(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap, const real qw_r3_N, const real rg)
{
    const real rs = rg + ap->r0;
    return qw_r3_N / (rg * cube(rs));
}

ALWAYS_INLINE HOT CONST_F OLD_GCC_EXTERNINLINE
inline real h_prob_slow(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap, const real qw_r3_N, const real rg)
{
    const real rs = rg + ap->r0;
    return qw_r3_N / (mw_powr(rg, ap->alpha) * mw_powr(rs, ap->alpha_delta3));
}

ALWAYS_INLINE HOT CONST_F OLD_GCC_EXTERNINLINE
inline LB_TRIG lb_trig(LB lb)
{
    LB_TRIG lbt;
    mw_sincos(d2r(LB_L(lb)), &lbt.lsin, &lbt.lcos);
    mw_sincos(d2r(LB_B(lb)), &lbt.bsin, &lbt.bcos);
    return lbt;
}

ALWAYS_INLINE OLD_GCC_EXTERNINLINE
inline void zero_st_probs(real* st_probs, const unsigned int nstream)
{
    unsigned int i;

    for (i = 0; i < nstream; ++i)
        st_probs[i] = 0.0;
}

ALWAYS_INLINE HOT OLD_GCC_EXTERNINLINE
inline void stream_sums(real* st_probs,
                        __MW_CONSTANT STREAM_CONSTANTS* sc,
                        const mwvector xyz,
                        const real qw_r3_N,
                        const unsigned int nstreams)
{
    unsigned int i;

    for (i = 0; i < nstreams; ++i)
        st_probs[i] += calc_st_prob_inc(&sc[i], xyz, qw_r3_N);
}

ALWAYS_INLINE OLD_GCC_EXTERNINLINE
inline void sum_probs(KAHAN* probs,
                      const real* st_probs,
                      const real V_reff_xr_rp3,
                      const unsigned int nstream)
{
    unsigned int i;

    for (i = 0; i < nstream; ++i)
        KAHAN_ADD(probs[i], V_reff_xr_rp3 * st_probs[i]);
}

ALWAYS_INLINE OLD_GCC_EXTERNINLINE
inline void mult_probs(real* st_probs, const real V_reff_xr_rp3, const unsigned int n_stream)
{
    unsigned int i;

    for (i = 0; i < n_stream; ++i)
        st_probs[i] *= V_reff_xr_rp3;
}

#ifdef __cplusplus
}
#endif

#endif /* _INTEGRALS_COMMON_H_ */

