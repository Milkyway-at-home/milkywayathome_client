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
#include "separation_types.h"
#include "coordinates.h"
#include "integrals_likelihood.h"
#include "r_points.h"


ALWAYS_INLINE
inline void zero_st_probs(real* st_probs, const unsigned int nstream)
{
    unsigned int i;

    for (i = 0; i < nstream; ++i)
        st_probs[i] = 0.0;
}

ALWAYS_INLINE HOT
inline void stream_sums(real* st_probs,
                        __MW_CONSTANT STREAM_CONSTANTS* sc,
                        const vector xyz,
                        const real qw_r3_N,
                        const unsigned int nstreams)
{
    unsigned int i;
    real dotted, xyz_norm;
    vector xyzs;

    for (i = 0; i < nstreams; ++i)
    {
        if (sc[i].large_sigma)
        {
            SUBV(xyzs, xyz, sc[i].c);
            DOTVP(dotted, sc[i].a, xyzs);
            INCSUBVMS(xyzs, dotted, sc[i].a);
            SQRV(xyz_norm, xyzs);

            st_probs[i] += qw_r3_N * mw_exp(-xyz_norm / sc[i].sigma_sq2);
        }
    }
}

ALWAYS_INLINE
inline void sum_probs(ST_PROBS* probs,
                      const real* st_probs,
                      const real V_reff_xr_rp3,
                      const unsigned int nstream)
{
    unsigned int i;
    for (i = 0; i < nstream; ++i)
        KAHAN_ADD(probs[i].st_prob_int, V_reff_xr_rp3 * st_probs[i], probs[i].st_prob_int_c);
}

ALWAYS_INLINE HOT
inline real aux_prob(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                     const real qw_r3_N,
                     const real r_in_mag,
                     const real r_in_mag2)
{
    return qw_r3_N * (ap->bg_a * r_in_mag2 + ap->bg_b * r_in_mag + ap->bg_b);
}

ALWAYS_INLINE HOT CONST_F
inline real rg_calc(const vector xyz, const real q)
{
    return mw_sqrt(sqr(X(xyz)) + sqr(Y(xyz)) + sqr(Z(xyz)) / sqr(q));
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


/* FIXME: I don't know what these do enough to name it properly */
ALWAYS_INLINE HOT
inline real sub_bg_probability1(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                                __MW_CONSTANT STREAM_CONSTANTS* sc,
                                __MW_CONSTANT STREAM_GAUSS* sg,
                                const LB integral_point,
                                const real gPrime,
                                const int aux_bg_profile,
                                const unsigned int convolve,
                                const R_POINTS* r_pts,
                                real* st_probs)
{
    unsigned int i;
    real h_prob;
    real rg, rs;
    real lsin, lcos;
    real bsin, bcos;
    vector xyz;
    real bg_prob = 0.0;

    mw_sincos(d2r(LB_L(integral_point)), &lsin, &lcos);
    mw_sincos(d2r(LB_B(integral_point)), &bsin, &bcos);

    for (i = 0; i < convolve; ++i)
    {
        lbr2xyz_2(xyz, r_pts[i].r_point, bsin, bcos, lsin, lcos);

        rg = rg_calc(xyz, ap->q);
        rs = rg + ap->r0;

        h_prob = h_prob_fast(r_pts[i].qw_r3_N, rg, rs);
        /* the Hernquist profile includes a quadratic term in g */
        if (aux_bg_profile)
            h_prob += aux_prob(ap, r_pts[i].qw_r3_N, r_pts[i].r_in_mag, r_pts[i].r_in_mag2);
        bg_prob += h_prob;

        stream_sums(st_probs, sc, xyz, r_pts[i].qw_r3_N, ap->number_streams);
    }

    return bg_prob;
}

ALWAYS_INLINE
inline real sub_bg_probability2(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                                __MW_CONSTANT STREAM_CONSTANTS* sc,
                                __MW_CONSTANT STREAM_GAUSS* sg,
                                const LB integral_point,
                                const real gPrime,
                                const unsigned int convolve,
                                const R_POINTS* r_pts,
                                real* st_probs)
{
    unsigned int i;
    real rg;
    real lsin, lcos;
    real bsin, bcos;
    vector xyz;
    R_POINTS r_pt;
    real bg_prob = 0.0;

    mw_sincos(d2r(LB_L(integral_point)), &lsin, &lcos);
    mw_sincos(d2r(LB_B(integral_point)), &bsin, &bcos);

    for (i = 0; i < convolve; ++i)
    {
        lbr2xyz_2(xyz, r_pts[i].r_point, bsin, bcos, lsin, lcos);

        rg = rg_calc(xyz, ap->q);

        bg_prob += h_prob_slow(ap, r_pts[i].qw_r3_N, rg);
        stream_sums(st_probs, sc, xyz, r_pts[i].qw_r3_N, ap->number_streams);
    }

    return bg_prob;
}

ALWAYS_INLINE HOT
inline real bg_probability(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                           __MW_CONSTANT STREAM_CONSTANTS* sc,
                           __MW_CONSTANT STREAM_GAUSS* sg,
                           const LB integral_point,
                           const real gPrime,
                           const real reff_xr_rp3,
                           const real V,
                           const R_POINTS* r_pts,
                           real* st_probs,
                           ST_PROBS* probs)


{
    real bg_prob;

    /* if q is 0, there is no probability */
    if (ap->q == 0)
        return -1.0;

    zero_st_probs(st_probs, ap->number_streams);

    if (ap->fast_h_prob)
    {
        bg_prob = sub_bg_probability1(ap,
                                      sc,
                                      sg,
                                      integral_point,
                                      gPrime,
                                      ap->aux_bg_profile,
                                      ap->convolve,
                                      r_pts,
                                      st_probs);
    }
    else
    {
        bg_prob = sub_bg_probability2(ap,
                                      sc,
                                      sg,
                                      integral_point,
                                      gPrime,
                                      ap->convolve,
                                      r_pts,
                                      st_probs);
    }

    sum_probs(probs, st_probs, V * reff_xr_rp3, ap->number_streams);

    bg_prob *= reff_xr_rp3;

    return bg_prob;
}


#ifdef __cplusplus
}
#endif

#endif /* _INTEGRALS_COMMON_H_ */

