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
#include "milkyway_vectors.h"
#include "separation_cl.h"

/* Used in innermost loops of integrals and likelihood calculation,
 * and we want it inlined. C99 inlining is annoying and sort of forces
 * you to have functions in headers which is unfortunate. */

/* FIXME: Better name? */
__attribute__ ((always_inline, pure, hot))
inline double probabilities_convolve(__MW_CONSTANT const STREAM_CONSTANTS* sc,
                                     __MW_LOCAL const R_POINTS* r_pts,
                                     __MW_LOCAL vector* const xyz,
                                     const unsigned int convolve)
{
    unsigned int i;
    double dotted, xyz_norm;
    vector xyzs;

    double st_prob = 0.0;

    for (i = 0; i < convolve; i++)
    {
        SUBV(xyzs, xyz[i], sc->c);
        DOTVP(dotted, sc->a, xyzs);
        INCSUBVMS(xyzs, dotted, sc->a);
        SQRV(xyz_norm, xyzs);

        st_prob += r_pts[i].qw_r3_N * exp(-xyz_norm / sc->sigma_sq2);
    }

    return st_prob;
}

/* FIXME: I don't know what these do enough to name it properly */
__attribute__ ((always_inline, hot))
inline double sub_bg_probability1(__MW_PRIVATE const ASTRONOMY_PARAMETERS* ap,
                                  __MW_LOCAL const R_POINTS* r_pts,
                                  __MW_LOCAL vector* const xyz,
                                  const LB integral_point,
                                  const int aux_bg_profile,
                                  const unsigned int convolve)
{
    unsigned int i;
    double h_prob, aux_prob;
    double rg, rs, zp;
    double lsin, lcos;
    double bsin, bcos;
    double bg_prob = 0.0;

    mw_sincos(d2r(LB_L(integral_point)), &lsin, &lcos);
    mw_sincos(d2r(LB_B(integral_point)), &bsin, &bcos);

    for (i = 0; i < convolve; ++i)
    {
        Z(xyz[i]) = r_pts[i].r_point * bsin;
        zp = r_pts[i].r_point * bcos;
        X(xyz[i]) = zp * lcos - sun_r0;
        Y(xyz[i]) = zp * lsin;

        rg = sqrt(sqr(X(xyz[i])) + sqr(Y(xyz[i])) + sqr(Z(xyz[i])) / sqr(ap->q));
        rs = rg + ap->r0;

        h_prob = r_pts[i].qw_r3_N / (rg * cube(rs));

        //the hernquist profile includes a quadratic term in g
        if (aux_bg_profile)
        {
            aux_prob = r_pts[i].qw_r3_N * (  ap->bg_a * r_pts[i].r_in_mag2
                                           + ap->bg_b * r_pts[i].r_in_mag
                                           + ap->bg_c );
            h_prob += aux_prob;
        }

        bg_prob += h_prob;
    }

    return bg_prob;
}

__attribute__ ((always_inline))
inline double sub_bg_probability2(__MW_PRIVATE const ASTRONOMY_PARAMETERS* ap,
                                  __MW_LOCAL const R_POINTS* r_pts,
                                  __MW_LOCAL vector* const xyz,
                                  const LB integral_point,
                                  const unsigned int convolve)
{
    unsigned int i;
    double rg, zp;
    double lsin, lcos;
    double bsin, bcos;
    double bg_prob = 0.0;

    mw_sincos(d2r(LB_L(integral_point)), &lsin, &lcos);
    mw_sincos(d2r(LB_B(integral_point)), &bsin, &bcos);

    for (i = 0; i < convolve; ++i)
    {
        Z(xyz[i]) = r_pts[i].r_point * bsin;
        zp = r_pts[i].r_point * bcos;
        X(xyz[i]) = zp * lcos - sun_r0;
        Y(xyz[i]) = zp * lsin;

        rg = sqrt(sqr(X(xyz[i])) + sqr(Y(xyz[i])) + sqr(Z(xyz[i])) / sqr(ap->q));

        bg_prob += r_pts[i].qw_r3_N / (pow(rg, ap->alpha) * pow(rg + ap->r0, ap->alpha_delta3));
    }

    return bg_prob;
}

__attribute__ ((always_inline, hot))
inline double bg_probability(__MW_PRIVATE const ASTRONOMY_PARAMETERS* ap,
                             __MW_LOCAL const R_POINTS* r_pts,
                             __MW_LOCAL vector* const xyz,
                             const LB integral_point,
                             const double reff_xr_rp3)
{
    double bg_prob;

    /* if q is 0, there is no probability */
    if (ap->q == 0)
        bg_prob = -1.0;
    else
    {
        if (ap->alpha == 1 && ap->delta == 1)
            bg_prob = sub_bg_probability1(ap, r_pts, xyz, integral_point, ap->aux_bg_profile, ap->convolve);
        else
            bg_prob = sub_bg_probability2(ap, r_pts, xyz, integral_point, ap->convolve);

        bg_prob *= reff_xr_rp3;
    }

    return bg_prob;
}

#ifdef __cplusplus
}
#endif

#endif /* _INTEGRALS_LIKELIHOOD_H_ */

