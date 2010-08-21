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
#include "milkyway_math.h"
#include "coordinates.h"
#include "integrals_likelihood.h"

__attribute__ ((always_inline))
inline static void probabilities(__MW_CONSTANT const STREAM_CONSTANTS* sc,
                                 __MW_LOCAL const R_POINTS* r_pts,
                                 __MW_LOCAL vector* const xyz,
                                 const real V,
                                 const real reff_xr_rp3,
                                 const unsigned int number_streams,
                                 const unsigned int nconvolve,
                                 __MW_LOCAL ST_PROBS* probs)
{
    unsigned int i;
    real st_prob;

    for (i = 0; i < number_streams; ++i)
    {
        if (sc[i].large_sigma)
            st_prob = V * reff_xr_rp3 * probabilities_convolve(&sc[i], r_pts, xyz, nconvolve);
        else
            st_prob = 0.0;


        real _tmp = probs[i].st_prob_int;
        probs[i].st_prob_int += st_prob;
        probs[i].st_prob_int_c += (st_prob) - ((probs[i].st_prob_int) - _tmp);

        //KAHAN_ADD(probs[i].st_prob_int, st_prob, probs[i].st_prob_int_c);
    }
}

__attribute__ ((always_inline, const))
inline static real distance_magnitude(const real m)
{
    return mw_powr(10.0, (m - 14.2) / 5.0);
}

/* Sum over mu steps using Kahan summation */
__attribute__ ((always_inline, hot))
inline static BG_PROB mu_sum(__MW_PRIVATE const ASTRONOMY_PARAMETERS* ap,
                             __MW_CONSTANT const STREAM_CONSTANTS* sc,
                             __MW_LOCAL const R_POINTS* r_pts,
                             const real irv,             /* r constants */
                             const real reff_xr_rp3,
                             const real nu_consts_id,    /* nu constants */
                             const real nu_consts_nu,
                             const unsigned int mu_steps,
                             const real mu_step_size,
                             const real mu_min,
                             __MW_LOCAL ST_PROBS* probs,
                             __MW_LOCAL vector* xyz)
{
    unsigned int mu_step_current;
    real mu, V;
    real bg_prob;
    BG_PROB bg_prob_int = ZERO_BG_PROB; /* for Kahan summation */
    LB lb;

    for (mu_step_current = 0; mu_step_current < mu_steps; ++mu_step_current)
    {
        mu = mu_min + (((real) mu_step_current + 0.5) * mu_step_size);

        lb = gc2lb(ap->wedge, mu, nu_consts_nu);

        bg_prob = bg_probability(ap, r_pts, xyz, lb, reff_xr_rp3);

        V = irv * nu_consts_id;
        bg_prob *= V;

        KAHAN_ADD(bg_prob_int.bg_int, bg_prob, bg_prob_int.correction);

        probabilities(sc, r_pts, xyz, V, reff_xr_rp3, ap->number_streams, ap->convolve, probs);
    }

    return bg_prob_int;
}

#ifdef __cplusplus
}
#endif

#endif /* _INTEGRALS_COMMON_H_ */

