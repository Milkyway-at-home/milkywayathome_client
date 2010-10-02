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
            st_probs[i] += calc_st_prob_inc(&sc[i], xyz, qw_r3_N);
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

ALWAYS_INLINE
inline NU_ID calc_nu_step(__MW_CONSTANT INTEGRAL_AREA* ia, const unsigned int nu_step)
{
    NU_ID nuid;
    real tmp1, tmp2;

    nuid.nu = ia->nu_min + (nu_step * ia->nu_step_size);

    tmp1 = d2r(90.0 - nuid.nu - ia->nu_step_size);
    tmp2 = d2r(90.0 - nuid.nu);

    nuid.id = mw_cos(tmp1) - mw_cos(tmp2);
    nuid.nu += 0.5 * ia->nu_step_size;

    return nuid;
}

#ifdef __cplusplus
}
#endif

#endif /* _INTEGRALS_COMMON_H_ */

