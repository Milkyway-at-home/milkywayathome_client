/*
Copyright 2008-2010 Travis Desell, Dave Przybylo, Nathan Cole, Matthew
Arsenault, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

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

#include "evaluation_state.h"
#include "integrals.h"
#include "integrals_likelihood.h"
#include "integrals_common.h"
#include "coordinates.h"
#include "r_points.h"
#include "milkyway_util.h"


#if BOINC_APPLICATION

ALWAYS_INLINE
static inline real progress(const EVALUATION_STATE* es,
                            const INTEGRAL_AREA* ia,
                            real total_calc_probs)
{
    /* This integral's progress */
    /* When checkpointing is done, ia->mu_step would always be 0 */
    unsigned int i_prog =  (es->r_step * ia->nu_steps * ia->mu_steps)
                         + (es->nu_step * ia->mu_steps); /* + es->mu_step */

    return (real)(i_prog + es->current_calc_probs) / total_calc_probs;
}

ALWAYS_INLINE
static inline void do_boinc_checkpoint(const EVALUATION_STATE* es,
                                       const INTEGRAL_AREA* ia,
                                       real total_calc_probs)
{
    if (boinc_time_to_checkpoint())
    {
        if (write_checkpoint(es))
            fail("Write checkpoint failed\n");
        boinc_checkpoint_completed();
    }

    boinc_fraction_done(progress(es, ia, total_calc_probs));
}

#else

#define do_boinc_checkpoint(es, ia, total_calc_probs)

#endif /* BOINC_APPLICATION */

ALWAYS_INLINE
static inline void nu_sum(const ASTRONOMY_PARAMETERS* ap,
                          const STREAM_CONSTANTS* sc,
                          const INTEGRAL_AREA* ia,
                          const real irv,
                          const real reff_xr_rp3,
                          const R_POINTS* r_pts,
                          const NU_CONSTANTS* nu_consts,
                          real* st_probs,
                          ST_PROBS* probs,
                          EVALUATION_STATE* es)
{
    BG_PROB mu_result;

    const unsigned int nu_steps = ia->nu_steps;
    const unsigned int mu_steps = ia->mu_steps;
    const real mu_min = ia->mu_min;
    const real mu_step_size = ia->mu_step_size;

    for ( ; es->nu_step < nu_steps; es->nu_step++)
    {
        do_boinc_checkpoint(es, ia, ap->total_calc_probs);

        mu_result = mu_sum(ap,
                           sc,
                           r_pts,
                           irv,
                           reff_xr_rp3,
                           nu_consts[es->nu_step].id,
                           nu_consts[es->nu_step].nu,
                           mu_steps,
                           mu_step_size,
                           mu_min,
                           st_probs,
                           probs);

        INCADD_BG_PROB(es->nu_acc, mu_result);
    }

    es->nu_step = 0;
}

static real r_sum(const ASTRONOMY_PARAMETERS* ap,
                  const INTEGRAL_AREA* ia,
                  const STREAM_CONSTANTS* sc,
                  const STREAM_GAUSS* sg,
                  const NU_CONSTANTS* nu_consts,
                  R_POINTS* r_pts,
                  real* st_probs,
                  ST_PROBS* probs,
                  EVALUATION_STATE* es)
{
    real reff_xr_rp3;
    R_PRIME rp;

    const unsigned int r_steps = ia->r_steps;

    for ( ; es->r_step < r_steps; es->r_step++)
    {
        rp = calcRPrime(ia, es->r_step);

        set_r_points(ap, sg, ap->convolve, rp.rPrime, r_pts);
        reff_xr_rp3 = calcReffXrRp3(rp.rPrime);

        nu_sum(ap, sc, ia, rp.irv, reff_xr_rp3, r_pts, nu_consts, st_probs, probs, es);

        INCADD_BG_PROB(es->r_acc, es->nu_acc);
        CLEAR_BG_PROB(es->nu_acc);
    }

    return es->r_acc.bg_int + es->r_acc.correction;
}

/* returns background integral */
real integrate(const ASTRONOMY_PARAMETERS* ap,
               const INTEGRAL_AREA* ia,
               const STREAM_CONSTANTS* sc,
               const STREAM_GAUSS* sg,
               ST_PROBS* probs,
               EVALUATION_STATE* es)
{
    real result;
    NU_CONSTANTS* nu_consts;
    R_POINTS* r_pts;
    real* st_probs;

    nu_consts = (NU_CONSTANTS*) prepare_nu_constants(ia->nu_steps, ia->nu_step_size, ia->nu_min);
    r_pts = (R_POINTS*) mallocSafe(sizeof(R_POINTS) * ap->convolve);
    st_probs = (real*) mallocSafe(sizeof(real) * ap->number_streams);

    result = r_sum(ap, ia, sc, sg, nu_consts, r_pts, st_probs, probs, es);
    es->r_step = 0;

    free(nu_consts);
    free(r_pts);
    free(st_probs);

    return result;
}


