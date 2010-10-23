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
#include "integrals_common.h"
#include "coordinates.h"
#include "r_points.h"
#include "milkyway_util.h"
#include "calculated_constants.h"


#if BOINC_APPLICATION

ALWAYS_INLINE
static inline real progress(const EVALUATION_STATE* es,
                            const INTEGRAL_AREA* ia,
                            real total_calc_probs)
{
    /* This integral's progress */
    /* When checkpointing is done, ia->mu_step would always be 0 */
    unsigned int i_prog =  (es->nu_step * ia->mu_steps * ia->r_steps)
                         + (es->mu_step * ia->r_steps); /* + es->r_step */

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

HOT
real bg_probability_fast_hprob(const ASTRONOMY_PARAMETERS* ap,
                               const STREAM_CONSTANTS* sc,
                               const STREAM_GAUSS* sg,
                               const LB_TRIG lbt,
                               const real gPrime,
                               const int aux_bg_profile,
                               const unsigned int convolve,
                               const R_POINTS* r_pts,
                               real* st_probs)
{
    unsigned int i;
    real h_prob, g, rg;
    mwvector xyz;
    real bg_prob = 0.0;

    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(r_pts[i].r_point, lbt);

        rg = rg_calc(xyz, ap->q_inv_sqr);


        h_prob = h_prob_fast(ap, r_pts[i].qw_r3_N, rg);

        /* Add a quadratic term in g to the the Hernquist profile */
        if (aux_bg_profile)
        {
            g = gPrime + sg[i].dx;
            h_prob += aux_prob(ap, r_pts[i].qw_r3_N, g);
        }
        bg_prob += h_prob;

        stream_sums(st_probs, sc, xyz, r_pts[i].qw_r3_N, ap->number_streams);
    }

    return bg_prob;
}


HOT
real bg_probability_slow_hprob(const ASTRONOMY_PARAMETERS* ap,
                               const STREAM_CONSTANTS* sc,
                               const STREAM_GAUSS* sg,
                               const LB_TRIG lbt,
                               const real gPrime,
                               const int aux_bg_profile,
                               const unsigned int convolve,
                               const R_POINTS* r_pts,
                               real* st_probs)
{
    unsigned int i;
    real rg, g;
    mwvector xyz;
    real bg_prob = 0.0;

    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(r_pts[i].r_point, lbt);

        rg = rg_calc(xyz, ap->q_inv_sqr);

        bg_prob += h_prob_slow(ap, r_pts[i].qw_r3_N, rg);
        if (aux_bg_profile)
        {
            g = gPrime + sg[i].dx;
            bg_prob += aux_prob(ap, r_pts[i].qw_r3_N, g);
        }

        stream_sums(st_probs, sc, xyz, r_pts[i].qw_r3_N, ap->number_streams);
    }

    return bg_prob;
}

ALWAYS_INLINE HOT
static inline real bg_probability(const ASTRONOMY_PARAMETERS* ap,
                                  const STREAM_CONSTANTS* sc,
                                  const STREAM_GAUSS* sg,
                                  const LB_TRIG lbt, /* integral point */
                                  const R_CONSTS rc,
                                  const real V,
                                  const R_POINTS* r_pts,
                                  real* st_probs,
                                  KAHAN* probs)
{
    real bg_prob;

    /* if q is 0, there is no probability */
    if (ap->zero_q)
        return -1.0;

    zero_st_probs(st_probs, ap->number_streams);

    bg_prob = ap->bg_prob_func(ap, sc, sg,
                               lbt,
                               rc.gPrime,
                               ap->aux_bg_profile,
                               ap->convolve,
                               r_pts,
                               st_probs);

    sum_probs(probs, st_probs, V * rc.reff_xr_rp3, ap->number_streams);
    bg_prob *= rc.reff_xr_rp3;

    return bg_prob;
}

HOT
static void r_sum(const ASTRONOMY_PARAMETERS* ap,
                  const INTEGRAL_AREA* ia,
                  const STREAM_CONSTANTS* sc,
                  const STREAM_GAUSS* sg,
                  const LB_TRIG lbt,
                  const real id,
                  real* st_probs,
                  KAHAN* probs,
                  EVALUATION_STATE* es,
                  const R_POINTS* r_pts,
                  const R_CONSTS* rc,
                  const unsigned int r_steps)
{
    unsigned int r_step;
    real V;
    real bg_prob;

    for (r_step = 0; r_step < r_steps; ++r_step)
    {
        V = id * rc[r_step].irv;

        bg_prob = V * bg_probability(ap, sc, sg, lbt, rc[r_step], V, &r_pts[r_step * ap->convolve], st_probs, probs);

        KAHAN_ADD(es->sum, bg_prob);
    }
}

ALWAYS_INLINE HOT
static inline void mu_sum(const ASTRONOMY_PARAMETERS* ap,
                          const INTEGRAL_AREA* ia,
                          const STREAM_CONSTANTS* sc,
                          const STREAM_GAUSS* sg,
                          const NU_ID nuid,
                          real* st_probs,
                          KAHAN* probs,
                          const R_POINTS* r_pts,
                          const R_CONSTS* rc,
                          EVALUATION_STATE* es)
{
    real mu;
    LB lb;
    LB_TRIG lbt;

    const real mu_step_size = ia->mu_step_size;
    const real mu_min = ia->mu_min;

    for (; es->mu_step < ia->mu_steps; es->mu_step++)
    {
        do_boinc_checkpoint(es, ia, ap->total_calc_probs);

        mu = mu_min + (((real) es->mu_step + 0.5) * mu_step_size);

        lb = gc2lb(ap->wedge, mu, nuid.nu); /* integral point */
        lbt = lb_trig(lb);

        r_sum(ap, ia, sc, sg, lbt, nuid.id, st_probs, probs, es, r_pts, rc, ia->r_steps);
    }

    es->mu_step = 0;
}

static void nu_sum(const ASTRONOMY_PARAMETERS* ap,
                   const INTEGRAL_AREA* ia,
                   const STREAM_CONSTANTS* sc,
                   const STREAM_GAUSS* sg,
                   real* st_probs,
                   KAHAN* probs,
                   const R_POINTS* r_pts,
                   const R_CONSTS* rc,
                   EVALUATION_STATE* es)
{
    NU_ID nuid;

    for ( ; es->nu_step < ia->nu_steps; es->nu_step++)
    {
        nuid = calc_nu_step(ia, es->nu_step);

        mu_sum(ap, ia, sc, sg,
               nuid,
               st_probs, probs, r_pts, rc, es);
    }
}


/* returns background integral */
real integrate(const ASTRONOMY_PARAMETERS* ap,
               const INTEGRAL_AREA* ia,
               const STREAM_CONSTANTS* sc,
               const STREAM_GAUSS* sg,
               KAHAN* probs,
               EVALUATION_STATE* es)
{
    real result;
    real* st_probs;
    R_POINTS* r_pts;
    R_CONSTS* rc;

    st_probs = (real*) mwMallocAligned(sizeof(real) * ap->number_streams, 2 * sizeof(real));
    r_pts = precalculate_r_pts(ap, ia, sg, &rc);

    nu_sum(ap, ia, sc, sg, st_probs, probs, r_pts, rc, es);
    es->nu_step = 0;

    result = es->sum.sum + es->sum.correction;

    mwAlignedFree(st_probs);
    mwAlignedFree(r_pts);
    mwAlignedFree(rc);

    return result;
}


