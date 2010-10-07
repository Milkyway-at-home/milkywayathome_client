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

/* FIXME: I don't know what these do enough to name it properly */
ALWAYS_INLINE HOT
static inline real sub_bg_probability1(const ASTRONOMY_PARAMETERS* ap,
                                       const STREAM_CONSTANTS* sc,
                                       const STREAM_GAUSS* sg,
                                       const LB_TRIG lbt,
                                       const int aux_bg_profile,
                                       const unsigned int convolve,
                                       const R_POINTS* r_pts,
                                       real* st_probs)
{
    unsigned int i;
    real h_prob, rg, rs;
    vector xyz;
    real bg_prob = 0.0;

    for (i = 0; i < convolve; ++i)
    {
        lbr2xyz_2(xyz, r_pts[i].r_point, lbt);

        rg = rg_calc(xyz, ap->q_inv_sqr);
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


static real sub_bg_probability2(const ASTRONOMY_PARAMETERS* ap,
                                const STREAM_CONSTANTS* sc,
                                const STREAM_GAUSS* sg,
                                const LB_TRIG lbt,
                                const unsigned int convolve,
                                const R_POINTS* r_pts,
                                real* st_probs)
{
    unsigned int i;
    real rg;
    vector xyz;
    real bg_prob = 0.0;

    for (i = 0; i < convolve; ++i)
    {
        lbr2xyz_2(xyz, r_pts[i].r_point, lbt);

        rg = rg_calc(xyz, ap->q_inv_sqr);

        bg_prob += h_prob_slow(ap, r_pts[i].qw_r3_N, rg);
        stream_sums(st_probs, sc, xyz, r_pts[i].qw_r3_N, ap->number_streams);
    }

    return bg_prob;
}

ALWAYS_INLINE HOT
static inline real bg_probability(const ASTRONOMY_PARAMETERS* ap,
                                  const STREAM_CONSTANTS* sc,
                                  const STREAM_GAUSS* sg,
                                  const LB_TRIG lbt, /* integral point */
                                  const real reff_xr_rp3,
                                  const real V,
                                  const R_POINTS* r_pts,
                                  real* st_probs,
                                  ST_PROBS* probs)
{
    real bg_prob;

    /* if q is 0, there is no probability */
    if (ap->zero_q)
        return -1.0;

    zero_st_probs(st_probs, ap->number_streams);

    if (ap->fast_h_prob)
    {
        bg_prob = sub_bg_probability1(ap,
                                      sc,
                                      sg,
                                      lbt,
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
                                      lbt,
                                      ap->convolve,
                                      r_pts,
                                      st_probs);
    }

    sum_probs(probs, st_probs, V * reff_xr_rp3, ap->number_streams);
    bg_prob *= reff_xr_rp3;

    return bg_prob;
}

HOT
static BG_PROB r_sum(const ASTRONOMY_PARAMETERS* ap,
                     const INTEGRAL_AREA* ia,
                     const STREAM_CONSTANTS* sc,
                     const STREAM_GAUSS* sg,
                     const LB_TRIG lbt,
                     const real id,
                     real* st_probs,
                     ST_PROBS* probs,
                     const R_POINTS* r_pts,
                     const R_CONSTS* rc,
                     const unsigned int r_steps)
{
    unsigned int r_step;
    real V;
    real bg_prob;
    BG_PROB bg_prob_int = ZERO_BG_PROB; /* for Kahan summation */

    for (r_step = 0; r_step < r_steps; ++r_step)
    {
        V = id * rc[r_step].irv;

        bg_prob = V * bg_probability(ap, sc, sg, lbt, rc[r_step].reff_xr_rp3, V, &r_pts[r_step * ap->convolve], st_probs, probs);

        KAHAN_ADD(bg_prob_int.bg_int, bg_prob, bg_prob_int.correction);
    }

    return bg_prob_int;
}

/* Sum over mu steps using Kahan summation */
ALWAYS_INLINE HOT
static inline void mu_sum(const ASTRONOMY_PARAMETERS* ap,
                          const INTEGRAL_AREA* ia,
                          const STREAM_CONSTANTS* sc,
                          const STREAM_GAUSS* sg,
                          const NU_ID nuid,
                          real* st_probs,
                          ST_PROBS* probs,
                          const R_POINTS* r_pts,
                          const R_CONSTS* rc,
                          EVALUATION_STATE* es)
{
    real mu;
    LB lb;
    LB_TRIG lbt;
    BG_PROB r_result;

    const unsigned int mu_steps = ia->mu_steps;
    const real mu_step_size = ia->mu_step_size;
    const real mu_min = ia->mu_min;

    for (; es->mu_step < mu_steps; es->mu_step++)
    {
        do_boinc_checkpoint(es, ia, ap->total_calc_probs);

        mu = mu_min + (((real) es->mu_step + 0.5) * mu_step_size);

        lb = gc2lb(ap->wedge, mu, nuid.nu); /* integral point */
        lbt = lb_trig(lb);

        r_result = r_sum(ap, ia, sc, sg, lbt, nuid.id, st_probs, probs, r_pts, rc, ia->r_steps);

        INCADD_BG_PROB(es->mu_acc, r_result);
    }

    es->mu_step = 0;
}

static real nu_sum(const ASTRONOMY_PARAMETERS* ap,
                   const INTEGRAL_AREA* ia,
                   const STREAM_CONSTANTS* sc,
                   const STREAM_GAUSS* sg,
                   real* st_probs,
                   ST_PROBS* probs,
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

        INCADD_BG_PROB(es->nu_acc, es->mu_acc);
        CLEAR_BG_PROB(es->mu_acc);
    }

    return es->nu_acc.bg_int + es->nu_acc.correction;
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
    real* st_probs;
    R_POINTS* r_pts;
    R_CONSTS* rc;

    st_probs = (real*) mwMallocAligned(sizeof(real) * ap->number_streams, 2 * sizeof(real));
    r_pts = precalculate_r_pts(ap, ia, sg, &rc);

    result = nu_sum(ap, ia, sc, sg, st_probs, probs, r_pts, rc, es);
    es->nu_step = 0;

    mwAlignedFree(st_probs);
    mwAlignedFree(r_pts);
    mwAlignedFree(rc);

    return result;
}


