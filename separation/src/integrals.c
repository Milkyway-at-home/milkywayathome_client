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
static inline real progress(const EvaluationState* es,
                            const IntegralArea* ia,
                            real total_calc_probs)
{
    /* This integral's progress */
    /* When checkpointing is done, ia->mu_step would always be 0 */
    unsigned int i_prog =  (es->nu_step * ia->mu_steps * ia->r_steps)
                         + (es->mu_step * ia->r_steps); /* + es->r_step */

    return (real)(i_prog + es->current_calc_probs) / total_calc_probs;
}

ALWAYS_INLINE
static inline void doBoincCheckpoint(const EvaluationState* es,
                                     const IntegralArea* ia,
                                     real total_calc_probs)
{
    if (boinc_time_to_checkpoint())
    {
        if (writeCheckpoint(es))
            fail("Write checkpoint failed\n");
        boinc_checkpoint_completed();
    }

    boinc_fraction_done(progress(es, ia, total_calc_probs));
}

#else

#define doBoincCheckpoint(es, ia, total_calc_probs)

#endif /* BOINC_APPLICATION */


HOT
real bg_probability_fast_hprob(const AstronomyParameters* ap,
                               const StreamConstants* sc,
                               const RPoints* r_pts,
                               const real* sg_dx,
                               const LBTrig lbt,
                               const real gPrime,
                               const int aux_bg_profile,
                               const unsigned int convolve,
                               real* st_probs)
{
    unsigned int i;
    real h_prob, g, rg;
    mwvector xyz;
    real bg_prob = 0.0;

    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(ap, r_pts[i].r_point, lbt);

        rg = rg_calc(ap, xyz);

        h_prob = h_prob_fast(ap, r_pts[i].qw_r3_N, rg);

        /* Add a quadratic term in g to the the Hernquist profile */
        if (aux_bg_profile)
        {
            g = gPrime + sg_dx[i];
            h_prob += aux_prob(ap, r_pts[i].qw_r3_N, g);
        }
        bg_prob += h_prob;

        stream_sums(st_probs, sc, xyz, r_pts[i].qw_r3_N, ap->number_streams);
    }

    return bg_prob;
}


HOT
real bg_probability_slow_hprob(const AstronomyParameters* ap,
                               const StreamConstants* sc,
                               const RPoints* r_pts,
                               const real* sg_dx,
                               const LBTrig lbt,
                               const real gPrime,
                               const int aux_bg_profile,
                               const unsigned int convolve,
                               real* st_probs)
{
    unsigned int i;
    real rg, g;
    mwvector xyz;
    real bg_prob = 0.0;

    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(ap, r_pts[i].r_point, lbt);

        rg = rg_calc(ap, xyz);

        bg_prob += h_prob_slow(ap, r_pts[i].qw_r3_N, rg);
        if (aux_bg_profile)
        {
            g = gPrime + sg_dx[i];
            bg_prob += aux_prob(ap, r_pts[i].qw_r3_N, g);
        }

        stream_sums(st_probs, sc, xyz, r_pts[i].qw_r3_N, ap->number_streams);
    }

    return bg_prob;
}

ALWAYS_INLINE HOT
static inline real bg_probability(const AstronomyParameters* ap,
                                  const StreamConstants* sc,
                                  const RPoints* r_pts,
                                  const real* sg_dx,
                                  const real gPrime,
                                  const LBTrig lbt, /* integral point */
                                  real* st_probs,
                                  Kahan* probs)
{
    real bg_prob;

    zero_st_probs(st_probs, ap->number_streams);

    bg_prob = ap->bg_prob_func(ap, sc, r_pts, sg_dx,
                               lbt,
                               gPrime,
                               ap->aux_bg_profile,
                               ap->convolve,
                               st_probs);

    return bg_prob;
}

HOT
static inline void r_sum(const AstronomyParameters* ap,
                         const StreamConstants* sc,
                         const real* sg_dx,
                         const LBTrig lbt,
                         const real id,
                         real* st_probs,
                         Kahan* probs,
                         EvaluationState* es,
                         const RPoints* r_pts,
                         const RConsts* rc,
                         const unsigned int r_steps)
{
    unsigned int r_step;
    real V_reff_xr_rp3;
    real bg_prob;

    for (r_step = 0; r_step < r_steps; ++r_step)
    {
        bg_prob = bg_probability(ap, sc,
                                 &r_pts[r_step * ap->convolve], sg_dx,
                                 rc[r_step].gPrime, lbt, st_probs, probs);

        V_reff_xr_rp3 = id * rc[r_step].irv_reff_xr_rp3;
        bg_prob *= V_reff_xr_rp3;
        sum_probs(probs, st_probs, V_reff_xr_rp3, ap->number_streams);

        KAHAN_ADD(es->sum, bg_prob);
    }
}

ALWAYS_INLINE HOT
static inline void mu_sum(const AstronomyParameters* ap,
                          const IntegralArea* ia,
                          const StreamConstants* sc,
                          const RConsts* rc,
                          const RPoints* r_pts,
                          const real* sg_dx,
                          const NuId nuid,
                          real* st_probs,
                          Kahan* probs,
                          EvaluationState* es)
{
    real mu;
    LB lb;
    LBTrig lbt;

    const real mu_step_size = ia->mu_step_size;
    const real mu_min = ia->mu_min;

    for (; es->mu_step < ia->mu_steps; es->mu_step++)
    {
        doBoincCheckpoint(es, ia, ap->total_calc_probs);

        mu = mu_min + (((real) es->mu_step + 0.5) * mu_step_size);

        lb = gc2lb(ap->wedge, mu, nuid.nu); /* integral point */
        lbt = lb_trig(lb);

        r_sum(ap, sc, sg_dx, lbt, nuid.id, st_probs, probs, es, r_pts, rc, ia->r_steps);
    }

    es->mu_step = 0;
}

static void nuSum(const AstronomyParameters* ap,
                  const IntegralArea* ia,
                  const StreamConstants* sc,
                  const RConsts* rc,
                  const RPoints* r_pts,
                  const real* sg_dx,
                  real* st_probs,
                  Kahan* probs,
                  EvaluationState* es)
{
    NuId nuid;

    for ( ; es->nu_step < ia->nu_steps; es->nu_step++)
    {
        nuid = calcNuStep(ia, es->nu_step);

        mu_sum(ap, ia, sc, rc, r_pts, sg_dx, nuid, st_probs, probs, es);
    }
}

/* returns background integral */
real integrate(const AstronomyParameters* ap,
               const IntegralArea* ia,
               const StreamConstants* sc,
               const StreamGauss sg,
               real* probs,
               Kahan* probs_sum,
               EvaluationState* es)
{
    real result;
    real* st_probs;
    RPoints* r_pts;
    RConsts* rc;
    unsigned int i;

    if (ap->q == 0.0)
    {
        /* if q is 0, there is no probability */
        /* Short circuit the entire integral rather than add up -1 many times. */
        warn("q is 0.0\n");
        return -1.0 * ia->nu_steps * ia->mu_steps * ia->r_steps;
    }

    st_probs = (real*) mwMallocA(sizeof(real) * ap->number_streams);
    r_pts = precalculateRPts(ap, ia, sg, &rc, 0);

    nuSum(ap, ia, sc, rc, r_pts, sg.dx, st_probs, probs_sum, es);
    es->nu_step = 0;

    result = es->sum.sum + es->sum.correction;
    for (i  = 0; i < ap->number_streams; ++i)
        probs[i] = probs_sum[i].sum + probs_sum[i].correction;

    mwFreeA(st_probs);
    mwFreeA(r_pts);
    mwFreeA(rc);

    return result;
}


