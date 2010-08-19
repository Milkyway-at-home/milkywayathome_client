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
#include "coordinates.h"
#include "r_points.h"


#if BOINC_APPLICATION

__attribute ((always_inline))
inline static double progress(const EVALUATION_STATE* es,
                              const INTEGRAL_AREA* ia,
                              unsigned int total_calc_probs)
{
    /* This integral's progress */
    /* When checkpointing is done, ia->mu_step would always be 0 */
    unsigned int i_prog =  (es->r_step * ia->nu_steps * ia->mu_steps)
                         + (es->nu_step * ia->mu_steps); /* + es->mu_step */

    return (double)(i_prog + es->current_calc_probs) / total_calc_probs;
}

__attribute ((always_inline))
inline static void do_boinc_checkpoint(const EVALUATION_STATE* es,
                                       const INTEGRAL_AREA* ia,
                                       unsigned int total_calc_probs)
{
    if (boinc_time_to_checkpoint())
    {
        if (write_checkpoint(es))
        {
            warn("Write checkpoint failed\n");
            return;
        }
        boinc_checkpoint_completed();
    }

    boinc_fraction_done(progress(es, ia, total_calc_probs));
}

#else

#define do_boinc_checkpoint(es, ia, total_calc_probs)

#endif /* BOINC_APPLICATION */

__attribute__ ((always_inline))
inline static void probabilities(const ASTRONOMY_PARAMETERS* ap,
                                 const STREAM_CONSTANTS* sc,
                                 const R_POINTS* r_pts,
                                 const double reff_xr_rp3,
                                 const double V,
                                 vector* const xyz,
                                 ST_PROBS* probs)
{
    unsigned int i;
    double st_prob;

    for (i = 0; i < ap->number_streams; ++i)
    {
        if (sc[i].large_sigma)
            st_prob = V * reff_xr_rp3 * probabilities_convolve(&sc[i], r_pts, xyz, ap->convolve);
        else
            st_prob = 0.0;

        KAHAN_ADD(probs[i].st_prob_int, st_prob, probs[i].st_prob_int_c);
    }
}

/* Sum over mu steps using Kahan summation */
__attribute__ ((always_inline, hot))
inline static BG_PROB mu_sum(const ASTRONOMY_PARAMETERS* ap,
                             const STREAM_CONSTANTS* sc,
                             const R_POINTS* r_pts,
                             const double irv,             /* r constants */
                             const double reff_xr_rp3,
                             const double nu_consts_id,    /* nu constants */
                             const double nu_consts_nu,
                             const unsigned int mu_steps,
                             const double mu_step_size,
                             const double mu_min,
                             ST_PROBS* probs,
                             vector* xyz)
{
    unsigned int mu_step_current;
    double mu, V;
    double bg_prob;
    BG_PROB bg_prob_int = ZERO_BG_PROB; /* for Kahan summation */
    LB lb;

    for (mu_step_current = 0; mu_step_current < mu_steps; ++mu_step_current)
    {
        mu = mu_min + (((double) mu_step_current + 0.5) * mu_step_size);

        lb = gc2lb(ap->wedge, mu, nu_consts_nu);

        bg_prob = bg_probability(ap, r_pts, reff_xr_rp3, lb, xyz);

        V = irv * nu_consts_id;
        bg_prob *= V;

        KAHAN_ADD(bg_prob_int.bg_int, bg_prob, bg_prob_int.correction);

        probabilities(ap, sc, r_pts, reff_xr_rp3, V, xyz, probs);
    }

    return bg_prob_int;
}

__attribute__ ((always_inline))
inline static void nu_sum(const ASTRONOMY_PARAMETERS* ap,
                          const STREAM_CONSTANTS* sc,
                          const INTEGRAL_AREA* ia,
                          const double irv,
                          const double reff_xr_rp3,
                          const R_POINTS* r_pts,
                          const NU_CONSTANTS* nu_consts,
                          ST_PROBS* probs,
                          vector* xyz,
                          EVALUATION_STATE* es)
{
    BG_PROB mu_result;

    const unsigned int nu_steps = ia->nu_steps;
    const unsigned int mu_steps = ia->mu_steps;
    const double mu_min = ia->mu_min;
    const double mu_step_size = ia->mu_step_size;

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
                           probs,
                           xyz);

        INCADD_BG_PROB(es->nu_acc, mu_result);
    }

    es->nu_step = 0;
}

double r_sum(const ASTRONOMY_PARAMETERS* ap,
             const INTEGRAL_AREA* ia,
             const STREAM_CONSTANTS* sc,
             const STREAM_GAUSS* sg,
             const NU_CONSTANTS* nu_consts,
             R_POINTS* r_pts,
             ST_PROBS* probs,
             vector* xyz,
             EVALUATION_STATE* es)
{
    double r, next_r, rPrime;
    double irv, reff_xr_rp3;

  #ifdef USE_KPC
    const double r_max           = ia->r_min + ia->r_step_size * r_steps;
    const double r_min_kpc       = pow(10.0, ((ia->r_min - 14.2) / 5.0));
    const double r_max_kpc       = pow(10.0, ((ia->r_max - 14.2) / 5.0));
    const double r_step_size_kpc = (r_max_kpc - r_min_kpc) / r_steps;
  #endif

    const unsigned int r_steps = ia->r_steps;

    for ( ; es->r_step < r_steps; es->r_step++)
    {
      #ifdef USE_KPC
        r = r_min_kpc + (es->r_step * r_step_size_kpc);
        next_r = r + r_step_size_kpc;
      #else
        double log_r = ia->r_min + (es->r_step * ia->r_step_size);
        r = pow(10.0, (log_r - 14.2) / 5.0);
        next_r = pow(10.0, (log_r + ia->r_step_size - 14.2) / 5.0);
      #endif

        irv = d2r(((cube(next_r) - cube(r)) / 3.0) * ia->mu_step_size);
        rPrime = (next_r + r) / 2.0;

        reff_xr_rp3 = set_r_points(ap, sg, ap->convolve, rPrime, r_pts);

        nu_sum(ap, sc, ia, irv, reff_xr_rp3, r_pts, nu_consts, probs, xyz, es);

        INCADD_BG_PROB(es->r_acc, es->nu_acc);
        CLEAR_BG_PROB(es->nu_acc);
    }

    return es->r_acc.bg_int + es->r_acc.correction;
}

