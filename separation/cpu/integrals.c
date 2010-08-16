/*
Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

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

#include "separation.h"
#include "evaluation.h"
#include "evaluation_state.h"
#include "integral_constants.h"
#include "integrals_likelihood.h"
#include "integrals.h"
#include "setup_cl.h"

inline static double progress(const EVALUATION_STATE* es,
                              const INTEGRAL_AREA* ia,
                              unsigned int total_calc_probs)
{
    /* This integral's progress */
    /* When checkpointing is done, ia->r_step would always be 0 */
    unsigned int i_prog =  (es->mu_step * ia->nu_steps * ia->r_steps)
                              + (es->nu_step * ia->r_steps); /* + es->r_step */

    return (double)(i_prog + es->current_calc_probs) / total_calc_probs;
}

#if BOINC_APPLICATION
inline static void do_boinc_checkpoint(const EVALUATION_STATE* es,
                                       const INTEGRAL_AREA* ia,
                                       unsigned int total_calc_probs)
{
    if (boinc_time_to_checkpoint())
    {
        if (write_checkpoint(es))
        {
            fprintf(stderr, "APP: write checkpoint failed\n");
            return;
        }
        boinc_checkpoint_completed();
    }

    boinc_fraction_done(progress(es, ia, total_calc_probs));
}

#else

inline static void do_boinc_checkpoint(const EVALUATION_STATE* es,
                                       const INTEGRAL_AREA* ia,
                                       unsigned int total_calc_probs)
{
  #pragma unused(ia)
  #pragma unused(es)
  #pragma unused(total_calc_probs)
}

#endif /* BOINC_APPLICATION */


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


/* Sum over r steps using Kahan summation */
inline static BG_PROB r_sum(const ASTRONOMY_PARAMETERS* ap,
                            const STREAM_CONSTANTS* sc,
                            const R_POINTS* r_pts,
                            const R_CONSTANTS* r_consts,
                            const double nu_consts_id,
                            const unsigned int r_steps,
                            ST_PROBS* probs,
                            vector* xyz,
                            const vector integral_point)
{
    unsigned int r_step_current;
    double V;
    double bg_prob;
    BG_PROB bg_prob_int = ZERO_BG_PROB; /* for Kahan summation */

    for (r_step_current = 0; r_step_current < r_steps; ++r_step_current)
    {
        bg_prob = bg_probability(ap,
                                 &r_pts[r_step_current * ap->convolve],
                                 r_consts[r_step_current].reff_xr_rp3,
                                 integral_point,
                                 xyz);

        V = r_consts[r_step_current].irv * nu_consts_id;
        bg_prob *= V;

        KAHAN_ADD(bg_prob_int.bg_int, bg_prob, bg_prob_int.correction);

        probabilities(ap,
                      sc,
                      &r_pts[r_step_current * ap->convolve],
                      r_consts[r_step_current].reff_xr_rp3,
                      V,
                      xyz,
                      probs);
    }

    return bg_prob_int;
}

inline static void nu_sum(const ASTRONOMY_PARAMETERS* ap,
                          const STREAM_CONSTANTS* sc,
                          const INTEGRAL_AREA* ia,
                          const R_CONSTANTS* r_consts,
                          const R_POINTS* r_pts,
                          const NU_CONSTANTS* nu_consts,
                          const double mu,
                          ST_PROBS* probs,
                          vector* xyz,
                          EVALUATION_STATE* es)
{
    vector integral_point;
    BG_PROB r_result;

    const unsigned int nu_steps = ia->nu_steps;
    const unsigned int r_steps = ia->r_steps;

    for ( ; es->nu_step < nu_steps; es->nu_step++)
    {
        do_boinc_checkpoint(es, ia, ap->total_calc_probs);

        ap->sgr_conversion(ap->wedge,
                           mu + 0.5 * ia->mu_step_size,
                           nu_consts[es->nu_step].nu,
                           &L(integral_point),
                           &B(integral_point));

        r_result = r_sum(ap,
                         sc,
                         r_pts,
                         r_consts,
                         nu_consts[es->nu_step].id,
                         r_steps,
                         probs,
                         xyz,
                         integral_point);

        INCADD_BG_PROB(es->nu_acc, r_result);
    }

    es->nu_step = 0;
}

/* returns background integral */
static double integrate(const ASTRONOMY_PARAMETERS* ap,
                        const STREAM_CONSTANTS* sc,
                        const R_CONSTANTS* r_consts,
                        const R_POINTS* r_pts,
                        const NU_CONSTANTS* nu_consts,
                        const INTEGRAL_AREA* ia,
                        ST_PROBS* probs,
                        vector* xyz,
                        EVALUATION_STATE* es)
{
    double mu;
    const unsigned int mu_steps = ia->mu_steps;

    for ( ; es->mu_step < mu_steps; es->mu_step++)
    {
        mu = ia->mu_min + (es->mu_step * ia->mu_step_size);

        nu_sum(ap, sc, ia, r_consts, r_pts, nu_consts, mu, probs, xyz, es);
        INCADD_BG_PROB(es->mu_acc, es->nu_acc);
        CLEAR_BG_PROB(es->nu_acc);
    }

    es->mu_step = 0;

    return es->mu_acc.bg_int + es->mu_acc.correction;
}

inline static void calculate_stream_integrals(const ST_PROBS* probs,
                                              double* stream_integrals,
                                              const unsigned int number_streams)
{
    unsigned int i;

    for (i = 0; i < number_streams; ++i)
        stream_integrals[i] = probs[i].st_prob_int + probs[i].st_prob_int_c;
}

/* Add up completed integrals for progress reporting */
inline static double completed_integral_progress(const ASTRONOMY_PARAMETERS* ap,
                                                 const EVALUATION_STATE* es)
{
    INTEGRAL_AREA* ia;
    unsigned int i, current_calc_probs = 0;

    for (i = 0; i < es->current_integral; ++i)
    {
        ia = &ap->integral[i];
        current_calc_probs += ia->r_steps * ia->mu_steps * ia->nu_steps;
    }

    return current_calc_probs;
}

void calculate_integrals(const ASTRONOMY_PARAMETERS* ap,
                         const STREAM_CONSTANTS* sc,
                         const STREAM_GAUSS* sg,
                         vector* xyz,
                         EVALUATION_STATE* es)
{
    unsigned int i;
    INTEGRAL_CONSTANTS ic;
    INTEGRAL* integral;
    INTEGRAL_AREA* ia;

    for (; es->current_integral < ap->number_integrals; es->current_integral++)
    {
        integral = &es->integrals[es->current_integral];
        ia = &ap->integral[es->current_integral];
        es->current_calc_probs = completed_integral_progress(ap, es);

        prepare_integral_constants(ap, sg, ia, &ic);

        /* FIXME: This will only work for 1 integral for now */
        //setupSeparationCL(ap, sc, ic.r_step_consts, ic.r_pts, ic.nu_consts, ia);

        //printf("CL Setup\n");
        //mw_finish(EXIT_SUCCESS);

        integral->background_integral = integrate(ap, sc,
                                                  ic.r_step_consts, ic.r_pts, ic.nu_consts,
                                                  ia, integral->probs, xyz, es);

        calculate_stream_integrals(integral->probs, integral->stream_integrals, ap->number_streams);

        free_integral_constants(&ic);
        CLEAR_BG_PROB(es->mu_acc);
    }
}

