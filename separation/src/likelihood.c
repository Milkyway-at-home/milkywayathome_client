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

#include <stdio.h>

#include "separation_types.h"
#include "likelihood.h"
#include "integrals_common.h"
#include "r_points.h"
#include "milkyway_util.h"


/* FIXME: Excessive duplication with stuff used in integrals which I
 * was too lazy to also fix here */


static inline real likelihood_bg_probability_main(const ASTRONOMY_PARAMETERS* ap,
                                                  const STREAM_CONSTANTS* sc,
                                                  const R_POINTS* r_pts,
                                                  const real* sg_dx,
                                                  const LB_TRIG lbt,
                                                  const R_CONSTS rc,
                                                  const unsigned int convolve,
                                                  real* st_probs)
{
    unsigned int i;
    real h_prob, g, rg;
    mwvector xyz;
    real bg_prob = 0.0;

    zero_st_probs(st_probs, ap->number_streams);

    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(ap, r_pts[i].r_point, lbt);
        rg = rg_calc(xyz, ap->q_inv_sqr);

        /* CHECKME: Not having quadratic term on slow one looks like a bug but I'm not sure */
        if (ap->fast_h_prob)
        {
            h_prob = h_prob_fast(ap, r_pts[i].qw_r3_N, rg);
            /* the Hernquist profile includes a quadratic term in g */
            if (ap->aux_bg_profile)
            {
                g = rc.gPrime + sg_dx[i];
                h_prob += aux_prob(ap, r_pts[i].qw_r3_N, g);
            }
        }
        else
        {
            h_prob = h_prob_slow(ap, r_pts[i].qw_r3_N, rg);
        }

        stream_sums(st_probs, sc, xyz, r_pts[i].qw_r3_N, ap->number_streams);

        bg_prob += h_prob;
    }

    mult_probs(st_probs, rc.reff_xr_rp3, ap->number_streams);

    return bg_prob;
}

static inline real likelihood_bg_probability(const ASTRONOMY_PARAMETERS* ap,
                                             const STREAM_CONSTANTS* sc,
                                             const R_POINTS* r_pts,
                                             const real* sg_dx,
                                             const LB_TRIG lbt,
                                             const R_CONSTS rc,
                                             real* st_probs)
{
    real bg_prob;

    /* if q is 0, there is no probability */
    if (ap->zero_q)
        return -1.0;

    bg_prob = likelihood_bg_probability_main(ap, sc, r_pts, sg_dx, lbt, rc, ap->convolve, st_probs);
    bg_prob *= rc.reff_xr_rp3;

    return bg_prob;
}

static inline real stream_sum(const unsigned int number_streams,
                              const FINAL_STREAM_INTEGRALS* fsi,
                              real* st_prob,
                              KAHAN* st_only_sum,
                              const real* exp_stream_weights,
                              const real sum_exp_weights,
                              real bg_only)
{
    unsigned int i;
    real st_only;
    real star_prob = bg_only;

    for (i = 0; i < number_streams; ++i)
    {
        st_only = st_prob[i] / fsi->stream_integrals[i] * exp_stream_weights[i];
        star_prob += st_only;

        if (st_only == 0.0)
            st_only = -238.0;
        else
            st_only = mw_log10(st_only / sum_exp_weights);

        KAHAN_ADD(st_only_sum[i], st_only);
    }
    star_prob /= sum_exp_weights;

    return star_prob;
}

/* Populates exp_stream_weights, and returns the sum */
static inline real get_exp_stream_weights(real* exp_stream_weights,
                                          const STREAMS* streams,
                                          real exp_background_weight)
{
    unsigned int i;
    real sum_exp_weights = exp_background_weight;
    for (i = 0; i < streams->number_streams; i++)
    {
        exp_stream_weights[i] = mw_exp(streams->stream_weight[i].weight);
        sum_exp_weights += exp_stream_weights[i];
    }

    sum_exp_weights *= 0.001;

    return sum_exp_weights;
}

static inline void get_stream_only_likelihood(KAHAN* st_only_sum,
                                              const unsigned int number_stars,
                                              const unsigned int number_streams)
{
    unsigned int i;
    fprintf(stderr, "<stream_only_likelihood>");
    for (i = 0; i < number_streams; i++)
    {
        st_only_sum[i].sum += st_only_sum[i].correction;
        st_only_sum[i].sum /= number_stars;

        fprintf(stderr, " %.20lf", st_only_sum[i].sum - 3.0);
    }
    fprintf(stderr, " </stream_only_likelihood>\n");
}

static real likelihood_sum(const ASTRONOMY_PARAMETERS* ap,
                           const STAR_POINTS* sp,
                           const STREAM_CONSTANTS* sc,
                           const STREAMS* streams,
                           const FINAL_STREAM_INTEGRALS* fsi,
                           const STREAM_GAUSS sg,
                           R_POINTS* r_pts,
                           KAHAN* st_sum,
                           real* st_prob,
                           const real* exp_stream_weights,
                           const real sum_exp_weights,
                           const real exp_background_weight)
{
    KAHAN prob = ZERO_KAHAN;
    KAHAN bg_only = ZERO_KAHAN;

    unsigned int current_star_point;
    real star_prob;
    real bg_prob, bg;
    LB lb;
    LB_TRIG lbt;
    R_CONSTS rc;

    unsigned int num_zero = 0;
    unsigned int bad_jacobians = 0;

    for (current_star_point = 0; current_star_point < sp->number_stars; ++current_star_point)
    {
        rc.gPrime = calcGPrime(Z(sp->stars[current_star_point]));
        set_r_points(ap, sg, ap->convolve, rc.gPrime, r_pts);
        rc.reff_xr_rp3 = calcReffXrRp3(Z(sp->stars[current_star_point]), rc.gPrime);

        LB_L(lb) = L(sp->stars[current_star_point]);
        LB_B(lb) = B(sp->stars[current_star_point]);

        lbt = lb_trig(lb);

        bg_prob = likelihood_bg_probability(ap, sc, r_pts, sg.dx, lbt, rc, st_prob);

        bg = (bg_prob / fsi->background_integral) * exp_background_weight;

        star_prob = stream_sum(streams->number_streams,
                               fsi,
                               st_prob,
                               st_sum,
                               exp_stream_weights,
                               sum_exp_weights,
                               bg);

        if (star_prob != 0.0)
        {
            star_prob = mw_log10(star_prob);
            KAHAN_ADD(prob, star_prob);
        }
        else
        {
            ++num_zero;
            prob.sum -= 238.0;
        }

        if (bg == 0.0)
            bg = -238.0;
        else
            bg = mw_log10(bg / sum_exp_weights);

        KAHAN_ADD(bg_only, bg);
    }

    prob.sum += prob.correction;
    bg_only.sum += bg_only.correction;
    bg_only.sum /= sp->number_stars;
    bg_only.sum -= 3.0;

    fprintf(stderr, "<background_only_likelihood> %.20lf </background_only_likelihood>\n", bg_only.sum);

    /*  log10(x * 0.001) = log10(x) - 3.0 */
    return (prob.sum / (sp->number_stars - bad_jacobians)) - 3.0;
}


real likelihood(const ASTRONOMY_PARAMETERS* ap,
                const STAR_POINTS* sp,
                const STREAM_CONSTANTS* sc,
                const STREAMS* streams,
                const FINAL_STREAM_INTEGRALS* fsi,
                const STREAM_GAUSS sg)

{
    real* st_prob;
    R_POINTS* r_pts;
    KAHAN* st_sum;
    real* exp_stream_weights;
    real sum_exp_weights;
    real likelihood_val;

    const real exp_background_weight = mw_exp(ap->background_weight);

    st_prob = (real*) mwMallocAligned(sizeof(real) * streams->number_streams, 2 * sizeof(real));
    r_pts = (R_POINTS*) mwMallocAligned(sizeof(R_POINTS) * ap->convolve, sizeof(R_POINTS));
    st_sum = (KAHAN*) mwCallocAligned(sizeof(KAHAN), streams->number_streams, sizeof(KAHAN));
    exp_stream_weights = (real*) mwMallocAligned(sizeof(real) * streams->number_streams, 2 * sizeof(real));

    sum_exp_weights = get_exp_stream_weights(exp_stream_weights, streams, exp_background_weight);

    likelihood_val = likelihood_sum(ap, sp, sc, streams, fsi,
                                    sg, r_pts,
                                    st_sum, st_prob,
                                    exp_stream_weights,
                                    sum_exp_weights,
                                    exp_background_weight);

    get_stream_only_likelihood(st_sum, sp->number_stars, streams->number_streams);

    mwAlignedFree(st_prob);
    mwAlignedFree(r_pts);
    mwAlignedFree(st_sum);
    mwAlignedFree(exp_stream_weights);

    return likelihood_val;
}

