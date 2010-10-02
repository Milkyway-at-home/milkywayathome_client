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
#include "integrals_likelihood.h"
#include "r_points.h"
#include "milkyway_util.h"


/* FIXME: Excessive duplication with stuff used in integrals which I
 * was too lazy to also fix here */

static inline real probabilities_convolve(__MW_CONSTANT STREAM_CONSTANTS* sc,
                                          __MW_LOCAL const R_POINTS* r_pts,
                                          __MW_LOCAL vector* const xyz,
                                          const unsigned int convolve)
{
    unsigned int i;
    real st_prob = 0.0;

    for (i = 0; i < convolve; ++i)
        st_prob += calc_st_prob_inc(sc, xyz[i], r_pts[i].qw_r3_N);

    return st_prob;
}

static inline void likelihood_probabilities(const ASTRONOMY_PARAMETERS* ap,
                                            const STREAM_CONSTANTS* sc,
                                            const R_POINTS* r_pts,
                                            const real reff_xr_rp3,
                                            vector* const xyz,
                                            real* probs)
{
    unsigned int i;

    for (i = 0; i < ap->number_streams; ++i)
    {
        if (sc[i].large_sigma)
            probs[i] = reff_xr_rp3 * probabilities_convolve(&sc[i], r_pts, xyz, ap->convolve);
        else
            probs[i] = 0.0;
    }
}

static inline real likelihood_bg_probability_main(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                                                  __MW_LOCAL const R_POINTS* r_pts,
                                                  __MW_LOCAL vector* const xyz,
                                                  const LB integral_point,
                                                  const unsigned int convolve)
{
    unsigned int i;
    real h_prob;
    real rg, rs;
    real lsin, lcos;
    real bsin, bcos;
    real bg_prob = 0.0;

    mw_sincos(d2r(LB_L(integral_point)), &lsin, &lcos);
    mw_sincos(d2r(LB_B(integral_point)), &bsin, &bcos);

    for (i = 0; i < convolve; ++i)
    {
        lbr2xyz_2(xyz[i], r_pts[i].r_point, bsin, bcos, lsin, lcos);

        rg = rg_calc(xyz[i], ap->q_inv_sqr);

        rs = rg + ap->r0;

        /* CHECKME: Not having quadratic term on slow one looks like a bug but I'm not sure */
        if (ap->fast_h_prob)
        {
            h_prob = h_prob_fast(r_pts[i].qw_r3_N, rg, rs);

            /* The Hernquist profile includes a quadratic term in g */
            if (ap->aux_bg_profile)
                h_prob += aux_prob(ap, r_pts[i].qw_r3_N, r_pts[i].r_in_mag, r_pts[i].r_in_mag2);
        }
        else
        {
            h_prob = h_prob_slow(ap, r_pts[i].qw_r3_N, rg);
        }

        bg_prob += h_prob;
    }

    return bg_prob;
}

static inline real likelihood_bg_probability(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                                             __MW_LOCAL const R_POINTS* r_pts,
                                             __MW_LOCAL vector* const xyz,
                                             const LB integral_point,
                                             const real reff_xr_rp3)
{
    real bg_prob;

    /* if q is 0, there is no probability */
    if (ap->zero_q)
        return -1.0;

    bg_prob = likelihood_bg_probability_main(ap, r_pts, xyz, integral_point, ap->convolve);
    bg_prob *= reff_xr_rp3;

    return bg_prob;
}

static inline real stream_sum(const unsigned int number_streams,
                              const FINAL_STREAM_INTEGRALS* fsi,
                              real* st_prob,
                              ST_SUM* st_sum,
                              const real* exp_stream_weights,
                              const real sum_exp_weights,
                              real bg_only)
{
    unsigned int current_stream;
    real st_only;
    real star_prob = bg_only;

    for (current_stream = 0; current_stream < number_streams; current_stream++)
    {
        st_only = st_prob[current_stream] / fsi->stream_integrals[current_stream] * exp_stream_weights[current_stream];
        star_prob += st_only;

        if (st_only == 0.0)
            st_only = -238.0;
        else
            st_only = mw_log10(st_only / sum_exp_weights);

        KAHAN_ADD(st_sum[current_stream].st_only_sum, st_only, st_sum[current_stream].st_only_sum_c);
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

static inline void get_stream_only_likelihood(ST_SUM* st_sum,
                                              const unsigned int number_stars,
                                              const unsigned int number_streams)
{
    unsigned int i;
    fprintf(stderr, "<stream_only_likelihood>");
    for (i = 0; i < number_streams; i++)
    {
        st_sum[i].st_only_sum += st_sum[i].st_only_sum_c;
        st_sum[i].st_only_sum /= number_stars;

        fprintf(stderr, " %.20lf", st_sum[i].st_only_sum - 3.0);
    }
    fprintf(stderr, " </stream_only_likelihood>\n");
}

static real likelihood_sum(const ASTRONOMY_PARAMETERS* ap,
                           const STAR_POINTS* sp,
                           const STREAM_CONSTANTS* sc,
                           const STREAMS* streams,
                           const FINAL_STREAM_INTEGRALS* fsi,
                           const STREAM_GAUSS* sg,
                           R_POINTS* r_pts,
                           ST_SUM* st_sum,
                           vector* xyz,
                           real* st_prob,
                           const real* exp_stream_weights,
                           const real sum_exp_weights,
                           const real exp_background_weight)
{
    PROB_SUM prob = ZERO_PROB_SUM;
    PROB_SUM bg_only = ZERO_PROB_SUM;

    unsigned int current_star_point;
    real star_prob;
    real bg_prob, bg, reff_xr_rp3;
    LB lb;
    real gPrime;

    unsigned int num_zero = 0;
    unsigned int bad_jacobians = 0;

    for (current_star_point = 0; current_star_point < sp->number_stars; ++current_star_point)
    {
        gPrime = calcGPrime(Z(sp->stars[current_star_point]));
        set_r_points(ap, sg, ap->convolve, gPrime, r_pts);
        reff_xr_rp3 = calcReffXrRp3(Z(sp->stars[current_star_point]), gPrime);

        LB_L(lb) = L(sp->stars[current_star_point]);
        LB_B(lb) = B(sp->stars[current_star_point]);

        bg_prob = likelihood_bg_probability(ap, r_pts, xyz, lb, reff_xr_rp3);

        bg = (bg_prob / fsi->background_integral) * exp_background_weight;

        likelihood_probabilities(ap, sc, r_pts, reff_xr_rp3, xyz, st_prob);

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
            KAHAN_ADD(prob.sum, star_prob, prob.correction);
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

        KAHAN_ADD(bg_only.sum, bg, bg_only.correction);
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
                const STREAM_GAUSS* sg)

{
    real* st_prob;
    R_POINTS* r_pts;
    ST_SUM* st_sum;
    real* exp_stream_weights;
    vector* xyzs;
    real sum_exp_weights;
    real likelihood_val;

    const real exp_background_weight = mw_exp(ap->background_weight);

    st_prob = (real*) mallocSafe(sizeof(real) * streams->number_streams);
    r_pts = (R_POINTS*) mallocSafe(sizeof(R_POINTS) * ap->convolve);
    st_sum = (ST_SUM*) callocSafe(sizeof(ST_SUM), streams->number_streams);
    exp_stream_weights = (real*) mallocSafe(sizeof(real) * streams->number_streams);
    xyzs = (vector*) callocSafe(sizeof(vector), ap->convolve);

    sum_exp_weights = get_exp_stream_weights(exp_stream_weights, streams, exp_background_weight);

    likelihood_val = likelihood_sum(ap, sp, sc, streams, fsi,
                                    sg, r_pts,
                                    st_sum, xyzs, st_prob,
                                    exp_stream_weights,
                                    sum_exp_weights,
                                    exp_background_weight);

    get_stream_only_likelihood(st_sum, sp->number_stars, streams->number_streams);

    free(st_prob);
    free(r_pts);
    free(st_sum);
    free(exp_stream_weights);
    free(xyzs);

    return likelihood_val;
}

