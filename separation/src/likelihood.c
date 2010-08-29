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

inline static real probabilities_convolve(__MW_CONSTANT STREAM_CONSTANTS* sc,
                                          __MW_LOCAL const R_POINTS* r_pts,
                                          __MW_LOCAL vector* const xyz,
                                          const unsigned int convolve)
{
    unsigned int i;
    real dotted, xyz_norm;
    vector xyzs;

    real st_prob = 0.0;

    for (i = 0; i < convolve; ++i)
    {
        SUBV(xyzs, xyz[i], sc->c);
        DOTVP(dotted, sc->a, xyzs);
        INCSUBVMS(xyzs, dotted, sc->a);
        SQRV(xyz_norm, xyzs);

        st_prob += r_pts[i].qw_r3_N * mw_exp(-xyz_norm / sc->sigma_sq2);
    }

    return st_prob;
}

inline static void likelihood_probabilities(const ASTRONOMY_PARAMETERS* ap,
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

inline static real likelihood_bg_probability_main(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                                                  __MW_LOCAL const R_POINTS* r_pts,
                                                  __MW_LOCAL vector* const xyz,
                                                  const LB integral_point,
                                                  const int aux_bg_profile,
                                                  const unsigned int convolve)
{
    unsigned int i;
    real h_prob, aux_prob;
    real rg, rs;
    real lsin, lcos;
    real bsin, bcos;
    real bg_prob = 0.0;

    mw_sincos(d2r(LB_L(integral_point)), &lsin, &lcos);
    mw_sincos(d2r(LB_B(integral_point)), &bsin, &bcos);

    for (i = 0; i < convolve; ++i)
    {
        lbr2xyz_2(xyz[i], r_pts[i].r_point, bsin, bcos, lsin, lcos);

        rg = mw_sqrt(sqr(X(xyz[i])) + sqr(Y(xyz[i])) + sqr(Z(xyz[i])) / sqr(ap->q));
        rs = rg + ap->r0;

        h_prob = r_pts[i].qw_r3_N / (rg * cube(rs));

        //the hernquist profile includes a quadratic term in g
        if (aux_bg_profile)
        {
            aux_prob = r_pts[i].qw_r3_N * (  ap->bg_a * r_pts[i].r_in_mag2
                                           + ap->bg_b * r_pts[i].r_in_mag
                                           + ap->bg_c );
            h_prob += aux_prob;
        }
        bg_prob += h_prob;
    }

    return bg_prob;
}

inline static real likelihood_bg_probability_full(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                                                  __MW_LOCAL const R_POINTS* r_pts,
                                                  __MW_LOCAL vector* const xyz,
                                                  const LB integral_point,
                                                  const unsigned int convolve)
{
    unsigned int i;
    real rg;
    real lsin, lcos;
    real bsin, bcos;
    real bg_prob = 0.0;

    mw_sincos(d2r(LB_L(integral_point)), &lsin, &lcos);
    mw_sincos(d2r(LB_B(integral_point)), &bsin, &bcos);

    for (i = 0; i < convolve; ++i)
    {
        lbr2xyz_2(xyz[i], r_pts[i].r_point, bsin, bcos, lsin, lcos);

        rg = mw_sqrt(sqr(X(xyz[i])) + sqr(Y(xyz[i])) + sqr(Z(xyz[i])) / sqr(ap->q));

        bg_prob += r_pts[i].qw_r3_N / (mw_powr(rg, ap->alpha) * mw_powr(rg + ap->r0, ap->alpha_delta3));
    }

    return bg_prob;
}

inline static real likelihood_bg_probability(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                                             __MW_LOCAL const R_POINTS* r_pts,
                                             __MW_LOCAL vector* const xyz,
                                             const LB integral_point,
                                             const real reff_xr_rp3)
{
    real bg_prob;

    /* if q is 0, there is no probability */
    if (ap->q == 0)
        bg_prob = -1.0;
    else
    {
        if (ap->alpha == 1 && ap->delta == 1)
        {
            bg_prob = likelihood_bg_probability_main(ap,
                                                     r_pts,
                                                     xyz,
                                                     integral_point,
                                                     ap->aux_bg_profile,
                                                     ap->convolve);
        }
        else
            bg_prob = likelihood_bg_probability_full(ap, r_pts, xyz, integral_point, ap->convolve);

        bg_prob *= reff_xr_rp3;
    }

    return bg_prob;
}

inline static real stream_sum(const unsigned int number_streams,
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
inline static real get_exp_stream_weights(real* exp_stream_weights,
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

inline static void get_stream_only_likelihood(ST_SUM* st_sum,
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

    unsigned int num_zero = 0;
    unsigned int bad_jacobians = 0;

    for (current_star_point = 0; current_star_point < sp->number_stars; ++current_star_point)
    {
        reff_xr_rp3 = set_r_points(ap, sg, ap->convolve, ZN(sp, current_star_point), r_pts);

        LB_L(lb) = LN(sp, current_star_point);
        LB_B(lb) = BN(sp, current_star_point);

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
    real* st_prob = mallocSafe(sizeof(real) * streams->number_streams);
    R_POINTS* r_pts = mallocSafe(sizeof(R_POINTS) * ap->convolve);
    ST_SUM* st_sum = callocSafe(sizeof(ST_SUM), streams->number_streams);
    real* exp_stream_weights = mallocSafe(sizeof(real) * streams->number_streams);
    vector* xyzs = callocSafe(sizeof(vector), ap->convolve);

    const real exp_background_weight = mw_exp(ap->background_weight);
    real sum_exp_weights = get_exp_stream_weights(exp_stream_weights, streams, exp_background_weight);

    real likelihood_val = likelihood_sum(ap, sp, sc, streams, fsi,
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

