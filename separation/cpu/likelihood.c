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

#include "likelihood.h"
#include "integrals_likelihood.h"
#include "r_points.h"
#include "milkyway_util.h"

inline static void likelihood_probabilities(const ASTRONOMY_PARAMETERS* ap,
                                            const STREAM_CONSTANTS* sc,
                                            const R_POINTS* r_pts,
                                            const double reff_xr_rp3,
                                            vector* const xyz,
                                            double* probs)
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

inline static double stream_sum(const unsigned int number_streams,
                                const FINAL_STREAM_INTEGRALS* fsi,
                                double* st_prob,
                                ST_SUM* st_sum,
                                const double* exp_stream_weights,
                                const double sum_exp_weights,
                                double bg_only)
{
    unsigned int current_stream;
    double st_only;
    double star_prob = bg_only;

    for (current_stream = 0; current_stream < number_streams; current_stream++)
    {
        st_only = st_prob[current_stream] / fsi->stream_integrals[current_stream] * exp_stream_weights[current_stream];
        star_prob += st_only;

        if (st_only == 0.0)
            st_only = -238.0;
        else
            st_only = log10(st_only / sum_exp_weights);

        KAHAN_ADD(st_sum[current_stream].st_only_sum, st_only, st_sum[current_stream].st_only_sum_c);
    }
    star_prob /= sum_exp_weights;

    return star_prob;
}

/* Populates exp_stream_weights, and returns the sum */
inline static double get_exp_stream_weights(double* exp_stream_weights,
                                            const STREAMS* streams,
                                            double exp_background_weight)
{
    unsigned int i;
    double sum_exp_weights = exp_background_weight;
    for (i = 0; i < streams->number_streams; i++)
    {
        exp_stream_weights[i] = exp(streams->stream_weight[i].weight);
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


static double likelihood_sum(const ASTRONOMY_PARAMETERS* ap,
                             const STAR_POINTS* sp,
                             const STREAM_CONSTANTS* sc,
                             const STREAMS* streams,
                             const FINAL_STREAM_INTEGRALS* fsi,
                             R_POINTS* r_pts,
                             STREAM_GAUSS* sg,
                             ST_SUM* st_sum,
                             vector* xyz,
                             double* st_prob,
                             const double* exp_stream_weights,
                             const double sum_exp_weights,
                             const double exp_background_weight)
{
    PROB_SUM prob = ZERO_PROB_SUM;
    PROB_SUM bg_only = ZERO_PROB_SUM;

    unsigned int current_star_point;
    double star_prob;
    double bg_prob, bg, reff_xr_rp3;
    LB lb;

    unsigned int num_zero = 0;
    unsigned int bad_jacobians = 0;



    for (current_star_point = 0; current_star_point < sp->number_stars; ++current_star_point)
    {
        reff_xr_rp3 = set_r_points(ap, sg, ap->convolve, ZN(sp, current_star_point), r_pts);

        LB_L(lb) = LN(sp, current_star_point);
        LB_B(lb) = BN(sp, current_star_point);

        bg_prob = bg_probability(ap, r_pts, reff_xr_rp3, lb, xyz);

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
            star_prob = log10(star_prob);
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
            bg = log10(bg / sum_exp_weights);

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


double likelihood(const ASTRONOMY_PARAMETERS* ap,
                  const STAR_POINTS* sp,
                  const STREAM_CONSTANTS* sc,
                  const STREAMS* streams,
                  const FINAL_STREAM_INTEGRALS* fsi,
                  STREAM_GAUSS* sg)

{
    double* st_prob = mallocSafe(sizeof(double) * streams->number_streams);
    R_POINTS* r_pts = mallocSafe(sizeof(R_POINTS) * ap->convolve);
    ST_SUM* st_sum = callocSafe(sizeof(ST_SUM), streams->number_streams);
    double* exp_stream_weights = mallocSafe(sizeof(double) * streams->number_streams);
    vector* xyzs = mallocSafe(sizeof(vector) * ap->convolve);

    const double exp_background_weight = exp(ap->background_weight);
    double sum_exp_weights = get_exp_stream_weights(exp_stream_weights, streams, exp_background_weight);

    double likelihood_val = likelihood_sum(ap, sp, sc, streams, fsi,
                                           r_pts, sg,
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

