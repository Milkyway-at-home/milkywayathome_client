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

#include "setup_cl.h"

#define stdev 0.6
#define xr (3.0 * stdev)
#define absm 4.2
#define SIGMA_LIMIT 0.0001


inline static double progress(const ASTRONOMY_PARAMETERS* ap, const EVALUATION_STATE* es)
{
    unsigned int i;
    unsigned int current_calc_probs = 0;
    const INTEGRAL_AREA* ia;

    /* Add up completed integrals */
    for (i = 0; i < es->current_integral; i++)
    {
        ia = &ap->integral[i];
        current_calc_probs += ia->r_steps * ia->mu_steps * ia->nu_steps;
    }

    ia = &ap->integral[es->current_integral];

    /* When checkpointing is done, ia->r_step would always be 0 */
    current_calc_probs +=   (es->mu_step * ia->nu_steps * ia->r_steps)
                          + (es->nu_step * ia->r_steps); /* + es->r_step */

    return (double)current_calc_probs / ap->total_calc_probs;
}

#if BOINC_APPLICATION
inline static void do_boinc_checkpoint(const ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es)
{
    double frac;

    static unsigned int i = 0;

    i = (i + 1) % 1000;

    //if (boinc_time_to_checkpoint())
    if ( i == 100 )
    {
        if (write_checkpoint(es))
        {
            fprintf(stderr, "APP: write checkpoint failed\n");
            return;
        }
        boinc_checkpoint_completed();
    }

    frac = progress(ap, es);
    boinc_fraction_done(frac);
}

#else

inline static void do_boinc_checkpoint(const ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es)
{
  #pragma unused(ap)
  #pragma unused(es)
}

#endif /* BOINC_APPLICATION */


/* FIXME: I don't know what these do enough to name it properly */
inline static double sub_bg_probability1(const ASTRONOMY_PARAMETERS* ap,
                                         const R_POINTS* r_pts,
                                         const unsigned int convolve,
                                         const int aux_bg_profile,
                                         const vector integral_point,
                                         vector* const xyz)
{
    unsigned int i;
    double h_prob, aux_prob;
    double rg, rs;
    double zp;
    double bg_prob = 0.0;

    const double lsin = sin(d2r(L(integral_point)));
    const double lcos = cos(d2r(L(integral_point)));
    const double bsin = sin(d2r(B(integral_point)));
    const double bcos = cos(d2r(B(integral_point)));

    for (i = 0; i < convolve; ++i)
    {
        Z(xyz[i]) = r_pts[i].r_point * bsin;
        zp = r_pts[i].r_point * bcos;
        X(xyz[i]) = zp * lcos - sun_r0;
        Y(xyz[i]) = zp * lsin;

        rg = sqrt(sqr(X(xyz[i])) + sqr(Y(xyz[i])) + sqr(Z(xyz[i])) / sqr(ap->q));
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

inline static double sub_bg_probability2(const ASTRONOMY_PARAMETERS* ap,
                                         const R_POINTS* r_pts,
                                         const unsigned int convolve,
                                         const vector integral_point,
                                         vector* const xyz)
{
    unsigned int i;
    double bg_prob = 0.0;
    double rg;
    double zp;

    const double lsin = sin(d2r(L(integral_point)));
    const double lcos = cos(d2r(L(integral_point)));
    const double bsin = sin(d2r(B(integral_point)));
    const double bcos = cos(d2r(B(integral_point)));

    for (i = 0; i < convolve; ++i)
    {
        Z(xyz[i]) = r_pts[i].r_point * bsin;
        zp = r_pts[i].r_point * bcos;
        X(xyz[i]) = zp * lcos - sun_r0;
        Y(xyz[i]) = zp * lsin;

        rg = sqrt(sqr(X(xyz[i])) + sqr(Y(xyz[i])) + sqr(Z(xyz[i])) / sqr(ap->q));

        bg_prob += r_pts[i].qw_r3_N / (pow(rg, ap->alpha) * pow(rg + ap->r0, ap->alpha_delta3));
    }

    return bg_prob;
}

inline static double bg_probability(const ASTRONOMY_PARAMETERS* ap,
                                    const R_POINTS* r_pts,
                                    const double reff_xr_rp3,
                                    const vector integral_point,
                                    vector* xyz)
{
    double bg_prob;

    /* if q is 0, there is no probability */
    if (ap->q == 0)
        bg_prob = -1.0;
    else
    {
        if (ap->alpha == 1 && ap->delta == 1)
            bg_prob = sub_bg_probability1(ap, r_pts, ap->convolve, ap->aux_bg_profile, integral_point, xyz);
        else
            bg_prob = sub_bg_probability2(ap, r_pts, ap->convolve, integral_point, xyz);

        bg_prob *= reff_xr_rp3;
    }

    return bg_prob;
}

/* FIXME: Better name? */
inline static double probabilities_convolve(const STREAM_CONSTANTS* sc,
                                            const R_POINTS* r_pts,
                                            const unsigned int convolve,
                                            vector* const xyz)
{
    unsigned int i;
    double dotted, xyz_norm;
    vector xyzs;

    double st_prob = 0.0;

    for (i = 0; i < convolve; i++)
    {
        X(xyzs) = X(xyz[i]) - X(sc->c);
        Y(xyzs) = Y(xyz[i]) - Y(sc->c);
        Z(xyzs) = Z(xyz[i]) - Z(sc->c);

        dotted = X(sc->a) * X(xyzs)
               + Y(sc->a) * Y(xyzs)
               + Z(sc->a) * Z(xyzs);

        X(xyzs) -= dotted * X(sc->a);
        Y(xyzs) -= dotted * Y(sc->a);
        Z(xyzs) -= dotted * Z(sc->a);

        xyz_norm = sqr(X(xyzs)) + sqr(Y(xyzs)) + sqr(Z(xyzs));

        st_prob += r_pts[i].qw_r3_N * exp(-xyz_norm / sc->sigma_sq2);
    }

    return st_prob;
}

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
            st_prob = V * reff_xr_rp3 * probabilities_convolve(&sc[i], r_pts, ap->convolve, xyz);
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
        do_boinc_checkpoint(ap, es);

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

void calculate_integrals(const ASTRONOMY_PARAMETERS* ap,
                         const STREAM_CONSTANTS* sc,
                         const STREAM_GAUSS* sg,
                         EVALUATION_STATE* es,
                         vector* xyz)
{
    INTEGRAL_CONSTANTS ic;
    INTEGRAL* integral;
    INTEGRAL_AREA* ia;

    for (; es->current_integral < ap->number_integrals; es->current_integral++)
    {
        integral = &es->integrals[es->current_integral];
        ia = &ap->integral[es->current_integral];

        prepare_integral_constants(ap, sg, ia, &ic);

        setupSeparationCL(ap, sc, ic.r_step_consts, ic.r_pts, ic.nu_consts, ia);
        printf("CL Setup\n");
        mw_finish(EXIT_SUCCESS);

        integral->background_integral = integrate(ap, sc,
                                                  ic.r_step_consts, ic.r_pts, ic.nu_consts,
                                                  ia, integral->probs, xyz, es);

        calculate_stream_integrals(integral->probs, integral->stream_integrals, ap->number_streams);

        free_integral_constants(&ic);
        CLEAR_BG_PROB(es->mu_acc);
    }
}

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
            probs[i] = reff_xr_rp3 * probabilities_convolve(&sc[i], r_pts, ap->convolve, xyz);
        else
            probs[i] = 0.0;
    }
}

/* Used in likelihood calculation */
inline static double stream_sum(const unsigned int number_streams,
                                EVALUATION_STATE* es,
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
        st_only = st_prob[current_stream] / es->stream_integrals[current_stream] * exp_stream_weights[current_stream];
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
                             const STREAM_CONSTANTS* sc,
                             const STAR_POINTS* sp,
                             const STREAMS* streams,
                             R_POINTS* r_pts,
                             EVALUATION_STATE* es,
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

    unsigned int num_zero = 0;
    unsigned int bad_jacobians = 0;

    for (current_star_point = 0; current_star_point < sp->number_stars; ++current_star_point)
    {
        reff_xr_rp3 = set_prob_consts(ap, &sg[0], ap->convolve, ZN(sp, current_star_point), r_pts);

        bg_prob = bg_probability(ap, r_pts,
                                 reff_xr_rp3, &VN(sp, current_star_point), xyz);

        bg = (bg_prob / es->background_integral) * exp_background_weight;

        likelihood_probabilities(ap, sc, r_pts, reff_xr_rp3, xyz, st_prob);

        star_prob = stream_sum(streams->number_streams,
                               es,
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
                  const STREAM_CONSTANTS* sc,
                  const STREAMS* streams,
                  EVALUATION_STATE* es,
                  STREAM_GAUSS* sg,
                  vector* xyz,
                  const STAR_POINTS* sp)
{
    double sum_exp_weights;

    double* st_prob = malloc(sizeof(double) * streams->number_streams);
    R_POINTS* r_pts = malloc(sizeof(R_POINTS) * ap->convolve);
    ST_SUM* st_sum = calloc(sizeof(ST_SUM), streams->number_streams);
    double* exp_stream_weights = malloc(sizeof(double) * streams->number_streams);

    const double exp_background_weight = exp(ap->background_weight);
    sum_exp_weights = get_exp_stream_weights(exp_stream_weights, streams, exp_background_weight);

    double likelihood_val = likelihood_sum(ap, sc, sp, streams,
                                           r_pts, es, sg,
                                           st_sum, xyz, st_prob,
                                           exp_stream_weights, sum_exp_weights, exp_background_weight);

    get_stream_only_likelihood(st_sum, sp->number_stars, streams->number_streams);

    free(exp_stream_weights);
    free(st_prob);
    free(r_pts);
    free(st_sum);

    return likelihood_val;
}

