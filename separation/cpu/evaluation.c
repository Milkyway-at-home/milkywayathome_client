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

#include "milkyway.h"
#include "milkyway_priv.h"

#include "evaluation.h"
#include "evaluation_state.h"

#define stdev 0.6
#define xr (3.0 * stdev)
#define absm 4.2
#define SIGMA_LIMIT 0.0001

static const double sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };

STREAM_CONSTANTS* init_constants(ASTRONOMY_PARAMETERS* ap,
                                 const BACKGROUND_PARAMETERS* bgp,
                                 const STREAMS* streams)
{
    unsigned int i;
    vector lbr;

    STREAM_CONSTANTS* sc = malloc(sizeof(STREAM_CONSTANTS) * streams->number_streams);

    ap->alpha = bgp->parameters[0];
    ap->q     = bgp->parameters[1];
    ap->r0    = bgp->parameters[2];
    ap->delta = bgp->parameters[3];

    if (ap->aux_bg_profile == 0)
    {
        ap->bg_a = 0;
        ap->bg_b = 0;
        ap->bg_c = 0;
    }
    else if (ap->aux_bg_profile == 1)
    {
        ap->bg_a = bgp->parameters[4];
        ap->bg_b = bgp->parameters[5];
        ap->bg_c = bgp->parameters[6];
    }
    else
    {
        fprintf(stderr, "Error: aux_bg_profile invalid");
    }

    ap->coeff = 1.0 / (stdev * sqrt(2.0 * pi));
    ap->alpha_delta3 = 3.0 - ap->alpha + ap->delta;

    for (i = 0; i < streams->number_streams; i++)
    {
        double stream_sigma = streams->parameters[i].stream_parameters[4];
        sc[i].large_sigma = (stream_sigma > SIGMA_LIMIT || stream_sigma < -SIGMA_LIMIT);
        sc[i].sigma_sq2 = 2.0 * sqr(stream_sigma);

        if (ap->sgr_coordinates == 0)
            ap->sgr_conversion = (SGRConversion) gc2lb;
        else if (ap->sgr_coordinates == 1)
        {
            fprintf(stderr, "gc2sgr probably broken right now, so refusing to run\n");
            ap->sgr_conversion = (SGRConversion) gc2sgr;
            mw_finish(EXIT_FAILURE);
        }
        else
        {
            fprintf(stderr, "Error: sgr_coordinates not valid");
            mw_finish(EXIT_FAILURE);
        }

        ap->sgr_conversion(ap->wedge,
                           streams->parameters[i].stream_parameters[0],
                           0,
                           &L(lbr),
                           &B(lbr));

        R(lbr) = streams->parameters[i].stream_parameters[1];
        lbr2xyz(lbr, sc[i].c);

        X(sc[i].a) =   sin(streams->parameters[i].stream_parameters[2])
                            * cos(streams->parameters[i].stream_parameters[3]);

        Y(sc[i].a) =   sin(streams->parameters[i].stream_parameters[2])
                            * sin(streams->parameters[i].stream_parameters[3]);

        Z(sc[i].a) = cos(streams->parameters[i].stream_parameters[2]);
    }

    return sc;
}

static void get_stream_gauss(const unsigned int convolve, STREAM_GAUSS* sg)
{
    unsigned int i;
    double* qgaus_X = malloc(sizeof(double) * convolve);

    sg->qgaus_W = (double*) malloc(sizeof(double) * convolve);
    sg->dx      = (double*) malloc(sizeof(double) * convolve);

    gaussLegendre(-1.0, 1.0, qgaus_X, sg->qgaus_W, convolve);

    for (i = 0; i < convolve; i++)
        sg->dx[i] = 3.0 * stdev * qgaus_X[i];

    free(qgaus_X);
}

static double set_prob_consts(const ASTRONOMY_PARAMETERS* ap,
                              const STREAM_GAUSS* sg,
                              const unsigned int n_convolve,
                              const double coords,
                              R_POINTS* rss)
{
    double g, exponent, r3, N;
    double reff_xr_rp3;
    unsigned int i;

    //R2MAG
    const double gPrime = 5.0 * (log10(coords * 1000.0) - 1.0) + absm;

    //REFF
    const double exp_result = exp(sigmoid_curve_params[1] * (gPrime - sigmoid_curve_params[2]));
    const double reff_value = sigmoid_curve_params[0] / (exp_result + 1.0);
    const double rPrime3 = cube(coords);

    for (i = 0; i < n_convolve; i++)
    {
        g = gPrime + sg->dx[i];

        //MAG2R
        rss[i].r_in_mag = g;
        rss[i].r_in_mag2 = sqr(g);
        rss[i].r_point = pow(10.0, (g - absm) / 5.0 + 1.0) / 1000.0;

        r3 = cube(rss[i].r_point);
        exponent = sqr(g - gPrime) / (2.0 * sqr(stdev));
        N = ap->coeff * exp(-exponent);
        rss[i].qw_r3_N = sg->qgaus_W[i] * r3 * N;
    }

    reff_xr_rp3 = reff_value * xr / rPrime3;
    return reff_xr_rp3;
}

/* FIXME: I don't know what these do enough to name it properly */
inline static double sub_bg_probability1(const ASTRONOMY_PARAMETERS* ap,
                                         const R_POINTS* rss,
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
        Z(xyz[i]) = rss[i].r_point * bsin;
        zp = rss[i].r_point * bcos;
        X(xyz[i]) = zp * lcos - sun_r0;
        Y(xyz[i]) = zp * lsin;

        rg = sqrt(sqr(X(xyz[i])) + sqr(Y(xyz[i])) + sqr(Z(xyz[i])) / sqr(ap->q));
        rs = rg + ap->r0;

        h_prob = rss[i].qw_r3_N / (rg * cube(rs));

        //the hernquist profile includes a quadratic term in g
        if (aux_bg_profile)
        {
            aux_prob = rss[i].qw_r3_N * (  ap->bg_a * rss[i].r_in_mag2
                                         + ap->bg_b * rss[i].r_in_mag
                                         + ap->bg_c );
            h_prob += aux_prob;
        }

        bg_prob += h_prob;
    }

    return bg_prob;
}

inline static double sub_bg_probability2(const ASTRONOMY_PARAMETERS* ap,
                                         const R_POINTS* rss,
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
        Z(xyz[i]) = rss[i].r_point * bsin;
        zp = rss[i].r_point * bcos;
        X(xyz[i]) = zp * lcos - sun_r0;
        Y(xyz[i]) = zp * lsin;

        rg = sqrt(sqr(X(xyz[i])) + sqr(Y(xyz[i])) + sqr(Z(xyz[i])) / sqr(ap->q));

        bg_prob += rss[i].qw_r3_N / (pow(rg, ap->alpha) * pow(rg + ap->r0, ap->alpha_delta3));
    }

    return bg_prob;
}

inline static double bg_probability(const ASTRONOMY_PARAMETERS* ap,
                                    const R_POINTS* rss,
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
            bg_prob = sub_bg_probability1(ap, rss, ap->convolve, ap->aux_bg_profile, integral_point, xyz);
        else
            bg_prob = sub_bg_probability2(ap, rss, ap->convolve, integral_point, xyz);

        bg_prob *= reff_xr_rp3;
    }

    return bg_prob;
}

/* FIXME: Better name? */
inline static double probabilities_convolve(const STREAM_CONSTANTS* sc,
                                            const R_POINTS* rss,
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

        st_prob += rss[i].qw_r3_N * exp(-xyz_norm / sc->sigma_sq2);
    }

    return st_prob;
}

inline static void probabilities(const ASTRONOMY_PARAMETERS* ap,
                                 const STREAM_CONSTANTS* sc,
                                 const R_POINTS* rss,
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
            st_prob = V * reff_xr_rp3 * probabilities_convolve(&sc[i], rss, ap->convolve, xyz);
        else
            st_prob = 0.0;

        KAHAN_ADD(probs[i].st_prob_int, st_prob, probs[i].st_prob_int_c);
    }
}

inline static void likelihood_probabilities(const ASTRONOMY_PARAMETERS* ap,
                                            const STREAM_CONSTANTS* sc,
                                            const R_POINTS* rss,
                                            const double reff_xr_rp3,
                                            vector* const xyz,
                                            double* probs)
{
    unsigned int i;

    for (i = 0; i < ap->number_streams; ++i)
    {
        if (sc[i].large_sigma)
            probs[i] = reff_xr_rp3 * probabilities_convolve(&sc[i], rss, ap->convolve, xyz);
        else
            probs[i] = 0.0;
    }
}

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


static void prepare_nu_constants(NU_CONSTANTS* nu_st,
                                 const unsigned int nu_steps,
                                 double nu_step_size,
                                 double nu_min)
{
    unsigned int i;
    double tmp1, tmp2;

    for (i = 0; i < nu_steps; i++)
    {
        nu_st[i].nu = nu_min + (i * nu_step_size);

        tmp1 = d2r(90.0 - nu_st[i].nu - nu_step_size);
        tmp2 = d2r(90.0 - nu_st[i].nu);

        nu_st[i].id = cos(tmp1) - cos(tmp2);
        nu_st[i].nu += 0.5 * nu_step_size;
    }
}

static R_CONSTANTS* prepare_r_constants(const ASTRONOMY_PARAMETERS* ap,
                                        const STREAM_GAUSS* sg,
                                        const unsigned int n_convolve,
                                        const unsigned int r_steps,
                                        const double r_min,
                                        const double r_step_size,
                                        const double mu_step_size,
                                        R_POINTS* rss)
{
    unsigned int i;
    double r, next_r, rPrime;
    R_CONSTANTS* r_step_consts = malloc(sizeof(R_CONSTANTS) * r_steps);

//vickej2_kpc edits to make volumes even in kpc rather than g
//vickej2_kpc        double log_r, r, next_r, rPrime;

  #ifdef USE_KPC
    const double r_max           = r_min + r_step_size * r_steps;
    const double r_min_kpc       = pow(10.0, ((r_min - 14.2) / 5.0));
    const double r_max_kpc       = pow(10.0, ((r_max - 14.2) / 5.0));
    const double r_step_size_kpc = (r_max_kpc - r_min_kpc) / r_steps;
  #endif

    for (i = 0; i < r_steps; i++)
    {
      #ifdef USE_KPC
        r = r_min_kpc + (i * r_step_size_kpc);
        next_r = r + r_step_size_kpc;
      #else
        double log_r = r_min + (i * r_step_size);
        r = pow(10.0, (log_r - 14.2) / 5.0);
        next_r = pow(10.0, (log_r + r_step_size - 14.2) / 5.0);
      #endif

        r_step_consts[i].irv = d2r(((cube(next_r) - cube(r)) / 3.0) * mu_step_size);
        rPrime = (next_r + r) / 2.0;

        r_step_consts[i].reff_xr_rp3 = set_prob_consts(ap, sg, n_convolve, rPrime, &rss[i * n_convolve]);
    }

    return r_step_consts;
}

static void prepare_integral_constants(const ASTRONOMY_PARAMETERS* ap,
                                       const STREAM_GAUSS* sg,
                                       const INTEGRAL_AREA* ia,
                                       INTEGRAL_CONSTANTS* ic)
{

    /* 2D block, ia->r_steps = rows, ap->convolve = columns */
    ic->rss = malloc(sizeof(R_POINTS) * ia->r_steps * ap->convolve);
    ic->r_step_consts = prepare_r_constants(ap,
                                            sg,
                                            ap->convolve,
                                            ia->r_steps,
                                            ia->r_min,
                                            ia->r_step_size,
                                            ia->mu_step_size,
                                            ic->rss);

    ic->nu_st = malloc(sizeof(NU_CONSTANTS) * ia->nu_steps);
    prepare_nu_constants(ic->nu_st, ia->nu_steps, ia->nu_step_size, ia->nu_min);
}

static void free_integral_constants(INTEGRAL_CONSTANTS* ic)
{
    free(ic->r_step_consts);
    free(ic->rss);
    free(ic->nu_st);
}

/* Sum over r steps using Kahan summation */
inline static BG_PROB r_sum(const ASTRONOMY_PARAMETERS* ap,
                            const STREAM_CONSTANTS* sc,
                            const INTEGRAL_CONSTANTS* ic,
                            const unsigned int r_steps,
                            vector* xyz,
                            ST_PROBS* probs,
                            const vector integral_point,
                            const unsigned int nu_step_current)

{
    unsigned int r_step_current;
    double V;
    double bg_prob;
    BG_PROB bg_prob_int = ZERO_BG_PROB; /* for Kahan summation */

    for (r_step_current = 0; r_step_current < r_steps; ++r_step_current)
    {
        bg_prob = bg_probability(ap,
                                 &ic->rss[r_step_current * ap->convolve],
                                 ic->r_step_consts[r_step_current].reff_xr_rp3,
                                 integral_point,
                                 xyz);

        V = ic->r_step_consts[r_step_current].irv * ic->nu_st[nu_step_current].id;
        bg_prob *= V;

        KAHAN_ADD(bg_prob_int.bg_int, bg_prob, bg_prob_int.correction);

        probabilities(ap,
                      sc,
                      &ic->rss[r_step_current * ap->convolve],
                      ic->r_step_consts[r_step_current].reff_xr_rp3,
                      V,
                      xyz,
                      probs);
    }

    return bg_prob_int;
}

inline static void nu_sum(const ASTRONOMY_PARAMETERS* ap,
                          const STREAM_CONSTANTS* sc,
                          const INTEGRAL_AREA* ia,
                          const INTEGRAL_CONSTANTS* ic,
                          EVALUATION_STATE* es,
                          vector* xyz,
                          ST_PROBS* probs,
                          const unsigned int mu_step_current)
{
    vector integral_point;
    BG_PROB r_result;

    const double mu = ia->mu_min + (mu_step_current * ia->mu_step_size);
    const unsigned int nu_steps = ia->nu_steps;
    const unsigned int r_steps = ia->r_steps;

    for ( ; es->nu_step < nu_steps; es->nu_step++)
    {
        do_boinc_checkpoint(ap, es);

        ap->sgr_conversion(ap->wedge,
                           mu + 0.5 * ia->mu_step_size,
                           ic->nu_st[es->nu_step].nu,
                           &L(integral_point),
                           &B(integral_point));

        r_result = r_sum(ap,
                         sc,
                         ic,
                         r_steps,
                         xyz,
                         probs,
                         integral_point,
                         es->nu_step);

        INCADD_BG_PROB(es->nu_acc, r_result);
    }

    es->nu_step = 0;
}

/* returns background integral */
static double integrate(const ASTRONOMY_PARAMETERS* ap,
                        const STREAM_CONSTANTS* sc,
                        const INTEGRAL_CONSTANTS* ic,
                        const INTEGRAL_AREA* ia,
                        vector* xyz,
                        EVALUATION_STATE* es,
                        ST_PROBS* probs)
{
    const unsigned int mu_steps = ia->mu_steps;

    for ( ; es->mu_step < mu_steps; es->mu_step++)
    {
        nu_sum(ap, sc, ia, ic, es, xyz, probs, es->mu_step);
        INCADD_BG_PROB(es->mu_acc, es->nu_acc);
        CLEAR_BG_PROB(es->nu_acc)
    }

    es->mu_step = 0;

    return es->mu_acc.bg_int + es->mu_acc.correction;
}

static void print_stream_integrals(EVALUATION_STATE* es, const unsigned int number_streams)
{
    unsigned int i;
    fprintf(stderr, "<background_integral> %.20lf </background_integral>\n", es->background_integral);
    fprintf(stderr, "<stream_integrals>");
    for (i = 0; i < number_streams; i++)
        fprintf(stderr, " %.20lf", es->stream_integrals[i]);
    fprintf(stderr, " </stream_integrals>\n");
}

static void final_stream_integrals(EVALUATION_STATE* es,
                                   const unsigned int number_streams,
                                   const unsigned int number_integrals)
{
    unsigned int i, j;

    es->background_integral = es->integrals[0].background_integral;
    for (i = 0; i < number_streams; ++i)
        es->stream_integrals[i] = es->integrals[0].stream_integrals[i];

    for (i = 1; i < number_integrals; ++i)
    {
        es->background_integral -= es->integrals[i].background_integral;
        for (j = 0; j < number_streams; j++)
            es->stream_integrals[j] -= es->integrals[i].stream_integrals[j];
    }

}

inline static void calculate_stream_integrals(const ST_PROBS* probs,
                                              double* stream_integrals,
                                              const unsigned int number_streams)
{
    unsigned int i;

    for (i = 0; i < number_streams; ++i)
        stream_integrals[i] = probs[i].st_prob_int + probs[i].st_prob_int_c;
}

static void calculate_integrals(const ASTRONOMY_PARAMETERS* ap,
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

        integral->background_integral = integrate(ap, sc, &ic, ia, xyz, es, integral->probs);
        calculate_stream_integrals(integral->probs, integral->stream_integrals, ap->number_streams);

        free_integral_constants(&ic);
        CLEAR_BG_PROB(es->mu_acc);
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

static double likelihood(const ASTRONOMY_PARAMETERS* ap,
                         const STREAM_CONSTANTS* sc,
                         const STREAMS* streams,
                         EVALUATION_STATE* es,
                         STREAM_GAUSS* sg,
                         vector* xyz,
                         const STAR_POINTS* sp)
{
    double star_prob;
    double bg_prob, bg, reff_xr_rp3;
    double exp_background_weight, sum_exp_weights;

    PROB_SUM prob = ZERO_PROB_SUM;
    PROB_SUM bg_only = ZERO_PROB_SUM;

    unsigned int current_star_point;
    unsigned int num_zero = 0;
    unsigned int bad_jacobians = 0;

    double* st_prob = malloc(sizeof(double) * streams->number_streams);

    R_POINTS* rss = malloc(sizeof(R_POINTS) * ap->convolve);
    ST_SUM* st_sum = calloc(sizeof(ST_SUM), streams->number_streams);

    double* exp_stream_weights = malloc(sizeof(double) * streams->number_streams);
    exp_background_weight = exp(ap->background_weight);
    sum_exp_weights = get_exp_stream_weights(exp_stream_weights, streams, exp_background_weight);


    for (current_star_point = 0; current_star_point < sp->number_stars; ++current_star_point)
    {
        reff_xr_rp3 = set_prob_consts(ap, &sg[0], ap->convolve, ZN(sp, current_star_point), rss);

        bg_prob = bg_probability(ap, rss,
                                 reff_xr_rp3, &VN(sp, current_star_point), xyz);

        bg = (bg_prob / es->background_integral) * exp_background_weight;

        likelihood_probabilities(ap, sc, rss, reff_xr_rp3, xyz, st_prob);

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

    get_stream_only_likelihood(st_sum, sp->number_stars, streams->number_streams);

    free(exp_stream_weights);
    free(st_prob);
    free(rss);
    free(st_sum);

    /*  log10(x * 0.001) = log10(x) - 3.0 */
    return (prob.sum / (sp->number_stars - bad_jacobians)) - 3.0;
}

static void free_stream_gauss(STREAM_GAUSS* sg)
{
    free(sg->dx);
    free(sg->qgaus_W);
}

double cpu_evaluate(const ASTRONOMY_PARAMETERS* ap,
                    const STAR_POINTS* sp,
                    const STREAMS* streams,
                    const STREAM_CONSTANTS* sc)
{
    double likelihood_val;
    EVALUATION_STATE es = EMPTY_EVALUATION_STATE;
    STREAM_GAUSS sg;

    initialize_state(ap, &es);
    get_stream_gauss(ap->convolve, &sg);

#if BOINC_APPLICATION
    if (boinc_file_exists(CHECKPOINT_FILE))
    {
        fprintf(stderr, "Checkpoint exists. Attempting to resume from it\n");

        if (read_checkpoint(&es))
        {
            fprintf(stderr, "Reading checkpoint failed\n");
            boinc_delete_file(CHECKPOINT_FILE);
            mw_finish(EXIT_FAILURE);
        }
    }
#endif

    vector* xyz = malloc(sizeof(vector) * ap->convolve);

    calculate_integrals(ap, sc, &sg, &es, xyz);

    /* FIXME: Force a checkpoint, don't wait for time. */
    /* Final checkpoint */
    do_boinc_checkpoint(ap, &es);

    final_stream_integrals(&es, ap->number_streams, ap->number_integrals);
    print_stream_integrals(&es, ap->number_streams);

    likelihood_val = likelihood(ap, sc, streams, &es, &sg, xyz, sp);

    free(xyz);
    free_evaluation_state(&es);
    free_stream_gauss(&sg);

  #if BOINC_APPLICATION
    boinc_delete_file(CHECKPOINT_FILE);
  #endif

    return likelihood_val;
}


