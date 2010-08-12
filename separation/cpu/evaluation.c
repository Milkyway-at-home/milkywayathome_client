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

#define CHECKPOINT_FILE "astronomy_checkpoint"

#include "milkyway.h"
#include "milkyway_priv.h"

#define stdev 0.6
#define xr (3.0 * stdev)
#define absm 4.2

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
        sc[i].stream_sigma = streams->parameters[i].stream_parameters[4];
        sc[i].stream_sigma_sq2 = 2.0 * sqr(sc[i].stream_sigma);

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
        lbr2xyz(lbr, sc[i].stream_c);

        X(sc[i].stream_a) =   sin(streams->parameters[i].stream_parameters[2])
                            * cos(streams->parameters[i].stream_parameters[3]);

        Y(sc[i].stream_a) =   sin(streams->parameters[i].stream_parameters[2])
                            * sin(streams->parameters[i].stream_parameters[3]);

        Z(sc[i].stream_a) = cos(streams->parameters[i].stream_parameters[2]);
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
                              R_STEP_STATE* rss)
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
        rss[i].r_in_mag2 = g * g;
        rss[i].r_point = pow(10.0, (g - absm) / 5.0 + 1.0) / 1000.0;

        r3 = rss[i].r_point * rss[i].r_point * rss[i].r_point;
        exponent = sqr(g - gPrime) / (2.0 * sqr(stdev));
        N = ap->coeff * exp(-exponent);
        rss[i].qw_r3_N = sg->qgaus_W[i] * r3 * N;
    }

    reff_xr_rp3 = reff_value * xr / rPrime3;
    return reff_xr_rp3;
}

/* FIXME: I don't know what these do enough to name it properly */
inline static double sub_bg_probability1(const ASTRONOMY_PARAMETERS* ap,
                                         const R_STEP_STATE* rss,
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
        xyz[i][2] = rss[i].r_point * bsin;
        zp = rss[i].r_point * bcos;
        xyz[i][0] = zp * lcos - sun_r0;
        xyz[i][1] = zp * lsin;

        rg = sqrt( sqr(xyz[i][0]) + sqr(xyz[i][1]) + sqr(xyz[i][2]) / sqr(ap->q));
        rs = rg + ap->r0;

        //the hernquist profile includes a quadratic term in g
        if (aux_bg_profile == 1)
        {
            h_prob = rss[i].qw_r3_N / (rg * cube(rs));
            aux_prob = rss[i].qw_r3_N * (ap->bg_a * rss[i].r_in_mag2 + ap->bg_b * rss[i].r_in_mag + ap->bg_c );
            bg_prob += h_prob + aux_prob;
        }
        else if (aux_bg_profile == 0)
        {
            bg_prob += rss[i].qw_r3_N / (rg * cube(rs));
        }
        else
        {
            fprintf(stderr, "Error: aux_bg_profile invalid");
        }
    }

    return bg_prob;
}

inline static double sub_bg_probability2(const ASTRONOMY_PARAMETERS* ap,
                                         const R_STEP_STATE* rss,
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
        xyz[i][2] = rss[i].r_point * bsin;
        zp = rss[i].r_point * bcos;
        xyz[i][0] = zp * lcos - sun_r0;
        xyz[i][1] = zp * lsin;

        rg = sqrt( sqr(xyz[i][0]) + sqr(xyz[i][1]) + sqr(xyz[i][2]) / sqr(ap->q));

        bg_prob += rss[i].qw_r3_N / (pow(rg, ap->alpha) * pow(rg + ap->r0, ap->alpha_delta3));
    }

    return bg_prob;
}

inline static double bg_probability(const ASTRONOMY_PARAMETERS* ap,
                                    const R_STEP_STATE* rss,
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
inline static void probabilities_convolve(const STREAM_CONSTANTS* sc,
                                          const R_STEP_STATE* rss,
                                          ST_PROBS* probs,
                                          const unsigned int convolve,
                                          const double reff_xr_rp3,
                                          vector* const xyz)
{
    unsigned int i;
    double dotted, xyz_norm;
    vector xyzs;

    for (i = 0; i < convolve; i++)
    {
        X(xyzs) = X(xyz[i]) - X(sc->stream_c);
        Y(xyzs) = Y(xyz[i]) - Y(sc->stream_c);
        Z(xyzs) = Z(xyz[i]) - Z(sc->stream_c);

        dotted = X(sc->stream_a) * X(xyzs)
               + Y(sc->stream_a) * Y(xyzs)
               + Z(sc->stream_a) * Z(xyzs);

        X(xyzs) = X(xyzs) - dotted * X(sc->stream_a);
        Y(xyzs) = Y(xyzs) - dotted * Y(sc->stream_a);
        Z(xyzs) = Z(xyzs) - dotted * Z(sc->stream_a);

        xyz_norm =  sqr(X(xyzs)) + sqr(Y(xyzs)) + sqr(Z(xyzs));

        probs->st_prob += rss[i].qw_r3_N * exp(-xyz_norm / sc->stream_sigma_sq2);
    }

    probs->st_prob *= reff_xr_rp3;
}

static void probabilities(const ASTRONOMY_PARAMETERS* ap,
                          const STREAM_CONSTANTS* sc,
                          const R_STEP_STATE* rss,
                          const double reff_xr_rp3,
                          vector* const xyz,
                          ST_PROBS* probs)
{
    unsigned int i;

    for (i = 0; i < ap->number_streams; i++)
    {
        probs[i].st_prob = 0.0;
        if (sc[i].stream_sigma > -0.0001 && sc[i].stream_sigma < 0.0001)
            continue;

        probabilities_convolve(&sc[i], rss, &probs[i], ap->convolve, reff_xr_rp3, xyz);
    }
}

inline static double progress(const ASTRONOMY_PARAMETERS* ap,
                              const EVALUATION_STATE* es,
                              unsigned int mu_step_current,
                              unsigned int nu_step_current)
{
    unsigned int i;
    const INTEGRAL_AREA* ia;

    unsigned int current_probs;
    unsigned int total_calc_probs = 0;
    unsigned int current_calc_probs = 0;

    for (i = 0; i < es->number_integrals; i++)
    {
        ia = &ap->integral[i];

        current_probs = ia->r_steps * ia->mu_steps * ia->nu_steps;
        total_calc_probs += current_probs;
        if (i < es->current_integral)
        {
            current_calc_probs += current_probs;
        }
        else if (i == es->current_integral)
        {
            /* When checkpointing is done, ia->r_step would always be 0 */
            current_calc_probs +=   (mu_step_current * ia->nu_steps * ia->r_steps)
                                  + (nu_step_current * ia->r_steps); /* + ia->r_step */
        }
    }

    total_calc_probs += es->total_stars;
    current_calc_probs += es->current_star_point;

    return (double)current_calc_probs / (double)total_calc_probs;
}

#if BOINC_APPLICATION
inline static void do_boinc_checkpoint(const ASTRONOMY_PARAMETERS* ap,
                                       EVALUATION_STATE* es,
                                       unsigned int mu_step_current,
                                       unsigned int nu_step_current)
{
    double frac;

    if (boinc_time_to_checkpoint())
    {
        /* FIXME: Make checkpointing make more sense, then we won't need this */
        es->integrals[es->current_integral].mu_step = mu_step_current;
        es->integrals[es->current_integral].nu_step = nu_step_current;

        int retval = write_checkpoint(es);
        if (retval)
        {
            fprintf(stderr, "APP: astronomy checkpoint failed %d\n", retval);
            return;
        }
        boinc_checkpoint_completed();
    }

    frac = progress(ap, es, mu_step_current, nu_step_current);
    //printf("progress: %.10f\n", frac);
    boinc_fraction_done(frac);
}

#else

inline static void do_boinc_checkpoint(const ASTRONOMY_PARAMETERS* ap,
                                       EVALUATION_STATE* es,
                                       unsigned int mu_step_current,
                                       unsigned int nu_step_current)
{
}

#endif /* BOINC_APPLICATION */


static void prepare_nu_constants(NU_STATE* nu_st,
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

static R_STEP_CONSTANTS* prepare_r_constants(const ASTRONOMY_PARAMETERS* ap,
                                             const STREAM_GAUSS* sg,
                                             const unsigned int n_convolve,
                                             const unsigned int r_steps,
                                             const double r_min,
                                             const double r_step_size,
                                             const double mu_step_size,
                                             R_STEP_STATE* rss)
{
    unsigned int i;
    double r, next_r, rPrime;
    R_STEP_CONSTANTS* r_step_consts = malloc(sizeof(R_STEP_CONSTANTS) * r_steps);


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

static void prepare_integral_state(const ASTRONOMY_PARAMETERS* ap,
                                   const STREAM_GAUSS* sg,
                                   INTEGRAL_AREA* ia,
                                   INTEGRAL_STATE* st)
{

    st->probs = (ST_PROBS*) malloc(sizeof(ST_PROBS) * ap->number_streams);

    /* 2D block, ia->r_steps = rows, ap->convolve = columns */
    st->rss = malloc(sizeof(R_STEP_STATE) * ia->r_steps * ap->convolve);
    st->r_step_consts = prepare_r_constants(ap,
                                            sg,
                                            ap->convolve,
                                            ia->r_steps,
                                            ia->r_min,
                                            ia->r_step_size,
                                            ia->mu_step_size,
                                            st->rss);

    st->nu_st = malloc(sizeof(NU_STATE) * ia->nu_steps);
    prepare_nu_constants(st->nu_st, ia->nu_steps, ia->nu_step_size, ia->nu_min);
}

static void free_integral_state(INTEGRAL_STATE* st)
{
    free(st->r_step_consts);
    free(st->probs);
    free(st->rss);
    free(st->nu_st);
}

inline static void update_probs(ST_PROBS* probs, const unsigned int n_streams, const double V)
{
    unsigned int i;
    double tmp;

    for (i = 0; i < n_streams; ++i)
    {
        probs[i].st_prob *= V;
        tmp = probs[i].st_prob_int;
        probs[i].st_prob_int += probs[i].st_prob;
        probs[i].st_prob_int_c += probs[i].st_prob - (probs[i].st_prob_int - tmp);
    }
}

/* Sum over r steps using Kahan summation */
inline static BG_PROB r_sum(const ASTRONOMY_PARAMETERS* ap,
                            const STREAM_CONSTANTS* sc,
                            const unsigned int r_steps,
                            INTEGRAL_STATE* st,
                            vector* xyz,
                            const vector integral_point,
                            const unsigned int nu_step_current)

{
    unsigned int r_step_current;
    double V, tmp;
    double bg_prob;
    BG_PROB bg_prob_int = { 0.0, 0.0 }; /* for Kahan summation */

    const unsigned int n_streams = ap->number_streams;

    for (r_step_current = 0; r_step_current < r_steps; ++r_step_current)
    {
        V = st->r_step_consts[r_step_current].irv * st->nu_st[nu_step_current].id;

        bg_prob = bg_probability(ap,
                                 &st->rss[r_step_current * ap->convolve],
                                 st->r_step_consts[r_step_current].reff_xr_rp3,
                                 integral_point,
                                 xyz);

        probabilities(ap,
                      sc,
                      &st->rss[r_step_current * ap->convolve],
                      st->r_step_consts[r_step_current].reff_xr_rp3,
                      xyz,
                      st->probs);

        bg_prob *= V;

        tmp = bg_prob_int.bg_int;
        bg_prob_int.bg_int += bg_prob;
        bg_prob_int.correction += bg_prob - (bg_prob_int.bg_int - tmp);

        update_probs(st->probs, n_streams, V);
    }

    return bg_prob_int;
}

inline static void apply_correction(const unsigned int number_streams,
                                    INTEGRAL* integral,
                                    ST_PROBS* probs,
                                    BG_PROB bg_prob_int)
{
    unsigned int i;
    integral->background_integral = bg_prob_int.bg_int + bg_prob_int.correction;
    for (i = 0; i < number_streams; i++)
        integral->stream_integrals[i] = probs[i].st_prob_int + probs[i].st_prob_int_c;
}

inline static BG_PROB nu_sum(const ASTRONOMY_PARAMETERS* ap,
                             const STREAM_CONSTANTS* sc,
                             INTEGRAL_AREA* ia,
                             EVALUATION_STATE* es,
                             INTEGRAL_STATE* st,
                             vector* xyz,
                             const unsigned int mu_step_current)
{
    unsigned int nu_step_current;
    vector integral_point;
    BG_PROB r_result;
    BG_PROB bg_prob_int = { 0.0, 0.0 };

    INTEGRAL* integral = &es->integrals[es->current_integral];

    const double mu = ia->mu_min + (mu_step_current * ia->mu_step_size);
    const unsigned int nu_steps = ia->nu_steps;
    const unsigned int r_steps = ia->r_steps;

    for (nu_step_current = 0; nu_step_current < nu_steps; ++nu_step_current)
    {
        apply_correction(ap->number_streams, integral, st->probs, bg_prob_int);

        do_boinc_checkpoint(ap, es, mu_step_current, nu_step_current);

        ap->sgr_conversion(ap->wedge,
                           mu + 0.5 * ia->mu_step_size,
                           st->nu_st[nu_step_current].nu,
                           &L(integral_point),
                           &B(integral_point));

        r_result = r_sum(ap,
                         sc,
                         r_steps,
                         st,
                         xyz,
                         integral_point,
                         nu_step_current);

        bg_prob_int.bg_int += r_result.bg_int;
        bg_prob_int.correction += r_result.correction;
    }

    return bg_prob_int;
}

inline static void init_st_probs(ST_PROBS* probs,
                                 const double* stream_integrals,
                                 const unsigned int n_streams)
{
    unsigned int i;
    for (i = 0; i < n_streams; i++)
    {
        probs[i].st_prob = 0.0;
        probs[i].st_prob_int = stream_integrals[i];
        probs[i].st_prob_int_c = 0.0;
    }
}

static void integrate(const ASTRONOMY_PARAMETERS* ap,
                      const STREAM_CONSTANTS* sc,
                      vector* xyz,
                      EVALUATION_STATE* es,
                      INTEGRAL_STATE* st)
{
    unsigned int mu_step_current;
    BG_PROB nu_result;

    INTEGRAL_AREA* ia = &ap->integral[es->current_integral];
    INTEGRAL* integral = &es->integrals[es->current_integral];
    BG_PROB bg_prob_int = { integral->background_integral, 0.0 };

    init_st_probs(st->probs, integral->stream_integrals, ap->number_streams);

    const unsigned int mu_steps = ia->mu_steps;

    for (mu_step_current = 0; mu_step_current < mu_steps; mu_step_current++)
    {
        nu_result = nu_sum(ap, sc, ia, es, st, xyz, mu_step_current);

        bg_prob_int.bg_int += nu_result.bg_int;
        bg_prob_int.correction += nu_result.correction;
    }

    apply_correction(ap->number_streams, integral, st->probs, bg_prob_int);
}

static void print_stream_integrals(const ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es)
{
    unsigned int i;
    fprintf(stderr, "<background_integral> %.20lf </background_integral>\n", es->background_integral);
    fprintf(stderr, "<stream_integrals>");
    for (i = 0; i < ap->number_streams; i++)
        fprintf(stderr, " %.20lf", es->stream_integrals[i]);
    fprintf(stderr, " </stream_integrals>\n");
}

static int calculate_integrals(const ASTRONOMY_PARAMETERS* ap,
                               const STREAM_CONSTANTS* sc,
                               const STREAM_GAUSS* sg,
                               EVALUATION_STATE* es,
                               vector* xyz)
{
    unsigned int i, j;
    INTEGRAL_STATE st;

  #if BOINC_APPLICATION
    read_checkpoint(es);
  #endif

    for (; es->current_integral < ap->number_integrals; es->current_integral++)
    {
        prepare_integral_state(ap, sg, &ap->integral[es->current_integral], &st);
        integrate(ap, sc, xyz, es, &st);
        free_integral_state(&st);
    }

    es->background_integral = es->integrals[0].background_integral;
    for (i = 0; i < ap->number_streams; i++)
        es->stream_integrals[i] = es->integrals[0].stream_integrals[i];

    for (i = 1; i < ap->number_integrals; i++)
    {
        es->background_integral -= es->integrals[i].background_integral;
        for (j = 0; j < ap->number_streams; j++)
            es->stream_integrals[j] -= es->integrals[i].stream_integrals[j];
    }

    print_stream_integrals(ap, es);

    return 0;
}

typedef struct
{
    double st_only_sum;
    double st_only_sum_c;
} ST_SUM;

/* Used in likelihood calculation */
inline static double stream_sum(const unsigned int number_streams,
                                EVALUATION_STATE* es,
                                ST_PROBS* st_prob,
                                ST_SUM* st_sum,
                                const double* exp_stream_weights,
                                const double sum_exp_weights,
                                double bg_only)
{
    unsigned int current_stream;
    double st_only, tmp;
    double star_prob = bg_only;

    for (current_stream = 0; current_stream < number_streams; current_stream++)
    {
        st_only = st_prob[current_stream].st_prob / es->stream_integrals[current_stream] * exp_stream_weights[current_stream];
        star_prob += st_only;

        if (st_only == 0.0)
            st_only = -238.0;
        else
            st_only = log10(st_only / sum_exp_weights);

        tmp = st_sum[current_stream].st_only_sum;
        st_sum[current_stream].st_only_sum += st_only;
        st_sum[current_stream].st_only_sum_c += st_only - (st_sum[current_stream].st_only_sum - tmp);
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
    double bg_prob;
    double prob_sum, prob_sum_c, temp;  // for Kahan summation
    double exp_background_weight, sum_exp_weights;
    double reff_xr_rp3;

    double bg_only, bg_only_sum, bg_only_sum_c;

    /* The correction terms aren't used here since this isn't the sum? */
    ST_PROBS* st_prob = (ST_PROBS*) malloc(sizeof(ST_PROBS) * streams->number_streams);
    R_STEP_STATE* rss = malloc(sizeof(R_STEP_STATE) * ap->convolve);
    ST_SUM* st_sum = calloc(sizeof(ST_SUM), streams->number_streams);

    double* exp_stream_weights = malloc(sizeof(double) * streams->number_streams);
    exp_background_weight = exp(ap->background_weight);
    sum_exp_weights = get_exp_stream_weights(exp_stream_weights, streams, exp_background_weight);

    do_boinc_checkpoint(ap, es, 0, 0); /* CHECKME: Steps? */

    prob_sum = 0.0;
    prob_sum_c = 0.0;

    bg_only_sum = 0.0;
    bg_only_sum_c = 0.0;

    for (; es->current_star_point < sp->number_stars; es->current_star_point++)
    {
        double star_prob;

        reff_xr_rp3 = set_prob_consts(ap, &sg[0], ap->convolve, ZN(sp, es->current_star_point), rss);

        bg_prob = bg_probability(ap, rss,
                                 reff_xr_rp3, &VN(sp, es->current_star_point), xyz);

        bg_only = (bg_prob / es->background_integral) * exp_background_weight;

        probabilities(ap, sc, rss, reff_xr_rp3, xyz, st_prob);

        star_prob = stream_sum(streams->number_streams,
                               es,
                               st_prob,
                               st_sum,
                               exp_stream_weights,
                               sum_exp_weights,
                               bg_only);

        if (star_prob != 0.0)
        {
            star_prob = log10(star_prob);
            temp = prob_sum;
            prob_sum += star_prob;
            prob_sum_c += star_prob - (prob_sum - temp);
        }
        else
        {
            es->num_zero++;
            prob_sum -= 238.0;
        }

        if (bg_only == 0.0)
            bg_only = -238.0;
        else
            bg_only = log10(bg_only / sum_exp_weights);

        temp = bg_only_sum;
        bg_only_sum += bg_only;
        bg_only_sum_c += bg_only - (bg_only_sum - temp);
    }
    es->prob_sum = prob_sum + prob_sum_c;
    bg_only_sum += bg_only_sum_c;
    bg_only_sum /= sp->number_stars;

    fprintf(stderr, "<background_only_likelihood> %.20lf </background_only_likelihood>\n", bg_only_sum - 3.0);

    get_stream_only_likelihood(st_sum, sp->number_stars, streams->number_streams);

    free(exp_stream_weights);
    free(st_prob);
    free(rss);
    free(st_sum);

    /*  log10(x * 0.001) = log10(x) - 3.0 */
    return (es->prob_sum / (sp->number_stars - es->bad_jacobians)) - 3.0;
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
    int retval;
    EVALUATION_STATE es = EMPTY_EVALUATION_STATE;
    STREAM_GAUSS sg;

    initialize_state(ap, sp, &es);
    get_stream_gauss(ap->convolve, &sg);

    reset_evaluation_state(&es);

    vector* xyz = malloc(sizeof(vector) * ap->convolve);

    retval = calculate_integrals(ap, sc, &sg, &es, xyz);
    if (retval)
    {
        fprintf(stderr, "APP: error calculating integrals: %d\n", retval);
        mw_finish(retval);
    }

    double likelihood_val = likelihood(ap, sc, streams, &es, &sg, xyz, sp);

    free(xyz);
    free_evaluation_state(&es);
    free_stream_gauss(&sg);

    return likelihood_val;
}


