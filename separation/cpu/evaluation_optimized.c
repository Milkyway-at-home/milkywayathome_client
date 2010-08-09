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

#include <math.h>
#include <time.h>

#include "milkyway.h"
#include "milkyway_priv.h"

#include "evaluation_optimized.h"
#include "parameters.h"
#include "probability.h"
#include "atSurveyGeometry.h"
#include "stCoords.h"
#include "star_points.h"
#include "numericalIntegration.h"
#include "evaluation_state.h"


#define stdev 0.6
#define xr 3.0 * stdev
#define lbr_r 8.5
#define absm 4.2

static const double sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };

STREAM_CONSTANTS* init_constants(ASTRONOMY_PARAMETERS* ap, STREAM_NUMS* sn)
{
    unsigned int i;
    STREAM_CONSTANTS* sc = malloc(sizeof(STREAM_CONSTANTS) * ap->number_streams);

    sn->alpha = ap->background_parameters[0];
    sn->q     = ap->background_parameters[1];
    sn->r0    = ap->background_parameters[2];
    sn->delta = ap->background_parameters[3];

    if (ap->aux_bg_profile == 0)
    {
        sn->bg_a = 0;
        sn->bg_b = 0;
        sn->bg_c = 0;
    }
    else if (ap->aux_bg_profile == 1)
    {
        sn->bg_a = ap->background_parameters[4];
        sn->bg_b = ap->background_parameters[5];
        sn->bg_c = ap->background_parameters[6];
    }
    else
    {
        fprintf(stderr, "Error: aux_bg_profile invalid");
    }

    sn->coeff = 1.0 / (stdev * sqrt(2.0 * pi));
    sn->alpha_delta3 = 3.0 - sn->alpha + sn->delta;

    for (i = 0; i < ap->number_streams; i++)
    {
        RA_DEC radec;
        double lamda, beta, l, b;
        vector lbr;

        sc[i].stream_sigma = ap->parameters[i].stream_parameters[4];
        sc[i].stream_sigma_sq2 = 2.0 * sc[i].stream_sigma * sc[i].stream_sigma;

        if (ap->sgr_coordinates == 0)
        {
            radec = atGCToEq(ap->parameters[i].stream_parameters[0],
                             0, wedge_incl(ap->wedge));
            atEqToGal(radec.ra, radec.dec, &l, &b);
        }
        else if (ap->sgr_coordinates == 1)
        {
            gcToSgr(ap->parameters[i].stream_parameters[0], 0, ap->wedge, &lamda, &beta);
            sgrToGal(lamda, beta, &l, &b);
        }
        else
        {
            fprintf(stderr, "Error: sgr_coordinates not valid");
        }

        L(lbr) = l;
        B(lbr) = b;
        R(lbr) = ap->parameters[i].stream_parameters[1];
        lbr2xyz(lbr, sc[i].stream_c);

        X(sc[i].stream_a) =   sin(ap->parameters[i].stream_parameters[2])
                            * cos(ap->parameters[i].stream_parameters[3]);

        Y(sc[i].stream_a) =   sin(ap->parameters[i].stream_parameters[2])
                            * sin(ap->parameters[i].stream_parameters[3]);

        Z(sc[i].stream_a) = cos(ap->parameters[i].stream_parameters[2]);
    }

    return sc;
}

static void get_stream_gauss(const ASTRONOMY_PARAMETERS* ap, STREAM_GAUSS* sg)
{
    unsigned int i;
    double* qgaus_X = malloc(sizeof(double) * ap->convolve);

    sg->qgaus_W = (double*) malloc(sizeof(double) * ap->convolve);
    sg->dx      = (double*) malloc(sizeof(double) * ap->convolve);

    gaussLegendre(-1.0, 1.0, qgaus_X, sg->qgaus_W, ap->convolve);

    for (i = 0; i < ap->convolve; i++)
        sg->dx[i] = 3.0 * stdev * qgaus_X[i];

    free(qgaus_X);
}

void free_constants(ASTRONOMY_PARAMETERS* ap)
{
}

static double set_probability_constants(const STREAM_NUMS* sn,
                                        const STREAM_GAUSS* sg,
                                        const unsigned int n_convolve,
                                        double coords,
                                        R_STEP_STATE* rss,
                                        unsigned int r_step)
{
    double gPrime, exp_result, g, exponent, r3, N, reff_value, rPrime3;
    double reff_xr_rp3;
    unsigned int i, idx;
    const unsigned int row_base = n_convolve * r_step;

    //R2MAG
    gPrime = 5.0 * (log10(coords * 1000.0) - 1.0) + absm;

    //REFF
    exp_result = exp(sigmoid_curve_params[1] * (gPrime - sigmoid_curve_params[2]));
    reff_value = sigmoid_curve_params[0] / (exp_result + 1);
    rPrime3 = coords * coords * coords;

    for (i = 0; i < n_convolve; i++)
    {
        idx = row_base + i;
        g = gPrime + sg->dx[i];

        //MAG2R
        rss[idx].r_in_mag = g;
        rss[idx].r_in_mag2 = g * g;
        rss[idx].r_point = pow(10.0, (g - absm) / 5.0 + 1.0) / 1000.0;

        r3 = rss[idx].r_point * rss[idx].r_point * rss[idx].r_point;
        exponent = (g - gPrime) * (g - gPrime) / (2.0 * stdev * stdev);
        N = sn->coeff * exp(-exponent);
        rss[idx].qw_r3_N = sg->qgaus_W[i] * r3 * N;
    }

    reff_xr_rp3 = reff_value * xr / rPrime3;
    return reff_xr_rp3;
}

static double calculate_bg_probability(const ASTRONOMY_PARAMETERS* ap,
                                       const STREAM_NUMS* sn,
                                       const R_STEP_STATE* rss,
                                       const unsigned int r_step_current,
                                       const double reff_xr_rp3,
                                       const vector integral_point,
                                       vector* xyz)
{
    double zp;
    double rg, rs;
    double h_prob, aux_prob;
    double bg_prob;
    unsigned int i, idx;

    const unsigned int row_base = r_step_current * ap->convolve;

    const double bsin = sin(integral_point[1] / deg);
    const double lsin = sin(integral_point[0] / deg);
    const double bcos = cos(integral_point[1] / deg);
    const double lcos = cos(integral_point[0] / deg);

    /* if q is 0, there is no probability */
    if (sn->q == 0)
    {
        bg_prob = -1.0;
    }
    else
    {
        bg_prob = 0.0;
        if (sn->alpha == 1 && sn->delta == 1)
        {
            for (i = 0; i < ap->convolve; i++)
            {
                /* Index into 2d array */
                idx = row_base + i;

                xyz[i][2] = rss[idx].r_point * bsin;
                zp = rss[idx].r_point * bcos;
                xyz[i][0] = zp * lcos - lbr_r;
                xyz[i][1] = zp * lsin;

                rg = sqrt( xyz[i][0] * xyz[i][0] + xyz[i][1] * xyz[i][1] + (xyz[i][2] * xyz[i][2])
                                   /
                            (sn->q * sn->q)
                         );
                rs = rg + sn->r0;

                //the hernquist profile includes a quadratic term in g
                if (ap->aux_bg_profile == 1)
                {
                    h_prob = rss[idx].qw_r3_N / (rg * rs * rs * rs);
                    aux_prob = rss[idx].qw_r3_N * (sn->bg_a * rss[idx].r_in_mag2 + sn->bg_b * rss[idx].r_in_mag + sn->bg_c );
                    bg_prob += h_prob + aux_prob;
                }
                else if (ap->aux_bg_profile == 0)
                {
                    bg_prob += rss[idx].qw_r3_N / (rg * rs * rs * rs);
                }
                else
                {
                    fprintf(stderr, "Error: aux_bg_profile invalid");
                }
            }
        }
        else
        {
            for (i = 0; i < ap->convolve; i++)
            {
                idx = row_base + i;
                xyz[i][2] = rss[idx].r_point * bsin;
                zp = rss[idx].r_point * bcos;
                xyz[i][0] = zp * lcos - lbr_r;
                xyz[i][1] = zp * lsin;

                rg = sqrt(xyz[i][0] * xyz[i][0] + xyz[i][1] * xyz[i][1] + (xyz[i][2] * xyz[i][2])
                          / (sn->q * sn->q));

                bg_prob += rss[idx].qw_r3_N / (pow(rg, sn->alpha) * pow(rg + sn->r0, sn->alpha_delta3));
            }
        }
        bg_prob *= reff_xr_rp3;
    }

    return bg_prob;
}

static void calculate_probabilities(const ASTRONOMY_PARAMETERS* ap,
                                    const STREAM_CONSTANTS* sc,
                                    const R_STEP_STATE* rss,
                                    const unsigned int r_step_current,
                                    const double reff_xr_rp3,
                                    vector* const xyz,
                                    ST_PROBS* probs)
{
    unsigned int i, j, idx;
    double dotted, xyz_norm;
    vector xyzs;
    const unsigned int row_base = r_step_current * ap->convolve;

    for (i = 0; i < ap->number_streams; i++)
    {
        probs[i].st_prob = 0.0;
        if (sc[i].stream_sigma > -0.0001 && sc[i].stream_sigma < 0.0001)
            continue;
        for (j = 0; j < ap->convolve; j++)
        {
            X(xyzs) = X(xyz[j]) - X(sc[i].stream_c);
            Y(xyzs) = Y(xyz[j]) - Y(sc[i].stream_c);
            Z(xyzs) = Z(xyz[j]) - Z(sc[i].stream_c);

            dotted = X(sc[i].stream_a) * X(xyzs)
                   + Y(sc[i].stream_a) * Y(xyzs)
                   + Z(sc[i].stream_a) * Z(xyzs);

            X(xyzs) = X(xyzs) - dotted * X(sc[i].stream_a);
            Y(xyzs) = Y(xyzs) - dotted * Y(sc[i].stream_a);
            Z(xyzs) = Z(xyzs) - dotted * Z(sc[i].stream_a);

            xyz_norm = X(xyzs) * X(xyzs)
                     + Y(xyzs) * Y(xyzs)
                     + Z(xyzs) * Z(xyzs);

            idx = row_base + j;
            probs[i].st_prob += rss[idx].qw_r3_N * exp(-xyz_norm / sc[i].stream_sigma_sq2);
        }
        probs[i].st_prob *= reff_xr_rp3;
    }
}

inline static double calculate_progress(const EVALUATION_STATE* es,
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
        ia = &es->integrals[i];

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
                                  + (nu_step_current * ia->r_steps); /* + ia->r_step*/
        }
    }

    total_calc_probs += es->total_stars;
    current_calc_probs += es->current_star_point;

    return (double)current_calc_probs / (double)total_calc_probs;
}


inline static void do_boinc_checkpoint(EVALUATION_STATE* es,
                                       unsigned int mu_step_current,
                                       unsigned int nu_step_current)
{
    double progress;

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

    progress = calculate_progress(es, mu_step_current, nu_step_current);
    //printf("progress: %.10f\n", progress);
    boinc_fraction_done(progress);
}


static void prepare_nu_constants(NU_STATE* nu_st,
                                 unsigned int nu_steps,
                                 double nu_step_size,
                                 double nu_min)
{
    unsigned int i;

    for (i = 0; i < nu_steps; i++)
    {
        nu_st[i].nus = nu_min + (i * nu_step_size);
        nu_st[i].ids = cos((90.0 - nu_st[i].nus - nu_step_size) / deg) - cos((90.0 - nu_st[i].nus) / deg);
        nu_st[i].nus += 0.5 * nu_step_size;
    }
}

static R_STEP_CONSTANTS* prepare_r_constants(const STREAM_NUMS* sn,
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

        r_step_consts[i].irv = (((next_r * next_r * next_r) - (r * r * r)) / 3.0) * mu_step_size / deg;
        rPrime = (next_r + r) / 2.0;

        r_step_consts[i].reff_xr_rp3 = set_probability_constants(sn,
                                                                 sg,
                                                                 n_convolve,
                                                                 rPrime,
                                                                 rss,
                                                                 i);
    }

    return r_step_consts;
}

static void prepare_integral_state(const ASTRONOMY_PARAMETERS* ap,
                                   const STREAM_NUMS* sn,
                                   const STREAM_GAUSS* sg,
                                   INTEGRAL_AREA* ia,
                                   INTEGRAL_STATE* st)
{

    st->probs = (ST_PROBS*) malloc(sizeof(ST_PROBS) * ap->number_streams);

    /* 2D block, ia->r_steps = rows, ap->convolve = columns */
    st->rss = malloc(sizeof(R_STEP_STATE) * ia->r_steps * ap->convolve);
    st->r_step_consts = prepare_r_constants(sn,
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

/* Sum over r steps using Kahan summation */
inline static BG_PROB r_sum(const ASTRONOMY_PARAMETERS* ap,
                            const STREAM_CONSTANTS* sc,
                            const STREAM_NUMS* sn,
                            const unsigned int r_steps,
                            INTEGRAL_STATE* st,
                            vector* xyz,
                            const vector integral_point,
                            const unsigned int nu_step_current)

{
    unsigned int i, r_step_current;
    double V, temp;
    double bg_prob;
    BG_PROB bg_prob_int = { 0.0, 0.0 }; /* for Kahan summation */

    for (r_step_current = 0; r_step_current < r_steps; ++r_step_current)
    {
        V = st->r_step_consts[r_step_current].irv * st->nu_st[nu_step_current].ids;

        bg_prob = calculate_bg_probability(ap,
                                           sn,
                                           st->rss,
                                           r_step_current,
                                           st->r_step_consts[r_step_current].reff_xr_rp3,
                                           integral_point,
                                           xyz);

        calculate_probabilities(ap,
                                sc,
                                st->rss,
                                r_step_current,
                                st->r_step_consts[r_step_current].reff_xr_rp3,
                                xyz,
                                st->probs);

        bg_prob *= V;

        temp = bg_prob_int.bg_int;
        bg_prob_int.bg_int += bg_prob;
        bg_prob_int.correction += bg_prob - (bg_prob_int.bg_int - temp);

        for (i = 0; i < ap->number_streams; i++)
        {
            st->probs[i].st_prob *= V;
            temp = st->probs[i].st_prob_int;
            st->probs[i].st_prob_int += st->probs[i].st_prob;
            st->probs[i].st_prob_int_c += st->probs[i].st_prob - (st->probs[i].st_prob_int - temp);
        }
    }

    return bg_prob_int;
}

inline static void apply_correction(const unsigned int number_streams,
                                    INTEGRAL_AREA* ia,
                                    INTEGRAL_STATE* st,
                                    BG_PROB bg_prob_int)
{
    unsigned int i;
    ia->background_integral = bg_prob_int.bg_int + bg_prob_int.correction;
    for (i = 0; i < number_streams; i++)
        ia->stream_integrals[i] = st->probs[i].st_prob_int + st->probs[i].st_prob_int_c;
}

inline static BG_PROB nu_sum(const ASTRONOMY_PARAMETERS* ap,
                             const STREAM_CONSTANTS* sc,
                             const STREAM_NUMS* sn,
                             INTEGRAL_AREA* ia,
                             EVALUATION_STATE* es,
                             INTEGRAL_STATE* st,
                             vector* xyz,
                             vector integral_point,
                             const unsigned int mu_step_current)
{
    BG_PROB bg_prob_int = { 0.0, 0.0 };
    BG_PROB r_result;
    unsigned int nu_step_current;

    double mu = ia->mu_min + (mu_step_current * ia->mu_step_size);

    const unsigned int nu_steps = ia->nu_steps;

    for (nu_step_current = 0; nu_step_current < nu_steps; ++nu_step_current)
    {
        apply_correction(ap->number_streams, ia, st, bg_prob_int);

        do_boinc_checkpoint(es, mu_step_current, nu_step_current);

        if (ap->sgr_coordinates == 0)
        {
            gc2lb(ap->wedge,
                  mu + 0.5 * ia->mu_step_size,
                  st->nu_st[nu_step_current].nus,
                  &L(integral_point),
                  &B(integral_point));
        }
        else if (ap->sgr_coordinates == 1)
        {
            gc2sgr(ap->wedge,
                   mu + 0.5 * ia->mu_step_size,
                   st->nu_st[nu_step_current].nus,
                   &L(integral_point),
                   &B(integral_point));
        }
        else
        {
            fprintf(stderr, "Error: ap->sgr_coordinates not valid");
        }

        r_result = r_sum(ap,
                         sc,
                         sn,
                         ia->r_steps,
                         st,
                         xyz,
                         integral_point,
                         nu_step_current);


        bg_prob_int.bg_int += r_result.bg_int;
        bg_prob_int.correction += r_result.correction;
    }

    return bg_prob_int;
}

static void calculate_integral(const ASTRONOMY_PARAMETERS* ap,
                               const STREAM_CONSTANTS* sc,
                               const STREAM_NUMS* sn,
                               vector* xyz,
                               EVALUATION_STATE* es,
                               INTEGRAL_STATE* st)
{
    unsigned int i, mu_step_current;
    vector integral_point;
    BG_PROB bg_prob_int;    /* for Kahan summation */
    BG_PROB nu_result;
    INTEGRAL_AREA* ia = &es->integrals[es->current_integral];

    bg_prob_int.bg_int = ia->background_integral;
    bg_prob_int.correction = 0.0;
    for (i = 0; i < ap->number_streams; i++)
    {
        st->probs[i].st_prob_int = ia->stream_integrals[i];
        st->probs[i].st_prob_int_c = 0.0;
    }

    const unsigned int mu_steps = ia->mu_steps;

    for (mu_step_current = 0; mu_step_current < mu_steps; mu_step_current++)
    {
        nu_result = nu_sum(ap,
                           sc,
                           sn,
                           ia,
                           es,
                           st,
                           xyz,
                           integral_point,
                           mu_step_current);

        bg_prob_int.bg_int += nu_result.bg_int;
        bg_prob_int.correction += nu_result.correction;

    }

    ia->mu_step = 0;

    apply_correction(ap->number_streams, ia, st, bg_prob_int);
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
                               const STREAM_NUMS* sn,
                               const STREAM_GAUSS* sg,
                               EVALUATION_STATE* es,
                               vector* xyz)
{
    unsigned int i, j;
    INTEGRAL_STATE st;

#ifdef MILKYWAY
    read_checkpoint(es);
#endif

    for (; es->current_integral < ap->number_integrals; es->current_integral++)
    {
        prepare_integral_state(ap, sn, sg, &es->integrals[es->current_integral], &st);
        calculate_integral(ap, sc, sn, xyz, es, &st);
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

#ifdef MILKYWAY
    print_stream_integrals(ap, es);
#endif

    return 0;
}

typedef struct
{
    double st_only_sum;
    double st_only_sum_c;
} ST_SUM;

static int calculate_likelihood(const ASTRONOMY_PARAMETERS* ap,
                                const STREAM_CONSTANTS* sc,
                                const STREAM_NUMS* sn,
                                EVALUATION_STATE* es,
                                STREAM_GAUSS* sg,
                                vector* xyz,
                                const STAR_POINTS* sp)
{
    unsigned int i, current_stream;
    double bg_prob;
    double prob_sum, prob_sum_c, temp;  // for Kahan summation
    double exp_background_weight, sum_exp_weights;
    double reff_xr_rp3;

    double bg_only, bg_only_sum, bg_only_sum_c;
    double st_only;

#ifdef MW_ENABLE_DEBUG
    time_t start_time, finish_time;
    time (&start_time);
#endif

    /* The correction terms aren't used here since this isn't the sum? */
    ST_PROBS* st_prob = (ST_PROBS*) malloc(sizeof(ST_PROBS) * ap->number_streams);
    double* exp_stream_weights = malloc(sizeof(double) * ap->number_streams);
    R_STEP_STATE* rss = malloc(sizeof(R_STEP_STATE) * ap->convolve);
    ST_SUM* st_sum = calloc(sizeof(ST_SUM), ap->number_streams);

    exp_background_weight = exp(ap->background_weight);
    sum_exp_weights = exp_background_weight;
    for (i = 0; i < ap->number_streams; i++)
    {
        exp_stream_weights[i] = exp(ap->stream[i].weights);
        sum_exp_weights += exp(ap->stream[i].weights);
    }
    sum_exp_weights *= 0.001;

    do_boinc_checkpoint(es, 0, 0); /* CHECKME: Steps? */

    prob_sum = 0.0;
    prob_sum_c = 0.0;

    bg_only_sum = 0.0;
    bg_only_sum_c = 0.0;

    for (; es->current_star_point < sp->number_stars; es->current_star_point++)
    {
        double star_prob;

        reff_xr_rp3 = set_probability_constants(sn,
                                                sg,
                                                ap->convolve,
                                                ZN(sp, es->current_star_point),
                                                rss,
                                                0);

        bg_prob = calculate_bg_probability(ap,
                                           sn,
                                           rss,
                                           0, /* Would be indexing the 2D block used by integration */
                                           reff_xr_rp3,
                                           &VN(sp, es->current_star_point),
                                           xyz);

        calculate_probabilities(ap,
                                sc,
                                rss,
                                0, /* Would be indexing the 2D block used by integration */
                                reff_xr_rp3,
                                xyz,
                                st_prob);

        bg_only = (bg_prob / es->background_integral) * exp_background_weight;
        star_prob = bg_only;

        if (bg_only == 0.0)
            bg_only = -238.0;
        else
            bg_only = log10(bg_only / sum_exp_weights);

        temp = bg_only_sum;
        bg_only_sum += bg_only;
        bg_only_sum_c += bg_only - (bg_only_sum - temp);


        for (current_stream = 0; current_stream < ap->number_streams; current_stream++)
        {
            st_only = st_prob[current_stream].st_prob / es->stream_integrals[current_stream] * exp_stream_weights[current_stream];
            star_prob += st_only;

            if (st_only == 0.0)
                st_only = -238.0;
            else
                st_only = log10(st_only / sum_exp_weights);

            temp = st_sum[current_stream].st_only_sum;
            st_sum[current_stream].st_only_sum += st_only;
            st_sum[current_stream].st_only_sum_c += st_only - (st_sum[current_stream].st_only_sum - temp);
        }
        star_prob /= sum_exp_weights;

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
    }
    es->prob_sum = prob_sum + prob_sum_c;
    bg_only_sum += bg_only_sum_c;
    bg_only_sum /= sp->number_stars;

#ifdef MILKYWAY
    fprintf(stderr, "<background_only_likelihood> %.20lf </background_only_likelihood>\n", bg_only_sum - 3.0);
    fprintf(stderr, "<stream_only_likelihood>");
    for (i = 0; i < ap->number_streams; i++)
    {
        st_sum[i].st_only_sum += st_sum[i].st_only_sum_c;
        st_sum[i].st_only_sum /= sp->number_stars;

        fprintf(stderr, " %.20lf", st_sum[i].st_only_sum - 3.0);
    }
    fprintf(stderr, " </stream_only_likelihood>\n");
#endif

    free(exp_stream_weights);
    free(st_prob);
    free(rss);

    free(st_sum);

#ifdef MW_ENABLE_DEBUG
    time(&finish_time);
    MW_DEBUG("likelihood calculated in: %lf\n", (double)finish_time - (double)start_time);
#endif

    return 0;
}


double cpu_evaluate(const ASTRONOMY_PARAMETERS* ap,
                    const STAR_POINTS* sp,
                    const STREAM_CONSTANTS* sc,
                    const STREAM_NUMS* sn)
{
    int retval;
    EVALUATION_STATE es = EMPTY_EVALUATION_STATE;
    STREAM_GAUSS sg;

    initialize_state(ap, sp, &es);
    get_stream_gauss(ap, &sg);

    reset_evaluation_state(&es);

    vector* xyz = malloc(sizeof(vector) * ap->convolve);

    retval = calculate_integrals(ap, sc, sn, &sg, &es, xyz);
    if (retval)
    {
        fprintf(stderr, "APP: error calculating integrals: %d\n", retval);
        mw_finish(retval);
    }

    retval = calculate_likelihood(ap, sc, sn, &es, &sg, xyz, sp);
    if (retval)
    {
        fprintf(stderr, "APP: error calculating likelihood: %d\n", retval);
        mw_finish(retval);
    }

    free(xyz);
    free(sg.dx);
    free(sg.qgaus_W);

    /*  log10(x * 0.001) = log10(x) - 3.0 */
    return (es.prob_sum / (sp->number_stars - es.bad_jacobians)) - 3.0;
}


