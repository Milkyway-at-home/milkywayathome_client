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

double sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };

double* dx;
vector* xyz;
double* qgaus_W;

STREAM_CONSTANTS* init_constants(ASTRONOMY_PARAMETERS* ap, STREAM_NUMS* sn)
{
    unsigned int i;
    double* qgaus_X;
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
        double ra, dec, lamda, beta, l, b;
        vector lbr;

        sc[i].stream_sigma = STREAM_PARAM_N(ap, i).stream_parameters[4];
        sc[i].stream_sigma_sq2 = 2.0 * sc[i].stream_sigma * sc[i].stream_sigma;

        if (ap->sgr_coordinates == 0)
        {
            atGCToEq(STREAM_PARAM_N(ap, i).stream_parameters[0],
                     0, &ra, &dec, get_node(), wedge_incl(ap->wedge));
            atEqToGal(ra, dec, &l, &b);
        }
        else if (ap->sgr_coordinates == 1)
        {
            gcToSgr(STREAM_PARAM_N(ap, i).stream_parameters[0], 0, ap->wedge, &lamda, &beta);
            sgrToGal(lamda, beta, &l, &b);
        }
        else
        {
            fprintf(stderr, "Error: sgr_coordinates not valid");
        }

        L(lbr) = l;
        B(lbr) = b;
        R(lbr) = STREAM_PARAM_N(ap, i).stream_parameters[1];
        lbr2xyz(lbr, sc[i].stream_c);

        X(sc[i].stream_a) =   sin(STREAM_PARAM_N(ap, i).stream_parameters[2])
                            * cos(STREAM_PARAM_N(ap, i).stream_parameters[3]);

        Y(sc[i].stream_a) =   sin(STREAM_PARAM_N(ap, i).stream_parameters[2])
                            * sin(STREAM_PARAM_N(ap, i).stream_parameters[3]);

        Z(sc[i].stream_a) = cos(STREAM_PARAM_N(ap, i).stream_parameters[2]);
    }

    qgaus_X = (double*) malloc(sizeof(double) * ap->convolve);

    xyz     = (vector*) malloc(sizeof(vector) * ap->convolve);
    qgaus_W = (double*) malloc(sizeof(double) * ap->convolve);
    dx      = (double*) malloc(sizeof(double) * ap->convolve);

    gaussLegendre(-1.0, 1.0, qgaus_X, qgaus_W, ap->convolve);

    for (i = 0; i < ap->convolve; i++)
        dx[i] = 3.0 * stdev * qgaus_X[i];

    free(qgaus_X);

    return sc;
}

void free_constants(ASTRONOMY_PARAMETERS* ap, STREAM_CONSTANTS* sc)
{
    //free(qgaus_X);
    //free(qgaus_W);
    free(dx);
    free(xyz);
    free(sc);
}

void set_probability_constants(const STREAM_NUMS* sn,
                               unsigned int n_convolve,
                               double coords,
                               R_STEP_STATE* rss,
                               unsigned int r_step,
                               double* reff_xr_rp3)
{
    double gPrime, exp_result, g, exponent, r3, N, reff_value, rPrime3;
    unsigned int i, idx;
    const unsigned int row_base = n_convolve * r_step;

    //R2MAG
    gPrime = 5.0 * (log10(coords * 1000.0) - 1.0) + absm;

    //REFF
    exp_result = exp(sigmoid_curve_params[1] * (gPrime - sigmoid_curve_params[2]));
    reff_value = sigmoid_curve_params[0] / (exp_result + 1);
    rPrime3 = coords * coords * coords;

    *reff_xr_rp3 = reff_value * xr / rPrime3;

    for (i = 0; i < n_convolve; i++)
    {
        idx = row_base + i;
        g = gPrime + dx[i];

        //MAG2R
        rss[idx].r_in_mag = g;
        rss[idx].r_in_mag2 = g * g;
        rss[idx].r_point = pow(10.0, (g - absm) / 5.0 + 1.0) / 1000.0;

        r3 = rss[idx].r_point * rss[idx].r_point * rss[idx].r_point;
        exponent = (g - gPrime) * (g - gPrime) / (2.0 * stdev * stdev);
        N = sn->coeff * exp(-exponent);
        rss[idx].qw_r3_N = qgaus_W[i] * r3 * N;
    }
}

void calculate_probabilities(const ASTRONOMY_PARAMETERS* ap,
                             const STREAM_CONSTANTS* cs,
                             const STREAM_NUMS* sn,
                             R_STEP_STATE* rss,
                             unsigned int r_step_current,
                             unsigned int r_steps,
                             double reff_xr_rp3,
                             double* integral_point,
                             double* bg_prob,
                             ST_PROBS* probs)
{
    double bsin, lsin, bcos, lcos, zp;
    double rg, rs, xyzs[3], dotted, xyz_norm;
    double h_prob, aux_prob;
    unsigned int i, j, idx;
    const unsigned int row_base = r_step_current * ap->convolve;

    bsin = sin(integral_point[1] / deg);
    lsin = sin(integral_point[0] / deg);
    bcos = cos(integral_point[1] / deg);
    lcos = cos(integral_point[0] / deg);

    /* if q is 0, there is no probability */
    if (sn->q == 0)
    {
        *bg_prob = -1;
    }
    else
    {
        *bg_prob = 0;
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

//vickej2_bg changing the hernquist profile to include a quadratic term in g

                if (ap->aux_bg_profile == 1)
                {
                    h_prob = rss[idx].qw_r3_N / (rg * rs * rs * rs);
                    aux_prob = rss[idx].qw_r3_N * (sn->bg_a * rss[idx].r_in_mag2 + sn->bg_b * rss[idx].r_in_mag + sn->bg_c );
                    *bg_prob += h_prob + aux_prob;
                }
                else if (ap->aux_bg_profile == 0)
                {
                    *bg_prob += rss[idx].qw_r3_N / (rg * rs * rs * rs);
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

                *bg_prob += rss[idx].qw_r3_N / (pow(rg, sn->alpha) * pow(rg + sn->r0, sn->alpha_delta3));
            }
        }
        *bg_prob *= reff_xr_rp3;
    }

    for (i = 0; i < ap->number_streams; i++)
    {
        probs[i].st_prob = 0.0;
        if (cs[i].stream_sigma > -0.0001 && cs[i].stream_sigma < 0.0001)
            continue;
        for (j = 0; j < ap->convolve; j++)
        {
            X(xyzs) = X(xyz[j]) - X(cs[i].stream_c);
            Y(xyzs) = Y(xyz[j]) - Y(cs[i].stream_c);
            Z(xyzs) = Z(xyz[j]) - Z(cs[i].stream_c);

            dotted = X(cs[i].stream_a) * X(xyzs)
                   + Y(cs[i].stream_a) * Y(xyzs)
                   + Z(cs[i].stream_a) * Z(xyzs);

            X(xyzs) = X(xyzs) - dotted * X(cs[i].stream_a);
            Y(xyzs) = Y(xyzs) - dotted * Y(cs[i].stream_a);
            Z(xyzs) = Z(xyzs) - dotted * Z(cs[i].stream_a);

            xyz_norm = X(xyzs) * X(xyzs)
                     + Y(xyzs) * Y(xyzs)
                     + Z(xyzs) * Z(xyzs);

            idx = row_base + j;
            probs[i].st_prob += rss[idx].qw_r3_N * exp(-xyz_norm / cs[i].stream_sigma_sq2);
        }
        probs[i].st_prob *= reff_xr_rp3;
    }
}

double calculate_progress(EVALUATION_STATE* es)
{
    double total_calc_probs, current_calc_probs, current_probs;
    unsigned int i, mu_step_current, nu_step_current, r_step_current;
    INTEGRAL_AREA* ia;

    total_calc_probs = 0;
    current_calc_probs = 0;

    for (i = 0; i < es->number_integrals; i++)
    {
        ia = &es->integrals[i];

        get_steps(ia, &mu_step_current, &nu_step_current, &r_step_current);

        current_probs = ia->r_steps * ia->mu_steps * ia->nu_steps;
        total_calc_probs += current_probs;
        if (i < es->current_integral)
        {
            current_calc_probs += current_probs;
        }
        else if (i == es->current_integral)
        {
            current_calc_probs += r_step_current + (nu_step_current * ia->r_steps) + (mu_step_current * ia->nu_steps * ia->r_steps);
        }

    }

    total_calc_probs += es->total_stars;
    current_calc_probs += es->current_star_point;

    return (double)current_calc_probs / (double)total_calc_probs;
}

#ifdef MILKYWAY
void do_boinc_checkpoint(EVALUATION_STATE* es)
{
    double progress;

    if (boinc_time_to_checkpoint())
    {
        int retval = write_checkpoint(es);
        if (retval)
        {
            fprintf(stderr, "APP: astronomy checkpoint failed %d\n", retval);
            return;
        }
        boinc_checkpoint_completed();
    }

    progress = calculate_progress(es);
//      printf("progress: %.10f\n", progress);
    boinc_fraction_done(progress);
}
#endif

void cpu__r_constants(const STREAM_NUMS* sn,
                      unsigned int n_convolve,
                      unsigned int r_steps,
                      double r_min,
                      double r_step_size,
                      double mu_step_size,
                      unsigned int nu_steps,
                      double nu_min,
                      double nu_step_size,
                      double* irv,
                      R_STEP_STATE* rss,
                      double* reff_xr_rp3,
                      double* nus,
                      double* ids)
{
    unsigned int i;

//vickej2_kpc edits to make volumes even in kpc rather than g
//vickej2_kpc        double log_r, r, next_r, rPrime;

    double r, next_r, rPrime, r_min_kpc, r_max_kpc, r_step_size_kpc, r_max;

    r_max = r_min + r_step_size * r_steps;

    r_min_kpc = pow(10.0, ((r_min - 14.2) / 5.0));
    r_max_kpc = pow(10.0, ((r_max - 14.2) / 5.0));
    r_step_size_kpc = (r_max_kpc - r_min_kpc) / r_steps;


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

        irv[i] = (((next_r * next_r * next_r) - (r * r * r)) / 3.0) * mu_step_size / deg;
        rPrime = (next_r + r) / 2.0;

        set_probability_constants(sn,
                                  n_convolve,
                                  rPrime,
                                  rss,
                                  i,
                                  &reff_xr_rp3[i]);
    }

    for (i = 0; i < nu_steps; i++)
    {
        nus[i] = nu_min + (i * nu_step_size);
        ids[i] = cos((90.0 - nus[i] - nu_step_size) / deg) - cos((90.0 - nus[i]) / deg);
        nus[i] += 0.5 * nu_step_size;
    }
}

void prepare_integral_state(const ASTRONOMY_PARAMETERS* ap,
                            const STREAM_NUMS* sn,
                            INTEGRAL_AREA* ia,
                            INTEGRAL_STATE* st)
{

    st->probs = (ST_PROBS*) malloc(sizeof(ST_PROBS) * ap->number_streams);

    st->irv         = (double*)malloc(sizeof(double) * ia->r_steps);
    st->reff_xr_rp3 = (double*)malloc(sizeof(double) * ia->r_steps);


    /* 2D block, ia->r_steps = rows, ap->convolve = columns */
    st->rss = (R_STEP_STATE*) malloc(sizeof(R_STEP_STATE) * ia->r_steps * ap->convolve);

    st->ids     = (double*)malloc(sizeof(double) * ia->nu_steps);
    st->nus     = (double*)malloc(sizeof(double) * ia->nu_steps);
    cpu__r_constants(sn,
                     ap->convolve, ia->r_steps, ia->r_min, ia->r_step_size,
                     ia->mu_step_size,
                     ia->nu_steps, ia->nu_min, ia->nu_step_size,
                     st->irv,
                     st->rss, st->reff_xr_rp3, st->nus, st->ids);

}

void free_integral_state(INTEGRAL_AREA* ia, INTEGRAL_STATE* st)
{
    free(st->nus);
    free(st->ids);
    free(st->irv);
    free(st->probs);
    free(st->reff_xr_rp3);
    free(st->rss);
}

void calculate_integral(const ASTRONOMY_PARAMETERS* ap,
                        const STREAM_CONSTANTS* sc,
                        const STREAM_NUMS* sn,
                        INTEGRAL_AREA* ia,
                        EVALUATION_STATE* es,
                        INTEGRAL_STATE* st)
{
    unsigned int i, mu_step_current, nu_step_current, r_step_current;
    double integral_point[3];
    double temp;
    double V;
    double bg_prob;
    double bg_prob_int, bg_prob_int_c; /* for kahan summation */


    get_steps(ia, &mu_step_current, &nu_step_current, &r_step_current);

    bg_prob_int = ia->background_integral;
    bg_prob_int_c = 0.0;
    for (i = 0; i < ap->number_streams; i++)
    {
        st->probs[i].st_prob_int = ia->stream_integrals[i];
        st->probs[i].st_prob_int_c = 0.0;
    }

    for (; mu_step_current < ia->mu_steps; mu_step_current++)
    {
        double mu = ia->mu_min + (mu_step_current * ia->mu_step_size);

        for (; nu_step_current < ia->nu_steps; nu_step_current++)
        {
#ifdef MILKYWAY
            ia->background_integral = bg_prob_int + bg_prob_int_c;  // apply correction
            for (i = 0; i < ap->number_streams; i++)
            {
                ia->stream_integrals[i] = st->probs[i].st_prob_int + st->probs[i].st_prob_int_c;  // apply correction
            }
            ia->mu_step = mu_step_current;
            ia->nu_step = nu_step_current;
            ia->r_step = r_step_current;

            do_boinc_checkpoint(es);

//              bg_prob_int_c = 0;
//              for (i = 0; i < ap->number_streams; i++) st_probs_int_c[i] = 0;
#endif

            if (ap->sgr_coordinates == 0)
            {
                double ra, dec;
                atGCToEq(mu + 0.5 * ia->mu_step_size, st->nus[nu_step_current], &ra, &dec, get_node(), wedge_incl(ap->wedge));
                atEqToGal(ra, dec, &integral_point[0], &integral_point[1]);
            }
            else if (ap->sgr_coordinates == 1)
            {
                double lamda, beta;
                gcToSgr(mu + 0.5 * ia->mu_step_size, st->nus[nu_step_current], ap->wedge, &lamda, &beta);
                sgrToGal(lamda, beta, &integral_point[0], &integral_point[1]);
            }
            else
            {
                printf("Error: ap->sgr_coordinates not valid");
            }

            for (; r_step_current < ia->r_steps; r_step_current++)
            {
                V = st->irv[r_step_current] * st->ids[nu_step_current];

                calculate_probabilities(ap,
                                        sc,
                                        sn,
                                        st->rss,
                                        r_step_current,
                                        ia->r_steps,
                                        st->reff_xr_rp3[r_step_current],
                                        integral_point,
                                        &bg_prob,
                                        st->probs);

                bg_prob *= V;

                temp = bg_prob_int;
                bg_prob_int += bg_prob;
                bg_prob_int_c += bg_prob - (bg_prob_int - temp);

//              ia->background_integral += bg_prob;
                for (i = 0; i < ap->number_streams; i++)
                {
                    st->probs[i].st_prob *= V;
                    temp = st->probs[i].st_prob_int;
                    st->probs[i].st_prob_int += st->probs[i].st_prob;
                    st->probs[i].st_prob_int_c += st->probs[i].st_prob - (st->probs[i].st_prob_int - temp);

//                  ia->stream_integrals[i] += st_probs[i] * V;
                }

#ifndef MILKYWAY
                ia->current_calculation++;
                if (ia->current_calculation >= ia->max_calculation)
                    break;
#endif
            }
#ifndef MILKYWAY
            if (ia->current_calculation >= ia->max_calculation)
                break;
#endif
            r_step_current = 0;
        }
#ifndef MILKYWAY
        if (ia->current_calculation >= ia->max_calculation)
            break;
#endif
        nu_step_current = 0;
    }
    mu_step_current = 0;

    ia->background_integral = bg_prob_int + bg_prob_int_c;  // apply correction
    for (i = 0; i < ap->number_streams; i++)
    {
        ia->stream_integrals[i] = st->probs[i].st_prob_int + st->probs[i].st_prob_int_c;  // apply correction
    }

//  printf("bg_int: %.15lf ", ia->background_integral);
//  for (i = 0; i < ap->number_streams; i++) printf("st_int[%d]: %.15lf ", i, ia->stream_integrals[i]);
//  printf("\n");

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

int calculate_integrals(const ASTRONOMY_PARAMETERS* ap,
                        const STREAM_CONSTANTS* sc,
                        const STREAM_NUMS* sn,
                        EVALUATION_STATE* es)
{
    unsigned int i, j;

#ifdef MILKYWAY
    read_checkpoint(es);
#endif

    INTEGRAL_STATE st;
    /* FIXME: the integral area not actually needed here, for some
     * reason they all carry the same information which never
     * changes. */
    prepare_integral_state(ap, sn, &es->integrals[0], &st);

    for (; es->current_integral < ap->number_integrals; es->current_integral++)
        calculate_integral(ap, sc, sn, &es->integrals[es->current_integral], es, &st);

    free_integral_state(&es->integrals[0], &st);

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

int calculate_likelihood(const ASTRONOMY_PARAMETERS* ap,
                         const STREAM_CONSTANTS* sc,
                         const STREAM_NUMS* sn,
                         EVALUATION_STATE* es,
                         const STAR_POINTS* sp)
{
    unsigned int i, current_stream;
    double bg_prob;
    double prob_sum, prob_sum_c, temp;  // for Kahan summation
    double exp_background_weight, sum_exp_weights, *exp_stream_weights;
    double reff_xr_rp3;

    double bg_only, bg_only_sum, bg_only_sum_c;
    double st_only, *st_only_sum, *st_only_sum_c;

#ifdef MW_ENABLE_DEBUG
    time_t start_time, finish_time;
    time (&start_time);
#endif

    /* The correction terms aren't used here since this isn't the sum? */
    ST_PROBS* st_prob = (ST_PROBS*) malloc(sizeof(ST_PROBS) * ap->number_streams);

    exp_stream_weights = (double*)malloc(sizeof(double) * ap->number_streams);

    R_STEP_STATE* rss = malloc(sizeof(R_STEP_STATE) * ap->convolve);

    st_only_sum = (double*)malloc(sizeof(double) * ap->number_streams);
    st_only_sum_c = (double*)malloc(sizeof(double) * ap->number_streams);

    exp_background_weight = exp(ap->background_weight);
    sum_exp_weights = exp_background_weight;
    for (i = 0; i < ap->number_streams; i++)
    {
        exp_stream_weights[i] = exp(STREAM_N(ap, i).weights);
        sum_exp_weights += exp(STREAM_N(ap, i).weights);
    }
    sum_exp_weights *= 0.001;

#ifdef MILKYWAY
    do_boinc_checkpoint(es);
#endif

    prob_sum = 0;
    prob_sum_c = 0;

    bg_only_sum = 0;
    bg_only_sum_c = 0;

    for (i = 0; i < ap->number_streams; i++)
    {
        st_only_sum[i] = 0;
        st_only_sum_c[i] = 0;
    }


    for (; es->current_star_point < sp->number_stars; es->current_star_point++)
    {
        double star_prob;

        set_probability_constants(sn,
                                  ap->convolve,
                                  ZN(sp, es->current_star_point),
                                  rss,
                                  0,
                                  &reff_xr_rp3);

        calculate_probabilities(ap,
                                sc,
                                sn,
                                rss,
                                0,   /* Would be for indexing the 2D block used by the integration */
                                0,
                                reff_xr_rp3,
                                &VN(sp, es->current_star_point),
                                &bg_prob,
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

            temp = st_only_sum[current_stream];
            st_only_sum[current_stream] += st_only;
            st_only_sum_c[current_stream] += st_only - (st_only_sum[current_stream] - temp);
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
        st_only_sum[i] += st_only_sum_c[i];
        st_only_sum[i] /= sp->number_stars;

        fprintf(stderr, " %.20lf", st_only_sum[i] - 3.0);
    }
    fprintf(stderr, " </stream_only_likelihood>\n");
#endif

    free(exp_stream_weights);
    free(st_prob);
    free(rss);

    free(st_only_sum);
    free(st_only_sum_c);

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

    initialize_state(ap, sp, &es);

    reset_evaluation_state(&es);

    retval = calculate_integrals(ap, sc, sn, &es);
    if (retval)
    {
        fprintf(stderr, "APP: error calculating integrals: %d\n", retval);
        mw_finish(retval);
    }

    retval = calculate_likelihood(ap, sc, sn, &es, sp);
    if (retval)
    {
        fprintf(stderr, "APP: error calculating likelihood: %d\n", retval);
        mw_finish(retval);
    }

    /*  log10(x * 0.001) = log10(x) - 3.0 */
    return (es.prob_sum / (sp->number_stars - es.bad_jacobians)) - 3.0;
}


