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

    ap->coeff = 1.0 / (stdev * sqrt(M_2PI));
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

void get_stream_gauss(const unsigned int convolve, STREAM_GAUSS* sg)
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

double set_prob_consts(const ASTRONOMY_PARAMETERS* ap,
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

void prepare_nu_constants(NU_CONSTANTS* nu_st,
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

R_CONSTANTS* prepare_r_constants(const ASTRONOMY_PARAMETERS* ap,
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

void prepare_integral_constants(const ASTRONOMY_PARAMETERS* ap,
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

void free_integral_constants(INTEGRAL_CONSTANTS* ic)
{
    free(ic->r_step_consts);
    free(ic->rss);
    free(ic->nu_st);
}

void free_stream_gauss(STREAM_GAUSS* sg)
{
    free(sg->dx);
    free(sg->qgaus_W);
}

