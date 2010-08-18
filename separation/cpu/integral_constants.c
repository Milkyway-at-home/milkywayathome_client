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

static const double sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };

STREAM_CONSTANTS* init_constants(ASTRONOMY_PARAMETERS* ap,
                                 const BACKGROUND_PARAMETERS* bgp,
                                 const STREAMS* streams)
{
    unsigned int i;
    vector lbr;

    STREAM_CONSTANTS* sc = mallocSafe(sizeof(STREAM_CONSTANTS) * streams->number_streams);

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

        if (ap->sgr_coordinates)
        {
            fprintf(stderr, "gc2sgr probably broken right now, so refusing to run\n");
            mw_finish(EXIT_FAILURE);
        }

        gc2lb(ap->wedge, streams->parameters[i].stream_parameters[0], 0, &L(lbr), &B(lbr));

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

STREAM_GAUSS* get_stream_gauss(const unsigned int convolve)
{
    unsigned int i;
    STREAM_GAUSS* sg;

    double* qgaus_X = mallocSafe(sizeof(double) * convolve);
    double* qgaus_W = mallocSafe(sizeof(double) * convolve);

    gaussLegendre(-1.0, 1.0, qgaus_X, qgaus_W, convolve);

    sg = mallocSafe(sizeof(STREAM_GAUSS) * convolve);

    /* Use separate buffers at first since that's what the gaussLegendre takes,
       but then pack them into a more coherent struct */
    for (i = 0; i < convolve; ++i)
    {
        sg[i].dx = 3.0 * stdev * qgaus_X[i];
        sg[i].qgaus_W = qgaus_W[i];
    }

    free(qgaus_X);
    free(qgaus_W);

    return sg;

}

double set_r_points(const ASTRONOMY_PARAMETERS* ap,
                    const STREAM_GAUSS* sg,
                    const unsigned int n_convolve,
                    const double coords,
                    R_POINTS* r_pts)
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

    for (i = 0; i < n_convolve; ++i)
    {
        g = gPrime + sg[i].dx;

        //MAG2R
        r_pts[i].r_in_mag = g;
        r_pts[i].r_in_mag2 = sqr(g);
        r_pts[i].r_point = pow(10.0, (g - absm) / 5.0 + 1.0) / 1000.0;

        r3 = cube(r_pts[i].r_point);
        exponent = sqr(g - gPrime) / (2.0 * sqr(stdev));
        N = ap->coeff * exp(-exponent);
        r_pts[i].qw_r3_N = sg[i].qgaus_W * r3 * N;
    }

    reff_xr_rp3 = reff_value * xr / rPrime3;
    return reff_xr_rp3;
}

NU_CONSTANTS* prepare_nu_constants(const unsigned int nu_steps,
                                   const double nu_step_size,
                                   const double nu_min)
{
    unsigned int i;
    double tmp1, tmp2;

    NU_CONSTANTS* nu_consts = mallocSafe(sizeof(NU_CONSTANTS) * nu_steps);

    for (i = 0; i < nu_steps; i++)
    {
        nu_consts[i].nu = nu_min + (i * nu_step_size);

        tmp1 = d2r(90.0 - nu_consts[i].nu - nu_step_size);
        tmp2 = d2r(90.0 - nu_consts[i].nu);

        nu_consts[i].id = cos(tmp1) - cos(tmp2);
        nu_consts[i].nu += 0.5 * nu_step_size;
    }

    return nu_consts;
}

