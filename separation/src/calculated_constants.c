/*
Copyright 2008-2010 Travis Desell, Dave Przybylo, Nathan Cole, Matthew Arsenault,
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

#include "separation_types.h"
#include "separation_constants.h"
#include "calculated_constants.h"
#include "milkyway_util.h"
#include "coordinates.h"
#include "gauss_legendre.h"

/* Convert sun-centered lbr (degrees) into galactic xyz coordinates. */
static void lbr2xyz(const vector lbr, vector xyz)
{
    real zp, d;
/* TODO: Use radians to begin with */

    real lsin, lcos;
    real bsin, bcos;

    mw_sincos(d2r(B(lbr)), &bsin, &bcos);
    mw_sincos(d2r(L(lbr)), &lsin, &lcos);

    Z(xyz) = R(lbr) * bsin;
    zp = R(lbr) * bcos;
    d = mw_sqrt(sqr(sun_r0) + sqr(zp) - 2.0 * sun_r0 * zp * lcos);
    X(xyz) = (sqr(zp) - sqr(sun_r0) - sqr(d)) / (2.0 * sun_r0);
    Y(xyz) = zp * lsin;
}

static inline void stream_a(vector a, real* parameters)
{
    X(a) = mw_sin(parameters[2]) * mw_cos(parameters[3]);
    Y(a) = mw_sin(parameters[2]) * mw_sin(parameters[3]);
    Z(a) = mw_cos(parameters[2]);
}

static inline void stream_c(vector c, int wedge, real mu, real r)
{
    LB lb;
    vector lbr;

    lb = gc2lb(wedge, mu, 0.0);

    L(lbr) = LB_L(lb);
    B(lbr) = LB_B(lb);
    R(lbr) = r;
    lbr2xyz(lbr, c);
}
int setAstronomyParameters(ASTRONOMY_PARAMETERS* ap, const BACKGROUND_PARAMETERS* bgp)

{
    STREAM_CONSTANTS* sc;

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
        warn("Error: aux_bg_profile invalid");
        return 1;
    }

    if (ap->sgr_coordinates)
    {
        warn("gc2sgr not implemented\n");
        return 1;
    }

    ap->coeff = 1.0 / (stdev * SQRT_2PI);
    ap->alpha_delta3 = 3.0 - ap->alpha + ap->delta;

    ap->fast_h_prob = (ap->alpha == 1 && ap->delta == 1);

    return 0;
}

STREAM_CONSTANTS* getStreamConstants(const ASTRONOMY_PARAMETERS* ap, const STREAMS* streams)

{
    unsigned int i;
    STREAM_CONSTANTS* sc;
    real stream_sigma;

    sc = (STREAM_CONSTANTS*) mwMallocAligned(sizeof(STREAM_CONSTANTS) * streams->number_streams,
                                             sizeof(STREAM_CONSTANTS));

    for (i = 0; i < streams->number_streams; i++)
    {
        stream_sigma = streams->parameters[i].stream_parameters[4];
        sc[i].large_sigma = (stream_sigma > SIGMA_LIMIT || stream_sigma < -SIGMA_LIMIT);
        sc[i].sigma_sq2 = 2.0 * sqr(stream_sigma);

        stream_a(sc[i].a, streams->parameters[i].stream_parameters);
        stream_c(sc[i].c,
                 ap->wedge,
                 streams->parameters[i].stream_parameters[0],
                 streams->parameters[i].stream_parameters[1]);
    }

    return sc;
}

STREAM_GAUSS* get_stream_gauss(const unsigned int convolve)
{
    unsigned int i;
    STREAM_GAUSS* sg;

    real* qgaus_X = (real*) mallocSafe(sizeof(real) * convolve);
    real* qgaus_W = (real*) mallocSafe(sizeof(real) * convolve);

    gaussLegendre(-1.0, 1.0, qgaus_X, qgaus_W, convolve);

    sg = (STREAM_GAUSS*) mallocSafe(sizeof(STREAM_GAUSS) * convolve);

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

NU_CONSTANTS* prepare_nu_constants(const unsigned int nu_steps,
                                   const real nu_step_size,
                                   const real nu_min)
{
    unsigned int i;
    real tmp1, tmp2;
    NU_CONSTANTS* nu_consts;

    nu_consts = (NU_CONSTANTS*) mallocSafe(sizeof(NU_CONSTANTS) * nu_steps);

    for (i = 0; i < nu_steps; ++i)
    {
        nu_consts[i].nu = nu_min + (i * nu_step_size);

        tmp1 = d2r(90.0 - nu_consts[i].nu - nu_step_size);
        tmp2 = d2r(90.0 - nu_consts[i].nu);

        nu_consts[i].id = mw_cos(tmp1) - mw_cos(tmp2);
        nu_consts[i].nu += 0.5 * nu_step_size;
    }

    return nu_consts;
}

