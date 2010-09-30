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

#include "separation_types.h"
#include "milkyway_cl.h"
#include "milkyway_extra.h"
#include "r_points.h"

#ifndef NDEBUG
extern real calcReffXrRp3(const real coords);
extern real calcGPrime(const real coords);
extern R_PRIME calcRPrime(__MW_CONSTANT INTEGRAL_AREA* ia, const unsigned int r_step);
#endif


static inline R_POINTS calc_r_point(__MW_CONSTANT STREAM_GAUSS* sg, const real gPrime, const real coeff)
{
    R_POINTS r_pt;
    real g, exponent, r3, N;

    g = gPrime + sg->dx;

    /* MAG2R */
    r_pt.r_in_mag = g;
    r_pt.r_in_mag2 = sqr(g);
    r_pt.r_point = mw_powr(RL10, (g - absm) / RL5 + RL1) / RL1000;

    r3 = cube(r_pt.r_point);
    exponent = sqr(g - gPrime) / (RL2 * sqr(stdev));
    N = coeff * mw_exp(-exponent);
    r_pt.qw_r3_N = sg->qgaus_W * r3 * N;

    return r_pt;
}


void set_r_points(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                  __MW_CONSTANT STREAM_GAUSS* sg,
                  const unsigned int n_convolve,
                  const real coords,
                  __MW_LOCAL R_POINTS* r_pts)
{
    unsigned int i;

    const real gPrime = calcGPrime(coords);

    for (i = 0; i < n_convolve; ++i)
        r_pts[i] = calc_r_point(&sg[i], gPrime, ap->coeff);
}

