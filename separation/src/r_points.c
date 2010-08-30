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

__attribute__ ((always_inline, const))
inline _MW_STATIC real distance_magnitude(const real m)
{
    return mw_powr(R10, (m - (real) 14.2) / R5);
}

R_PRIME calcRPrime(__MW_CONSTANT INTEGRAL_AREA* ia, const unsigned int r_step)
{
    real r, next_r, log_r;
    R_PRIME ret;

    log_r = ia->r_min + (r_step * ia->r_step_size);
    r = distance_magnitude(log_r);
    next_r = distance_magnitude(log_r + ia->r_step_size);

    ret.irv = d2r(((cube(next_r) - cube(r)) / R3) * ia->mu_step_size);
    ret.rPrime = (next_r + r) / R2;

    return ret;
}

void set_r_points(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                  __MW_CONSTANT STREAM_GAUSS* sg,
                  const unsigned int n_convolve,
                  const real coords,
                  __MW_LOCAL R_POINTS* r_pts)
{
    real g, exponent, r3, N;
    unsigned int i;

    const real gPrime = calcGPrime(coords);

    for (i = 0; i < n_convolve; ++i)
    {
        g = gPrime + sg[i].dx;

        /* MAG2R */
        r_pts[i].r_in_mag = g;
        r_pts[i].r_in_mag2 = sqr(g);
        r_pts[i].r_point = mw_powr(R10, (g - absm) / R5 + R1) / R1000;

        r3 = cube(r_pts[i].r_point);
        exponent = sqr(g - gPrime) / (R2 * sqr(stdev));
        N = ap->coeff * mw_exp(-exponent);
        r_pts[i].qw_r3_N = sg[i].qgaus_W * r3 * N;
    }
}

