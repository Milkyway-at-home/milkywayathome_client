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
#include "separation_constants.h"
#include "milkyway_cl.h"
#include "milkyway_extra.h"
#include "r_points.h"

/* Literals are assumed to be doubles by default, and the
 * -cl-single-precision-constant flag seems to not be working when
 * trying to use float */
#define R1 ((real) 1.0)
#define R2 ((real) 2.0)
#define R5 ((real) 5.0)
#define R10 ((real) 10.0)
#define R1000 ((real) 1000.0)

real set_r_points(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                  __MW_CONSTANT STREAM_GAUSS* sg,
                  const unsigned int n_convolve,
                  const real coords,
                  __MW_LOCAL R_POINTS* r_pts)
{
    real g, exponent, r3, N;
    unsigned int i;

    /* R2MAG */
    _MW_STATIC const real sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };
    const real gPrime = R5 * (mw_log10(coords * R1000) - R1) + absm;

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

    /* REFF */
    const real exp_result = mw_exp(sigmoid_curve_params[1] * (gPrime - sigmoid_curve_params[2]));
    const real reff_value = sigmoid_curve_params[0] / (exp_result + R1);
    const real rPrime3 = cube(coords);
    const real reff_xr_rp3 = reff_value * xr / rPrime3;
    return reff_xr_rp3;
}

