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
#include "separation_cl.h"
#include "milkyway_math.h"
#include "r_points.h"

static const real sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };

real set_r_points(__MW_PRIVATE const ASTRONOMY_PARAMETERS* ap,
                  __MW_CONSTANT const STREAM_GAUSS* sg,
                  const unsigned int n_convolve,
                  const real coords,
                  __MW_LOCAL R_POINTS* r_pts)
{
    real g, exponent, r3, N;
    real reff_xr_rp3;
    unsigned int i;

    /* R2MAG */
    const real gPrime = 5.0 * (mw_log10(coords * 1000.0) - 1.0) + absm;

    /* REFF */
    const real exp_result = mw_exp(sigmoid_curve_params[1] * (gPrime - sigmoid_curve_params[2]));
    const real reff_value = sigmoid_curve_params[0] / (exp_result + 1.0);
    const real rPrime3 = cube(coords);

    for (i = 0; i < n_convolve; ++i)
    {
        g = gPrime + sg[i].dx;

        /* MAG2R */
        r_pts[i].r_in_mag = g;
        r_pts[i].r_in_mag2 = sqr(g);
        r_pts[i].r_point = mw_powr(10.0, (g - absm) / 5.0 + 1.0) / 1000.0;

        r3 = cube(r_pts[i].r_point);
        exponent = sqr(g - gPrime) / (2.0 * sqr(stdev));
        N = ap->coeff * mw_exp(-exponent);
        r_pts[i].qw_r3_N = sg[i].qgaus_W * r3 * N;
    }

    reff_xr_rp3 = reff_value * xr / rPrime3;
    return reff_xr_rp3;
}

