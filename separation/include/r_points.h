/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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

#ifndef _R_POINTS_H_
#define _R_POINTS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "separation_types.h"
#include "separation_constants.h"
#include "milkyway_cl.h"
#include "milkyway_extra.h"


/* Literals are assumed to be doubles by default, and the
 * -cl-single-precision-constant flag seems to not be working when
 * trying to use float */
#define R1 ((real) 1.0)
#define R2 ((real) 2.0)
#define R3 ((real) 3.0)
#define R5 ((real) 5.0)
#define R10 ((real) 10.0)
#define R1000 ((real) 1000.0)


inline real calcGPrime(const real coords)
{
    return R5 * (mw_log10(coords * R1000) - R1) + absm;
}

inline real calcReffXrRp3(const real coords)
{
    _MW_STATIC const real sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };
    const real gPrime = calcGPrime(coords);

    /* REFF */
    const real exp_result = mw_exp(sigmoid_curve_params[1] * (gPrime - sigmoid_curve_params[2]));
    const real reff_value = sigmoid_curve_params[0] / (exp_result + R1);
    const real rPrime3 = cube(coords);
    const real reff_xr_rp3 = reff_value * xr / rPrime3;
    return reff_xr_rp3;
}

R_PRIME calcRPrime(__MW_CONSTANT INTEGRAL_AREA* ia, const unsigned int r_step);

void set_r_points(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                  __MW_CONSTANT STREAM_GAUSS* sg,
                  const unsigned int n_convolve,
                  const real coords,
                  __MW_LOCAL R_POINTS* r_pts);

#ifdef __cplusplus
}
#endif

#endif /* _INTEGRAL_CONSTANTS_H_ */

