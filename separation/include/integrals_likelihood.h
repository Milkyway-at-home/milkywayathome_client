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

#ifndef _INTEGRALS_LIKELIHOOD_H_
#define _INTEGRALS_LIKELIHOOD_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "separation_constants.h"
#include "separation_types.h"

#define lbr2xyz_2(xyz, r_point, bsin, bcos, lsin, lcos) \
    {                                                   \
        real zp = r_point * bcos;                       \
        X(xyz) = zp * lcos - sun_r0;                    \
        Y(xyz) = zp * lsin;                             \
        Z(xyz) = r_point * bsin;                        \
    }

ALWAYS_INLINE HOT
inline real calc_st_prob_inc(__MW_CONSTANT STREAM_CONSTANTS* sc, const vector xyz, const real qw_r3_N)
{
    vector xyzs;
    real xyz_norm, dotted;

    SUBV(xyzs, xyz, sc->c);
    DOTVP(dotted, sc->a, xyzs);
    INCSUBVMS(xyzs, dotted, sc->a);
    SQRV(xyz_norm, xyzs);

    return qw_r3_N * mw_exp(-xyz_norm * sc->sigma_sq2_inv);
}


#ifdef __cplusplus
}
#endif

#endif /* _INTEGRALS_LIKELIHOOD_H_ */

