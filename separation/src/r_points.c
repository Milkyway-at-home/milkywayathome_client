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

#include "milkyway_util.h"

#ifndef NDEBUG
extern real calcReffXrRp3(const real coords, const real gPrime);
extern real calcGPrime(const real coords);
extern R_PRIME calcRPrime(__MW_CONSTANT INTEGRAL_AREA* ia, const unsigned int r_step);
extern R_POINTS calc_r_point(__MW_CONSTANT STREAM_GAUSS* sg, const real gPrime, const real coeff);
#endif


void set_r_points(const ASTRONOMY_PARAMETERS* ap,
                  const STREAM_GAUSS* sg,
                  const unsigned int n_convolve,
                  const real gPrime,
                  R_POINTS* r_pts)
{
    unsigned int i;

    for (i = 0; i < n_convolve; ++i)
        r_pts[i] = calc_r_point(&sg[i], gPrime, ap->coeff);
}

R_POINTS* precalculate_r_pts(const ASTRONOMY_PARAMETERS* ap,
                             const INTEGRAL_AREA* ia,
                             const STREAM_GAUSS* sg,
                             R_CONSTS** rc_out)
{
    unsigned int i;
    R_POINTS* r_pts;
    R_PRIME rp;
    R_CONSTS* rc;

    size_t rPtsSize = sizeof(R_POINTS) * ap->convolve * ia->r_steps;
    size_t rConstsSize = sizeof(R_CONSTS) * ia->r_steps;

    r_pts = (R_POINTS*) mwMallocAligned(rPtsSize, sizeof(R_POINTS));
    rc = (R_CONSTS*) mwMallocAligned(rConstsSize, sizeof(R_CONSTS));

    for (i = 0; i < ia->r_steps; ++i)
    {
        rp = calcRPrime(ia, i);
        rc[i] = calcRConsts(rp);
        set_r_points(ap, sg, ap->convolve, rc[i].gPrime, &r_pts[i * ap->convolve]);
    }

    *rc_out = rc;
    return r_pts;
}

