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

#ifndef _EVALUATION_H_
#define _EVALUATION_H_

#include "parameters.h"
#include "star_points.h"
#include "integral_constants.h"

typedef struct
{
    double st_prob_int;    /* for Kahan summation */
    double st_prob_int_c;
} ST_PROBS;

typedef struct
{
    double bg_int;
    double correction;   /* Correction for Kahan summation */
} BG_PROB;

/* TODO: All these tuples of doubles really serve the same
 * purpose. Fix having all of them. */
typedef struct
{
    double sum;
    double correction;
} PROB_SUM;

#define ZERO_PROB_SUM { 0.0, 0.0 }

#define CLEAR_BG_PROB(bgp) { (bgp).bg_int = 0.0; (bgp).correction = 0.0; }

/* Add b to a */
#define INCADD_BG_PROB(a, b) { (a).bg_int += (b).bg_int; (a).correction += (b).correction; }

#define ZERO_BG_PROB { 0.0, 0.0 }

typedef struct
{
    double st_only_sum;
    double st_only_sum_c;
} ST_SUM;

#define KAHAN_ADD(sum, item, correction)                                \
    {                                                                   \
        double _tmp = sum;                                              \
        sum += item;                                                    \
        correction +=  item - (sum - _tmp);                             \
    }

double evaluate(const ASTRONOMY_PARAMETERS* ap,
                const STAR_POINTS* sp,
                const STREAMS* streams,
                const STREAM_CONSTANTS* sc);

#endif /* _EVALUATION_H_ */

