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

#include "separation_types.h"
#include "milkyway_util.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef real (*ProbabilityFunc)(const AstronomyParameters* ap,
                                const StreamConstants* sc,
                                const real* RESTRICT sg_dx,
                                const real* RESTRICT r_point,
                                const real* RESTRICT qw_r3_N,
                                LBTrig lbt,
                                real gPrime,
                                real reff_xr_rp3,
                                real* RESTRICT streamTmps);


extern ProbabilityFunc probabilityFunc;


int evaluate(SeparationResults* results,
             const AstronomyParameters* ap,
             const IntegralArea* ias,
             const Streams* streams,
             const StreamConstants* sc,
             const char* star_points_file,
             const CLRequest* clr,
             int do_separation,
             int ignoreCheckpoint,
             const char* separation_outfile);

#ifdef __cplusplus
}
#endif

#endif /* _EVALUATION_H_ */

