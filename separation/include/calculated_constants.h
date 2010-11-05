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

#ifndef _CALCULATED_CONSTANTS_H_
#define _CALCULATED_CONSTANTS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "separation_types.h"

STREAM_CONSTANTS* getStreamConstants(const ASTRONOMY_PARAMETERS* ap, const STREAMS* streams);

int setAstronomyParameters(ASTRONOMY_PARAMETERS* ap, const BACKGROUND_PARAMETERS* bgp);

STREAM_GAUSS get_stream_gauss(const unsigned int convolve);
void free_stream_gauss(STREAM_GAUSS sg);

NU_CONSTANTS* prepare_nu_constants(const unsigned int nu_steps,
                                   const real nu_step_size,
                                   const real nu_min);


#ifdef __cplusplus
}
#endif

#endif /* _CALCULATED_CONSTANTS_H_ */

