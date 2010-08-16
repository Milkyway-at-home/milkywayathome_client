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

#ifndef _INTEGRAL_CONSTANTS_H_
#define _INTEGRAL_CONSTANTS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "separation_types.h"

STREAM_CONSTANTS* init_constants(ASTRONOMY_PARAMETERS* ap,
                                 const BACKGROUND_PARAMETERS* bgp,
                                 const STREAMS* streams);

void get_stream_gauss(STREAM_GAUSS* sg, const unsigned int convolve);

double set_prob_consts(const ASTRONOMY_PARAMETERS* ap,
                       const STREAM_GAUSS* sg,
                       const unsigned int n_convolve,
                       const double coords,
                       R_POINTS* r_pts);

void prepare_nu_constants(NU_CONSTANTS* nu_consts,
                          const unsigned int nu_steps,
                          double nu_step_size,
                          double nu_min);

R_CONSTANTS* prepare_r_constants(const ASTRONOMY_PARAMETERS* ap,
                                 const STREAM_GAUSS* sg,
                                 const unsigned int n_convolve,
                                 const unsigned int r_steps,
                                 const double r_min,
                                 const double r_step_size,
                                 const double mu_step_size,
                                 R_POINTS* r_pts);

void prepare_integral_constants(const ASTRONOMY_PARAMETERS* ap,
                                const STREAM_GAUSS* sg,
                                const INTEGRAL_AREA* ia,
                                INTEGRAL_CONSTANTS* ic);

void free_integral_constants(INTEGRAL_CONSTANTS* ic);
void free_stream_gauss(STREAM_GAUSS* sg);


#ifdef __cplusplus
}
#endif

#endif /* _INTEGRAL_CONSTANTS_H_ */

