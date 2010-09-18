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

#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include "separation_types.h"

void free_background_parameters(BACKGROUND_PARAMETERS* bgp);
void free_streams(STREAMS* streams);
void free_stream_parameters(STREAM_PARAMETERS* p);

unsigned int get_optimized_parameter_count(ASTRONOMY_PARAMETERS* ap,
                                           BACKGROUND_PARAMETERS* bgp,
                                           STREAMS* streams);

INTEGRAL_AREA* read_parameters(const char* file,
                               ASTRONOMY_PARAMETERS* ap,
                               BACKGROUND_PARAMETERS* bgp,
                               STREAMS* streams);

int write_parameters(const char* file,
                     ASTRONOMY_PARAMETERS* ap,
                     INTEGRAL_AREA* integral,
                     BACKGROUND_PARAMETERS* bgp,
                     STREAMS* streams);

void set_parameters(ASTRONOMY_PARAMETERS* ap,
                    BACKGROUND_PARAMETERS* bgp,
                    STREAMS* streams,
                    const real* parameters);


#endif /* _PARAMETERS_H_ */

