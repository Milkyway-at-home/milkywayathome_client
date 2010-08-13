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


#ifndef _EVALUATION_STATE_H_
#define _EVALUATION_STATE_H_

#include <stdio.h>
#include "parameters.h"
#include "star_points.h"

typedef struct
{
    /* State for integral calculation. */
    INTEGRAL* integrals;
    double* stream_integrals;
    unsigned int number_integrals;
    unsigned int current_integral;
    unsigned int number_streams;

    double background_integral;
} EVALUATION_STATE;

#define EMPTY_EVALUATION_STATE { NULL, NULL, 0, 0, 0, 0.0 }

void initialize_state(const ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es);
void free_evaluation_state(EVALUATION_STATE* es);
void reset_evaluation_state(EVALUATION_STATE* es);

int write_checkpoint(EVALUATION_STATE* es);
int read_checkpoint(EVALUATION_STATE* es);


#endif /* _EVALUATION_STATE_H_ */

