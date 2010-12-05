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

#include "separation_types.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Completed integral state */
typedef struct
{
    real background_integral;
    real* stream_integrals;
    Kahan* probs;
} Integral;

typedef struct
{
    /* State for integral calculation. */
    Integral* integrals;
    unsigned int nu_step, mu_step;   /* r_steps aren't checkpointed */
    Kahan sum;

    unsigned int current_calc_probs; /* progress of completed integrals */

    unsigned int current_integral;

    unsigned int number_integrals;
    unsigned int number_streams;
} EvaluationState;

EvaluationState* newEvaluationState(const AstronomyParameters* ap);
void freeEvaluationState(EvaluationState* es);
void copyEvaluationState(EvaluationState* esDest, const EvaluationState* esSrc);

int writeCheckpoint(const EvaluationState* es);
int readCheckpoint(EvaluationState* es);
int resolveCheckpoint();
int maybeResume(EvaluationState* es);

#ifdef __cplusplus
}
#endif


#endif /* _EVALUATION_STATE_H_ */

