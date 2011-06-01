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
    real bgIntegral;       /* Background integral */
    real* streamIntegrals;
} Cut;

typedef struct
{
    /* State for integral calculation. */
    Cut* cuts;
    Cut* cut;                        /* es->cuts[es->currentCut] */
    unsigned int nu_step, mu_step;   /* r_steps aren't checkpointed */
    Kahan bgSum;
    Kahan* streamSums;

    real bgTmp;
    real* streamTmps;

    unsigned int lastCheckpointNuStep; /* Nu step of last checkpointed (only used by GPU) */
    unsigned int current_calc_probs; /* progress of completed cuts */

    unsigned int currentCut;

    unsigned int numberCuts;
    unsigned int numberStreams;
} EvaluationState;

EvaluationState* newEvaluationState(const AstronomyParameters* ap);
void freeEvaluationState(EvaluationState* es);
void copyEvaluationState(EvaluationState* esDest, const EvaluationState* esSrc);
void clearEvaluationStateTmpSums(EvaluationState* es);
void printEvaluationState(const EvaluationState* es);

void addTmpSums(EvaluationState* es);

int writeCheckpoint(EvaluationState* es);
int readCheckpoint(EvaluationState* es);
int resolveCheckpoint();
int maybeResume(EvaluationState* es);
int timeToCheckpointGPU(const EvaluationState* es, const IntegralArea* ia);

#ifdef __cplusplus
}
#endif


#endif /* _EVALUATION_STATE_H_ */

