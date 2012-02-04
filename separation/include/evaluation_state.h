/*
 *  Copyright (c) 2008-2010 Travis Desell, Nathan Cole, Dave Przybylo
 *  Copyright (c) 2008-2010 Boleslaw Szymanski, Heidi Newberg
 *  Copyright (c) 2008-2010 Carlos Varela, Malik Magdon-Ismail
 *  Copyright (c) 2008-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _EVALUATION_STATE_H_
#define _EVALUATION_STATE_H_

#include "separation_types.h"

#ifdef __cplusplus
extern "C" {
#endif

EvaluationState* newEvaluationState(const AstronomyParameters* ap);
void freeEvaluationState(EvaluationState* es);
void copyEvaluationState(EvaluationState* esDest, const EvaluationState* esSrc);
void clearEvaluationStateTmpSums(EvaluationState* es);
void printEvaluationState(const EvaluationState* es);
int integralsAreDone(const EvaluationState* es);

void addTmpCheckpointSums(EvaluationState* es);
int writeCheckpoint(EvaluationState* es);
int readCheckpoint(EvaluationState* es);
int resolveCheckpoint(void);
int maybeResume(EvaluationState* es);
int timeToCheckpointGPU(const EvaluationState* es, const IntegralArea* ia);

#ifdef __cplusplus
}
#endif


#endif /* _EVALUATION_STATE_H_ */

