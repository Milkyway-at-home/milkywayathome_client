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

#ifndef _PROBABILITIES_H_
#define _PROBABILITIES_H_

#include "separation_types.h"
#include "evaluation_state.h"
#include "evaluation.h"

#ifdef __cplusplus
extern "C" {
#endif

/* probabilities will be rebuilt for each SSE level */
#if MW_IS_X86
  #if defined(__SSE3__)
    #define INIT_PROBABILITIES initProbabilities_SSE3
  #elif defined(__SSE2__)
    #define INIT_PROBABILITIES initProbabilities_SSE2
  #else
    #define INIT_PROBABILITIES initProbabilities
  #endif
#else
  #define INIT_PROBABILITIES initProbabilities
#endif /* MW_IS_X86 */


#define DEFINE_INIT_PROBABILITIES(level) ProbabilityFunc initProbabilities##level(const AstronomyParameters* ap, int useIntrinsics)

#if MW_IS_X86
DEFINE_INIT_PROBABILITIES();
DEFINE_INIT_PROBABILITIES(_SSE2);
DEFINE_INIT_PROBABILITIES(_SSE3);
#endif /* MW_IS_X86 */


#ifdef __cplusplus
}
#endif

#endif /* _PROBABILITIES_H_ */

