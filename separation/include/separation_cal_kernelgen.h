/*
Copyright (C) 2010  Matthew Arsenault

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

#ifndef _SEPARATION_CAL_KERNELGEN_H_
#define _SEPARATION_CAL_KERNELGEN_H_

#include "separation_types.h"
#include "evaluation_state.h"
#include "milkyway_util.h"


#ifdef __cplusplus
extern "C" {
#endif

char* separationIntegralKernelSrc(const AstronomyParameters* ap,
                                  const StreamConstants* sc,
                                  CALuint device);

char* separationLikelihoodKernelSrc(const AstronomyParameters* ap,
                                    const StreamConstants* sc,
                                    CALuint device);


#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
std::string createSeparationIntegralKernel(const AstronomyParameters* ap,
                                           const StreamConstants* sc,
                                           CALuint device);

std::string createSeparationLikelihoodKernel(const AstronomyParameters* ap,
                                             const StreamConstants* sc,
                                             CALuint device);

#endif




#endif /* _SEPARATION_CAL_KERNELGEN_H_ */

