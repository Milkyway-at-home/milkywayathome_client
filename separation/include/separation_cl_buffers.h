/*
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
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

#ifndef _SEPARATION_CL_BUFFERS_H_
#define _SEPARATION_CL_BUFFERS_H_

#include "milkyway_cl.h"
#include "mw_cl.h"
#include "separation_types.h"
#include "setup_cl.h"


#ifdef __cplusplus
extern "C" {
#endif

cl_int createSeparationBuffers(CLInfo* ci,
                               SeparationCLMem* cm,
                               const AstronomyParameters* ap,
                               const IntegralArea* ia,
                               const StreamConstants* sc,
                               const StreamGauss sg,
                               const SeparationSizes* sizes);

void releaseSeparationBuffers(SeparationCLMem* cm);

void calculateSizes(SeparationSizes* sizes, const AstronomyParameters* ap, const IntegralArea* ia);


#ifdef __cplusplus
}
#endif

#endif /* _SEPARATION_CL_BUFFERS_H_ */

