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

#ifndef _SEPARATION_CAL_SETUP_H_
#define _SEPARATION_CAL_SETUP_H_

#include "milkyway_util.h"
#include "separation_types.h"
#include "evaluation_state.h"
#include "separation_cal_types.h"

#include <CAL/cal.h>
#include <CAL/cal_ext.h>

#ifdef __cplusplus
extern "C" {
#endif

#define cal_warn(msg, err, ...) fprintf(stderr, msg ": %s (%s)\n", ##__VA_ARGS__, calGetErrorString(), showCALresult(err))

CALresult separationCALInit(MWCALInfo* ci, const CLRequest* clr);

CALresult separationLoadKernel(MWCALInfo* ci,
                               const AstronomyParameters* ap,
                               const StreamConstants* sc);


CALresult createConstantBuffer1D(MWMemRes* mr,
                                 MWCALInfo* ci,
                                 const CALvoid* src,
                                 CALformat format,
                                 CALuint width);

CALresult mwUnloadKernel(MWCALInfo* ci);
CALresult mwCALShutdown(MWCALInfo* ci);

CALresult getModuleNames(MWCALInfo* ci, SeparationCALNames* cn, CALuint numberStreams);
void destroyModuleNames(SeparationCALNames* cn);

CALresult mapMWMemRes(MWMemRes* mr, CALvoid** pPtr, CALuint* pitch);
CALresult unmapMWMemRes(MWMemRes* mr);

CALresult setKernelArguments(MWCALInfo* ci, SeparationCALMem* cm, SeparationCALNames* cn);

CALresult createSeparationBuffers(MWCALInfo* ci,
                                  SeparationCALMem* cm,
                                  const AstronomyParameters* ap,
                                  const IntegralArea* ia,
                                  const StreamGauss sg,
                                  const CALSeparationSizes* sizes);

CALresult releaseSeparationBuffers(MWCALInfo* ci, SeparationCALMem* cm);

typedef CALresult (CALAPIENTRYP PFNCALCTXWAITFOREVENTS)(CALcontext ctx, CALevent *event, CALuint num, CALuint flags);
CALAPI PFNCALCTXWAITFOREVENTS mw_calCtxWaitForEvents;


#ifdef __cplusplus
}
#endif

#endif /* _SEPARATION_CAL_SETUP_H_ */

