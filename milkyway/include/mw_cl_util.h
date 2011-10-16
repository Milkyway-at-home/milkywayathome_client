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

#ifndef _MW_CL_UTIL_H_
#define _MW_CL_UTIL_H_

#include "mw_cl_setup.h"
#include "mw_cl_types.h"

#ifdef __cplusplus
extern "C" {
#endif


cl_int mwGetWorkGroupInfo(cl_kernel kern, const CLInfo* ci, WGInfo* wgi);
void mwPrintWorkGroupInfo(const WGInfo* wgi);

cl_ulong mwEventTimeNS(cl_event ev);
double mwEventTimeMS(cl_event ev);
double mwEventTime(cl_event ev);

cl_int mwWaitReleaseEvent(cl_event* ev);

#ifdef CL_VERSION_1_1
cl_event mwCreateEvent(CLInfo* ci);
cl_int mwFinishEvent(cl_event ev);
#endif /* CL_VERSION_1_1 */

cl_mem mwCreateZeroReadWriteBuffer(CLInfo* ci, size_t size);
cl_mem mwDuplicateBuffer(CLInfo* ci, cl_mem buf);

/* Print a message with the name of a cl_int error the end of the line */
void mwPerrorCL(cl_int err, const char* fmt, ...);

#ifdef __cplusplus
}
#endif

#endif /* _MW_CL_UTIL_H_ */

