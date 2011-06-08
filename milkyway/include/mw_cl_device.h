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

#if !defined(_MILKYWAY_CL_H_INSIDE_) && !defined(MILKYWAY_CL_COMPILATION)
  #error "Only milkyway_cl.h can be included directly."
#endif


#ifndef _MW_CL_DEVICE_H_
#define _MW_CL_DEVICE_H_

#include "mw_cl_types.h"
#include "milkyway_util.h"


#ifdef __cplusplus
extern "C" {
#endif

cl_int mwSelectDevice(CLInfo* ci, const cl_device_id* devs, const CLRequest* clr, const cl_uint nDev);
cl_int mwGetDevInfo(DevInfo* di, cl_device_id dev);
void mwPrintDevInfoShort(const DevInfo* di);
void mwPrintDevInfo(const DevInfo* di);
void mwPrintPlatforms(cl_platform_id* platforms, cl_uint n_platforms);

MWDoubleExts mwGetDoubleExts(const char* extensions);
cl_bool mwSupportsDoubles(const DevInfo* di);

cl_device_id* mwGetAllDevices(cl_platform_id platform, cl_uint* numDevOut);
cl_platform_id* mwGetAllPlatformIDs(CLInfo* ci, cl_uint* n_platforms_out);


size_t mwFindGroupSize(const DevInfo* di);
cl_uint mwFindGroupsPerCU(const DevInfo* di);
cl_uint mwBlockSize(const DevInfo* di);


/* Nvidia specific functions */
cl_bool minComputeCapabilityCheck(const DevInfo* di, cl_uint major, cl_uint minor);
cl_bool computeCapabilityIs(const DevInfo* di, cl_uint major, cl_uint minor);
cl_uint cudaCoresPerComputeUnit(const DevInfo* di);
cl_double cudaEstimateGFLOPs(const DevInfo* di);
cl_uint cudaEstimateDoubleFrac(const DevInfo* di);

cl_double referenceGFLOPsGTX285(cl_bool doubleprec);
cl_double referenceGFLOPsGTX480(cl_bool doubleprec);
cl_double referenceGFLOPsRadeon5870(cl_bool doubleprec);

#ifdef __cplusplus
}
#endif

#endif /* _MW_CL_DEVICE_H_ */

