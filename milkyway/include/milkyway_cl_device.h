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

#if !defined(_MILKYWAY_CL_H_INSIDE_) && !defined(MILKYWAY_CL_COMPILATION)
  #error "Only milkyway_cl.h can be included directly."
#endif


#ifndef _MILKYWAY_CL_DEVICE_H_
#define _MILKYWAY_CL_DEVICE_H_

#include "milkyway_cl_types.h"
#include "milkyway_util.h"


#ifdef __cplusplus
extern "C" {
#endif

cl_int mwSelectDevice(CLInfo* ci, const cl_device_id* devs, const CLRequest* clr, const cl_uint nDev);
cl_int mwGetDevInfo(DevInfo* di, cl_device_id dev);
void mwPrintDevInfoShort(const DevInfo* di);
void mwPrintDevInfo(const DevInfo* di);
void mwPrintPlatforms(cl_platform_id* platforms, cl_uint n_platforms);

cl_bool mwPlatformSupportsAMDOfflineDevices(const CLInfo* ci);

MWDoubleExts mwGetDoubleExts(const char* extensions);
cl_bool mwSupportsDoubles(const DevInfo* di);

cl_bool mwDeviceHasDenormals(const DevInfo* di, cl_bool doublePrec);
cl_bool mwDeviceHasFMA(const DevInfo* di, cl_bool doublePrec);
cl_bool mwDeviceHasInfNan(const DevInfo* di, cl_bool doublePrec);
cl_bool mwDeviceHasRTN(const DevInfo* di, cl_bool doublePrec);


cl_device_id* mwGetAllDevices(cl_platform_id platform, cl_uint* numDevOut);
cl_platform_id* mwGetAllPlatformIDs(cl_uint* n_platforms_out);

cl_double mwDeviceEstimateGFLOPs(const DevInfo* di, cl_bool useDouble);

cl_bool mwIsNvidiaGPUDevice(const DevInfo* di);
cl_bool mwIsAMDGPUDevice(const DevInfo* di);


/* AMD specific functions */
cl_double mwAMDEstimateGFLOPs(const DevInfo* di, cl_bool useDouble);
cl_bool mwDeviceVendorIsAMD(const DevInfo* di);
cl_int mwUAVIdFromMWCALtargetEnum(MWCALtargetEnum x);

/* Nvidia specific functions */
cl_bool mwMinComputeCapabilityCheck(const DevInfo* di, cl_uint major, cl_uint minor);
cl_bool mwComputeCapabilityIs(const DevInfo* di, cl_uint major, cl_uint minor);
cl_double mwCUDAEstimateGFLOPs(const DevInfo* di, cl_bool useDouble);
cl_bool mwHasNvidiaCompilerFlags(const DevInfo* di);
cl_bool mwDeviceVendorIsNvidia(const DevInfo* di);
cl_bool mwNvidiaDriverVersionGreaterEqual(const DevInfo* di, cl_uint minMajor, cl_uint minMinor);
cl_bool mwNvidiaInlinePTXAvailable(cl_platform_id platform);

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_CL_DEVICE_H_ */

