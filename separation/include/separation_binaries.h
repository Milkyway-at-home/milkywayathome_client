/*
 *  Copyright (c) 2010 Matthew Arsenault
 *  Copyright (c) 2010 Rensselaer Polytechnic Institute.
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

#ifndef _SEPARATION_BINARIES_H_
#define _SEPARATION_BINARIES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "milkyway_cl.h"

typedef struct
{
    cl_uint versionMajor;   /* Separation version information */
    cl_uint versionMinor;
    cl_int doublePrec;

    cl_int number_streams;   /* Constants compiled into kernel */
    cl_int fast_h_prob;
    cl_int aux_bg_profile;
    cl_int zero_q;

    size_t binSize;          /* Size of program binary */

    cl_device_type devType;  /* Check device/platform/driver versions */
    cl_uint vendorID;
    char devName[128];
    char deviceVersion[128];
    char driverVersion[128];
    char _reserved[512];
} SeparationBinaryHeader;

unsigned char* separationLoadBinary(const AstronomyParameters* ap,
                                    const DevInfo* di,
                                    const char* filename,
                                    size_t* binSizeOut);


cl_bool separationSaveBinary(const AstronomyParameters* ap,
                             const DevInfo* di,
                             const unsigned char* bin,
                             const size_t binSize,
                             const char* filename);

#ifdef __cplusplus
}
#endif

#endif /* _SEPARATION_BINRIES_H_ */

