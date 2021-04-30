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

#ifndef _MILKYWAY_CL_PROGRAM_H_
#define _MILKYWAY_CL_PROGRAM_H_

#include "milkyway_cl_types.h"

#ifdef __cplusplus
extern "C" {
#endif

unsigned char* mwGetProgramBinary(cl_program program, size_t* binSizeOut);
int mwSaveProgramBinaryToFile(cl_program program, const char* filename);

cl_program mwCreateProgramFromBin(CLInfo* ci, const unsigned char* bin, size_t binSize);
cl_program mwCreateProgramFromSrc(CLInfo* ci,
                                  cl_uint srcCount,
                                  const char** src,
                                  const size_t* lengths,
                                  const char* compileDefs);
cl_int mwBuildProgram(cl_program program, cl_device_id device, const char* options);

cl_kernel mwCreateKernel(cl_program program, const char* name);

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_CL_PROGRAM_H_ */

