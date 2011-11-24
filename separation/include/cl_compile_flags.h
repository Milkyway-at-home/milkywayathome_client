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

#ifndef _CL_COMPILER_FLAGS_H_
#define _CL_COMPILER_FLAGS_H_

#include "milkyway_cl.h"
#include "separation_types.h"


#ifdef __cplusplus
extern "C" {
#endif

char* getCompilerFlags(const CLInfo* ci, const AstronomyParameters* ap, cl_bool useILernel);

#ifdef __cplusplus
}
#endif

#endif /* _CL_COMPILER_FLAGS_H_ */

