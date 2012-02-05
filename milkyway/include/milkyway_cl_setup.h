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


#ifndef _MILKYWAY_CL_SETUP_H_
#define _MILKYWAY_CL_SETUP_H_

#include "milkyway_cl_types.h"
#include "milkyway_util.h"

#ifdef __cplusplus
extern "C" {
#endif

cl_int mwSetupCL(CLInfo* ci, const CLRequest* clr);
cl_int mwDestroyCLInfo(CLInfo* ci);

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_CL_SETUP_H_ */

