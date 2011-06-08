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

#ifndef _MW_CL_H_
#define _MW_CL_H_

#include "milkyway_config.h"

#ifdef __APPLE__
  #include <OpenCL/cl.h>
  #include <OpenCL/cl_platform.h>
  #include <OpenCL/cl_ext.h>
#else
  #include <CL/cl.h>
  #include <CL/cl_platform.h>
  #include <CL/cl_ext.h>
#endif /* __APPLE__ */

#endif /* _MW_CL_H_ */

