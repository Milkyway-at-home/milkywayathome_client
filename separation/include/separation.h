/* Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

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

#ifndef _SEPARATION_H_
#define _SEPARATION_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "separation_config.h"
#include "milkyway_math.h"
#include "milkyway_util.h"
#include "separation_types.h"
#include "evaluation.h"
#include "parameters.h"
#include "star_points.h"
#include "calculated_constants.h"
#include "separation_constants.h"
#include "r_points.h"
#include "gauss_legendre.h"
#include "io_util.h"
#include "coordinates.h"
#include "integrals.h"
#include "likelihood.h"

#if BOINC_APPLICATION
  #include <boinc_api.h>
  #include <filesys.h>

  #if BOINC_DEBUG
    #include <diagnostics.h>
  #endif /* BOINC_DEBUG */

  #if BOINC_APP_GRAPHICS
    #include <graphics_api.h>
    #include <graphics_lib.h>
  #endif /* BOINC_APP_GRAPHICS */

  #ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
    #define VC_EXTRALEAN
    #include <windows.h>
  #endif /* _WIN32 */
#endif /* BOINC_APPLICATION */

#ifdef __cplusplus
}
#endif


#endif /* _SEPARATION_H_ */

