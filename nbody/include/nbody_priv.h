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

#ifndef _NBODY_PRIV_H_
#define _NBODY_PRIV_H_

#define _GNU_SOURCE

#ifdef __cplusplus
extern "C" {
#endif

#include "nbody_config.h"
#include "nbody_types.h"
#include "stdinc.h"
#include "vectmath.h"
#include "real.h"
#include "nbody_util.h"
#include "show.h"
#include "io.h"
#include "grav.h"
#include "chisq.h"
#include "load.h"
#include "orbitintegrator.h"
#include "plummer.h"

#if BOINC_APPLICATION
  #include <boinc_api.h>
  #if BOINC_DEBUG
    #include <diagnostics.h>
  #endif /* BOINC_DEBUG */
#endif /* BOINC_APPLICATION */

#endif /* _NBODY_PRIV_H_ */

#include <stdlib.h>

#ifdef __cplusplus
}
#endif
