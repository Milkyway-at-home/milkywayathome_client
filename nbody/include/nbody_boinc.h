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

#ifndef _NBODY_BOINC_H_
#define _NBODY_BOINC_H_

#include "nbody_config.h"

/* TODO: Naming, and sharing with separation */
#if BOINC_APPLICATION
  #define WIN32_LEAN_AND_MEAN
  #define VC_EXTRALEAN
  #if BOINC_DEBUG
    #include <boinc/diagnostics.h>
  #endif /* BOINC_DEBUG */

  #include <stdio.h>
  #include <filesys.h>
  #include <boinc_api.h>

  #define nbody_finish(x) boinc_finish(x)
  #define nbody_fopen(x,y) boinc_fopen((x),(y))
  #define nbody_remove(x) boinc_delete_file((x))
#else
  #define nbody_finish(x) exit(x)
  #define nbody_fopen(x,y) fopen((x),(y))
  #define nbody_remove(x) remove((x))
#endif /* BOINC_APPLICATION */

#endif /* _NBODY_BOINC_H_ */

