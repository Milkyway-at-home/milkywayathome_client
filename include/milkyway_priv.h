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

#ifndef _MILKYWAY_PRIV_H_
#define _MILKYWAY_PRIV_H_

#include <stdio.h>
#include <assert.h>

#include "config.h"

#ifdef MW_ENABLE_DEBUG
    /* convenient functions for printing debugging stuffs */
    #define MW_DEBUG(msg, ...) printf("%s():%d: ", __func__, __LINE__);\
                                 printf(msg, __VA_ARGS__);
    #define MW_DEBUGMSG(msg) puts(msg)

#else
    #define MW_DEBUG(msg, ...) ((void) 0)
    #define MW_DEBUGMSG(msg, ...) ((void) 0)
#endif

#if BOINC_APPLICATION
  #define mw_finish(x) boinc_finish(x)
#else
  #define mw_finish(x) exit(x)
#endif /* BOINC_APPLICATION */

#endif /* _MILKYWAY_PRIV_H_ */

