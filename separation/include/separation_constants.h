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

#ifndef _SEPARATION_CONSTANTS_H_
#define _SEPARATION_CONSTANTS_H_

#if DOUBLEPREC
  #define EPS (3.0e-11)
#else
  #define EPS (5.96e-08)
#endif /* DOUBLEPREC */

#define stdev (0.6)
#define xr (3.0 * stdev)
#define absm (4.2)
#define SIGMA_LIMIT (0.0001)

#define stripeSeparation  (2.5)
#define surveyCenterRa (185.0)
#define surveyCenterDec (32.5)

/* The node of the GC coordinates used in the survey. */
#define NODE_GC_COORDS (surveyCenterRa - 90.0)
#define NODE_GC_COORDS_RAD d2r(NODE_GC_COORDS)

#define sun_r0 (8.5)

#define CHECKPOINT_FILE "separation_checkpoint"
#define CHECKPOINT_FILE_TMP "separation_checkpoint_tmp"

#endif /* _SEPARATION_CONSTANTS_H_ */

