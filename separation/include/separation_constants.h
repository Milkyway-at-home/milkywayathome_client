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
  #define EPS ((real) 3.0e-11)
#else
  #define EPS ((real) 5.96e-08)
#endif /* DOUBLEPREC */

#if DOUBLEPREC
  #define SEPARATION_EPS (1.0e-15)
#else
  #define SEPARATION_EPS (1.0e-7)
#endif /* DOUBLEPREC */

#define stdev ((real) 0.6)
#define xr ((real) 3.0 * stdev)
#define absm ((real) 4.2)
#define SIGMA_LIMIT ((real) 0.0001)

#define stripeSeparation  ((real) 2.5)
#define surveyCenterRa ((real) 185.0)
#define surveyCenterDec ((real) 32.5)
#define surveyCenterDec_rad (d2r(surveyCenterDec))

/* The node of the GC coordinates used in the survey. */
#define NODE_GC_COORDS (surveyCenterRa - (real) 90.0)
#define NODE_GC_COORDS_RAD ((real) d2r(NODE_GC_COORDS))

#define const_sun_r0 ((real) 8.5)

#define CHECKPOINT_FILE "separation_checkpoint"
#define CHECKPOINT_FILE_TMP "separation_checkpoint_tmp"

#define MAX_CONVOLVE 256

#endif /* _SEPARATION_CONSTANTS_H_ */

