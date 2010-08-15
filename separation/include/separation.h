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


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "milkyway_util.h"
#include "evaluation.h"
#include "parameters.h"
#include "probability.h"
#include "star_points.h"
#include "integral_constants.h"
#include "numerical_integration.h"
#include "../util/io_util.h"
#include "coordinates.h"
#include "integrals_likelihood.h"


#if BOINC_APPLICATION
  #include <boinc_api.h>
  #include <filesys.h>
#endif

#define EPS 3.0e-11

#define dmod(A,B) ((B) != 0.0 ? ((A)*(B) > 0.0 ? (A) - (B) * floor((A)/(B))\
                             :(A) + (B)*floor(-(A)/(B))):(A))
#define dsign(A,B) ((B) < 0.0 ? -(A) : (A))

#define stripeSeparation  (2.5)
#define surveyCenterRa (185.0)
#define surveyCenterDec (32.5)

/* The node of the GC coordinates used in the survey. */
#define NODE_GC_COORDS (surveyCenterRa - 90.0)
#define sun_r0 8.5

#define CHECKPOINT_FILE "separation_checkpoint"
#define CHECKPOINT_FILE_TMP "separation_checkpoint_tmp"

#endif /* _SEPARATION_H_ */

