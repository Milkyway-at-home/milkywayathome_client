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

#ifndef _MILKYWAY_AT_HOME_
#define _MILKYWAY_AT_HOME_


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"

#include "evaluation_optimized.h"
#include "simple_evaluator.h"
#include "parameters.h"
#include "probability.h"
#include "atSurveyGeometry.h"
#include "star_points.h"
#include "numericalIntegration.h"
#include "../util/io_util.h"

#if BOINC_APPLICATION
  #include <boinc_api.h>
  #include <filesys.h>
#endif


#ifndef M_PI
	#define pi 3.1415926535897932384626433832795028841971693993751
#else
	#define pi M_PI
#endif

#define deg (180.0/pi)

#define EPS 3.0e-11
#define PI (double) 3.1415926535897932384626433832795028841971693993751
#define D2PI (double) 6.2831853071795864769252867665590057683943387987502
#define DPI (double) 3.1415926535897932384626433832795028841971693993751

#define PI_4_3 4.18879020478639098461685784437267051226289253250014109463325945641042\
1875048278664837379767122822757

#define PI_2_3 2.09439510239319549230842892218633525613144626625007054731662972820521\
0937524139332418689883561411379

#endif /* _MILKYWAY_AT_HOME_ */

