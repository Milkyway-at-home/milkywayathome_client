/*
 *  Copyright (c) 2010, 2011 Matthew Arsenault
 *  Copyright (c) 2010, 2011 Rensselaer Polytechnic Institute
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MILKYWAY_MATH_H_
#define _MILKYWAY_MATH_H_

#define _MILKYWAY_MATH_H_INSIDE_

#define _USE_MATH_DEFINES

#include "milkyway_extra.h"


#ifndef DOUBLEPREC
  #error DOUBLEPREC not defined
#endif

#if HAVE_FLOAT_H
  #include <float.h>
#endif

#ifndef __cplusplus
  #include <math.h>
#else
  #include <cmath>
#endif /* __cplusplus */

/* crlibm is a math.h supplement.*/
#if ENABLE_CRLIBM
  #include <crlibm.h>
#endif /* ENABLE_CRLIBM */


#if DOUBLEPREC
  #include "milkyway_math_double.h"
#else
  #include "milkyway_math_float.h"
#endif

#include "milkyway_math_supplemental.h"
#include "milkyway_vectors.h"


#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#undef _MILKYWAY_MATH_H_INSIDE_

#endif /* _MILKYWAY_MATH_H_ */
