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

#ifndef _MILKYWAY_EXTRA_H_
#define _MILKYWAY_EXTRA_H_

#ifndef _WIN32
  #include <inttypes.h>
#endif

/* FIXME: OpenCL bool size */
typedef short int mwbool;

#ifndef TRUE
  #define TRUE  1
#endif
#ifndef FALSE
  #define FALSE 0
#endif

#ifdef _MSC_VER
  #define strdup _strdup
  #define isnan _isnan
  #define isfinite _finite
  #define copysign _copysign
  #define access _access
  #define snprintf _snprintf
  #define getcwd _getcwd
  #define strncasecmp(a, b, n) _strnicmp(a, b, n)
  #define strcasecmp(a, b) _stricmp(a, b)
#endif /* _MSC_VER */

/* Horrible workaround for lack of C99 in MSVCRT and it being
   impossible to print size_t correctly and standardly */
#ifdef _WIN32
  #define ZU "%Iu"
  #define LLU "%I64u"
#else
  #define ZU "%zu"
  #define LLU "%"PRIu64
#endif /* _WIN32 */

#endif /* _MILKYWAY_EXTRA_H_ */

