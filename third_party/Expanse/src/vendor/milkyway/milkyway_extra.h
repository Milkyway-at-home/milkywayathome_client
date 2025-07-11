/*
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
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

#ifndef _MILKYWAY_EXTRA_H_
#define _MILKYWAY_EXTRA_H_

#include "milkyway_config.h"

#if HAVE_INTTYPES_H
  #include <inttypes.h>
#endif

#if HAVE_STDINT_H
  #include <stdint.h>
#endif

#if HAVE_FLOAT_H
  #include <float.h>
#endif

#if HAVE_PROCESS_H
  #include <process.h>
#endif

typedef char mwbool;

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
  #define getpid _getpid
  #define open _open
  #define close _close
  #define read _read
  #define lseek _lseek
  #define mkdir _mkdir

  #if _MSC_VER > 1600
    #define write _write
  #endif

  #define isinf(x) (!_finite(x) && !_isnan(x))
#endif /* _MSC_VER */

/* Horrible workaround for lack of C99 in MSVCRT and it being
   impossible to print size_t correctly and standardly */
#ifdef _WIN32
  #define ZU "%Iu"
  #define LLU "%I64u"
#else
  #define ZU "%zu"
  #define LLU "%" PRIu64
#endif /* _WIN32 */

#endif /* _MILKYWAY_EXTRA_H_ */
