/*
Copyright (C) 2011  Matthew Arsenault

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

#include <stdlib.h>
#include <stdio.h>

#include "mw_asprintf.h"

#ifndef _WIN32
  /* Random internet person says this should work but it seems to not */
  #define _vscprintf(f,a) vsnprintf(NULL, (size_t) 0, f, a)
#endif

#if defined(_WIN32) && !HAVE_ASPRINTF

/* Intended for use on windows where asprintf is missing */
int _mw_asprintf(char** bufOut, const char* format, ...)
{
    int size, rc;
    char* buf;
    va_list args;

    va_start(args, format);

    size = _vscprintf(format, args);
    if (size < 0)
    {
        *bufOut = NULL;
        return size;
    }

    buf = (char*) malloc((size + 1) * sizeof(char));  /* Include null terminator */
    if (!buf)
    {
        *bufOut = NULL;
        return -1;
    }

    rc = vsprintf(buf, format, args);
    va_end(args);

    *bufOut = buf;
    return rc;
}

#endif /* defined(_WIN32) && !HAVE_ASPRINTF */

