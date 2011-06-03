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

#ifndef _MW_ASPRINTF_H_
#define _MW_ASPRINTF_H_

#include <stdarg.h>
#include "milkyway_config.h"

#ifdef __cplusplus
extern "C" {
#endif


#if defined(_WIN32) && !HAVE_ASPRINTF
int _mw_asprintf(char** buf, const char* format, ...);
#define asprintf(buf, fmt, ...) _mw_asprintf(buf, fmt, ##__VA_ARGS__)

#endif

#ifdef __cplusplus
}
#endif


#endif /* _MW_ASPRINTF_H_ */

