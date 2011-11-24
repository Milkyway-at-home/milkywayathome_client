/*
 *  Copyright (c) 2011 Matthew Arsenault
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

#ifndef _MILKYWAY_ASPRINTF_H_
#define _MILKYWAY_ASPRINTF_H_

#include "milkyway_config.h"

#ifdef __cplusplus
extern "C" {
#endif


#if defined(_WIN32) && !HAVE_ASPRINTF
int asprintf(char** buf, const char* format, ...);
#endif

#ifdef __cplusplus
}
#endif


#endif /* _MILKYWAY_ASPRINTF_H_ */

