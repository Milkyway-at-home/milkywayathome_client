/*
Copyright (C) 2010  Matthew Arsenault

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

#ifndef _SHOW_CAL_TYPES_H_
#define _SHOW_CAL_TYPES_H_

#include <CAL/cal.h>
#include <CAL/calcl.h>

#ifdef __cplusplus
extern "C" {
#endif

const char* showCALboolean(CALboolean x);
const char* showCALresult(CALresult x);
const char* showCALtargetEnum(enum CALtargetEnum x);

#ifdef __cplusplus
}
#endif

#endif /* _SHOW_CAL_TYPES_H_ */

