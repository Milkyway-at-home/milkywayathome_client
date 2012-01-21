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

#ifndef _MILKYWAY_TIMING_H_
#define _MILKYWAY_TIMING_H_

#include "milkyway_config.h"
#include "milkyway_extra.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct
{
    uint64_t sec;
    uint64_t nSec; /* Nanoseconds */
} MWHighResTime;


int mwGetHighResTime_RealTime(MWHighResTime* timer);
MWHighResTime mwDiffMWHighResTime(const MWHighResTime* a, const MWHighResTime* b);
void mwIncMWHighResTime(MWHighResTime* a, const MWHighResTime* b);


double mwGetTime(void);
double mwGetTimeMilli(void);

#ifdef _WIN32
  #define mwMilliSleep(x) Sleep((DWORD) (x))
  /* The usleep() in MinGW tries to round up to avoid sleeping for 0 */
  #define mwMicroSleep(x) Sleep(((DWORD) (x) + 999) / 1000)
#else
  #define mwMilliSleep(x) usleep((useconds_t) 1000 * (x))
  #define mwMicroSleep(x) usleep((useconds_t)(x))
#endif /* _WIN32 */

int mwSetTimerMinResolution(void);
int mwResetTimerResolution(void);


#ifdef __cplusplus
}
#endif


#endif /* _MILKYWAY_UTIL_H_ */

