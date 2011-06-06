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

#ifndef _MILKYWAY_CPP_UTIL_H_
#define _MILKYWAY_CPP_UTIL_H_


#if BOINC_APPLICATION
  #include <boinc_api.h>
#endif /* BOINC_APPLICATION */

#ifdef __cplusplus
extern "C" {
#endif

/* Stuff that might be actually useful */
typedef struct
{
    int majorVersion;
    int minorVersion;
    int release;
    int appVersion;
    double checkpointPeriod;
    char wuName[256];
    char projectDir[256];
    char boincDir[256];
    char projectPrefs[4096];
} MWAppInitData;

#if BOINC_APPLICATION
int mwGetMWAppInitData(MWAppInitData* mwaid);
void mwGetBoincOptionsDefault(BOINC_OPTIONS* options);
void* mw_graphics_make_shmem(const char* x, int y);
void* mw_graphics_get_shmem(const char* x);

#endif /* BOINC_APPLICATION */

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_CPP_UTIL_H_ */

