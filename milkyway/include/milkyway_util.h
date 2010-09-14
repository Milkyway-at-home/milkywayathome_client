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

#ifndef _MILKYWAY_UTIL_H_
#define _MILKYWAY_UTIL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "milkyway_config.h"

#include <stdlib.h>
#include <stdio.h>

#if BOINC_APPLICATION
  #include <boinc_api.h>
  #include <filesys.h>

  #define WIN32_LEAN_AND_MEAN
  #define VC_EXTRALEAN
  #if BOINC_DEBUG
    #include <boinc/diagnostics.h>
  #endif /* BOINC_DEBUG */
#endif /* BOINC_APPLICATION */

#ifdef MW_ENABLE_DEBUG
    /* convenient functions for printing debugging stuffs */
    #define MW_DEBUG(msg, ...) printf("%s():%d: ", __func__, __LINE__);\
                                 printf(msg, __VA_ARGS__);
    #define MW_DEBUGMSG(msg) puts(msg)
#else
    #define MW_DEBUG(msg, ...) ((void) 0)
    #define MW_DEBUGMSG(msg, ...) ((void) 0)
#endif

#if BOINC_APPLICATION
  #define mw_finish(x) boinc_finish(x)
  #define mw_fopen(x,y) boinc_fopen((x),(y))
  #define mw_remove(x) boinc_delete_file((x))
  #define mw_rename(x, y) boinc_rename((x), (y))
#else
  #define mw_finish(x) exit(x)
  #define mw_fopen(x,y) fopen((x),(y))
  #define mw_remove(x) remove((x))
  #define mw_rename(x, y) rename((x), (y))
#endif /* BOINC_APPLICATION */

void* mallocSafe(size_t size);
void* callocSafe(size_t count, size_t size);


#define warn(msg, ...) fprintf(stderr, msg, ##__VA_ARGS__)
#define fail(msg, ...) { fprintf(stderr, msg, ##__VA_ARGS__);  \
                         mw_finish(EXIT_FAILURE); }

char* mwReadFile(const char* filename);
FILE* mwOpenResolved(const char* filename, const char* mode);

int mwRename(const char* oldf, const char* newf);

double get_time();

#if defined(__SSE__) && DISABLE_DENORMALS
int mwDisableDenormalsSSE();
#endif /* defined(__SSE__) && DISABLE_DENORMALS */

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_UTIL_H_ */

