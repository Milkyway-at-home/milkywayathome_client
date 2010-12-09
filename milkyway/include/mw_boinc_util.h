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

#ifndef _MW_BOINC_UTIL_H_
#define _MW_BOINC_UTIL_H_

#include "milkyway_config.h"

#include <stdio.h>

#if BOINC_APPLICATION
  #include <boinc/boinc_api.h>
  #include <boinc/filesys.h>

  #if BOINC_APP_GRAPHICS
    #include <boinc/graphics_api.h>
    #include <boinc/graphics_lib.h>
  #endif /* BOINC_APP_GRAPHICS */
#endif /* BOINC_APPLICATION */

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#if BOINC_APPLICATION
  #define mw_boinc_print(f, msg, ...) fprintf(f, msg, ##__VA_ARGS__)
  #define mw_finish(x) boinc_finish(x)
  #define mw_fopen(x,y) boinc_fopen((x),(y))
  #define mw_remove(x) boinc_delete_file((x))
  #define mw_rename(x, y) boinc_rename((x), (y))
#else
  #define mw_boinc_print(f, msg, ...)
  #define mw_finish(x) exit(x)
  #define mw_fopen(x,y) fopen((x),(y))
  #define mw_remove(x) remove((x))
  #define mw_rename(x, y) rename((x), (y))
#endif /* BOINC_APPLICATION */

int mwBoincInit(const char* appname, int useDebug);
char* mwReadFileResolved(const char* filename);
FILE* mwOpenResolved(const char* filename, const char* mode);
int mwRename(const char* oldf, const char* newf);


#ifdef __cplusplus
}
#endif

#endif /* _MW_BOINC_UTIL_H_ */

