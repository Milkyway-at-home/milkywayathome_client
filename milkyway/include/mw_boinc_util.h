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

#if BOINC_APPLICATION
  #ifndef _WIN32
    /* Workaround: Old version of BOINC libraries missed including this */
    #include <sys/types.h>
  #endif
  #include <boinc/boinc_api.h>
  #include <boinc/filesys.h>
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
  #define mw_begin_critical_section() boinc_begin_critical_section()
  #define mw_end_critical_section() boinc_end_critical_section()
#else
  #define mw_boinc_print(f, msg, ...)
  #define mw_finish(x) exit(x)
  #define mw_fopen(x,y) fopen((x),(y))
  #define mw_remove(x) remove((x))
  #define mw_begin_critical_section()
  #define mw_end_critical_section()
#endif /* BOINC_APPLICATION */

typedef enum
{
    MW_PLAIN       = 1 << 0,
    MW_MULTITHREAD = 1 << 1,
    MW_CAL         = 1 << 2,
    MW_OPENCL      = 1 << 3,
    MW_DEBUG       = 1 << 4,
    MW_GRAPHICS    = 1 << 5
} MWInitType;

int mwBoincInit(MWInitType type);
char* mwReadFileResolved(const char* filename);
FILE* mwOpenResolved(const char* filename, const char* mode);
int mw_rename(const char* oldf, const char* newf);

int mw_resolve_filename(const char* filename, char* buf, size_t bufSize);
int mw_file_exists(const char* file);


#ifdef __cplusplus
}
#endif

#endif /* _MW_BOINC_UTIL_H_ */

