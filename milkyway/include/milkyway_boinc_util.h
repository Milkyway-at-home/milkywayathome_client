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

#ifndef _MILKYWAY_BOINC_UTIL_H_
#define _MILKYWAY_BOINC_UTIL_H_

#include "milkyway_config.h"


#if BOINC_APPLICATION
  #if HAVE_SYS_TYPES_H
    /* Workaround: Old version of BOINC libraries missed including this */
    #include <sys/types.h>
  #endif
  #include <boinc_api.h>
  #include <filesys.h>
#endif /* BOINC_APPLICATION */

#include <stdio.h>



typedef enum
{
    MW_PLAIN       = 0,
    MW_MULTITHREAD = 1 << 0,
    MW_CAL         = 1 << 1,
    MW_OPENCL      = 1 << 2,
    MW_DEBUG       = 1 << 3,
    MW_GRAPHICS    = 1 << 4
} MWInitType;

typedef enum
{
    MW_PREF_NONE,
    MW_PREF_DOUBLE,
    MW_PREF_BOOL,
    MW_PREF_INT,
    MW_PREF_STRING
} MWPrefType;

typedef struct
{
    const char* name;
    MWPrefType type;
    int found;
    void* value;
}  MWProjectPrefs;

#define END_MW_PROJECT_PREFS { NULL, MW_PREF_NONE, FALSE, NULL }


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
  #define mw_fraction_done(x) boinc_fraction_done(x)
  #define mw_time_to_checkpoint() boinc_time_to_checkpoint()
  #define mw_checkpoint_completed() boinc_checkpoint_completed()
  #define mw_is_standalone() boinc_is_standalone()
#else
  #define mw_boinc_print(f, msg, ...)
  #define mw_finish(x) exit(x)
  #define mw_fopen(x,y) fopen((x),(y))
  #define mw_remove(x) remove((x))
  #define mw_begin_critical_section()
  #define mw_end_critical_section()
  #define mw_fraction_done(x) ((void) (x))
  #define mw_time_to_checkpoint() (FALSE)
  #define mw_checkpoint_completed()
  #define mw_is_standalone() (TRUE)
#endif /* BOINC_APPLICATION */


int mwBoincInit(MWInitType type);

/* BOINC file function wrappers */
char* mwReadFileResolved(const char* filename);
FILE* mwOpenResolved(const char* filename, const char* mode);

int mw_resolve_filename(const char* filename, char* buf, size_t bufSize);
int mw_file_exists(const char* file);


/* Things for checkinb BOINC status */
int mwGetAppInitData(void);
int mwReadProjectPrefs(MWProjectPrefs* prefs, const char* prefConfig);
const char* mwGuessPreferredPlatform(const char* progName);


int mwIsFirstRun(void);
const char* mwGetProjectPrefs(void);

int mwGetBoincNumCPU(void);
int mwGetBoincOpenCLDeviceIndex(void);
const char* mwGetBoincOpenCLPlatformVendor(void);


#if BOINC_APPLICATION
void mw_boinc_sleep(double t);
void* mw_graphics_make_shmem(const char* x, int y);
void* mw_graphics_get_shmem(const char* x);
#endif /* BOINC_APPLICATION */

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_BOINC_UTIL_H_ */

