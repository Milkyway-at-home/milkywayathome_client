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


#include "milkyway_config.h"
#include "milkyway_extra.h"
#include "milkyway_cpp_util.h"
#include <string.h>

#if BOINC_APPLICATION
  #include <graphics2.h>
#endif

/* Work around areas broken in the BOINC libraries which make you use
 * C++ */

#if BOINC_APPLICATION

int mwGetMWAppInitData(MWAppInitData* mwaid)
{
    struct APP_INIT_DATA aid;

    /* The C API function, yet the headers are broken so you need C++ anyway */
    if (boinc_get_init_data_p(&aid))
        return 1;

    mwaid->majorVersion = aid.major_version;
    mwaid->minorVersion = aid.minor_version;
    mwaid->release = aid.release;
    mwaid->appVersion = aid.app_version;
    mwaid->checkpointPeriod = aid.checkpoint_period;

    strncpy(mwaid->wuName, aid.wu_name, sizeof(mwaid->wuName));
    strncpy(mwaid->projectDir, aid.project_dir, sizeof(mwaid->projectDir));
    strncpy(mwaid->boincDir, aid.boinc_dir, sizeof(mwaid->boincDir));
    strncpy(mwaid->projectPrefs, aid.project_preferences, sizeof(mwaid->projectPrefs));

    return 0;
}

void mwGetBoincOptionsDefault(BOINC_OPTIONS* options)
{
    boinc_options_defaults(*options);
}

/* The BOINC functions have C++ linkage for no reason */
void* mw_graphics_make_shmem(const char* x, int y)
{
    return boinc_graphics_make_shmem(x, y);
}

void* mw_graphics_get_shmem(const char* x)
{
    return boinc_graphics_get_shmem(x);
}

#endif /* BOINC_APPLICATION */

