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

#include "milkyway_config.h"

#if BOINC_APPLICATION
  #include <graphics2.h>
  #include <util.h>
  #include <diagnostics.h>
  #include <parse.h>
#endif

#if HAVE_SYS_TIME_H
  #include <sys/time.h>
#endif

#if HAVE_SYS_STAT_H
  #include <sys/stat.h>
#endif

#include <errno.h>
#include <limits.h>

#include "milkyway_util.h"
#include "milkyway_alloc.h"
#include "milkyway_boinc_util.h"


/* Work around areas broken in the BOINC libraries which make you use
 * C++ */

#if BOINC_APPLICATION

static APP_INIT_DATA mwAppInitData;
static int mwAppInitDataReady = FALSE;

int mwGetAppInitData(void)
{
    return (mwAppInitDataReady = boinc_get_init_data(mwAppInitData));
}

int mwIsFirstRun(void)
{
    if (!mwAppInitDataReady) /* If this fails we can't tell */
        return TRUE;

    return mwAppInitData.wu_cpu_time < 0.1;
}

const char* mwGetProjectPrefs(void)
{
    return mwAppInitData.project_preferences;
}

/* OpenCL stuff needs 6.13.7+ */
static bool mwBoincHasOpenCLData(void)
{
    const APP_INIT_DATA& aid = mwAppInitData;

    if (aid.major_version > 6)
        return true;

    if (aid.major_version == 6)
        return (aid.minor_version > 13 || (aid.minor_version == 13 && aid.release >= 1));

    return false;
}

int mwGetBoincOpenCLDeviceIndex(void)
{
    if (!mwAppInitDataReady || !mwBoincHasOpenCLData())
        return INT_MIN;

    return mwAppInitData.gpu_opencl_dev_index;
}

const char* mwGetBoincOpenCLPlatformVendor(void)
{
    if (!mwAppInitDataReady || !mwBoincHasOpenCLData())
        return NULL;

    return mwAppInitData.gpu_type;
}

static const int debugOptions = BOINC_DIAG_DUMPCALLSTACKENABLED
                              | BOINC_DIAG_HEAPCHECKENABLED
                              | BOINC_DIAG_MEMORYLEAKCHECKENABLED;

/* I don't understand why the graphics have a separate debug with
 * diagnostics API type stuff. */
static int mwBoincInitGraphics(int useDebug)
{
    return boinc_init_graphics_diagnostics(useDebug ? debugOptions : BOINC_DIAG_DEFAULTS);
}

static int mwBoincInitNormal(MWInitType type)
{
    int rc = 0;
    BOINC_OPTIONS options;

    if (type & MW_DEBUG)
    {
        rc |= boinc_init_diagnostics(debugOptions);
    }

    boinc_options_defaults(options);


    options.multi_thread = (type & MW_MULTITHREAD) > 0;
    options.normal_thread_priority = ((type & MW_CAL) || (type & MW_OPENCL)) > 0;

    rc |= boinc_init_options(&options);

    return rc;
}

int mwBoincInit(MWInitType type)
{
    int rc = 0;

    if (type & MW_GRAPHICS)
    {
        rc = mwBoincInitGraphics(type);
    }
    else
    {
        rc = mwBoincInitNormal(type);
    }

    if (rc)
    {
        mw_printf("Failed to init BOINC\n");
    }
    else
    {
        mwGetAppInitData();
    }

    return rc;
}

FILE* mwOpenResolved(const char* filename, const char* mode)
{
    int ret;
    char resolvedPath[1024];

    ret = boinc_resolve_filename(filename, resolvedPath, sizeof(resolvedPath));
    if (ret)
    {
        mw_printf("Error resolving file '%s': %d\n", filename, ret);
        return NULL;
    }

    return mw_fopen(resolvedPath, mode);
}

char* mwReadFileResolved(const char* filename)
{
    return mwFreadFile(mwOpenResolved(filename, "rb"), filename);
}

int mw_resolve_filename(const char* filename, char* buf, size_t bufSize)
{
    return boinc_resolve_filename(filename, buf, bufSize);
}

int mw_file_exists(const char* file)
{
    return boinc_file_exists(file);
}

/* The BOINC functions have C++ linkage for no reason */

void mw_boinc_sleep(double t)
{
    boinc_sleep(t);
}

void* mw_graphics_make_shmem(const char* x, int y)
{
    return boinc_graphics_make_shmem(x, y);
}

void* mw_graphics_get_shmem(const char* x)
{
    return boinc_graphics_get_shmem(x);
}

int mwReadProjectPrefs(MWProjectPrefs* prefs, const char* prefConfig)
{
    MIOFILE xmlPrefsFile;
    XML_PARSER parser(&xmlPrefsFile);

    if (!prefs || !prefConfig)
        return 1;

    xmlPrefsFile.init_buf_read(prefConfig);

    if (!parser.parse_start("project_preferences"))
    {
        mw_printf("Missing project_preferences start tag\n");
        return 1;
    }

    while (!parser.get_tag())
    {
        bool parseSuccess = false;
        MWProjectPrefs* p = prefs;

        if (!parser.is_tag)
        {
            mw_printf("Non-tag found in preferences: '%s'\n", parser.parsed_tag);
            continue;
        }

        if (parser.match_tag("/project_preferences"))
        {
            return 0;
        }

        while (p->name)
        {
            if (parser.match_tag(p->name))
            {
                bool tmp;

                switch (p->type)
                {
                    case MW_PREF_DOUBLE:
                        parseSuccess = parser.parse_double(p->name, *(double*) p->value);
                        break;

                    case MW_PREF_BOOL:
                        parseSuccess = parser.parse_bool(p->name, tmp);
                        *(int*) p->value = (int) tmp;
                        break;

                    case MW_PREF_INT:
                        parseSuccess = parser.parse_int(p->name, *(int*) p->value);
                        break;

                    case MW_PREF_STRING:
                    case MW_PREF_NONE:
                    default:
                        mw_panic("Unhandled preference type %d\n", p->type);
                }

                if (!parseSuccess)
                {
                    mw_printf("Error reading preference '%s'\n", p->name);
                }
                else
                {
                    p->found = TRUE;
                }

                break;
            }

            ++p;
        }

        if (!parseSuccess)
        {
            parser.skip_unexpected(true, "project preferences");
        }
    }

    mw_printf("Reading preferences ended prematurely\n");
    return 1;
}

#else /* !BOINC_APPLICATION */

int mwGetAppInitData(void)
{
    return 0;
}

int mwIsFirstRun(void)
{
    return TRUE;
}

const char* mwGetProjectPrefs(void)
{
    return NULL;
}

int mwBoincInit(MWInitType type)
{
    (void) type;
    return 0;
}

FILE* mwOpenResolved(const char* filename, const char* mode)
{
    return mw_fopen(filename, mode);
}

char* mwReadFileResolved(const char* filename)
{
    return mwReadFile(filename);
}

int mw_resolve_filename(const char* filename, char* buf, size_t bufSize)
{
    int rc;

    assert(buf != filename);
    rc = snprintf(buf, bufSize, "%s", filename);
    return (rc == -1) || ((size_t) rc == bufSize);
}

int mw_file_exists(const char* file)
{
    struct stat statBuf;
    return !stat(file, &statBuf);
}

int mwReadProjectPrefs(MWProjectPrefs* prefs, const char* prefConfig)
{
    return 0;
}

int mwGetBoincOpenCLDeviceIndex(void)
{
    return INT_MIN;
}

const char* mwGetBoincOpenCLPlatformVendor(void)
{
    return NULL;
}

#endif /* BOINC_APPLICATION */

