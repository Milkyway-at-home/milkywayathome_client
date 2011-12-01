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
#endif

#ifndef _WIN32
  #include <sys/time.h>
#endif

#if HAVE_SYS_STAT_H
  #include <sys/stat.h>
#endif

#include <errno.h>
#include <string.h>

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

#else /* !BOINC_APPLICATION */

int mwGetAppInitData(void)
{
    return 0;
}

int mwIsFirstRun(void)
{
    return FALSE;
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

#endif /* BOINC_APPLICATION */



#define MAX_TAG_NAME_LENGTH 256

static char* matchTagName(const char* prefs, const char* name)
{
    int rcOpen, rcClose;
    char openTag[MAX_TAG_NAME_LENGTH];
    char closeTag[MAX_TAG_NAME_LENGTH];

    const char* openPos = NULL;
    const char* closePos = NULL;
    const char* beginItem = NULL;
    size_t itemSize = 0;
    char* item = NULL;

    if (!prefs || !name)
        return NULL;

    rcOpen = snprintf(openTag, sizeof(openTag), "<%s>", name);
    rcClose = snprintf(closeTag, sizeof(closeTag), "</%s>", name);
    if ((size_t) rcOpen == sizeof(openTag) || (size_t) rcClose == sizeof(closeTag))
    {
        mw_panic("Tag buffers too small\n");
    }

    openPos = strstr(prefs, openTag);
    if (!openPos)
    {
        mw_printf("Didn't find opening tag for preference '%s'\n", name);
        return NULL;
    }

    beginItem = openPos + rcOpen;
    closePos = strstr(beginItem, closeTag);
    if (!closePos)
    {
        mw_printf("Didn't find close tag for preference '%s'\n", name);
        return NULL;
    }

    itemSize = closePos - beginItem;
    item = (char*) mwMalloc(itemSize + 1);
    item[itemSize] = '\0';

    /* no strndup() on windows */
    strncpy(item, beginItem, itemSize);
    return item;
}

static int mwReadPref(MWProjectPrefs* pref, const char* prefConfig)
{
    char* item = NULL;
    char* endP = NULL;
    double d = 0.0;
    int i = 0;
    int rc = 0;

    item = matchTagName(prefConfig, pref->name);
    if (!item)
    {
        return 1;
    }

    errno = 0;
    switch (pref->type)
    {
        case MW_PREF_DOUBLE:
            d = strtod(item, &endP);
            if (item == endP || errno != 0)
            {
                rc = 1;
                mwPerror("Error parsing preference double from item '%s'", item);
            }
            else
            {
                pref->found = TRUE;
                *(double*) pref->value = d;
            }
            break;

        case MW_PREF_BOOL:
        case MW_PREF_INT:
            i = (int) strtol(item, &endP, 10);
            if (item == endP || errno != 0)
            {
                rc = 1;
                mwPerror("Error parsing preference int or bool from item '%s'", item);
            }
            else
            {
                pref->found = TRUE;
                *(int*) pref->value = i;
            }
            break;

        case MW_PREF_STRING:
        case MW_PREF_NONE:
        default:
            mw_panic("Implement me\n");
    }

    free(item);

    return rc;
}

int mwReadProjectPrefs(MWProjectPrefs* prefs, const char* prefConfig)
{
    int rc = 0;
    MWProjectPrefs* p = prefs;

    if (!p || !prefConfig)
        return 1;

    while (p->name)
    {
        rc |= mwReadPref(p, prefConfig);
        ++p;
    }

    return rc;
}

