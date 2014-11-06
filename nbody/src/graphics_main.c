/*
 * Copyright (c) 2011-2012 Matthew Arsenault
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef _WIN32
  #pragma comment(linker, "/SUBSYSTEM:WINDOWS")
#endif

#include "nbody_config.h"
#include "nbody_gl.h"
#include "nbody_graphics.h"
#include "milkyway_util.h"

#include <signal.h>

#if USE_POSIX_SHMEM
  #include <sys/mman.h>
  #include <sys/stat.h>
  #include <fcntl.h>
  #include <errno.h>
#endif /* USE_POSIX_SHMEM */

typedef struct NBodyGLPrefs
{
    int eventPollPeriod;
    int blockSimulation;
    int updatePeriod;
    int floatView;
    double floatSpeed;
    double texturedPointSize;
    double pointPointSize;
    int untexturedPoints;
    int monochromatic;
    int originCentered;
    int showInfo;
    int showAxes;
    int showOrbitTrace;
} NBodyGLPrefs;

static void nbglPrintCopyright(void)
{
    mw_printf(
        "Milkyway@Home N-body graphics client %d.%d\n\n"
        "Copyright (c) 2012 Matthew Arsenault\n"
        "\n"
        "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        "This is free software: you are free to change and redistribute it.\n"
        "There is NO WARRANTY, to the extent permitted by law.\n"
        "\n"
        " Incorporates works covered by the following copyright and\n"
        " permission notices:\n"
        "\n"
        "--------------------------------------------------------------------------------\n"
        "Milkyway image courtesy NASA/JPL-Caltech.\n"
        "--------------------------------------------------------------------------------\n"
        "Font: Roboto Regular Copyright (c) 2008 The Android Open Source Project (Apache 2.0)\n"
        "--------------------------------------------------------------------------------\n"
        "This software contains source code provided by NVIDIA Corporation.\n"
        "--------------------------------------------------------------------------------\n"
        "GLM: Copyright (c) 2005 - 2012 G-Truc Creation\n"
        "glutil: Copyright (c) 2011 by Jason L. McKesson\n"
        "\n"
        "The MIT License\n"
        "\n"
        " Permission is hereby granted, free of charge, to any person obtaining a copy\n"
        " of this software and associated documentation files (the \"Software\"), to deal\n"
        " in the Software without restriction, including without limitation the rights\n"
        " to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
        " copies of the Software, and to permit persons to whom the Software is\n"
        " furnished to do so, subject to the following conditions:\n"
        "\n"
        " The above copyright notice and this permission notice shall be included in\n"
        " all copies or substantial portions of the Software.\n"
        "\n"
        " THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
        " IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
        " FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE\n"
        " AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
        " LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
        " OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
        " THE SOFTWARE.\n"
        "\n"
        "--------------------------------------------------------------------------------\n"
        "GLFW - Free, open source, portable framework for OpenGL application development.\n"
        "\n"
        "Copyright (c) 2002-2006 Marcus Geelnard\n"
        "Copyright (c) 2006-2010 Camilla Berglund <elmindreda@elmindreda.org>\n"
        "\n"
        "This software is provided 'as-is', without any express or implied\n"
        "warranty. In no event will the authors be held liable for any damages\n"
        "arising from the use of this software.\n"
        "\n"
        "Permission is granted to anyone to use this software for any purpose,\n"
        "including commercial applications, and to alter it and redistribute it\n"
        "freely, subject to the following restrictions:\n"
        "\n"
        "1. The origin of this software must not be misrepresented; you must not\n"
        "   claim that you wrote the original software. If you use this software\n"
        "   in a product, an acknowledgment in the product documentation would\n"
        "   be appreciated but is not required.\n"
        "\n"
        "2. Altered source versions must be plainly marked as such, and must not\n"
        "   be misrepresented as being the original software.\n"
        "\n"
        "3. This notice may not be removed or altered from any source\n"
        "   distribution.\n"
        "\n"
        "--------------------------------------------------------------------------------\n"
        "Freetype GL - A C OpenGL Freetype engine\n"
        "\n"
        "Copyright 2011, 2012 Nicolas P. Rougier. All rights reserved.\n"
        "Redistribution and use in source and binary forms, with or without\n"
        "modification, are permitted provided that the following conditions are met:\n"
        "\n"
        "  1. Redistributions of source code must retain the above copyright notice,\n"
        "    this list of conditions and the following disclaimer.\n"
        "\n"
        "  2. Redistributions in binary form must reproduce the above copyright\n"
        "    notice, this list of conditions and the following disclaimer in the\n"
        "    documentation and/or other materials provided with the distribution.\n"
        "\n"
        "THIS SOFTWARE IS PROVIDED BY NICOLAS P. ROUGIER ''AS IS'' AND ANY EXPRESS OR\n"
        "IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF\n"
        "MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO\n"
        "EVENT SHALL NICOLAS P. ROUGIER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,\n"
        "INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES\n"
        "(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\n"
        "LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND\n"
        "ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n"
        "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF\n"
        "THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n"
        " * The views and conclusions contained in the software and documentation are\n"
        "those of the authors and should not be interpreted as representing official\n"
        "policies, either expressed or implied, of Nicolas P. Rougier.\n"
        "--------------------------------------------------------------------------------\n"
        "\n",
        NBODY_VERSION_MAJOR,
        NBODY_VERSION_MINOR
        );
}

static const VisArgs defaultVisArgs =
{
    /* .fullscreen        */ FALSE,
    /* .plainFullscreen   */ FALSE,
    /* .width             */ 0,
    /* .height            */ 0,
    /* .eventPollPeriod   */ DEFAULT_EVENT_POLL_PERIOD,
    /* .printFrames       */ DEFAULT_PRINTFRAMES,

    /* .quitOnComplete    */ DEFAULT_QUIT_ON_COMPLETE,
    /* .blockSimulation   */ DEFAULT_BLOCK_SIMULATION,
    /* .updatePeriod      */ DEFAULT_UPDATE_PERIOD,
    /* .noFloat           */ FALSE,
    /* .floatSpeed        */ DEFAULT_FLOAT_SPEED,
    /* .texturedPointSize */ DEFAULT_TEXTURED_POINT_SIZE,
    /* .pointPointSize    */ DEFAULT_POINT_POINT_SIZE,
    /* .untexturedPoints  */ DEFAULT_UNTEXTURED_POINTS,
    /* .monochrome        */ DEFAULT_MONOCHROMATIC,
    /* .originCentered    */ DEFAULT_ORIGIN_CENTERED,
    /* .noDrawInfo        */ FALSE,
    /* .drawAxes          */ DEFAULT_SHOW_AXES,
    /* .drawOrbitTrace    */ DEFAULT_SHOW_ORBIT_TRACE,

    /* .pid               */ 0,
    /* .file              */ NULL,
    /* .instanceId        */ -1
};

static void freeVisArgs(VisArgs* args)
{
    free(args->file);
    args->file = NULL;
}

static const NBodyGLPrefs nbglDefaultPrefs =
{
    /* .eventPollPeriod   */ DEFAULT_EVENT_POLL_PERIOD,
    /* .blockSimulation   */ DEFAULT_BLOCK_SIMULATION,
    /* .updatePeriod      */ DEFAULT_UPDATE_PERIOD,
    /* .floatView         */ DEFAULT_FLOAT,
    /* .floatSpeed        */ DEFAULT_FLOAT_SPEED,
    /* .texturedPointSize */ DEFAULT_TEXTURED_POINT_SIZE,
    /* .pointPointSize    */ DEFAULT_POINT_POINT_SIZE,
    /* .untexturedPoints  */ DEFAULT_UNTEXTURED_POINTS,
    /* .monochromatic     */ DEFAULT_MONOCHROMATIC,
    /* .originCentered    */ DEFAULT_ORIGIN_CENTERED,
    /* .showInfo          */ DEFAULT_SHOW_INFO,
    /* .drawAxes          */ DEFAULT_SHOW_AXES,
    /* .drawOrbitTrace    */ DEFAULT_SHOW_ORBIT_TRACE
};

static void nbglReadPreferences(VisArgs* args)
{
    static NBodyGLPrefs prefs;

    static MWProjectPrefs nbglPrefs[] =
        {
            { "event_poll_period",   MW_PREF_INT,    FALSE, &prefs.eventPollPeriod   },
            { "block_simulation",    MW_PREF_BOOL,   FALSE, &prefs.blockSimulation   },
            { "update_period",       MW_PREF_INT,    FALSE, &prefs.updatePeriod      },
            { "float",               MW_PREF_BOOL,   FALSE, &prefs.floatView         },
            { "float_speed",         MW_PREF_DOUBLE, FALSE, &prefs.floatSpeed        },
            { "textured_point_size", MW_PREF_DOUBLE, FALSE, &prefs.texturedPointSize },
            { "point_point_size",    MW_PREF_DOUBLE, FALSE, &prefs.pointPointSize    },
            { "untextured_points",   MW_PREF_BOOL,   FALSE, &prefs.untexturedPoints  },
            { "monochromatic",       MW_PREF_BOOL,   FALSE, &prefs.monochromatic     },
            { "origin_centered",     MW_PREF_BOOL,   FALSE, &prefs.originCentered    },
            { "show_info",           MW_PREF_BOOL,   FALSE, &prefs.showInfo          },
            { "show_axes",           MW_PREF_BOOL,   FALSE, &prefs.showAxes          },
            { "show_orbit_trace",    MW_PREF_BOOL,   FALSE, &prefs.showOrbitTrace    },
            END_MW_PROJECT_PREFS
        };

    prefs = nbglDefaultPrefs;

    if (mwGetAppInitData())
    {
        mw_printf("Error reading app init data. Project preferences will not be used\n");
    }
    else
    {
        mwReadProjectPrefs(nbglPrefs, mwGetProjectPrefs());
    }

    /* Any successfully found setting will be used; otherwise it will get the default */
    args->eventPollPeriod   = prefs.eventPollPeriod;
    args->blockSimulation   = prefs.blockSimulation;
    args->updatePeriod      = prefs.updatePeriod;
    args->noFloat           = !prefs.floatView;
    args->floatSpeed        = (float) prefs.floatSpeed;
    args->texturedPointSize = (float) prefs.texturedPointSize;
    args->pointPointSize    = (float) prefs.pointPointSize;
    args->untexturedPoints  = prefs.untexturedPoints;
    args->monochromatic     = prefs.monochromatic;
    args->originCentered    = prefs.originCentered;
    args->noDrawInfo        = !prefs.showInfo;
    args->drawAxes          = prefs.showAxes;
    args->drawOrbitTrace    = prefs.showOrbitTrace;
}

static int nbglHandleVisArguments(int argc, const char** argv, VisArgs* visOut)
{
    poptContext context;
    int failed = FALSE;
    static int copyright = FALSE;
    static VisArgs visArgs;

    static const struct poptOption options[] =
    {
        {
            "fullscreen", 'f',
            POPT_ARG_NONE, &visArgs.fullscreen,
            0, "Start in screensaver mode", NULL
        },

        /* Since BOINC passes --fullscreen to start as a screensaver,
         * and we might want to start as fullscreen but without any
         * keys/moving quitting */
        {
            "plain-fullscreen", '\0',
            POPT_ARG_NONE, &visArgs.plainFullscreen,
            0, "Start as fullscreen, but without quitting on any motion", NULL
        },

        {
            "instance-id", 'i',
            POPT_ARG_INT, &visArgs.instanceId,
            0, "Instance id of main process to attach", NULL
        },

        {
            "static-input", 's',
            POPT_ARG_STRING, &visArgs.file,
            0, "Load from an output file (in Cartesian coordinates) and statically display", NULL
        },

        {
            "width", 'w',
            POPT_ARG_INT, &visArgs.width,
            0, "Starting width of window", NULL
        },

        {
            "height", 'h',
            POPT_ARG_INT, &visArgs.height,
            0, "Starting height of window", NULL
        },

        {
            "event-poll-period", 'e',
            POPT_ARG_INT, &visArgs.eventPollPeriod,
            0, "Period to poll for events / draw in milliseconds", NULL
        },

        {
            "block-simulation", 'b',
            POPT_ARG_NONE, &visArgs.blockSimulation,
            0, "Make simulation wait for graphics so every frame is drawn", NULL
        },

        {
            "update-period", 'p',
            POPT_ARG_INT, &visArgs.updatePeriod,
            0, "Interval between scene refreshes", NULL
        },

        {
            "quit-on-complete", 'q',
            POPT_ARG_NONE, &visArgs.quitOnComplete,
            0, "Graphics should quit when the simulation completes", NULL
        },

        {
            "no-float", 'r',
            POPT_ARG_NONE, &visArgs.noFloat,
            0, "By default do not float view around randomly", NULL
        },

        {
            "untextured-points", 'T',
            POPT_ARG_NONE, &visArgs.untexturedPoints,
            0, "Use faster but possibly uglier drawing of points", NULL
        },

        {
            "monochromatic", 'm',
            POPT_ARG_NONE, &visArgs.monochromatic,
            0, "All particles have same color", NULL
        },

        {
            "origin-centered", 'o',
            POPT_ARG_NONE, &visArgs.originCentered,
            0, "Focus on the galactic center instead of system's center of mass", NULL
        },

        {
            "no-show-info", 'I',
            POPT_ARG_NONE, &visArgs.noDrawInfo,
            0, "Disable displaying basic information about simulation", NULL
        },

        {
            "show-axes", 'a',
            POPT_ARG_NONE, &visArgs.drawAxes,
            0, "Display simple axes for reference", NULL
        },

        {
            "show-orbit-trace", 't',
            POPT_ARG_NONE, &visArgs.drawOrbitTrace,
            0, "Show path of center of mass", NULL
        },

        {
            "print-frames", 'y',
            POPT_ARG_STRING, &visArgs.printFrames,
            0, "Print frame images for making videos", NULL
        },


        {
            "copyright", '\0',
            POPT_ARG_NONE, &copyright,
            0, "Print copyright information and exit", NULL
        },

        POPT_AUTOHELP
        POPT_TABLEEND
    };

    /* TODO: Check project prefs */

    visArgs = defaultVisArgs;

    if (BOINC_APPLICATION)
    {
        nbglReadPreferences(&visArgs);
    }

    context = poptGetContext(argv[0], argc, argv, options, 0);

    if (mwReadArguments(context) < 0)
    {
        poptPrintHelp(context, stderr, 0);
        failed = TRUE;
    }
    poptFreeContext(context);

    if (visArgs.instanceId < 0) /* Default to first */
    {
        visArgs.instanceId = 0;
    }

    if (copyright)
    {
        nbglPrintCopyright();
        exit(0);
    }

    *visOut = visArgs;
    return failed;
}

static int nbglBoincGraphicsInit(int debug)
{
    if (BOINC_APPLICATION)
    {
        MWInitType type = MW_GRAPHICS | (debug ? MW_DEBUG : 0);
        if (mwBoincInit(type))
        {
            mw_printf("BOINC graphics init failed\n");
            return 1;
        }
    }

    return 0;
}


static int readCoordinateSystem = FALSE;
static int readHasGalaxy = FALSE;
static int readCenterOfMass = FALSE;
static int fileUsesCartesian = TRUE;

/* Return TRUE if succesfully matched a line */
static int nbglTryReadSceneItems(const char* lineBuf, scene_t* scene)
{

    if (!readCoordinateSystem)
    {
        if (sscanf(lineBuf,
                   " cartesian = %d \n",
                   &fileUsesCartesian) == 1)
        {
            readCoordinateSystem = TRUE;

            if (!fileUsesCartesian)
            {
                mw_printf("Warning: Noncartesian output files not implemented\n");
            }

            return TRUE;
        }
    }

    if (!readHasGalaxy)
    {
        if (sscanf(lineBuf,
                   " hasMilkyway = %d \n",
                   &scene->hasGalaxy) == 1)
        {
            readHasGalaxy = TRUE;
            return TRUE;
        }
    }

    if (!readCenterOfMass)
    {
        float* cmPos = scene->queue.info[0].rootCenterOfMass;

        if (sscanf(lineBuf,
                   " centerOfMass = %f , %f , %f \n",
                   &cmPos[0], &cmPos[1], &cmPos[2]) == 3)
        {
            readCenterOfMass = TRUE;
            return TRUE;
        }
        else
        {
            cmPos[0] = cmPos[1] = cmPos[2] = 0.0f;
        }
    }

    return FALSE;
}

static scene_t* nbglLoadStaticSceneFromFile(const char* filename)
{
    FILE* f;
    size_t lnCount; /* ~= nbody */
    size_t line = 0;
    int nbody = 0;
    char lnBuf[4096];
    int rc = 0;
    int ignore;
    double x, y, z;
    double vx, vy, vz;
    double lambda;
    FloatPos* r;
    int hasError = FALSE;
    scene_t* scene = NULL;

    f = fopen(filename, "r");
    if (!f)
    {
        mwPerror("Failed to open file '%s'", filename);
        return NULL;
    }

    lnCount = mwCountLinesInFile(f);
    if (lnCount == 0)
    {
        mw_printf("Error counting lines from file '%s'\n", filename);
        fclose(f);
        return NULL;
    }

    scene = mwCalloc(sizeof(scene_t) + 1 * (lnCount * sizeof(FloatPos)), sizeof(char));

    /* We don't neeqd the buffering features so just use first buffer slot */
    r = nbSceneGetQueueBuffer(scene, 0);
    scene->hasInfo = FALSE;
    scene->staticScene = TRUE;

    /* Make read data fake that we have 1 element in the queue */
    OPA_store_int(&scene->queue.head, 0);
    OPA_store_int(&scene->queue.tail, 1);


    readCoordinateSystem = FALSE;
    readHasGalaxy = FALSE;
    readCenterOfMass = FALSE;

    while (fgets(lnBuf, (int) sizeof(lnBuf), f) && line < lnCount)
    {
        ++line;

        if (strlen(lnBuf) + 1 >= sizeof(lnBuf))
        {
            mw_printf("Error reading histogram line "ZU" (Line buffer too small): %s", line, lnBuf);
            hasError = TRUE;
            break;
        }

        /* Skip comments and blank lines */
        if (lnBuf[0] == '#' || lnBuf[0] == '\n')
            continue;

        if (nbglTryReadSceneItems(lnBuf, scene))
            continue;

        rc = sscanf(lnBuf,
                    "%d , %lf , %lf , %lf , %lf , %lf , %lf ",
                    &ignore,
                    &x, &y, &z,
                    &vx, &vy, &vz);
        if (rc == 7)
        {
            /* May or may not be there */
            rc = sscanf(lnBuf, " , %lf \n", &lambda);
            if (rc != 1)
            {
                sscanf(lnBuf, " \n");
            }

            r[nbody].x = (float) x;
            r[nbody].y = (float) y;
            r[nbody].z = (float) z;

            ++nbody;
        }
        else
        {
            hasError = TRUE;
        }
    }

    if (fclose(f))
    {
        mwPerror("Failed to close file '%s'", filename);
    }

    if (!readCoordinateSystem || !readHasGalaxy || !readCenterOfMass)
    {
        mw_printf("Warning: Failed to read some scene info items\n");
    }

    scene->nbody = nbody;
    if (hasError)
    {
        free(scene);
        return NULL;
    }

    return scene;
}

static const char* nbGetShmemName(int instanceId)
{
    static char name[256];

    if (snprintf(name, sizeof(name), NBODY_SHMEM_NAME_FMT_STR, instanceId) == sizeof(name))
    {
        mw_panic("Name buffer too small for shared memory name\n");
    }

    return name;
}

#if USE_POSIX_SHMEM

static scene_t* nbglConnectSharedScene(int instanceId)
{
    int shmId;
    const int mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
    struct stat sb;
    const char* name;
    size_t calcSize;
    scene_t* scene = NULL;

    name = nbGetShmemName(instanceId);
    shmId = shm_open(name, O_RDWR, mode);
    if (shmId < 0)
    {
        mwPerror("Error getting shared memory");
        return NULL;
    }

    if (fstat(shmId, &sb) < 0)
    {
        mwPerror("shmem fstat");
        shm_unlink(name);
        return NULL;
    }

    scene = (scene_t*) mmap(NULL, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, shmId, 0);
    if (scene == MAP_FAILED)
    {
        mwPerror("mmap: Failed to mmap shared memory");
        if (shm_unlink(name) < 0)
        {
            mwPerror("Unlink shared memory");
        }

        return NULL;
    }

    if (   sb.st_size < (ssize_t) sizeof(scene_t)
        || sb.st_size < (ssize_t) scene->sceneSize
        || sb.st_size < (ssize_t) (calcSize = nbFindShmemSize(scene->nbody, scene->nSteps))
        || calcSize != scene->sceneSize)
    {
        mw_printf("Shared memory segment is impossibly small ("ZU")\n", (size_t) sb.st_size);
        if (shm_unlink(name) < 0)
        {
            mwPerror("Unlink shared memory");
        }

        munmap(scene, (size_t) sb.st_size);
        return NULL;
    }

    return scene;
}

static void nbglUnmapScene(scene_t* scene)
{
    munmap(scene, scene->sceneSize);
}

#elif USE_WIN32_SHARED_MAP

static scene_t* nbglConnectSharedScene(int instanceId)
{
    HANDLE mapFile;
    scene_t* scene;
    const char* name;
    size_t size;

    name = nbGetShmemName(instanceId);

    mapFile = OpenFileMapping(
        FILE_MAP_ALL_ACCESS,   /* read/write access */
        FALSE,                 /* do not inherit the name */
        name);                 /* name of mapping object */
    if (!mapFile)
    {
        mwPerrorW32("Could not open shared file mapping object '%s'", name);
        return NULL;
    }

    scene = (scene_t*) MapViewOfFile(mapFile,              /* handle to map object */
                                     FILE_MAP_ALL_ACCESS,  /* read/write permission */
                                     0,
                                     0,
                                     0); /* map entire file mapping */
    CloseHandle(mapFile);
    if (!scene)
    {
        mwPerrorW32("Could not map view of file mapping object '%s'", name);
        UnmapViewOfFile((LPCVOID) scene);
        return NULL;
    }

    /* Because this API sucks and doesn't give us a way to find the size of a
       paging file backed shared mapped file use the size we stored ourselves. */
    size = scene->sceneSize;
    if (size < sizeof(scene_t) || size != nbFindShmemSize(scene->nbody))
    {
        mw_printf("Shared memory segment '%s' is impossibly small (%u)\n", name, size);
        CloseHandle(mapFile);
        UnmapViewOfFile((LPCVOID) scene);
        return NULL;
    }

    return scene;
}

static void nbglUnmapScene(scene_t* scene)
{
    UnmapViewOfFile((LPCVOID) scene);
}

#elif USE_BOINC_SHMEM

static scene_t* nbglAttemptConnectSharedScene(void)
{
    scene_t* scene = (scene_t*) mw_graphics_get_shmem(NBODY_BIN_NAME);
    if (!scene)
    {
        mw_printf("Failed to connect to shared scene\n");
    }

    return scene;
}

#define MAX_TRIES 5
#define RETRY_INTERVAL 250

/* In case the main application isn't ready yet, try and wait for a while */
static scene_t* nbglConnectSharedScene(int instanceId)
{
    int tries = 0;
    scene_t* scene = NULL;

    (void) instanceId;

    while (tries < MAX_TRIES)
    {
        if ((scene = nbglAttemptConnectSharedScene()))
        {
            return scene; /* Error if something already attached */
        }

        mwMilliSleep(RETRY_INTERVAL);
        ++tries;
    }

    mw_printf("Could not attach to simulation after %d attempts\n", MAX_TRIES);
    return NULL;
}

static void nbglUnmapScene(scene_t* scene)
{

}

#else
  #error No shared memory method used
#endif /* USE_POSIX_SHMEM */

static int nbglCheckConnectedVersion(const scene_t* scene)
{
    if (   scene->nbodyMajorVersion != NBODY_VERSION_MAJOR
        || scene->nbodyMinorVersion != NBODY_VERSION_MINOR)
    {
        mw_printf("Graphics version (%d.%d) does not match application version (%d.%d)\n",
                  NBODY_VERSION_MAJOR,
                  NBODY_VERSION_MINOR,
                  scene->nbodyMajorVersion,
                  scene->nbodyMinorVersion);
        return 1;
    }

    return 0;
}

static void nbglReleaseSceneLocks(scene_t* scene)
{
    OPA_store_int(&scene->paused, 0);
    OPA_store_int(&scene->blockSimulationOnGraphics, 0);
    OPA_store_int(&scene->attachedLock, 0);

    /* Set this last */
    OPA_store_int(&scene->attachedPID, 0);
}

static int nbglGetExclusiveSceneAccess(scene_t* scene)
{
    int pid = (int) getpid();
    int oldPID = OPA_cas_int(&scene->attachedLock, 0, pid);
    if (oldPID != 0)
    {
        if (mwProcessIsAlive(oldPID))
        {
            mw_printf("Could not get exclusive access to simulation shared segment "
                      "(Owned by process %d)\n",
                      oldPID);
            return 1;
        }
        else
        {
            mw_printf("Simulation shared segment owned by dead process %d, stealing it\n",
                      oldPID);

            /* Process is dead, steal the lock */
            nbglReleaseSceneLocks(scene);
            return 0;
        }
    }
    else
    {
        OPA_store_int(&scene->attachedPID, 0);
        return 0;
    }
}

static scene_t* g_scene = NULL;

static void nbglCleanupAttached(void)
{
    if (g_scene)
    {
        nbglReleaseSceneLocks(g_scene);
        nbglUnmapScene(g_scene);
        g_scene = NULL;
    }
}

#if HAVE_SIGACTION

static void nbglSigInfoHandler(int sig, siginfo_t* siginfo, void* context)
{
    (void) siginfo, (void) context;

    nbglCleanupAttached();
    raise(sig);
}

static int nbglInstallExitHandlers(void)
{
    int rc = 0;
    struct sigaction action;

    memset(&action, 0, sizeof(action));

    action.sa_sigaction = nbglSigInfoHandler;
    sigemptyset(&action.sa_mask);
    action.sa_flags = SA_SIGINFO | SA_RESETHAND; /* Use sa_sigaction handler, reset handler on call */


    rc |= sigaction(SIGINT, &action, NULL);
    rc |= sigaction(SIGABRT, &action, NULL);
    rc |= sigaction(SIGFPE, &action, NULL);
    rc |= sigaction(SIGSEGV, &action, NULL);

    rc |= sigaction(SIGQUIT, &action, NULL);
    rc |= sigaction(SIGKILL, &action, NULL);
    rc |= sigaction(SIGUSR1, &action, NULL);

    rc |= atexit(nbglCleanupAttached);

    return rc;
}

#else

static void nbglSigHandler(int sig)
{
    nbglCleanupAttached();
    signal(sig, SIG_DFL);
    raise(sig);
}

static int nbglInstallExitHandlers(void)
{
    signal(SIGINT, nbglSigHandler);
    signal(SIGABRT, nbglSigHandler);
    signal(SIGFPE, nbglSigHandler);
    signal(SIGSEGV, nbglSigHandler);

  #ifndef _WIN32
    signal(SIGQUIT, nbglSigHandler);
    signal(SIGKILL, nbglSigHandler);
    signal(SIGUSR1, nbglSigHandler);
  #endif

    return atexit(nbglCleanupAttached);
}

#endif /* HAVE_SIGACTION */

#ifdef __APPLE__
#define main nbgl_main_apple
#endif

int main(int argc, const char* argv[])
{
    int rc;
    VisArgs flags;
    scene_t* scene = NULL;

    if (nbglBoincGraphicsInit(FALSE))
        return 1;

    if (nbglHandleVisArguments(argc, (const char**) argv, &flags))
    {
        freeVisArgs(&flags);
        return 1;
    }

    if (!flags.file)
    {
        scene = nbglConnectSharedScene(flags.instanceId);
        if (!scene)
        {
            freeVisArgs(&flags);
            return 1;
        }

        if (nbglCheckConnectedVersion(scene) || nbglGetExclusiveSceneAccess(scene))
        {
            freeVisArgs(&flags);
            return 1;
        }

        mw_report("Process %d acquired instance id %d\n", (int) getpid(), flags.instanceId);

        nbglInstallExitHandlers();
        g_scene = scene;
    }
    else
    {
        scene = nbglLoadStaticSceneFromFile(flags.file);
        if (!scene)
        {
            freeVisArgs(&flags);
            return 1;
        }
    }

    rc = nbglRunGraphics(scene, &flags);

    if (flags.file)
    {
        free(scene);
    }

    freeVisArgs(&flags);
    mw_finish(rc);

    return rc;
}

