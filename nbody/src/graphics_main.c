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

#include "nbody_config.h"
#include "nbody_gl.h"
#include "nbody_graphics.h"
#include "milkyway_util.h"

#include <signal.h>

#if USE_SHMEM
  #include <sys/mman.h>
  #include <sys/stat.h>
  #include <fcntl.h>
  #include <errno.h>
#endif /* USE_SHMEM */

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
    /* .fullscreen       */ FALSE,
    /* .plainFullscreen  */ FALSE,
    /* .width            */ 0,
    /* .height           */ 0,
    /* .blockSimulation  */ FALSE,
    /* .monochrome       */ FALSE,
    /* .untexturedPoints */ FALSE,
    /* .originCenter     */ FALSE,
    /* .drawAxes         */ FALSE,
    /* .drawOrbitTrace   */ FALSE,
    /* .noFloat          */ FALSE,
    /* .pid              */ 0,
    /* .file             */ NULL,
    /* .instanceId       */ -1
};

static void freeVisArgs(VisArgs* args)
{
    free(args->file);
    args->file = NULL;
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
            "block-simulation", 'b',
            POPT_ARG_NONE, &visArgs.blockSimulation,
            0, "Make simulation wait for graphics so every frame is drawn", NULL
        },

        {
            "monochromatic", 'm',
            POPT_ARG_NONE, &visArgs.monochrome,
            0, "All particles have same color", NULL
        },

        {
            "untextured-points", 'T',
            POPT_ARG_NONE, &visArgs.untexturedPoints,
            0, "Use faster but possibly uglier drawing of points", NULL
        },

        {
            "origin-center", 'o',
            POPT_ARG_NONE, &visArgs.originCenter,
            0, "Focus on the galactic center instead of system's center of mass", NULL
        },

        {
            "show-axes", 'a',
            POPT_ARG_NONE, &visArgs.drawAxes,
            0, "Draw simple axes for reference", NULL
        },

        {
            "show-orbit-trace", 't',
            POPT_ARG_NONE, &visArgs.drawOrbitTrace,
            0, "Show path of center of mass", NULL
        },

        {
            "no-float", 'r',
            POPT_ARG_NONE, &visArgs.noFloat,
            0, "By default do not float view around randomly", NULL
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
            "copyright", '\0',
            POPT_ARG_NONE, &copyright,
            0, "Print copyright information and exit", NULL
        },

        POPT_AUTOHELP
        POPT_TABLEEND
    };

    /* TODO: Check project prefs */

    visArgs = defaultVisArgs;
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

    /* We don't need the buffering features so just use first buffer slot */
    r = &scene->queue.bodyData[0];
    scene->hasInfo = FALSE;
    scene->staticScene = TRUE;


    /* Skip the 1st line with the # comment */
    fgets(lnBuf, sizeof(lnBuf), f);

    while (rc != EOF)
    {
        ++line;
        rc = fscanf(f,
                    "%d , %lf , %lf , %lf , %lf , %lf , %lf ",
                    &ignore,
                    &x, &y, &z,
                    &vx, &vy, &vz);
        if (rc == 7)
        {
            ++nbody;
            /* May or may not be there */
            rc = fscanf(f, " , %lf \n", &lambda);
            if (rc != 1)
            {
                fscanf(f, " \n");
            }
            assert(line < lnCount);

            r[line].x = (float) x;
            r[line].y = (float) y;
            r[line].z = (float) z;
        }
        else if (rc != EOF)
        {
            mw_printf("Error reading '%s' at line "ZU"\n", filename, line);
        }
    }

    if (rc != EOF)
    {
        fclose(f);
        free(scene);
        return NULL;
    }

    scene->nbody = nbody;

    if (fclose(f))
    {
        mwPerror("Failed to close file '%s'", filename);
    }

    return scene;
}


#if USE_SHMEM

/* FIXME: Duplicated in nbody_shmem.c */
static size_t nbFindShmemSize(int nbody)
{
    size_t snapshotSize = sizeof(NBodyCircularQueue) + nbody * sizeof(FloatPos);
    return sizeof(scene_t) + NBODY_CIRC_QUEUE_SIZE * snapshotSize;
}

static scene_t* nbglConnectSharedScene(int instanceId)
{
    int shmId;
    const int mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
    struct stat sb;
    char name[128];
    scene_t* scene = NULL;

    if (snprintf(name, sizeof(name), "/milkyway_nbody_%d", instanceId) == sizeof(name))
    {
        mw_panic("name buffer too small for shared memory name\n");
    }

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

    if (sb.st_size < (ssize_t) sizeof(scene_t) || sb.st_size < (ssize_t) nbFindShmemSize(scene->nbody))
    {
        mw_printf("Shared memory segment is impossibly small ("ZU")\n", (size_t) sb.st_size);
        if (shm_unlink(name) < 0)
        {
            mwPerror("Unlink shared memory");
        }

        return NULL;
    }

    return scene;
}

#else

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

#endif /* USE_SHMEM */

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
            OPA_store_int(&scene->attachedLock, pid);
            OPA_store_int(&scene->attachedPID, 0);
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
    printf("Cleanup\n");

    if (g_scene)
    {
        printf("Release scene\n");
        OPA_store_int(&g_scene->attachedPID, 0);
        OPA_store_int(&g_scene->attachedLock, 0);
        g_scene = NULL;
    }
}

static void nbglSigHandler(int sig)
{
    nbglCleanupAttached();
    signal(sig, SIG_DFL);
    raise(sig);
}

static void nbglInstallExitHandlers()
{
    /* TODO: Use sigaction() if available instead */
    signal(SIGINT, nbglSigHandler);
    signal(SIGQUIT, nbglSigHandler);
    signal(SIGINT, nbglSigHandler);
    signal(SIGABRT, nbglSigHandler);
    signal(SIGFPE, nbglSigHandler);
    signal(SIGKILL, nbglSigHandler);
    signal(SIGSEGV, nbglSigHandler);
    signal(SIGUSR1, nbglSigHandler);
    atexit(nbglCleanupAttached);
}

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
        if (scene)
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

