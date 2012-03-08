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

#include "nbody_gl.h"
#include "nbody_graphics.h"
#include "milkyway_util.h"

#include <signal.h>

static const VisArgs defaultVisArgs =
{
    /* .fullscreen      */ FALSE,
    /* .plainFullscreen */ FALSE,
    /* .width           */ 0,
    /* .height          */ 0,
    /* .monochrome      */ FALSE,
    /* .notUseGLPoints  */ FALSE,
    /* .originCenter    */ FALSE,
    /* .noFloat         */ FALSE,
    /* .pid             */ 0,
    /* .file            */ NULL,
    /* .instanceId      */ -1
};

static void freeVisArgs(VisArgs* args)
{
    free(args->file);
    args->file = NULL;
}

static int handleVisArguments(int argc, const char** argv, VisArgs* visOut)
{
    poptContext context;
    int failed = FALSE;
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
            "monochromatic", 'm',
            POPT_ARG_NONE, &visArgs.monochrome,
            0, "All particles have same color", NULL
        },

        {
            "not-use-gl-points", 'n',
            POPT_ARG_NONE, &visArgs.notUseGLPoints,
            0, "Use faster but possibly uglier drawing", NULL
        },

        {
            "origin-center", 'o',
            POPT_ARG_NONE, &visArgs.originCenter,
            0, "Focus on the galactic center instead of system's center of mass", NULL
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

    *visOut = visArgs;
    return failed;
}

static int nbodyGraphicsInit()
{
#if BOINC_APPLICATION
    if (boinc_parse_init_data_file())
    {
        mw_printf("Error parsing init data file\n");
    }

    if (mwBoincInit(MW_GRAPHICS))
        return 1;
#endif /* BOINC_APPLICATION */

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

    scene = mwCalloc(sizeof(scene_t) + lnCount * sizeof(FloatPos), sizeof(char));
    r = scene->rTrace;

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

static scene_t* g_scene = NULL;

static void cleanupAttachedCount(void)
{
    printf("Cleanup\n");

    if (g_scene)
    {
        printf("Decrement\n");
        OPA_decr_int(&g_scene->attachedCount);
        g_scene = NULL;
    }
}

static void sigHandler(int sig)
{
    cleanupAttachedCount();
    signal(sig, SIG_DFL);
    raise(sig);
}

static void installExitHandlers()
{
    /* TODO: Use sigaction() if available instead */
    signal(SIGINT, sigHandler);
    signal(SIGQUIT, sigHandler);
    signal(SIGINT, sigHandler);
	signal(SIGABRT, sigHandler);
	signal(SIGFPE, sigHandler);
	signal(SIGKILL, sigHandler);
	signal(SIGSEGV, sigHandler);
	signal(SIGUSR1, sigHandler);
    atexit(cleanupAttachedCount);
}

int main(int argc, const char* argv[])
{
    int rc;
    VisArgs flags;
    scene_t* scene = NULL;

    if (nbodyGraphicsInit())
        return 1;

    if (handleVisArguments(argc, (const char**) argv, &flags))
    {
        freeVisArgs(&flags);
        return 1;
    }

    if (!flags.file)
    {
        scene = nbConnectSharedScene(flags.instanceId);
        if (!scene)
        {
            freeVisArgs(&flags);
            return 1;
        }

        if (nbCheckConnectedVersion(scene))
        {
            freeVisArgs(&flags);
            return 1;
        }

        installExitHandlers();
        OPA_incr_int(&scene->attachedCount);
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

    rc = nbRunGraphics(scene, &flags);

    if (flags.file)
    {
        free(scene);
    }

    freeVisArgs(&flags);
    mw_finish(rc);

    return rc;
}

