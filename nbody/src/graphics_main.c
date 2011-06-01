/*
Copyright (C) 2011  Matthew Arsenault

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

#include "nbody_gl.h"
#include "nbody_graphics.h"
#include "milkyway_util.h"

static const VisArgs defaultVisArgs =
{
    /* .fullscreen */     FALSE,
    /* .width      */     0,
    /* .height     */     0,
    /* .monochrome */     FALSE,
    /* .notUseGLPoints */ FALSE
};


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
            0, "Start in fullscreen mode", NULL
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
            "not-use-gl-points", 'p',
            POPT_ARG_NONE, &visArgs.notUseGLPoints,
            0, "Use faster but possibly uglier drawing", NULL
        },

        POPT_AUTOHELP
        POPT_TABLEEND
    };

    /* TODO: Check project prefs */

    visArgs = defaultVisArgs;
    context = poptGetContext(argv[0], argc, argv, options, 0);

    if (mwReadArguments(context))
    {
        poptPrintHelp(context, stderr, 0);
        failed = TRUE;
    }

    poptFreeContext(context);

    *visOut = visArgs;
    return failed;
}

static int nbodyGraphicsInit()
{
#if BOINC_APPLICATION
    if (boinc_parse_init_data_file())
    {
        warn("Error parsing init data file\n");
    }

    if (mwBoincInit(MW_GRAPHICS))
        return 1;
#endif /* BOINC_APPLICATION */

    return 0;
}

int main(int argc, char* argv[])
{
    int rc;
    VisArgs flags;

    if (nbodyGraphicsInit())
        return 1;

    if (connectSharedScene())
        return 1;

    glutInit(&argc, argv);

    if (handleVisArguments(argc, (const char**) argv, &flags))
        return 1;

    rc = nbodyGLSetup(&flags);
    if (rc)
    {
        warn("Failed to setup GL\n");
        mw_finish(rc);
    }

    glutMainLoop();

    return 0;
}

