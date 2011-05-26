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
#include "milkyway_util.h"

typedef struct
{
    int fullscreen;
    int width;
    int height;

    /* pid_t pid */
    /* char* nbody bin */

} VisArgs;

static int handleVisArguments(int argc, const char* argv[], VisArgs* visOut)
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

        POPT_AUTOHELP
        POPT_TABLEEND
    };

    /* TODO: Check project prefs */

    context = poptGetContext(argv[0], argc, argv, options, 0);

    if (argc < 2)
    {
        poptPrintUsage(context, stderr, 0);
        poptFreeContext(context);
        return TRUE;
    }

    if (mwReadArguments(context))
    {
        poptPrintHelp(context, stderr, 0);
        failed = TRUE;
    }

    poptFreeContext(context);

    *visOut = visArgs;
    return failed;
}

int main(int argc, const char* argv[])
{
    int rc;

    rc = nbodyGLSetup(&argc, argv);
    if (rc)
    {
        warn("Failed to setup GL\n");
        mw_finish(rc);
    }

    nbodyInitDrawState();
    glutMainLoop();

    return 0;
}


