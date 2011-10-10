/*
 * Copyright (c) 2011 Matthew Arsenault
 * Copyright (c) 2011 Rensselaer Polytechnic Institute
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

#include "nbody_priv.h"
#include "nbody_curses.h"

#if ENABLE_CURSES
  #include <errno.h>
#endif


#if ENABLE_CURSES


static int stderrBack = -1;       /* Backup of stderr from dup */
static char stderrTmpFile[128];   /* Name of stderr dump file */
static mwbool cursesSetup = FALSE;


int nbSetupCursesOutput(void)
{
    int stderrTmp;
    int redirect; /* Will become stderr, let it "leak" */

    /* Since we use stderr for everything, use stdout for curses stuff.

       We don't want the printed stuff to disappear at the end after
       endwin() clears. Redirect stderr to a temporary file and print
       that at the end.

       Ideally we would also prevent the initial clearing, but I can't
       figure out how to do that. Internet says you can't without
       dropping to a lower level of terminal controls.
     */

    snprintf(stderrTmpFile, sizeof(stderrTmpFile), "stderr_tmp_%d", getpid());

    fflush(stderr);
    stderrBack = dup(fileno(stderr));
    if (stderrBack < 0)
    {
        mwPerror("dup() error on STDERR_FILENO");
        return errno;
    }

    stderrTmp = open(stderrTmpFile, O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (stderrTmp < 0)
    {
        mwPerror("Failed to get stderr dump file '%s'", stderrTmpFile);
        close(stderrBack);
        return errno;
    }

    redirect = dup2(stderrTmp, fileno(stderr));
    if (redirect < 0)
    {
        mwPerror("dup2() error to stderr redirect");
        close(stderrBack);
        return errno;
    }

    if (close(stderrTmp) < 0)
    {
        printf("Error closing stderrTmp: %s\n", strerror(errno));
        /* Try to restore stderr to console */
        dup2(stderrBack, fileno(stderr));

        close(stderrBack);
        return errno;
    }

    initscr();
    cursesSetup = TRUE;

    return 0;
}

void nbCleanupCursesOutput(void)
{
    FILE* f;
    size_t readSize;
    char buf[4096];

    if (!cursesSetup)
    {
        return;
    }

    if (endwin() != OK)
    {
        mw_printf("endwin() error\n");
    }

    /* Restore stderr to the console */
    fflush(stderr);
    if (dup2(stderrBack, fileno(stderr)) < 0)
    {
        printf("Error restoring stderr: %s\n", strerror(errno));
    }

    if (close(stderrBack) < 0)
    {
        mwPerror("Error closing stderrBack");
    }

    /* You can't seek on the redirected fd, so reopen the stderr tmp
     * file so we can print it. */
    f = fopen(stderrTmpFile, "r");
    if (!f)
    {
        mwPerror("Error reopening stderr log '%s'", stderrTmpFile);
        return;
    }

    /* Copy file contents to stderr */
    while (!feof(f))
    {
        readSize = fread(buf, 1, sizeof(buf), f);
        fwrite(buf, 1, readSize, stderr);
    }

    fclose(f);
    remove(stderrTmpFile);
}

#else

int nbSetupCursesOutput()
{
    return 0;
}

void nbCleanupCursesOutput()
{
}

#endif /* ENABLE_CURSES */


