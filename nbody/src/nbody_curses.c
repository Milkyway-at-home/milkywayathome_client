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
  #include <termios.h>
#endif


#if ENABLE_CURSES

static SCREEN* scr[2] = { NULL, NULL };

static struct termios oldTermSet, newTermSet;

/* We want newline -> CR-newline when regular printfs are called in curses */
static void fixTermNewlines()
{
    tcgetattr(0, &oldTermSet);
    newTermSet = oldTermSet;
    newTermSet.c_oflag = newTermSet.c_oflag | ONLCR;
    tcsetattr(0, TCSANOW, &newTermSet);
}

static void resetTermSettings()
{
    tcsetattr(0, TCSANOW, &oldTermSet);
}

int setupCursesOutput(void)
{
    /* Open a separate terminal for stdout and stderr. Then we must
       put the stderr one in shell mode, and then fix the terminal
       settings to translate newlines into newline + carriage return.

       This way we get the progress reported to stdout, and non-broken
       looking stderr printings also.

       The main reason to do this is at the end, the stderr output we
       want will still be in the terminal and not cleared like it
       normally is.

       I'm not entirely sure how / why this works. I had it preserving
       stderr at the end and printing stderr to the curses, but I
       broke it somehow. I think the 2 terminals is so that stderr can
       go to one in shell mode and stdout to one in curses mode.

       Ideally we would also prevent the initial clearing, but I can't
       figure out how to do that. Internet says you can't without
       dropping to a lower level of terminal controls.
     */

    scr[0] = newterm(NULL, stdout, stdin);
    if (!scr[0])
    {
        mw_printf("Failed to get term0\n");
        return 1;
    }
    def_prog_mode();
    //fixTermNewlines();

    scr[1] = newterm(NULL, stderr, stdin);
    if (!scr[1])
    {
        mw_printf("Failed to get term1\n");
        return 1;
    }
    def_prog_mode();
    fixTermNewlines(); /* Fix fprintf(stderr) output while the curses runs */
    mw_reset_prog_mode();


    set_term(scr[0]);       /* Switch to the curses mode terminal */
    mw_reset_shell_mode();


    /* FIXME: It would be good if the stderr went both to the
       non-cleared area and to the curses stuff
     */

    return 0;
}

void cleanupCursesOutput(void)
{
    if (endwin())
    {
        mw_printf("endwin() error\n");
    }

    if (endwin())
    {
        mw_printf("endwin() error\n");
    }

    delscreen(scr[0]);
    delscreen(scr[1]);

    resetTermSettings();
}

#else

int setupCursesOutput()
{
    return 0;
}

void cleanupCursesOutput()
{
}

#endif /* ENABLE_CURSES */


