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

#ifndef _NBODY_CURSES_H_
#define _NBODY_CURSES_H_

#include "nbody_config.h"


#if ENABLE_CURSES
  #include <ncurses.h>
#endif


#ifdef __cplusplus
extern "C" {
#endif


#if ENABLE_CURSES
  #define mw_mvprintw(y, x, ...) mvprintw(y, x, ##__VA_ARGS__)
  #define mw_printw(...) printw(__VA_ARGS__)
  #define mw_initscr() initscr()
  #define mw_refresh() refresh()
  #define mw_def_shell_mode() def_shell_mode()
  #define mw_def_prog_mode() def_prog_mode()
  #define mw_reset_prog_mode() reset_prog_mode()
  #define mw_reset_shell_mode() reset_shell_mode()
  #define mw_getch() getch()
#else
  #define mw_printw mw_printf
  #define mw_mvprintw(y, x, ...) printf(__VA_ARGS__)
  #define mw_initscr()
  #define mw_delscreen()
  #define mw_refresh()
  #define mw_def_shell_mode()
  #define mw_def_prog_mode()
  #define mw_reset_prog_mode()
  #define mw_reset_shell_mode()
  #define mw_getch()
#endif /* ENABLE_CURSES */


int nbSetupCursesOutput(void);
void nbCleanupCursesOutput(void);


#ifdef __cplusplus
}
#endif

#endif /* _NBODY_CURSES_H_ */

