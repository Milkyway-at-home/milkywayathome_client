/*
 * Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
 * Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
 * and Rensselaer Polytechnic Institute.
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
 * */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "search_arguments.h"


int argument_exists(const char* name, int argc, char** argv) {
	int i;
	for (i = 0; i < argc; i++) if (!strcmp(argv[i], name)) return 1;
	return 0;
}

int get_int_arg(const char* name, int argc, char** argv) {
	int i;
	for (i = 0; i < argc; i++) if (!strcmp(argv[i], name)) return atoi(argv[++i]);
	return -1;
}

long get_long_arg(const char* name, int argc, char** argv) {
	int i;
	for (i = 0; i < argc; i++) if (!strcmp(argv[i], name)) return atol(argv[++i]);
	return -1;
}

double get_double_arg(const char* name, int argc, char** argv) {
	int i;
	for (i = 0; i < argc; i++) if (!strcmp(argv[i], name)) return atof(argv[++i]);
	return -1.0;
}
