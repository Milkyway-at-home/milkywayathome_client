/*
Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "settings.h"

char working_directory[FILENAME_SIZE] = ".";

char* get_working_directory() {
	return working_directory;
}

void set_working_directory(char* wd) {
	sprintf(working_directory, "%s", wd);
	working_directory[strlen(wd)] = '\0';
}

void remove_arg(int *target, int *count, char ***values) {
	int i, current;
	char** new_values;

	new_values = (char**)malloc(sizeof(char*) * ((*count) - 1));
	current = 0;
	for (i = 0; i < (*count); i++) {
		if (i == (*target)) continue;
		new_values[current] = (*values)[i];
		current++;
	}

	printf("freeing %d of %d\n", (*target), (*count));
	printf("\tvalue: %s\n", (*values)[(*target)]);
	free( (*values)[(*target)] );
	free( (*values) );

	(*values) = new_values;
	(*target)--;
	(*count)--;
}
