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

#ifndef FGDO_SETTINGS_H
#define FGDO_SETTINGS_H

#define FILENAME_SIZE		2048
#define SEARCH_QUALIFIER_SIZE	64
#define SEARCH_NAME_SIZE	512
#define METADATA_SIZE		2048
#define BOINC_SEARCH_PATH "./boinc_testing/"

char* get_working_directory();
void set_working_directory(char* working_directory);

void remove_arg(int *target, int *count, char ***values);

#endif
