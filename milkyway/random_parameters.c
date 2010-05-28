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

#include <stdio.h>
#include <time.h>

#include "parameters.h"
#include "../searches/recombination.h"


int main(int argc, char** argv) {
	double *min_p, *max_p, *rand_p;
	int retval, number_parameters;
	time_t current_time;

	time(&current_time);
	srand48(current_time);

	ASTRONOMY_PARAMETERS *ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	retval = read_astronomy_parameters(argv[1], ap);
	if (retval) {
		fprintf(stderr, "ERROR reading parameters file: %s\n", argv[2]);
		exit(0);
	}
	number_parameters = get_optimized_parameter_count(ap);
	get_min_parameters(ap, &min_p);
	get_max_parameters(ap, &max_p);
	rand_p = (double*)malloc(sizeof(double) * number_parameters);
	random_recombination(min_p, max_p, number_parameters, rand_p);

	set_astronomy_parameters(ap, rand_p);
	retval = write_astronomy_parameters(argv[2], ap);
	if (retval) {
		fprintf(stderr, "ERROR writing parameters file: %s\n", argv[2]);
		exit(0);
	}

	return 0;
}
