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

#include "evaluator.h"
#include "mpi_search_manager.h"
#include "search_manager.h"
#include "../searches/search_parameters.h"


void start_mpi_search_manager(int argc, char** argv) {
	SEARCH_PARAMETERS **sp;
	int generated, generation_rate, i;

	printf("[mpi_search_manager]: init.\n");
	init_search_manager(argc, argv);
	printf("[mpi_search_manager]: successful.\n");

	generated = 1;
	generation_rate = get_generation_rate();
	sp = (SEARCH_PARAMETERS**)malloc(sizeof(SEARCH_PARAMETERS*) * generation_rate);

	while (generated > 0) {
		generated = generate_search_parameters(sp);
		for (i = 0; i < generated; i++) {
			sp[i]->fitness = evaluate(sp[i]->parameters);
			insert_search_parameters(sp[i]);
			free_search_parameters(sp[i]);
		}
	}
}
