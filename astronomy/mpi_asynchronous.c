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

#include "mpi.h"

/********
	*	Includes for astronomy
 ********/
#include "astronomy_worker.h"

/********
	*	Includes for FGDO
 ********/
#include "../searches/asynchronous_newton_method.h"
#include "../evaluation/evaluator.h"
#include "../evaluation/mpi_evaluator.h"
#include "../evaluation/mpi_search_manager.h"


int main(int number_arguments, char **arguments){
	double *point, *range, *min_bound, *max_bound;
	int number_parameters, integral_parameter_length, integral_results_length, likelihood_parameter_length, likelihood_results_length;
	int i;

	/********
		*	Initialize the mpi_evaluator
	 ********/
	mpi_evaluator__init(&number_arguments, &arguments);

	sprintf(star_points_file, "stars.txt");
	sprintf(astronomy_parameters_file, "parameters.txt");
	for (i = 0; i < number_arguments; i++) {
		if (!strcmp(arguments[i], "-stars")) {
			sprintf(star_points_file, "%s", arguments[++i]);
		} else if (!strcmp(arguments[i], "-parameters")) {
			sprintf(astronomy_parameters_file, "%s", arguments[++i]);
		}
	}

	mpi_evaluator__read_data(read_data);

	number_parameters = get_optimized_parameter_count(ap);
	integral_parameter_length = number_parameters;
	integral_results_length = 1 + ap->number_streams;
	likelihood_parameter_length = 1 + ap->number_streams;
	likelihood_results_length = 2;

	printf("number_parameters: %d\n", number_parameters);
	printf("integral_parameter_length: %d, integral_results_length: %d\n", integral_parameter_length, integral_results_length);
	printf("likelihood_parameter_length: %d, likelihood_results_length: %d\n", likelihood_parameter_length, likelihood_results_length);

	evaluator__init_integral(integral_f, integral_parameter_length, integral_compose, integral_results_length);
	evaluator__init_likelihood(likelihood_f, likelihood_parameter_length, likelihood_compose, likelihood_results_length);

	/********
		*	Start the mpi_evaluator
	 ********/

	printf("starting...\n");
	mpi_evaluator__start();

	get_min_parameters(ap, &min_bound);
	get_max_parameters(ap, &max_bound);
	get_search_parameters(ap, &point);
	get_step(ap, &range);

	/********
		*	Register searches and parse the search parameters
	 ********/
	register_search(get_asynchronous_newton_method());

	/********
		*	Start the search
	 ********/
	mpi_asynchronous_search(number_arguments, arguments, number_parameters, point, range, min_bound, max_bound);

	return 0;
}
