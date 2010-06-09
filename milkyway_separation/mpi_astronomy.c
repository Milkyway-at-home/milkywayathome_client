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
    *   Includes for astronomy
 ********/
#include "astronomy_worker.h"

/********
    *   Includes for FGDO
 ********/
#include "../searches/asynchronous_newton_method.h"
#include "../searches/asynchronous_genetic_search.h"
#include "../searches/asynchronous_particle_swarm.h"
#include "../searches/asynchronous_search.h"
#include "../searches/bounds.h"
#include "../searches/synchronous_gradient_descent.h"
#include "../searches/synchronous_newton_method.h"
#include "../searches/search_arguments.h"
#include "../evaluation/evaluator.h"
#include "../evaluation/mpi_evaluator.h"


int main(int number_arguments, char** arguments)
{
    double* point, *range, *min_bound, *max_bound;
    int number_parameters, integral_parameter_length, integral_results_length, likelihood_parameter_length, likelihood_results_length;
    int i;

    /********
        *   Initialize the mpi_evaluator
     ********/
    mpi_evaluator__init(&number_arguments, &arguments);

    sprintf(star_points_file, "stars.txt");
    sprintf(astronomy_parameters_file, "parameters.txt");
    for (i = 0; i < number_arguments; i++)
    {
        if (!strcmp(arguments[i], "-stars"))
        {
            sprintf(star_points_file, "%s", arguments[++i]);
        }
        else if (!strcmp(arguments[i], "-parameters"))
        {
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

    evaluator__init_integral(integral_f, integral_parameter_length, integral_compose, integral_results_length);
    evaluator__init_likelihood(likelihood_f, likelihood_parameter_length, likelihood_compose, likelihood_results_length);

    /********
        *   Start the mpi_evaluator
     ********/

    mpi_evaluator__start();

    get_min_parameters(ap, &min_bound);
    get_max_parameters(ap, &max_bound);
    get_search_parameters(ap, &point);
    get_step(ap, &range);

    /********
        *   Start the search
     ********/
    if (argument_exists("-asynch", number_arguments, arguments))
    {
        BOUNDS* bounds;
        if (argument_exists("-nm", number_arguments, arguments))    register_search(get_asynchronous_newton_method());
        else if (argument_exists("-gs", number_arguments, arguments))   register_search(get_asynchronous_genetic_search());
//      if (argument_exists("-de", number_arguments, arguments))    register_search(get_asynchronous_differential_evolution());
        else if (argument_exists("-ps", number_arguments, arguments))   register_search(get_asynchronous_particle_swarm());
//      if (argument_exists("-sx", number_arguments, arguments))    register_search(get_asynchronous_simplex());
        else
        {
            printf("Search not specified.\n");
            return 0;
        }
        new_bounds(&bounds, number_parameters, min_bound, max_bound, NULL);

        printf("doing asynchronous search\n");
        asynchronous_search(number_arguments, arguments, number_parameters, point, range, bounds);
    }
    else
    {
        if (argument_exists("-nm", number_arguments, arguments))    synchronous_newton_method(number_arguments, arguments, number_parameters, point, range);
        else if (argument_exists("-gd", number_arguments, arguments))   synchronous_gradient_descent(number_arguments, arguments, number_parameters, point, range);
        else if (argument_exists("-cgd", number_arguments, arguments))  synchronous_conjugate_gradient_descent(number_arguments, arguments, number_parameters, point, range);
        else
        {
            printf("Search not specified.\n");
            return 0;
        }
    }
    MPI_Abort(MPI_COMM_WORLD, 0);
    return 0;
}
