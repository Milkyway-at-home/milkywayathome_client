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

#ifndef FGDO_ASYNCHRONOUS_PSO_H
#define FGDO_ASYNCHRONOUS_PSO_H

#include <stdio.h>

#include "asynchronous_search.h"
#include "bounds.h"
#include "population.h"
#include "redundancies.h"

#define NEWTON_ERROR_RANGE 1
#define NEWTON_UPDATE_RANGE 2
#define NEWTON_LINE_SEARCH 3

#define NMS_HESSIAN 1
#define NMS_LINE_SEARCH 2

typedef struct particle_swarm_optimization {
	int current_particle, size;
	int number_parameters;
	double w, c0, c1, c2;
	double redundancy_rate;
	long analyzed;

	BOUNDS *bounds;

	double global_best_fitness, *global_best;
	POPULATION *particles;
	POPULATION *velocities;
	POPULATION *local_best;

	REDUNDANCIES *redundancies;
} PARTICLE_SWARM_OPTIMIZATION;

ASYNCHRONOUS_SEARCH* get_asynchronous_particle_swarm();

int create_particle_swarm(char* search_name, int number_arguments, char** arguments, int number_parameters, double *point, double *range, BOUNDS *bounds);
int read_particle_swarm(char* search_name, void** search_data);
int checkpoint_particle_swarm(char* search_name, void* search_data);
int pso_insert_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS *sp);
int pso_generate_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS *sp);
#endif
