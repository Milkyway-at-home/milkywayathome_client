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

#ifndef FGDO_POPULATION_H
#define FGDO_POPULATION_H

#include <stdio.h>
#include "../util/settings.h"

typedef struct population {
	int size;
	int max_size;

	int number_parameters;
	double** individuals;
	double* fitness;
	char** os_names;
	char** app_versions;
} POPULATION;


void fwrite_population_statistics(FILE* file, POPULATION* population);

int new_population(int max_size, int number_parameters, POPULATION** population);
void free_population(POPULATION* population);

double get_best_individual(POPULATION *p, double *best_individual);
void get_population_statistics(POPULATION *p, double *best_fitness, double *average_fitness, double *median_fitness, double *worst_fitness, double *st_dev);

int population_contains(POPULATION* population, double fitness, double *point);
int individual_exists(POPULATION* population, int position);

void insert_individual(POPULATION* population, int position, double* parameters, double fitness);
void insert_individual_info(POPULATION* population, int position, double* parameters, double fitness, char *os_name, char *app_version);
void insert_incremental(POPULATION* population, double* parameters, double fitness);
void insert_incremental_info(POPULATION* population, double* parameters, double fitness, char *os_name, char *app_version);
void insert_sorted(POPULATION* population, double* parameters, double fitness);
void insert_sorted_info(POPULATION* population, double* parameters, double fitness, char *os_name, char *app_version);

void remove_individual(POPULATION* population, int position);
void remove_incremental(POPULATION* population, int position);
void remove_sorted(POPULATION* population, int position);

void get_n_distinct(POPULATION *population, int number_parents, POPULATION **n_distinct);

int fread_population(FILE* file, POPULATION** p);
int read_population(char path[512], POPULATION** p);

int fwrite_individual(FILE* file, POPULATION* population, int position);
int fwrite_population(FILE* file, POPULATION* population);
int write_population(char path[512], POPULATION* population);

#endif
