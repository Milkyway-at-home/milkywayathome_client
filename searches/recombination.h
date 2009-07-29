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

#ifndef H_GEM_RECOMBINATION
#define H_GEM_RECOMBINATION

void mutate(double* parent, double* min_parameters, double* max_parameters, int number_parameters, double *parameters);

void random_recombination(int number_parameters, double* min_parameters, double* max_parameters, double* result);

void range_recombination(int number_parameters, double* point, double* range, double *result);
double random_linear_recombination(int number_parameters, double min, double max, double* initial, double* step, double *parameters);

void average_recombination(double** parents, int number_parents, int number_parameters, double *result);
void higher_recombination(double** parents, int number_parents, int number_parameters, double *result);
void lower_recombination(double** parents, int number_parents, int number_parameters, double *result);

double simplex_recombination(double** parents, double* fitness, int number_parents, int number_parameters, double ls_center, double ls_outside, double *result);

double* binomial_recombination(double** parents, int number_parents, int number_parameters, double crossover_rate, double crossover_scale);
double* exponential_recombination(double** parents, int number_parents, int number_parameters, double crossover_rate, double crossover_scale);

double* get_pair_sum(double **individuals, int number_individuals, int number_parameters, int number_pairs, double scale);
double* get_dir_sum(double **individuals, double *fitness, int number_individuals, int number_parameters, int number_pairs, double scale);

#endif
