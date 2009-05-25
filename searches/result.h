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

#ifndef FGDO_RESULT
#define FGDO_RESULT

#include <stdio.h>

int fwrite_gpu_result(FILE *file, int number_parameters, double **hessian, double *gradient, double initial_fitness, double *inital_parameters, double result_fitness, double *result_parameters, int number_evaluations, char *metadata);
int write_gpu_result(char *filename, int number_parameters, double **hessian, double *gradient, double initial_fitness, double *inital_parameters, double result_fitness, double *result_parameters, int number_evaluations, char *metadata);

int fread_gpu_result(FILE *file, int number_parameters, double **hessian, double *gradient, double *initial_fitness, double *initial_parameters, double *result_fitness, double *result_parameters, int *number_evaluations, char *metadata);
int fread_gpu_result(char *filename, int number_parameters, double **hessian, double *gradient, double *initial_fitness, double *initial_parameters, double *result_fitness, double *result_parameters, int *number_evaluations, char *metadata);


int fwrite_cpu_result(FILE *file, int number_parameters, double *parameters, double likelihood, char *metadata);
int write_cpu_result(char *filename, int number_parameters, double *parameters, double likelihood, char *metadata);

int fread_cpu_result(FILE *file, int number_parameters, double *parameters, double fitness, char *metadata);
int read_cpu_result(char *filename, int number_parameters, double *parameters, double fitness, char *metadata);

int fread_cpu_result__realloc(FILE *file, int *number_parameters, double **parameters, double *fitness, char *metadata);
int read_cpu_result__realloc(const char *filename, int *number_parameters, double **parameters, double *fitness, char *metadata);

#endif
