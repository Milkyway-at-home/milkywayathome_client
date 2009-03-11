#ifndef FGDO_REDUNDANCY_H
#define FGDO_REDUNDANCY_H

#include "stdio.h"

typedef struct redundancy REDUNDANCY;

struct redundancy {
	double *parameters;
	double *velocity;
	double fitness;

	REDUNDANCY *next;
};

int parameters_match(int number_parameters, double *p1, double *p2);
int fitness_match(double f1, double f2);

void free_redundancy(REDUNDANCY **r);
void new_redundancy(REDUNDANCY **r, double fitness, int number_parameters, double *parameters, double *velocity);

int fread_redundancy(FILE *file, REDUNDANCY **r, int number_parameters, int *particle);
int fwrite_redundancy(FILE *file, REDUNDANCY *r, int number_parameters, int particle);
#endif
