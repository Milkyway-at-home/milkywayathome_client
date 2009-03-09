#ifndef FGDO_OUTLIERS_H
#define FGDO_OUTLIERS_H

#include "population.h"

double distance_error_from(POPULATION *p, double fitness, double *parameters);

void get_distance_errors(POPULATION *p, double **errors);
void get_error_stats(double *errors, int error_size, double *min_error, double *max_error, double *median_error, double *average_error);
void population_error_stats(POPULATION *p, double *min_error, double *max_error, double *median_error, double *average_error);

double get_error2(POPULATION *p, int position, double fitness);
double get_distance2(POPULATION *p, int position, double *parameters);
double distance_error2(POPULATION *p, int position, double fitness, double *parameters);

int remove_outliers(POPULATION *p, double range);
int remove_outliers_incremental(POPULATION *p, double range);
int remove_outliers_sorted(POPULATION *p, double range);
#endif
