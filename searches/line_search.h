#ifndef FGDO_LINE_SEARCH_H
#define FGDO_LINE_SEARCH_H

extern const char *LS_STR[];
#define LS_SUCCESS	0
#define LS_LOOP1_NAN	1
#define LS_LOOP2_NAN	2
#define LS_LOOP3_NAN	3
#define LS_LOOP1_TOL	4
#define LS_LOOP1_MAX	5
#define LS_LOOP2_MAX	6


int line_search(double* point, double initial_fitness, double* direction, int number_parameters, double** minimum, double* fitness, int *evaluations_done);

#endif
