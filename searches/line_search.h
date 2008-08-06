#ifndef FGDO_LINE_SEARCH_H
#define FGDO_LINE_SEARCH_H

int synchronous_line_search(double* point, double initial_fitness, double* direction, int number_parameters, double** minimum, double* fitness);

#endif
