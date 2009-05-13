#ifndef CPU_COORDS_H
#define CPU_COORDS_H

#include "../astronomy/parameters.h"

void gc_to_gal(int wedge, double amu_rad, double anu_rad, double *glong, double *glat);
void cpu__gc_to_gal(int wedge, INTEGRAL *integral, double **cpu__sinb, double **cpu__sinl, double **cpu__cosb, double **cpu_cosl); 

#endif
