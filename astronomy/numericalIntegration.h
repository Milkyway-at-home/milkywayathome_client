#ifndef NUMERICALINTEGRATION_H
#define NUMERICALINTEGRATION_H

void gaussLegendre(double x1, double x2, double x[], double w[], int n);

void setWeights(double x[], double w[], int numpoints);

void freeWeights();

double qgaus(double (*func)(double, int), double xm, double xr, int wedge, int numpoints);

double qgaus_stream(double (*func)(double, int, int), double xm, double xr, int wedge, int numpoints, int sgr_coordinates);

#endif
