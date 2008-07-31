#ifndef NUMERICALINTEGRATION_H
#define NUMERICALINTEGRATION_H

void gaussLegendre(double x1, double x2, double x[], double w[], int n);

double qgaus(double (*func)(double, int), double a, double b, int n, int wedge);

#endif
