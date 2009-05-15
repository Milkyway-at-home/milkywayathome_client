/*
Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef NUMERICALINTEGRATION_H
#define NUMERICALINTEGRATION_H

void gaussLegendre__float(float x1, float x2, float x[], float w[], int n);
void gaussLegendre(double x1, double x2, double x[], double w[], int n);

void setWeights(int numpoints);

void freeWeights();

double qgaus(double (*func)(double, int), double xm, double xr, int wedge, int numpoints);

double qgaus_stream(double (*func)(double, int, int), double xm, double xr, int wedge, int numpoints, int sgr_coordinates);

#endif
