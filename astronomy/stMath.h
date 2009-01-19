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

#ifndef STMATH_H
#define STMATH_H

#include "stCnum.h"
#include "stVector.h"

double stCEval(cnum a, cnum b, cnum c, cnum x);
double stQEval(double a, double b, double c, double d, cnum x);
int DCubic(double a, double b, double c, cnum *roots, int verb);
int CnumCubic( cnum a, cnum b, cnum c, cnum *roots, int verb);
int quartic(double a, double b, double c, double d, cnum u[], int flag[], int verb);

int stRoot3(double a2, double a1, double a0,double *ans, int verb);
int stRoot4(double a3, double a2, double a1, double a0, double r[], int flag[], int verb);

double sum(double array[], int size);
int min(double array[], int flag[], int size, int numgood);

double fact(int d);
void Makepoints();

#endif
