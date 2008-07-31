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
