// src/basis.h
#ifndef BASIS_H
#define BASIS_H

// Calculates sin(m*phi) and cos(m*phi)
void basis_sincos(int mmax, double phi, double* c, double* s);

// Calculates both P_l^m(x) and their derivatives
void basis_legendre_deriv(int lmax, double x, double* p, double* dp);

// Calculates Gegenbauer polynomials C_n^alpha(x) for n=0 to nmax
void basis_gegenbauer(int nmax, double alpha, double x, double* C_out);

#endif // BASIS_H