#ifndef STPROB_H
#define STPROB_H

double r2mag(double r);
double mag2r(double g);
double reff(double kr);

double Jacob( const double* a, const double* b, double sint, double xp, int verb );

double stPbxConvolved(const double* coordpar, const double* bpars, int numpoints, int wedge);
double stPsgConvolved(const double* coordpar, const double* spars, int numpoints, int wedge, int sgr_coordinates);
double stPbxFunction( const double* coordpar, const double* bpars);
double stPsgFunction( const double* coordpar, const double* spars, int wedge, int sgr_coordiates);
double stPbx( const double* coordpar, const double* bpars);
double stPsg( const double* coordpar, const double* spars, int wedge, int sgr_coordinates);

double backgroundConvolve(double g, int wedge);
double streamConvolve(double g, int wedge, int sgr_coordinates);

#endif
