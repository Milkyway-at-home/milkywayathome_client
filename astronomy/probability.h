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

#ifndef STPROB_H
#define STPROB_H

double r2mag(double r);
double mag2r(double g);
double reff(double kr);

double Jacob( const double* a, const double* b, double sint, double xp, int verb );

double stPbxConvolved(const double* coordpar, const double* bpars, int wedge, int numpoints);
double stPsgConvolved(const double* coordpar, const double* spars, int wedge, int numpoints, int sgr_coordinates);
double stPbxFunction( const double* coordpar, const double* bpars);
double stPsgFunction( const double* coordpar, const double* spars, int wedge, int sgr_coordiates);
double stPbx( const double* coordpar, const double* bpars);
double stPsg( const double* coordpar, const double* spars, int wedge, int sgr_coordinates);

double backgroundConvolve(double g, int wedge);
double streamConvolve(double g, int wedge, int sgr_coordinates);

int prob_ok (int n, double* p);
void prob_ok_init ();

#endif
