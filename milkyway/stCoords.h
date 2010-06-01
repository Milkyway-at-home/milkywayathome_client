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

#ifndef STCOORDS_H
#define STCOORDS_H 1

void lbr2xyz( const double* lbr, double* xyz );
void xyz2lbr( const double* xyz, double* lbr );

int lbr2stream( const double* lbr, const double* spars, double* stream, int verb );
int xyz2stream( const double* xyz, const double* spars, double* stream, int verb );

void stream2lbr( const double* stream, const double* spars, double* lbr );
void stream2xyz( const double* stream, const double* spars, double* xyz );

double wedge_eta ( int wedge );
double wedge_incl( int wedge );
double get_node();

void gc2lb( int wedge, double mu, double nu, double* l, double* b );

void stripe_normal ( int wedge, double* xyz);

void lbToXyz ( double l, double b, double* xyz );

void lbToXyz(double l, double b, double* xyz);

void xyz_mag(double* point, double offset, double* logPoint);

void xyz2lbg(double* point, double offset, double* logPoint);

void sgr_stripe_normal(int wedge, double* xyz);

#endif /* STCOORDS_H */

