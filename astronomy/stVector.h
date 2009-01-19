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

#ifndef SPVECTOR_H
#define SPVECTOR_H 1

double** matalloc( int nrows, int ncols );
void matfree( double** mat, int nrows );

double norm( const double* vec );
void normalize( double* vec );

double dotp( const double* a, const double* b );
void crossp( const double* a, const double* b, double* prod );

double vecangle( const double* a, const double* b );

void get_transform( const double* f, const double* t, double** mat );

void do_transform( double* v, double* const* mat );

void transform_point( double* point, double** cmat, double* xsun, double* logPoint ); 

#endif /* SPVECTOR_H */
