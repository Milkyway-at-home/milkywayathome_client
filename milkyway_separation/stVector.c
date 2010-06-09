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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "stVector.h"
#include "stCoords.h"

/* Allocate a 2-dimensional matrix of doubles. */
double** matalloc( int nrows, int ncols )
{
    double** mat;
    int i;

    mat = (double**)malloc( nrows * sizeof(double*) );
    for ( i = 0; i < nrows; ++i )
        mat[i] = (double*)malloc( ncols * sizeof(double) );

    return mat;
}


/* Free memory returned from a call to list2dblmat() or matalloc(). */
void matfree( double** mat, int nrows )
{
    int i;

    for ( i = 0; i < nrows; ++i )
        free( mat[i] );

    free( mat );
}


/* Get norm of input vector */
double norm( const double* vec )
{
    return sqrt( vec[0] * vec[0] +
                 vec[1] * vec[1] +
                 vec[2] * vec[2] );
}


/* Normalize input vector */
void normalize( double* vec )
{
    double vnorm = norm( vec );

    vec[0] /= vnorm;
    vec[1] /= vnorm;
    vec[2] /= vnorm;
}


/* Dot product */
double dotp( const double* a, const double* b )
{
    return a[0] * b[0] +
           a[1] * b[1] +
           a[2] * b[2];
}


/* Cross product; stores result in prod */
void crossp( const double* a, const double* b, double* prod )
{
    prod[0] = a[1] * b[2] - a[2] * b[1];
    prod[1] = a[2] * b[0] - a[0] * b[2];
    prod[2] = a[0] * b[1] - a[1] * b[0];
}


/* Angle between two vectors, in the range [0,pi] */
double vecangle( const double* a, const double* b )
{
    double anorm, bnorm, dprod;

    anorm = norm( a );
    bnorm = norm( b );
    dprod = dotp( a, b );

    return acos( dprod / (anorm * bnorm) );
}


/* Determine the rotation matrix needed to transform f into t.  The result is an
   array of 9 elements representing a (flattened) 3x3 matrix.
   Adapted from information at http://www.flipcode.com/documents/matrfaq.html */
void get_transform( const double* f, const double* t, double** mat )
{
    double angle, sin_a;
    double x, y, z, w;
    double axis[3];

    crossp( f, t, axis );
    normalize( axis );

    angle = vecangle( f, t );
    sin_a = sin( angle / 2 );

    x = axis[0] * sin_a;
    y = axis[1] * sin_a;
    z = axis[2] * sin_a;
    w = cos( angle / 2 );

    mat[0][0] = 1 - 2 * (y * y + z * z);
    mat[0][1] =     2 * (x * y - z * w);
    mat[0][2] =     2 * (x * z + y * w);

    mat[1][0] =     2 * (x * y + z * w);
    mat[1][1] = 1 - 2 * (x * x + z * z);
    mat[1][2] =     2 * (y * z - x * w);

    mat[2][0] =     2 * (x * z - y * w);
    mat[2][1] =     2 * (y * z + x * w);
    mat[2][2] = 1 - 2 * (x * x + y * y);
}


/* Transform v by applying the rotation matrix mat */
void do_transform( double* v, double* const* mat )
{
    double newv[3];

    newv[0] = dotp( mat[0], v );
    newv[1] = dotp( mat[1], v );
    newv[2] = dotp( mat[2], v );

    v[0] = newv[0];
    v[1] = newv[1];
    v[2] = newv[2];
}

/* apply coordinate transformations to the given point */
void transform_point(double* point, double** cmat, double* xsun, double* logPoint)
{
    double mcutoff = 11.0;

    xyz_mag(point, mcutoff, logPoint);

    double newx = logPoint[0] - xsun[0];
    double newy = logPoint[1] - xsun[1];
    double newz = logPoint[2] - xsun[2];
    logPoint[0] = newx;
    logPoint[1] = newy;
    logPoint[2] = newz;

    do_transform(logPoint, cmat);
}

