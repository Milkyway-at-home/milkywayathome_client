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

#ifndef _STCOORDS_H_
#define _STCOORDS_H_

#include "milkyway.h"
#include "milkyway_priv.h"


#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
                             :(A)+(B)*floor(-(A)/(B))):(A))
#define dsign(A,B) ((B)<0.0?-(A):(A))

#define at_stripeSeparation  (2.5)
#define at_surveyCenterRa (185.0)
#define at_surveyCenterDec (32.5)
#define at_deg2Rad (DPI/180.0)
#define at_rad2Deg (180.0/DPI)

/* The node of the GC coordinates used in the survey. */
#define NODE_GC_COORDS (at_surveyCenterRa - 90.0)
#define sun_r0 8.5

typedef struct
{
    double ra;
    double dec;
} RA_DEC;

/* wedge, mu, nu, l, b. gc2lb or gc2sgr  */
typedef void (*SGRConversion)(int, double, double, double*, double*);

void gc2lb( int wedge, double mu, double nu, double* l, double* b );
void gc2sgr( int wedge, double mu, double nu, double* l, double* b );

void lbr2xyz( const double* lbr, double* xyz );


#endif /* _STCOORDS_H_ */

