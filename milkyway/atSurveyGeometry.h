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

#ifndef ATSURVEYGEOMETRY_H
#define ATSURVEYGEOMETRY_H

#include "milkyway.h"

#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
                             :(A)+(B)*floor(-(A)/(B))):(A))
#define dsign(A,B) ((B)<0.0?-(A):(A))

#define at_stripeSeparation  2.5;
#define at_surveyCenterRa 185.0;
#define at_surveyCenterDec 32.5;
#define at_deg2Rad DPI/180.0;
#define at_rad2Deg 180.0/DPI;

void atGCToEq (
    double amu,  /* IN */
    double anu,  /* IN */
    double* ra,  /* OUT */
    double* dec, /* OUT */
    double anode,    /* IN */
    double ainc  /* IN */
);

void atEqToGal (
    double ra,  /* IN */
    double dec, /* IN */
    double* glong,  /* OUT: Galactic longitude */
    double* glat    /* OUT: Galactic latitude */
);

void atBound (
    double* angle,    /* MODIFIED -- the angle to bound */
    double min,   /* IN -- inclusive minimum value */
    double max    /* IN -- exclusive maximum value */
);

void atBound2(
    double* theta,    /* MODIFIED -- the -90 to 90 angle */
    double* phi   /* MODIFIED -- the 0 to 360 angle */
);

void slaDcc2s ( double v[3], double* a, double* b );

void slaDimxv ( double dm[3][3], double va[3], double vb[3] );

void slaDcs2c ( double a, double b, double v[3] );

void slaDmxv ( double dm[3][3], double va[3], double vb[3] );

double slaDrange ( double angle );

double slaDranrm ( double angle );

void slaEqgal ( double dr, double dd, double* dl, double* db );

void gcToSgr ( double mu, double nu, int wedge, double* lamda, double* beta );

void sgrToGal ( double lamda, double beta, double* l, double* b);

void atSurveyToEq ( double slong, double slat, double* ra, double* dec);

double atEtaFromStripeNumber ( int wedge );

#endif /* ATSURVEYGEOMETRY_H */

