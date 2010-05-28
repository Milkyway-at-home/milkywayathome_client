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

#ifndef ST_CNUM_H
#define ST_CNUM_H

typedef struct {
	double real;
	double imagine;
} cnum;

/* complex print function */
void PrintCnum(cnum *x, int n);

/* complex add functions */
cnum CnumAdd(cnum c1, cnum c2);
cnum CnumAddD(cnum c1, double d);

/* complex subtraction functions */
cnum CnumSub(cnum c1, cnum c2);
cnum CnumSubD(cnum c1, double d);
cnum CnumDSub(double d, cnum c1);

/* complex multiplication functions */
cnum CnumMult(cnum c1, cnum c2);
cnum CnumMultD(cnum c1, double d);

/* complex division functions */
cnum CnumDiv(cnum c1, cnum c2);
cnum CnumDivD(cnum c1, double d);
cnum CnumDDiv(double x, cnum c1);

/* complex square root and cube root functions */
cnum CnumSqrt(cnum c1);
cnum CnumCbrt(cnum c1, int verb);

#endif
