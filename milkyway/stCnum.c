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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "stCnum.h"
//#include "values.h"

#ifndef _WIN32
#define pi M_PI
#else
#define pi 3.14159265358979323846
#endif
#define deg (180.0/pi)

void PrintCnum(cnum* x, int n)
{
    double r,i;
    r = x->real;
    i = x->imagine;

    printf("%g + i * %g",r,i);
    if(n == 1)
        printf("\n");
}

cnum CnumAdd(cnum c1, cnum c2)
{
    cnum temp;
    temp.real = c1.real + c2.real;
    temp.imagine = c1.imagine + c2.imagine;

    return temp;
}

cnum CnumAddD(cnum c1, double d)
{
    cnum temp;
    temp.real = c1.real + d;
    temp.imagine = c1.imagine;

    return temp;
}

cnum CnumSub(cnum c1, cnum c2)
{
    cnum temp;
    temp.real = c1.real - c2.real;
    temp.imagine = c1.imagine - c2.imagine;

    return temp;
}

cnum CnumSubD(cnum c1, double d)
{
    cnum temp;
    temp.real = c1.real - d;
    temp.imagine = c1.imagine;

    return temp;
}

cnum CnumDSub(double d, cnum c1)
{
    cnum temp;
    temp.real = d - c1.real;
    temp.imagine = -c1.imagine;

    return temp;
}

cnum CnumMult(cnum c1, cnum c2)
{
    cnum temp;
    temp.real = c1.real * c2.real - c1.imagine * c2.imagine;
    temp.imagine = c1.real * c2.imagine + c1.imagine * c2.real;

    return temp;
}

cnum CnumMultD(cnum c1, double d)
{
    cnum temp;
    temp.real = c1.real * d;
    temp.imagine = c1.imagine * d;

    return temp;
}

cnum CnumDiv(cnum c1, cnum c2)
{
    double d;
    cnum temp;

    d = c2.real * c2.real + c2.imagine * c2.imagine;
    temp.real = (c1.real * c2.real - c1.imagine * c2.imagine)/d;
    temp.imagine = (c1.imagine * c2.real - c1.real * c2.imagine)/d;

    return temp;
}

cnum CnumDivD(cnum c1, double d)
{
    cnum temp;

    temp.real = c1.real / d;
    temp.imagine = c1.imagine / d;

    return temp;
}

cnum CnumDDiv(double x, cnum c1)
{
    double d;
    cnum temp;

    d = c1.real * c1.real + c1.imagine * c1.imagine;
    temp.real = x*c1.real/d;
    temp.imagine = -x*c1.imagine/d;

    return temp;
}

cnum CnumSqrt(cnum c1)
{
    double x,y;
    cnum temp;

    x = c1.real; y = c1.imagine;

    temp.real = (0.70710678118)*(sqrt(sqrt(x*x+y*y)+x));
    temp.imagine = sqrt(x*x+y*y)-x;
    if( fabs(temp.imagine) < 1E-3)
      temp.imagine = 0.0;
    temp.imagine = (0.70710678118)*(sqrt(temp.imagine));
    if(y < 0)
        temp.imagine = -temp.imagine;

    return temp;
}

cnum CnumCbrt(cnum c1, int verb)
{
    double r,theta;
    int index;
    double err,x,y;
    cnum roots[3],root;

    if (c1.real == 0.0 && c1.imagine == 0.0) {
        root.real = 0.0; root.imagine = 0.0;
        return root;
    }

    r = sqrt(c1.real * c1.real + c1.imagine * c1.imagine);
    r = pow(r,1.0/3.0);
    theta = atan(c1.imagine/c1.real);

    roots[0].real = r * cos(theta/3.0);
    roots[0].imagine = r * sin(theta/3.0);
    root = CnumMult(CnumMult(roots[0],roots[0]),roots[0]);
    x = root.real - c1.real; y = root.imagine - c1.imagine;
    err = sqrt(x*x + y*y); index = 0;

    roots[1].real = r * cos(theta/3.0 + 2.0*pi/3.0);
    roots[1].imagine = r * sin(theta/3.0 + 2.0*pi/3.0);
    root = CnumMult(CnumMult(roots[1],roots[1]),roots[1]);
    x = root.real - c1.real; y = root.imagine - c1.imagine;
    if( err > sqrt(x*x + y*y))
    {
        err = sqrt(x*x + y*y);
        index = 1;
    }

    roots[2].real = r * cos(theta/3.0 + 4.0*pi/3.0);
    roots[2].imagine = r * sin(theta/3.0 + 4.0 * pi/3.0);
    root = CnumMult(CnumMult(roots[2],roots[2]),roots[2]);
    x = root.real - c1.real; y = root.imagine - c1.imagine;
    if( err > sqrt(x*x + y*y))
    {
        err = sqrt(x*x + y*y);
        index = 2;
    }

    /* negative r */
    r = -1.0 * r;

    roots[0].real = r * cos(theta/3.0);
    roots[0].imagine = r * sin(theta/3.0);
    root = CnumMult(CnumMult(roots[0],roots[0]),roots[0]);
    x = root.real - c1.real; y = root.imagine - c1.imagine;
    if( err > sqrt(x*x + y*y))
    {
        err = sqrt(x*x + y*y);
        index = 3;
    }

    roots[1].real = r * cos(theta/3.0 + 2.0*pi/3.0);
    roots[1].imagine = r * sin(theta/3.0 + 2.0*pi/3.0);
    root = CnumMult(CnumMult(roots[1],roots[1]),roots[1]);
    x = root.real - c1.real; y = root.imagine - c1.imagine;
    if( err > sqrt(x*x + y*y))
    {
        err = sqrt(x*x + y*y);
        index = 4;
    }

    roots[2].real = r * cos(theta/3.0 + 4.0*pi/3.0);
    roots[2].imagine = r * sin(theta/3.0 + 4.0 * pi/3.0);
    root = CnumMult(CnumMult(roots[2],roots[2]),roots[2]);
    x = root.real - c1.real; y = root.imagine - c1.imagine;
    if( err > sqrt(x*x + y*y))
    {
        err = sqrt(x*x + y*y);
        index = 5;
    }

    if(index > 2)
        return roots[index-3];

    roots[index].real *= -1.0;
    roots[index].imagine *= -1.0;

    return roots[index];
}

/* EOF */
