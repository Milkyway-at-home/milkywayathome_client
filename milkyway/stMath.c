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

#include "milkyway_complex.h"
#include "stMath.h"


#define stSq(x) x * x



/* Evaluates a cubic: x^3 + ax^2 + bx + c */
/* returns y-value */
double stCEval(double complex a, double complex b, double complex c, double complex x)
{
    double complex y = 1.0 + 0.0 * I;

    y = x * y + a;
    y = x * y + b;
    y = x * y + c;

    return creal(y);
}


/* Evaluates the value of the quartic: x^4 + ax^3 + bx^2 + cx + d */
double stQEval(double a, double b, double c, double d, double complex x)
{
    double complex eval = 1.0 + 0.0 * I;

    eval *= x += a;
    eval *= x += b;
    eval *= x += c;
    eval *= x += d;

    return creal(eval);
}


/* stRoot3 - finds the real root of the cubic equation: */
/* x^3 + a2*x^2 + a1*x + a0 = 0 */
/* ans = max real root of the cubic */
/* returns 1 on success and -1 on error */
int stRoot3(double a2, double a1, double a0, double *ans, int verb)
{
    double complex roots[3];
    int numroots, i, index;
    const double lim = 0.001;

    numroots = DCubic(a2, a1, a0, roots, verb);

    if (numroots == -1)
        return -1;
    if (numroots == 0)
    {
        printf("Error: no real roots from cubic: (%g,%g,%g)\n", a2, a1, a0);
        ans[0] = -1;
        return 1;
    }

    index = -1;
    for (i = 0; i < 3; ++i)
    {
        if (fabs(cimag(roots[i])) < lim)
        {
            if (index == -1 || creal(roots[i]) > creal(roots[index]))
                index = i;
        }
    }

    if (index == -1)
        return -1;
    ans[0] = creal(roots[index]);
    return 1;
}


/* This is a wrapper function. It calls CnumCubic to find the roots */
/* of the cubic: x^3 + ax^2 + bx + c. Used when we only have real coeffs */
int DCubic(double a, double b, double c, double complex *roots, int verb)
{
    double complex A, B, C;

    A = a + 0.0 * I;
    B = b + 0.0 * I;
    C = c + 0.0 * I;

    return CnumCubic(A,B,C,roots,verb);
}


/* Finds roots of the cubic: x^3 + ax^2 + bx + c */
/* This version can handle complex coefficients */
int CnumCubic( double complex a, double complex b, double complex c, double complex *roots, int verb)
{
    double complex p, q, bq, cq;
    double complex temp;
    unsigned int i, index;
    const double lim = 0.000001;
    const double lim2 = 0.0001;

    if ( (cabs(a) < lim) && (cabs(b) < lim) )
    {
        /* if a and b are zero, then we have x^3 + c = 0 */
        /* so our first root is the cube root of -c. */
        roots[0] = ccbrt(-c);
    }
    else
    {
        /* p = (3.0*b - a*a)/3.0; */
        /* q = (27.0*c - 9.0*a*b + 2.0*a*a*a)/27.0; */
        p = b - (a*a / 3.0);
        q = (27.0*c - 9.0*a*b + 2.0*a*a*a)/27.0;

        if (verb)
        {
            printf("p = %g + i * %g\n",creal(p), cimag(p));
            printf("q = %g + i * %g\n",creal(q),cimag(q));
        }

        /* roots[0] = (-q + sqrt(q*q+ 4.0*p*p*p/27.0))/2.0; */
        /* roots[0] = std::pow(roots[0],1.0/3.0); */
        /* roots[0] = roots[0] - p/(3.0*roots[0]) - a/3.0; */

        temp = csqrt( q*q + (4.0 * p*p*p / 27.0))/2.0;

        roots[0] = (-q + temp) / 2.0;
        roots[1] = (-q - temp) / 2.0;

        if ( fabs(creal(roots[0])) < fabs(creal(roots[1])) )
            roots[0] = roots[1];

        if (verb)
            printf("z = %g + i * %g\n", creal(roots[0]), cimag(roots[0]));

        roots[0] = ccbrt(roots[0]);

        if (verb)
            printf("z^(1/3) = %g + i * %g\n",creal(roots[0]), cimag(roots[0]));

        //roots[0] -= (p / (3.0 * roots[0])) + a/3.0;
        roots[0] = roots[0] - p / (3.0 * roots[0]) - a/3.0;

        if (verb)
        {
            printf("Final Root: %g + i * %g\n", creal(roots[0]), cimag(roots[0]));
            printf("Evaluated: %g\n",stCEval(a,b,c,roots[0]));
        }

        if (stCEval(a,b,c,roots[0]) > lim2)
        {
            if (verb)
                printf("Recalculating root0...\n");

            roots[0] = (-q - csqrt( q*q + (4.0 * p*p*p) / 27.0)) / 2.0;

            if (verb)
                printf("z = %g + i * %g\n", creal(roots[0]), cimag(roots[0]));

            roots[0] = ccbrt(roots[0]);

            if (verb)
                printf("z^(1/3) = %g + i * %g\n",creal(roots[0]), cimag(roots[0]));

            roots[0] += - (p / (3.0 * roots[0])) - (a / 3.0);

            if (verb)
                printf("Final Root: %g + i * %g\n", creal(roots[0]),cimag(roots[0]));
        }
    }

    /* use a quadratic to find the other roots */
    /* bq = a+roots[0]; */
    /* cq = -c/roots[0]; */ /* or could be cq = b + roots[0]*bq */
    bq = a + roots[0];
    cq = b + bq * roots[0];

    /* roots[1] = (-bq + sqrt(bq*bq-4.0*cq))/2.0; */
    /* roots[2] = (-bq - sqrt(bq*bq-4.0*cq))/2.0; */

    temp = csqrt(bq*bq - 4.0 * cq);
    roots[1] = (-bq + temp) / 2.0;
    roots[2] = (-bq - temp) / 2.0;

    index = 0;

    for ( i = 0; i < 3; ++i)
    {
        if (cimag(roots[i]) < lim2)
            index++;

        if (verb)
            printf("Root#%d: %g + i * %g\n", i, creal(roots[i]), cimag(roots[i]));
    }

    return index;
}


/* wrapper function. Returns the same as quartic, except that the roots */
/* are returned as double data types instead of cnum structures */
int stRoot4(double a, double b, double c, double d, double r[], int flag[], int verb)
{
    double complex u[4];
    unsigned int numroots, i;
    const double lim = 0.1;

    /* Comment this out to allow debugging messages  */
       verb = 0;

    numroots = 0;
    quartic(a,b,c,d,u,flag,verb);

    for (i = 0; i < 4; ++i)
    {
        if ( fabs(cimag(u[i])) < lim)
        {
            ++numroots;
            flag[i] = 1;
        }
        else
            flag[i] = 0;

        r[i] = creal(u[i]);
    }

    return numroots;
}


/* finds the roots of the quartic: x^4 + ax^3 + bx^2 + cx + d */
/* u[] holds the roots and flag holds the flags for which are real(1) */
/* and which are complex(0). It returns the number of real roots found. */
int quartic(double a, double b, double c, double d, double complex u[], int flag[], int verb)
{
    const double lim = 0.0001;
    const double lim2 = 0.00001;
    int numroots, i, j, k;
    double complex A, B, C, Ap, Bp, Cp;
    double complex Csqrt;
    double complex roots[4];
    double complex z[3];
    double complex D, E, F, G;
    double complex Dp, Ep;
    double teval;
    double complex temp;

    if (fabs(d) < lim)
    {
        /* we know that one root is x = 0 */
        u[0] = COMPLEXZERO;

        /* the other roots are roots of the cubic: x^3 + ax^2 + bx + c */
        if ( DCubic(a,b,c,z,verb) == -1)
        {
            printf("The cubic: x^3 + %g*x^2 + %g*x + %g = 0\n",a,b,c);
            printf("has no real roots\n");
            return 0;
        }

        u[1] = z[0];
        u[2] = z[1];
        u[3] = z[2];

        numroots = 0;
        for (i = 0; i < 4; ++i)
        {
            if ( fabs(cimag(roots[i]) < lim) )
            {
                flag[i] = 1;
                ++numroots;
            }
            else
                flag[i] = 0;
        }

        return numroots;
    }

    /* find deflated coeffs */
    /* A = -3.0/8.0*a*a + b;
       B = a*a*a/8.0 - a*b/2.0 + c;
       C = -3.0*a*a*a*a/256.0 + a*a*b/16.0 - a*c/4.0 + d; */
    A = (-3.0/8.0*a*a + b) + 0.0 * I;
    B = (a*a*a/8.0 - a*b/2.0 + c) + 0.0 * I;
    C = (-3.0*a*a*a*a/256.0 + a*a*b/16.0 - a*c/4.0 + d) + 0.0 * I;
    Csqrt = csqrt(C);

    if (verb)
    {
        printf("A = %g + i * %g\n", creal(A), cimag(A));
        printf("B = %g + i * %g\n", creal(B), cimag(B));
        printf("C = %g + i * %g\n", creal(C), cimag(C));
    }

    /* find cubic coeffs */
    /* Ap = sqrt(C) * 3.0 - A/2.0; */
    /* Bp = 2.0*C - A * sqrt(C); */
    /* Cp = -B*B/8.0; */
    Ap = csqrt(C) * 2.0 - A/2.0;
    Bp = 2.0 * C - A * csqrt(C);
    Cp = (-B*B / 8.0) + 0.0 * I;

    /* find real root of cubic */
    if( CnumCubic(Ap,Bp,Cp,z,0) == -1)
    {
        printf("The cubic: x^3 + %g*x^2 + %g*x + %g = 0\n", creal(Ap), creal(Bp), creal(Cp));
        printf("has no real roots\n");
        return 0;
    }

    if (verb)
    {
        printf("Ap = %g + i * %g\n", creal(Ap), cimag(Ap));
        printf("Bp = %g + i * %g\n", creal(Bp), cimag(Bp));
        printf("Cp = %g + i * %g\n", creal(Cp), creal(Cp));
        printf("cubic real root: %g + i*%g\n", creal(z[0]), creal(z[0]));
    }

    teval = -1.0;
    for (i = 0; i < 3; ++i)
    {
        Dp = csqrt(2.0 * (z[i] + Csqrt) - A );
        Ep = csqrt( z[i] * (2.0 * Csqrt + z[i]));

        for (j = -1; j <= 1; j += 2)
        {
            for (k = -1; k <= 1; k += 2)
            {
                /* D = +/- sqrt(2.0*sqrt(C) - A + 2.0*z); */
                /* E = +/- sqrt(2.0*z * sqrt(C) + z*z) */
                /* temp = sqrt( (D*D) - 4.0*( sqrt(C) + z + E)) */
                /* F = (-D + temp)/2.0; */
                /* G = (-D - temp)/2.0; */

                D = Dp * (double)j;
                E = Ep * (double)k;

                temp = csqrt(D*D - 4.0 * (csqrt(C) + z[i] + E));
                F = (-D + temp) / 2.0;
                G = (-D - temp) / 2.0;

                /* adjustments */
                F -= a / 4.0;
                G -= a / 4.0;

                if(verb)
                {
                    printf("Cube root #%d   Sqrt D: %d  Sqrt E: %d\n",i,j,k);
                    printf("F: %g + i %g\n", creal(F), cimag(F));
                    printf("G: %g + i %g\n", creal(G), cimag(G));
                    printf("Eval(F): %g    Eval(G): %g\n", fabs(stQEval(a,b,c,d,F)),fabs(stQEval(a,b,c,d,G)));
                    printf("Eval: %g\n", (fabs(stQEval(a,b,c,d,F)) + fabs(stQEval(a,b,c,d,G)))/1.0 );
                }

                /* we need to exclude F and G if F == G */
                if ( F == G )
                {
                    if (verb) printf("F == G, skipped\n");
                    continue;
                }

                /* there shouldn't be a root of x = 0 */
                if ( fabs(creal(F)) < lim2 || fabs(creal(G)) < lim2 )
                {
                    if(verb) printf("x = 0, skipped\n");
                    continue;
                }

                /* if it's the first loop, or we find values that fit better to the roots */
                /* we take them */
                if ( (teval == -1.0) || ( ( (fabs(stQEval(a,b,c,d,F)) + fabs(stQEval(a,b,c,d,G)))/1.0 ) < teval) )
                {
                    /* average their evals and set roots */
                    teval = ( fabs(stQEval(a,b,c,d,F)) + fabs(stQEval(a,b,c,d,G)))/1.0;
                    roots[0] = F;
                    roots[1] = G;
                }
            }
        }
    }


    /* using the relation of the quartic coeffs to the roots,  */
    /* we can derive a quadratic equation to find the last two roots */
    if( (roots[0] == COMPLEXZERO) || (roots[1] == COMPLEXZERO) )
    {
        /* both roots are 0 */
        if (  (roots[0] == COMPLEXZERO) && (roots[1] == COMPLEXZERO) )
        {
            B = a + 0.0 * I;
            C = b + 0.0 * I;
        }
        else if ( roots[0] == COMPLEXZERO ) /* only root[0] is 0 */
        {
            B = roots[1] + a;
            C = -c / roots[1];
        }
        else /* only root[1] is 0 */
        {
            B = roots[0] + a;
            C = -c / roots[0];
        }
    }
    else
    {
        B = a + roots[0] + roots[1];
        C = d / (roots[0] * roots[1]);
    }

    temp = csqrt(B*B - 4.0 * C);
    roots[2] = (-B + temp) / 2.0;
    roots[3] = (-B - temp) / 2.0;

    /* lets see if -0.9999 or 0.9999 works better, use them */
    temp = 0.9999 + 0.0 * I;
    if ( fabs(stQEval(a,b,c,d,temp)) < fabs(stQEval(a,b,c,d,roots[2])))
        roots[2] = 0.9999 + 0.0 * I;

    temp = (-0.9999) + cimag(temp);
    if ( fabs(stQEval(a,b,c,d,temp)) < fabs(stQEval(a,b,c,d,roots[3])))
        roots[3] = 0.9999 + 0.0 * I;

    numroots = 1;
    for (i = 0; i < 4; ++i)
    {
        u[i] = roots[i];
        if (verb)
            printf("Root #%d: %g + i %g\n", i,creal(roots[i]), cimag(roots[i]));
    }

    return numroots;
}


/* Sum - returns the sum of all 'size' elements in array[] */
double sum(double array[], int size)
{
    int i;
    double total = 0.0;

    for (i = 0; i < size; ++i)
        total = total + array[i];

    return total;
}


/* Min - finds the minimum value in array[]. Only values that have
   a 1 in flag[] at the corresponding index are evalutated. array[]
   and flag[] have 'size' number of elements. flag[] has 'numgood'
   number of 1's. Returns the index to the min value */
int min(double array[], int flag[], int size, int numgood)
{
    int i, j, mini;
    double minimum;
    int *newindex;
    double *comp;

    newindex = (int*)malloc(sizeof(int)*numgood);
    comp = (double*)malloc(sizeof(double)*numgood);
    j = 0;

    for (i = 0; i < size; ++i)
    {
        if (flag[i] == 1)
        {
            comp[j] = array[i];
            newindex[j] = i;
            j++;
        }
    }
    if (j != numgood)
        printf("numgood != j\n");

    minimum = comp[0];
    mini = 0;
    for (i = 0; i<numgood; ++i)
    {
        if (comp[i]<minimum)
        {
            minimum = comp[i];
            mini = i;
        }
    }

    mini=newindex[mini];

    free(newindex);
    free(comp);

    return mini;
}


double fact(int d)
{
    if (d == 0)
        return 1.0;
    return (double)d*fact(d-1);
}


void MakePoints()
{
    double x;
    double y;

    FILE* ifile;

    ifile=fopen ("xy-points.txt", "w");
    for (x=-10; x<=10; x=x+0.5){
        y=-x*x+6*x-5;
        fprintf(ifile, "%g %g\n", x,y);
    }
    fclose (ifile);
}

