#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "stMath.h"


double stSq( double x )
{
    return x * x;
}


cnum complexPower(cnum base,double Power)
{
    double a,theta;
    cnum temp;

    a = sqrt( stSq(base.real) + stSq(base.imagine) );

    if( a )
    {
        if(base.real == 0.0 && base.imagine < 0.0)
            theta = -2.0*atan(1.0);
        else if(base.real == 0.0 && base.imagine > 0.0)
            theta = 2.0*atan(1.0);
        else
            theta = atan(base.imagine/base.real);

        if(base.imagine <= 0.0 && base.real <= 0.0)
            theta =-4*atan(1.0)+theta;
        else if(base.imagine >= 0.0 && base.real <= 0.0)
            theta =4*atan(1.0)+theta;
        
        temp.real = pow(a,Power)*cos(theta*Power);
        temp.imagine = pow(a,Power)*sin(theta*Power);
        return temp;
    }
    else
    {
        temp.real = 0; temp.imagine = 0;
        return temp;
    }
}


/* Evaluates a cubic: x^3 + ax^2 + bx + c */
/* returns y-value */
double stCEval(cnum a, cnum b, cnum c, cnum x)
{
    cnum y;

    y.real = 1.0; y.imagine = 0.0;
    y = CnumAdd(CnumMult(y,x),a);
    y = CnumAdd(CnumMult(y,x),b);
    y = CnumAdd(CnumMult(y,x),c);
    	
    return y.real;
}


/* Evaluates the value of the quartic: x^4 + ax^3 + bx^2 + cx + d */
double stQEval(double a, double b, double c, double d, cnum x)
{
    cnum eval;

    eval.real = 1.0; eval.imagine = 0.0;
    eval = CnumMult(eval,x); eval.real += a;
    eval = CnumMult(eval,x); eval.real += b;
    eval = CnumMult(eval,x); eval.real += c;
    eval = CnumMult(eval,x); eval.real += d;

    return eval.real;
}


/* stRoot3 - finds the real root of the cubic equation: */
/* x^3 + a2*x^2 + a1*x + a0 = 0 */
/* ans = max real root of the cubic */
/* returns 1 on success and -1 on error */
int stRoot3(double a2, double a1, double a0, double *ans, int verb)
{
    cnum roots[3];
    int numroots,i,index;

    numroots = DCubic(a2,a1,a0,roots,verb);
	
    if(numroots == -1)
        return -1;
    if(numroots == 0)
    {
        printf("Error: no real roots from cubic: (%g,%g,%g)\n",a2,a1,a0);
        ans[0] = -1;
        return 1;
    }

    index = -1;
    for(i = 0; i < 3; i++)
    {
        if(fabs(roots[i].imagine) < 0.001)
        {
            /* printf("Root#%d is a real root\n",i); */
            if(index == -1 || roots[i].real > roots[index].real)
                index = i;
        }
    }

    if(index == -1)
        return -1;

    /* printf("Chose Root #%d\n",index); */
	
    ans[0] = roots[index].real;
    return 1;
}


/* This is a wrapper function. It calls CnumCubic to find the roots */
/* of the cubic: x^3 + ax^2 + bx + c. Used when we only have real coeffs */
int DCubic(double a, double b, double c, cnum *roots, int verb)
{
    cnum A,B,C;

    A.real = a; A.imagine = 0.0;
    B.real = b; B.imagine = 0.0;
    C.real = c; C.imagine = 0.0;

    return CnumCubic(A,B,C,roots,verb);
}


/* Finds roots of the cubic: x^3 + ax^2 + bx + c */
/* This version can handle complex coefficients */
int CnumCubic( cnum a, cnum b, cnum c, cnum *roots, int verb)
{
    cnum p,q,bq,cq;
    cnum temp;
    int i,index;

    if( (fabs(a.real) < 0.000001) && (fabs(a.imagine) < 0.000001) && (fabs(b.real) < 0.000001) && (fabs(b.imagine) < 0.000001) )
    {
        /* if a and b are zero, then we have x^3 + c = 0 */
        /* so our first root is the cube root of -c. */
        roots[0] = CnumCbrt( CnumMultD(c,-1.0),verb);
    }
    else
    {
        /* p = (3.0*b - a*a)/3.0; */
        /* q = (27.0*c - 9.0*a*b + 2.0*a*a*a)/27.0; */
        p = CnumDivD(CnumSub(CnumMultD(b,3.0),CnumMult(a,a)),3.0);
        q = CnumDivD(CnumAdd(CnumSub(CnumMultD(c,27.0),
                                     CnumMult(a,CnumMultD(b,9.0))), 
                             CnumMultD(CnumMult(CnumMult(a,a),a),2.0)),27.0);

        if(verb)
        {
            printf("p = %g + i * %g\n",p.real,p.imagine);
            printf("q = %g + i * %g\n",q.real,q.imagine);
        }

        /* roots[0] = (-q + sqrt(q*q+ 4.0*p*p*p/27.0))/2.0; */
        /* roots[0] = std::pow(roots[0],1.0/3.0); */
        /* roots[0] = roots[0] - p/(3.0*roots[0]) - a/3.0; */
        temp = CnumSqrt( CnumAdd(CnumMult(q,q),CnumDivD( CnumMultD(CnumMult(p,CnumMult(p,p)),4.0) ,27.0)));
        roots[0] = CnumDivD(CnumAdd(CnumMultD(q,-1.0), temp),2.0);
        roots[1] = CnumDivD(CnumSub(CnumMultD(q,-1.0), temp),2.0);
		
        if( fabs(roots[0].real) < fabs(roots[1].real))
        {
            roots[0].real    = roots[1].real;
            roots[0].imagine = roots[1].imagine;
        }

        if(verb)
            printf("z = %g + i * %g\n",roots[0].real, roots[0].imagine);

        roots[0] = CnumCbrt(roots[0],verb);	

        if(verb)
            printf("z^(1/3) = %g + i * %g\n",roots[0].real, roots[0].imagine);

        roots[0] = CnumSub(CnumSub(roots[0],CnumDiv(p, CnumMultD(roots[0],3.0))),CnumDivD(a,3.0));

        if(verb)
        {
            printf("Final Root: %g + i * %g\n",roots[0].real,roots[0].imagine);		
            printf("Evaluated: %g\n",stCEval(a,b,c,roots[0]));
        }

        if( stCEval(a,b,c,roots[0]) > 0.0001)
        {
            if(verb)
                printf("Recalculating root0...\n");

            roots[0] = CnumDivD(CnumSub(CnumMultD(q,-1.0), CnumSqrt( CnumAdd(CnumMult(q,q),
                                                                             CnumDivD( CnumMultD(CnumMult(p,CnumMult(p,p)),4.0) ,27.0)))),2.0);
            if(verb)
                printf("z = %g + i * %g\n",roots[0].real, roots[0].imagine);

            roots[0] = CnumCbrt(roots[0],verb);	

            if(verb)
                printf("z^(1/3) = %g + i * %g\n",roots[0].real, roots[0].imagine);

            roots[0] = CnumSub(CnumSub(roots[0],CnumDiv(p, CnumMultD(roots[0],3.0))),CnumDivD(a,3.0));

            if(verb)
                printf("Final Root: %g + i * %g\n",roots[0].real,roots[0].imagine);
        }
    }

    /* use a quadratic to find the other roots */
    /* bq = a+roots[0]; */
    /* cq = -c/roots[0]; */ /* or could be cq = b + roots[0]*bq */
    bq = CnumAdd(a,roots[0]);
    cq = CnumAdd(b,CnumMult(roots[0],bq));
	
    /* roots[1] = (-bq + sqrt(bq*bq-4.0*cq))/2.0; */
    /* roots[2] = (-bq - sqrt(bq*bq-4.0*cq))/2.0; */
    temp = CnumSqrt(CnumSub(CnumMult(bq,bq),CnumMultD(cq,4.0)));
    roots[1] = CnumDivD(CnumAdd(CnumMultD(bq,-1.0), temp),2.0);
    roots[2] = CnumDivD(CnumSub(CnumMultD(bq,-1.0), temp),2.0);
	
    index = 0;
	
    for( i = 0; i < 3; i++)
    {
        if(roots[i].imagine < 0.0001)
            index++;

        if(verb)
            printf("Root#%d: %g + i * %g\n",i,roots[i].real,roots[i].imagine);		
    }
	
    return index;
}


/* wrapper function. Returns the same as quartic, except that the roots */
/* are returned as double data types instead of cnum structures */
int stRoot4(double a, double b, double c, double d, double r[], int flag[], int verb)
{
    cnum u[4];
    int numroots,i;
    
    /* Comment this out to allow debugging messages  */
       verb = 0;
    
    numroots = 0; 
    quartic(a,b,c,d,u,flag,verb);

    for(i = 0; i < 4; i++)
    {
        if( fabs(u[i].imagine) < 0.1)
        {
            numroots++;
            flag[i] = 1;
        }
        else
            flag[i] = 0;

        r[i] = u[i].real;
    }

    return numroots;
}


/* finds the roots of the quartic: x^4 + ax^3 + bx^2 + cx + d */
/* u[] holds the roots and flag holds the flags for which are real(1) */
/* and which are complex(0). It returns the number of real roots found. */
int quartic(double a, double b, double c, double d, cnum u[], int flag[], int verb)
{
    int numroots,i,j,k;
    cnum A,B,C,Ap,Bp,Cp;
    cnum Csqrt;
    cnum roots[4];
    cnum z[3];
    cnum D,E,F,G;
    cnum Dp,Ep;
    double teval;
    cnum temp;

    if(fabs(d) < 0.0001)
    {
        /* we know that one root is x = 0 */
        u[0].real = 0.0; u[0].imagine = 0.0;

        /* the other roots are roots of the cubic: x^3 + ax^2 + bx + c */
        if( DCubic(a,b,c,z,verb) == -1)
        {
            printf("The cubic: x^3 + %g*x^2 + %g*x + %g = 0\n",a,b,c);
            printf("has no real roots\n");
            return 0;
        }
        
        u[1] = z[0];
        u[2] = z[1];
        u[3] = z[2];
        
        numroots = 0;
        for(i = 0; i < 4; i++)
        {
            if( fabs(roots[i].imagine) < 0.0001)
            {
                flag[i] = 1;
                numroots++;
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
    A.real = -3.0/8.0*a*a + b;                               A.imagine = 0.0;
    B.real = a*a*a/8.0 - a*b/2.0 + c;                        B.imagine = 0.0;
    C.real = -3.0*a*a*a*a/256.0 + a*a*b/16.0 - a*c/4.0 + d;  C.imagine = 0.0;
    Csqrt = CnumSqrt(C);
    
    if(verb)
    {
        printf("A = %g + i * %g\n",A.real,A.imagine);
        printf("B = %g + i * %g\n",B.real,B.imagine);
               printf("C = %g + i * %g\n",C.real,C.imagine);
    }
    
    /* find cubic coeffs */
    /* Ap = sqrt(C) * 3.0 - A/2.0; */
    /* Bp = 2.0*C - A * sqrt(C); */
    /* Cp = -B*B/8.0; */
    Ap = CnumSub(CnumMultD(Csqrt,3.0), CnumDivD(A,2.0));
    Bp = CnumSub(CnumMultD(C,2.0),CnumMult(A,Csqrt));
    /* Cp = CnumDivD(CnumMult(CnumMultD(B,-1.0),B),8.0); */
    Cp.real = -(B.real*B.real)/8.0; Cp.imagine = 0.0;
    
    /* find real root of cubic */
    if( CnumCubic(Ap,Bp,Cp,z,0) == -1)
    {
        printf("The cubic: x^3 + %g*x^2 + %g*x + %g = 0\n",Ap.real,Bp.real,Cp.real);
        printf("has no real roots\n");
        return 0;
    }
    
    if(verb)
    {
        printf("Ap = %g + i * %g\n",Ap.real,Ap.imagine);
        printf("Bp = %g + i * %g\n",Bp.real,Bp.imagine);
        printf("Cp = %g + i * %g\n",Cp.real,Cp.imagine);
        printf("cubic real root: %g + i*%g\n",z[0].real,z[0].imagine);
    }
    
    teval = -1.0;
    for(i = 0; i < 3; i++)
    {
        Dp = CnumSqrt( CnumAdd(CnumSub(CnumMultD(Csqrt,2.0),A),CnumMultD(z[i],2.0)));
        Ep = CnumSqrt( CnumAdd(CnumMultD(CnumMult(z[i],Csqrt),2.0), CnumMult(z[i],z[i])));
        
        for(j = -1; j <= 1; j += 2)
        {
            for(k = -1; k <= 1; k += 2)
            {
                /* D = +/- sqrt(2.0*sqrt(C) - A + 2.0*z); */
                /* E = +/- sqrt(2.0*z * sqrt(C) + z*z) */
                /* temp = sqrt( (D*D) - 4.0*( sqrt(C) + z + E)) */
                /* F = (-D + temp)/2.0; */
                /* G = (-D - temp)/2.0; */
                D = CnumMultD(Dp,(double)j);
                E = CnumMultD(Ep,(double)k);
                
                temp = CnumSqrt(CnumSub(CnumMult(D,D), CnumMultD(CnumAdd(Csqrt, CnumAdd(z[i],E)),4.0)));				
                F = CnumDivD(CnumAdd(CnumMultD(D,-1.0),temp),2.0);
                G = CnumDivD(CnumSub(CnumMultD(D,-1.0),temp),2.0);
                
                /* adjustments */
                F = CnumSubD(F,a/4.0);
                G = CnumSubD(G,a/4.0);
                
                if(verb)
                {
                    printf("Cube root #%d   Sqrt D: %d  Sqrt E: %d\n",i,j,k);
                    printf("F: %g + i %g\n",F.real,F.imagine);
                    printf("G: %g + i %g\n",G.real,G.imagine);
                    printf("Eval(F): %g    Eval(G): %g\n", fabs(stQEval(a,b,c,d,F)),fabs(stQEval(a,b,c,d,G)));
                    printf("Eval: %g\n", (fabs(stQEval(a,b,c,d,F)) + fabs(stQEval(a,b,c,d,G)))/1.0 );
                }
                
                /* we need to exclude F and G if F == G */
                if( (F.real == G.real) && (F.imagine == G.imagine) )
                {
                    if(verb) printf("F == G, skipped\n");
                    continue;
                }
                
                /* there shouldn't be a root of x = 0 */
                if( fabs(F.real) < 0.00001 || fabs(G.real) < 0.00001 )
                {
                    if(verb) printf("x = 0, skipped\n");
                    continue;
                }

                /* if it's the first loop, or we find values that fit better to the roots */
                /* we take them */
                if( (teval == -1.0) || ( ( (fabs(stQEval(a,b,c,d,F)) + fabs(stQEval(a,b,c,d,G)))/1.0 ) < teval) )
                {
                    /* average their evals and set roots */
                    teval = ( fabs(stQEval(a,b,c,d,F)) + fabs(stQEval(a,b,c,d,G)))/1.0; 
                    roots[0] = F;  roots[1] = G;
                }				
            }
        }
    }	
    
    
    /* using the relation of the quartic coeffs to the roots,  */
    /* we can derive a quadratic equation to find the last two roots */
    if( (roots[0].real == 0.0 && roots[0].imagine == 0.0) || (roots[1].real == 0.0 && roots[1].imagine == 0.0) )
    {
        /* both roots are 0 */
        if( (roots[0].real == 0.0 && roots[0].imagine == 0.0) && (roots[1].real == 0.0 && roots[1].imagine == 0.0) )
        {
            /* B = complex<double>(a,0.0); */
            /* C = complex<double>(b,0.0); */
            B.real = a; B.imagine = 0.0;
            C.real = b; C.imagine = 0.0;			
        }
        else if( roots[0].real == 0.0 && roots[0].imagine == 0.0 ) /* only root[0] is 0 */
        {
            /* B = (roots[1] + a); */
            /* C = -c/roots[1]; */
            B = CnumAddD(roots[1],a);
            C = CnumDDiv(-c,roots[1]);
        }
        else /* only root[1] is 0 */
        {
            /* B = (roots[0] + a); */
            /* C = -c/roots[0]; */
            B = CnumAddD(roots[0],a);
            C = CnumDDiv(-c,roots[0]);
        }
    }
    else
    {
        /* B = a+roots[0]+roots[1]; */
        /* C = d/(roots[0]*roots[1]); */
        B = CnumAddD(CnumAdd(roots[0],roots[1]),a);
        C = CnumDDiv(d,CnumMult(roots[0],roots[1]));
    }
    
    /* temp = sqrt(B*B-4.0*C) */
    /* roots[2] = (-B + temp)/2.0; */
    /* roots[3] = (-B - temp)/2.0; */
    temp = CnumSqrt(CnumSub( CnumMult(B,B), CnumMultD(C,4.0)));
    roots[2] = CnumDivD(CnumAdd(CnumMultD(B,-1.0),temp),2.0);
    roots[3] = CnumDivD(CnumSub(CnumMultD(B,-1.0),temp),2.0);
    
    /* lets see if -0.9999 or 0.9999 works better, use them */
    temp.real = 0.9999; temp.imagine = 0.0;
    if( fabs(stQEval(a,b,c,d,temp)) < fabs(stQEval(a,b,c,d,roots[2])))
    {
        roots[2].real = 0.9999;
        roots[2].imagine = 0.0;
    }
    
    temp.real = -0.9999;
    if( fabs(stQEval(a,b,c,d,temp)) < fabs(stQEval(a,b,c,d,roots[3])))
    {
        roots[3].real = 0.9999;
        roots[3].imagine = 0.0;
    }
    
    numroots = 1;
    for(i = 0; i < 4; i++)
    {
        u[i] = roots[i];
        if(verb)
        {
            printf("Root #%d: %g + i %g\n",i,roots[i].real,roots[i].imagine);
        }
    }
    
    return numroots;	
}


/* Sum - returns the sum of all 'size' elements in array[] */
double sum(double array[], int size)
{
    int i;
    double total;
    
    total=0.0;

    for(i=0; i<size; i++)
        total=total+array[i];
    
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
    j=0;
    
    for(i=0; i<size; i++)
    {
        if(flag[i]==1)
        {
            comp[j]=array[i];
            newindex[j]=i;
            j++;
        }
    }
    if (j != numgood)
        printf("numgood != j\n");
    
    minimum=comp[0];
    mini=0;
    for (i=0; i<numgood; i++)
    {
        if (comp[i]<minimum)
        {
            minimum=comp[i];
            mini=i;
        }
    }
    
    mini=newindex[mini];

    free(newindex);
    free(comp);
    
    return mini;
}


double fact(int d)
{
    if(d==0)
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
