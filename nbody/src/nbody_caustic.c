/*
Copyright (C) 2014 Julie Dumas
Copyright (C) 2014 Jake Weiss

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
/* WARNING THIS CODE DOES NOT USE MILKYWAY@HOME LIBRARIES. BEWARE WHEN RUNNING IT ON MULTIPLE SYSTEMS*/
#include "nbody_caustic.h"

const double G = 1.0;
const double a_n[] = {1.0,40.1,20.1,13.6,10.4,8.4,7.0,6.1,5.3,4.8,4.3,4.0,3.7,3.4,3.2,3.0,2.8,2.7,2.5,2.4,2.3};
const double V_n[] = {1.0,517,523,523,523,522,521,521,520,517,515,512,510,507,505,503,501,499,497,496,494};
const double rate_n[] = {1.0,53,23,14,10,7.8,6.3,5.3,4.5,3.9,3.4,3.1,2.8,2.5,2.3,2.1,2.0,1.8,1.7,1.6,1.5};
const double p_n[] = {1.0,0.3,0.3,1.0,0.3,0.15,0.12,0.6,0.23,0.41,0.25,0.19,0.17,0.11,0.09,0.09,0.09,0.09,0.09,0.09,0.09};

#define ONE_THIRD ((double)1.0/(double)3.0)
#define c1 (double complex) 1.0
#define c2 (double complex) 2.0
#define c3 (double complex) 3.0
#define c4 (double complex) 4.0
#define c12 (double complex) 12.0
#define c05 (double complex) 0.5
#define II (double complex) 1.0*I
#define TOLERANCE ((double)0.00001)
#define TOLERANCE2 ((double)0.001)

/* Calculate gfield far away for nth flow */
static void gfield_far(double rho, double z, int n, double *rfield, double *zfield)
{

    /* simulation units are kpc, gyr, ms=222288.47*Ms */

    double r_squared = rho*rho+z*z;

    double shift = a_n[n]+p_n[n]/4.0;

    double A_n = (8.0*M_PI*G*rate_n[n]*4498.6589)/(V_n[n]*1.0226831);  //convert from M_sun/yr*sr to m_sim/gyr*sr with 4498 and km/s to kpc/gyr with 1.0226

    /* caustic radius shifted by 0.25 so that g goes to zero beyond a_n[n] */
    double s = hypot(r_squared-(shift)*(shift),2.0*(shift)*z);

    double factor = (-A_n)/(s*(2.0*(shift)*(shift)+s));

    /* result in kpc/gyr^2 */

    *rfield += factor*(r_squared*rho-(shift)*(shift)*rho);
    *zfield += factor*(r_squared*z+(shift)*(shift)*z);
}


/* now let's compute the values of T at which the various poles cross the real axis */
static inline double complex f1(double complex x)
{
    return (2.0+4.0*x)/3.0;
}

static inline double complex f2(double complex x, double complex z)
{
    return 4.0*((1.0-x)*(1.0-x)-3.0*z*z);
} 

static inline double complex f3(double complex x, double complex z)
{
    return -128.0*cpow(x-1.0,3.0)-1728.0*z*z-1152.0*(x-1.0)*z*z;
}

/* for convenience define f4, which appears frequently in the roots */
static double complex f4(double complex x, double complex z)
{
    return cpow( (f3(x,z)+csqrt(-256.0*f2(x,z)*f2(x,z)*f2(x,z)+f3(x,z)*f3(x,z)))/2.0, ONE_THIRD );
}

static inline double complex T1(double complex x, double complex z)
{
  
    if(creal(x) >= -0.125 && cabs(f4(x,z))<TOLERANCE2)
    {
        double complex temp = cpow((-1.0-6.0*x+15.0*x*x-8.0*x*x*x),ONE_THIRD)/3.0;
        return c05*(c1 - csqrt(f1(x)/2.0 + temp) - csqrt(f1(x) - temp - c2*x/csqrt(f1(x)/2.0+temp) ));
    }
    else
    {
        double complex temp = csqrt(f1(x)/2.0 + f2(x,z)/(c3*f4(x,z)) + f4(x,z)/c12);
        return c05*(c1 - temp - csqrt(f1(x) - f2(x,z)/(c3*f4(x,z)) - f4(x,z)/c12 - c2*x/temp ));
    }
}

static inline double complex T2(double complex x, double complex z)
{
  
    if(creal(x) >= -0.125 && cabs(f4(x,z))<TOLERANCE2)
    {
        double complex temp = cpow((-1.0-6.0*x+15.0*x*x-8.0*x*x*x),ONE_THIRD)/3.0;
        return c05*(c1 - csqrt(f1(x)/2.0 + temp) + csqrt(f1(x) - temp - c2*x/sqrt(f1(x)/2.0+temp) ));
    }
    else
    {
        double complex temp = csqrt(f1(x)/2.0 + f2(x,z)/(c3*f4(x,z)) + f4(x,z)/c12);
        return c05*(c1 - temp + csqrt(f1(x) - f2(x,z)/(c3*f4(x,z)) - f4(x,z)/c12 - c2*x/temp ));
    }
}

static inline double complex T3(double complex x, double complex z)
{
  
    if(creal(x) >= -0.125 && cabs(f4(x,z))<TOLERANCE2)
    {
        double complex temp = cpow((-1.0-6.0*x+15.0*x*x-8.0*x*x*x),ONE_THIRD)/3.0;
        return c05*(c1 + csqrt(f1(x)/2.0 + temp) - csqrt(f1(x) - temp + c2*x/csqrt(f1(x)/2.0+temp) ));
    }
    else
    {
        double complex temp = csqrt(f1(x)/2.0 + f2(x,z)/(c3*f4(x,z)) + f4(x,z)/c12);
        return c05*(c1 + temp - csqrt(f1(x) - f2(x,z)/(c3*f4(x,z)) - f4(x,z)/c12 + c2*x/temp ));
    }
}

static inline double complex T4(double complex x, double complex z)
{
  
    if(creal(x) >= -0.125 && cabs(f4(x,z))<TOLERANCE2)
    {
        double complex temp = cpow((-1.0-6.0*x+15.0*x*x-8.0*x*x*x),ONE_THIRD)/3.0;
        return c05*(c1 + csqrt(f1(x)/2.0 + temp) + csqrt(f1(x) - temp + c2*x/csqrt(f1(x)/2.0+temp) ));
    }
    else
    {
        double complex temp = csqrt(f1(x)/2.0 + f2(x,z)/(c3*f4(x,z)) + f4(x,z)/c12);
        return c05*(c1 + temp + csqrt(f1(x) - f2(x,z)/(c3*f4(x,z)) - f4(x,z)/c12 + c2*x/temp ));
    }
}

/* f6 my suspicion is that the gravitational field in r and z constitute the real and complex components of the individual terms in f5; so let's do a test to see if this is right */
static inline double complex f5(double complex x, double complex z, double complex T)
{
  return (csqrt(c2*T - c1 + x - II*z))/c2;
}

static void gfield_close(double rho, double z, int n, double *rfield, double *zfield)
{

    z = z/p_n[n];
    z = (double complex)z;
    double complex x = ((double complex)rho-(double complex)a_n[n])/(double complex)p_n[n];

    //two cases; within the caustic, all four roots are real, so we use T1 for d4 and T2,3,4 for d3
    //outside the caustic there should be two real roots and two complex roots
    //the vector intlimits hold the limits of integration; for the first case, intlimits is of length four; second, length two

    double im[4] = {cabs(cimag(T1(x,z))), cabs(cimag(T2(x,z))), cabs(cimag(T3(x,z))), cabs(cimag(T4(x,z)))};
    double re[4] = {creal(T1(x,z)), creal(T2(x,z)), creal(T3(x,z)), creal(T4(x,z))};

    int count=0, i=0, j=0;
    double temp=0.0;

    //count the roots
    for (i = 0; i < 4; i++)
    {
        if(im[i]<TOLERANCE)
        {
            count++;
        }  
    }

    if(count==4)
    {
        /* sort in ascending order */
        for(i=0;i<4;i++)
        {
            for(j=i+1;j<4;j++)
            {
                if(re[i]>re[j])
                {
                    temp=re[i];
                    re[i]=re[j];
                    re[j]=temp;
                } 
            } 
        }
    }

    int c=0;
    double re2[2]={0.0,0.0};
    if(count==2)
    {
        for (i = 0; i < 4; i++)
        {
            if(im[i]<TOLERANCE && c<3)
            {
                re2[c]=re[i];    
                c++;
            } 
        }

    /* sort in ascending order */
        for(i=0;i<2;i++)
        {
            for(j=i+1;j<2;j++)
            {
                if(re2[i]>re2[j])
                {
                    temp=re2[i];
                    re2[i]=re2[j];
                    re2[j]=temp;
                } 
            } 
        }  
    }

    /* roots stored in re[] and re2[] */


    double factor = -8.0*M_PI*G*rate_n[n]*4498.6589/(rho*V_n[n]*1.0226831);

    if(count==2)
    {
        /*cout << "At (X,Z)=(" << X << "," << Z << "), you are outside the caustic." << endl; */

        *rfield += factor*(creal(f5(x,z,re2[1]))-creal(f5(x,z,re2[0]))-0.5); 

        *zfield += factor*(cimag(f5(x,z,re2[1]))-cimag(f5(x,z,re2[0])));

    }
    else
    {
        /*there are four elements, so we are in the caustic
        * it seems that the answer is just off by 0.5, so we subtract it by hand*/

        *rfield += factor*(creal(f5(x,z,re[3]))-creal(f5(x,z,re[2])) + creal(f5(x,z,re[1]))-creal(f5(x,z,re[0]))-0.5);


        *zfield += factor*(cimag(f5(x,z,re[3]))-cimag(f5(x,z,re[2])) + cimag(f5(x,z,re[1]))-cimag(f5(x,z,re[0])));
    }  
}
 
mwvector causticHaloAccel(const Halo* h, mwvector pos, real r)
{

    mwvector accel;

/* 20070507 bwillett used hypot from math.h */


    real rho=0.0, z=0.0, rfield, zfield;
    real R, l, tr, tl;
    int n;

    rfield = 0.0;
    zfield = 0.0;



    rho = mw_sqrt(sqr(X(pos))+sqr(Y(pos)));

    for (n = 1; n <= 20; n++)
    {

        /* tricusp */
        R=(3.0-mw_sqrt(1.0+(8.0/p_n[n])*(rho-a_n[n])))/4.0;
        l=(3.0+mw_sqrt(1.0+(8.0/p_n[n])*(rho-a_n[n])))/4.0;
        tr=2.0*p_n[n]*mw_sqrt(cube(R)*(1.0-R));
        tl=2.0*p_n[n]*mw_sqrt(cube(l)*(1.0-l));


        if( (Z(pos)<=tr && Z(pos)>=0.0 && rho>=a_n[n] && rho<=a_n[n]+p_n[n]) || (Z(pos)>=tl && Z(pos)<=tr && rho>=(a_n[n]-p_n[n]/8.0) && rho<=a_n[n]) || (Z(pos)>=-tr && Z(pos)<=0.0 && rho>=a_n[n] && rho<=a_n[n]+p_n[n]) || (Z(pos)<=-tl && Z(pos)>=-tr && rho>=(a_n[n]-p_n[n]/8.0) && rho<=a_n[n]) )  //close
        {
            gfield_close(rho,Z(pos),n,&rfield,&zfield);
        }

        else /* far */
        {
            gfield_far(rho,Z(pos),n,&rfield,&zfield);        
        }

    }

    if(rho<0.000001)
    {
        X(accel) = 0.0;
        Y(accel) = 0.0;
    }
    else
    {
        X(accel) = (rfield*X(pos))/rho;
        Y(accel) = (rfield*Y(pos))/rho;
    }

    Z(accel) = zfield;

    //printf("%f, %f, %f, %f, %f, %f\n", X(pos), Y(pos), Z(pos), X(accel), Y(accel), Z(accel));
    
    return accel;
}