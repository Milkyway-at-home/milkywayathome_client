/* Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
     Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.    If not, see <http://www.gnu.org/licenses/>.

The minimum bracketing method is based on results from "Numerical Recipes,
3rd ed." and conforms to the authors' defintion of intellectual
property under their license, and as such does not conflict with
their copyright to their programs which execute similar algorithms.
*/

#include "nbody_priv.h"
#include "milkyway_util.h"
#include "milkyway_math.h"
#include "milkyway_lua.h"
#include "nbody_lua_types.h"
#include "nbody_isotropic.h"

/*Note: minusfivehalves(x) raises to x^-5/2 power and minushalf(x) is x^-1/2*/


/*      MODEL SPECIFIC FUNCTIONS       */
static inline real potential( real r, real * args, dsfmt_t* dsfmtState)
{
    /*Be Careful! this function returns the negative of the potential! this is the value of interest, psi*/
    real mass_l   = args[0];
    real mass_d   = args[1];
    real rscale_l = args[2];
    real rscale_d = args[3];

    real potential_result = -(mass_l/mw_sqrt(sqr(r) + sqr(rscale_l)) + mass_d/mw_sqrt(sqr(r) + sqr(rscale_d)) );

    return (-potential_result);
}

static inline real density( real r, real * args, dsfmt_t* dsfmtState)
{
    /*this is the density distribution function. Returns the density at a given radius.*/
    real mass_l   = args[0];
    real mass_d   = args[1];
    real rscale_l = args[2];
    real rscale_d = args[3];

    real rscale_lCube = cube(rscale_l); 
    real rscale_dCube = cube(rscale_d);
    real density_light = (mass_l/rscale_lCube) * (minusfivehalves( (1.0 + sqr(r)/sqr(rscale_l)) ) );
    real density_dark  = (mass_d/rscale_dCube) * (minusfivehalves( (1.0 + sqr(r)/sqr(rscale_d)) ) ); 
    real density_result = (3.0/(4.0 * (M_PI))) * ( density_light + density_dark );

    return density_result;
}

static inline real mass_en( real r, real mass, real scaleRad)
{
    /*BE CAREFUL! this function returns the mass enclosed in a single plummer sphere!*/
    real mass_enclosed = mass * cube(r) * minusthreehalves( ( sqr(r) + sqr(scaleRad) ) ) ;

    return mass_enclosed;
}

static inline real profile_rho(real r, real * args, dsfmt_t* dsfmtState)
{
    real result = r * r * density(r, args, dsfmtState);    
    return result;
}



/*      GENERAL PURPOSE DERIVATIVE, INTEGRATION, MAX FINDING, ROOT FINDING, AND ARRAY SHUFFLER FUNCTIONS        */
static inline real first_derivative(real (*func)(real, real *, dsfmt_t*), real x, real * funcargs, dsfmt_t* dsfmtState)
{
    /*yes, this does in fact use a 5-point stencil*/
    real h = 0.1;
    real deriv;
    real p1, p2, p3, p4, denom;
    
    p1 =   1.0 * (*func)( (x - 2.0 * h), funcargs, dsfmtState);
    p2 = - 8.0 * (*func)( (x - h)      , funcargs, dsfmtState);
    p3 = - 1.0 * (*func)( (x + 2.0 * h), funcargs, dsfmtState);
    p4 =   8.0 * (*func)( (x + h)      , funcargs, dsfmtState);
    denom = inv( 12.0 * h);
    deriv =   (p1 + p2 + p3 + p4) * denom;
    return deriv;
}

static inline real second_derivative(real (*func)(real, real *, dsfmt_t*), real x, real * funcargs, dsfmt_t* dsfmtState)
{
    /*yes, this also uses a five point stencil*/
    real h = 0.1;
    real deriv;
    real p1, p2, p3, p4, p5, denom;

    p1 = - 1.0 * (*func)( (x + 2.0 * h) , funcargs, dsfmtState);
    p2 =  16.0 * (*func)( (x + h)       , funcargs, dsfmtState);
    p3 = -30.0 * (*func)( (x)           , funcargs, dsfmtState);
    p4 =  16.0 * (*func)( (x - h)       , funcargs, dsfmtState);
    p5 = - 1.0 * (*func)( (x - 2.0 * h) , funcargs, dsfmtState);
    denom = inv( 12.0 * h * h);
    deriv =   (p1 + p2 + p3 + p4 + p5) * denom;
    return deriv;
}

static inline real gauss_quad(real (*func)(real, real *, dsfmt_t*), real lower, real upper, real benchmark, real * funcargs, dsfmt_t* dsfmtState)
{
    /*This is a guassian quadrature routine. It will test to always integrate from the lower to higher of the two limits.
     * If switching the order of the limits was needed to do this then the negative of the integral is returned.
     */
    real Ng,hg,lowerg, upperg;
    real intv = 0;//initial value of integral
    real coef1,coef2;//parameters for gaussian quad
    real c1,c2,c3;
    real x1,x2,x3;
    real x1n,x2n,x3n;
    real a, b;
    
    if(lower > upper)
    {
        a = upper;
        b = lower;
    }
    else
    {
        a = lower; 
        b = upper;
    }
    
    benchmark = 1.2 * a;
    Ng = 100.0;//integral resolution
    hg = (benchmark - a) / (Ng);
    lowerg = a;
    upperg = lowerg + hg;
    

    coef2 = (lowerg + upperg) / 2.0;//initializes the first coeff to change the function limits
    coef1 = (upperg - lowerg) / 2.0;//initializes the second coeff to change the function limits
    c1 = 5.0 / 9.0;
    c2 = 8.0 / 9.0;
    c3 = 5.0 / 9.0;
    x1 = -sqrt(3.0 / 5.0);
    x2 = 0.0;
    x3 = sqrt(3.0 / 5.0);
    x1n = (coef1 * x1 + coef2);
    x2n = (coef1 * x2 + coef2);
    x3n = (coef1 * x3 + coef2);
    int counter = 0;
    while (1)
    {
                //gauss quad
//         mw_printf("\nx1 = %f  x2 = %f x3 = %f low = %f  upp = %f upper = %f lower = %f\n", x1n, x2n, x3n, lowerg, upperg, upper, lower);
        intv = intv + c1 * (*func)(x1n, funcargs, dsfmtState) * coef1 +
                      c2 * (*func)(x2n, funcargs, dsfmtState) * coef1 + 
                      c3 * (*func)(x3n, funcargs, dsfmtState) * coef1;

        lowerg = upperg;
        upperg = upperg + hg;
        coef2 = (lowerg + upperg)/2.0;//initializes the first coeff to change the function limits
        coef1 = (upperg - lowerg)/2.0;

        x1n = ((coef1) * x1 + coef2);
        x2n = ((coef1) * x2 + coef2);
        x3n = ((coef1) * x3 + coef2);

        if(lowerg > benchmark)
        {
            Ng = 20.0;//integral resolution
            hg = (b - benchmark) / (Ng);
        }
            
        
        
        
        if(upper > lower)
        {
            if(lowerg >= upper)//loop termination clause
            {
                break;
            }
        }
        else if(lower > upper)
        {
            if(lowerg >= lower)//loop termination clause
            {
                break;
            }
        }
        
        if(counter > 100000)
        {
            break;
        }
        else
        {
            counter++;
        }
        
        
    }
    
    if(lower > upper)
    {
        a = upper;
        b = lower;
        intv *= -1.0;
    }
    else
    {
        a = lower; 
        b = upper;
    }
    
//     mw_printf("int = %0.20f\n", intv);
    return intv;
}

static real simpson( real (*func)(real, real *, dsfmt_t*), real lower_limit, real upper_limit, real * args, dsfmt_t* dsfmtState)
{
    real a,b,N,h;
    real intv = 0.0;
    real lower, upper, mid, coef;
    a = lower_limit;//lower limit 
    b = upper_limit;//upper limit
    N = 10;//number of intervals
    h = (b-a)/N;//step width
    lower = a;//initializes the temp lower limit
    upper = a + h;//initializes the temp upper limit
    mid = (lower + upper)/2;//initializes the temp limit midpoint
    coef = (upper - lower)/6;//initializes the coeff in the Simpson's rule formula
    while (1)
    {
    //simpsons rule
        intv = intv + coef * (   (*func)(lower, args, dsfmtState) + 
                             4 * (*func)(mid, args, dsfmtState) +
                                 (*func)(upper, args, dsfmtState) );//gets the integral value for the interval
        lower = upper;//iterates the lower limit
        upper = upper + h;//iterates the upper limit
        mid = (lower + upper)/2;//gets midpoint of new interval
        coef = (upper - lower)/6;//gets the coeff for the new interval

        if(b > a)
        {
            if(lower >= b)//loop termination clause
            {
                break;
            }
        }
        else if(a > b)
        {
            if(lower <= b)//loop termination clause
            {
                break;
            }
        }

    }
    
    return intv;
    
}

static inline real max_finder(real (*profile)(real , real*, dsfmt_t*), real* profileParams, real a, real b, real c, int limit, real tolerance, dsfmt_t* dsfmtState)
{
    /*this is a maxfinding routine to find the maximum of the density.
     * It uses Golden Section Search as outlined in Numerical Recipes 3rd edition
     */
    real RATIO = 0.61803399;
    real RATIO_COMPLEMENT = 1 - RATIO;
    int counter = 0;
    
    real profile_x1,profile_x2,x0,x1,x2,x3;
    x0 = a;
    x3 = c;
    
    if (mw_fabs(b - c) > mw_fabs(b - a))
    {
        x1 = b;
        x2 = b + (RATIO_COMPLEMENT * (c - b)); 
    }
    else
    {
        x2 = b;
        x1 = b - (RATIO_COMPLEMENT * (b - a));
    }

    profile_x1 = -(*profile)(x1, profileParams, dsfmtState);
    profile_x2 = -(*profile)(x2, profileParams, dsfmtState);
    
    while (mw_fabs(x3 - x0) > (tolerance * (mw_fabs(x1) + mw_fabs(x2)) ) )
    {
        counter++;
        if (profile_x2 < profile_x1)
        {
            x0 = x1;
            x1 = x2;
            x2 = RATIO * x2 + RATIO_COMPLEMENT * x3;
            profile_x1 = (real)profile_x2;
            profile_x2 = -(*profile)(x2,profileParams, dsfmtState);
        }
        else
        {
            x3 = x2;
            x2 = x1;
            x1 = RATIO * x1 + RATIO_COMPLEMENT * x0;
            profile_x2 = (real)profile_x1;
            profile_x1 = -(*profile)(x1, profileParams, dsfmtState);
        }
        
        if(counter > limit)
        {
            break;
        }
    }

    if (profile_x1 < profile_x2)
    {
        return (-profile_x1);
    }
    else
    {
        return (-profile_x2);
    }
}

static inline real root_finder(real (*rootFunc)(real, real*, dsfmt_t*), real* rootFuncParams, real funcValue, real lowBound, real upperBound, dsfmt_t* dsfmtState)
{
    //requires lowBound and upperBound to evaluate to opposite sign when rootFunc-funcValue
    if(rootFuncParams == NULL || rootFunc == NULL)
    {
        exit(-1);
    }
    unsigned int i = 0;
    /* Can find up to 20 roots, but may miss a root if they are too close together */
    int N = 4;
    unsigned int numSteps = N;
    real interval;
    real * values = mwCalloc(numSteps, sizeof(real));
    
    /* Divide the function area up into bins in case there is more than one root */
    /*numSteps+1 because you want to include the upperbound in the interval*/
    for(i = 0; i < numSteps + 1; i++)
    {
        interval = ((upperBound - lowBound) * (real)i) / (real)numSteps + lowBound;
        values[i] = (*rootFunc)(interval, rootFuncParams, dsfmtState) - funcValue;
    }
    
    real midPoint = 0;
    real midVal = 0;
    unsigned int nsteps = 0;
    real curUpper = 0;
    real curLower = 0;
    int rootsFound = 0;
    int q = 0;
    
    /* Find the roots using bisection because it was easy to code and good enough for our purposes 
     * this will hope around the different intervals until it checks all of them. This way it does not 
     * favor any root.
     */
    for(i = 0; i < numSteps; i++)
    {
        q = i;
        if((values[q] > 0 && values[q + 1] < 0) || (values[q] < 0 && values[q + 1] > 0))
        {
            if(values[q] < 0 && values[q + 1] > 0)
            {
                curLower = ((upperBound - lowBound) * (real)q)/(real)numSteps + lowBound;
                curUpper = ((upperBound - lowBound) * (real)(q + 1))/(real)numSteps + lowBound;
            }
            else if(values[q] > 0 && values[q+1] < 0)
            {
                curLower = ((upperBound - lowBound) * (real)(q + 1))/(real)numSteps + lowBound;
                curUpper = ((upperBound - lowBound) * (real)q)/(real)numSteps + lowBound;
            }
            else
            {
                continue;
            }
            midVal = 1;
            nsteps = 0;
            while(mw_fabs(midVal) > .0001 || nsteps >= 10000)
            {
                midPoint = (curLower + curUpper)/2.0;
                midVal = (*rootFunc)(midPoint, rootFuncParams, dsfmtState) - funcValue;
                
                if(midVal < 0.0)
                {
                    curLower = midPoint;
                }
                else
                {
                    curUpper = midPoint;
                }
                ++nsteps;
            }
            
            if(nsteps < 10000)
            {
                ++rootsFound;
            }
            else
            {
                return midPoint = 0.0;
            }
            
        }
        
        if(rootsFound != 0)
        {
            break;
        }
    }
//     mw_printf("rootsFound= %i\n", rootsFound);

    free(values);

    return midPoint;
}

static inline real check(real (*func)(real, real *, dsfmt_t*), real x, real * funcargs, dsfmt_t* dsfmtState)
{
    real h = 0.1;
    real funct;
    real p1, p2, p3, p4, p5, denom;

    p1 = - 1.0 * (*func)( (x + 2.0 * h) , funcargs, dsfmtState);
    p2 =  16.0 * (*func)( (x + h)       , funcargs, dsfmtState);
    p3 = -30.0 * (*func)( (x)           , funcargs, dsfmtState);
    p4 =  16.0 * (*func)( (x - h)       , funcargs, dsfmtState);
    p5 = - 1.0 * (*func)( (x - 2.0 * h) , funcargs, dsfmtState);
    denom = inv( 12.0 * h * h);
    funct = ( density(p1, funcargs, dsfmtState) + density(p2, funcargs, dsfmtState) + density(p3, funcargs, dsfmtState) + density(p4, funcargs, dsfmtState) + density(p5, funcargs, dsfmtState) ) * denom;
    
    return funct;
}


/*      VELOCITY DISTRIBUTION FUNCTION CALCULATION      */
static real fun(real ri, real * args, dsfmt_t* dsfmtState)
{
    
    real energy   = args[4];

    real first_deriv_psi;
    real second_deriv_psi;
    real first_deriv_density;
    real second_deriv_density;
    real dsqden_dpsisq;/*second derivative of density with respect to -potential (psi) */
    real dsqden_dpsisq1;
    real denominator; /*the demoninator of the distribution function: 1/sqrt(E-Psi)*/
    real func;

    first_deriv_psi      = first_derivative(potential, ri, args, dsfmtState);
    first_deriv_density  = first_derivative(density,   ri, args, dsfmtState);

    second_deriv_psi     = second_derivative(potential, ri, args, dsfmtState);
    second_deriv_density = second_derivative(density,   ri, args, dsfmtState);
    
    /*
    * Instead of calculating the second derivative of density with respect to -pot directly, 
    * did product rule since both density and pot are functions of radius. 
    */
    
    
    
    /*
     * we take the absolute value in the squareroot even though we shouldn't. We do this because there is a singularity in the
     * denom. After this occurs, the numbers in the squareroot become negative or: E-psi = neg because psi > E after the singularity.
     * we took the absolute value to avoid NANs from this. It is ok to do this because the alternative would be stopping the procedure
     * just before it goes to the singlularity. Either way, we over estimate or under estimate the denom by the same amount (the step size)
     */
    
    denominator = minushalf( mw_fabs(energy - potential(ri, args, dsfmtState) ) );
    
    if(first_deriv_psi != 0.0  && denominator != 0.0)
    {
            dsqden_dpsisq = second_deriv_density / first_deriv_psi - first_deriv_density * second_deriv_psi / (sqr(first_deriv_psi));
            dsqden_dpsisq *= inv(first_deriv_psi);
    }
    else
    {
        denominator = minushalf( mw_fabs(energy - potential(ri + 0.01, args, dsfmtState) ) );
    }
    
    
    func = first_deriv_psi * dsqden_dpsisq * denominator;
//     mw_printf("%f  %f\n", dsqden_dpsisq, dsqden_dpsisq1);
//     mw_printf("\nfunc = %0.20f f_psi = %0.20f s_psi = %0.20f  f_den = %0.20f s_den = %0.20f  denom = %f, dsdendpsi = %0.20f \n", func, first_deriv_psi, second_deriv_psi, first_deriv_density,second_deriv_density, denominator, dsqden_dpsisq); 
    
    return func;
        
}

static inline real find_upperlimit_r(real * args, real energy, dsfmt_t* dsfmtState, real search_range)
{
    real r = args[4];
    
    int counter = 0;
    real upperlimit_r = 0.0;
    
    do
    {
        upperlimit_r = root_finder(potential, args, energy, 0.0, search_range, dsfmtState); 
        
        if(isinf(upperlimit_r) == FALSE && upperlimit_r != 0.0 && isnan(upperlimit_r) == FALSE){break;}
        
        counter++;
        
        if(counter > 100)
        {
            upperlimit_r = r;
            break;
        }
        
    }while(1);
        
    return mw_fabs(upperlimit_r);
}
 
 
real distribution(real v, real r, real * dwarfargs, real * args, dsfmt_t* dsfmtState)
{
    real mass    = dwarfargs[0];
    real r_scale = dwarfargs[1];
    real coeff = 24.0 * mw_sqrt(2.0) * inv( 7.0 * cube(M_PI) );
    real energy = potential(r, args, dsfmtState) - 0.5 * sqr(v) ;
    real f = v * v * coeff * inv( mw_pow(mass, 4.0 )) * sqr(r_scale) * mw_pow(fabs(energy), 3.5);
    
    return f;
} 

static real dist_fun(real v, real * args, dsfmt_t* dsfmtState)
{
    /*This returns the value of the distribution function*/
    real mass_l   = args[0];
    real mass_d   = args[1];
    real rscale_l = args[2];
    real rscale_d = args[3];
    real ifmax    = args[4];
    real r        = args[5];
    real v_esc    = args[6];
//     mw_printf("ml= %f  md= %f  rl= %f   rd= %f  r= %f  ifmax=%f\n", mass_l, mass_d, rscale_l, rscale_d, r, ifmax);
    
    real distribution_function = 0.0;
    real c = inv( (mw_sqrt(8) * sqr(M_PI)) );
    real energy = 0;
    real upperlimit_r = 0.0;
    real lowerlimit_r = 0.0; 
    int counter = 0.0;
    real search_range = 0.0;   
    
    if(v == v_esc)
    {
        ifmax == 1;
    }
    
    if(ifmax == 1)
    {
        energy = potential(r, args, dsfmtState);
        upperlimit_r = r;
    }
    else if(ifmax == 0)
    {
        energy = potential(r, args, dsfmtState) - 0.5 * v * v; 
        search_range = 20.0 * mw_sqrt( mw_fabs( sqr(mass_d/energy) - sqr(rscale_d) ));
        
        /*dynamic search range*/
        while(potential(search_range, args, dsfmtState) > energy)
        {
            search_range = 100.0 * search_range;
            if(counter > 100)
            {
                search_range = 10 * (rscale_l + rscale_d);
                break;
            }
            counter++;
        }
        upperlimit_r = find_upperlimit_r(args, energy, dsfmtState, search_range);
//         mw_printf("\nenergy = %0.20f\t rnew = %0.20f\n", energy, upperlimit_r);
    }
    
    real funcargs[5] = {mass_l, mass_d, rscale_l, rscale_d, energy};
    lowerlimit_r = 5.0 * (upperlimit_r);
    
    /*This calls guassian quad to integrate the function for a given energy*/
    distribution_function = v * v * c * gauss_quad(fun, lowerlimit_r, upperlimit_r, rscale_d, funcargs, dsfmtState);
    return mw_fabs(distribution_function);
}



real test(real x, real * args, dsfmt_t* dsfmtState)
{
 real f = exp(x) * sin(x) * x;
 return f;
}

real calc_vesc(real r, real * args, dsfmt_t* dsfmtState)
{
    real vesc = mw_sqrt( mw_fabs(2.0 * potential( r, args, dsfmtState) ) );
    return vesc;
}

real calc_energy(real v, real r, real * args, dsfmt_t* dsfmtState)
{
    real e = potential(r, args, dsfmtState) - 0.5 * v * v;
    return e;
}


/*      SAMPLING FUNCTIONS      */
static inline mwvector angles(dsfmt_t* dsfmtState, real rad)
{
    /* assigns angles. Allows for non-circular orbits.*/
    mwvector vec;
    real phi, theta;
    
    /*defining some angles*/
    theta = mw_acos( mwXrandom(dsfmtState, -1.0, 1.0) );
    phi = mwXrandom( dsfmtState, 0.0, 1.0 ) * 2 * M_PI;

    /*this is standard formula for x,y,z components in spherical*/
    X(vec) = rad * sin( theta ) * cos( phi );        /*x component*/
    Y(vec) = rad * sin( theta ) * sin( phi );        /*y component*/
    Z(vec) = rad * cos( theta );                   /*z component*/

    return vec;
}

static inline real r_mag(dsfmt_t* dsfmtState, real * args, real rho_max)
{
    real rscale_l = args[2];
    real rscale_d = args[3];
    
    int counter = 0;
    real r;
    real u;
    real val;
    
    while (1)
    {
        r = (real)mwXrandom(dsfmtState, 0.0, 5.0 * (rscale_l + rscale_d) );
        u = (real)mwXrandom(dsfmtState, 0.0, 1.0);
        val = r * r * density(r, args, dsfmtState);

        if(val/rho_max > u || counter > 1000)
        {
            break;
        }
        else
        {
            counter++;
        }
    }
        
    return r;
}

static inline real vel_mag(dsfmt_t* dsfmtState, real r, real * args)
{
    
    /*
     * WE TOOK IN 
     * MASS IN SIMULATION UNITS, WHICH HAVE THE UNITS OF KPC^3/GY^2 
     * LENGTH IN KPC
     * AND TIME IN GY
     * THEREFORE, THE velocities ARE OUTPUTING IN KPC/GY
     * THIS IS EQUAL TO 0.977813107 KM/S
     * 
     */
    real mass_l   = args[0];
    real mass_d   = args[1];
    real rscale_l = args[2];
    real rscale_d = args[3];
    
    int counter = 0;
    real v, u, d;
    real ifmax = 1;
    real v_esc = mw_sqrt( mw_fabs(2.0 * potential( r, args, dsfmtState) ) );
    
    real parameters[7] = {mass_l, mass_d, rscale_l, rscale_d, ifmax, r, v_esc};
    real dist_max = dist_fun(0.5 * v_esc, parameters, dsfmtState);
    //max_finder(dist_fun, parameters, 0.0, .5*v_esc, v_esc, 10, 1e-2, dsfmtState);
    
    //
//     mw_printf("\n%f\n", dist_max);
    ifmax = 0;
    parameters[4] = ifmax;
   
//     mw_printf("rejection sampling...");
    while(1)
    {

        v = (real)mwXrandom(dsfmtState, 0.0, v_esc);
        u = (real)mwXrandom(dsfmtState, 0.0, 1.0);
        
        d = dist_fun(v, parameters, dsfmtState);
        if(mw_fabs(d/dist_max) > u)
        {
            break;
        }
        
        if(counter > 1000)
        {
            v = 0;
            break;
        }
        else
        {
            counter++;
        }
    }
    

    v *= 0.977813107;//changing from kpc/gy to km/s
    return mw_fabs(v); //km/s
}

static inline mwvector get_vec(dsfmt_t* dsfmtState, mwvector shift, real x)
{
    mwvector vec;

    vec = angles(dsfmtState, x);    /* pick scaled position */
    mw_incaddv(vec, shift);                             /* move the position */

    return vec;
}


/*      DWARF GENERATION        */
static int nbGenerateIsotropicCore(lua_State* luaSt, dsfmt_t* prng, unsigned int nbody, real mass1, real mass2, mwbool ignore, mwvector rShift, mwvector vShift, real radiusScale1, real radiusScale2)
{
    /* generatePlummer: generate Plummer model initial conditions for test
    * runs, scaled to units such that M = -4E = G = 1 (Henon, Heggie,
    * etc).    See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37,
    * 183.
    */
        unsigned int i;
        int table;
        Body b;
        real r, v;
        real mass_l   = mass1; /*mass of the light component*/
        real mass_d   = mass2; /*mass of the dark component*/
        real rscale_l = radiusScale1; /*scale radius of the light component*/
        real rscale_d = radiusScale2; /*scale radius of the dark component*/
        
        /*for normal*/
//         unsigned int half_bodies = nbody / 2;
//         real mass_light_particle = mass_l / (real)(0.5 * nbody);//half the particles are light matter
//         real mass_dark_particle = mass_d / (real)(0.5 * nbody);
        
        
        /*for all dark*/
        unsigned int half_bodies = 0; 
        mass_l = 0.0;
        real dwarfargs[2] = {mass_d, rscale_d};
        
        /*for all light*/
//         int half_bodies = nbody; 
//         mass_d = 0.0;
//         real dwarfargs[2] = {mass_l, rscale_l};
        
        
        real mass_light_particle = mass_l / (real)(nbody);//half the particles are light matter
        real mass_dark_particle = mass_d /  (real)(nbody);//half dark matter
        
        /*dark matter type is TRUE or 1. Light matter type is False, or 0*/
        mwbool isdark = TRUE;
        mwbool islight = FALSE;
        
        real * all_r = mwCalloc(nbody, sizeof(real));
       
        real args[4] = {mass_l, mass_d, rscale_l, rscale_d};
        real parameters_light[4] = {mass_l, 0.0, rscale_l, rscale_d};
        real parameters_dark[4]  = {0.0, mass_d, rscale_l, rscale_d};
        
        /*finding the max of the individual components*/
        real rho_max_light = max_finder(profile_rho, parameters_light, 0, rscale_l, 2.0 * (rscale_l), 20, 1e-4, prng );
        real rho_max_dark  = max_finder(profile_rho, parameters_dark, 0, rscale_d, 2.0 * (rscale_d), 20, 1e-4, prng );
        
//         mw_printf("light max= %f \t dark max= %f\n", rho_max_light, rho_max_dark);
        
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                /*DEBUGGING CODE*/
// {
//         
// //         real tst1 = root_finder(test, args, 4.0, 0.0, 5.0, prng);
// //         mw_printf("test=%f\n", tst1 );
//         
//         real tst1 = gauss_quad(test, 1.0, 5.0, 5.0, args, prng);
//         real tst2 = gauss_quad(test, 5.0, 1.0, 5.0, args, prng);
//         mw_printf("test = %f  %f\n", tst1, tst2);
// //       
// //         real dtst1 = first_derivative(test, 4.0, args, prng);
// //         real dtst2 = second_derivative(test, 4.0, args, prng);
// //         mw_printf("test = %f  %f\n", dtst1, dtst2);
// // 
// //         mw_printf("ml = %f md = %f  rl = %f  rd = %f\n", mass_l, mass_d, rscale_l, rscale_d);
// //         real tst1 = fun(3.0, args, prng);
// //         real tst2 = fun(321.310376, args, prng);
// //         mw_printf("test = %0.20f  %0.20f\n", tst1, tst2);
//         
//         ///////////////////////////////////////////////////////////
//         real r_pot = 0.0;
//         real breakrange = 100 * (rscale_l + rscale_d);
//         FILE * pot;
//         pot = fopen("pot.txt", "w");
//         real pt, pt2, pt3;
//         while(1)
//         {
//             pt  = potential(r_pot, parameters_light, prng);
//             pt2 = potential(r_pot, parameters_dark, prng);
//             pt3 = potential(r_pot, args, prng);
//             r_pot  = r_pot + 0.01;
//             fprintf(pot, "%f \t %f \t %f\t %f\n", pt, r_pot, pt2, pt3);
// //             mw_printf("\r printing density functions: %f %", r_pot / (breakrange) * 100);
//             if(r_pot > 0.1 * breakrange){break;}
//         }
//         fclose(pot);
//         ///////////////////////////////////////////////////////////
//         
//         /////////////////////////////////////////////////////////
//         real r_den = 0.0;
//         FILE * rho;
//         rho = fopen("rho.txt", "w");
//         real den, den2, den3;
//         while(1)
//         {
//             den  = r_den * r_den * density(r_den, parameters_light, prng);
//             den2 = r_den * r_den * density(r_den, parameters_dark, prng);
//             den3 = den + den2;//r_den * r_den * density(r_den, args, prng);
//             r_den   = r_den + 0.01;
//             fprintf(rho, "%f \t %f \t %f\t %f\n", den, r_den, den2, den3);
// //             mw_printf("\r printing density functions: %f %", r_den / (breakrange) * 100);
//             if(r_den > 0.1 * breakrange){break;}
//         }
//         fclose(rho);
//         /////////////////////////////////////////////////////////
//         
//         /////////////////////////////////////////////////////////
//         real r_dist = 0.1;
//         real v_dist = 0;
//         real v_dist_esc = 0.0;
//         int ifmax = 0;
//         real dsb, dsb_theory;
//         real dis_tst_b[7] = {mass_l, mass_d, rscale_l, rscale_d, ifmax, r_dist, v_dist_esc};
//         
//        
// /*         distribution function as a function of r        */
//         FILE * dist2;
//         dist2 = fopen("dist2.txt", "w");
//         
//         /*fixing v*/
//         r_dist = dwarfargs[1];
//         v_dist = 0.5 * calc_vesc(r_dist, args, prng); 
//         
//         /*resetting r*/
//         r_dist = 10 * rscale_d;
//         
//         while(1)
//         {
//             dis_tst_b[5] = r_dist;
//             
//             dsb  = dist_fun(v_dist, dis_tst_b, prng);
//             
//             dsb_theory = distribution(v_dist, r_dist, dwarfargs, args, prng);
//             
//             fprintf(dist2, "%f\t %f \t  %f\n", r_dist, dsb, dsb_theory);
//             
// //             mw_printf("\r printing density functions 1/2: %f %", (0.01 / r_dist ) * 100);
//             r_dist  -= .1;
//             if(r_dist < 0.1){break;}
//         }
//         fclose(dist2);
//         
//         
// /*         distribution function as a function of v        */
//         FILE * dist3;
//         dist3 = fopen("dist3.txt", "w");
//         
//         /*fixing r*/
//         r_dist = dwarfargs[1];
//         dis_tst_b[5] = r_dist;
//         v_dist_esc = calc_vesc(r_dist, args, prng); 
//         
//         v_dist = 0.0;
//         while(1)
//         {
//             dsb  = dist_fun(v_dist, dis_tst_b, prng);
//             
//             dsb_theory = distribution(v_dist, r_dist, dwarfargs, args, prng);
// 
//             fprintf(dist3, "%0.20f \t %2.20f \t %0.20f\n", v_dist, dsb, dsb_theory);
// 
//             v_dist += .001;
//                 if(v_dist > (v_dist_esc)){break;}
//                 
// //             mw_printf("\r printing density functions 2/2: %f %", (v_dist / v_dist_esc) * 100);
//         }
//         fclose(dist3);
//         
//         
//         /////////////////////////////////////////////////////////
//         
//         ///////////////////////////////////////////////////////
// /*         integrand as a function of r prime       */
//         real fun_tst;
//         real energy_fun = 0.0;
//         real funcargs[5] = {mass_l, mass_d, rscale_l, rscale_d, energy_fun};
//         
//         real r_fun = .5; //dwarfargs[1];
//         real v_fun = .5 * calc_vesc(r_fun, args, prng); 
//         energy_fun = calc_energy(v_fun, r_fun, args, prng);
//         real upperlimit_r = root_finder(potential, args, energy_fun, 0.0, 10, prng); 
// //         mw_printf("energy = %f  r= %f\n", energy_fun, upperlimit_r);
//         funcargs[4] = energy_fun;
//         
//         real r_fun_p = 50 * rscale_d;
//         
//         FILE * func;
//         func = fopen("fun.txt", "w");
//         mw_printf("\n printing density function integrand");
//         while(1)
//         {
//             fun_tst = fun(r_fun_p, funcargs, prng);
//             fprintf(func, "%f \t %f \n", r_fun_p, fun_tst);
//             r_fun_p -= 0.01;
//             
// //             if(potential(r_fun_p, args, prng) > energy_fun){break;}
//             if(r_fun_p < 0.01){break;}
//         }
//         fclose(func);
//         
//         
//         
// }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
        memset(&b, 0, sizeof(b));
        lua_createtable(luaSt, nbody, 0);
        table = lua_gettop(luaSt);      
        int counter = 0;
        
        /*getting the radii and velocities for the bodies*/
        for (i = 0; i < nbody; i++)
        {
            counter = 0;
//             mw_printf(" \r initalizing particle %i. ", i+1);
            do
            {
                
                if(i < half_bodies)
                {
                    r = r_mag(prng, parameters_light, rho_max_light);
                }
                else if(i >= half_bodies)
                {
                    r = r_mag(prng, parameters_dark, rho_max_dark);
                }
                /*to ensure that r is finite and nonzero*/
                if(isinf(r) == FALSE && r != 0.0 && isnan(r) == FALSE)
                {
                    break;
                }
                
                if(counter > 1000)
                {
                    exit(-1);
                }
                else
                {
                    counter++;
                }
                
            }while (1);

            /*storing r's*/
            all_r[i] = r;
        }
        

//         mw_printf("\n");
        
        /*this actually gets the position and velocity vectors and pushes table of bodies*/
        for (i = 0; i < nbody; i++)
        {
            r = all_r[i];
            if(i < half_bodies)
            {
                
                b.bodynode.mass = mass_light_particle;
                b.bodynode.type = BODY(islight);

            }
            else if(i >= half_bodies)
            {
                b.bodynode.mass = mass_dark_particle;
                b.bodynode.type = BODY(isdark);
            }
            
            mw_printf("\r velocity of particle %i", i+1);
            counter = 0;
            do
            {
                v = vel_mag(prng, r, args);
                if(isinf(v) == FALSE && v != 0.0 && isnan(v) == FALSE){break;}
                
                if(counter > 1000)
                {
                    exit(-1);
                }
                else
                {
                    counter++;
                }
                
            }while (1);
//                 mw_printf(" vel %i    \t v= %f \t r=%f \n", i, v, r);
            
            /*vShift and rShift are both zero here. They are meant to give the dwarf an initial position and vel*/
            b.vel = get_vec(prng, vShift, v);
            b.bodynode.pos = get_vec(prng, rShift, r);
            
            assert(nbPositionValid(b.bodynode.pos));
            pushBody(luaSt, &b);
            lua_rawseti(luaSt, table, i + 1);
        }
        free(all_r);
        return 1;             
             
}


int nbGenerateIsotropic(lua_State* luaSt)
{
        static dsfmt_t* prng;
        static const mwvector* position = NULL;
        static const mwvector* velocity = NULL;
        static mwbool ignore;
        static real mass1 = 0.0, nbodyf = 0.0, radiusScale1 = 0.0;
        static real mass2 = 0.0, radiusScale2 = 0.0;

        static const MWNamedArg argTable[] =
        {
            { "nbody",                LUA_TNUMBER,     NULL,                    TRUE,    &nbodyf            },
            { "mass1",                LUA_TNUMBER,     NULL,                    TRUE,    &mass1             },
            { "mass2",                LUA_TNUMBER,     NULL,                    TRUE,    &mass2             },
            { "scaleRadius1",         LUA_TNUMBER,     NULL,                    TRUE,    &radiusScale1      },
            { "scaleRadius2",         LUA_TNUMBER,     NULL,                    TRUE,    &radiusScale2      },
            { "position",             LUA_TUSERDATA,   MWVECTOR_TYPE,           TRUE,    &position          },
            { "velocity",             LUA_TUSERDATA,   MWVECTOR_TYPE,           TRUE,    &velocity          },
            { "ignore",               LUA_TBOOLEAN,    NULL,                    FALSE,   &ignore            },
            { "prng",                 LUA_TUSERDATA,   DSFMT_TYPE,              TRUE,    &prng              },
            END_MW_NAMED_ARG
            
        };

        if (lua_gettop(luaSt) != 1)
            return luaL_argerror(luaSt, 1, "Expected 1 arguments");
        
        handleNamedArgumentTable(luaSt, argTable, 1);
        
        return nbGenerateIsotropicCore(luaSt, prng, (unsigned int) nbodyf, mass1, mass2, ignore,
                                                                 *position, *velocity, radiusScale1, radiusScale2);
}

void registerGenerateIsotropic(lua_State* luaSt)
{
    lua_register(luaSt, "generateIsotropic", nbGenerateIsotropic);
}


