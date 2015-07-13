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

    real potential_result = -1.0 * (mass_l/mw_sqrt(sqr(r) + sqr(rscale_l)) + mass_d/mw_sqrt(sqr(r) + sqr(rscale_d)) );

    return (-1.0 * potential_result);
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



/*      GENERAL PURPOSE DERIVATIVE, INTEGRATION, MAX FINDING, ROOT FINDING, AND ARRAY SHUFFLER FUNCTIONS  */
static inline real first_derivative(real (*rootFunc)(real, real *, dsfmt_t*), real x, real * funcargs, dsfmt_t* dsfmtState)
{
    /*yes, this does in fact use a 5-point stencil*/
    real h = 0.01;
    real deriv;
    real p1, p2, p3, p4, denom;
    
    p1 =   1.0 * (*rootFunc)( (x - 2.0 * h), funcargs, dsfmtState);
    p2 = - 8.0 * (*rootFunc)( (x - h)      , funcargs, dsfmtState);
    p3 = - 1.0 * (*rootFunc)( (x + 2.0 * h), funcargs, dsfmtState);
    p4 =   8.0 * (*rootFunc)( (x + h)      , funcargs, dsfmtState);
    denom = inv( 12.0 * h);
    deriv =   (p1 + p2 + p3 + p4) * denom;
    return deriv;
}

static inline real second_derivative(real (*rootFunc)(real, real *, dsfmt_t*), real x, real * funcargs, dsfmt_t* dsfmtState)
{
    /*yes, this also uses a five point stencil*/
    real h = 0.01;
    real deriv;
    real p1, p2, p3, p4, p5, denom;

    p1 = - 1.0 * (*rootFunc)( (x + 2.0 * h) , funcargs, dsfmtState);
    p2 =  16.0 * (*rootFunc)( (x + h)       , funcargs, dsfmtState);
    p3 = -30.0 * (*rootFunc)( (x)           , funcargs, dsfmtState);
    p4 =  16.0 * (*rootFunc)( (x - h)       , funcargs, dsfmtState);
    p5 = - 1.0 * (*rootFunc)( (x - 2.0 * h) , funcargs, dsfmtState);
    denom = inv( 12.0 * h * h);
    deriv =   (p1 + p2 + p3 + p4 + p5) * denom;
    return deriv;
}

static real gauss_quad(real (*rootFunc)(real, real *, dsfmt_t*), real lower, real upper, real * funcargs, dsfmt_t* dsfmtState)
{
    /*This is a guassian quadrature routine. */
    real Ng,hg,lowerg, upperg;
    real intv;
    real coef1,coef2;//parameters for gaussian quad
    real c1,c2,c3;
    real x1,x2,x3;
    real x1n,x2n,x3n;
    
    real a = lower; 
    real b = upper;
    
    intv = 0;//initial value of integral
    Ng = 20.0;//integral resolution
    hg = (b-a)/(Ng);
    
    lowerg = a;
    upperg = lowerg+hg;
    

    coef2 = (lowerg+upperg)/2.0;//initializes the first coeff to change the function limits
    coef1 = (upperg-lowerg)/2.0;//initializes the second coeff to change the function limits
    c1 = 5.0/9.0;
    c2 = 8.0/9.0;
    c3 = 5.0/9.0;
    x1 = -sqrt(3.0/5.0);
    x2 = 0.0;
    x3 = sqrt(3.0/5.0);
    x1n = (coef1 * x1 + coef2);
    x2n = (coef1 * x2 + coef2);
    x3n = (coef1 * x3 + coef2);
    int counter=0;
    while (1)
    {
                //gauss quad
        intv= intv + c1 * (*rootFunc)(x1n, funcargs, dsfmtState) * coef1 +
                     c2 * (*rootFunc)(x2n, funcargs, dsfmtState) * coef1 + 
                     c3 * (*rootFunc)(x3n, funcargs, dsfmtState) * coef1;

        lowerg = upperg;
        upperg = upperg + hg;
        coef2 = (lowerg + upperg)/2.0;//initializes the first coeff to change the function limits
        coef1 = (upperg - lowerg)/2.0;

        x1n = ((coef1) * x1 + coef2);
        x2n = ((coef1) * x2 + coef2);
        x3n = ((coef1) * x3 + coef2);


        if(upper > lower)
        {
            if(lowerg >= upper)//loop termination clause
            {
                break;
            }
        }
        else if(lower > upper)
        {
            if(lowerg <= upper)//loop termination clause
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
    
    
    return intv;
}

static real simpson( real (*rootFunc)(real, real *, dsfmt_t*), real lower_limit, real upper_limit, real * args, dsfmt_t* dsfmtState)
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
        intv = intv + coef * (   (*rootFunc)(lower, args, dsfmtState) + 
                             4 * (*rootFunc)(mid, args, dsfmtState) +
                                 (*rootFunc)(upper, args, dsfmtState) );//gets the integral value for the interval
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
static inline void shuffle(real * init_vec, int length,dsfmt_t* dsfmtState )
{
    /*a quick array shuffling algorithm */
    real * test = mwCalloc(length, sizeof(real));/*new array to return*/
    int * store = mwCalloc(length, sizeof(real));/*stores array coordinates*/
    mwbool good_q = FALSE;
    mwbool dropped_the_cards = FALSE;
    int matches = 0;
    int q;
    int counter = 0;
    /*picks a random array element and assigns that array element of the input
     * array into the new array. To avoid repetition, the assigned array coordinate is stored. 
     * All future coordinates are tested to make sure they were not already assigned.*/
    for(int i = 0; i < length; i++)
        {
            /*The first array element obviously would not have been assigned already*/
            if(i == 0)
            {
                q = (int)(mwXrandom(dsfmtState,0.0,length));
                store[i] = q;
                test[i] = init_vec[q];
            }
            else 
            {
                /*Chooses an array coordinate and checks to see if it has been used already*/
                while(1)
                {
                    q = (int)(mwXrandom(dsfmtState,0.0,length));
                    matches = 0;
                    for(int j = 0; j < i; j++)
                    {
                        if(store[j] == q)
                        {
                            matches++;
                            break;
                        }
                    }
                    
                    /*if that coordinate was not already used, then its good*/
                    if(matches == 0)
                    {
                        good_q = TRUE;
                        break;
                    }
                    /*I do not like infinite loops*/
                    if(counter > 10000)
                    {
                        dropped_the_cards = TRUE;
                        break;
                    }
                    else
                    {
                        counter++;
                    }
                }
                /*stores the choice of coordinate and assigns that element to the new array*/
                store[i] = q;
                test[i] = init_vec[q];
            }
        }
    
        for(int k = 0; k < length; k++)
        {
            /*if loop counter was hit, then shuffling did not work. In that case return original array*/
            if(dropped_the_cards == FALSE)
            {
                init_vec[k] = test[k];
            }
        }
        free(test);
        free(store);
}

static real root_finder(real (*rootFunc)(real, real*, dsfmt_t*), real* rootFuncParams, real funcValue, real lowBound, real upperBound, dsfmt_t* dsfmtState)
{
    //requires lowBound and upperBound to evaluate to opposite sign when rootFunc-funcValue
    if(rootFuncParams == NULL || rootFunc == NULL)
    {
        exit(-1);
    }
    unsigned int i = 0;
    /* Can find up to 20 roots, but may miss a root if they are too close together */
    int N=4;
    unsigned int numSteps = N;
    real interval;
    real * values = mwCalloc(numSteps, sizeof(real));
    int * deck = mwCalloc(numSteps, sizeof(real));
    
    for(int k = 0; k < N; k++)
    {
        deck[k] = k;
    }
    shuffle(deck, numSteps, dsfmtState);
    
    
    /* Divide the function area up into bins in case there is more than one root */
    /*numSteps+1 because you want to include the upperbound in the interval*/
    for(i = 0; i < numSteps+1; i++)
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
        q = deck[i];
        if((values[q] > 0 && values[q+1] < 0) || (values[q] < 0 && values[q+1] > 0))
        {
            if(values[q] < 0 && values[q+1] > 0)
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
    free(deck);
    return midPoint;
}



/*      VELOCITY DISTRIBUTION FUNCTION CALCULATION      */
static real fun(real ri, real * args, dsfmt_t* dsfmtState)
{
    
    real mass_l   = args[0];
    real mass_d   = args[1];
    real rscale_l = args[2];
    real rscale_d = args[3];
    real energy   = args[4];

    real first_deriv_psi;
    real second_deriv_psi;
    real first_deriv_density;
    real second_deriv_density;
    real dsqden_dpsisq;/*second derivative of density with respect to -potential (psi) */
    real denominator; /*the demoninator of the distribution function: 1/sqrt(E-Psi)*/
    real func;


    first_deriv_psi      = first_derivative(potential, ri, args, dsfmtState);
    first_deriv_density  = first_derivative(density, ri, args, dsfmtState);


    second_deriv_psi     = second_derivative(potential, ri, args, dsfmtState);
    second_deriv_density = second_derivative(density, ri, args, dsfmtState);
    
    

    /*
    * Instead of calculating the second derivative of density with respect to -pot directly, 
    * did product rule since both density and pot are functions of radius. 
    */
    if(first_deriv_psi != 0.0)
    {
            dsqden_dpsisq = second_deriv_density / first_deriv_psi - first_deriv_density * second_deriv_psi / (sqr(first_deriv_psi));
    }
    else
    {
            dsqden_dpsisq = 0.0;
    }
//     mw_printf("%0.10f \t %0.10f\t%0.10f\t%0.10f\n", first_deriv_psi, first_deriv_density, second_deriv_density, second_deriv_psi);
    denominator = minushalf( mw_fabs(energy - potential(ri,args,dsfmtState) ) );
    func = first_deriv_psi * dsqden_dpsisq * denominator;
    return func;
        
}

static inline real find_upperlimit_r(real * args, real energy, dsfmt_t* dsfmtState, real v, real search_range)
{
    real mass_l   = args[0];
    real mass_d   = args[1];
    real rscale_l = args[2];
    real rscale_d = args[3];
    real r        = args[4];
    real v_esc = mw_sqrt( mw_fabs(2.0 * potential( r, args, dsfmtState) ) );
    real energy_max = potential(r, args,dsfmtState);
    int counter = 0;
    mwbool limit_reached = FALSE;
    real upperlimit_r = 0.0;
    do
    {
//         mw_printf("\t \t fetching root...\n");
        upperlimit_r = root_finder(potential, args, energy, 0.0, search_range, dsfmtState); 
        
        if(isinf(upperlimit_r) == FALSE && upperlimit_r != 0.0 && isnan(upperlimit_r) == FALSE){break;}
        counter++;
        if(counter > 100)
        {
//             mw_printf("this ran \n");
//             mw_printf("\t \t  v=%f v_esc=%f root=%f energy=%f en_max=%f\n",v, v_esc, upperlimit_r, energy, energy_max);
            upperlimit_r = r;
//             mw_printf("potential range = %f \n ", potential(10 * rscale_d, args,dsfmtState));
//             mw_printf("\t \t  v=%f v_esc=%f root=%f energy=%f en_max=%f\n",v, v_esc, upperlimit_r, energy, energy_max);
            break;
        }
    }while(1);
        
//     mw_printf("upperlimit_r= %f energy=%f\n", upperlimit_r, energy);
    upperlimit_r = mw_fabs(upperlimit_r);
    return upperlimit_r;
}

 
static real dist_fun(real v, real * args, dsfmt_t* dsfmtState)
{
    /*This returns the value of the distribution function*/
    real mass_l   = args[0];
    real mass_d   = args[1];
    real rscale_l = args[2];
    real rscale_d = args[3];
    real r        = args[4];
    real ifmax    = args[5];
    
    real distribution_function = 0.0;
    real c = inv( (mw_sqrt(8) * sqr(M_PI)) );
    real energy = 0;
    real upperlimit_r = 0.0;
    real lowerlimit_r = 0.0; 
    int counter = 0.0;
    real search_range = 0.0;   

    if(ifmax == 1)
    {
        energy = potential(r, args,dsfmtState);
        upperlimit_r = r;
    }
    else if(ifmax == 0)
    {
        energy = potential(r, args,dsfmtState) - 0.5 * v * v; 
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
        upperlimit_r = find_upperlimit_r(args, energy, dsfmtState, v, search_range);
    }
    
    real funcargs[5] = {mass_l, mass_d, rscale_l, rscale_d, energy};
    lowerlimit_r = 50.0 * (rscale_d);
   
//     mw_printf("upper: %f \t lower: %f\n", upperlimit_r, lowerlimit_r);
    /*This calls guassian quad to integrate the function for a given energy*/
    distribution_function = v * v * c * gauss_quad(fun, lowerlimit_r, upperlimit_r, funcargs, dsfmtState);
        
    return mw_fabs(distribution_function);
}



real test(real x, real * args, dsfmt_t* dsfmtState)
{
 real f = exp(x) * sin(x) * x;
 return f;
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
    real mass_l   = args[0];
    real mass_d   = args[1];
    real rscale_l = args[2];
    real rscale_d = args[3];
    
    int counter = 0;
    real r;
    real u;
    real val;
    real bound;
    
    mwbool GOOD_RADIUS = 0;
    if(mass_l = 0.0)
    {
        bound = 5.0 * ( rscale_d);
    }
    else if(mass_d = 0.0)
    {
        bound = 5.0 * (rscale_l );
    }
    while (GOOD_RADIUS != 1)
    {
        r = (real)mwXrandom(dsfmtState, 0.0, 5.0 * (rscale_l + rscale_d) );
        u = (real)mwXrandom(dsfmtState, 0.0, 1.0);
        val = r * r * density(r, args, dsfmtState);

        if(val/rho_max > u || counter > 1000)
        {
            GOOD_RADIUS = 1;
        }
        else
        {
            counter++;
        }
    }
        
    return r;
}

static inline real vel_mag(dsfmt_t* dsfmtState,real r, real * args)
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
    real v,u, d;
    real ifmax = 1;
    real v_esc = mw_sqrt( mw_fabs(2.0 * potential( r, args, dsfmtState) ) );
    
    real parameters[6] = {mass_l, mass_d, rscale_l, rscale_d, r, ifmax};
    real dist_max = dist_fun(v_esc, parameters, dsfmtState);
    //max_finder(dist_fun, parameters, 0.0, .5*v_esc, v_esc, 10, 1e-2, dsfmtState);
    
    
    ifmax = 0;
    parameters[5] = ifmax;
    
//     mw_printf("rejection sampling...");
    while(1)
    {

        v = (real)mwXrandom(dsfmtState, 0.0, v_esc);
        u = (real)mwXrandom(dsfmtState, 0.0, 1.0);
        
        d = dist_fun(v, parameters, dsfmtState);
        
        if(mw_fabs(d/dist_max) > u || counter > 1000)
        {
            break;
        }
        else
        {
            counter++;
        }
    }

    v *= 0.977813107;//changing from kpc/gy to km/s
    return fabs(v); //km/s
}

static inline mwvector r_vec(dsfmt_t* dsfmtState, mwvector rshift, real r)
{
    mwvector pos;

    pos = angles(dsfmtState, r);    /* pick scaled position */
    mw_incaddv(pos, rshift);                             /* move the position */

    return pos;
}


static inline mwvector vel_vec(dsfmt_t* dsfmtState, mwvector vshift, real v)
{
    mwvector vel;

    vel = angles(dsfmtState, v);     /* pick scaled velocity */
    mw_incaddv(vel, vshift);                                /* move the velocity */

    return vel;
}


/*      DWARF GENERATION        */
/* generatePlummer: generate Plummer model initial conditions for test
 * runs, scaled to units such that M = -4E = G = 1 (Henon, Heggie,
 * etc).    See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37,
 * 183.
 */

static int nbGenerateIsotropicCore(lua_State* luaSt, dsfmt_t* prng, unsigned int nbody, 
                                   real mass1, real mass2, 
                                   mwbool ignore, mwvector rShift, mwvector vShift, 
                                   real radiusScale1, real radiusScale2)
{
        unsigned int i;
        int table;
        Body b;
        real r, v;
        real mass_l   = mass1;
        real mass_d   = mass2;
        real rscale_l = radiusScale1;
        real rscale_d = radiusScale2;
        
        real half_bodies = 0.5 * nbody;
        real mass_light_particle = mass_l / (half_bodies);//half the particles are light matter
        real mass_dark_particle = mass_d / (half_bodies);//half dark matter
        
        /*dark matter type is TRUE or 1. Light matter type is False, or 0*/
        mwbool isdark = TRUE;
        mwbool islight = FALSE;
        int dark = 1;
        int light = 0;
        
        real * dark_r = mwCalloc(half_bodies, sizeof(real));
        real * light_r = mwCalloc(half_bodies, sizeof(real));
       
        real args[4] = {mass_l, mass_d, rscale_l, rscale_d};
        real parameters_light[4] = {mass_l, 0.0, rscale_l, rscale_d};
        real parameters_dark[4] = {0.0, mass_d, rscale_l, rscale_d};
        
//         real tst1 = root_finder(test, args, 4.0, 0.0, 5.0, prng);
//         mw_printf("test=%f\n", tst1 );
        
        real tst1 = gauss_quad(test, 1.0, 5.0, args, prng);
        real tst2 = gauss_quad(test, 5.0, 1.0, args, prng);
        mw_printf("test = %f  %f\n", tst1, tst2);
        real dtst1 = first_derivative(test, 4.0, args, prng);
        real dtst2 = second_derivative(test, 4.0, args, prng);
        mw_printf("test = %f  %f\n", dtst1, dtst2);
        
        /*finding the max of the individual components*/
        real rho_max_light = max_finder(profile_rho, parameters_light, 0, rscale_l, 2.0*(rscale_l), 20, 1e-4, prng );
        real rho_max_dark =  max_finder(profile_rho, parameters_dark, 0, rscale_d, 2.0*(rscale_d), 20, 1e-4, prng );
        real rho_max = max_finder(profile_rho, args, 0, rscale_l, (rscale_d), 20, 1e-4, prng );
//         mw_printf("light max= %f \t dark max= %f\n", rho_max_light, rho_max_dark);
        
        
//         ///////////////////////////////////////////////////////////
//         real w1 = 0.0;
//         FILE * pot;
//         pot = fopen("pot.txt", "w");
//         real pt, pt2, pt3;
//         while(1)
//         {
//             pt  = potential(w1, parameters_light, prng);
//             pt2 = potential(w1, parameters_dark, prng);
//             pt3 = potential(w1, args, prng);
//             w1  = w1 + 0.01;
//             fprintf(pot, "%f \t %f \t %f\t %f\n", pt, w1, pt2, pt3);
// //             mw_printf("\r printing density functions: %f %", w/(5*(rscale_l+rscale_d))*100);
//             if(w1 > 5 * (rscale_l+rscale_d)){break;}
//         }
//         fclose(pot);
//         ///////////////////////////////////////////////////////////
//         
//         /////////////////////////////////////////////////////////
//         real w = 0.0;
//         FILE * rho;
//         rho = fopen("rho.txt", "w");
//         real de, de2, de3;
//         while(1)
//         {
//             de  = w * w * density(w, parameters_light, prng);
//             de2 = w * w * density(w, parameters_dark, prng);
//             de3 = w * w * density(w, args, prng);
//             w   = w + 0.01;
//             fprintf(rho, "%f \t %f \t %f\t %f\n", de, w, de2, de3);
// //             mw_printf("\r printing density functions: %f %", w/(5*(rscale_l+rscale_d))*100);
//             if(w > 5 * (rscale_l+rscale_d)){break;}
//         }
//         fclose(rho);
//         /////////////////////////////////////////////////////////
        
//         /////////////////////////////////////////////////////////
//         real w2 = 0.1;
//         FILE * dist;
//         dist = fopen("dist.txt", "w");
//         real dsl, dsd, dsb;
//         real breakrange = 100 * (rscale_l+rscale_d);
//         int ifmax = 0;
//         real tst_l[6] = {mass_l, 0.0,    rscale_l, rscale_d, w2, ifmax};
//         real tst_d[6] = {0.0,    mass_d, rscale_l, rscale_d, w2, ifmax};
//         real tst_b[6] = {mass_l, mass_d, rscale_l, rscale_d, w2, ifmax};
//         real vtst = 0;
//         real vsc_l, vsc_d, vsc_b;
//         
//         while(1)
//         {
//             vsc_l = mw_sqrt( mw_fabs(2.0 * potential( w2, tst_l, prng) ) );
//             vsc_d = mw_sqrt( mw_fabs(2.0 * potential( w2, tst_d, prng) ) );
//             vsc_b = mw_sqrt( mw_fabs(2.0 * potential( w2, tst_b, prng) ) );
//             vtst = 0.0;
//             while(1)
//             {
//                 tst_l[4] = w2;
//                 tst_d[4] = w2;
//                 tst_b[4] = w2;
//                 dsl = dist_fun(vtst, tst_l, prng);
//                 dsd = dist_fun(vtst, tst_d, prng);
//                 dsb = dist_fun(vtst, tst_b, prng);
//                 fprintf(dist, "%f \t %f \t %f \t %f \t %f \n", w2, vtst,  dsl, dsd, dsb);
//                 vtst = vtst + .1;
//                 
//                  if(vtst > (vsc_b)){break;}
//             }
//             mw_printf("\r printing density functions: %f %", (w2 / breakrange) * 100);
//             w2  = w2 + .01;
//             if(w2 > breakrange){break;}
//         }
//         fclose(dist);
//         /////////////////////////////////////////////////////////
        
        
        
        
      
        memset(&b, 0, sizeof(b));
        lua_createtable(luaSt, nbody, 0);
        table = lua_gettop(luaSt);      
        int counter = 0;
        int j=0;
        /*getting the radii and velocities for the bodies*/
        for (i = 0; i < nbody; i++)
        {
            counter = 0;
//             mw_printf(" \r initalizing particle %i. ", j+1);
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
//                 mw_printf(" particle %i    \t r= %f \n", i, r);

            /*storing r's, mass type*/
            if(i < half_bodies)
            {
                light_r[i] = r;
            }
            else if(i >= half_bodies)
            {
                dark_r[j] = r;
                j++;
            }
        }
        

//         mw_printf("\n");
        
        /*this actually gets the position and velocity vectors and pushes table of bodies*/
        j=0;
//         mw_printf("\n");
        for (i = 0; i < nbody; i++)
        {
            if(i < half_bodies)
            {
                r = light_r[i];
                b.bodynode.mass = mass_light_particle;
                b.bodynode.type = BODY(islight);
            }
            else if(i >= half_bodies)
            {
                r = dark_r[j];
                b.bodynode.mass = mass_dark_particle;
                b.bodynode.type = BODY(isdark);
                j++;
            }
            
//             mw_printf("\r velocity of  particle %i", i+1);
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
            
            b.vel = vel_vec(prng, vShift, v);
            b.bodynode.pos = r_vec(prng, rShift, r);
//             mw_printf("%f\n", mw_acos(b.bodynode.pos.z/r));
            assert(nbPositionValid(b.bodynode.pos));
            pushBody(luaSt, &b);
            lua_rawseti(luaSt, table, i + 1);
        }
        free(light_r);
        free(dark_r);
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
            { "scaleRadius1", LUA_TNUMBER,     NULL,                    TRUE,    &radiusScale1},
            { "scaleRadius2", LUA_TNUMBER,     NULL,                    TRUE,    &radiusScale2},
            { "position",         LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,    &position        },
            { "velocity",         LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,    &velocity        },
            { "ignore",             LUA_TBOOLEAN,    NULL,                    FALSE, &ignore            },
            { "prng",                 LUA_TUSERDATA, DSFMT_TYPE,        TRUE,    &prng                },
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


