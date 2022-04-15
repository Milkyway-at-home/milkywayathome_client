/* Copyright (c) 2016 Siddhartha Shelton
  
  Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
  Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
   Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
   Rensselaer Polytechnic Institute.

   Copyright (c) 2016-2018 Siddhartha Shelton

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
    //-------------------------------
    real mass_l   = args[0];
    real mass_d   = args[1];
    real rscale_l = args[2];
    real rscale_d = args[3];
    //-------------------------------
    real potential_light  = mass_l / mw_sqrt(sqr(r) + sqr(rscale_l));
    real potential_dark   = mass_d / mw_sqrt(sqr(r) + sqr(rscale_d));
    real potential_result = -(potential_light + potential_dark);

    return (-potential_result);
}

static inline real density( real r, real * args, dsfmt_t* dsfmtState)
{
    /*this is the density distribution function. Returns the density at a given radius.*/
    //-------------------------------
    real mass_l   = args[0];
    real mass_d   = args[1];
    real rscale_l = args[2];
    real rscale_d = args[3];
    //-------------------------------
    
    real rscale_lCube = cube(rscale_l); 
    real rscale_dCube = cube(rscale_d);
    real density_light = (mass_l / rscale_lCube) * (minusfivehalves( (1.0 + sqr(r)/sqr(rscale_l)) ) );
    real density_dark  = (mass_d / rscale_dCube) * (minusfivehalves( (1.0 + sqr(r)/sqr(rscale_d)) ) ); 
    real density_result = (3.0 / (4.0 * M_PI)) * ( density_light + density_dark );

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
    real h = 0.001;
    real deriv;
    real p1, p2, p3, p4, denom;
    
    p1 =   1.0 * (*func)( (x - 2.0 * h), funcargs, dsfmtState);
    p2 = - 8.0 * (*func)( (x - h)      , funcargs, dsfmtState);
    p3 = - 1.0 * (*func)( (x + 2.0 * h), funcargs, dsfmtState);
    p4 =   8.0 * (*func)( (x + h)      , funcargs, dsfmtState);
    denom = inv( 12.0 * h);
    deriv = (p1 + p2 + p3 + p4) * denom;
    return deriv;
}

static inline real second_derivative(real (*func)(real, real *, dsfmt_t*), real x, real * funcargs, dsfmt_t* dsfmtState)
{
    /*yes, this also uses a five point stencil*/
    real h = 0.001;
    real deriv;
    real p1, p2, p3, p4, p5, denom;

    p1 = - 1.0 * (*func)( (x + 2.0 * h) , funcargs, dsfmtState);
    p2 =  16.0 * (*func)( (x + h)       , funcargs, dsfmtState);
    p3 = -30.0 * (*func)( (x)           , funcargs, dsfmtState);
    p4 =  16.0 * (*func)( (x - h)       , funcargs, dsfmtState);
    p5 = - 1.0 * (*func)( (x - 2.0 * h) , funcargs, dsfmtState);
    denom = inv( 12.0 * h * h);
    deriv = (p1 + p2 + p3 + p4 + p5) * denom;
    return deriv;
}

static real gauss_quad(real (*func)(real, real *, dsfmt_t*), real lower, real upper, real * funcargs, dsfmt_t* dsfmtState)
{
    /*This is a guassian quadrature routine. It will test to always integrate from the lower to higher of the two limits.
     * If switching the order of the limits was needed to do this then the negative of the integral is returned.
     */
    real Ng, hg, lowerg, upperg;
    real intv = 0.0;//initial value of integral
    real coef1, coef2;//parameters for gaussian quad
    real c1, c2, c3;
    real x1, x2, x3;
    real x1n, x2n, x3n;
    real a, b;
    real benchmark;
    
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
    c1 = 0.55555555555; //5.0 / 9.0;
    c2 = 0.88888888888; //8.0 / 9.0;
    c3 = 0.55555555555; //5.0 / 9.0;
    x1 = -0.77459666924;//-sqrt(3.0 / 5.0);
    x2 = 0.00000000000;
    x3 = 0.77459666924; //sqrt(3.0 / 5.0);
    x1n = (coef1 * x1 + coef2);
    /*should be: x2n = (coef1 * x2 + coef2);*/
    x2n = (coef2);
    x3n = (coef1 * x3 + coef2);
    int counter = 0;
    while (1)
    {
                //gauss quad
        intv = intv + c1 * (*func)(x1n, funcargs, dsfmtState) * coef1 +
                      c2 * (*func)(x2n, funcargs, dsfmtState) * coef1 + 
                      c3 * (*func)(x3n, funcargs, dsfmtState) * coef1;

        lowerg = upperg;
        upperg = upperg + hg;
        coef2 = (lowerg + upperg) / 2.0;//initializes the first coeff to change the function limits
        coef1 = (upperg - lowerg) / 2.0;

        x1n = ((coef1) * x1 + coef2);
        /*should be: x2n = (coef1 * x2 + coef2);*/
        x2n = (coef2);
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
        intv *= -1.0;
    }
    
    return intv;
}

static inline real max_finder(real (*profile)(real , real*, dsfmt_t*), real* profileParams, real a, real b, real c, int limit, real tolerance, dsfmt_t* dsfmtState)
{
    /*this is a maxfinding routine to find the maximum of the density.
     * It uses Golden Section Search as outlined in Numerical Recipes 3rd edition
     */
    real RATIO = 0.61803399;
    real RATIO_COMPLEMENT = 1.0 - RATIO;
    int counter = 0;
    
    real profile_x1, profile_x2, x0, x1, x2, x3;
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
            profile_x2 = -(*profile)(x2, profileParams, dsfmtState);
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


static inline real root_finder(real (*func)(real, real*, dsfmt_t*), real* function_parameters, real function_value, real lower_bound, real upper_bound, dsfmt_t* dsfmtState)
{
    //requires lower_bound and upper_bound to evaluate to opposite sign when func-function_value
    if(function_parameters == NULL || func == NULL)
    {
        exit(-1);
    }
    unsigned int i = 0;

    int N = 4;
    unsigned int intervals = N;
    real interval_bound;

    /*interval + 1 because for N intervals there are N + 1 values*/
    real * values = mwCalloc(intervals + 1, sizeof(real));
    real * interval_bounds = mwCalloc(intervals + 1, sizeof(real));
    /*intervals+1 because you want to include the upperbound in the interval*/
    for(i = 0; i < intervals + 1; i++)
    {
        /* breaking up the range between bounds into smaller intervals*/
        interval_bound = ((upper_bound - lower_bound) * (real)i) / (real)intervals + lower_bound;
        interval_bounds[i] = interval_bound;
        /*function value at those intervals*/
        values[i] = (*func)(interval_bound, function_parameters, dsfmtState) - function_value;
    }
    
    real mid_point = 0;
    real mid_point_funcval = 0;
    unsigned int counter = 0;
    real new_upper_bound = 0;
    real new_lower_bound = 0;
    int roots_found = 0;
    int q = 0;
    
    /* Find the roots using bisection because it was easy to code and good enough for our purposes 
     * this will hope around the different intervals until it checks all of them. This way it does not 
     * favor any root.
     */
    for(i = 0; i < intervals; i++)
    {
        q = i;
        if((values[q] > 0 && values[q + 1] < 0) || (values[q] < 0 && values[q + 1] > 0))
        {
            if(values[q] < 0 && values[q + 1] > 0)
            {
                
                new_lower_bound = interval_bounds[q];
                new_upper_bound = interval_bounds[q + 1];
            }
            else if(values[q] > 0 && values[q + 1] < 0)
            {
                
                new_lower_bound = interval_bounds[q + 1];
                new_upper_bound = interval_bounds[q];
            }
            else
            {
                continue;
            }
            
            mid_point_funcval = 1;
            counter = 0;
            while(mw_fabs(mid_point_funcval) > .0001)
            {
                mid_point = (new_lower_bound + new_upper_bound) / 2.0;
                mid_point_funcval = (*func)(mid_point, function_parameters, dsfmtState) - function_value;
                
                if(mid_point_funcval < 0.0)
                {
                    new_lower_bound = mid_point;
                }
                else
                {
                    new_upper_bound = mid_point;
                }
                counter++;
                
                if(counter > 10000)
                {
                    break;
                }
            }
            
            /* If it found a sign change, then the root finder definitly got close. So it will always say it found one. */
            roots_found++;
            
        }
        
        if(roots_found != 0)
        {
            break;
        }
    }

    if(roots_found == 0)
    {
        mid_point = 0.0;
    }
    
    free(values);
    free(interval_bounds);

    return mid_point;
}

/*      VELOCITY DISTRIBUTION FUNCTION CALCULATION      */
real fun(real ri, real * args, dsfmt_t* dsfmtState)
{
    //-------------------------------    
    real energy   = args[4];
    //-------------------------------
    
    real first_deriv_psi;
    real second_deriv_psi;
    real first_deriv_density;
    real second_deriv_density;
    real dsqden_dpsisq;/*second derivative of density with respect to -potential (psi) */
    real denominator; /*the demoninator of the distribution function: 1/sqrt(E-Psi)*/
    real diff;
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
    
    
    
    /*just in case*/
    if(first_deriv_psi == 0.0)
    {
        first_deriv_psi = 1.0e-6;//this should be small enough
    }
    
    dsqden_dpsisq = second_deriv_density * inv(first_deriv_psi) - first_deriv_density * second_deriv_psi * inv(sqr(first_deriv_psi));
    diff = mw_fabs(energy - potential(ri, args, dsfmtState));
    
    
    /*we don't want to have a 0 in the demon*/
    if(diff != 0.0)
    {
        denominator = minushalf( mw_fabs(energy - potential(ri, args, dsfmtState) ) );
    }
    else
    {
        /*if the r is exactly at the singularity then move it a small amount.*/
        denominator = minushalf( mw_fabs(energy - potential(ri + 0.0001, args, dsfmtState) ) );
    }
    
    
    /*
     * the second derivative term should be divided by the first derivate of psi. 
     * However, from changing from dpsi to dr we multiply by first derivative of psi. 
     * Since these undo each other we left them out completely.
     */
    
    func = dsqden_dpsisq * denominator; 
    
    return func;
        
}

static inline real find_upperlimit_r(dsfmt_t* dsfmtState, real * args, real energy, real search_range, real r)
{
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
 
static inline real dist_fun(real v, real * args, dsfmt_t* dsfmtState)
{
    /*This returns the value of the distribution function*/
    
    //-------------------------------
    real mass_l   = args[0];
    real mass_d   = args[1];
    real rscale_l = args[2];
    real rscale_d = args[3];
    real r        = args[4];
    //-------------------------------
    
    
    real distribution_function = 0.0;
//     real c = inv( (mw_sqrt(8.0) * sqr(M_PI)) );
    real c = 0.03582244801567226;
    real energy = 0.0;
    real upperlimit_r = 0.0;
    real lowerlimit_r = 0.0; 
    int counter = 0;
    real search_range = 0.0;   
    
    /*energy as defined in binney*/
    energy = potential(r, args, dsfmtState) - 0.5 * v * v; 
    
    /*this starting point is 20 times where the dark matter component is equal to the energy, since the dark matter dominates*/
    search_range = 20.0 * mw_sqrt( mw_fabs( sqr(mass_d / energy) - sqr(rscale_d) ));
    
    /*dynamic search range*/
    /* This is done this way because we are searching for the r' where:
     * psi(r') = energy = psi(r) - .5 v^2
     * since psi is a positive quantity, the right hand side is always less than/equal to psi(r),
     * this corresponds to larger r (smaller psi). Therefore,
     * as long as the psi(r') > energy we continue to expand the search range
     * in order to have that energy inside the search range,
     * we want to be able to find a root within the search range, so we make sure that the range includes the root.
     * By this, we mean that we want to find a root within a range (r1, r2), where 
     * psi(r1) > energy and psi(r2) < energy
     */
    while(potential(search_range, args, dsfmtState) > energy)
    {
        search_range = 100.0 * search_range;
        if(counter > 100)
        {
            search_range = 100.0 * (rscale_l + rscale_d);//default
            break;
        }
        counter++;
    }
    
    upperlimit_r = find_upperlimit_r(dsfmtState, args, energy, search_range, r);
    
    
    real funcargs[5] = {mass_l, mass_d, rscale_l, rscale_d, energy};
    
    /*This lowerlimit should be good enough. In the important case where the upperlimit is small (close to the singularity in the integrand)
     * then 5 times it is already where the integrand is close to 0 since it goes to 0 quickly. 
     */
    lowerlimit_r = 5.0 * (upperlimit_r);
    
    /*This calls guassian quad to integrate the function for a given energy*/
    distribution_function = v * v * c * gauss_quad(fun, lowerlimit_r, upperlimit_r, funcargs, dsfmtState);
    
    return distribution_function;
}


/*      SAMPLING FUNCTIONS      */
static inline real r_mag(dsfmt_t* dsfmtState, real * args, real rho_max, real bound)
{
    int counter = 0;
    real r, u, val;
    
    while (1)
    {
        r = (real)mwXrandom(dsfmtState, 0.0, bound);
        u = (real)mwXrandom(dsfmtState, 0.0, 1.0);
        val = r * r * density(r, args, dsfmtState);

        if(val / rho_max > u)
        {
            break;
        }
        
        if(counter > 1000)
        {
            r = 0;
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
     * WE TOOK IN MASS IN SIMULATION UNITS, WHICH HAVE THE UNITS OF KPC^3/GY^2 
     * LENGTH IN KPC AND TIME IN GY THEREFORE, THE velocities ARE OUTPUTING IN KPC/GY
     * THIS IS EQUAL TO 0.977813107 KM/S
     */
    
    //-------------------------------
    real mass_l   = args[0];
    real mass_d   = args[1];
    real rscale_l = args[2];
    real rscale_d = args[3];
    //-------------------------------
    
    int counter = 0;
    real v, u, d;
    real v_esc = mw_sqrt( mw_fabs(2.0 * potential( r, args, dsfmtState) ) );
    
    real parameters[5] = {mass_l, mass_d, rscale_l, rscale_d, r};
    real dist_max = max_finder(dist_fun, parameters, 0.0, 0.5 * v_esc, v_esc, 10, 1.0e-2, dsfmtState);
   
    while(1)
    {

        v = (real)mwXrandom(dsfmtState, 0.0, 1.0) * v_esc;
        u = (real)mwXrandom(dsfmtState, 0.0, 1.0);
        
        d = dist_fun(v, parameters, dsfmtState);
        if(mw_fabs(d / dist_max) > u)
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
    

//     v *= 0.977813107;//changing from kpc/gy to km/s
    return v; //km/s
}

static inline mwvector get_components(dsfmt_t* dsfmtState, real rad)
{
    /* assigns angles. Allows for non-circular orbits.*/
    mwvector vec;
    real phi, theta;
    
    /*defining some angles*/
    theta = mw_acos( mwXrandom(dsfmtState, -1.0, 1.0) );
    phi = mwXrandom( dsfmtState, 0.0, 1.0 ) * 2.0 * M_PI;

    /*this is standard formula for x,y,z components in spherical*/
    X(vec) = rad * mw_sin( theta ) * mw_cos( phi );        /*x component*/
    Y(vec) = rad * mw_sin( theta ) * mw_sin( phi );        /*y component*/
    Z(vec) = rad * mw_cos( theta );                   /*z component*/

    return vec;
}


static int cm_correction(real * x, real * y, real * z, real * vx, real * vy, real * vz, real * mass, mwvector rShift, mwvector vShift, real dwarf_mass, int nbody)
{
    /*  
     * This function takes the table of bodies produced and zeroes the center of mass 
     * and center of momentum. It then shifts the center of mass and center of momentum
     * to the expected value for its position in the orbit.
     */
    real cm_x = 0.0;
    real cm_y = 0.0;
    real cm_z = 0.0;
    real cm_vx = 0.0;
    real cm_vy = 0.0;
    real cm_vz = 0.0;
    int i;
    for(i = 0; i < nbody; i++)
    {
        cm_x += mass[i] * x[i];
        cm_y += mass[i] * y[i];
        cm_z += mass[i] * z[i];
        
        cm_vx += mass[i] * vx[i];
        cm_vy += mass[i] * vy[i];
        cm_vz += mass[i] * vz[i];
    }
     
    cm_x = cm_x / (dwarf_mass);
    cm_y = cm_y / (dwarf_mass);
    cm_z = cm_z / (dwarf_mass);
    
    cm_vx = cm_vx / (dwarf_mass);
    cm_vy = cm_vy / (dwarf_mass);
    cm_vz = cm_vz / (dwarf_mass);
    
    for(i = 0; i < nbody; i++)
    {
        x[i] = x[i] - cm_x + rShift.x;
        y[i] = y[i] - cm_y + rShift.y;
        z[i] = z[i] - cm_z + rShift.z;
        
        vx[i] = vx[i] - cm_vx + vShift.x;
        vy[i] = vy[i] - cm_vy + vShift.y;
        vz[i] = vz[i] - cm_vz + vShift.z;
    }
    

    return 1;
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
 
        real * x  = mwCalloc(nbody, sizeof(real));
        real * y  = mwCalloc(nbody, sizeof(real));
        real * z  = mwCalloc(nbody, sizeof(real));
        real * vx = mwCalloc(nbody, sizeof(real));
        real * vy = mwCalloc(nbody, sizeof(real));
        real * vz = mwCalloc(nbody, sizeof(real));
        real * masses = mwCalloc(nbody, sizeof(real));
        
        mwvector vec;
        real dwarf_mass = mass1 + mass2;
        
        
        real mass_l   = mass1; /*mass of the light component*/
        real mass_d   = mass2; /*mass of the dark component*/
        real rscale_l = radiusScale1; /*scale radius of the light component*/
        real rscale_d = radiusScale2; /*scale radius of the dark component*/
        
        real bound = 50.0 * (rscale_l + rscale_d);

    //---------------------------------------------------------------------------------------------------        
        /*for normal*/
        unsigned int half_bodies = nbody / 2;
        real mass_light_particle = mass_l / (real)(0.5 * (real) nbody);//half the particles are light matter
        real mass_dark_particle = mass_d / (real)(0.5 * (real) nbody);
    //----------------------------------------------------------------------------------------------------

        /*dark matter type is TRUE or 1. Light matter type is False, or 0*/
        mwbool isdark = TRUE;
        mwbool islight = FALSE;
        
       /*since the potential and density are a sum of the two components, setting the mass of one 
        * component to zero effectively gives a single component potential and density. 
        */
        real args[4] = {mass_l, mass_d, rscale_l, rscale_d};
        real parameters_light[4] = {mass_l, 0.0, rscale_l, rscale_d};
        real parameters_dark[4]  = {0.0, mass_d, rscale_l, rscale_d};

        /*finding the max of the individual components*/
        real rho_max_light = max_finder(profile_rho, parameters_light, 0, rscale_l, 2.0 * (rscale_l), 20, 1e-4, prng );
        real rho_max_dark  = max_finder(profile_rho, parameters_dark, 0, rscale_d, 2.0 * (rscale_d), 20, 1e-4, prng );
        
     
     /*initializing particles:*/
        memset(&b, 0, sizeof(b));
        lua_createtable(luaSt, nbody, 0);
        table = lua_gettop(luaSt);      
        int counter = 0;
        
        /*getting the radii and velocities for the bodies*/
        for (i = 0; i < nbody; i++)
        {
            counter = 0;
            do
            {
                
                if(i < half_bodies)
                {
                    r = r_mag(prng, parameters_light, rho_max_light, bound);
                    masses[i] = mass_light_particle;
                }
                else if(i >= half_bodies)
                {
                    r = r_mag(prng, parameters_dark, rho_max_dark, bound);
                    masses[i] = mass_dark_particle;
                }
                /*to ensure that r is finite and nonzero*/
                if(isinf(r) == FALSE && r != 0.0 && isnan(r) == FALSE){break;}
                
                if(counter > 1000)
                {
                    exit(-1);
                }
                else
                {
                    counter++;
                }
                
            }while (1);
            
            
            
//             mw_printf("\r velocity of particle %i", i+1);
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

            vec = get_components(prng, v);   
            vx[i] = vec.x;
            vy[i] = vec.y;
            vz[i] = vec.z;
            
            vec = get_components(prng, r);  
            x[i] = vec.x;
            y[i] = vec.y;
            z[i] = vec.z;
        }
        
        /* getting the center of mass and momentum correction */
        cm_correction(x, y, z, vx, vy, vz, masses, rShift, vShift, dwarf_mass, nbody);


        /* pushing the bodies */
        for (i = 0; i < nbody; i++)
        {
            b.bodynode.id = i + 1;
            if(i < half_bodies)
            {
                b.bodynode.type = BODY(islight);
            }
            else if(i >= half_bodies)
            {
                b.bodynode.type = BODY(isdark);
            }
            
            b.bodynode.mass = masses[i];
            /*this actually gets the position and velocity vectors and pushes table of bodies*/
            /*They are meant to give the dwarf an initial position and vel*/
            /* you have to work for your bodynode */
            b.bodynode.pos.x = x[i];
            b.bodynode.pos.y = y[i];
            b.bodynode.pos.z = z[i];
            
            b.vel.x = vx[i];
            b.vel.y = vy[i];
            b.vel.z = vz[i];
            
            assert(nbPositionValid(b.bodynode.pos));
            pushBody(luaSt, &b);
            lua_rawseti(luaSt, table, i + 1);
        }
        
        /* go now and be free!*/
        free(x);
        free(y);
        free(z);
        free(vx);
        free(vy);
        free(vz);
        free(masses);
        
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


// As this code runs, know that it is running on the rotting corpses of a thousand bugs.