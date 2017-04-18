/* Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
     Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

Copyright (c) 2016 Siddhartha Shelton
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
#include "nbody_dwarf_potential.h"
#include "nbody_mixeddwarf.h"
#include "nbody_types.h"
#include "nbody_potential_types.h"

/*Note: minusfivehalves(x) raises to x^-5/2 power and minushalf(x) is x^-1/2*/


/*      MODEL SPECIFIC FUNCTIONS       */
static inline real potential( real r, const Dwarf* comp1, const Dwarf* comp2)
{
    /*Be Careful! this function returns the negative of the potential! this is the value of interest, psi*/
    real potential_light  = get_potential(comp1, r);
    real potential_dark   = get_potential(comp2, r);
    real potential_result = (potential_light + potential_dark);

    return (potential_result);
}

static inline real density( real r, const Dwarf* comp1, const Dwarf* comp2)
{
    /*this is the density distribution function. Returns the density at a given radius.*/
    
    real density_light = get_density(comp1, r);
    real density_dark  = get_density(comp2, r);
    real density_result = (density_light + density_dark );

    return density_result;
}

static inline real profile_rho(real r, real r_placeholder, 
                               const Dwarf* comp, const Dwarf* comp_placeholder)
{
    //there are two place holders because max finder is also used for dis_fun
    real result = r * r * get_density(comp, r);    
    return result;
}

/*      GENERAL PURPOSE DERIVATIVE, INTEGRATION, MAX FINDING, ROOT FINDING, AND ARRAY SHUFFLER FUNCTIONS        */
static inline real first_derivative(real (*func)(real, const Dwarf*, const Dwarf*), real x, const Dwarf* comp1, const Dwarf* comp2)
{
    /*yes, this does in fact use a 5-point stencil*/
    real h = 0.001;
    real deriv;
    real p1, p2, p3, p4, denom;
    
    p1 =   1.0 * (*func)( (x - 2.0 * h), comp1, comp2);
    p2 = - 8.0 * (*func)( (x - h)      , comp1, comp2);
    p3 = - 1.0 * (*func)( (x + 2.0 * h), comp1, comp2);
    p4 =   8.0 * (*func)( (x + h)      , comp1, comp2);
    denom = inv( 12.0 * h);
    deriv = (p1 + p2 + p3 + p4) * denom;
    return deriv;
}

static inline real second_derivative(real (*func)(real, const Dwarf*, const Dwarf*), real x, const Dwarf* comp1, const Dwarf* comp2)
{
    /*yes, this also uses a five point stencil*/
    real h = 0.001;
    real deriv;
    real p1, p2, p3, p4, p5, denom;

    p1 = - 1.0 * (*func)( (x + 2.0 * h) , comp1, comp2);
    p2 =  16.0 * (*func)( (x + h)       , comp1, comp2);
    p3 = -30.0 * (*func)( (x)           , comp1, comp2);
    p4 =  16.0 * (*func)( (x - h)       , comp1, comp2);
    p5 = - 1.0 * (*func)( (x - 2.0 * h) , comp1, comp2);
    denom = inv( 12.0 * h * h);
    deriv = (p1 + p2 + p3 + p4 + p5) * denom;
    return deriv;
}

static real gauss_quad(real (*func)(real, const Dwarf*, const Dwarf*, real), real lower, real upper, const Dwarf* comp1, const Dwarf* comp2, real energy)
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
    
    benchmark = 1.5 * a;
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
        intv = intv + c1 * (*func)(x1n, comp1, comp2, energy) * coef1 +
                      c2 * (*func)(x2n, comp1, comp2, energy) * coef1 + 
                      c3 * (*func)(x3n, comp1, comp2, energy) * coef1;

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
            Ng = 10.0;//integral resolution
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

static inline real max_finder(real (*profile)(real , real , const Dwarf*, const Dwarf*), real r, const Dwarf* comp1, const Dwarf* comp2, real a, real b, real c, int limit, real tolerance)
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

    profile_x1 = -(*profile)(x1, r, comp1, comp2);
    profile_x2 = -(*profile)(x2, r, comp1, comp2);
    
    while (mw_fabs(x3 - x0) > (tolerance * (mw_fabs(x1) + mw_fabs(x2)) ) )
    {
        counter++;
        if (profile_x2 < profile_x1)
        {
            x0 = x1;
            x1 = x2;
            x2 = RATIO * x2 + RATIO_COMPLEMENT * x3;
            profile_x1 = (real)profile_x2;
            profile_x2 = -(*profile)(x2, r, comp1, comp2);
        }
        else
        {
            x3 = x2;
            x2 = x1;
            x1 = RATIO * x1 + RATIO_COMPLEMENT * x0;
            profile_x2 = (real)profile_x1;
            profile_x1 = -(*profile)(x1, r, comp1, comp2);
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


static inline real root_finder(real (*func)(real, const Dwarf*, const Dwarf*), const Dwarf* comp1, const Dwarf* comp2, real function_value, real lower_bound, real upper_bound)
{
    //requires lower_bound and upper_bound to evaluate to opposite sign when func-function_value
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
        values[i] = (*func)(interval_bound, comp1, comp2) - function_value;
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
                mid_point_funcval = (*func)(mid_point, comp1, comp2) - function_value;
                
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
static real fun(real ri, const Dwarf* comp1, const Dwarf* comp2, real energy)
{
    
    real first_deriv_psi;
    real second_deriv_psi;
    real first_deriv_density;
    real second_deriv_density;
    real dsqden_dpsisq;/*second derivative of density with respect to -potential (psi) */
    real denominator; /*the demoninator of the distribution function: 1/sqrt(E-Psi)*/
    real diff;
    real func;

    first_deriv_psi      = first_derivative(potential, ri, comp1, comp2);
    first_deriv_density  = first_derivative(density,   ri, comp1, comp2);

    second_deriv_psi     = second_derivative(potential, ri, comp1, comp2);
    second_deriv_density = second_derivative(density,   ri, comp1, comp2);
    
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
    diff = mw_fabs(energy - potential(ri, comp1, comp2));
    
    /*we don't want to have a 0 in the demon*/
    if(diff != 0.0)
    {
        denominator = minushalf( diff );
    }
    else
    {
        /*if the r is exactly at the singularity then move it a small amount.*/
        denominator = minushalf( mw_fabs(energy - potential(ri + 0.0001, comp1, comp2) ) );
    }
    
    
    /*
     * the second derivative term should be divided by the first derivate of psi. 
     * However, from changing from dpsi to dr we multiply by first derivative of psi. 
     * Since these undo each other we left them out completely.
     */
    
    func = dsqden_dpsisq * denominator; 
    
    return func;
        
}

static inline real find_upperlimit_r(const Dwarf* comp1, const Dwarf* comp2, real energy, real search_range, real r)
{
    int counter = 0;
    real upperlimit_r = 0.0;

    do
    {
        upperlimit_r = root_finder(potential, comp1, comp2, energy, 0.0, search_range); 

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
 
static inline real dist_fun(real v, real r, const Dwarf* comp1, const Dwarf* comp2)
{
    /*This returns the value of the distribution function*/
    
    //-------------------------------
    real mass_l   = comp1->mass; //comp1[0]; /*mass of the light component*/
    real mass_d   = comp2->mass; //comp2[0]; /*mass of the dark component*/
    real rscale_l = comp1->scaleLength; //comp1[1]; /*scale radius of the light component*/
    real rscale_d = comp2->scaleLength; //comp2[1]; /*scale radius of the dark component*/
    //-------------------------------
    
    
    real distribution_function = 0.0;
//     real cons = inv( (mw_sqrt(8.0) * sqr(M_PI)) );
    real cons = 0.03582244801567226;
    real energy = 0.0;
    real upperlimit_r = 0.0;
    real lowerlimit_r = 0.0; 
    int counter = 0;
    real search_range = 0.0;   
    
    /*energy as defined in binney*/
    energy = potential(r, comp1, comp2) - 0.5 * v * v; 
    
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
    
    
    while(potential(search_range, comp1, comp2) > energy)
    {
        search_range = 100.0 * search_range;
        
        if(counter > 100)
        {
            search_range = 100.0 * (rscale_l + rscale_d);//default
            break;
        }
        counter++;
    }
    upperlimit_r = find_upperlimit_r(comp1, comp2, energy, search_range, r);
    /* This lowerlimit should be good enough. In the important case where the upperlimit is small (close to the singularity in the integrand)
     * then 5 times it is already where the integrand is close to 0 since it goes to 0 quickly. 
     */
    lowerlimit_r = 10.0 * (upperlimit_r);

    /*This calls guassian quad to integrate the function for a given energy*/
    distribution_function = v * v * cons * gauss_quad(fun, lowerlimit_r, upperlimit_r, comp1, comp2, energy);
    return distribution_function;
}

/*      SAMPLING FUNCTIONS      */
static inline real r_mag(dsfmt_t* dsfmtState, const Dwarf* comp, real rho_max, real bound)
{
    int counter = 0;
    real r, u, val;
    
    /*this technically calls the massless density but that is fine because
    * the masses would cancel in the denom since 
    * we are sampling the one component model.
    */
    
    /* the sampling is protected from r = 0. if profiles have a singularity there they would return inf or NANs
     * this would not satisfy the break conidition so it would choose another r.
     * if counter limit is reached r = 0 is returned which isn't accepted in the calling function so sampling is redone.
     */
    while (1)
    {
        r = (real)mwXrandom(dsfmtState, 0.0, 1.0) * bound;
        u = (real)mwXrandom(dsfmtState, 0.0, 1.0);
        val = r * r * get_density(comp, r);

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

static inline real vel_mag(real r, const Dwarf* comp1, const Dwarf* comp2, dsfmt_t* dsfmtState)
{
    
    /*
     * WE TOOK IN MASS IN SIMULATION UNITS, WHICH HAVE THE UNITS OF KPC^3/GY^2 
     * LENGTH IN KPC AND TIME IN GY THEREFORE, THE velocities ARE OUTPUTING IN KPC/GY
     * THIS IS EQUAL TO 0.977813107 KM/S
     */
    
    
    int counter = 0;
    real v, u, d;
    
    /* having the upper limit as exactly v_esc is bad since the dist fun seems to blow up there for small r. */
    real v_esc = 0.99 * mw_sqrt( mw_fabs(2.0 * potential( r, comp1, comp2) ) );
    
    real dist_max = max_finder(dist_fun, r, comp1, comp2, 0.0, 0.5 * v_esc, v_esc, 10, 1.0e-2);
    while(1)
    {

        v = (real)mwXrandom(dsfmtState, 0.0, 1.0) * v_esc;
        u = (real)mwXrandom(dsfmtState, 0.0, 1.0);

        d = dist_fun(v, r, comp1, comp2);
        
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
    v *= 0.977813107;//changing from kpc/gy to km/s
    return v; //km/s
}

static inline mwvector get_components(dsfmt_t* dsfmtState, real rad)
{
    /* assigns angles. Allows for non-circular orbits.*/
    /* have to sample in this way because sampling angles and then converting
     * to xyz leads to strong dependence on the rad, which could lead to bunching 
     * at the poles.
     */
    real r_sq, r_scaling;
    mwvector vec;

    do                                       
    {
        vec = mwRandomUnitPoint(dsfmtState); /* pick point in NDIM-space */
        r_sq = mw_sqrv(vec);                 /* compute radius squared */
    }
    while (r_sq > 1.0);                      /* reject if outside unit sphere */

    r_scaling = rad / mw_sqrt(r_sq);         /* compute scaling factor */
    mw_incmulvs(vec, r_scaling);             /* rescale to radius given */
    
    /* this is r * (u_vec / |u|). 
     * the r gives the magnitude, rad.
     * u_vec, which is the unit point original picked, vec, 
     * divided by the magnitude |u|, which is sqrt(r_sq),
     * gives it a direction (unit vector).
     */
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


static inline void set_p0(Dwarf* comp)
{
    /*this is only used for the nfw but it is technically valid for all the profiles. easier to have it here*/
    /* this is the pcrit * delta_crit from the nfw 1997 paper or just p0 from binney */
    //as defined in Binney and Tremaine 2nd ed:
    //the r200 is now used for all potentials to provide the bounds for density sampling
    real mass = comp->mass; 
    real rscale = comp->scaleLength;
    real r200 = mw_cbrt(mass / (vol_pcrit));//vol_pcrit = 200.0 * pcrit * PI_4_3
    real c = r200 / rscale; //halo concentration
    real term = mw_log(1.0 + c) - c / (1.0 + c);
    real p0 = inv(1.0) * 200.0 * cube(c) * pcrit / (3.0 * term); //rho_0 as defined in Navarro et. al. 1997
    comp->r200 = r200;
    comp->p0 = p0;
}


static inline void get_extra_nfw_mass(Dwarf* comp, real bound)
{
    real rs = comp->scaleLength;
    real r = bound;
    real m = 4.0 * M_PI * comp->p0 * cube(rs) * (mw_log( (rs + r) / rs) - r / (rs + r));
    comp->mass = m;
    
}


/*      DWARF GENERATION        */
static int nbGenerateMixedDwarfCore(lua_State* luaSt, dsfmt_t* prng, unsigned int nbody, 
                                     Dwarf* comp1,  Dwarf* comp2, 
                                    mwbool ignore, mwvector rShift, mwvector vShift)
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
        real rscale_l = comp1->scaleLength; //comp1[1]; /*scale radius of the light component*/
        real rscale_d = comp2->scaleLength; //comp2[1]; /*scale radius of the dark component*/
        set_p0(comp1);
        set_p0(comp2);
        real bound1 ;
        real bound2 ;
        
        /*
         * NOTE: the nfw technically works but it does not output in virial equilibrium. 
         * will return to it at a later date. The Plummer and hern work.
         */
        
        switch(comp1->type)
        {
            case Plummer:
                bound1 =  50.0 * (rscale_l + rscale_d);
                break;
            case NFW:
                bound1 = 5 * comp1->r200;
                get_extra_nfw_mass(comp1, bound1);
                break;
            case General_Hernquist:
                bound1 =  50.0 * (rscale_l + rscale_d);
                break;
        }

        switch(comp2->type)
        {
            case Plummer:
                bound2 =  50.0 * (rscale_l + rscale_d);
                break;
            case NFW:
                bound2 = 5 * comp2->r200;
                get_extra_nfw_mass(comp2, bound2);
                break;
            case General_Hernquist:
                bound2 =  50.0 * (rscale_l + rscale_d);
                break;
        }
        
        
        real mass_l   = comp1->mass; //comp1[0]; /*mass of the light component*/
        real mass_d   = comp2->mass; //comp2[0]; /*mass of the dark component*/
        real dwarf_mass = mass_l + mass_d;

    //---------------------------------------------------------------------------------------------------        
        unsigned int half_bodies = nbody / 2;
        real mass_light_particle = mass_l / (real)(0.5 * (real) nbody);//half the particles are light matter
        real mass_dark_particle = mass_d / (real)(0.5 * (real) nbody);
    //----------------------------------------------------------------------------------------------------

        /* dark matter type is TRUE or 1. Light matter type is False, or 0*/
        mwbool isdark = TRUE;
        mwbool islight = FALSE;
       
        
        /*finding the max of the individual components*/
        int place_holder = 0;
        real rho_max_light = max_finder(profile_rho, place_holder, comp1, comp1, 0, rscale_l, 2.0 * (rscale_l), 20, 1e-4 );
        real rho_max_dark  = max_finder(profile_rho, place_holder, comp2, comp2, 0, rscale_d, 2.0 * (rscale_d), 20, 1e-4 );

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
                    r = r_mag(prng, comp1, rho_max_light, bound1);
                    masses[i] = mass_light_particle;
                }
                else if(i >= half_bodies)
                {
                    r = r_mag(prng, comp2, rho_max_dark, bound2);
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
            
//             mw_printf("\rvelocity of particle %i", i + 1);
            counter = 0;
            do
            {
                v = vel_mag(r, comp1, comp2, prng);
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

int nbGenerateMixedDwarf(lua_State* luaSt)
{
        static dsfmt_t* prng;
        static const mwvector* position = NULL;
        static const mwvector* velocity = NULL;
        static mwbool ignore;
        static real nbodyf = 0.0;
        static Dwarf* comp1 = NULL;
        static Dwarf* comp2 = NULL;
        static const MWNamedArg argTable[] =
        {
            { "nbody",                LUA_TNUMBER,     NULL,                    TRUE,    &nbodyf            },
            { "comp1",                LUA_TUSERDATA,   DWARF_TYPE,              TRUE,    &comp1             },
            { "comp2",                LUA_TUSERDATA,   DWARF_TYPE,              TRUE,    &comp2             },
            { "position",             LUA_TUSERDATA,   MWVECTOR_TYPE,           TRUE,    &position          },
            { "velocity",             LUA_TUSERDATA,   MWVECTOR_TYPE,           TRUE,    &velocity          },
            { "ignore",               LUA_TBOOLEAN,    NULL,                    FALSE,   &ignore            },
            { "prng",                 LUA_TUSERDATA,   DSFMT_TYPE,              TRUE,    &prng              },
            END_MW_NAMED_ARG
            
        };

        if (lua_gettop(luaSt) != 1)
            return luaL_argerror(luaSt, 1, "Expected 1 arguments");
        
        handleNamedArgumentTable(luaSt, argTable, 1);
        
        return nbGenerateMixedDwarfCore(luaSt, prng, (unsigned int) nbodyf, comp1, comp2, ignore,
                                                                 *position, *velocity);
}


void registerGenerateMixedDwarf(lua_State* luaSt)
{
    lua_register(luaSt, "generatemixeddwarf", nbGenerateMixedDwarf);
}


// As this code runs, know that it is running on the rotting corpses of a thousand bugs.