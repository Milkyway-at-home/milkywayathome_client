/* Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
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
#include "nbody_dwarf_potential.h"
#include "nbody_mixeddwarf.h"
#include "nbody_types.h"
#include "nbody_potential_types.h"

#include "nbody_autodiff.h"
#include <stdio.h>
#include <time.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

/*Note: minusfivehalves(x) raises to x^-5/2 power and minushalf(x) is x^-1/2*/

#define H_STEPSIZE ((real_0) 0.0001) //To approximate initial derivatives in AUTODIFF

/*      MODEL SPECIFIC FUNCTIONS       */
static inline real_0 potential( real_0 r, const Dwarf* comp1, const Dwarf* comp2)
{
    /*Be Careful! this function returns the negative of the potential! this is the value of interest, psi*/
    real_0 potential_light  = get_potential(comp1, r);
    real_0 potential_dark   = get_potential(comp2, r);
    real_0 potential_result = (potential_light + potential_dark);

    return (potential_result);
}

static inline real_0 density( real_0 r, const Dwarf* comp1, const Dwarf* comp2)
{
    /*this is the density distribution function. Returns the density at a given radius.*/
    
    real_0 density_light = get_density(comp1, r);
    real_0 density_dark  = get_density(comp2, r);
    real_0 density_result = (density_light + density_dark );

    return density_result;
}


/*      GENERAL PURPOSE DERIVATIVE, INTEGRATION, MAX FINDING, ROOT FINDING, AND ARRAY SHUFFLER FUNCTIONS        */
static inline real_0 first_derivative(real_0 (*func)(real_0, const Dwarf*, const Dwarf*), real_0 x, const Dwarf* comp1, const Dwarf* comp2)
{
    /*yes, this does in fact use a 5-point stencil*/
    real_0 h = 0.001;
    real_0 deriv;
    real_0 p1, p2, p3, p4, denom;
    
    p1 =   1.0 * (*func)( (x - 2.0 * h), comp1, comp2);
    p2 = - 8.0 * (*func)( (x - h)      , comp1, comp2);
    p3 = - 1.0 * (*func)( (x + 2.0 * h), comp1, comp2);
    p4 =   8.0 * (*func)( (x + h)      , comp1, comp2);
    denom = inv_0( 12.0 * h);
    deriv = (p1 + p2 + p3 + p4) * denom;
    return deriv;
}

static inline real_0 second_derivative(real_0 (*func)(real_0, const Dwarf*, const Dwarf*), real_0 x, const Dwarf* comp1, const Dwarf* comp2)
{
    /*yes, this also uses a five point stencil*/
    real_0 h = 0.001;
    real_0 deriv;
    real_0 p1, p2, p3, p4, p5, denom;

    p1 = - 1.0 * (*func)( (x + 2.0 * h) , comp1, comp2);
    p2 =  16.0 * (*func)( (x + h)       , comp1, comp2);
    p3 = -30.0 * (*func)( (x)           , comp1, comp2);
    p4 =  16.0 * (*func)( (x - h)       , comp1, comp2);
    p5 = - 1.0 * (*func)( (x - 2.0 * h) , comp1, comp2);
    denom = inv_0( 12.0 * h * h);
    deriv = (p1 + p2 + p3 + p4 + p5) * denom;
    return deriv;
}

static real_0 gauss_quad(real_0 (*func)(real_0, const Dwarf*, const Dwarf*, real_0), real_0 lower, real_0 upper, const Dwarf* comp1, const Dwarf* comp2, real_0 energy)
{
    /*This is a guassian quadrature routine. It will test to always integrate from the lower to higher of the two limits.
     * If switching the order of the limits was needed to do this then the negative of the integral is returned.
     */
    real_0 Ng, hg, lowerg, upperg;
    real_0 intv = 0.0;//initial value of integral
    real_0 coef1, coef2;//parameters for gaussian quad
    real_0 c1, c2, c3;
    real_0 x1, x2, x3;
    real_0 x1n, x2n, x3n;
    real_0 a, b;
    real_0 benchmark;
    
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

static inline real_0 max_finder(real_0 (*profile)(real_0 , real_0 , const Dwarf*, const Dwarf*), real_0 r, const Dwarf* comp1, const Dwarf* comp2, real_0 a, real_0 b, real_0 c, int limit, real_0 tolerance)
{
    /*this is a maxfinding routine to find the maximum of the density.
     * It uses Golden Section Search as outlined in Numerical Recipes 3rd edition
     */
    real_0 RATIO = 0.61803399;
    real_0 RATIO_COMPLEMENT = 1.0 - RATIO;
    int counter = 0;
    
    real_0 profile_x1, profile_x2, x0, x1, x2, x3;
    x0 = a;
    x3 = c;
    
    if (mw_fabs_0(b - c) > mw_fabs_0(b - a))
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
    
    while (mw_fabs_0(x3 - x0) > (tolerance * (mw_fabs_0(x1) + mw_fabs_0(x2)) ) )
    {
        counter++;
        if (profile_x2 < profile_x1)
        {
            x0 = x1;
            x1 = x2;
            x2 = RATIO * x2 + RATIO_COMPLEMENT * x3;
            profile_x1 = (real_0)profile_x2;
            profile_x2 = -(*profile)(x2, r, comp1, comp2);
        }
        else
        {
            x3 = x2;
            x2 = x1;
            x1 = RATIO * x1 + RATIO_COMPLEMENT * x0;
            profile_x2 = (real_0)profile_x1;
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


static inline real_0 root_finder(real_0 (*func)(real_0, const Dwarf*, const Dwarf*), const Dwarf* comp1, const Dwarf* comp2, real_0 function_value, real_0 lower_bound, real_0 upper_bound)
{
    //requires lower_bound and upper_bound to evaluate to opposite sign when func-function_value
    unsigned int i = 0;

    int N = 4;
    unsigned int intervals = N;
    real_0 interval_bound;

    /*interval + 1 because for N intervals there are N + 1 values*/
    real_0 * values = mwCalloc(intervals + 1, sizeof(real_0));
    real_0 * interval_bounds = mwCalloc(intervals + 1, sizeof(real_0));
    /*intervals+1 because you want to include the upperbound in the interval*/
    for(i = 0; i < intervals + 1; i++)
    {
        /* breaking up the range between bounds into smaller intervals*/
        interval_bound = ((upper_bound - lower_bound) * (real_0)i) / (real_0)intervals + lower_bound;
        /*function value at those intervals*/
        values[i] = (*func)(interval_bound, comp1, comp2) - function_value;
	if(isnan(values[i]))
	{	
		/*If the interval bound is at the singularity, shift it slightly to prevent a nan so that the root can still be found;
		* added specifically to prevent issues with finding the root for an NFW profile*/
		interval_bound += 0.0001;
		values[i] = (*func)(interval_bound, comp1, comp2) - function_value;
	}
	interval_bounds[i] = interval_bound;
    }
    
    real_0 mid_point = 0;
    real_0 mid_point_funcval = 0;
    unsigned int counter = 0;
    real_0 new_upper_bound = 0;
    real_0 new_lower_bound = 0;
    int roots_found = 0;
    int q = 0;
    
    /* Find the roots using bisection because it was easy to code and good enough for our purposes 
     * this will hop around the different intervals until it checks all of them. This way it does not 
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
            while(mw_fabs_0(mid_point_funcval) > .0001)
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
static real_0 fun(real_0 ri, const Dwarf* comp1, const Dwarf* comp2, real_0 energy)
{
    
    real_0 first_deriv_psi;
    real_0 second_deriv_psi;
    real_0 first_deriv_density;
    real_0 second_deriv_density;
    real_0 dsqden_dpsisq;/*second derivative of density with respect to -potential (psi) */
    real_0 denominator; /*the demoninator of the distribution function: 1/sqrt(E-Psi)*/
    real_0 diff;
    real_0 func;

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
    
    dsqden_dpsisq = second_deriv_density * inv_0(first_deriv_psi) - first_deriv_density * second_deriv_psi * inv_0(sqr_0(first_deriv_psi));
    diff = mw_fabs_0(energy - potential(ri, comp1, comp2));
    
    /*we don't want to have a 0 in the demon*/
    if(diff != 0.0)
    {
        denominator = minushalf_0( diff );
    }
    else
    {
        /*if the r is exactly at the singularity then move it a small amount.*/
        denominator = minushalf_0( mw_fabs_0(energy - potential(ri + 0.0001, comp1, comp2) ) );
    }
    
    
    /*
     * the second derivative term should be divided by the first derivate of psi. 
     * However, from changing from dpsi to dr we multiply by first derivative of psi. 
     * Since these undo each other we left them out completely.
     */
    
    func = dsqden_dpsisq * denominator;
    //printf("radius: %1f, energy: %1f, numerator: %1f, denom: %1f, func: %1f\n", ri, energy, dsqden_dpsisq, 1.0 / denominator, func);
    
    return func;
        
}

static inline real_0 find_upperlimit_r(const Dwarf* comp1, const Dwarf* comp2, real_0 energy, real_0 search_range, real_0 r)
{
    int counter = 0;
    real_0 upperlimit_r = 0.0;

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
        
    return mw_fabs_0(upperlimit_r);
}
 
static inline real_0 dist_fun(real_0 v, real_0 r, const Dwarf* comp1, const Dwarf* comp2)
{
    /*This returns the value of the distribution function*/
    
    //-------------------------------
    real_0 mass_l   = comp1->mass; //comp1[0]; /*mass of the light component*/
    real_0 mass_d   = comp2->mass; //comp2[0]; /*mass of the dark component*/
    real_0 rscale_l = comp1->scaleLength; //comp1[1]; /*scale radius of the light component*/
    real_0 rscale_d = comp2->scaleLength; //comp2[1]; /*scale radius of the dark component*/
    //-------------------------------
    
    
    real_0 distribution_function = 0.0;
//     real_0 cons = inv( (mw_sqrt(8.0) * sqr(M_PI)) );
    real_0 cons = 0.03582244801567226;
    real_0 energy = 0.0;
    real_0 upperlimit_r = 0.0;
    real_0 lowerlimit_r = 0.0; 
    int counter = 0;
    real_0 search_range = 0.0;   
    
    /*energy as defined in binney*/
    energy = potential(r, comp1, comp2) - 0.5 * v * v; 
    
    /*this starting point is 20 times where the dark matter component is equal to the energy, since the dark matter dominates*/
    search_range = 20.0 * mw_sqrt_0( mw_fabs_0( sqr_0(mass_d / energy) - sqr_0(rscale_d) ));
    
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

static inline real_0 radial_dist(real_0 r, const Dwarf* comp1, const Dwarf* comp2, mwbool isLight) //Normalized radial distribution
{
    if(isLight)
    {
        return 4 * M_PI * r * r * get_density(comp1, r) / comp1->mass;
    }
    else
    {
        return 4 * M_PI * r * r * get_density(comp2, r) / comp2->mass;
    }
}

static inline real_0 vel_dist(real_0 v, real_0 r, const Dwarf* comp1, const Dwarf* comp2) //Normalized velocity distribution
{
    //Calculate normalization for velocity distribution function
    real_0 a = 0.0;
    real_0 b = 0.99 * mw_sqrt_0( mw_fabs_0(2.0 * potential( r, comp1, comp2) ) ); //v_esc
    real_0 numBins = 10.0;
    real_0 width = (b-a)/numBins;
    real_0 a_tmp;
    real_0 b_tmp;
    real_0 c1 = 0.55555555555; //5.0 / 9.0;
    real_0 c2 = 0.88888888888; //8.0 / 9.0;
    real_0 c3 = 0.55555555555; //5.0 / 9.0;
    real_0 x1 = -0.77459666924;//-sqrt(3.0 / 5.0);
    real_0 x2 = 0.00000000000;
    real_0 x3 = 0.77459666924; //sqrt(3.0 / 5.0);
    real_0 v1, v2, v3, integrand1, integrand2, integrand3;

    real_0 integral = 0;
    for(int k = 0; k < numBins; k++)
    {
        a_tmp = k*width;
        b_tmp = (k+1)*width;
        v1 = (b_tmp-a_tmp)*x1/2.0 + (b_tmp+a_tmp)/2.0;
        v2 = (b_tmp-a_tmp)*x2/2.0 + (b_tmp+a_tmp)/2.0;
        v3 = (b_tmp-a_tmp)*x3/2.0 + (b_tmp+a_tmp)/2.0;

        integrand1 = dist_fun(v1, r, comp1, comp2);
        integrand2 = dist_fun(v2, r, comp1, comp2);
        integrand3 = dist_fun(v3, r, comp1, comp2);

        integral += (b_tmp-a_tmp)/2.0*c1*( integrand1 )
                  + (b_tmp-a_tmp)/2.0*c2*( integrand2 )
                  + (b_tmp-a_tmp)/2.0*c3*( integrand3 );
    }

    //Report normalized distribution value
    return dist_fun(v, r, comp1, comp2) / integral;    
}

/* These functions prep the Dwarf structs for rejection sampling */
static inline void set_p0(Dwarf* comp)
{
    /*this is only used for the nfw but it is technically valid for all the profiles. easier to have it here*/
    /* this is the pcrit * delta_crit from the nfw 1997 paper or just p0 from binney */
    //as defined in Binney and Tremaine 2nd ed:
    //the r200 is now used for all potentials to provide the bounds for density sampling
    real_0 mass = comp->mass;
    comp->originmass = comp->mass; //need to store original mass value for nfw AUTODIFF
    real_0 rscale = comp->scaleLength;
    real_0 r200 = mw_cbrt_0(mass / (vol_pcrit));//vol_pcrit = 200.0 * pcrit * PI_4_3
    real_0 c = r200 / rscale; //halo concentration
    real_0 term = mw_log_0(1.0 + c) - c / (1.0 + c);
    real_0 p0 = 200.0 * cube_0(c) * pcrit / (3.0 * term); //rho_0 as defined in Navarro et. al. 1997
    comp->r200 = r200;
    comp->p0 = p0;
}

static inline void get_extra_nfw_mass(Dwarf* comp, real_0 bound)
{
    /* The mass inputted is taken to be the M200 mass (mass within radius r200).*/
    /* If the sampling boundary goes above or below r200, this function resets the mass of the component.*/
    real_0 rs = comp->scaleLength;
    real_0 r = bound;
    real_0 m = 4.0 * M_PI * comp->p0 * cube_0(rs) * (mw_log_0( (rs + r) / rs) - r / (rs + r));
    comp->mass = m;
}

/* AUTODIFF METHODS FOR CALCULATING INITIAL DERIVATIVES */
#if AUTODIFF
  static inline Dwarf getShiftedDwarf(const Dwarf* comp, real_0 a_b, real_0 xi_R, real_0 M_b, real_0 xi_M, mwbool isLight)
  {
      Dwarf shiftedDwarf = *comp;
      real_0 a_d = (1.0/xi_R - 1.0)*a_b;
      real_0 M_d = (1.0/xi_M - 1.0)*M_b;
      if(isLight)
      {
          shiftedDwarf.scaleLength = a_b;
          shiftedDwarf.mass = M_b;
      }
      else
      {
          shiftedDwarf.scaleLength = a_d;
          shiftedDwarf.mass = M_d;
      }
      set_p0(&shiftedDwarf);
      switch(shiftedDwarf.type)
      {
          case NFW:
              get_extra_nfw_mass(&shiftedDwarf, 5.0 * shiftedDwarf.r200);
      }
      return shiftedDwarf;
  }

  /*---------------------------------------------------METHODS TO CALCULATE DERIVATIVES FOR RADIUS---------------------------------------------------*/
  static inline real_0 gradIntegrand_rad(real_0 r, const Dwarf* comp1, const Dwarf* comp2, mwbool isLight, const Dwarf* shiftedDwarf1, const Dwarf* shiftedDwarf2)
  {
      return (radial_dist(r, shiftedDwarf1, shiftedDwarf2, isLight) - radial_dist(r, comp1, comp2, isLight))/H_STEPSIZE;
  }

  static inline real_0 getGradIntegral_rad(real* r, const Dwarf* comp1, const Dwarf* comp2, mwbool isLight, int i)
  {
      real_0 scale_b = comp1->scaleLength;
      real_0 mass_b = comp1->originmass;
      real_0 xi_scale = comp1->scaleLength / (comp1->scaleLength + comp2->scaleLength);
      real_0 xi_mass = comp1->originmass / (comp1->originmass + comp2->originmass);
      Dwarf shiftedDwarf1, shiftedDwarf2;
      real_0 numBins = 10.0;
      real_0 integral;

      real_0 a = 0.0;
      real_0 b = r->value;
      real_0 width = (b-a)/numBins;
      real_0 a_tmp;
      real_0 b_tmp;
      real_0 c1 = 0.55555555555; //5.0 / 9.0;
      real_0 c2 = 0.88888888888; //8.0 / 9.0;
      real_0 c3 = 0.55555555555; //5.0 / 9.0;
      real_0 x1 = -0.77459666924;//-sqrt(3.0 / 5.0);
      real_0 x2 = 0.00000000000;
      real_0 x3 = 0.77459666924; //sqrt(3.0 / 5.0);
      real_0 r1, r2, r3, integrand1, integrand2, integrand3;

      switch(i)
      {
          case(BARYON_RADIUS_POS):
              scale_b += H_STEPSIZE;
          case(RADIUS_RATIO_POS):
              xi_scale += H_STEPSIZE;
          case(BARYON_MASS_POS):
              mass_b += H_STEPSIZE;
          case(MASS_RATIO_POS):
              xi_mass += H_STEPSIZE;
      }

      shiftedDwarf1 = getShiftedDwarf(comp1, scale_b, xi_scale, mass_b, xi_mass, TRUE);
      shiftedDwarf2 = getShiftedDwarf(comp2, scale_b, xi_scale, mass_b, xi_mass, FALSE);
      integral = 0;
      for(int k = 0; k < numBins; k++)
      {
          a_tmp = k*width;
          b_tmp = (k+1)*width;
          r1 = (b_tmp-a_tmp)*x1/2.0 + (b_tmp+a_tmp)/2.0;
          r2 = (b_tmp-a_tmp)*x2/2.0 + (b_tmp+a_tmp)/2.0;
          r3 = (b_tmp-a_tmp)*x3/2.0 + (b_tmp+a_tmp)/2.0;

          integrand1 = gradIntegrand_rad(r1, comp1, comp2, isLight, &shiftedDwarf1, &shiftedDwarf2);
          integrand2 = gradIntegrand_rad(r2, comp1, comp2, isLight, &shiftedDwarf1, &shiftedDwarf2);
          integrand3 = gradIntegrand_rad(r3, comp1, comp2, isLight, &shiftedDwarf1, &shiftedDwarf2);

          integral += (b_tmp-a_tmp)/2.0*c1*( integrand1 )
                    + (b_tmp-a_tmp)/2.0*c2*( integrand2 )
                    + (b_tmp-a_tmp)/2.0*c3*( integrand3 );
      }
      integral = -integral / radial_dist(r->value, comp1, comp2, isLight);
      return integral;
  }

  static inline void getGradientInfo_rad(real* r, const Dwarf* comp1, const Dwarf* comp2, mwbool isLight)
  {
      real_0 integral;

      //a_b derivative
      integral = getGradIntegral_rad(r, comp1, comp2, isLight, BARYON_RADIUS_POS);
      setRealGradient(r, integral, BARYON_RADIUS_POS);

      //xi_R derivative
      integral = getGradIntegral_rad(r, comp1, comp2, isLight, RADIUS_RATIO_POS);
      setRealGradient(r, integral, RADIUS_RATIO_POS);

      //M_b derivative
      integral = getGradIntegral_rad(r, comp1, comp2, isLight, BARYON_MASS_POS);
      setRealGradient(r, integral, BARYON_MASS_POS);

      //xi_M derivative
      integral = getGradIntegral_rad(r, comp1, comp2, isLight, MASS_RATIO_POS);
      setRealGradient(r, integral, MASS_RATIO_POS);
  }

  static inline real_0 hessIntegrand_rad(real_0 r, const Dwarf* comp1, const Dwarf* comp2, mwbool isLight, const Dwarf* shiftedDwarf1_i, const Dwarf* shiftedDwarf1_j, const Dwarf* shiftedDwarf1_d, const Dwarf* shiftedDwarf2_i, const Dwarf* shiftedDwarf2_j, const Dwarf* shiftedDwarf2_d, int i, int j)
  {
      real_0 integrand = (  radial_dist(r, shiftedDwarf1_d, shiftedDwarf2_d, isLight)
                          - radial_dist(r, shiftedDwarf1_i, shiftedDwarf2_i, isLight)
                          - radial_dist(r, shiftedDwarf1_j, shiftedDwarf2_j, isLight)
                          + radial_dist(r, comp1, comp2, isLight)                    )/H_STEPSIZE/H_STEPSIZE;

      return integrand;
  }

  static inline real_0 getHessIntegral_rad(real* r, const Dwarf* comp1, const Dwarf* comp2, mwbool isLight, int i, int j)
  {
      real_0 scale_b = comp1->scaleLength;
      real_0 mass_b = comp1->originmass;
      real_0 xi_scale = comp1->scaleLength / (comp1->scaleLength + comp2->scaleLength);
      real_0 xi_mass = comp1->originmass / (comp1->originmass + comp2->originmass);
      Dwarf shiftedDwarf1_i, shiftedDwarf1_j, shiftedDwarf1_d, shiftedDwarf2_i, shiftedDwarf2_j, shiftedDwarf2_d;
      real_0 integral;
      real_0 numBins = 10.0;

      real_0 a = 0.0;
      real_0 b = r->value;
      real_0 width = (b-a)/numBins;
      real_0 a_tmp;
      real_0 b_tmp;
      real_0 c1 = 0.55555555555; //5.0 / 9.0;
      real_0 c2 = 0.88888888888; //8.0 / 9.0;
      real_0 c3 = 0.55555555555; //5.0 / 9.0;
      real_0 x1 = -0.77459666924;//-sqrt(3.0 / 5.0);
      real_0 x2 = 0.00000000000;
      real_0 x3 = 0.77459666924; //sqrt(3.0 / 5.0);
      real_0 r1, r2, r3, integrand1, integrand2, integrand3;

      real_0 scale_b_i = scale_b;
      real_0 scale_b_j = scale_b;
      real_0 scale_b_d = scale_b;

      real_0 mass_b_i = mass_b;
      real_0 mass_b_j = mass_b;
      real_0 mass_b_d = mass_b;

      real_0 xi_scale_i = xi_scale;
      real_0 xi_scale_j = xi_scale;
      real_0 xi_scale_d = xi_scale;

      real_0 xi_mass_i = xi_mass;
      real_0 xi_mass_j = xi_mass;
      real_0 xi_mass_d = xi_mass;

      switch(i)
      {
          case(BARYON_RADIUS_POS):
              scale_b_i += H_STEPSIZE;
              scale_b_d += H_STEPSIZE;
          case(RADIUS_RATIO_POS):
              xi_scale_i += H_STEPSIZE;
              xi_scale_d += H_STEPSIZE;
          case(BARYON_MASS_POS):
              mass_b_i += H_STEPSIZE;
              mass_b_d += H_STEPSIZE;
          case(MASS_RATIO_POS):
              xi_mass_i += H_STEPSIZE;
              xi_mass_d += H_STEPSIZE;
      }

      switch(j)
      {
          case(BARYON_RADIUS_POS):
              scale_b_j += H_STEPSIZE;
              scale_b_d += H_STEPSIZE;
          case(RADIUS_RATIO_POS):
              xi_scale_j += H_STEPSIZE;
              xi_scale_d += H_STEPSIZE;
          case(BARYON_MASS_POS):
              mass_b_j += H_STEPSIZE;
              mass_b_d += H_STEPSIZE;
          case(MASS_RATIO_POS):
              xi_mass_j += H_STEPSIZE;
              xi_mass_d += H_STEPSIZE;
      }

      shiftedDwarf1_i = getShiftedDwarf(comp1, scale_b_i, xi_scale_i, mass_b_i, xi_mass_i, TRUE);
      shiftedDwarf1_j = getShiftedDwarf(comp1, scale_b_j, xi_scale_j, mass_b_j, xi_mass_j, TRUE);
      shiftedDwarf1_d = getShiftedDwarf(comp1, scale_b_d, xi_scale_d, mass_b_d, xi_mass_d, TRUE);
      shiftedDwarf2_i = getShiftedDwarf(comp2, scale_b_i, xi_scale_i, mass_b_i, xi_mass_i, FALSE);
      shiftedDwarf2_j = getShiftedDwarf(comp2, scale_b_j, xi_scale_j, mass_b_j, xi_mass_j, FALSE);
      shiftedDwarf2_d = getShiftedDwarf(comp2, scale_b_d, xi_scale_d, mass_b_d, xi_mass_d, FALSE);

      real_0 df_dx_i = gradIntegrand_rad(r->value, comp1, comp2, isLight, &shiftedDwarf1_i, &shiftedDwarf2_i);
      real_0 df_dx_j = gradIntegrand_rad(r->value, comp1, comp2, isLight, &shiftedDwarf1_j, &shiftedDwarf2_j);
      real_0 df_de = (radial_dist(r->value + H_STEPSIZE, comp1, comp2, isLight) - radial_dist(r->value, comp1, comp2, isLight))/H_STEPSIZE;

      integral = 0.0;
      for(int k = 0; k < numBins; k++)
      {
          a_tmp = k*width;
          b_tmp = (k+1)*width;
          r1 = (b_tmp-a_tmp)*x1/2.0 + (b_tmp+a_tmp)/2.0;
          r2 = (b_tmp-a_tmp)*x2/2.0 + (b_tmp+a_tmp)/2.0;
          r3 = (b_tmp-a_tmp)*x3/2.0 + (b_tmp+a_tmp)/2.0;

          integrand1 = hessIntegrand_rad(r1, comp1, comp2, isLight, &shiftedDwarf1_i, &shiftedDwarf1_j, &shiftedDwarf1_d, &shiftedDwarf2_i, &shiftedDwarf2_j, &shiftedDwarf2_d, i, j);
          integrand2 = hessIntegrand_rad(r2, comp1, comp2, isLight, &shiftedDwarf1_i, &shiftedDwarf1_j, &shiftedDwarf1_d, &shiftedDwarf2_i, &shiftedDwarf2_j, &shiftedDwarf2_d, i, j);
          integrand3 = hessIntegrand_rad(r3, comp1, comp2, isLight, &shiftedDwarf1_i, &shiftedDwarf1_j, &shiftedDwarf1_d, &shiftedDwarf2_i, &shiftedDwarf2_j, &shiftedDwarf2_d, i, j);

          integral += (b_tmp-a_tmp)/2.0*c1*( integrand1 )
                    + (b_tmp-a_tmp)/2.0*c2*( integrand2 )
                    + (b_tmp-a_tmp)/2.0*c3*( integrand3 );
      }

      integral += df_dx_i * r->gradient[j] + df_dx_j * r->gradient[i] + df_de * r->gradient[i] * r->gradient[j]; // Other terms of the Hessian
      integral = -integral / radial_dist(r->value, comp1, comp2, isLight);
      return integral;
  }

  static inline void getHessianInfo_rad(real* r, const Dwarf* comp1, const Dwarf* comp2, mwbool isLight)
  {
      real_0 integral;

      //hessian (a_b, a_b)
      integral = getHessIntegral_rad(r, comp1, comp2, isLight, BARYON_RADIUS_POS, BARYON_RADIUS_POS);
      setRealHessian(r, integral, BARYON_RADIUS_POS, BARYON_RADIUS_POS);

      //hessian (a_b, xi_R)
      integral = getHessIntegral_rad(r, comp1, comp2, isLight, BARYON_RADIUS_POS, RADIUS_RATIO_POS);
      setRealHessian(r, integral, BARYON_RADIUS_POS, RADIUS_RATIO_POS);

      //hessian (a_b, M_b)
      integral = getHessIntegral_rad(r, comp1, comp2, isLight, BARYON_RADIUS_POS, BARYON_MASS_POS);
      setRealHessian(r, integral, BARYON_RADIUS_POS, BARYON_MASS_POS);

      //hessian (a_b, xi_M)
      integral = getHessIntegral_rad(r, comp1, comp2, isLight, BARYON_RADIUS_POS, MASS_RATIO_POS);
      setRealHessian(r, integral, BARYON_RADIUS_POS, MASS_RATIO_POS);

      //hessian (xi_R, xi_R)
      integral = getHessIntegral_rad(r, comp1, comp2, isLight, RADIUS_RATIO_POS, RADIUS_RATIO_POS);
      setRealHessian(r, integral, RADIUS_RATIO_POS, RADIUS_RATIO_POS);

      //hessian (xi_R, M_b)
      integral = getHessIntegral_rad(r, comp1, comp2, isLight, RADIUS_RATIO_POS, BARYON_MASS_POS);
      setRealHessian(r, integral, RADIUS_RATIO_POS, BARYON_MASS_POS);

      //hessian (xi_R, xi_M)
      integral = getHessIntegral_rad(r, comp1, comp2, isLight, RADIUS_RATIO_POS, MASS_RATIO_POS);
      setRealHessian(r, integral, RADIUS_RATIO_POS, MASS_RATIO_POS);

      //hessian (M_b, M_b)
      integral = getHessIntegral_rad(r, comp1, comp2, isLight, BARYON_MASS_POS, BARYON_MASS_POS);
      setRealHessian(r, integral, BARYON_MASS_POS, BARYON_MASS_POS);

      //hessian (M_b, xi_M)
      integral = getHessIntegral_rad(r, comp1, comp2, isLight, BARYON_MASS_POS, MASS_RATIO_POS);
      setRealHessian(r, integral, BARYON_MASS_POS, MASS_RATIO_POS);

      //hessian (xi_M, xi_M)
      integral = getHessIntegral_rad(r, comp1, comp2, isLight, MASS_RATIO_POS, MASS_RATIO_POS);
      setRealHessian(r, integral, MASS_RATIO_POS, MASS_RATIO_POS);      
  }

  static inline void getDwarfDerivativeInfo_rad(real* r, const Dwarf* comp1, const Dwarf* comp2, mwbool isLight)
  {
      //Calculate gradient
      getGradientInfo_rad(r, comp1, comp2, isLight);

      //Calculate hessian
      getHessianInfo_rad(r, comp1, comp2, isLight);
  }

  /*---------------------------------------------------METHODS TO CALCULATE DERIVATIVES FOR VELOCITY---------------------------------------------------*/
  static inline real_0 gradIntegrand_vel(real_0 v, real* r, const Dwarf* comp1, const Dwarf* comp2, const Dwarf* shiftedDwarf1, const Dwarf* shiftedDwarf2, int i)
  {
      real_0 shifted_r = r->value + r->gradient[i]*H_STEPSIZE; //r has derivative information, so we also apply a shift to r
      return (vel_dist(v, shifted_r, shiftedDwarf1, shiftedDwarf2) - vel_dist(v, r->value, comp1, comp2))/H_STEPSIZE;
  }

  static inline real_0 getGradIntegral_vel(real* v, real* r, const Dwarf* comp1, const Dwarf* comp2, int i)
  {
      real_0 scale_b = comp1->scaleLength;
      real_0 mass_b = comp1->originmass;
      real_0 xi_scale = comp1->scaleLength / (comp1->scaleLength + comp2->scaleLength);
      real_0 xi_mass = comp1->originmass / (comp1->originmass + comp2->originmass);
      Dwarf shiftedDwarf1, shiftedDwarf2;
      real_0 numBins = 10.0;
      real_0 integral;

      real_0 a = 0.0;
      real_0 b = v->value;
      real_0 width = (b-a)/numBins;
      real_0 a_tmp;
      real_0 b_tmp;
      real_0 c1 = 0.55555555555; //5.0 / 9.0;
      real_0 c2 = 0.88888888888; //8.0 / 9.0;
      real_0 c3 = 0.55555555555; //5.0 / 9.0;
      real_0 x1 = -0.77459666924;//-sqrt(3.0 / 5.0);
      real_0 x2 = 0.00000000000;
      real_0 x3 = 0.77459666924; //sqrt(3.0 / 5.0);
      real_0 v1, v2, v3, integrand1, integrand2, integrand3;

      switch(i)
      {
          case(BARYON_RADIUS_POS):
              scale_b += H_STEPSIZE;
          case(RADIUS_RATIO_POS):
              xi_scale += H_STEPSIZE;
          case(BARYON_MASS_POS):
              mass_b += H_STEPSIZE;
          case(MASS_RATIO_POS):
              xi_mass += H_STEPSIZE;
      }

      shiftedDwarf1 = getShiftedDwarf(comp1, scale_b, xi_scale, mass_b, xi_mass, TRUE);
      shiftedDwarf2 = getShiftedDwarf(comp2, scale_b, xi_scale, mass_b, xi_mass, FALSE);
      integral = 0;
      for(int k = 0; k < numBins; k++)
      {
          a_tmp = k*width;
          b_tmp = (k+1)*width;
          v1 = (b_tmp-a_tmp)*x1/2.0 + (b_tmp+a_tmp)/2.0;
          v2 = (b_tmp-a_tmp)*x2/2.0 + (b_tmp+a_tmp)/2.0;
          v3 = (b_tmp-a_tmp)*x3/2.0 + (b_tmp+a_tmp)/2.0;

          integrand1 = gradIntegrand_vel(v1, r, comp1, comp2, &shiftedDwarf1, &shiftedDwarf2, i);
          integrand2 = gradIntegrand_vel(v2, r, comp1, comp2, &shiftedDwarf1, &shiftedDwarf2, i);
          integrand3 = gradIntegrand_vel(v3, r, comp1, comp2, &shiftedDwarf1, &shiftedDwarf2, i);

          integral += (b_tmp-a_tmp)/2.0*c1*( integrand1 )
                    + (b_tmp-a_tmp)/2.0*c2*( integrand2 )
                    + (b_tmp-a_tmp)/2.0*c3*( integrand3 );
      }
      integral = -integral / vel_dist(v->value, r->value, comp1, comp2);
      return integral;
  }

  static inline void getGradientInfo_vel(real* v, real* r, const Dwarf* comp1, const Dwarf* comp2)
  {
      real_0 integral;

      //a_b derivative
      integral = getGradIntegral_vel(v, r, comp1, comp2, BARYON_RADIUS_POS);
      setRealGradient(v, integral, BARYON_RADIUS_POS);

      //xi_R derivative
      integral = getGradIntegral_vel(v, r, comp1, comp2, RADIUS_RATIO_POS);
      setRealGradient(v, integral, RADIUS_RATIO_POS);

      //M_b derivative
      integral = getGradIntegral_vel(v, r, comp1, comp2, BARYON_MASS_POS);
      setRealGradient(v, integral, BARYON_MASS_POS);

      //xi_M derivative
      integral = getGradIntegral_vel(v, r, comp1, comp2, MASS_RATIO_POS);
      setRealGradient(v, integral, MASS_RATIO_POS);
  }

  static inline real_0 hessIntegrand_vel(real_0 v, real* r, const Dwarf* comp1, const Dwarf* comp2, const Dwarf* shiftedDwarf1_i, const Dwarf* shiftedDwarf1_j, const Dwarf* shiftedDwarf1_d, const Dwarf* shiftedDwarf2_i, const Dwarf* shiftedDwarf2_j, const Dwarf* shiftedDwarf2_d, int i, int j)
  {
      //r has derivative information, so we also apply shifts to r
      real_0 shifted_r_i = r->value + r->gradient[i]*H_STEPSIZE;
      real_0 shifted_r_j = r->value + r->gradient[j]*H_STEPSIZE;
      real_0 shifted_r_d = r->value + r->gradient[i]*H_STEPSIZE + r->gradient[j]*H_STEPSIZE;

      real_0 integrand =  (  mw_log_0(vel_dist(v, shifted_r_d, shiftedDwarf1_d, shiftedDwarf2_d))
                           - mw_log_0(vel_dist(v, shifted_r_i, shiftedDwarf1_i, shiftedDwarf2_i))
                           - mw_log_0(vel_dist(v, shifted_r_j, shiftedDwarf1_j, shiftedDwarf2_j))
                           + mw_log_0(vel_dist(v, r->value, comp1, comp2))                                    )/H_STEPSIZE/H_STEPSIZE;

      return integrand;
  }

  static inline real_0 getHessIntegral_vel(real* v, real* r, const Dwarf* comp1, const Dwarf* comp2, int i, int j)
  {
      real_0 scale_b = comp1->scaleLength;
      real_0 mass_b = comp1->originmass;
      real_0 xi_scale = comp1->scaleLength / (comp1->scaleLength + comp2->scaleLength);
      real_0 xi_mass = comp1->originmass / (comp1->originmass + comp2->originmass);
      Dwarf shiftedDwarf1_i, shiftedDwarf1_j, shiftedDwarf1_d, shiftedDwarf2_i, shiftedDwarf2_j, shiftedDwarf2_d;
      real_0 integral;
      real_0 numBins = 10.0;

      real_0 a = 0.0;
      real_0 b = v->value;
      real_0 width = (b-a)/numBins;
      real_0 a_tmp;
      real_0 b_tmp;
      real_0 c1 = 0.55555555555; //5.0 / 9.0;
      real_0 c2 = 0.88888888888; //8.0 / 9.0;
      real_0 c3 = 0.55555555555; //5.0 / 9.0;
      real_0 x1 = -0.77459666924;//-sqrt(3.0 / 5.0);
      real_0 x2 = 0.00000000000;
      real_0 x3 = 0.77459666924; //sqrt(3.0 / 5.0);
      real_0 v1, v2, v3, integrand1, integrand2, integrand3;

      real_0 scale_b_i = scale_b;
      real_0 scale_b_j = scale_b;
      real_0 scale_b_d = scale_b;

      real_0 mass_b_i = mass_b;
      real_0 mass_b_j = mass_b;
      real_0 mass_b_d = mass_b;

      real_0 xi_scale_i = xi_scale;
      real_0 xi_scale_j = xi_scale;
      real_0 xi_scale_d = xi_scale;

      real_0 xi_mass_i = xi_mass;
      real_0 xi_mass_j = xi_mass;
      real_0 xi_mass_d = xi_mass;

      switch(i)
      {
          case(BARYON_RADIUS_POS):
              scale_b_i += H_STEPSIZE;
              scale_b_d += H_STEPSIZE;
          case(RADIUS_RATIO_POS):
              xi_scale_i += H_STEPSIZE;
              xi_scale_d += H_STEPSIZE;
          case(BARYON_MASS_POS):
              mass_b_i += H_STEPSIZE;
              mass_b_d += H_STEPSIZE;
          case(MASS_RATIO_POS):
              xi_mass_i += H_STEPSIZE;
              xi_mass_d += H_STEPSIZE;
      }

      switch(j)
      {
          case(BARYON_RADIUS_POS):
              scale_b_j += H_STEPSIZE;
              scale_b_d += H_STEPSIZE;
          case(RADIUS_RATIO_POS):
              xi_scale_j += H_STEPSIZE;
              xi_scale_d += H_STEPSIZE;
          case(BARYON_MASS_POS):
              mass_b_j += H_STEPSIZE;
              mass_b_d += H_STEPSIZE;
          case(MASS_RATIO_POS):
              xi_mass_j += H_STEPSIZE;
              xi_mass_d += H_STEPSIZE;
      }

      shiftedDwarf1_i = getShiftedDwarf(comp1, scale_b_i, xi_scale_i, mass_b_i, xi_mass_i, TRUE);
      shiftedDwarf1_j = getShiftedDwarf(comp1, scale_b_j, xi_scale_j, mass_b_j, xi_mass_j, TRUE);
      shiftedDwarf1_d = getShiftedDwarf(comp1, scale_b_d, xi_scale_d, mass_b_d, xi_mass_d, TRUE);
      shiftedDwarf2_i = getShiftedDwarf(comp2, scale_b_i, xi_scale_i, mass_b_i, xi_mass_i, FALSE);
      shiftedDwarf2_j = getShiftedDwarf(comp2, scale_b_j, xi_scale_j, mass_b_j, xi_mass_j, FALSE);
      shiftedDwarf2_d = getShiftedDwarf(comp2, scale_b_d, xi_scale_d, mass_b_d, xi_mass_d, FALSE);

      real_0 df_dx_i = gradIntegrand_vel(v->value, r, comp1, comp2, &shiftedDwarf1_i, &shiftedDwarf2_i, i);
      real_0 df_dx_j = gradIntegrand_vel(v->value, r, comp1, comp2, &shiftedDwarf1_j, &shiftedDwarf2_j, j);
      real_0 df_de = (vel_dist(v->value + H_STEPSIZE, r->value, comp1, comp2) - vel_dist(v->value, r->value, comp1, comp2))/H_STEPSIZE;

      integral = 0.0;
      for(int k = 0; k < numBins; k++)
      {
          a_tmp = k*width;
          b_tmp = (k+1)*width;
          v1 = (b_tmp-a_tmp)*x1/2.0 + (b_tmp+a_tmp)/2.0;
          v2 = (b_tmp-a_tmp)*x2/2.0 + (b_tmp+a_tmp)/2.0;
          v3 = (b_tmp-a_tmp)*x3/2.0 + (b_tmp+a_tmp)/2.0;

          integrand1 = hessIntegrand_vel(v1, r, comp1, comp2, &shiftedDwarf1_i, &shiftedDwarf1_j, &shiftedDwarf1_d, &shiftedDwarf2_i, &shiftedDwarf2_j, &shiftedDwarf2_d, i, j);
          integrand2 = hessIntegrand_vel(v2, r, comp1, comp2, &shiftedDwarf1_i, &shiftedDwarf1_j, &shiftedDwarf1_d, &shiftedDwarf2_i, &shiftedDwarf2_j, &shiftedDwarf2_d, i, j);
          integrand3 = hessIntegrand_vel(v3, r, comp1, comp2, &shiftedDwarf1_i, &shiftedDwarf1_j, &shiftedDwarf1_d, &shiftedDwarf2_i, &shiftedDwarf2_j, &shiftedDwarf2_d, i, j);

          integral += (b_tmp-a_tmp)/2.0*c1*( integrand1 )
                    + (b_tmp-a_tmp)/2.0*c2*( integrand2 )
                    + (b_tmp-a_tmp)/2.0*c3*( integrand3 );
      }

      integral += df_dx_i * v->gradient[j] + df_dx_j * v->gradient[i] + df_de * v->gradient[i] * v->gradient[j]; // Other terms of the Hessian
      integral = -integral / vel_dist(v->value, r->value, comp1, comp2);
      return integral;
  }

  static inline void getHessianInfo_vel(real* v, real* r, const Dwarf* comp1, const Dwarf* comp2)
  {
      real_0 integral;

      //hessian (a_b, a_b)
      integral = getHessIntegral_vel(v, r, comp1, comp2, BARYON_RADIUS_POS, BARYON_RADIUS_POS);
      setRealHessian(v, integral, BARYON_RADIUS_POS, BARYON_RADIUS_POS);

      //hessian (a_b, xi_R)
      integral = getHessIntegral_vel(v, r, comp1, comp2, BARYON_RADIUS_POS, RADIUS_RATIO_POS);
      setRealHessian(v, integral, BARYON_RADIUS_POS, RADIUS_RATIO_POS);

      //hessian (a_b, M_b)
      integral = getHessIntegral_vel(v, r, comp1, comp2, BARYON_RADIUS_POS, BARYON_MASS_POS);
      setRealHessian(v, integral, BARYON_RADIUS_POS, BARYON_MASS_POS);

      //hessian (a_b, xi_M)
      integral = getHessIntegral_vel(v, r, comp1, comp2, BARYON_RADIUS_POS, MASS_RATIO_POS);
      setRealHessian(v, integral, BARYON_RADIUS_POS, MASS_RATIO_POS);

      //hessian (xi_R, xi_R)
      integral = getHessIntegral_vel(v, r, comp1, comp2, RADIUS_RATIO_POS, RADIUS_RATIO_POS);
      setRealHessian(v, integral, RADIUS_RATIO_POS, RADIUS_RATIO_POS);

      //hessian (xi_R, M_b)
      integral = getHessIntegral_vel(v, r, comp1, comp2, RADIUS_RATIO_POS, BARYON_MASS_POS);
      setRealHessian(v, integral, RADIUS_RATIO_POS, BARYON_MASS_POS);

      //hessian (xi_R, xi_M)
      integral = getHessIntegral_vel(v, r, comp1, comp2, RADIUS_RATIO_POS, MASS_RATIO_POS);
      setRealHessian(v, integral, RADIUS_RATIO_POS, MASS_RATIO_POS);

      //hessian (M_b, M_b)
      integral = getHessIntegral_vel(v, r, comp1, comp2, BARYON_MASS_POS, BARYON_MASS_POS);
      setRealHessian(v, integral, BARYON_MASS_POS, BARYON_MASS_POS);

      //hessian (M_b, xi_M)
      integral = getHessIntegral_vel(v, r, comp1, comp2, BARYON_MASS_POS, MASS_RATIO_POS);
      setRealHessian(v, integral, BARYON_MASS_POS, MASS_RATIO_POS);

      //hessian (xi_M, xi_M)
      integral = getHessIntegral_vel(v, r, comp1, comp2, MASS_RATIO_POS, MASS_RATIO_POS);
      setRealHessian(v, integral, MASS_RATIO_POS, MASS_RATIO_POS);  
  }

  static inline void getDwarfDerivativeInfo_vel(real* v, real* r, const Dwarf* comp1, const Dwarf* comp2)
  {
      //Calculate gradient
      getGradientInfo_vel(v, r, comp1, comp2);

      //Calculate hessian
      getHessianInfo_vel(v, r, comp1, comp2);
  }

#endif //AUTODIFF

/*      SAMPLING FUNCTIONS      */
static inline real r_mag(dsfmt_t* dsfmtState, const Dwarf* comp1, const Dwarf* comp2, mwbool isLight)
{
    int counter = 0;
    real_0 r, u, val, bound1, bound2, rho_max;
    /*finding the max of the individual components*/
    real_0 rho_max_light = 0;
    real_0 rho_max_dark  = 0;

    switch(comp1->type)
    {
        case Plummer:
            bound1 =  50.0 * (comp1->scaleLength + comp2->scaleLength);
            rho_max_light = mw_sqrt_0(2.0 / 3.0) * comp1->scaleLength;
            break;
        case NFW:
            bound1 = 5.0 * comp1->r200;
            get_extra_nfw_mass(comp1, bound1);
            rho_max_light = comp1->scaleLength;
            break;
        case General_Hernquist:
            bound1 =  50.0 * (comp1->scaleLength + comp2->scaleLength);
            rho_max_light = comp1->scaleLength / 2.0;
            break;
    }

    switch(comp2->type)
    {
        case Plummer:
            bound2 =  50.0 * (comp1->scaleLength + comp2->scaleLength);
            rho_max_dark = mw_sqrt_0(2.0 / 3.0) * comp2->scaleLength;
            break;
        case NFW:
            bound2 = 5.0 * comp2->r200;
            get_extra_nfw_mass(comp2, bound2);
            rho_max_dark = comp2->scaleLength;
            break;
        case General_Hernquist:
            bound2 =  50.0 * (comp1->scaleLength + comp2->scaleLength);
            rho_max_dark = comp2->scaleLength / 2.0;
            break;
    }

    rho_max_light = 4 * M_PI * sqr_0(rho_max_light) * get_density(comp1, rho_max_light) / comp1->mass;
    rho_max_dark  = 4 * M_PI * sqr_0(rho_max_dark)  * get_density(comp2, rho_max_dark) / comp2->mass;
    
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
        if(isLight)
        {
            r = (real_0)mwXrandom(dsfmtState, 0.0, 1.0) * bound1;
            rho_max = rho_max_light;
        }
        else
        {
            r = (real_0)mwXrandom(dsfmtState, 0.0, 1.0) * bound2;
            rho_max = rho_max_dark;
        }
        u = (real_0)mwXrandom(dsfmtState, 0.0, 1.0);
        val = radial_dist(r, comp1, comp2, isLight);

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
    real r_val = mw_real_const(r);

    // Report radius
    return r_val;
}

static inline real vel_mag(real* r, const Dwarf* comp1, const Dwarf* comp2, dsfmt_t* dsfmtState)
{    
    int counter = 0;
    real_0 v, u, d;
    
    /* having the upper limit as exactly v_esc is bad since the dist fun seems to blow up there for small r. */
    real_0 v_esc = 0.99 * mw_sqrt_0( mw_fabs_0(2.0 * potential( showRealValue(r), comp1, comp2) ) );
    
    real_0 dist_max = max_finder(dist_fun, showRealValue(r), comp1, comp2, 0.0, 0.5 * v_esc, v_esc, 10, 1.0e-2);
    while(1)
    {

        v = (real_0)mwXrandom(dsfmtState, 0.0, 1.0) * v_esc;
        u = (real_0)mwXrandom(dsfmtState, 0.0, 1.0);

        d = dist_fun(v, showRealValue(r), comp1, comp2);
        

        if(mw_fabs_0(d / dist_max) > u)
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
    real v_val = mw_real_const(v);

    //Report veloctiy
    return v_val;
}

static inline mwvector get_components(dsfmt_t* dsfmtState, real* rad)
{
    /* assigns angles. Allows for non-circular orbits.*/
    /* have to sample in this way because sampling angles and then converting
     * to xyz leads to strong dependence on the rad, which could lead to bunching 
     * at the poles.
     */
    real r_sq, r_scaling, tmp;
    mwvector vec;

    do                                       
    {
        vec = mwRandomUnitPoint(dsfmtState); /* pick point in NDIM-space */
        r_sq = mw_sqrv(&vec);                 /* compute radius squared */
    }
    while (showRealValue(&r_sq) > 1.0);                      /* reject if outside unit sphere */

    tmp = mw_sqrt(&r_sq);
    r_scaling = mw_div(rad, &tmp);         /* compute scaling factor */
    vec.x = mw_mul(&vec.x, &r_scaling);    /* rescale to radius given */
    vec.y = mw_mul(&vec.y, &r_scaling);
    vec.z = mw_mul(&vec.z, &r_scaling);
    
    /* this is r * (u_vec / |u|). 
     * the r gives the magnitude, rad.
     * u_vec, which is the unit point original picked, vec, 
     * divided by the magnitude |u|, which is sqrt(r_sq),
     * gives it a direction (unit vector).
     */
    return vec;
}

static int cm_correction_by_comp(real * x, real * y, real * z, real * vx, real * vy, real * vz, real * mass, 
								mwvector* rShift, mwvector* vShift, real_0 comp_mass, unsigned int compStart, unsigned int compEnd)
{
    /*  
     * This function takes the table of bodies produced and zeroes the center of mass 
     * and center of momentum. It then shifts the center of mass and center of momentum
     * to the expected value for its position in the orbit.
     */
    real cm_x = ZERO_REAL;
    real cm_y = ZERO_REAL;
    real cm_z = ZERO_REAL;
    real cm_vx = ZERO_REAL;
    real cm_vy = ZERO_REAL;
    real cm_vz = ZERO_REAL;
    real tmp;
    unsigned int i;
    for(i = compStart; i < compEnd; i++)
    {
        cm_x = mw_mad(&mass[i], &x[i], &cm_x);
        cm_y = mw_mad(&mass[i], &y[i], &cm_y);
        cm_z = mw_mad(&mass[i], &z[i], &cm_z);
        
        cm_vx = mw_mad(&mass[i], &vx[i], &cm_vx);
        cm_vy = mw_mad(&mass[i], &vy[i], &cm_vy);
        cm_vz = mw_mad(&mass[i], &vz[i], &cm_vz);
    }
     
    cm_x = mw_mul_s(&cm_x, inv_0(comp_mass));
    cm_y = mw_mul_s(&cm_y, inv_0(comp_mass));
    cm_z = mw_mul_s(&cm_z, inv_0(comp_mass));
    
    cm_vx = mw_mul_s(&cm_vx, inv_0(comp_mass));
    cm_vy = mw_mul_s(&cm_vy, inv_0(comp_mass));
    cm_vz = mw_mul_s(&cm_vz, inv_0(comp_mass));

    for(i = compStart; i < compEnd; i++)
    {
        tmp = mw_sub(&x[i], &cm_x);
        x[i] = mw_add(&tmp, &rShift->x);
        tmp = mw_sub(&y[i], &cm_y);
        y[i] = mw_add(&tmp, &rShift->y);
        tmp = mw_sub(&z[i], &cm_z);
        z[i] = mw_add(&tmp, &rShift->z);
        
        tmp = mw_sub(&vx[i], &cm_vx);
        vx[i] = mw_add(&tmp, &vShift->x);
        tmp = mw_sub(&vy[i], &cm_vy);
        vy[i] = mw_add(&tmp, &vShift->y);
        tmp = mw_sub(&vz[i], &cm_vz);
        vz[i] = mw_add(&tmp, &vShift->z);
    }
    return 1;
}


/*      DWARF GENERATION        */
static int nbGenerateMixedDwarfCore(lua_State* luaSt, dsfmt_t* prng, unsigned int nbody, 
                                     Dwarf* comp1,  Dwarf* comp2, 
                                    mwbool ignore, mwvector* rShift, mwvector* vShift)
{
    /* generatePlummer: generate Plummer model initial conditions for test
    * runs, scaled to units such that M = -4E = G = 1 (Henon, Heggie,
    * etc).    See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37,
    * 183.
    */
        mw_printf("Generating Mixed Dwarf...\n");
        unsigned int i;
        int table;
        Body b;
        real r, v, tmp;
 
        real * x  = mwCalloc(nbody, sizeof(real));
        real * y  = mwCalloc(nbody, sizeof(real));
        real * z  = mwCalloc(nbody, sizeof(real));
        real * vx = mwCalloc(nbody, sizeof(real));
        real * vy = mwCalloc(nbody, sizeof(real));
        real * vz = mwCalloc(nbody, sizeof(real));
        real * masses = mwCalloc(nbody, sizeof(real));
        
        
        mwvector vec;
        real_0 rscale_l = comp1->scaleLength; //comp1[1]; /*scale radius of the light component*/
        real_0 rscale_d = comp2->scaleLength; //comp2[1]; /*scale radius of the dark component*/
        set_p0(comp1);
        set_p0(comp2);
        real_0 bound1 ;
        real_0 bound2 ;
        
        
        real_0 mass_l   = comp1->mass; //comp1[0]; /*mass of the light component*/
        real_0 mass_d   = comp2->mass; //comp2[0]; /*mass of the dark component*/
        real_0 dwarf_mass = mass_l + mass_d;

        real mass_light = mw_real_var(mass_l, BARYON_MASS_POS);            /*Setting dwarf light mass to position 4 of gradient*/
        real mass_xi    = mw_real_var(mass_l/dwarf_mass, MASS_RATIO_POS);  /*Setting dwarf mass ratio to position 5 of gradient*/
        tmp = inv(&mass_xi);
        tmp = mw_add_s(&tmp, -1.0);
        real mass_dark  = mw_mul(&mass_light, &tmp);


    //---------------------------------------------------------------------------------------------------        
        unsigned int half_bodies = nbody / 2;
        real mass_light_particle = mw_mul_s(&mass_light, inv_0((real_0) half_bodies));//half the particles are light matter
        real mass_dark_particle  = mw_mul_s(&mass_dark , inv_0((real_0) half_bodies));
    //----------------------------------------------------------------------------------------------------

	
        /* dark matter type is TRUE or 1. Light matter type is False, or 0*/
        mwbool isdark = TRUE;
        mwbool islight = FALSE;
        
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
                    r = r_mag(prng, comp1, comp2, TRUE);
                    masses[i] = mass_light_particle;
                }
                else if(i >= half_bodies)
                {
                    r = r_mag(prng, comp1, comp2, FALSE);
                    masses[i] = mass_dark_particle;
                }
                /*to ensure that r is finite and nonzero*/
                if(isinf(showRealValue(&r)) == FALSE && showRealValue(&r) != 0.0 && isnan(showRealValue(&r)) == FALSE){break;}
                
                if(counter > 1000)
                {
                    exit(-1);
                }
                else
                {
                    counter++;
                }
                
            }while (1);
            //printReal(&r,"Body Rad");

            counter = 0;
            do
            {
                v = vel_mag(&r, comp1, comp2, prng);
                //mw_printf("VELOCTIY ATTEMPT\n");
                if(isinf(showRealValue(&v)) == FALSE && showRealValue(&v) != 0.0 && isnan(showRealValue(&v)) == FALSE){break;}
                
                if(counter > 1000)
                {
                    exit(-1);
                }
                else
                {
                    counter++;
                }
                
            }while (1);
            //printReal(&v,"Body Vel");

            //mw_printf("BODY VEL = %.15f\n", showRealValue(&v));
            vec = get_components(prng, &v); 
            vx[i] = vec.x;
            vy[i] = vec.y;
            vz[i] = vec.z;
            //mw_printf("BODY VEL = [%.15f, %.15f, %.15f]\n", showRealValue(&vec.x), showRealValue(&vec.y), showRealValue(&vec.z));

            //mw_printf("BODY POS = %.15f\n", showRealValue(&r));
            vec = get_components(prng, &r);  
            x[i] = vec.x;
            y[i] = vec.y;
            z[i] = vec.z;
            //mw_printf("BODY POS = [%.15f, %.15f, %.15f]\n", showRealValue(&vec.x), showRealValue(&vec.y), showRealValue(&vec.z));
        }

      #if AUTODIFF // Calculating the initial dwarf derivatives is time consuming, so we are multithreading this process alone
          mw_printf("    Calculating initial dwarf derivatives...\n");
          real r_val, v_val;
          real_0 x_frac, y_frac, z_frac, vx_frac, vy_frac, vz_frac;
        #ifdef _OPENMP
          real_0 t = omp_get_wtime();
          #pragma omp parallel for private(i, r_val, v_val, x_frac, y_frac, z_frac, vx_frac, vy_frac, vz_frac) shared(x, y, z, vx, vy, vz, comp1, comp2, half_bodies) schedule(dynamic)
        #else
          clock_t t;
          t = clock();
        #endif
          for (i = 0; i < nbody; i++)
          {      
              // Get dwarf derivative information for radius     
              r_val = mw_hypot(&x[i], &y[i]);
              r_val = mw_hypot(&r_val, &z[i]);
              x_frac = x[i].value / r_val.value;
              y_frac = y[i].value / r_val.value;
              z_frac = z[i].value / r_val.value;
              if (i < half_bodies)
              {
                  getDwarfDerivativeInfo_rad(&r_val, comp1, comp2, TRUE);
              }
              else if(i >= half_bodies)
              {
                  getDwarfDerivativeInfo_rad(&r_val, comp1, comp2, FALSE);
              }
              x[i] = mw_mul_s(&r_val, x_frac);
              y[i] = mw_mul_s(&r_val, y_frac);
              z[i] = mw_mul_s(&r_val, z_frac);

              // Get dwarf derivative information for velocity     
              v_val = mw_hypot(&vx[i], &vy[i]);
              v_val = mw_hypot(&v_val, &vz[i]);
              vx_frac = vx[i].value / v_val.value;
              vy_frac = vy[i].value / v_val.value;
              vz_frac = vz[i].value / v_val.value;
              getDwarfDerivativeInfo_vel(&v_val, &r_val, comp1, comp2);
              vx[i] = mw_mul_s(&v_val, vx_frac);
              vy[i] = mw_mul_s(&v_val, vy_frac);
              vz[i] = mw_mul_s(&v_val, vz_frac);
              //mw_printf("      Particle %u logged\n", i+1);
              //printReal(&r_val, "Radius");
              //printReal(&v_val, "Velocity");
          }
        #ifdef _OPENMP
          real_0 time_taken = omp_get_wtime() - t;
        #else
          t = clock() - t;
          real_0 time_taken = ((real_0)t)/CLOCKS_PER_SEC; // in seconds
        #endif
          mw_printf("    Derivative calculation took %.3f seconds\n", time_taken);
      #endif


        /* getting the center of mass and momentum correction */
		cm_correction_by_comp(x, y, z, vx, vy, vz, masses, rShift, vShift, mass_l, 0, half_bodies); //corrects light component
		cm_correction_by_comp(x, y, z, vx, vy, vz, masses, rShift, vShift, mass_d, half_bodies, nbody); //corrects dark component
        //cm_correction(x, y, z, vx, vy, vz, masses, rShift, vShift, dwarf_mass, nbody);


        /* pushing the bodies */
        for (i = 0; i < nbody; i++)
        {
            b.bodynode.id   = i + 1;
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

            //mw_printf("BODY POS  = [%.15f, %.15f, %.15f]\n", showRealValue(&b.bodynode.pos.x), showRealValue(&b.bodynode.pos.y), showRealValue(&b.bodynode.pos.z));
            
            b.vel.x = vx[i];
            b.vel.y = vy[i];
            b.vel.z = vz[i];

            //mw_printf(" BODY VEL = [%.15f, %.15f, %.15f]\n", showRealValue(&b.vel.x), showRealValue(&b.vel.y), showRealValue(&b.vel.z));
            
            assert(nbPositionValid(&b.bodynode.pos));
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

        mw_printf("Mixed Dwarf Generated!\n");
        
        return 1;             
        
}


int nbGenerateMixedDwarfCore_TESTVER(mwvector* pos, mwvector* vel, real* bodyMasses, dsfmt_t* prng, unsigned int nbody, 
                                     Dwarf* comp1,  Dwarf* comp2, mwvector* rShift, mwvector* vShift)
{
    /* NOTE: unction is designed to mimic the above function, but bypass the need for the
	* for the lua state. It is used in the test unit for constructing multi-component
	* dwarf galaxies. Any changes to the above function not pertaining to pushing the bodies
	* should be made here
    */
        unsigned int i;
        int table;
        Body b;
        real r, v, tmp;
 
        real * x  = mwCalloc(nbody, sizeof(real));
        real * y  = mwCalloc(nbody, sizeof(real));
        real * z  = mwCalloc(nbody, sizeof(real));
        real * vx = mwCalloc(nbody, sizeof(real));
        real * vy = mwCalloc(nbody, sizeof(real));
        real * vz = mwCalloc(nbody, sizeof(real));
        real * masses = mwCalloc(nbody, sizeof(real));
        
        
        mwvector vec;
        real_0 rscale_l = comp1->scaleLength; //comp1[1]; /*scale radius of the light component*/
        real_0 rscale_d = comp2->scaleLength; //comp2[1]; /*scale radius of the dark component*/
        set_p0(comp1);
        set_p0(comp2);
        real_0 bound1 ;
        real_0 bound2 ;        
        
        real_0 mass_l   = comp1->mass; //comp1[0]; /*mass of the light component*/
        real_0 mass_d   = comp2->mass; //comp2[0]; /*mass of the dark component*/
        real_0 dwarf_mass = mass_l + mass_d;

        real mass_light = mw_real_var(mass_l, BARYON_MASS_POS);            /*Setting dwarf light mass to position 4 of gradient*/
        real mass_xi    = mw_real_var(mass_l/dwarf_mass, MASS_RATIO_POS); /*Setting dwarf mass ratio to position 5 of gradient*/
        tmp = inv(&mass_xi);
        tmp = mw_add_s(&tmp, -1.0);
        real mass_dark  = mw_mul(&mass_light, &tmp);


    //---------------------------------------------------------------------------------------------------        
        unsigned int half_bodies = nbody / 2;
        real mass_light_particle = mw_mul_s(&mass_light, inv_0((real_0) half_bodies));//half the particles are light matter
        real mass_dark_particle  = mw_mul_s(&mass_dark , inv_0((real_0) half_bodies));
    //----------------------------------------------------------------------------------------------------

	
        /* dark matter type is TRUE or 1. Light matter type is False, or 0*/
        mwbool isdark = TRUE;
        mwbool islight = FALSE;
        
    	/*initializing particles:*/
        //memset(&b, 0, sizeof(b));
        //lua_createtable(luaSt, nbody, 0);
        //table = lua_gettop(luaSt);      
        int counter = 0;
        

        /*getting the radii and velocities for the bodies*/
        for (i = 0; i < nbody; i++)
        {
            counter = 0;
            do
            {
                
                if(i < half_bodies)
                {
                    r = r_mag(prng, comp1, comp2, TRUE);
                    masses[i] = mass_light_particle;
                }
                else if(i >= half_bodies)
                {
                    r = r_mag(prng, comp1, comp2, FALSE);
                    masses[i] = mass_dark_particle;
                }
                /*to ensure that r is finite and nonzero*/
                if(isinf(showRealValue(&r)) == FALSE && showRealValue(&r) != 0.0 && isnan(showRealValue(&r)) == FALSE){break;}
                
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
                v = vel_mag(&r, comp1, comp2, prng);
                if(isinf(showRealValue(&v)) == FALSE && showRealValue(&v) != 0.0 && isnan(showRealValue(&v)) == FALSE){break;}
                
                if(counter > 1000)
                {
                    exit(-1);
                }
                else
                {
                    counter++;
                }
                
            }while (1);
			
            vec   = get_components(prng, &v);   
            vx[i] = vec.x;
            vy[i] = vec.y;
            vz[i] = vec.z;
            vec   = get_components(prng, &r);  
            x[i] = vec.x;
            y[i] = vec.y;
            z[i] = vec.z;
        }


        /* getting the center of mass and momentum correction */
		cm_correction_by_comp(x, y, z, vx, vy, vz, masses, rShift, vShift, mass_l, 0, half_bodies); //corrects light component
		cm_correction_by_comp(x, y, z, vx, vy, vz, masses, rShift, vShift, mass_d, half_bodies, nbody); //corrects dark component
        //cm_correction(x, y, z, vx, vy, vz, masses, rShift, vShift, dwarf_mass, nbody);


        /* pushing the bodies */
        for (i = 0; i < nbody; i++)
        {

			pos[i].x  = x[i];
			pos[i].y  = y[i];
			pos[i].z  = z[i];

			vel[i].x = vx[i];
			vel[i].y = vy[i];
			vel[i].z = vz[i];

			bodyMasses[i] = masses[i];
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
        static real_0 nbodyf = 0.0;
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
                                                                 position, velocity);
}


void registerGenerateMixedDwarf(lua_State* luaSt)
{
    lua_register(luaSt, "generatemixeddwarf", nbGenerateMixedDwarf);
}


// As this code runs, know that it is running on the rotting corpses of a thousand bugs.
